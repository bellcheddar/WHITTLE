"""Interface metrics for the co-folded complexes. Spec 6.4 items 4, 5 and 6.

For each complex:

* **Buried surface area**, computed with FreeSASA as the difference between the
  complex and its separated chains, halved. Reported **per subunit and in
  total** for the TNFa and TL1A trimers, so a candidate binding across the
  subunit cleft is distinguishable from one binding a single face (spec 5.2).
* **Interface residues**, any antigen residue with a heavy atom within
  `CONTACT_CUTOFF` of a nanobody heavy atom, and the antigen chains they sit on.
* **PRODIGY interface delta G** (`prodigy-prot`, Apache-2.0), where it is
  installed. Absent rather than faked if it is not.
* **Epitope binning** across candidates for a target: cluster by interface
  residue set, so "these six candidates hit the same face" is visible rather
  than implied.

Confidences (ipTM, pTM, pLDDT, PAE) come straight out of boltz's own
`confidence_*.json`; nothing here recomputes them.

    pixi run python -m pipeline.interface
"""
from __future__ import annotations

import argparse
import json
import math
import re
import sys
from collections import defaultdict
from pathlib import Path

from pipeline.common import (RUNS, STRUCTURES, TARGETS, connect, log,
                             manifest_read, manifest_write, targets_config)

CONTACT_CUTOFF = 5.0          # heavy-atom distance, angstroms
EPITOPE_JACCARD = 0.5         # two candidates share a bin above this overlap
NANOBODY_CHAIN = "N"
# A construct is a cut out of a larger protein, so its first and last residues
# are a knife edge rather than a surface. Anything within this many residues of
# either terminus is on that edge.
EDGE_RESIDUES = 6
# Above this share of the interface sitting on cut edges, the epitope is an
# artefact of the construct and could not exist in the intact protein.
EDGE_ARTEFACT_FRACTION = 0.30


# --------------------------------------------------------------- parsing

def read_complex(path: Path) -> dict[str, list[tuple]]:
    """{chain: [(resnum, resname, x, y, z), ...]} for heavy atoms."""
    chains: dict[str, list[tuple]] = defaultdict(list)
    for line in path.read_text().splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        element = (line[76:78].strip() or line[12:16].strip()[:1])
        if element == "H":
            continue
        chain = line[21].strip() or "A"
        chains[chain].append((int(line[22:26]), line[17:20].strip(),
                              float(line[30:38]), float(line[38:46]),
                              float(line[46:54])))
    return dict(chains)


def chain_plddt(path: Path, antigen_chains: set[str] | list[str]
                ) -> tuple[float | None, float | None]:
    """Mean CA pLDDT for the nanobody chain and for the antigen chains.

    boltz stores pLDDT in the B-factor column and writes no per-chain summary,
    so the split is read off the coordinates. It is worth reading because
    `complex_plddt` is an average over the whole complex and the antigen is
    several times the nanobody's size, which dilutes anything that happens to
    the nanobody by roughly the ratio of their lengths. Pooled over the paired
    controls the whole-complex figure separates a candidate from the control
    shuffled from it by 0.018 and the nanobody chain alone by 5.1 pLDDT, and
    the antigen chain, which is the same sequence from the same MSA inside a
    pair, sits at chance. See `tools/plddt_locus.py` and DECISIONS 49.

    CA atoms only: every residue has exactly one, so the mean is per residue
    rather than weighted by how many atoms a side chain happens to have.
    """
    antigen = set(antigen_chains)
    nanobody: list[float] = []
    antigen_values: list[float] = []
    for line in path.read_text().splitlines():
        if not line.startswith("ATOM") or line[12:16].strip() != "CA":
            continue
        try:
            value = float(line[60:66])
        except ValueError:
            continue
        chain = line[21].strip() or "A"
        (antigen_values if chain in antigen else nanobody).append(value)
    mean = lambda xs: sum(xs) / len(xs) if xs else None  # noqa: E731
    return mean(nanobody), mean(antigen_values)


CONFIDENCE_KEYS = ("confidence_score", "ptm", "iptm", "complex_plddt",
                   "complex_iplddt", "complex_pde", "complex_ipde",
                   "ligand_iptm", "protein_iptm")


def confidence_for(model_path: Path) -> dict:
    """The confidence numbers for **this** diffusion sample, not an average.

    boltz writes `<stem>_model_<n>.pdb` beside
    `confidence_<stem>_model_<n>.json`, so each sample has its own numbers.
    Averaging them across samples and storing the same value on every row
    would throw away precisely the sampling variability the five samples exist
    to show, and the spread between samples is itself a confidence statement.
    """
    path = model_path.with_name(f"confidence_{model_path.stem}.json")
    if not path.exists():
        return {}
    try:
        payload = json.loads(path.read_text())
    except json.JSONDecodeError:
        return {}
    out = {k: round(float(payload[k]), 4) for k in CONFIDENCE_KEYS
           if isinstance(payload.get(k), (int, float))}
    # Per-chain and pairwise interface confidences, which is what the spec
    # means by "the per-interface confidences".
    pair = payload.get("pair_chains_iptm") or {}
    try:
        values = [float(v) for row in pair.values() for k, v in row.items()
                  if k != list(pair).index]
        out["pair_iptm_min"] = round(min(values), 4) if values else None
    except Exception:
        pass
    return out


MODEL_INDEX = re.compile(r"_model_(\d+)$")


def structures(run_dir: Path) -> list[Path]:
    """The diffusion samples for one complex, in sample order.

    The index is matched rather than split off the end, because the run
    directories live under iCloud-synced Documents and a directory rewritten
    while a sync is in flight leaves the superseded copy behind as
    `<stem>_model_4 2.pdb`. Those files are real PDBs from an earlier run of
    the same candidate, so a plain glob picks them up as extra samples, and
    `int("4 2")` then raises ValueError and takes the whole stage down. This
    cost the first interface pass of 2026-11-14: 250 co-folds finished and the
    refresh died on the first candidate that had one.

    A file whose index is not an integer is skipped rather than repaired. It
    is a duplicate of a sample this function is already returning, so dropping
    it loses nothing, and silently averaging a stale pose into a live run
    would be worse than either crashing or skipping.
    """
    indexed: list[tuple[int, Path]] = []
    for path in run_dir.rglob("*_model_*.pdb"):
        if "pre_affinity" in path.name:
            continue
        found = MODEL_INDEX.search(path.stem)
        if not found:
            continue
        indexed.append((int(found.group(1)), path))
    return [path for _, path in sorted(indexed, key=lambda pair: pair[0])]


# ------------------------------------------------------------------- metrics

def interface_residues(chains: dict[str, list[tuple]], antigen_chains: list[str]
                       ) -> dict:
    """Antigen residues contacting the nanobody, and which chains they are on.

    Also reports how much of the interface sits on the construct's cut edges.
    Every antigen here is a slice of a larger protein, so its first and last
    residues are a knife edge that does not exist as a surface in the real
    molecule: for GLP1R the N-terminal edge of the ECD is where the 7TM bundle
    continues, and a nanobody modelled binding there has found a hole in the
    construct rather than an epitope. The smoke co-fold put three of eighteen
    contacts on residues 1 to 3, which is what prompted the check.
    """
    nano = chains.get(NANOBODY_CHAIN) or []
    if not nano:
        return {}
    cutoff2 = CONTACT_CUTOFF ** 2
    contacts: dict[str, set[int]] = {c: set() for c in antigen_chains}
    nano_contacts: set[int] = set()
    for chain in antigen_chains:
        for resnum, _name, ax, ay, az in chains.get(chain, []):
            for nres, _n, bx, by, bz in nano:
                if (ax - bx) ** 2 + (ay - by) ** 2 + (az - bz) ** 2 <= cutoff2:
                    contacts[chain].add(resnum)
                    nano_contacts.add(nres)
                    break
    touched = [c for c, residues in contacts.items() if residues]

    # Cut-edge analysis, per chain, using each chain's own residue range.
    edge_hits = 0
    total_hits = 0
    for chain in antigen_chains:
        numbers = [resnum for resnum, *_ in chains.get(chain, [])]
        if not numbers:
            continue
        lo, hi = min(numbers), max(numbers)
        for resnum in contacts.get(chain, ()):  # already a set
            total_hits += 1
            if resnum - lo < EDGE_RESIDUES or hi - resnum < EDGE_RESIDUES:
                edge_hits += 1
    edge_fraction = (edge_hits / total_hits) if total_hits else 0.0

    return {
        "per_chain": {c: sorted(r) for c, r in contacts.items() if r},
        "n_antigen_residues": sum(len(r) for r in contacts.values()),
        # The paratope. Needed to ask whether the CDRs are doing the binding,
        # which the antigen side of the interface cannot answer on its own.
        "nanobody_residues": sorted(nano_contacts),
        "n_nanobody_residues": len(nano_contacts),
        "n_subunits_contacted": len(touched),
        "chains_contacted": touched,
        "n_edge_contacts": edge_hits,
        "edge_fraction": round(edge_fraction, 3),
        # A flag, not a score: the UI says so and nothing in the ranking reads
        # it. What it means is that this pose should not be believed as an
        # epitope on the intact protein.
        "edge_artefact": 1 if edge_fraction > EDGE_ARTEFACT_FRACTION else 0,
        # The identity of an epitope, for binning. Residue numbers only, with
        # the chain dropped: in a homotrimer the same face on chain A and on
        # chain B is the same epitope.
        "epitope_key": sorted({r for residues in contacts.values()
                               for r in residues}),
    }


_CONTROL_REGIONS: dict[str, dict] = {}


def _control_regions(row: dict) -> dict | None:
    """CDRs for a control, numbered on demand and remembered."""
    cid = row["id"]
    if cid in _CONTROL_REGIONS:
        return _CONTROL_REGIONS[cid]
    try:
        from pipeline import numbering as nb

        regions = nb.number_many([(cid, row["sequence"])])[cid].regions or {}
    except Exception:
        regions = {}
    out = {"sequence": row["sequence"], **{k: regions.get(k) for k in
                                           ("cdr1", "cdr2", "cdr3")}}
    _CONTROL_REGIONS[cid] = out
    return out


def cdr_spans(sequence: str, cdrs: dict[str, str]) -> dict[str, tuple[int, int]]:
    """Residue ranges of each CDR in chain N, 1-based and inclusive.

    Chain N is numbered 1..len(sequence) and carries the candidate sequence
    exactly, checked across the run, so a substring search is enough. The CDRs
    occur in order, so when a short CDR1 or CDR2 matches in more than one place
    the search continues from the end of the previous one. One candidate in
    2500 needs that tie-break.
    """
    spans: dict[str, tuple[int, int]] = {}
    cursor = 0
    for name in ("cdr1", "cdr2", "cdr3"):
        motif = (cdrs.get(name) or "").strip()
        if not motif:
            continue
        at = sequence.find(motif, cursor)
        if at < 0:                      # out of order: fall back to anywhere
            at = sequence.find(motif)
        if at < 0:
            continue
        spans[name] = (at + 1, at + len(motif))
        cursor = at + len(motif)
    return spans


def paratope_composition(nanobody_residues: list[int],
                         spans: dict[str, tuple[int, int]]) -> dict:
    """How much of the paratope is CDR rather than framework.

    A nanobody modelled binding through its framework is either a modelling
    artefact or a contact that would not be specific, because the framework is
    the conserved part. Scale-free by construction: a fraction, not a count, so
    it does not simply restate interface size the way the existing affinity
    components do.
    """
    if not nanobody_residues or not spans:
        return {}
    detail = {name: 0 for name in spans}
    detail["framework"] = 0
    for resnum in nanobody_residues:
        for name, (lo, hi) in spans.items():
            if lo <= resnum <= hi:
                detail[name] += 1
                break
        else:
            detail["framework"] += 1
    in_cdr = sum(v for k, v in detail.items() if k != "framework")
    return {
        "cdr_interface_fraction": round(in_cdr / len(nanobody_residues), 3),
        "cdr_interface_detail": detail,
    }


def buried_area(path: Path, antigen_chains: list[str]) -> dict:
    """BSA from FreeSASA: (complex - separated) / 2, total and per subunit."""
    import tempfile

    import freesasa

    freesasa.setVerbosity(freesasa.silent)

    def sasa_of(lines: list[str]) -> float:
        with tempfile.NamedTemporaryFile("w", suffix=".pdb", delete=False) as fh:
            fh.write("".join(lines))
            tmp = fh.name
        try:
            structure = freesasa.Structure(tmp)
            return float(freesasa.calc(structure).totalArea())
        except Exception:
            return float("nan")
        finally:
            Path(tmp).unlink(missing_ok=True)

    text = path.read_text().splitlines(keepends=True)
    atoms = [ln for ln in text if ln.startswith("ATOM")]
    nano = [ln for ln in atoms if ln[21] == NANOBODY_CHAIN]
    if not nano:
        return {}
    whole = sasa_of(atoms)
    sasa_nano = sasa_of(nano)
    antigen_all = [ln for ln in atoms if ln[21] in antigen_chains]
    sasa_antigen = sasa_of(antigen_all)
    if any(math.isnan(v) for v in (whole, sasa_nano, sasa_antigen)):
        return {}
    total = (sasa_nano + sasa_antigen - whole) / 2.0

    per_subunit: dict[str, float] = {}
    if len(antigen_chains) > 1:
        for chain in antigen_chains:
            subunit = [ln for ln in atoms if ln[21] == chain]
            if not subunit:
                continue
            pair = sasa_of(subunit + nano)
            value = (sasa_of(subunit) + sasa_nano - pair) / 2.0
            if not math.isnan(value) and value > 1.0:
                per_subunit[chain] = round(value, 1)
    return {"bsa_total": round(total, 1), "bsa_per_subunit": per_subunit}


def prodigy_dg(path: Path, antigen_chains: list[str]) -> float | None:
    """PRODIGY interface delta G, if prodigy-prot is installed.

    Prodigy takes a Bio.PDB *Model*, not a Structure, and its `selection` is a
    list of chain groups: the antigen chains as one comma-joined group against
    the nanobody chain. Checked against the installed 2.4.0 signature rather
    than assumed, because an older import path (`prodigy_prot.predict_IC`) is
    what most examples show and it no longer exists.
    """
    try:
        from prodigy_prot.modules.prodigy import Prodigy
    except Exception:
        return None
    try:
        from Bio.PDB import PDBParser

        parser = PDBParser(QUIET=True)
        model = parser.get_structure("c", str(path))[0]
        prodigy = Prodigy(model, selection=[",".join(antigen_chains),
                                            NANOBODY_CHAIN], temp=25.0)
        prodigy.predict()
        return round(float(prodigy.ba_val), 3)
    except Exception:
        return None


# ------------------------------------------------------------ epitope binning

def pose_consensus(samples: list[dict]) -> dict:
    """Do the independent diffusion samples agree about where this binds?

    Mean pairwise Jaccard over the epitope of each seed. This is the one
    interface measure that costs nothing and is genuinely orthogonal to
    interface size: a candidate can bury a lot of surface in a different place
    every time, and that is a pose the ranking should not trust. A single
    good-looking sample out of five is a lucky draw, not evidence.

    Returned with the ipTM spread, which says the same thing in the model's
    own confidence terms.
    """
    keys = [set(s.get("epitope_key") or ()) for s in samples]
    keys = [k for k in keys if k]
    out: dict[str, float | None] = {"pose_consensus": None, "iptm_spread": None}
    if len(keys) >= 2:
        scores = []
        for i in range(len(keys)):
            for j in range(i + 1, len(keys)):
                union = keys[i] | keys[j]
                scores.append(len(keys[i] & keys[j]) / len(union) if union else 0.0)
        if scores:
            out["pose_consensus"] = round(sum(scores) / len(scores), 3)
    iptms = [s["iptm"] for s in samples if s.get("iptm") is not None]
    if len(iptms) >= 2:
        out["iptm_spread"] = round(max(iptms) - min(iptms), 4)
    return out


def publish_complex(candidate_id: str, target: str, model: Path) -> str | None:
    """Copy the best co-folded complex where the drawer can load it.

    Tier 2's whole point is showing where a candidate binds, and the complexes
    live under `data/runs/`, which is deliberately not deployed: it holds every
    diffusion sample, the MSA and boltz's processed inputs, and is far too
    large to ship. Without this step the site would carry ipTM, interface
    energies and epitope bins with no picture to attach them to.

    One file per candidate, the highest-ipTM sample, gzipped beside the Tier 1
    model so the same nginx rule and the same deploy already handle it.
    """
    import gzip

    out = STRUCTURES / target / f"{candidate_id}.complex.pdb.gz"
    try:
        out.parent.mkdir(parents=True, exist_ok=True)
        raw = model.read_bytes()
        with gzip.open(out, "wb", compresslevel=9) as fh:
            fh.write(b"REMARK  CO-FOLDED COMPLEX: CHAIN N IS THE NANOBODY, "
                     b"A (AND B, C) THE ANTIGEN\n")
            fh.write(raw)
        return out.name
    except Exception:
        return None


def bin_epitopes(records: list[dict]) -> dict[str, int]:
    """Single-linkage clustering of candidates by interface residue overlap."""
    keys = {r["candidate_id"]: set(r.get("epitope_key") or []) for r in records}
    ids = [cid for cid, k in keys.items() if k]
    parent = {cid: cid for cid in ids}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for i, a in enumerate(ids):
        for b in ids[i + 1:]:
            union = keys[a] | keys[b]
            if not union:
                continue
            if len(keys[a] & keys[b]) / len(union) >= EPITOPE_JACCARD:
                ra, rb = find(a), find(b)
                if ra != rb:
                    parent[rb] = ra
    labels: dict[str, int] = {}
    numbering: dict[str, int] = {}
    for cid in ids:
        root = find(cid)
        if root not in numbering:
            numbering[root] = len(numbering) + 1
        labels[cid] = numbering[root]
    return labels


# ----------------------------------------------------------------------- main

def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--target", action="append", default=None)
    args = parser.parse_args()

    targets = args.target or TARGETS
    cfg = targets_config()["targets"]
    conn = connect()
    done_runs = {str(r["id"]): r for r in manifest_read("cofold")
                 if r.get("status") == "ok"}

    for target in targets:
        stoich = int(cfg[target]["stoichiometry"])
        antigen_chains = ["A", "B", "C"][:stoich]
        # CDR spans come from the candidate's own sequence: chain N is numbered
        # 1..len(sequence) and carries that sequence exactly.
        seqs = {r["id"]: dict(r) for r in conn.execute(
            "SELECT id, sequence, cdr1, cdr2, cdr3 FROM candidates WHERE target = ?",
            (target,))}
        # Scrambled controls go through this same code (spec 8.2). They carry
        # no CDR columns, so their regions are numbered here; the CDRs of a
        # shuffled sequence are still positionally defined, which is what the
        # paratope fraction needs.
        controls = {r["id"]: dict(r) for r in conn.execute(
            "SELECT id, sequence FROM controls WHERE target = ? "
            "AND sequence IS NOT NULL", (target,))}
        records: list[dict] = []
        for cid, run in done_runs.items():
            if run.get("target") != target:
                continue
            run_dir = Path(run["path"])
            if not run_dir.exists():
                continue
            models = structures(run_dir)
            if not models:
                manifest_write("interface", {"id": cid, "target": target,
                                             "status": "failed:no structure"})
                continue
            for sample, model in enumerate(models):
                conf = confidence_for(model)
                nb_plddt, ag_plddt = chain_plddt(model, antigen_chains)
                chains = read_complex(model)
                iface = interface_residues(chains, antigen_chains)
                areas = buried_area(model, antigen_chains)
                record = {
                    "candidate_id": cid, "target": target, "seed": sample,
                    "iptm": conf.get("iptm"), "ptm": conf.get("ptm"),
                    "plddt": conf.get("complex_plddt"),
                    "nanobody_plddt": nb_plddt,
                    "antigen_plddt": ag_plddt,
                    "interface_plddt": conf.get("complex_iplddt"),
                    "pae": conf.get("complex_pde"),
                    "complex_iptm": conf.get("protein_iptm") or conf.get("iptm"),
                    "prodigy_dg": prodigy_dg(model, antigen_chains),
                    "bsa_total": areas.get("bsa_total"),
                    "bsa_per_subunit_json": json.dumps(
                        areas.get("bsa_per_subunit") or {}),
                    "n_subunits_contacted": iface.get("n_subunits_contacted"),
                    "edge_fraction": iface.get("edge_fraction"),
                    "edge_artefact": iface.get("edge_artefact"),
                    "interface_residues_json": json.dumps(
                        iface.get("per_chain") or {}),
                    "nanobody_residues_json": json.dumps(
                        iface.get("nanobody_residues") or []),
                    "epitope_key": iface.get("epitope_key") or [],
                    "path": str(model),
                    "status": "ok",
                }
                meta = seqs.get(cid)
                if meta is None and cid in controls:
                    meta = _control_regions(controls[cid])
                record["is_control"] = 1 if cid in controls else 0
                if meta:
                    spans = cdr_spans(meta["sequence"], meta)
                    comp = paratope_composition(
                        iface.get("nanobody_residues") or [], spans)
                    record["cdr_interface_fraction"] = comp.get(
                        "cdr_interface_fraction")
                    record["cdr_interface_detail_json"] = json.dumps(
                        comp.get("cdr_interface_detail") or {})
                records.append(record)
        if not records:
            log(f"interface: no co-folds for {target} yet")
            continue

        # Epitope bins are per target, over the best sample of each candidate.
        best: dict[str, dict] = {}
        for rec in records:
            current = best.get(rec["candidate_id"])
            if current is None or (rec.get("iptm") or 0) > (current.get("iptm") or 0):
                best[rec["candidate_id"]] = rec
        bins = bin_epitopes(list(best.values()))
        by_candidate: dict[str, list[dict]] = defaultdict(list)
        for rec in records:
            by_candidate[rec["candidate_id"]].append(rec)
        consensus = {cid: pose_consensus(samples)
                     for cid, samples in by_candidate.items()}
        for rec in records:
            rec["epitope_bin"] = bins.get(rec["candidate_id"])
            # Per candidate, repeated on each of its rows, like epitope_bin.
            rec.update(consensus.get(rec["candidate_id"]) or {})

        published = 0
        for cid, rec in best.items():
            # Controls are scored, not shown. They have no candidate page, so
            # publishing their complexes would put scrambled-sequence models
            # on a public site with nothing to explain what they are.
            if rec.get("is_control"):
                continue
            if rec.get("path") and publish_complex(cid, target, Path(rec["path"])):
                published += 1
        if published:
            log(f"interface: {target} published {published} co-folded complexes "
                f"for the drawer, the highest-ipTM sample of each")

        columns = ["candidate_id", "target", "seed", "iptm", "ptm", "plddt",
                   "pae", "complex_iptm", "prodigy_dg", "bsa_total",
                   "bsa_per_subunit_json", "n_subunits_contacted",
                   "edge_fraction", "edge_artefact",
                   "epitope_bin", "interface_residues_json",
                   "nanobody_residues_json", "cdr_interface_fraction",
                   "cdr_interface_detail_json", "pose_consensus",
                   "iptm_spread", "nanobody_plddt", "antigen_plddt",
                   "interface_plddt", "is_control", "path", "status"]
        conn.executemany(
            f"INSERT OR REPLACE INTO cofolds ({', '.join(columns)}) VALUES "
            f"({', '.join(':' + c for c in columns)})",
            [{c: rec.get(c) for c in columns} for rec in records])
        conn.commit()
        multi = sum(1 for r in best.values()
                    if (r.get("n_subunits_contacted") or 0) > 1)
        edge = sum(1 for r in best.values() if r.get("edge_artefact"))
        if edge:
            log(f"interface: {target} {edge} of {len(best)} candidates place more "
                f"than {EDGE_ARTEFACT_FRACTION:.0%} of their interface on the "
                f"construct's cut edges. Those poses are artefacts of the "
                f"truncation and cannot be epitopes on the intact protein")
        log(f"interface: {target} {len(records)} samples across {len(best)} "
            f"candidates, {len(set(bins.values()))} epitope bins"
            + (f", {multi} of {len(best)} contact more than one subunit"
               if stoich > 1 else ""))
        manifest_write("interface", {"id": target, "status": "ok",
                                     "samples": len(records),
                                     "candidates": len(best),
                                     "bins": len(set(bins.values())),
                                     "multi_subunit": multi,
                                     "edge_artefacts": edge})
    conn.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
