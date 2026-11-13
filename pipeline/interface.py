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
import sys
from collections import defaultdict
from pathlib import Path

from pipeline.common import (RUNS, TARGETS, connect, log, manifest_read,
                             manifest_write, targets_config)

CONTACT_CUTOFF = 5.0          # heavy-atom distance, angstroms
EPITOPE_JACCARD = 0.5         # two candidates share a bin above this overlap
NANOBODY_CHAIN = "N"


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


def confidence(run_dir: Path) -> dict:
    """boltz's own confidence numbers, averaged over the diffusion samples."""
    files = sorted(run_dir.rglob("confidence_*.json"))
    if not files:
        return {}
    values: dict[str, list[float]] = defaultdict(list)
    for path in files:
        try:
            payload = json.loads(path.read_text())
        except json.JSONDecodeError:
            continue
        for key in ("confidence_score", "ptm", "iptm", "complex_plddt",
                    "complex_iplddt", "complex_pde", "complex_ipde",
                    "ligand_iptm", "protein_iptm"):
            if isinstance(payload.get(key), (int, float)):
                values[key].append(float(payload[key]))
    out = {k: round(sum(v) / len(v), 4) for k, v in values.items() if v}
    out["n_samples"] = len(files)
    return out


def structures(run_dir: Path) -> list[Path]:
    return sorted(p for p in run_dir.rglob("*.pdb")
                  if "pre_affinity" not in p.name)


# ------------------------------------------------------------------- metrics

def interface_residues(chains: dict[str, list[tuple]], antigen_chains: list[str]
                       ) -> dict:
    """Antigen residues contacting the nanobody, and which chains they are on."""
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
    return {
        "per_chain": {c: sorted(r) for c, r in contacts.items() if r},
        "n_antigen_residues": sum(len(r) for r in contacts.values()),
        "n_nanobody_residues": len(nano_contacts),
        "n_subunits_contacted": len(touched),
        "chains_contacted": touched,
        # The identity of an epitope, for binning. Residue numbers only, with
        # the chain dropped: in a homotrimer the same face on chain A and on
        # chain B is the same epitope.
        "epitope_key": sorted({r for residues in contacts.values()
                               for r in residues}),
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
        records: list[dict] = []
        for cid, run in done_runs.items():
            if run.get("target") != target:
                continue
            run_dir = Path(run["path"])
            if not run_dir.exists():
                continue
            conf = confidence(run_dir)
            models = structures(run_dir)
            if not models:
                manifest_write("interface", {"id": cid, "target": target,
                                             "status": "failed:no structure"})
                continue
            for sample, model in enumerate(models):
                chains = read_complex(model)
                iface = interface_residues(chains, antigen_chains)
                areas = buried_area(model, antigen_chains)
                record = {
                    "candidate_id": cid, "target": target, "seed": sample,
                    "iptm": conf.get("iptm"), "ptm": conf.get("ptm"),
                    "plddt": conf.get("complex_plddt"),
                    "pae": conf.get("complex_pde"),
                    "complex_iptm": conf.get("protein_iptm") or conf.get("iptm"),
                    "prodigy_dg": prodigy_dg(model, antigen_chains),
                    "bsa_total": areas.get("bsa_total"),
                    "bsa_per_subunit_json": json.dumps(
                        areas.get("bsa_per_subunit") or {}),
                    "n_subunits_contacted": iface.get("n_subunits_contacted"),
                    "interface_residues_json": json.dumps(
                        iface.get("per_chain") or {}),
                    "epitope_key": iface.get("epitope_key") or [],
                    "path": str(model),
                    "status": "ok",
                }
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
        for rec in records:
            rec["epitope_bin"] = bins.get(rec["candidate_id"])

        columns = ["candidate_id", "target", "seed", "iptm", "ptm", "plddt",
                   "pae", "complex_iptm", "prodigy_dg", "bsa_total",
                   "bsa_per_subunit_json", "n_subunits_contacted",
                   "epitope_bin", "interface_residues_json", "path", "status"]
        conn.executemany(
            f"INSERT OR REPLACE INTO cofolds ({', '.join(columns)}) VALUES "
            f"({', '.join(':' + c for c in columns)})",
            [{c: rec.get(c) for c in columns} for rec in records])
        conn.commit()
        multi = sum(1 for r in best.values()
                    if (r.get("n_subunits_contacted") or 0) > 1)
        log(f"interface: {target} {len(records)} samples across {len(best)} "
            f"candidates, {len(set(bins.values()))} epitope bins"
            + (f", {multi} of {len(best)} contact more than one subunit"
               if stoich > 1 else ""))
        manifest_write("interface", {"id": target, "status": "ok",
                                     "samples": len(records),
                                     "candidates": len(best),
                                     "bins": len(set(bins.values())),
                                     "multi_subunit": multi})
    conn.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
