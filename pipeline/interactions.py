"""Non-covalent interactions across the co-folded interface. Spec 6.4.

PLIP (Salentin 2015, Adasme 2021, GPLv2) profiles the contacts between the
nanobody and the antigen: hydrogen bonds, hydrophobic contacts, salt bridges,
pi-stacking and pi-cation. PLIP is built for protein-ligand complexes, so the
nanobody chain is declared a peptide chain (`config.PEPTIDES`), which is the
route its own authors provide for protein-protein and protein-peptide work.

Two things are produced, and they are deliberately different objects:

* **Counts and densities**, stored per sample in `cofolds`. The counts on their
  own are close to a restatement of interface size, which the affinity axis
  already has three of (ipTM, PRODIGY delta G, BSA). What is not a restatement
  is the **density** and the **polar fraction**: hydrogen bonds per 100 square
  angstroms separates a well-packed polar interface from a merely large one.
* **The interactions themselves**, written beside the published complex as
  `<id>.interactions.json`, so the drawer can list them and the viewer can
  draw them. Both endpoints carry coordinates, so the viewer needs no atom
  index mapping and the picture cannot disagree with the list.

Everything here is conditioned on a predicted pose. A detailed contact
inventory of a wrong pose is false precision, not information, so the caller
gates on `edge_artefact` and on pose confidence before any of it is ranked.

    pixi run python -m pipeline.interactions
"""
from __future__ import annotations

import argparse
import json
import sys
import warnings
from pathlib import Path

from pipeline.common import (STRUCTURES, TARGETS, connect, log, manifest_write,
                             targets_config)

NANOBODY_CHAIN = "N"

# What a published interaction record is keyed on. `antigen` and `nanobody`
# rather than PLIP's protein/ligand, because "ligand" for a 111-residue VHH
# would read as an error in the drawer.
POLAR_TYPES = ("hbond", "saltbridge", "waterbridge")


def _atom_name(atom) -> str | None:
    """The PDB atom name, not OpenBabel's atom type.

    PLIP's `atype`/`dtype` are OpenBabel perception types, so an interaction
    would read as "O3 to O3" in the drawer where a reader expects "OH to OG".
    """
    try:
        ob = atom.OBAtom
        residue = ob.GetResidue()
        if residue is None:
            return None
        name = residue.GetAtomID(ob).strip()
        return name or None
    except Exception:
        return None


def _xyz(value) -> list[float] | None:
    try:
        return [round(float(v), 3) for v in value]
    except Exception:
        return None


def _centre(group) -> list[float] | None:
    for attr in ("center", "centroid"):
        got = getattr(group, attr, None)
        if got is not None:
            out = _xyz(got)
            if out:
                return out
    atoms = getattr(group, "atoms", None)
    if atoms:
        pts = [a.coords for a in atoms]
        if pts:
            return [round(sum(p[i] for p in pts) / len(pts), 3) for i in range(3)]
    return None


def _record(kind: str, antigen: dict, nanobody: dict, dist: float | None,
            a_xyz, n_xyz, extra: dict | None = None) -> dict:
    rec = {
        "type": kind,
        "antigen": antigen,
        "nanobody": nanobody,
        "distance": round(float(dist), 2) if dist is not None else None,
        "a_xyz": a_xyz,
        "n_xyz": n_xyz,
    }
    if extra:
        rec.update(extra)
    return rec


def _side(restype, resnr, reschain, atom=None) -> dict:
    out = {"resname": str(restype), "resnum": int(resnr), "chain": str(reschain)}
    if atom is not None:
        out["atom"] = str(atom)
    return out


def _dedupe(records: list[dict]) -> list[dict]:
    """One physical contact, one record.

    Two hydroxyls can each donate to the other, so PLIP reports that pair
    twice, once per direction, with a different hydrogen placement and angle.
    Counting both would inflate the hydrogen bond count and the polar fraction
    on exactly the polar interfaces the density metric is meant to identify.
    The pair with the better donor angle is kept.
    """
    best: dict[tuple, dict] = {}
    for rec in records:
        key = (rec["type"],
               rec["antigen"]["chain"], rec["antigen"]["resnum"],
               rec["antigen"].get("atom"),
               rec["nanobody"]["chain"], rec["nanobody"]["resnum"],
               rec["nanobody"].get("atom"))
        current = best.get(key)
        if current is None or (rec.get("angle") or 0) > (current.get("angle") or 0):
            best[key] = rec
    return list(best.values())


def profile(model: Path) -> dict:
    """Every nanobody-antigen interaction in one co-folded model."""
    warnings.filterwarnings("ignore")
    from plip.basic import config

    # Declare the nanobody a peptide chain. Without this PLIP finds no ligand
    # in an all-protein complex and reports nothing at all.
    config.PEPTIDES = [NANOBODY_CHAIN]
    from plip.structure.preparation import PDBComplex

    complex_ = PDBComplex()
    complex_.load_pdb(str(model))
    complex_.analyze()
    if not complex_.interaction_sets:
        return {"records": [], "counts": {}}

    records: list[dict] = []
    for site in complex_.interaction_sets.values():
        for hb in list(getattr(site, "hbonds_pdon", [])) + \
                  list(getattr(site, "hbonds_ldon", [])):
            # protisdon says which side donates, which decides which atom
            # belongs to which chain.
            donor, acceptor = _xyz(hb.d.coords), _xyz(hb.a.coords)
            prot_is_donor = bool(getattr(hb, "protisdon", True))
            d_name, a_name = _atom_name(hb.d), _atom_name(hb.a)
            records.append(_record(
                "hbond",
                _side(hb.restype, hb.resnr, hb.reschain,
                      d_name if prot_is_donor else a_name),
                _side(hb.restype_l, hb.resnr_l, hb.reschain_l,
                      a_name if prot_is_donor else d_name),
                hb.distance_ad,
                donor if prot_is_donor else acceptor,
                acceptor if prot_is_donor else donor,
                {"angle": round(float(hb.angle), 1),
                 "sidechain": bool(getattr(hb, "sidechain", False)),
                 "protein_donates": prot_is_donor}))

        for hc in getattr(site, "hydrophobic_contacts", []):
            records.append(_record(
                "hydrophobic",
                _side(hc.restype, hc.resnr, hc.reschain,
                      _atom_name(hc.bsatom)),
                _side(hc.restype_l, hc.resnr_l, hc.reschain_l,
                      _atom_name(hc.ligatom)),
                hc.distance, _xyz(hc.bsatom.coords), _xyz(hc.ligatom.coords)))

        for sb in list(getattr(site, "saltbridge_lneg", [])) + \
                  list(getattr(site, "saltbridge_pneg", [])):
            pos, neg = getattr(sb, "positive", None), getattr(sb, "negative", None)
            # The protein side is whichever group is not on the nanobody.
            a_pt, n_pt = _centre(pos), _centre(neg)
            if str(getattr(sb, "reschain", "")) == NANOBODY_CHAIN:
                a_pt, n_pt = n_pt, a_pt
            records.append(_record(
                "saltbridge",
                _side(sb.restype, sb.resnr, sb.reschain),
                _side(sb.restype_l, sb.resnr_l, sb.reschain_l),
                getattr(sb, "distance", None), a_pt, n_pt))

        for ps in getattr(site, "pistacking", []):
            records.append(_record(
                "pistacking",
                _side(ps.restype, ps.resnr, ps.reschain),
                _side(ps.restype_l, ps.resnr_l, ps.reschain_l),
                getattr(ps, "distance", None),
                _centre(getattr(ps, "proteinring", None)),
                _centre(getattr(ps, "ligandring", None)),
                {"stack_type": str(getattr(ps, "type", ""))}))

        for pc in list(getattr(site, "pication_laro", [])) + \
                  list(getattr(site, "pication_paro", [])):
            records.append(_record(
                "pication",
                _side(pc.restype, pc.resnr, pc.reschain),
                _side(pc.restype_l, pc.resnr_l, pc.reschain_l),
                getattr(pc, "distance", None),
                _centre(getattr(pc, "ring", None)),
                _centre(getattr(pc, "charge", None))))

    records = _dedupe(records)
    counts: dict[str, int] = {}
    for rec in records:
        counts[rec["type"]] = counts.get(rec["type"], 0) + 1
    # Design intent for the two His-enriched libraries: did a CDR histidine
    # actually reach the interface and make a contact? A pH-switch library
    # whose histidines never touch the antigen has not done its job.
    counts["his_contacts"] = sum(
        1 for r in records if r["nanobody"]["resname"] == "HIS")
    counts["his_residues"] = len({
        r["nanobody"]["resnum"] for r in records
        if r["nanobody"]["resname"] == "HIS"})
    return {"records": records, "counts": counts}


def derive(counts: dict, bsa_total: float | None) -> dict:
    """Densities, which is the part that is not a restatement of size."""
    total = sum(counts.get(k, 0) for k in
                ("hbond", "hydrophobic", "saltbridge", "pistacking", "pication"))
    polar = sum(counts.get(k, 0) for k in POLAR_TYPES)
    out = {
        "n_hbond": counts.get("hbond", 0),
        "n_hydrophobic": counts.get("hydrophobic", 0),
        "n_saltbridge": counts.get("saltbridge", 0),
        "n_pistack": counts.get("pistacking", 0) + counts.get("pication", 0),
        "n_interactions": total,
        "n_his_contacts": counts.get("his_contacts", 0),
        "polar_fraction": round(polar / total, 3) if total else None,
        "hbond_per_100a2": None,
        "interactions_per_100a2": None,
    }
    if bsa_total and bsa_total > 0:
        out["hbond_per_100a2"] = round(
            counts.get("hbond", 0) / (bsa_total / 100.0), 3)
        out["interactions_per_100a2"] = round(total / (bsa_total / 100.0), 3)
    return out


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--target", action="append", default=None)
    parser.add_argument("--limit", type=int, default=None)
    args = parser.parse_args()

    conn = connect()
    targets = args.target or TARGETS
    cfg = targets_config()["targets"]
    total_done = 0

    for target in targets:
        rows = [dict(r) for r in conn.execute(
            "SELECT candidate_id, seed, path, bsa_total, iptm, "
            "COALESCE(is_control, 0) AS is_control FROM cofolds "
            "WHERE target = ? AND status = 'ok' ORDER BY candidate_id, seed",
            (target,))]
        if not rows:
            continue
        if args.limit:
            rows = rows[:args.limit]

        best: dict[str, tuple[float, dict, list]] = {}
        updates: list[dict] = []
        for row in rows:
            model = Path(row["path"] or "")
            if not model.exists():
                continue
            try:
                out = profile(model)
            except Exception as exc:
                log(f"interactions: {row['candidate_id']} seed {row['seed']} "
                    f"failed: {exc}")
                continue
            metrics = derive(out["counts"], row.get("bsa_total"))
            updates.append({**metrics, "candidate_id": row["candidate_id"],
                            "target": target, "seed": row["seed"]})
            # The published complex is the highest-ipTM sample, so the
            # interaction list has to come from that same sample or the
            # drawer would describe a structure it is not showing.
            iptm = row.get("iptm") or 0.0
            if row["candidate_id"] not in best or iptm > best[row["candidate_id"]][0]:
                best[row["candidate_id"]] = (iptm, metrics, out["records"])

        if updates:
            conn.executemany(
                "UPDATE cofolds SET n_hbond = :n_hbond, "
                "n_hydrophobic = :n_hydrophobic, n_saltbridge = :n_saltbridge, "
                "n_pistack = :n_pistack, n_interactions = :n_interactions, "
                "n_his_contacts = :n_his_contacts, "
                "polar_fraction = :polar_fraction, "
                "hbond_per_100a2 = :hbond_per_100a2, "
                "interactions_per_100a2 = :interactions_per_100a2 "
                "WHERE candidate_id = :candidate_id AND target = :target "
                "AND seed = :seed", updates)
            conn.commit()

        # Controls are profiled for the calibration but never published:
        # there is no page to show them on.
        control_ids = {r["candidate_id"] for r in rows if r["is_control"]}
        out_dir = STRUCTURES / target
        out_dir.mkdir(parents=True, exist_ok=True)
        for cid, (iptm, metrics, records) in best.items():
            if cid in control_ids:
                continue
            payload = {
                "candidate_id": cid, "target": target,
                "iptm": iptm, "counts": metrics, "interactions": records,
                "note": ("Interactions are computed on a predicted co-folded "
                         "pose, not an experimental structure."),
            }
            (out_dir / f"{cid}.interactions.json").write_text(
                json.dumps(payload, separators=(",", ":")))
        total_done += len(best)
        log(f"interactions: {target} profiled {len(updates)} samples across "
            f"{len(best)} candidates, "
            f"median {sorted(u['n_interactions'] for u in updates)[len(updates)//2]} "
            f"interactions per sample")

    manifest_write("interactions", {"id": "interactions", "status": "ok",
                                    "candidates": total_done})
    conn.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
