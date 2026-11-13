"""Surface hydrophobic and charged patches from the predicted model. Spec 6.2 item 3.

FreeSASA gives per-residue solvent accessibility; the patch definition on top of
it is stated here in full, because an unstated patch definition is an
uninterpretable number:

1. A residue is **exposed** when its relative accessibility exceeds
   `EXPOSURE_THRESHOLD` (0.25), against the Tien et al. (2013) theoretical
   maximum accessible areas.
2. Exposed residues are grouped by class: hydrophobic (A, V, L, I, M, F, W, Y,
   C), positive (K, R, and H, which carries partial charge at pH 7.4) and
   negative (D, E).
3. The **patch area** for a class is continuous, in the manner of Spatial
   Aggregation Propensity (Chennamsetty et al. (2009) PNAS 106:11937): for each
   exposed residue of that class, sum the exposed area of every same-class
   residue with a side-chain heavy atom within `NEIGHBOURHOOD_RADIUS`
   (7.5 angstroms, the radius the antibody developability literature uses), and
   take the largest such neighbourhood. That is the biggest contiguous surface
   of that character on the molecule.
4. The total exposed area of the class is reported beside it, and a
   single-linkage cluster count at the same radius is kept as a shape
   descriptor.

A hard "cluster of at least three" rule was tried first and rejected: on a
125-residue VHH it returns zero for almost every model, which is a threshold
artefact rather than a finding (see DECISIONS.md).

This is a geometric measurement of the model, not a trained predictor. The
trained developability axis is TNP; these columns exist so the drawer can say
*where* on the surface the hydrophobicity is, and so the expression axis has a
structural term that is not simply a function of composition.
"""
from __future__ import annotations

import gzip
import math
from collections import defaultdict
from pathlib import Path

EXPOSURE_THRESHOLD = 0.25
NEIGHBOURHOOD_RADIUS = 7.5
MIN_CLUSTER_RESIDUES = 3

HYDROPHOBIC = set("AVLIMFWYC")
POSITIVE = set("KRH")
NEGATIVE = set("DE")

# Tien, Meyer, Spielman, Wilke and Marcotte (2013) PLoS ONE 8:e80635,
# theoretical maximum accessible surface area per residue, square angstroms.
MAX_ASA = {
    "A": 129.0, "R": 274.0, "N": 195.0, "D": 193.0, "C": 167.0, "E": 223.0,
    "Q": 225.0, "G": 104.0, "H": 224.0, "I": 197.0, "L": 201.0, "K": 236.0,
    "M": 224.0, "F": 240.0, "P": 159.0, "S": 155.0, "T": 172.0, "W": 285.0,
    "Y": 263.0, "V": 174.0,
}

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}

BACKBONE = {"N", "CA", "C", "O", "OXT"}


def read_model(path: Path) -> tuple[str, list[dict]]:
    """(pdb text, residues). Residues carry their side-chain heavy atoms."""
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        text = fh.read()
    residues: dict[int, dict] = {}
    for line in text.splitlines():
        if not line.startswith("ATOM"):
            continue
        name = line[12:16].strip()
        element = line[76:78].strip() or name[:1]
        if element == "H":
            continue
        resnum = int(line[22:26])
        resname = line[17:20].strip()
        rec = residues.setdefault(resnum, {
            "num": resnum, "aa": THREE_TO_ONE.get(resname, "X"), "atoms": [],
            "sidechain": [],
        })
        xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        rec["atoms"].append(xyz)
        # Glycine has no side chain; use CA so it can still take part in a
        # contact graph rather than silently vanishing from it.
        if name not in BACKBONE or (rec["aa"] == "G" and name == "CA"):
            rec["sidechain"].append(xyz)
    return text, [residues[k] for k in sorted(residues)]


def residue_sasa(path: Path) -> dict[int, float]:
    """Per-residue total SASA, square angstroms, from FreeSASA."""
    import freesasa

    freesasa.setVerbosity(freesasa.silent)
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        text = fh.read()
    # FreeSASA reads a file; write the decompressed model beside it in a temp
    # file rather than keeping a second gzip-aware parser in step with this one.
    import tempfile

    with tempfile.NamedTemporaryFile("w", suffix=".pdb", delete=False) as tmp:
        tmp.write(text)
        tmp_path = tmp.name
    try:
        structure = freesasa.Structure(tmp_path)
        result = freesasa.calc(structure)
        areas = result.residueAreas()
    finally:
        Path(tmp_path).unlink(missing_ok=True)
    out: dict[int, float] = {}
    for _chain, residues in areas.items():
        for number, area in residues.items():
            try:
                out[int(number)] = float(area.total)
            except (TypeError, ValueError):
                continue
    return out


def _clusters(members: list[dict]) -> list[list[dict]]:
    """Single-linkage clustering on side-chain heavy-atom contact."""
    n = len(members)
    parent = list(range(n))

    def find(i: int) -> int:
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    def union(i: int, j: int) -> None:
        ri, rj = find(i), find(j)
        if ri != rj:
            parent[rj] = ri

    cutoff2 = NEIGHBOURHOOD_RADIUS ** 2
    for i in range(n):
        for j in range(i + 1, n):
            for ax, ay, az in members[i]["sidechain"]:
                hit = False
                for bx, by, bz in members[j]["sidechain"]:
                    if (ax - bx) ** 2 + (ay - by) ** 2 + (az - bz) ** 2 <= cutoff2:
                        union(i, j)
                        hit = True
                        break
                if hit:
                    break
    groups: dict[int, list[dict]] = defaultdict(list)
    for i in range(n):
        groups[find(i)].append(members[i])
    return list(groups.values())


def patches(path: Path, cdr_positions: set[int] | None = None) -> dict:
    """Patch areas for one model. Returns NaN-free floats, or an empty dict."""
    try:
        _text, residues = read_model(path)
        sasa = residue_sasa(path)
    except Exception:
        return {}
    if not residues or not sasa:
        return {}

    exposed: dict[str, list[dict]] = {"hydrophobic": [], "positive": [], "negative": []}
    total_sasa = 0.0
    for res in residues:
        area = sasa.get(res["num"])
        if area is None:
            continue
        total_sasa += area
        maximum = MAX_ASA.get(res["aa"])
        if not maximum:
            continue
        relative = area / maximum
        res["sasa"] = area
        res["rsa"] = relative
        if relative < EXPOSURE_THRESHOLD or not res["sidechain"]:
            continue
        if res["aa"] in HYDROPHOBIC:
            exposed["hydrophobic"].append(res)
        elif res["aa"] in POSITIVE:
            exposed["positive"].append(res)
        elif res["aa"] in NEGATIVE:
            exposed["negative"].append(res)

    out: dict[str, float] = {"sasa_total": round(total_sasa, 1)}
    radius2 = NEIGHBOURHOOD_RADIUS ** 2
    for label, members in exposed.items():
        total = sum(r["sasa"] for r in members)
        largest = 0.0
        for i, centre in enumerate(members):
            neighbourhood = centre["sasa"]
            for j, other in enumerate(members):
                if i == j:
                    continue
                close = False
                for ax, ay, az in centre["sidechain"]:
                    for bx, by, bz in other["sidechain"]:
                        if (ax - bx) ** 2 + (ay - by) ** 2 + (az - bz) ** 2 <= radius2:
                            close = True
                            break
                    if close:
                        break
                if close:
                    neighbourhood += other["sasa"]
            largest = max(largest, neighbourhood)
        n_clusters = sum(1 for c in _clusters(members) if len(c) >= MIN_CLUSTER_RESIDUES)
        out[f"{label}_patch_area"] = round(largest, 1)
        out[f"{label}_exposed_area"] = round(total, 1)
        out[f"{label}_patch_count"] = n_clusters
        out[f"{label}_exposed_residues"] = len(members)

    out["charged_patch_area"] = round(
        max(out.get("positive_patch_area", 0.0), out.get("negative_patch_area", 0.0)), 1)
    out["charged_exposed_area"] = round(
        out.get("positive_exposed_area", 0.0) + out.get("negative_exposed_area", 0.0), 1)

    if cdr_positions:
        cdr_exposed = sum(r.get("sasa", 0.0) for r in residues
                          if r["num"] in cdr_positions
                          and r.get("rsa", 0) >= EXPOSURE_THRESHOLD
                          and r["aa"] in HYDROPHOBIC)
        out["cdr_hydrophobic_exposed"] = round(cdr_exposed, 1)
    return out


def definition() -> dict:
    """The patch definition, for the About tab. Stated, not implied."""
    return {
        "exposure_threshold": EXPOSURE_THRESHOLD,
        "neighbourhood_radius_angstroms": NEIGHBOURHOOD_RADIUS,
        "minimum_cluster_residues": MIN_CLUSTER_RESIDUES,
        "hydrophobic": "".join(sorted(HYDROPHOBIC)),
        "positive": "".join(sorted(POSITIVE)),
        "negative": "".join(sorted(NEGATIVE)),
        "max_asa_reference": "Tien et al. (2013) PLoS ONE 8:e80635",
        "sasa_engine": "FreeSASA, Lee-Richards",
        "patch_definition": "Chennamsetty et al. (2009) PNAS 106:11937, in spirit",
        "note": ("Relative accessibility is measured against the theoretical "
                 "maximum per residue. The patch area for a class is the largest "
                 "neighbourhood sum: for each exposed residue of that class, the "
                 "exposed area of every same-class residue within 7.5 angstroms, "
                 "maximised over residues. Total exposed area and a cluster count "
                 "are reported beside it."),
    }


if __name__ == "__main__":
    import json
    import sys

    for arg in sys.argv[1:]:
        print(arg, json.dumps(patches(Path(arg)), indent=2))
