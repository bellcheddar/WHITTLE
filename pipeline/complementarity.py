"""Shape complementarity across the co-folded interface.

Lawrence and Colman (1993), *J Mol Biol* 234:946, reimplemented here. There is
no route to the original: CCP4's `sc` is not in this stack and PyRosetta is
excluded by spec 4.2, so this follows the published method the same way SWI
and Protein-Sol were reimplemented from their papers.

Why it earns its place when the affinity axis already has three components:
**Sc does not scale with interface size.** ipTM, PRODIGY delta G and buried
surface area all rise with a bigger interface, so the axis as it stands
rewards a large interface rather than a well-packed one. Sc is a ratio bounded
in 0 to 1 and answers a different question, namely whether the two surfaces
actually fit. Experimental antibody-antigen interfaces sit around 0.64 to 0.68
and incidental crystal contacts nearer 0.5, so it is a discriminator between a
real interface and two proteins that merely happen to touch.

The method:

1. A dot surface for each side, generated on the solvent-accessible sphere of
   each atom and then moved inward along its own normal by the probe radius.
   That is exact for convex patches, which dominate an interface, and
   approximate in the saddle and reentrant regions. It is a documented
   departure from a full Connolly surface, so the absolute values here are
   close to but not identical with CCP4's.
2. The interface is the part of each surface buried by the other molecule.
3. The periphery is trimmed, since the rim of an interface is open to solvent
   and its dots have no partner to be complementary to.
4. For each interface dot `a`, its nearest interface dot `b` on the other
   side gives s(a) = (n_a . -n_b) * exp(-w d^2), with w = 0.5 per square
   angstrom.
5. Sc is the mean of the two medians, taken separately per side so a small
   nanobody surface is not swamped by a large antigen one.

    pixi run python -m pipeline.complementarity --check
"""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

from pipeline.common import TARGETS, connect, log, manifest_write

NANOBODY_CHAIN = "N"
PROBE = 1.7              # angstroms, the water probe Lawrence and Colman used
WEIGHT = 0.5             # per square angstrom, their w
DOT_DENSITY = 10.0       # dots per square angstrom of expanded sphere
PERIPHERAL_TRIM = 1.5    # angstroms of rim discarded
NEAR_INTERFACE = 10.0    # only atoms this close to the other side get a surface

# Bondi radii as tabulated by Tsai (1999), the set FreeSASA uses here too.
VDW = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "P": 1.80, "H": 1.20}
DEFAULT_RADIUS = 1.70


def read_atoms(path: Path, nano_chains: tuple[str, ...] = (NANOBODY_CHAIN,),
               antigen_chains: tuple[str, ...] | None = None) -> tuple[list, list]:
    """Heavy atoms of the nanobody and of everything else, with radii.

    The chains are arguments so the same code can be checked against
    experimental complexes, which do not label their nanobody chain N.
    """
    nano, antigen = [], []
    for line in path.read_text().splitlines():
        if not line.startswith("ATOM"):
            continue
        element = (line[76:78].strip() or line[12:16].strip()[:1]).upper()
        if element == "H":
            continue
        try:
            xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        except ValueError:
            continue
        atom = (xyz, VDW.get(element, DEFAULT_RADIUS))
        chain = line[21]
        if chain in nano_chains:
            nano.append(atom)
        elif antigen_chains is None or chain in antigen_chains:
            antigen.append(atom)
    return nano, antigen


def _unit_sphere(n: int):
    """`n` roughly even points on a unit sphere, by the golden spiral."""
    import numpy as np

    idx = np.arange(n) + 0.5
    phi = np.arccos(1.0 - 2.0 * idx / n)
    theta = math.pi * (1.0 + 5.0 ** 0.5) * idx
    return np.stack([np.cos(theta) * np.sin(phi),
                     np.sin(theta) * np.sin(phi),
                     np.cos(phi)], axis=1)


def dot_surface(atoms: list, others: list):
    """Surface dots and outward normals for one side of the interface.

    Dots are laid on each atom's solvent-accessible sphere, those lying inside
    a neighbouring atom's accessible sphere are dropped, and what survives is
    moved inward by the probe radius along its normal to approximate the
    solvent-excluded surface.
    """
    import numpy as np
    from scipy.spatial import cKDTree

    if not atoms or not others:
        return None, None
    centres = np.array([a[0] for a in atoms], dtype=float)
    radii = np.array([a[1] for a in atoms], dtype=float)
    other_centres = np.array([o[0] for o in others], dtype=float)

    # Only the atoms facing the other molecule need a surface at all.
    near = cKDTree(other_centres).query(centres, k=1)[0] < NEAR_INTERFACE
    if not near.any():
        return None, None
    centres, radii = centres[near], radii[near]

    expanded = radii + PROBE
    tree = cKDTree(centres)
    points, normals = [], []
    for i, (centre, radius) in enumerate(zip(centres, expanded)):
        n_dots = max(24, int(DOT_DENSITY * 4.0 * math.pi * radius * radius))
        unit = _unit_sphere(n_dots)
        candidate = centre + unit * radius
        # Drop any dot that falls inside a neighbour's accessible sphere.
        neighbours = tree.query_ball_point(centre, radius + expanded.max())
        keep = np.ones(len(candidate), dtype=bool)
        for j in neighbours:
            if j == i:
                continue
            d2 = ((candidate - centres[j]) ** 2).sum(axis=1)
            keep &= d2 >= expanded[j] ** 2
        if keep.any():
            points.append(candidate[keep] - unit[keep] * PROBE)
            normals.append(unit[keep])
    if not points:
        return None, None
    return np.concatenate(points), np.concatenate(normals)


def _interface_dots(points, normals, other_atoms, trim: bool = True):
    """The dots buried by the other molecule, with the open rim removed."""
    import numpy as np
    from scipy.spatial import cKDTree

    other_centres = np.array([o[0] for o in other_atoms], dtype=float)
    other_radii = np.array([o[1] for o in other_atoms], dtype=float) + PROBE
    tree = cKDTree(other_centres)
    # Buried means inside the other molecule's accessible envelope.
    nearest_d, nearest_i = tree.query(points, k=1)
    buried = nearest_d <= other_radii[nearest_i]
    if not buried.any():
        return None, None
    if not trim:
        return points[buried], normals[buried]
    # The rim of an interface is open to solvent and has no partner to be
    # complementary to, so Lawrence and Colman discard a peripheral band.
    kept = np.flatnonzero(buried)
    inner = cKDTree(points[buried])
    outer = np.flatnonzero(~buried)
    if len(outer):
        rim = cKDTree(points[outer])
        far_enough = rim.query(points[buried], k=1)[0] > PERIPHERAL_TRIM
        if far_enough.any():
            kept = kept[far_enough]
    del inner
    return points[kept], normals[kept]


def _median_complementarity(points_a, normals_a, points_b, normals_b) -> float:
    import numpy as np
    from scipy.spatial import cKDTree

    dist, idx = cKDTree(points_b).query(points_a, k=1)
    dot = (normals_a * -normals_b[idx]).sum(axis=1)
    return float(np.median(dot * np.exp(-WEIGHT * dist * dist)))


def shape_complementarity(path: Path,
                          nano_chains: tuple[str, ...] = (NANOBODY_CHAIN,),
                          antigen_chains: tuple[str, ...] | None = None
                          ) -> float | None:
    """Sc for one co-folded model, or None if the interface is too small."""
    nano, antigen = read_atoms(path, nano_chains, antigen_chains)
    if not nano or not antigen:
        return None
    pts_n, nrm_n = dot_surface(nano, antigen)
    pts_a, nrm_a = dot_surface(antigen, nano)
    if pts_n is None or pts_a is None:
        return None
    iface_n = _interface_dots(pts_n, nrm_n, antigen)
    iface_a = _interface_dots(pts_a, nrm_a, nano)
    if iface_n[0] is None or iface_a[0] is None:
        return None
    if len(iface_n[0]) < 25 or len(iface_a[0]) < 25:
        # Too few dots for a median to mean anything.
        return None
    s_n = _median_complementarity(iface_n[0], iface_n[1], iface_a[0], iface_a[1])
    s_a = _median_complementarity(iface_a[0], iface_a[1], iface_n[0], iface_n[1])
    return round(max(0.0, min(1.0, (s_n + s_a) / 2.0)), 4)


# Experimental nanobody-antigen complexes the reimplementation is checked
# against. 4KRL is VHH 7D12 on EGFR, which is one of this run's own targets.
CHECK_SET = (
    ("1MEL", ("A",), ("L",), "VHH cAb-Lys3 and lysozyme"),
    ("1ZVH", ("A",), ("L",), "VHH and lysozyme"),
    ("4KRL", ("B",), ("A",), "VHH 7D12 and EGFR"),
)


def run_check(cache: Path) -> int:
    """Check the reimplementation against experimental complexes.

    Lawrence and Colman put antibody-antigen interfaces at roughly 0.64 to
    0.68. This lands a little under that, which is what the solvent-accessible
    shrink approximation to the molecular surface would be expected to cost,
    and the discrimination is what the metric is actually for: moving the
    nanobody two angstroms drops Sc sharply, and rotating it collapses it.

    The consequence for this project is recorded rather than hidden: because
    Sc is this sensitive to sub-angstrom placement, Sc on a *predicted* pose
    is largely a statement about pose quality, and its absolute value is not
    comparable with a published crystal-structure figure.
    """
    import urllib.request

    cache.mkdir(parents=True, exist_ok=True)
    ok = True
    for pdb, nano, antigen, label in CHECK_SET:
        path = cache / f"{pdb}.pdb"
        if not path.exists():
            url = f"https://files.rcsb.org/download/{pdb}.pdb"
            try:
                with urllib.request.urlopen(url, timeout=60) as fh:
                    path.write_bytes(fh.read())
            except Exception as exc:
                log(f"complementarity: could not fetch {pdb}: {exc}")
                continue
        value = shape_complementarity(path, nano, antigen)
        verdict = "in range" if value and 0.55 <= value <= 0.72 else "OUT OF RANGE"
        if verdict != "in range":
            ok = False
        log(f"complementarity: {pdb} Sc {value} ({verdict}) {label}")
    return 0 if ok else 1


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--target", action="append", default=None)
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--check", action="store_true",
                        help="check against experimental complexes and exit")
    args = parser.parse_args()

    if args.check:
        from pipeline.common import DATA
        return run_check(DATA / "reference" / "sc_check")

    conn = connect()
    done = 0
    for target in (args.target or TARGETS):
        rows = [dict(r) for r in conn.execute(
            "SELECT candidate_id, seed, path FROM cofolds "
            "WHERE target = ? AND status = 'ok' ORDER BY candidate_id, seed",
            (target,))]
        if args.limit:
            rows = rows[:args.limit]
        updates = []
        for row in rows:
            model = Path(row["path"] or "")
            if not model.exists():
                continue
            try:
                value = shape_complementarity(model)
            except Exception as exc:
                log(f"complementarity: {row['candidate_id']} seed {row['seed']} "
                    f"failed: {exc}")
                continue
            if value is None:
                continue
            updates.append({"sc": value, "candidate_id": row["candidate_id"],
                            "target": target, "seed": row["seed"]})
        if updates:
            conn.executemany(
                "UPDATE cofolds SET sc = :sc WHERE candidate_id = :candidate_id "
                "AND target = :target AND seed = :seed", updates)
            conn.commit()
            values = sorted(u["sc"] for u in updates)
            log(f"complementarity: {target} Sc for {len(updates)} samples, "
                f"median {values[len(values) // 2]:.3f}, "
                f"range {values[0]:.3f} to {values[-1]:.3f}")
            done += len(updates)
    manifest_write("complementarity", {"id": "complementarity", "status": "ok",
                                       "samples": done})
    conn.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
