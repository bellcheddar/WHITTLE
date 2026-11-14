#!/usr/bin/env python3
"""Where does the pLDDT gap between a candidate and its own shuffled twin live?

    ./tools/plddt_locus.py

This exists because `plddt` won the affinity calibration and is not a component
of the affinity axis. Pooled over the 50 paired controls, a candidate's complex
carries a higher `complex_plddt` than the complex of the control shuffled from
it on 41 of 50 pairs (p = 5.6e-06), separating real from shuffled better than
every metric the axis actually weights. That is either a missing component or a
measurement of something other than binding, and the two look identical in the
calibration output, so the question has to be settled somewhere else.

`complex_plddt` is an average over the whole complex, which for these targets is
mostly antigen: the nanobody is about 125 residues against 130 for the GLP1R ECD
and roughly 450 for the TNFa trimer. This splits the average by chain, and the
split is interpretable because of what a CDR shuffle holds fixed. Within a pair
the antigen is the identical sequence folded from the identical cached MSA, and
the nanobody differs only in the order of its CDR residues: same length, same
composition, same framework, same seeds, same five diffusion samples. So the
antigen chain is a negative control on the measurement itself. If it moves, the
pipeline is doing something systematic to one arm of each pair and no reading of
the nanobody chain can be trusted.

Boltz writes no per-chain pLDDT, so the chain means are read from the B-factor
column of the predicted PDBs, which is where pLDDT is stored. `complex_iplddt`
is taken from the confidence JSON beside each sample, and is reported because
an interface-local confidence is the one form of pLDDT that could honestly
belong on an affinity axis.

The result, 2026-11-14, is recorded in DECISIONS 49: the gap is in the nanobody
chain (48 of 50 pairs, p = 2.3e-12, mean 5.1 pLDDT), the antigen chain is at
chance (31 of 50, p = 0.12), so `plddt` is reporting how predictable a VHH's own
loops are rather than anything about the interface. Nothing here is wired into
scoring; it is a diagnostic, and the one conclusion it licences is a negative.
"""
import json
import re
import sqlite3
import statistics
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tools"))

from pipeline.interface import structures  # noqa: E402
from vacuous_tests import sign_test  # noqa: E402

DB = ROOT / "data/whittle.sqlite"
RUNS = ROOT / "data/runs"
# Higher is better for every one of these, so "parent higher" is "parent wins".
FIELDS = ["complex_plddt", "nanobody_plddt", "antigen_plddt", "interface_plddt"]
NANOBODY_CHAIN = "N"


def pairs(conn) -> list[tuple[str, str, str]]:
    """Each scrambled control with the candidate it was shuffled from.

    The parent is recovered from the control's own note rather than a join,
    because that note is what `pipeline/controls.py` wrote when it made
    the shuffle and is the only record of which candidate it came from.
    """
    folded = {r["candidate_id"] for r in
              conn.execute("SELECT DISTINCT candidate_id FROM cofolds")}
    out = []
    for row in conn.execute("SELECT id, target, notes FROM controls "
                            "WHERE kind = 'scrambled' AND sequence IS NOT NULL"):
        found = re.search(r"derived from ([A-Za-z0-9]+-\d+)", row["notes"] or "")
        if not found:
            continue
        parent = found.group(1)
        if row["id"] in folded and parent in folded:
            out.append((row["id"], parent, row["target"]))
    return out


def run_dir(target: str, cid: str) -> Path | None:
    for cand in (RUNS / target / cid, RUNS / "controls" / target / cid):
        if cand.is_dir():
            return cand
    hits = [p for p in RUNS.rglob(f"*{cid}*") if p.is_dir()]
    return hits[0] if hits else None


def chain_means(pdb: Path) -> tuple[float | None, float | None]:
    """Mean CA B-factor for the nanobody chain and for the antigen chains."""
    nb: list[float] = []
    ag: list[float] = []
    with pdb.open() as fh:
        for line in fh:
            if not line.startswith("ATOM") or line[12:16].strip() != "CA":
                continue
            (nb if line[21] == NANOBODY_CHAIN else ag).append(float(line[60:66]))
    return (statistics.fmean(nb) if nb else None,
            statistics.fmean(ag) if ag else None)


def measure(target: str, cid: str) -> dict[str, float] | None:
    """One complex, averaged over its diffusion samples."""
    directory = run_dir(target, cid)
    if directory is None:
        return None
    nb, ag, interface, complexes = [], [], [], []
    for pdb in structures(directory):
        chain_nb, chain_ag = chain_means(pdb)
        if chain_nb is not None:
            nb.append(chain_nb)
        if chain_ag is not None:
            ag.append(chain_ag)
        conf = pdb.parent / f"confidence_{pdb.stem}.json"
        if conf.is_file():
            body = json.loads(conf.read_text())
            for key, sink in (("complex_iplddt", interface),
                              ("complex_plddt", complexes)):
                if body.get(key) is not None:
                    sink.append(float(body[key]))
    if not nb:
        return None
    mean = statistics.fmean
    return {"nanobody_plddt": mean(nb),
            "antigen_plddt": mean(ag) if ag else None,
            "interface_plddt": mean(interface) if interface else None,
            "complex_plddt": mean(complexes) if complexes else None}


def main() -> int:
    conn = sqlite3.connect(f"file:{DB}?mode=ro", uri=True)
    conn.row_factory = sqlite3.Row
    todo = pairs(conn)

    wins = {f: [0, 0] for f in FIELDS}
    deltas: dict[str, list[float]] = {f: [] for f in FIELDS}
    by_target: dict[str, dict[str, list[int]]] = {}
    unmeasurable = 0

    for cid, parent, target in todo:
        control, real = measure(target, cid), measure(target, parent)
        if control is None or real is None:
            unmeasurable += 1
            continue
        for field in FIELDS:
            if control.get(field) is None or real.get(field) is None:
                continue
            won = real[field] > control[field]
            wins[field][1] += 1
            wins[field][0] += won
            deltas[field].append(real[field] - control[field])
            tally = by_target.setdefault(target, {}).setdefault(field, [0, 0])
            tally[1] += 1
            tally[0] += won

    print(f"pairs: {len(todo)}, unmeasurable: {unmeasurable}")
    print("\nparent beats the control shuffled from it:\n")
    print(f"{'pLDDT over':18} {'parent higher':>14} {'frac':>6} {'p':>9} {'mean delta':>12}")
    for field in FIELDS:
        won, n = wins[field]
        if not n:
            print(f"{field:18} {'no data':>14}")
            continue
        delta = statistics.fmean(deltas[field]) if deltas[field] else float("nan")
        print(f"{field:18} {won:>8}/{n:<5} {won / n:6.2f} "
              f"{sign_test(won, n):9.2g} {delta:+12.4f}")

    print("\nper target, parent higher of decided pairs")
    print(f"{'target':8}" + "".join(
        f"{f.replace('_plddt', ''):>12}" for f in FIELDS))
    for target in sorted(by_target):
        row = f"{target:8}"
        for field in FIELDS:
            won, n = by_target[target].get(field, [0, 0])
            row += f"{f'{won}/{n}' if n else '-':>12}"
        print(row)

    nb_won, nb_n = wins["nanobody_plddt"]
    ag_won, ag_n = wins["antigen_plddt"]
    if nb_n and ag_n:
        print()
        ag_p = sign_test(ag_won, ag_n)
        if ag_p is not None and ag_p < 0.05:
            print("The antigen chain is also separating. Within a pair it is the "
                  "same sequence and the same MSA, so something systematic is "
                  "acting on one arm of each pair and the nanobody chain figure "
                  "below cannot be read as a property of the shuffle.")
        else:
            print(f"The antigen chain is at chance ({ag_won}/{ag_n}, p = {ag_p:.2g}), "
                  "which is what an identical sequence folded from an identical "
                  "MSA should do, so the measurement carries no systematic bias "
                  "between the two arms of a pair.")
        print(f"The gap is in the nanobody chain ({nb_won}/{nb_n}, "
              f"mean {statistics.fmean(deltas['nanobody_plddt']):+.2f} pLDDT). A "
              "shuffled VHH is a less predictable protein, which is foldability "
              "and not binding, so pLDDT does not belong on the affinity axis.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
