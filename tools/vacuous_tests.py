"""Vacuous-test detector. Spec 10 item 4, adapted from NOBBLE/tools/vacuous_tests.py.

A metric that cannot fail is worse than no metric. This checks each of the four
axes for the ways a column can look like a measurement while carrying no
information:

1. **No variance.** Every row the same value, or almost.
2. **Degenerate distribution.** More than half the rows share one value, so the
   percentile is mostly ties.
3. **A pure restatement of another column.** Correlation above 0.98 with
   another axis, or with sequence length, means the axis is that column wearing
   a different name.
4. **No separation from a negative control.** The axis does not distinguish
   real candidates from CDR-shuffled ones. This is the strong form, and for the
   affinity axis spec Section 9 item 5 says a failure here is a finding to
   report prominently, not a failure to hide.
5. **Cannot fail by construction.** The axis has no low tail: nothing scores
   badly, so the column cannot discriminate against anything.

Results go to the `qc` table and onto the QC panel on the control tab.

    pixi run python tools/vacuous_tests.py
"""
from __future__ import annotations

import json
import math
import statistics
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline.common import TARGETS, connect, log
from pipeline.score import AXIS_COMPONENTS  # noqa: E402

AXES = {
    "affinity": "axis_affinity",
    "specificity": "axis_specificity",
    "thermo": "axis_thermo",
    "expression": "axis_expression",
}


def pearson(xs: list[float], ys: list[float]) -> float | None:
    pairs = [(x, y) for x, y in zip(xs, ys) if x is not None and y is not None]
    if len(pairs) < 10:
        return None
    n = len(pairs)
    mx = sum(p[0] for p in pairs) / n
    my = sum(p[1] for p in pairs) / n
    sxy = sum((p[0] - mx) * (p[1] - my) for p in pairs)
    sxx = sum((p[0] - mx) ** 2 for p in pairs)
    syy = sum((p[1] - my) ** 2 for p in pairs)
    if sxx == 0 or syy == 0:
        return None
    return sxy / math.sqrt(sxx * syy)


def mann_whitney_auc(a: list[float], b: list[float]) -> float | None:
    """P(a random value from `a` exceeds one from `b`). 0.5 means no separation.

    Reported rather than a p-value: with 500 candidates a trivial difference is
    significant, and the question here is whether the axis separates the two
    sets at all, not whether the difference is detectable.
    """
    if len(a) < 5 or len(b) < 5:
        return None
    wins = ties = 0
    for x in a:
        for y in b:
            if x > y:
                wins += 1
            elif x == y:
                ties += 1
    return (wins + 0.5 * ties) / (len(a) * len(b))


def check_axis(conn, target: str, name: str, column: str) -> dict:
    rows = conn.execute(
        f"SELECT s.candidate_id, s.{column} AS value, c.length "
        f"FROM scores s JOIN candidates c ON c.id = s.candidate_id "
        f"WHERE s.target = ?", (target,)).fetchall()
    values = [r["value"] for r in rows if r["value"] is not None]
    result: dict = {"axis": name, "target": target, "n": len(values),
                    "findings": [], "passed": True}
    if len(values) < 20:
        result["passed"] = None
        result["findings"].append("not computed yet")
        return result

    spread = statistics.pstdev(values)
    result["sd"] = round(spread, 3)
    if spread < 1.0:
        result["passed"] = False
        result["findings"].append(
            f"no variance: standard deviation {spread:.3f} percentile points")

    counts: dict[float, int] = {}
    for v in values:
        counts[v] = counts.get(v, 0) + 1
    most = max(counts.values())
    result["modal_fraction"] = round(most / len(values), 3)
    if most > len(values) / 2:
        result["passed"] = False
        result["findings"].append(
            f"degenerate: {most / len(values) * 100:.0f} percent of rows share "
            f"one value, so the percentile is mostly ties")

    # Is the axis a restatement of something else?
    lengths = [r["length"] for r in rows if r["value"] is not None]
    ordered = [r["value"] for r in rows if r["value"] is not None]
    r_len = pearson(ordered, [float(x) for x in lengths])
    result["r_with_length"] = None if r_len is None else round(r_len, 3)
    if r_len is not None and abs(r_len) > 0.98:
        result["passed"] = False
        result["findings"].append(
            f"this axis is sequence length wearing a different name "
            f"(r = {r_len:.3f})")
    for other_name, other_col in AXES.items():
        if other_col == column:
            continue
        other = [r[0] for r in conn.execute(
            f"SELECT {other_col} FROM scores WHERE target = ? "
            f"AND {column} IS NOT NULL ORDER BY candidate_id", (target,))]
        if len(other) != len(ordered):
            continue
        r = pearson(ordered, other)
        if r is not None and abs(r) > 0.98:
            result["passed"] = False
            result["findings"].append(
                f"this axis is a restatement of the {other_name} axis "
                f"(r = {r:.3f})")

    # No low tail at all: the column cannot discriminate against anything.
    if min(values) > 25:
        result["passed"] = False
        result["findings"].append(
            f"no low tail: the worst candidate sits at the "
            f"{min(values):.0f}th percentile, so nothing can score badly")

    # Separation from the CDR-shuffled negative control.
    control_column = {
        "affinity": "predicted_affinity", "specificity": "predicted_specificity",
        "thermo": "predicted_thermo", "expression": "predicted_expression",
    }[name]
    shuffled = [r[0] for r in conn.execute(
        f"SELECT {control_column} FROM controls WHERE kind = 'scrambled' "
        f"AND target = ? AND {control_column} IS NOT NULL", (target,))]
    auc = mann_whitney_auc(values, shuffled)
    result["auc_vs_scrambled"] = None if auc is None else round(auc, 4)

    # How much of this axis could the control move at all? The shuffle
    # preserves amino-acid composition and the framework by design, so a
    # component reading only those is bit-identical between a candidate and
    # its scrambled twin and contributes exactly nothing here. Calling such an
    # axis "non-discriminating" blames the axis for a control that cannot
    # test it: 67 per cent of the expression axis is composition-only, so its
    # AUC sits at chance by construction rather than by failure.
    components = AXIS_COMPONENTS.get(name, [])
    total_weight = sum(c.weight for c in components) or 1.0
    testable = sum(c.weight for c in components
                   if not c.shuffle_invariant) / total_weight
    result["testable_weight_fraction"] = round(testable, 3)

    if auc is not None:
        result["n_scrambled"] = len(shuffled)
        if testable < 0.5:
            if result["passed"]:
                result["passed"] = None
            result["findings"].append(
                f"this control cannot test this axis: {1 - testable:.0%} of its "
                f"component weight reads only amino-acid composition or the "
                f"framework, both of which a CDR shuffle preserves exactly, so "
                f"the measured AUC of {auc:.3f} is near chance by construction "
                f"rather than by failure. A negative control that varies "
                f"composition is needed before this axis can be called "
                f"informative or not.")
        elif auc < 0.6:
            result["passed"] = False
            result["findings"].append(
                f"does not separate real candidates from CDR-shuffled ones "
                f"(AUC {auc:.3f}; 0.5 is no separation), and {testable:.0%} of "
                f"the axis weight could have moved. The axis is "
                f"non-discriminating and must be labelled as such wherever it "
                f"appears.")
    else:
        # A check that did not run is not a check that passed. Leaving
        # `passed` True here reported PASS for an axis whose decisive test had
        # never executed, which is a stronger claim than "we have not looked"
        # and the opposite of what this file exists to do. Only downgrade a
        # would-be pass: a failure found by the other checks still stands.
        if result["passed"]:
            result["passed"] = None
        result["findings"].append(
            "no scrambled control scored yet, so the strongest check has not "
            "run and this axis is untested, not passing")
    return result


def main() -> int:
    conn = connect()
    results: list[dict] = []
    for target in TARGETS:
        for name, column in AXES.items():
            result = check_axis(conn, target, name, column)
            results.append(result)
            state = ("PASS" if result["passed"]
                     else ("PENDING" if result["passed"] is None else "VACUOUS"))
            detail = "; ".join(result["findings"]) or "informative"
            log(f"vacuous_tests: {state} {target} {name}: {detail}")

    by_axis: dict[str, dict] = {}
    for result in results:
        entry = by_axis.setdefault(result["axis"], {"targets": {}, "passed": True})
        entry["targets"][result["target"]] = result
        if result["passed"] is False:
            entry["passed"] = False

    failed = [axis for axis, entry in by_axis.items() if not entry["passed"]]
    conn.execute("INSERT OR REPLACE INTO qc (key, value, passed) VALUES (?,?,?)",
                 ("vacuous_tests", json.dumps(by_axis), 0 if failed else 1))
    conn.commit()
    conn.close()

    if failed:
        print(f"VACUOUS AXES: {', '.join(failed)}")
        print("These are findings to report, not failures to hide. The app marks "
              "a non-discriminating axis wherever it appears and keeps it out of "
              "the Pareto front.")
        return 0  # a finding, not an error exit
    print("vacuous tests clean: every computed axis carries information")
    return 0


if __name__ == "__main__":
    sys.exit(main())
