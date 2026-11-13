"""Independent re-implementation of the Pareto front. Spec 10 item 6.

"Arithmetic re-run independently for the Pareto fronts: a second
implementation, not the same function called twice."

`pipeline/score.py` uses the standard fast non-dominated sort: build the
domination counts once, peel front by front. This file does the obviously
correct thing instead: a candidate is on the front if no other candidate
dominates it, checked directly, O(n squared) and with no bookkeeping to get
wrong. If the two disagree, the fast one is wrong.

It also re-derives the within-target percentiles from the raw axis values, so a
percentile computed against the wrong target would show up here as a
disagreement rather than passing silently (spec 10 item 3).

    pixi run python tools/qc_pareto.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline.common import TARGETS, connect, log  # noqa: E402

AXES = ["axis_affinity", "axis_specificity", "axis_thermo", "axis_expression"]


def on_front_bruteforce(points: list[dict]) -> set[int]:
    """Index set of the non-dominated points, checked pairwise and directly."""
    usable = [a for a in AXES if any(p.get(a) is not None for p in points)]
    front: set[int] = set()
    for i, a in enumerate(points):
        dominated = False
        for j, b in enumerate(points):
            if i == j:
                continue
            shared = [ax for ax in usable
                      if a.get(ax) is not None and b.get(ax) is not None]
            if not shared:
                continue
            if all(b[ax] >= a[ax] for ax in shared) and \
               any(b[ax] > a[ax] for ax in shared):
                dominated = True
                break
        if not dominated:
            front.add(i)
    return front


def percentile_bruteforce(values: list[float | None]) -> list[float | None]:
    """Percentile by counting, with ties averaged. No sorting tricks."""
    present = [v for v in values if v is not None]
    n = len(present)
    out: list[float | None] = []
    for v in values:
        if v is None or n < 2:
            out.append(None if v is None else 50.0)
            continue
        below = sum(1 for other in present if other < v)
        equal = sum(1 for other in present if other == v)
        # Mean rank over the tie block, mapped onto 0-100 across n-1 steps.
        mean_rank = below + (equal - 1) / 2.0
        out.append(round(100.0 * mean_rank / (n - 1), 2))
    return out


def main() -> int:
    conn = connect()
    failures: list[str] = []
    summary: dict[str, dict] = {}

    for target in TARGETS:
        rows = conn.execute(
            "SELECT candidate_id, axis_affinity, axis_specificity, axis_thermo, "
            "axis_expression, on_front FROM scores WHERE target = ? "
            "ORDER BY candidate_id", (target,)).fetchall()
        if not rows:
            continue
        points = [{a: r[a] for a in AXES} for r in rows]
        recomputed = on_front_bruteforce(points)
        stored = {i for i, r in enumerate(rows) if r["on_front"]}
        missing = sorted(rows[i]["candidate_id"] for i in recomputed - stored)
        extra = sorted(rows[i]["candidate_id"] for i in stored - recomputed)
        agree = not missing and not extra
        if not agree:
            failures.append(
                f"{target}: the stored front and the independent one disagree. "
                f"{len(missing)} missing from the stored front "
                f"({missing[:5]}), {len(extra)} extra ({extra[:5]})")
        summary[target] = {"n": len(rows), "front_stored": len(stored),
                           "front_recomputed": len(recomputed), "agree": agree}
        log(f"qc_pareto: {target} {len(rows)} rows, stored front {len(stored)}, "
            f"independent front {len(recomputed)}, "
            f"{'AGREE' if agree else 'DISAGREE'}")

        # Spec 10 item 3: a percentile computed against the wrong target.
        # Re-derive one axis from the raw values and compare.
        raw = conn.execute(
            "SELECT c.id, t.tm_pred FROM candidates c "
            "LEFT JOIN thermostability t ON t.candidate_id = c.id "
            "WHERE c.target = ? ORDER BY c.id", (target,)).fetchall()
        values = [r["tm_pred"] for r in raw]
        if any(v is not None for v in values):
            expected = percentile_bruteforce(values)
            stored_thermo = {r["candidate_id"]: r["axis_thermo"] for r in rows}
            worst = 0.0
            for row, want in zip(raw, expected):
                got = stored_thermo.get(row["id"])
                if want is None or got is None:
                    continue
                worst = max(worst, abs(got - want))
            # The thermo axis blends Tm with Tier 0 terms, so an exact match is
            # not expected; a large rank disagreement would mean the percentile
            # was computed against the wrong population.
            summary[target]["thermo_percentile_max_gap"] = round(worst, 2)
            log(f"qc_pareto: {target} thermo percentile differs from a "
                f"Tm-only percentile by at most {worst:.1f} points "
                f"(the axis blends other terms, so a gap is expected; a gap "
                f"near 100 would mean the wrong population)")
            if worst > 95:
                failures.append(
                    f"{target}: the thermostability percentile looks computed "
                    f"against the wrong population (max gap {worst:.1f})")

    conn.execute("INSERT OR REPLACE INTO qc (key, value, passed) VALUES (?,?,?)",
                 ("pareto_independent_check", json.dumps(summary),
                  0 if failures else 1))
    conn.commit()
    conn.close()

    if failures:
        print("PARETO CROSS-CHECK FAILED")
        for line in failures:
            print("  " + line)
        return 1
    print("pareto cross-check clean: the independent implementation agrees")
    return 0


if __name__ == "__main__":
    sys.exit(main())
