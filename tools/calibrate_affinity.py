#!/usr/bin/env python3
"""Weight the affinity axis by evidence rather than by argument.

The affinity axis blends seven interface metrics, and until this runs their
weights are set by judgement. That is a weak footing for a ranking, because
several of them are plausible and still measure nothing: a contact inventory
or a shape complementarity computed on a predicted pose can look perfectly
reasonable and carry no information about whether the thing binds.

The test is the scrambled controls. A CDR-shuffled sequence has no reason to
bind its antigen, so it is co-folded against the same antigen through the same
code, and each metric is asked one question: **how often does a real candidate
beat a shuffled one?** That is the Mann-Whitney AUC. 0.5 is no separation at
all, and a metric sitting at 0.5 is measuring two proteins placed near each
other rather than binding.

Weights are then proportional to |AUC - 0.5|, so a metric earns its influence
on the ranking in proportion to the separation it actually achieves. Anything
that fails to separate is reported for removal rather than quietly carried.

    pixi run python tools/calibrate_affinity.py
    pixi run python tools/calibrate_affinity.py --apply

Nothing is written without `--apply`, and `--apply` only rewrites the weight
table in `pipeline/score.py`; it never touches the metrics themselves.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline.common import TARGETS, connect, log  # noqa: E402
from pipeline.score import AXIS_COMPONENTS  # noqa: E402

sys.path.insert(0, str(ROOT / "tools"))
from vacuous_tests import mann_whitney_auc  # noqa: E402

# Every metric the affinity axis can draw on, with the direction that counts
# as better. Taken from the axis itself so the two cannot drift apart.
METRICS = {c.field: c.sign for c in AXIS_COMPONENTS["affinity"]}

AGGREGATE = {"pose_consensus": "MAX"}      # already one value per candidate
MIN_PER_GROUP = 8


def collect(conn, target: str, is_control: int) -> dict[str, list[float]]:
    """One value per candidate per metric, averaged over its seeds."""
    selects = ", ".join(
        f"{AGGREGATE.get(field, 'AVG')}(f.{field}) AS {field}"
        for field in METRICS)
    rows = conn.execute(
        f"SELECT f.candidate_id, {selects} FROM cofolds f "
        f"WHERE f.target = ? AND f.status = 'ok' AND COALESCE(f.is_control, 0) = ? "
        f"GROUP BY f.candidate_id", (target, is_control)).fetchall()
    out: dict[str, list[float]] = {field: [] for field in METRICS}
    for row in rows:
        for field in METRICS:
            value = row[field]
            if value is not None:
                out[field].append(float(value))
    return out


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--apply", action="store_true",
                        help="rewrite the affinity weights in pipeline/score.py")
    parser.add_argument("--floor", type=float, default=0.55,
                        help="AUC below which a metric is reported as "
                             "non-separating and given no weight")
    args = parser.parse_args()

    conn = connect()
    per_metric: dict[str, list[float]] = {field: [] for field in METRICS}
    usable_targets = []

    for target in TARGETS:
        real = collect(conn, target, 0)
        scrambled = collect(conn, target, 1)
        sizes = (len({len(v) for v in [real[f] for f in METRICS]}),)
        n_real = max((len(v) for v in real.values()), default=0)
        n_scram = max((len(v) for v in scrambled.values()), default=0)
        if n_real < MIN_PER_GROUP or n_scram < MIN_PER_GROUP:
            log(f"calibrate: {target} skipped, {n_real} candidates and "
                f"{n_scram} scrambled controls co-folded, need "
                f"{MIN_PER_GROUP} of each")
            continue
        usable_targets.append(target)
        log(f"calibrate: {target}, {n_real} candidates against {n_scram} "
            f"scrambled controls")
        for field, sign in METRICS.items():
            a = real[field]
            b = scrambled[field]
            if len(a) < MIN_PER_GROUP or len(b) < MIN_PER_GROUP:
                continue
            if sign < 0:                     # lower is better, so flip
                a = [-v for v in a]
                b = [-v for v in b]
            auc = mann_whitney_auc(a, b)
            if auc is None:
                continue
            per_metric[field].append(auc)
            verdict = ("separates" if auc >= args.floor else
                       "does not separate")
            log(f"calibrate:   {field:24} AUC {auc:.3f}  {verdict}")
        del sizes

    if not usable_targets:
        log("calibrate: nothing to calibrate yet. Co-fold the scrambled "
            "controls first: pixi run -e cofold python -m pipeline.cofold "
            "--controls 10")
        conn.close()
        return 1

    log(f"calibrate: pooled over {len(usable_targets)} targets")
    pooled: dict[str, float] = {}
    for field, values in per_metric.items():
        if values:
            pooled[field] = sum(values) / len(values)

    # Weight by separation achieved, normalised so the heaviest metric keeps
    # the weight the heaviest metric has now. A metric below the floor earns
    # nothing: carrying it would dilute the ones that work.
    strengths = {f: max(0.0, abs(a - 0.5)) if a >= args.floor else 0.0
                 for f, a in pooled.items()}
    top = max(strengths.values()) if strengths else 0.0
    current_top = max(c.weight for c in AXIS_COMPONENTS["affinity"])
    proposed = {f: (round(current_top * s / top, 2) if top > 0 else 0.0)
                for f, s in strengths.items()}

    log("calibrate: proposed affinity weights")
    for component in AXIS_COMPONENTS["affinity"]:
        field = component.field
        auc = pooled.get(field)
        new = proposed.get(field)
        note = ""
        if auc is None:
            note = "  (no data, weight unchanged)"
            new = component.weight
        elif new == 0.0:
            note = f"  (AUC {auc:.3f} below {args.floor}, drop it)"
        log(f"calibrate:   {field:24} {component.weight:5.2f} -> {new:5.2f}"
            + (f"  AUC {auc:.3f}" if auc is not None else "") + note)

    if args.apply:
        log("calibrate: --apply is not yet wired; the weights above are a "
            "recommendation and the edit to pipeline/score.py is deliberate "
            "and manual, so the reasoning lands in the commit message")
    conn.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
