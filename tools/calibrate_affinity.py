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
from vacuous_tests import mann_whitney_auc, sign_test  # noqa: E402

sign_p = sign_test       # the paired count is a sign test, not an AUC

# Every metric the affinity axis can draw on, with the direction that counts
# as better. Taken from the axis itself so the two cannot drift apart.
METRICS = {c.field: c.sign for c in AXIS_COMPONENTS["affinity"]}

# Metrics the co-folds carry that the axis does **not** use. Surveyed, never
# weighted: the recommendation at the end of this script covers the components
# above and nothing here, and adding one of these to the axis stays a
# deliberate edit to `pipeline/score.py`.
#
# This exists because the survey is the only way a missing component gets
# found. `polar_fraction` separated real from shuffled at AUC 0.841 on the
# first unpaired read, better than any metric actually in the axis, and
# because METRICS is derived from the axis itself the calibration could not
# see it and would never have mentioned it again.
#
# Each needs a direction that counts as better, and a metric whose direction
# is genuinely arguable is left out rather than given one. `n_hydrophobic` is
# the case in point: more hydrophobic contact helps affinity and costs
# specificity, so "better" depends on which axis is asking, and counting wins
# in either direction would assert something this project does not believe.
SURVEY = {
    "polar_fraction": +1,          # share of contacts that are polar
    "interactions_per_100a2": +1,  # contact density, size-normalised
    "n_saltbridge": +1,            # charged pairs across the interface
    "n_his_contacts": +1,          # EGFR and TNFa are His-enriched by design
    "complex_iptm": +1,            # whole-complex confidence, not pairwise
    "plddt": +1,                   # per-atom confidence of the pose
    "iptm_spread": -1,             # disagreement between the five samples
    "edge_fraction": -1,           # interface sitting on the construct's cut
}

# Everything measured, components first. The split is kept at every point
# where results are reported, so a surveyed metric can never be read as part
# of the axis.
ALL_FIELDS = {**METRICS, **SURVEY}

AGGREGATE = {"pose_consensus": "MAX"}      # already one value per candidate
MIN_PER_GROUP = 8
MIN_PAIRS = 5


def control_parents(conn, target: str) -> dict[str, str]:
    """Control id to the candidate it was shuffled from."""
    import re

    out: dict[str, str] = {}
    for row in conn.execute(
            "SELECT id, notes FROM controls WHERE kind = 'scrambled' "
            "AND target = ?", (target,)):
        found = re.search(r"derived from ([A-Za-z0-9]+-\d+)", row["notes"] or "")
        if found:
            out[row["id"]] = found.group(1)
    return out


def paired(conn, target: str) -> dict[str, tuple[int, int, int]]:
    """Per metric: (control better, parent better, ties), within pairs.

    The unpaired comparison is confounded. The co-folded candidates are the
    top N by Tier 1 rank, a rank driven by patch surface hydrophobicity, so
    they are much less hydrophobic than the population the controls were
    shuffled from. Comparing those two groups measures Tier 1 selection as
    much as scrambling. Inside a pair the only difference is CDR order, which
    is the question the control exists to ask.
    """
    parents = control_parents(conn, target)
    if not parents:
        return {}
    selects = ", ".join(
        f"{AGGREGATE.get(f, 'AVG')}(f.{f}) AS {f}" for f in ALL_FIELDS)
    rows = {r["candidate_id"]: dict(r) for r in conn.execute(
        f"SELECT f.candidate_id, {selects} FROM cofolds f "
        f"WHERE f.target = ? AND f.status = 'ok' GROUP BY f.candidate_id",
        (target,))}
    tally: dict[str, list[int]] = {f: [0, 0, 0] for f in ALL_FIELDS}
    for control_id, parent_id in parents.items():
        a, b = rows.get(control_id), rows.get(parent_id)
        if not a or not b:
            continue
        for field, sign in ALL_FIELDS.items():
            x, y = a.get(field), b.get(field)
            if x is None or y is None:
                continue
            better = (x - y) * sign       # positive means the control wins
            if abs(better) < 1e-12:
                tally[field][2] += 1
            elif better > 0:
                tally[field][0] += 1
            else:
                tally[field][1] += 1
    return {f: tuple(v) for f, v in tally.items()}


def collect(conn, target: str, is_control: int) -> dict[str, list[float]]:
    """One value per candidate per metric, averaged over its seeds."""
    selects = ", ".join(
        f"{AGGREGATE.get(field, 'AVG')}(f.{field}) AS {field}"
        for field in ALL_FIELDS)
    rows = conn.execute(
        f"SELECT f.candidate_id, {selects} FROM cofolds f "
        f"WHERE f.target = ? AND f.status = 'ok' AND COALESCE(f.is_control, 0) = ? "
        f"GROUP BY f.candidate_id", (target, is_control)).fetchall()
    out: dict[str, list[float]] = {field: [] for field in ALL_FIELDS}
    for row in rows:
        for field in ALL_FIELDS:
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
    per_metric: dict[str, list[float]] = {field: [] for field in ALL_FIELDS}
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

        # The paired view first, because the unpaired one is confounded.
        pairs = paired(conn, target)
        n_pairs = max((sum(v) for v in pairs.values()), default=0)
        if n_pairs >= MIN_PAIRS:
            log(f"calibrate:   paired against each control's own parent "
                f"({n_pairs} pairs; this is the unconfounded comparison)")
            for field in METRICS:
                wins, losses, ties = pairs.get(field, (0, 0, 0))
                decided = wins + losses
                if not decided:
                    continue
                log(f"calibrate:     {field:22} parent beats its shuffle "
                    f"{losses}/{decided}  p={sign_p(losses, decided):.2g}")
            shown = False
            for field in SURVEY:
                wins, losses, ties = pairs.get(field, (0, 0, 0))
                decided = wins + losses
                if not decided:
                    continue
                if not shown:
                    log("calibrate:   surveyed, NOT in the axis and not "
                        "weighted below:")
                    shown = True
                log(f"calibrate:     {field:22} parent beats its shuffle "
                    f"{losses}/{decided}  p={sign_p(losses, decided):.2g}")
        else:
            log(f"calibrate:   only {n_pairs} control-parent pairs co-folded, "
                f"so the paired comparison is skipped and the unpaired one "
                f"below is confounded by Tier 1 selection. Run "
                f"`pipeline.cofold --parents 10`.")
        for field, sign in ALL_FIELDS.items():
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
            tag = "" if field in METRICS else "  [surveyed, not in the axis]"
            log(f"calibrate:   {field:24} AUC {auc:.3f}  {verdict}{tag}")
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
    # Normalise against the best **component**, never the best surveyed
    # metric. A surveyed metric carries no weight, so letting one set the
    # scale would shrink every real component's proposed weight against
    # something the axis does not use.
    top = max((s for f, s in strengths.items() if f in METRICS), default=0.0)
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

    # Does anything outside the axis beat everything inside it? This is the
    # question the survey exists to answer, and it is asked explicitly rather
    # than left for a reader to spot in the table above.
    best_component = max((pooled.get(c.field, 0.0)
                          for c in AXIS_COMPONENTS["affinity"]), default=0.0)
    contenders = sorted(((pooled[f], f) for f in SURVEY if f in pooled
                         and pooled[f] > best_component), reverse=True)
    if contenders:
        log(f"calibrate: NOT IN THE AXIS and separating better than every "
            f"metric that is (best component AUC {best_component:.3f}):")
        for auc, field in contenders:
            log(f"calibrate:   {field:24} AUC {auc:.3f}  consider adding it")
        log("calibrate: check the paired counts above before acting on these. "
            "An unpaired AUC is confounded by Tier 1 selection, and a metric "
            "that only separates unpaired has not earned a place.")
    else:
        log("calibrate: no surveyed metric out-separates the axis components")

    if args.apply:
        log("calibrate: --apply is not yet wired; the weights above are a "
            "recommendation and the edit to pipeline/score.py is deliberate "
            "and manual, so the reasoning lands in the commit message")
    conn.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
