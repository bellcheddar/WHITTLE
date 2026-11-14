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

   This is asked twice, because the two forms answer different questions and
   on this data they disagree. **Unpaired** compares every scrambled control
   against all 500 candidates: it asks where a shuffled sequence would rank,
   and it is the number spec Section 9 asks for. **Paired** compares each
   control against the one candidate it was shuffled from, which is the same
   sequence with the same composition, length and framework, so CDR order is
   the only difference left inside a pair. Between-sequence variance drops out
   of the paired form, and with 50 pairs per target it is far the more
   powerful of the two. The paired form therefore decides the verdict and the
   unpaired AUC is still reported beside it.
5. **Cannot fail by construction.** The axis has no low tail: nothing scores
   badly, so the column cannot discriminate against anything.

Results go to the `qc` table and onto the QC panel on the control tab.

    pixi run python tools/vacuous_tests.py
"""
from __future__ import annotations

import json
import math
import re
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


def sign_test(wins: int, n: int) -> float | None:
    """Two-sided exact binomial test against p = 0.5. No scipy on the droplet.

    A p-value is the right statistic here where an AUC is not, because the
    paired question is a count of 50 and not a comparison of two large
    samples: 30 of 50 looks like separation and is what a fair coin does one
    time in five.

    Both tails are taken, and the smaller one is doubled. Taking only the upper
    tail, as this did until 2026-11-14, makes the test blind to the one result
    it most needs to report: a metric that runs backwards. `bsa_total` has the
    candidate beating the control shuffled from it on 10 of 50 pairs, which is
    p = 2.4e-05 that the shuffle is better, and the upper tail alone scored it
    p = 1 and read as no evidence of anything. A metric pointing the wrong way
    is worse than a metric carrying no signal, because its weight actively
    moves the ranking away from the candidates it should favour, so reporting
    the two the same way was the dangerous direction to be wrong in.
    """
    if n == 0:
        return None
    upper = sum(math.comb(n, k) for k in range(wins, n + 1)) / 2 ** n
    lower = sum(math.comb(n, k) for k in range(0, wins + 1)) / 2 ** n
    return min(1.0, 2 * min(upper, lower))


def paired_against_parents(conn, target: str, control_column: str,
                           axis_column: str) -> dict:
    """Does a candidate beat the control that was shuffled from it?

    The parent is recovered from the control's own note rather than stored as a
    column, because `pipeline.controls` writes the provenance there and
    duplicating it in the schema would let the two disagree. A control whose
    parent is missing or unscored is dropped from the count rather than
    counted as a loss.

    Both sides are on the same scale: `pipeline.controls_run` ranks a control
    against its target's candidate distribution using scales fitted on the
    candidates alone, which is exactly the scale `scores.axis_*` carries.
    """
    wins = losses = ties = 0
    for row in conn.execute(
            f"SELECT notes, {control_column} AS value FROM controls "
            "WHERE kind = 'scrambled' AND target = ?", (target,)):
        if row["value"] is None:
            continue
        found = re.search(r"derived from ([A-Za-z0-9]+-\d+)", row["notes"] or "")
        if not found:
            continue
        parent = conn.execute(
            f"SELECT {axis_column} AS value FROM scores WHERE candidate_id = ?",
            (found.group(1),)).fetchone()
        if parent is None or parent["value"] is None:
            continue
        if parent["value"] > row["value"]:
            wins += 1
        elif parent["value"] < row["value"]:
            losses += 1
        else:
            ties += 1
    n = wins + losses + ties
    out = {"paired_n": n, "paired_wins": wins, "paired_ties": ties}
    if n < 10:
        out["paired_fraction"] = None
        out["paired_p"] = None
        return out
    out["paired_fraction"] = round(wins / n, 3)
    out["paired_p"] = sign_test(wins, n)
    return out


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

    # The same question asked pairwise, which is the powerful form.
    result.update(paired_against_parents(conn, target, control_column, column))

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

    result["n_scrambled"] = len(shuffled)
    # The separation verdict is NOT decided here. It is decided once per axis
    # over all five targets by `pooled_separation`, because fifty pairs cannot
    # resolve the effect these axes actually have: at a true 60 per cent rate
    # the expected p-value at n = 50 is 0.20, so a per-target verdict is a coin
    # toss dressed as a finding. It published specificity as non-discriminating
    # on EGFR and TNFa, unconfirmed on GLP1R and TL1A and informative on PRLR,
    # which is one 60 per cent effect wearing three labels. The per-target
    # counts are kept and displayed; what they no longer do is decide.
    #
    # The structural checks above stay per target, because no-variance, modal
    # degeneracy and restatement of another axis are genuinely properties of
    # one target's distribution and need no pooling to be seen.
    if auc is None:
        # A check that did not run is not a check that passed. Leaving `passed`
        # True here reported PASS for an axis whose decisive test had never
        # executed, which is a stronger claim than "we have not looked".
        if result["passed"]:
            result["passed"] = None
        result["findings"].append(
            "no scrambled control scored yet, so the strongest check has not "
            "run and this axis is untested, not passing")
    return result


def pooled_separation(axis: str, results: list[dict]) -> dict:
    """One separation verdict per axis, decided over all five targets.

    Fifty pairs per target cannot resolve the effect sizes these axes have.
    At a true 60 per cent win rate the expected sign-test p-value at n = 50 is
    0.20, so asking each target separately returns noise: specificity came
    back non-discriminating on EGFR (28 of 46) and TNFa (24 of 49), weak on
    GLP1R (30 of 47) and TL1A (30 of 50) and informative on PRLR (32 of 48),
    which is a single 60 per cent effect labelled three different ways by luck.
    Pooled it is 144 of 240, p = 0.0024: weak, and real. Marc's call,
    2026-11-14 (DECISIONS 53).

    Pooling settles only whether the axis carries information. It does not
    pool the rankings, which stay within-target percentiles throughout, and it
    does not hide the per-target counts, which travel with the verdict so a
    reader can see both the answer and how thin the evidence is per target.

    Pooling across targets is sound here because both sides of every pair are
    already within-target percentiles on the same 0 to 100 scale, so a win is
    the same event on EGFR as on TNFa. What it assumes is that the axis
    behaves alike across targets, which n = 50 per target cannot confirm; the
    assumption is stated rather than buried, and the fix is more controls.
    """
    wins = sum(r.get("paired_wins") or 0 for r in results
               if r.get("paired_p") is not None)
    n = sum(r.get("paired_n") or 0 for r in results
            if r.get("paired_p") is not None)
    aucs = [r["auc_vs_scrambled"] for r in results
            if r.get("auc_vs_scrambled") is not None]
    testable = next((r["testable_weight_fraction"] for r in results
                     if r.get("testable_weight_fraction") is not None), 1.0)
    per_target = {r["target"]: f"{r.get('paired_wins')}/{r.get('paired_n')}"
                  for r in results if r.get("paired_p") is not None}

    out = {"pooled_wins": wins, "pooled_n": n, "per_target": per_target,
           "testable_weight_fraction": testable,
           "auc_range": [min(aucs), max(aucs)] if aucs else None}
    if not n:
        out["verdict"] = "pending"
        out["finding"] = ("no scrambled control has been scored on this axis, "
                          "so its strongest check has not run")
        return out

    fraction = wins / n
    p = sign_test(wins, n)
    out["pooled_fraction"] = round(fraction, 3)
    out["pooled_p"] = p
    # Only say the per-target counts are unresolvable when they are. Thermo
    # clears p < 0.05 on every target by itself; specificity clears it on none.
    resolved = [r["target"] for r in results
                if r.get("paired_p") is not None and r["paired_p"] < 0.05]
    if not resolved:
        tail = ", none of which reaches significance alone at this size."
    elif len(resolved) == len(per_target):
        tail = ", every one of which is significant on its own."
    else:
        tail = (f", of which {', '.join(sorted(resolved))} "
                f"{'is' if len(resolved) == 1 else 'are'} significant alone.")
    spread = (" Per target: "
              + ", ".join(f"{t} {c}" for t, c in sorted(per_target.items()))
              + tail)

    if testable < 0.5:
        out["verdict"] = "untestable"
        out["finding"] = (
            f"this control cannot test this axis: {1 - testable:.0%} of its "
            f"component weight reads only amino-acid composition or the "
            f"framework, both of which a CDR shuffle preserves exactly, so a "
            f"result near chance is expected by construction rather than by "
            f"failure. No synthetic control settles it either, because the two "
            f"heaviest components are defined from composition and a "
            f"deliberately insoluble sequence would score badly by arithmetic "
            f"and not by evidence. Deciding this axis needs measured "
            f"expression or solubility on real nanobodies, of which this "
            f"control set has none, so it stays untested rather than "
            f"informative or non-discriminating. Pooled over all five targets "
            f"it is {wins} of {n} ({fraction:.0%}, p = {p:.2g}), which is what "
            f"a control that cannot move two thirds of the weight looks like.")
    elif p < 0.05 and fraction > 0.5:
        # "weak" is a verdict of its own rather than a modifier on
        # "informative", so that one value drives every label in the app. A
        # 60 per cent axis and an 84 per cent axis both separate, and telling a
        # reader to weigh them equally would be the overstatement this file
        # exists to prevent.
        clear = fraction >= 0.70
        strength = "clearly" if clear else "weakly but measurably"
        out["verdict"] = "informative" if clear else "weak"
        out["strength"] = "clear" if clear else "weak"
        out["finding"] = (
            f"{strength} separates a candidate from its own CDR-shuffled twin: "
            f"pooled over all five targets, {wins} of {n} pairs "
            f"({fraction:.0%}, sign test p = {p:.2g}).{spread}")
    elif p < 0.05 and fraction < 0.5:
        # Reachable only since the sign test started reading both tails
        # (DECISIONS 51). An inverted axis is worse than an empty one, because
        # its weight moves the ranking away from the candidates it should
        # favour, so it gets a label of its own rather than sharing one with
        # "carries no information".
        out["verdict"] = "inverted"
        out["finding"] = (
            f"runs backwards: the CDR-shuffled control beats the candidate it "
            f"was shuffled from on {n - wins} of {n} pairs (p = {p:.2g}). This "
            f"is not an absence of signal, it is signal with the wrong sign, "
            f"and nothing should be ranked on it until the cause is "
            f"understood.{spread}")
    else:
        out["verdict"] = "vacuous"
        out["finding"] = (
            f"does not separate real candidates from CDR-shuffled ones: "
            f"pooled over all five targets, {wins} of {n} pairs "
            f"({fraction:.0%}, sign test p = {p:.2g}), with {testable:.0%} of "
            f"the axis weight free to move. At this pooled size the test would "
            f"have found a 60 per cent effect, so this is a real negative "
            f"rather than a shortage of controls. The axis is "
            f"non-discriminating and is labelled as such wherever it "
            f"appears.{spread}")
    return out


# "weak" passes: it separates, by a small margin that the label states.
PASSES = {"informative", "weak"}
FAILS = {"vacuous", "inverted"}


def main() -> int:
    conn = connect()
    results: list[dict] = []
    for target in TARGETS:
        for name, column in AXES.items():
            results.append(check_axis(conn, target, name, column))

    by_axis: dict[str, dict] = {}
    for result in results:
        entry = by_axis.setdefault(result["axis"], {"targets": {}})
        entry["targets"][result["target"]] = result

    # One separation verdict per axis, over all five targets. The per-target
    # structural checks still stand on their own and can fail a single target
    # by themselves; what pooling replaces is the coin toss that fifty pairs
    # made of the separation question. See `pooled_separation`.
    for axis, entry in by_axis.items():
        per = list(entry["targets"].values())
        pooled = pooled_separation(axis, per)
        entry["pooled"] = pooled
        verdict = pooled["verdict"]
        entry["verdict"] = verdict
        entry["passed"] = (True if verdict in PASSES
                           else False if verdict in FAILS else None)
        # A structural failure on one target is not overridden by a pooled
        # pass: it is a different fault and it still stands.
        structural = [t for t, r in entry["targets"].items()
                      if r["passed"] is False]
        if structural:
            entry["passed"] = False
            entry["structural_failures"] = structural
        for result in per:
            result["verdict"] = verdict
            result["pooled_finding"] = pooled["finding"]
        # Kept for the app, which renders a label per target: they now agree
        # by construction, which is the point.
        entry["verdicts"] = {t: verdict for t in entry["targets"]}

    for axis, entry in sorted(by_axis.items()):
        pooled = entry["pooled"]
        state = ("PASS" if entry["passed"]
                 else "PENDING" if entry["passed"] is None else "FAIL")
        counts = (f"pooled {pooled['pooled_wins']}/{pooled['pooled_n']}, "
                  f"p = {pooled['pooled_p']:.2g}"
                  if pooled.get("pooled_p") is not None else "no controls")
        log(f"vacuous_tests: {state} {axis} [{counts}] "
            f"{pooled['verdict'].upper()}: {pooled['finding']}")
        for target, result in sorted(entry["targets"].items()):
            own = [f for f in result["findings"]]
            if own:
                log(f"vacuous_tests:   {target} also: {'; '.join(own)}")

    failed = [axis for axis, entry in by_axis.items() if entry["passed"] is False]
    unresolved = [axis for axis, entry in by_axis.items() if entry["passed"] is None]
    conn.execute("INSERT OR REPLACE INTO qc (key, value, passed) VALUES (?,?,?)",
                 ("vacuous_tests", json.dumps(by_axis), 0 if failed else 1))
    conn.commit()
    conn.close()

    if failed:
        print(f"VACUOUS AXES: {', '.join(failed)}")
        print("These are findings to report, not failures to hide. The app marks "
              "a non-discriminating axis wherever it appears.")
        print("An axis that fails is excluded from the ranking, which is what "
              "run.tier2_reorders_shortlist does for affinity. The verdict is "
              "now one per axis rather than one per target (DECISIONS 53), so "
              "there is no longer a case where an axis would be ranked on for "
              "some targets and not others.")
        return 0  # a finding, not an error exit
    if unresolved:
        print(f"UNRESOLVED: {', '.join(unresolved)} were tested and the control "
              f"set could not rule on them. Not the same as passing.")
        return 0
    print("vacuous tests clean: every computed axis carries information")
    return 0


if __name__ == "__main__":
    sys.exit(main())
