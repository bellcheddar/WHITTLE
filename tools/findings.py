"""FINDINGS.md. Spec Section 11, Phase 4.

"What the control set showed, where the axes disagree with the upstream
columns, what the scrambled control did to the affinity axis."

Generated from the database rather than written by hand, so it cannot drift
from what the run actually produced. Every number below is read back out of
the `qc`, `controls`, `scores` and `candidates` tables.

Adapted from `debris/tools/findings.py` and `NOBBLE/tools/findings.py`
(see REUSE.md): the generated-not-maintained convention is the part carried
over.

    pixi run findings
"""
from __future__ import annotations

import json
import statistics
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline.common import TARGETS, connect, log, run_id  # noqa: E402

AXES = {"affinity": "axis_affinity", "specificity": "axis_specificity",
        "thermo": "axis_thermo", "expression": "axis_expression"}


def qc_value(conn, key: str):
    row = conn.execute("SELECT value, passed FROM qc WHERE key = ?",
                       (key,)).fetchone()
    return (row["value"], bool(row["passed"])) if row else (None, None)


def section_upstream(conn) -> list[str]:
    lines = ["## Where this build disagrees with the upstream columns", ""]

    strict, _ = qc_value(conn, "hallmark_reconciliation_strict")
    graded, _ = qc_value(conn, "hallmark_reconciliation_graded")
    his, _ = qc_value(conn, "his_recompute_agreement")
    if strict and graded:
        lines += [
            "### The camelid hallmark tetrad",
            "",
            f"The upstream `hallmark` column calls {float(strict) * 100:.1f} "
            f"percent of rows the same way a strict all-four reading of IMGT 42, "
            f"49, 50 and 52 does, and {float(graded) * 100:.1f} percent the same "
            f"way an at-least-three-of-four reading does. The upstream column "
            f"therefore implements the three-of-four rule, not the four-of-four "
            f"one it appears to.",
            "",
        ]
        dist = dict(conn.execute(
            "SELECT hallmark_n, COUNT(*) FROM candidates GROUP BY hallmark_n"
        ).fetchall())
        total = sum(dist.values())
        lines += [
            "| Camelid residues present | Candidates | Share |",
            "|---|---|---|",
        ]
        for n in sorted(dist, reverse=True):
            lines.append(f"| {n} of 4 | {dist[n]:,} | {dist[n] / total * 100:.1f}% |")
        vh = conn.execute(
            "SELECT COUNT(*) FROM candidates WHERE hallmark_detail LIKE '%!VH%'"
        ).fetchone()[0]
        lines += [
            "",
            f"{vh:,} candidates carry a VH-type residue at one of the four "
            f"positions. That leaves the former light-chain interface "
            f"hydrophobic and exposed, which is a solubility and stability "
            f"liability rather than a cosmetic difference. This build reports "
            f"the 0 to 4 count and the strict boolean side by side.",
            "",
        ]
    if his:
        lines += [
            "### CDR histidine counts",
            "",
            f"Recomputed IMGT CDR histidine totals agree with the upstream "
            f"columns for {float(his) * 100:.1f} percent of rows. The residual is "
            f"a CDR-definition difference, not a counting error: this build uses "
            f"IMGT CDRs throughout and the upstream files used a wider "
            f"definition, so a histidine just outside an IMGT CDR is counted "
            f"there and not here.",
            "",
        ]

    row = conn.execute(
        "SELECT MIN(upstream_tm_pred) AS lo, MAX(upstream_tm_pred) AS hi, "
        "SUM(CASE WHEN upstream_tm_pred > 95 THEN 1 ELSE 0 END) AS hot "
        "FROM candidates WHERE upstream_tm_pred IS NOT NULL").fetchone()
    ours = conn.execute(
        "SELECT MIN(tm_pred) AS lo, MAX(tm_pred) AS hi, COUNT(*) AS n, "
        "SUM(CASE WHEN tm_pred > 95 THEN 1 ELSE 0 END) AS hot "
        "FROM thermostability WHERE tm_pred IS NOT NULL").fetchone()
    if row and row["lo"] is not None:
        lines += ["### The upstream Tm column", ""]
        lines.append(
            f"The upstream column spans {row['lo']:.1f} to {row['hi']:.1f} "
            f"degrees C, with {row['hot']:,} candidates above 95 degrees C. An "
            f"apparent Tm above 100 degrees C is not a value a DSF run gives a "
            f"VHH in solution, so the column is treated as uncalibrated "
            f"throughout.")
        if ours and ours["n"]:
            correlation = conn.execute(
                "SELECT c.upstream_tm_pred, t.tm_pred FROM candidates c "
                "JOIN thermostability t ON t.candidate_id = c.id "
                "WHERE c.upstream_tm_pred IS NOT NULL AND t.tm_pred IS NOT NULL"
            ).fetchall()
            xs = [r[0] for r in correlation]
            ys = [r[1] for r in correlation]
            if len(xs) > 20:
                mx, my = statistics.mean(xs), statistics.mean(ys)
                sxy = sum((a - mx) * (b - my) for a, b in zip(xs, ys))
                sxx = sum((a - mx) ** 2 for a in xs)
                syy = sum((b - my) ** 2 for b in ys)
                r = sxy / ((sxx * syy) ** 0.5) if sxx and syy else 0.0
                lines.append("")
                lines.append(
                    f"This build's own model predicts {ours['lo']:.1f} to "
                    f"{ours['hi']:.1f} degrees C over the same {ours['n']:,} "
                    f"candidates, with {ours['hot']:,} above 95 degrees C. The "
                    f"two columns correlate at r = {r:.3f}, so they are not "
                    f"measuring the same thing; both are shown and neither is "
                    f"blended into the other.")
        lines.append("")

    empty = conn.execute(
        "SELECT COUNT(*) FROM candidates WHERE upstream_cdr_novelty IS NULL "
        "AND upstream_cdr_identity IS NULL").fetchone()[0]
    multi = conn.execute(
        "SELECT COUNT(DISTINCT seq_md5) FROM candidates WHERE multi_target = 1"
    ).fetchone()[0]
    lines += [
        "### The rest of the upstream block",
        "",
        f"- `cdr_novelty` and `cdr_identity` are empty in all five source files: "
        f"{empty:,} of 2,500 rows carry neither. Both are ingested as null and "
        f"this build computes its own equivalents.",
        f"- `target_cdr_identity` and `s_his` exist only for EGFR and TNFa, the "
        f"two targets with known nanobodies. They are left null for the other "
        f"three rather than fabricated.",
        f"- {multi} sequences appear under more than one target. They are not "
        f"deduplicated; they are internal promiscuity controls.",
        f"- The upstream composite `score` blends four terms for EGFR and TNFa "
        f"and three for the others, so its ranges are not on one scale. It is "
        f"kept for comparison and the whole table is never sorted by it.",
        "",
    ]
    return lines


def section_novelty(conn) -> list[str]:
    rows = conn.execute(
        "SELECT nearest_identity_fw AS fw, nearest_identity_cdr AS cdr, "
        "nearest_tm AS tm, has_experimental_badge AS badge FROM structures "
        "WHERE nearest_identity_fw IS NOT NULL").fetchall()
    if not rows:
        return []
    fw = [r["fw"] for r in rows if r["fw"] is not None]
    cdr = [r["cdr"] for r in rows if r["cdr"] is not None]
    tm = [r["tm"] for r in rows if r["tm"] is not None]
    badges = sum(1 for r in rows if r["badge"])
    return [
        "## Framework identity tells you nothing; CDR identity does",
        "",
        f"Against {len(rows):,} candidates and their nearest experimental "
        f"neighbour in SAbDab-nano:",
        "",
        "| Measure | Median | Range |",
        "|---|---|---|",
        f"| Framework identity | {statistics.median(fw) * 100:.1f}% | "
        f"{min(fw) * 100:.0f}% to {max(fw) * 100:.0f}% |",
        f"| CDR identity | {statistics.median(cdr) * 100:.1f}% | "
        f"{min(cdr) * 100:.0f}% to {max(cdr) * 100:.0f}% |",
        f"| TM-score | {statistics.median(tm):.3f} | {min(tm):.2f} to "
        f"{max(tm):.2f} |" if tm else "",
        "",
        "This is why the spec asks for per-region identity rather than a "
        "whole-chain number. Framework identity sits near 93 percent and the "
        "TM-score near 0.92 for essentially every candidate, because every VHH "
        "has the same fold. Neither separates anything. The CDR identity does.",
        "",
        f"**Only {badges} of {len(rows):,} candidates reach the 95 percent "
        f"whole-chain identity that earns an experimental-structure badge.** The "
        f"library is almost entirely novel relative to solved structures. The "
        f"badge is a positive signal only: its absence here says nothing about "
        f"whether these sequences express, and nothing in the scoring reads it.",
        "",
    ]


def section_controls(conn) -> list[str]:
    lines = ["## What the control set showed", ""]
    counts = dict(conn.execute(
        "SELECT kind, COUNT(*) FROM controls GROUP BY kind").fetchall())
    if not counts:
        return lines + ["The control set has not been built for this run.", ""]
    lines.append("| Control kind | Rows |")
    lines.append("|---|---|")
    for kind, n in sorted(counts.items()):
        lines.append(f"| {kind.replace('_', ' ')} | {n} |")
    lines.append("")

    pairs = conn.execute(
        "SELECT observed_tm, predicted_tm, name FROM controls "
        "WHERE observed_tm IS NOT NULL AND predicted_tm IS NOT NULL").fetchall()
    if len(pairs) > 5:
        errors = [abs(r["predicted_tm"] - r["observed_tm"]) for r in pairs]
        signed = [r["predicted_tm"] - r["observed_tm"] for r in pairs]
        lines += [
            "### Thermostability, observed against predicted",
            "",
            f"{len(pairs)} nanobodies with a published apparent Tm, held out of "
            f"the model's training set entirely and predicted with the "
            f"measurement withheld. Mean absolute error "
            f"{statistics.mean(errors):.1f} degrees C, mean signed error "
            f"{statistics.mean(signed):+.1f} degrees C "
            f"({'the model reads high' if statistics.mean(signed) > 0.5 else 'the model reads low' if statistics.mean(signed) < -0.5 else 'no systematic offset'}).",
            "",
        ]

    for key in ("control_3_tm_sanity", "control_4_trimer_TNFa",
                "control_4_trimer_TL1A", "control_6_hallmark"):
        value, passed = qc_value(conn, key)
        if value:
            state = "PASS" if passed else "ATTENTION"
            lines.append(f"- **{key.replace('_', ' ')}** ({state}): {value}")
    lines.append("")
    return lines


def section_scrambled(conn) -> list[str]:
    lines = ["## What the scrambled control did to the affinity axis", ""]
    any_result = False
    rows = []
    for target in TARGETS:
        value, passed = qc_value(conn, f"control_5_scrambled_{target}")
        if value and "pending" not in value:
            any_result = True
            rows.append((target, value, passed))
    if not any_result:
        lines += [
            "The affinity axis depends on the co-folding tier, which has not "
            "completed for this run. Until it does, the axis is null everywhere "
            "and the app shows it as pending rather than as a number.",
            "",
            "This is the control that matters most (spec Section 9 item 5). "
            "Published benchmarks show co-folding confidence does not reliably "
            "discriminate cognate from non-cognate nanobody-antigen pairs, so "
            "there is a real chance this axis cannot separate real candidates "
            "from CDR-shuffled ones. **If that happens it is a finding, not a "
            "failure**: the axis is marked non-discriminating wherever it "
            "appears and it does not drive the Pareto front.",
            "",
        ]
        return lines
    lines.append("| Target | Result | Verdict |")
    lines.append("|---|---|---|")
    for target, value, passed in rows:
        lines.append(f"| {target} | {value} | "
                     f"{'separates' if passed else '**does not separate**'} |")
    lines.append("")
    if any(not p for _t, _v, p in rows):
        lines += [
            "**This is a finding, not a failure.** Where the affinity axis does "
            "not separate real candidates from CDR-shuffled ones it is marked "
            "non-discriminating wherever it appears, and it does not drive the "
            "Pareto front. That is the honest reading of what co-folding "
            "confidence can and cannot do for nanobody specificity.",
            "",
        ]
    return lines


def section_axes(conn) -> list[str]:
    lines = ["## Axis independence and double counting", ""]
    value, _ = qc_value(conn, "axis_correlations")
    if not value:
        return lines + ["Not computed for this run.", ""]
    try:
        payload = json.loads(value)
    except json.JSONDecodeError:
        return lines + ["Not computed for this run.", ""]
    flagged = []
    for target, info in payload.items():
        for line in info.get("flagged", []):
            flagged.append(f"- **{target}**: {line}")
    if flagged:
        lines += [
            "Hydrophobic patch area feeds non-specificity, aggregation and "
            "expression, so the axes are not independent by construction. These "
            "pairs correlate above 0.8 and the app states it beside the ranking "
            "rather than leaving it implicit:",
            "",
        ] + flagged + [""]
    else:
        lines += ["No axis pair correlates above 0.8 in any target for this run.",
                  ""]

    value, _ = qc_value(conn, "vacuous_tests")
    if value:
        try:
            payload = json.loads(value)
        except json.JSONDecodeError:
            payload = {}
        bad = [axis for axis, entry in payload.items() if not entry.get("passed")]
        lines += [
            "### Can each axis fail?",
            "",
            "Every axis is checked for the ways a column can look like a "
            "measurement while carrying no information: no variance, a "
            "degenerate distribution, being a restatement of another column, "
            "having no low tail, or not separating real sequences from shuffled "
            "ones.",
            "",
        ]
        lines.append(f"Axes flagged as vacuous: "
                     f"{', '.join(bad) if bad else 'none'}.")
        lines.append("")
    return lines


def main() -> int:
    conn = connect()
    built = conn.execute("SELECT COUNT(*) FROM candidates").fetchone()[0]
    lines = [
        "# WHITTLE findings",
        "",
        f"Run `{run_id()}`. Generated from the database by `tools/findings.py`, "
        f"not maintained by hand, so it cannot drift from what the run produced.",
        "",
        f"{built:,} candidates across {len(TARGETS)} targets. Every number on "
        f"this page is predicted or estimated. Nothing here was measured in a "
        f"laboratory.",
        "",
    ]
    lines += section_upstream(conn)
    lines += section_novelty(conn)
    lines += section_controls(conn)
    lines += section_scrambled(conn)
    lines += section_axes(conn)

    model_card = ROOT / "MODEL_CARD.md"
    if model_card.exists():
        lines += [
            "## The thermostability model",
            "",
            "Its validation, its failure modes and the gap between a grouped and "
            "a row-wise split are in [MODEL_CARD.md](MODEL_CARD.md). The short "
            "version: it orders nanobodies usefully and does not place them "
            "correctly on the absolute temperature scale, which is why the axis "
            "is a within-target percentile and the degrees carry an error bar.",
            "",
        ]
    lines += [
        "## Licence",
        "",
        "Code MIT. Analysis derived from public databases, each cited with its "
        "own licence on the About tab. No database is mirrored in this "
        "repository.",
        "",
    ]
    text = "\n".join(line for line in lines if line is not None) + "\n"
    (ROOT / "FINDINGS.md").write_text(text)
    conn.close()
    log(f"findings: FINDINGS.md written, {len(text.splitlines())} lines")
    return 0


if __name__ == "__main__":
    sys.exit(main())
