# WHITTLE findings

Run `whittle-2026-11`. Generated from the database by `tools/findings.py`, not maintained by hand, so it cannot drift from what the run produced.

2,500 candidates across 5 targets. Every number on this page is predicted or estimated. Nothing here was measured in a laboratory.

## Where this build disagrees with the upstream columns

### The camelid hallmark tetrad

The upstream `hallmark` column calls 92.8 percent of rows the same way a strict all-four reading of IMGT 42, 49, 50 and 52 does, and 99.3 percent the same way an at-least-three-of-four reading does. The upstream column therefore implements the three-of-four rule, not the four-of-four one it appears to.

| Camelid residues present | Candidates | Share |
|---|---|---|
| 4 of 4 | 2,321 | 92.8% |
| 3 of 4 | 165 | 6.6% |
| 2 of 4 | 14 | 0.6% |

93 candidates carry a VH-type residue at one of the four positions. That leaves the former light-chain interface hydrophobic and exposed, which is a solubility and stability liability rather than a cosmetic difference. This build reports the 0 to 4 count and the strict boolean side by side.

### CDR histidine counts

Recomputed IMGT CDR histidine totals agree with the upstream columns for 90.9 percent of rows. The residual is a CDR-definition difference, not a counting error: this build uses IMGT CDRs throughout and the upstream files used a wider definition, so a histidine just outside an IMGT CDR is counted there and not here.

### The upstream Tm column

The upstream column spans 71.2 to 101.1 degrees C, with 38 candidates above 95 degrees C. An apparent Tm above 100 degrees C is not a value a DSF run gives a VHH in solution, so the column is treated as uncalibrated throughout.

This build's own model predicts 58.8 to 91.1 degrees C over the same 2,500 candidates, with 0 above 95 degrees C. The two columns correlate at r = 0.135, so they are not measuring the same thing; both are shown and neither is blended into the other.

### The rest of the upstream block

- `cdr_novelty` and `cdr_identity` are empty in all five source files: 2,500 of 2,500 rows carry neither. Both are ingested as null and this build computes its own equivalents.
- `target_cdr_identity` and `s_his` exist only for EGFR and TNFa, the two targets with known nanobodies. They are left null for the other three rather than fabricated.
- 10 sequences appear under more than one target. They are not deduplicated; they are internal promiscuity controls.
- The upstream composite `score` blends four terms for EGFR and TNFa and three for the others, so its ranges are not on one scale. It is kept for comparison and the whole table is never sorted by it.

## Framework identity tells you nothing; CDR identity does

Against 2,500 candidates and their nearest experimental neighbour in SAbDab-nano:

| Measure | Median | Range |
|---|---|---|
| Framework identity | 93.4% | 74% to 100% |
| CDR identity | 60.0% | 17% to 100% |
| TM-score | 0.918 | 0.70 to 0.98 |

This is why the spec asks for per-region identity rather than a whole-chain number. Framework identity sits near 93 percent and the TM-score near 0.92 for essentially every candidate, because every VHH has the same fold. Neither separates anything. The CDR identity does.

**Only 2 of 2,500 candidates reach the 95 percent whole-chain identity that earns an experimental-structure badge.** The library is almost entirely novel relative to solved structures. The badge is a positive signal only: its absence here says nothing about whether these sequences express, and nothing in the scoring reads it.

## What the control set showed

| Control kind | Rows |
|---|---|
| known | 94 |
| multi target | 23 |
| scrambled | 1000 |

### Does each axis separate a candidate from its own shuffled twin?

Each scrambled control keeps its parent's amino-acid composition,
length and framework exactly, so CDR order is the only difference
left. The paired sign test below compares each control with the one
candidate it was shuffled from. The verdict is pooled over all five
targets and the per-target counts travel with it.

| Axis | Pooled | Rate | p | Verdict | Per target |
|---|---|---|---|---|---|
| affinity | pending | | | pending | |
| specificity | 679 of 839 | 81% | 9e-77 | informative | EGFR 130/169, GLP1R 141/162, PRLR 140/171, TL1A 137/169, TNFa 131/168 |
| thermo | 837 of 1000 | 84% | 1e-109 | informative | EGFR 164/200, GLP1R 162/200, PRLR 174/200, TL1A 167/200, TNFa 170/200 |
| expression | 982 of 1000 | 98% | 3e-263 | untestable | EGFR 195/200, GLP1R 199/200, PRLR 199/200, TL1A 195/200, TNFa 194/200 |

The expression row is the one to read carefully. Its 98 per cent is
not evidence that the axis works: two thirds of its weight is
amino-acid composition, which the shuffle preserves exactly, so the
comparison is largely a sequence against itself and a high rate is
what that looks like. It is reported untestable for that reason,
not ranked above thermostability.

An earlier build reported the specificity axis at 60 per cent and
weak. That was withdrawn: the control set had been scored before
the non-specific-binding predictor had run on it, so the heaviest
component of the axis was missing from the control side and present
on the candidate side, and a third of the axis was being ranked
against all of it. `pipeline.score.score_external` now refuses to
place an external row unless it used exactly the components the
candidates were scored on. See DECISIONS 59.

### Thermostability, observed against predicted

77 nanobodies with a published apparent Tm, held out of the model's training set entirely and predicted with the measurement withheld. Mean absolute error 7.0 degrees C, mean signed error +3.5 degrees C (the model reads high).

- **control 3 tm sanity** (PASS): 0 of 2500 predictions above 95 deg C; hottest 91.1 | Spec 9 item 3. A prediction above 95 deg C means the model is extrapolating and its error bar must reflect that.
- **control 4 trimer TNFa** (PASS): 311 of 340 co-folds contact more than one subunit | Spec 9 item 4. If every interface sits on a single chain, suspect the trimer was never assembled.
- **control 4 trimer TL1A** (PASS): 304 of 340 co-folds contact more than one subunit | Spec 9 item 4. If every interface sits on a single chain, suspect the trimer was never assembled.
- **control 6 hallmark** (PASS): 0.993 agreement under the three-of-four rule the upstream column implements | Spec 9 item 6. The strict all-four call agrees for 0.928 of rows; both are reported.

## What the scrambled control did to the affinity axis

The affinity axis depends on the co-folding tier, which has not completed for this run. Until it does, the axis is null everywhere and the app shows it as pending rather than as a number.

This is the control that matters most (spec Section 9 item 5). Published benchmarks show co-folding confidence does not reliably discriminate cognate from non-cognate nanobody-antigen pairs, so there is a real chance this axis cannot separate real candidates from CDR-shuffled ones. **If that happens it is a finding, not a failure**: the axis is marked non-discriminating wherever it appears and it does not drive the Pareto front.

## Axis independence and double counting

No axis pair correlates above 0.8 in any target for this run.

### Can each axis fail?

Every axis is checked for the ways a column can look like a measurement while carrying no information: no variance, a degenerate distribution, being a restatement of another column, having no low tail, or not separating real sequences from shuffled ones.

Axes flagged as vacuous: affinity, expression.

## The thermostability model

Its validation, its failure modes and the gap between a grouped and a row-wise split are in [MODEL_CARD.md](MODEL_CARD.md). The short version: it orders nanobodies usefully and does not place them correctly on the absolute temperature scale, which is why the axis is a within-target percentile and the degrees carry an error bar.

## Licence

Code MIT. Analysis derived from public databases, each cited with its own licence on the About tab. No database is mirrored in this repository.

