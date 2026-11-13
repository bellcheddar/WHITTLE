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
| multi target | 10 |
| scrambled | 250 |

- **control 3 tm sanity** (ATTENTION): no Tm predictions yet | Spec 9 item 3, pending the thermostability model.
- **control 4 trimer TNFa** (ATTENTION): no co-folds yet | Spec 9 item 4, pending the co-folding tier.
- **control 4 trimer TL1A** (ATTENTION): no co-folds yet | Spec 9 item 4, pending the co-folding tier.
- **control 6 hallmark** (PASS): 0.993 agreement under the three-of-four rule the upstream column implements | Spec 9 item 6. The strict all-four call agrees for 0.928 of rows; both are reported.

## What the scrambled control did to the affinity axis

The affinity axis depends on the co-folding tier, which has not completed for this run. Until it does, the axis is null everywhere and the app shows it as pending rather than as a number.

This is the control that matters most (spec Section 9 item 5). Published benchmarks show co-folding confidence does not reliably discriminate cognate from non-cognate nanobody-antigen pairs, so there is a real chance this axis cannot separate real candidates from CDR-shuffled ones. **If that happens it is a finding, not a failure**: the axis is marked non-discriminating wherever it appears and it does not drive the Pareto front.

## Axis independence and double counting

No axis pair correlates above 0.8 in any target for this run.

### Can each axis fail?

Every axis is checked for the ways a column can look like a measurement while carrying no information: no variance, a degenerate distribution, being a restatement of another column, having no low tail, or not separating real sequences from shuffled ones.

Axes flagged as vacuous: none.

## The thermostability model

Its validation, its failure modes and the gap between a grouped and a row-wise split are in [MODEL_CARD.md](MODEL_CARD.md). The short version: it orders nanobodies usefully and does not place them correctly on the absolute temperature scale, which is why the axis is a within-target percentile and the degrees carry an error bar.

## Licence

Code MIT. Analysis derived from public databases, each cited with its own licence on the About tab. No database is mirrored in this repository.

