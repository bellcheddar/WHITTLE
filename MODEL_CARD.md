# WHITTLE thermostability model, `whittle_tm_v1`

An **estimate of apparent melting temperature**, not a measurement. Apparent Tm depends on the method, the buffer and the protein concentration, and nothing on this site was measured in a laboratory.

Built 2026-11-13 17:05.

## Training set

- 650 unique nanobody sequences with a measured apparent Tm, spanning 26.64 to 98.23 degrees C.
- 377 clusters at 90 percent CDR3 identity. The cross-validation splits on these clusters, never on rows.
- Sources:
  - **nanomelt**: Supplementary Dataset 1, 640 curated apparent-Tm nanobodies (doi:10.1080/19420862.2024.2442750, CC BY, as the article).
  - **nbthermo**: NbThermo, 548 entries; AHo sequences with gaps stripped (doi:10.1093/database/baad021, the repository carries no licence file; cited to the CC-BY article, derived values only).

### How good are the labels

486 sequences carry more than one measurement, and 135 of those disagree between measurements, by a mean of 3.63 and at most 30.5 degrees C. That disagreement is the floor on any honest error bar, and the calibration below is floored at it.

## Features

One base Gaussian Process and one Support Vector regressor per block, stacked under a Ridge meta-learner (`sklearn.ensemble.StackingRegressor`).

| Block | Dimensions | What it is |
|---|---|---|
| one-hot | 106 | amino-acid composition of the whole chain and of each IMGT region, plus lengths |
| VHSE | 80 | the eight VHSE descriptors (Mei et al. 2005), mean and standard deviation per region |
| ESM-2 | 1280 | mean-pooled `esm2_t33_650M_UR50D` embeddings, float32 |

## Validation

| Split | n | R squared | RMSE, deg C | MAE, deg C | Spearman |
|---|---|---|---|---|---|
| **Grouped by 90 percent CDR3 identity** | 650 | **0.390** | 8.40 | 6.21 | 0.646 |
| Row-wise, for contrast | 650 | 0.697 | 5.92 | 4.36 | 0.799 |
| **External held-out set** | 77 | **-0.058** | 9.16 | 6.91 | 0.673 |

The grouped figure is the one to quote. The row-wise figure is +0.307 R squared higher, and that gap is the near-duplicate leak a row-wise split hides: it is reported here rather than used.

## Error estimate

Per row, from the cosine distance to the nearest training neighbour in the ESM-2 block. Linear fit of the median absolute residual against the cosine distance to the nearest training neighbour, with a floor at the 2.4 degree c mean disagreement between repeat measurements reported for this label set. The quoted interval covers 76.0 percent of the out-of-fold residuals and 66 percent of the external-set residuals
.

A candidate far from the training manifold gets a wide bar, not a crisp number, and the app draws the bar rather than only the number.

## Known failure modes

- **Away from the training manifold the model regresses towards the training mean.** The error bar widens with distance, but a wide bar on a central prediction is still a weak statement.
- **Apparent Tm is method-dependent.** The labels mix nanoDSF, DSF, circular dichroism and DSC. Differences of a few degrees between methods for the same nanobody are ordinary.
- **The training set is small.** A few hundred clusters cannot cover the sequence space of 2,500 designed candidates.
- **Long or stapled CDR3s are under-represented** in the labels relative to the candidate set.
- **The model predicts apparent Tm, not stability in any particular formulation**, and says nothing about aggregation on storage.

## Positive control

Spec Section 9 item 3 asserts that no predicted apparent Tm in this model's output exceeds 95 degrees C. If one appears, the model is extrapolating and the error bar must reflect it. The check runs in `tools/controls.py` and its result appears on the QC panel.

## Licence and citation

Model code MIT, with the rest of WHITTLE. scikit-learn is BSD-3 and ESM-2 is MIT. The labels are cited to their CC-BY articles above; no label set is redistributed in this repository.
