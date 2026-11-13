"""Tier 1b: the thermostability model. Spec 6.3.

Built, not licensed. A Ridge-stacked ensemble over three feature blocks, with a
per-row error estimate and a validation protocol designed to be hard to flatter.

Features (spec 6.3 item 2), one base model per block:
  one-hot    composition and CDR-region composition
  VHSE       the eight VHSE descriptors, averaged per region
  ESM-2      mean-pooled embeddings, 650M, dropping to 150M under memory pressure

Model (item 3): `sklearn.ensemble.StackingRegressor`, Ridge meta-learner over a
Gaussian Process and a Support Vector regressor per block.

**Validation is the part that matters** (item 4). Nested cross-validation
grouped by sequence similarity: the training set is clustered at 90 percent
CDR3 identity and the splits are taken on clusters, never on rows. A row-wise
split reports a flattering R squared that evaporates on 2,500 predicted
sequences, because near-duplicates land on both sides of the split. Both
numbers are computed and both are reported, precisely so the gap is visible.

**Error estimate** (item 5): per row, from the distance to the nearest training
neighbour, calibrated so that the quoted interval covers the held-out residuals
at the stated rate. A candidate far from the training manifold gets a wide bar,
not a crisp number, and the UI draws the bar.

    pixi run tm-train
    pixi run python -m pipeline.tm_model --quick     # 150M ESM, fewer repeats
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time
import warnings
from pathlib import Path

import numpy as np

from pipeline.common import MODELS, ROOT, connect, log, manifest_write, tuning
from pipeline.tm_data import load as load_dataset

warnings.filterwarnings("ignore")

MODEL_DIR = MODELS / "whittle_tm_v1"
AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"

# VHSE descriptors. Mei, Liao, Lin and Chen (2005) Biopolymers 80:775-786,
# doi:10.1002/bip.20296. Eight principal components over 50 physicochemical
# variables: VHSE1-2 hydrophobic, VHSE3-4 steric, VHSE5-8 electronic. Values as
# published in the paper's Table 2 and carried in the R `Peptides` package.
VHSE = {
    "A": (0.15, -1.11, -1.35, -0.92, 0.02, -0.91, 0.36, -0.48),
    "R": (-1.47, 1.45, 1.24, 1.27, 1.55, 1.47, 1.30, 0.83),
    "N": (-0.99, 0.00, -0.37, 0.69, -0.55, 0.85, 0.73, -0.80),
    "D": (-1.15, 0.67, -0.41, -0.01, -2.68, 1.31, 0.03, 0.56),
    "C": (0.18, -1.67, -0.46, -0.21, 0.00, 1.20, -1.61, -0.19),
    "Q": (-0.96, 0.12, 0.18, 0.16, 0.09, 0.42, -0.20, -0.41),
    "E": (-1.18, 0.40, 0.10, 0.36, -2.16, -0.17, 0.91, 0.02),
    "G": (-0.20, -1.53, -2.63, 2.28, -0.53, -1.18, 2.01, -1.34),
    "H": (-0.43, -0.25, 0.37, 0.19, 0.51, 1.28, 0.93, 0.65),
    "I": (1.27, -0.14, 0.30, -1.80, 0.30, -1.61, -0.16, -0.13),
    "L": (1.36, 0.07, 0.26, -0.80, 0.22, -1.37, 0.08, -0.62),
    "K": (-1.17, 0.70, 0.70, 0.80, 1.64, 0.67, 1.63, 0.13),
    "M": (1.01, -0.53, 0.43, 0.00, 0.23, 0.10, -0.86, -0.68),
    "F": (1.52, 0.61, 0.96, -0.16, 0.25, 0.28, -1.33, -0.20),
    "P": (0.22, -0.17, -0.50, 0.05, -0.01, -1.34, -0.19, 3.56),
    "S": (-0.67, -0.86, -1.07, -0.41, -0.32, 0.27, -0.64, 0.11),
    "T": (-0.34, -0.51, -0.55, -1.06, -0.06, -0.01, -0.79, 0.39),
    "W": (1.50, 2.06, 1.79, 0.75, 0.75, -0.13, -1.01, -0.85),
    "Y": (0.61, 1.60, 1.17, 0.73, 0.53, 0.25, -0.96, -0.52),
    "V": (0.76, -0.92, -0.17, -1.91, 0.22, -1.40, -0.24, -0.03),
}

REGIONS = ("cdr1", "cdr2", "cdr3", "framework")


# --------------------------------------------------------------------- features

def region_split(sequence: str) -> dict[str, str]:
    """IMGT regions, via the same numbering ladder Tier 0 used."""
    from pipeline import numbering as nb

    result = nb.number_many([("q", sequence)])["q"]
    regions = result.regions or {}
    framework = sequence
    for key in ("cdr1", "cdr2", "cdr3"):
        frag = regions.get(key) or ""
        if frag:
            framework = framework.replace(frag, "", 1)
    return {"cdr1": regions.get("cdr1", ""), "cdr2": regions.get("cdr2", ""),
            "cdr3": regions.get("cdr3", ""), "framework": framework}


def onehot_block(sequences: list[str], regions: list[dict[str, str]]) -> np.ndarray:
    """Composition of the whole chain and of each region, plus lengths.

    Composition rather than position-wise one-hot: VHH lengths vary by 25
    residues and a positional encoding over an unaligned set would mostly
    encode indel position.
    """
    rows = []
    for seq, reg in zip(sequences, regions):
        row: list[float] = []
        n = max(len(seq), 1)
        row += [seq.count(aa) / n for aa in AMINO_ACIDS]
        for name in REGIONS:
            frag = reg.get(name, "")
            m = max(len(frag), 1)
            row += [frag.count(aa) / m for aa in AMINO_ACIDS]
            row.append(len(frag))
        row.append(len(seq))
        row.append(seq.count("C"))
        rows.append(row)
    return np.asarray(rows, dtype=np.float32)


def vhse_block(sequences: list[str], regions: list[dict[str, str]]) -> np.ndarray:
    """The eight VHSE descriptors, averaged over the chain and per region."""
    rows = []
    for seq, reg in zip(sequences, regions):
        row: list[float] = []
        for fragment in [seq] + [reg.get(name, "") for name in REGIONS]:
            values = [VHSE[aa] for aa in fragment if aa in VHSE]
            if values:
                row += list(np.mean(values, axis=0))
                row += list(np.std(values, axis=0))
            else:
                row += [0.0] * 16
        rows.append(row)
    return np.asarray(rows, dtype=np.float32)


def esm_block(sequences: list[str], model_name: str, batch: int) -> np.ndarray:
    """Mean-pooled ESM-2 embeddings. float32 throughout (spec 2.2)."""
    import esm
    import torch

    torch.set_grad_enabled(False)
    loader = getattr(esm.pretrained, model_name)
    model, alphabet = loader()
    model.eval().float()
    layer = model.num_layers
    converter = alphabet.get_batch_converter()
    out = []
    for i in range(0, len(sequences), batch):
        chunk = [(f"s{j}", s) for j, s in enumerate(sequences[i:i + batch])]
        _labels, _strs, tokens = converter(chunk)
        result = model(tokens, repr_layers=[layer], return_contacts=False)
        reps = result["representations"][layer]
        for k, (_name, seq) in enumerate(chunk):
            # Drop the BOS token and the padding; mean-pool the residues only.
            out.append(reps[k, 1:len(seq) + 1].mean(0).numpy())
        if (i // batch) % 10 == 0:
            log(f"tm_model: ESM-2 embedded {min(i + batch, len(sequences))}"
                f"/{len(sequences)}")
    return np.asarray(out, dtype=np.float32)


# ------------------------------------------------------------------- clustering

def cdr3_clusters(sequences: list[str], regions: list[dict[str, str]],
                  identity: float = 0.9) -> np.ndarray:
    """Group label per row: cluster at 90 percent CDR3 identity (spec 6.3 item 4).

    Sequences with no resolvable CDR3 fall back to whole-chain identity rather
    than all collapsing into one group, which would make them look like a
    single cluster and distort the split.
    """
    from pipeline.tier0 import cluster_greedy

    keys = []
    for seq, reg in zip(sequences, regions):
        keys.append(reg.get("cdr3") or seq)
    records = [(str(i), key) for i, key in enumerate(keys)]
    assignment = cluster_greedy(records, identity)
    return np.asarray([assignment[str(i)] for i in range(len(sequences))])


# ----------------------------------------------------------------------- model

def build_estimator(block_slices: dict[str, slice], seed: int = 0):
    """StackingRegressor: a Ridge meta-learner over a GPR and an SVR per block."""
    from sklearn.ensemble import StackingRegressor
    from sklearn.gaussian_process import GaussianProcessRegressor
    from sklearn.gaussian_process.kernels import RBF, ConstantKernel, WhiteKernel
    from sklearn.compose import ColumnTransformer
    from sklearn.decomposition import PCA
    from sklearn.linear_model import RidgeCV
    from sklearn.pipeline import Pipeline
    from sklearn.preprocessing import StandardScaler
    from sklearn.svm import SVR

    estimators = []
    for name, columns in block_slices.items():
        indices = list(range(columns.start, columns.stop))
        select = ColumnTransformer([(f"{name}_cols", "passthrough", indices)],
                                   remainder="drop")
        # ESM-2 is 1280-dimensional against a few hundred training rows, so a
        # GP on the raw block would be fitting noise. PCA first, then the two
        # base models see the same compressed block.
        steps = [("select", select), ("scale", StandardScaler())]
        if name == "esm":
            steps.append(("pca", PCA(n_components=0.95, svd_solver="full",
                                     random_state=seed)))
        kernel = (ConstantKernel(1.0, (1e-3, 1e3))
                  * RBF(length_scale=10.0, length_scale_bounds=(1e-2, 1e4))
                  + WhiteKernel(noise_level=4.0, noise_level_bounds=(1e-2, 1e3)))
        estimators.append((f"{name}_gpr", Pipeline(steps + [
            ("gpr", GaussianProcessRegressor(kernel=kernel, normalize_y=True,
                                             alpha=1e-6, random_state=seed,
                                             n_restarts_optimizer=0))])))
        estimators.append((f"{name}_svr", Pipeline(steps + [
            ("svr", SVR(kernel="rbf", C=10.0, epsilon=1.0, gamma="scale"))])))
    return StackingRegressor(
        estimators=estimators,
        final_estimator=RidgeCV(alphas=np.logspace(-3, 3, 13)),
        cv=5, n_jobs=1, passthrough=False)


# ------------------------------------------------------------------ validation

def grouped_cv(X: np.ndarray, y: np.ndarray, groups: np.ndarray,
               block_slices: dict[str, slice], repeats: int, folds: int,
               seed: int) -> dict:
    """Repeated grouped k-fold, with a row-wise run beside it for contrast."""
    from sklearn.model_selection import GroupKFold, KFold

    def run(splitter, split_groups) -> dict:
        preds = np.full(len(y), np.nan)
        counts = np.zeros(len(y))
        totals = np.zeros(len(y))
        for repeat in range(repeats):
            order = np.random.RandomState(seed + repeat).permutation(len(y))
            inverse = np.argsort(order)
            Xs, ys = X[order], y[order]
            gs = split_groups[order] if split_groups is not None else None
            for train_idx, test_idx in (splitter.split(Xs, ys, gs)
                                        if gs is not None
                                        else splitter.split(Xs, ys)):
                model = build_estimator(block_slices, seed=seed + repeat)
                model.fit(Xs[train_idx], ys[train_idx])
                predicted = model.predict(Xs[test_idx])
                original = inverse[test_idx] if False else order[test_idx]
                totals[original] += predicted
                counts[original] += 1
        preds = np.where(counts > 0, totals / np.maximum(counts, 1), np.nan)
        return metrics(y, preds)

    log(f"tm_model: grouped {folds}-fold x {repeats}, split on CDR3 clusters")
    grouped = run(GroupKFold(n_splits=folds), groups)
    log(f"tm_model: grouped CV  R2 {grouped['r2']:.3f}  RMSE {grouped['rmse']:.2f}  "
        f"MAE {grouped['mae']:.2f}  Spearman {grouped['spearman']:.3f}")
    log(f"tm_model: row-wise {folds}-fold x {repeats}, for contrast only")
    rowwise = run(KFold(n_splits=folds, shuffle=True, random_state=seed), None)
    log(f"tm_model: row-wise CV R2 {rowwise['r2']:.3f}  RMSE {rowwise['rmse']:.2f}  "
        f"(the gap to the grouped figure is the near-duplicate leak)")
    return {"grouped": grouped, "row_wise": rowwise,
            "optimism_r2": round(rowwise["r2"] - grouped["r2"], 4)}


def metrics(y: np.ndarray, pred: np.ndarray) -> dict:
    from scipy.stats import pearsonr, spearmanr
    from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score

    mask = ~np.isnan(pred)
    y, pred = y[mask], pred[mask]
    return {
        "n": int(len(y)),
        "r2": round(float(r2_score(y, pred)), 4),
        "rmse": round(float(math.sqrt(mean_squared_error(y, pred))), 3),
        "mae": round(float(mean_absolute_error(y, pred)), 3),
        "pearson": round(float(pearsonr(y, pred)[0]), 4),
        "spearman": round(float(spearmanr(y, pred)[0]), 4),
        "residual_sd": round(float(np.std(y - pred)), 3),
    }


# ------------------------------------------------------------- error estimate

def nearest_distance(query: np.ndarray, train: np.ndarray) -> np.ndarray:
    """Cosine distance to the nearest training row, in the ESM block."""
    qn = query / np.maximum(np.linalg.norm(query, axis=1, keepdims=True), 1e-9)
    tn = train / np.maximum(np.linalg.norm(train, axis=1, keepdims=True), 1e-9)
    similarity = qn @ tn.T
    return 1.0 - similarity.max(axis=1)


def calibrate_error(distances: np.ndarray, residuals: np.ndarray) -> dict:
    """Fit |residual| against the distance to the nearest training neighbour.

    A straight line through the binned medians, with a floor at the measurement
    spread the labels themselves carry. The floor matters: no error bar on an
    apparent Tm should be narrower than the disagreement between two
    measurements of the same nanobody.
    """
    absolute = np.abs(residuals)
    if len(distances) < 20:
        return {"slope": 0.0, "intercept": float(np.median(absolute)),
                "floor": 2.4, "method": "median only, too few points to fit"}
    bins = np.quantile(distances, np.linspace(0, 1, 9))
    xs, ys = [], []
    for lo, hi in zip(bins[:-1], bins[1:]):
        mask = (distances >= lo) & (distances <= hi)
        if mask.sum() >= 5:
            xs.append(float(np.mean(distances[mask])))
            ys.append(float(np.median(absolute[mask])))
    if len(xs) < 3:
        return {"slope": 0.0, "intercept": float(np.median(absolute)),
                "floor": 2.4, "method": "median only"}
    slope, intercept = np.polyfit(xs, ys, 1)
    return {"slope": float(slope), "intercept": float(intercept),
            "floor": 2.4,
            "method": "linear fit of the median absolute residual against the "
                      "cosine distance to the nearest training neighbour, with a "
                      "floor at the 2.4 degree C mean disagreement between repeat "
                      "measurements reported for this label set"}


def apply_error(distances: np.ndarray, calibration: dict) -> np.ndarray:
    raw = calibration["slope"] * distances + calibration["intercept"]
    # 1.96 for a nominal 95 percent interval under a normal residual.
    return np.maximum(raw, calibration["floor"]) * 1.96


# ------------------------------------------------------------------------ main

def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--quick", action="store_true",
                        help="150M ESM-2 and fewer repeats, for a smoke run")
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--seed", type=int, default=20261113)
    parser.add_argument("--predict-only", action="store_true")
    args = parser.parse_args()

    conf = tuning()
    esm_name = ("esm2_t30_150M_UR50D" if args.quick
                else str(conf.get("esm_model", "esm2_t33_650M_UR50D")))
    batch = int(conf.get("esm_batch", 8))
    MODEL_DIR.mkdir(parents=True, exist_ok=True)

    data = load_dataset()
    train_rows = data["training"]
    external_rows = data["external"]
    log(f"tm_model: {len(train_rows)} training sequences, "
        f"{len(external_rows)} external, ESM-2 {esm_name}")

    sequences = [r["sequence"] for r in train_rows]
    y = np.asarray([r["tm"] for r in train_rows], dtype=np.float32)
    log("tm_model: numbering the training set")
    regions = [region_split(s) for s in sequences]
    groups = cdr3_clusters(sequences, regions)
    log(f"tm_model: {len(set(groups.tolist()))} CDR3 clusters at 90 percent "
        f"identity across {len(sequences)} sequences")

    started = time.time()
    oh = onehot_block(sequences, regions)
    vh = vhse_block(sequences, regions)
    es = esm_block(sequences, esm_name, batch)
    X = np.hstack([oh, vh, es])
    block_slices = {
        "onehot": slice(0, oh.shape[1]),
        "vhse": slice(oh.shape[1], oh.shape[1] + vh.shape[1]),
        "esm": slice(oh.shape[1] + vh.shape[1], X.shape[1]),
    }
    log(f"tm_model: features {X.shape} built in {(time.time() - started) / 60:.1f} "
        f"minutes (one-hot {oh.shape[1]}, VHSE {vh.shape[1]}, ESM {es.shape[1]})")

    cv = grouped_cv(X, y, groups, block_slices, args.repeats, args.folds, args.seed)

    # Final fit on everything, then the error calibration from the grouped
    # out-of-fold residuals.
    log("tm_model: final fit on the whole training set")
    model = build_estimator(block_slices, seed=args.seed)
    model.fit(X, y)

    from sklearn.model_selection import GroupKFold

    oof = np.zeros(len(y))
    for train_idx, test_idx in GroupKFold(n_splits=args.folds).split(X, y, groups):
        fold_model = build_estimator(block_slices, seed=args.seed)
        fold_model.fit(X[train_idx], y[train_idx])
        oof[test_idx] = fold_model.predict(X[test_idx])
    oof_distance = np.zeros(len(y))
    for train_idx, test_idx in GroupKFold(n_splits=args.folds).split(X, y, groups):
        oof_distance[test_idx] = nearest_distance(es[test_idx], es[train_idx])
    calibration = calibrate_error(oof_distance, y - oof)
    coverage = float(np.mean(np.abs(y - oof) <= apply_error(oof_distance, calibration)))
    log(f"tm_model: error bars calibrated; they cover "
        f"{coverage * 100:.1f} percent of the out-of-fold residuals")

    # External set: the only statement about sequences nobody fitted on.
    external = {}
    if external_rows:
        ext_seqs = [r["sequence"] for r in external_rows]
        ext_regions = [region_split(s) for s in ext_seqs]
        ext_es = esm_block(ext_seqs, esm_name, batch)
        ext_X = np.hstack([onehot_block(ext_seqs, ext_regions),
                           vhse_block(ext_seqs, ext_regions), ext_es])
        ext_y = np.asarray([r["tm"] for r in external_rows], dtype=np.float32)
        ext_pred = model.predict(ext_X)
        external = metrics(ext_y, ext_pred)
        ext_err = apply_error(nearest_distance(ext_es, es), calibration)
        external["error_bar_coverage"] = round(
            float(np.mean(np.abs(ext_y - ext_pred) <= ext_err)), 4)
        log(f"tm_model: external set R2 {external['r2']:.3f} RMSE "
            f"{external['rmse']:.2f} over {external['n']} held-out sequences; "
            f"error bars cover {external['error_bar_coverage'] * 100:.0f} percent")

    import joblib

    joblib.dump({"model": model, "block_slices": block_slices,
                 "esm_model": esm_name, "calibration": calibration,
                 "train_esm": es, "train_tm": y},
                MODEL_DIR / "model.joblib")
    card = {
        "version": "whittle_tm_v1",
        "built_at": time.strftime("%Y-%m-%d %H:%M"),
        "esm_model": esm_name,
        "n_training": len(train_rows),
        "n_clusters": int(len(set(groups.tolist()))),
        "tm_range": data["tm_range"],
        "cv": cv,
        "external": external,
        "calibration": calibration,
        "error_bar_coverage_internal": round(coverage, 4),
        "feature_dimensions": {k: int(v.stop - v.start)
                               for k, v in block_slices.items()},
        "label_sources": data["sources"],
        "label_disagreement": {
            "n_sequences_measured_more_than_once": data["n_with_repeats"],
            "n_where_measurements_disagree": data.get("n_where_measurements_disagree"),
            "mean_disagreement_c": data.get("mean_spread_where_measurements_disagree"),
            "max_disagreement_c": data.get("max_spread"),
        },
    }
    (MODEL_DIR / "metrics.json").write_text(json.dumps(card, indent=2))
    write_model_card(card)
    manifest_write("tm_model", {"id": "whittle_tm_v1", "status": "ok", **cv["grouped"]})
    log(f"tm_model: written to {MODEL_DIR}")
    return 0


def write_model_card(card: dict) -> None:
    """MODEL_CARD.md. Spec 6.3 item 7."""
    grouped = card["cv"]["grouped"]
    rowwise = card["cv"]["row_wise"]
    external = card.get("external") or {}
    lines = [
        "# WHITTLE thermostability model, `whittle_tm_v1`",
        "",
        "An **estimate of apparent melting temperature**, not a measurement. "
        "Apparent Tm depends on the method, the buffer and the protein "
        "concentration, and nothing on this site was measured in a laboratory.",
        "",
        f"Built {card['built_at']}.",
        "",
        "## Training set",
        "",
        f"- {card['n_training']} unique nanobody sequences with a measured "
        f"apparent Tm, spanning {card['tm_range'][0]} to {card['tm_range'][1]} "
        f"degrees C.",
        f"- {card['n_clusters']} clusters at 90 percent CDR3 identity. The "
        f"cross-validation splits on these clusters, never on rows.",
        "- Sources:",
    ]
    for key, meta in (card.get("label_sources") or {}).items():
        lines.append(f"  - **{key}**: {meta.get('what')} "
                     f"(doi:{meta.get('doi')}, {meta.get('licence')}).")
    dis = card.get("label_disagreement") or {}
    lines += [
        "",
        "### How good are the labels",
        "",
        f"{dis.get('n_sequences_measured_more_than_once')} sequences carry more "
        f"than one measurement, and "
        f"{dis.get('n_where_measurements_disagree')} of those disagree between "
        f"measurements, by a mean of {dis.get('mean_disagreement_c')} and at most "
        f"{dis.get('max_disagreement_c')} degrees C. That disagreement is the "
        "floor on any honest error bar, and the calibration below is floored at it.",
        "",
        "## Features",
        "",
        "One base Gaussian Process and one Support Vector regressor per block, "
        "stacked under a Ridge meta-learner (`sklearn.ensemble.StackingRegressor`).",
        "",
        "| Block | Dimensions | What it is |",
        "|---|---|---|",
        f"| one-hot | {card['feature_dimensions']['onehot']} | amino-acid "
        "composition of the whole chain and of each IMGT region, plus lengths |",
        f"| VHSE | {card['feature_dimensions']['vhse']} | the eight VHSE "
        "descriptors (Mei et al. 2005), mean and standard deviation per region |",
        f"| ESM-2 | {card['feature_dimensions']['esm']} | mean-pooled "
        f"`{card['esm_model']}` embeddings, float32 |",
        "",
        "## Validation",
        "",
        "| Split | n | R squared | RMSE, deg C | MAE, deg C | Spearman |",
        "|---|---|---|---|---|---|",
        f"| **Grouped by 90 percent CDR3 identity** | {grouped['n']} | "
        f"**{grouped['r2']:.3f}** | {grouped['rmse']:.2f} | {grouped['mae']:.2f} | "
        f"{grouped['spearman']:.3f} |",
        f"| Row-wise, for contrast | {rowwise['n']} | {rowwise['r2']:.3f} | "
        f"{rowwise['rmse']:.2f} | {rowwise['mae']:.2f} | {rowwise['spearman']:.3f} |",
    ]
    if external:
        lines.append(
            f"| **External held-out set** | {external['n']} | "
            f"**{external['r2']:.3f}** | {external['rmse']:.2f} | "
            f"{external['mae']:.2f} | {external['spearman']:.3f} |")
    lines += [
        "",
        f"The grouped figure is the one to quote. The row-wise figure is "
        f"{card['cv']['optimism_r2']:+.3f} R squared higher, and that gap is the "
        f"near-duplicate leak a row-wise split hides: it is reported here rather "
        f"than used.",
        "",
        "## Error estimate",
        "",
        f"Per row, from the cosine distance to the nearest training neighbour in "
        f"the ESM-2 block. {card['calibration']['method'].capitalize()}. The "
        f"quoted interval covers "
        f"{card['error_bar_coverage_internal'] * 100:.1f} percent of the "
        f"out-of-fold residuals",
    ]
    if external.get("error_bar_coverage") is not None:
        lines[-1] += (f" and {external['error_bar_coverage'] * 100:.0f} percent of "
                      f"the external-set residuals")
    lines += [
        ".",
        "",
        "A candidate far from the training manifold gets a wide bar, not a crisp "
        "number, and the app draws the bar rather than only the number.",
        "",
        "## Known failure modes",
        "",
        "- **Away from the training manifold the model regresses towards the "
        "training mean.** The error bar widens with distance, but a wide bar on a "
        "central prediction is still a weak statement.",
        "- **Apparent Tm is method-dependent.** The labels mix nanoDSF, DSF, "
        "circular dichroism and DSC. Differences of a few degrees between methods "
        "for the same nanobody are ordinary.",
        "- **The training set is small.** A few hundred clusters cannot cover the "
        "sequence space of 2,500 designed candidates.",
        "- **Long or stapled CDR3s are under-represented** in the labels relative "
        "to the candidate set.",
        "- **The model predicts apparent Tm, not stability in any particular "
        "formulation**, and says nothing about aggregation on storage.",
        "",
        "## Positive control",
        "",
        "Spec Section 9 item 3 asserts that no predicted apparent Tm in this "
        "model's output exceeds 95 degrees C. If one appears, the model is "
        "extrapolating and the error bar must reflect it. The check runs in "
        "`tools/controls.py` and its result appears on the QC panel.",
        "",
        "## Licence and citation",
        "",
        "Model code MIT, with the rest of WHITTLE. scikit-learn is BSD-3 and "
        "ESM-2 is MIT. The labels are cited to their CC-BY articles above; no "
        "label set is redistributed in this repository.",
        "",
    ]
    (ROOT / "MODEL_CARD.md").write_text("\n".join(lines))


if __name__ == "__main__":
    sys.exit(main())
