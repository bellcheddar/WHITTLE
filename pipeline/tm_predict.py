"""Apply the trained thermostability model. Spec 6.3 item 5.

Loads `models/whittle_tm_v1/model.joblib` and predicts apparent Tm for every
candidate, with a per-row error estimate from the distance to the nearest
training neighbour in the ESM-2 feature space. A candidate far from the
training manifold gets a wide bar, not a crisp number.

Also runs the TemStaPro cross-check where TemStaPro is present (spec 6.3 item
6). TemStaPro reports a probability of exceeding each of a ladder of
temperature thresholds rather than a temperature, so the two are reported side
by side and neither is blended into the other. If they disagree systematically
that is stated on the About tab rather than resolved by picking a winner.

    pixi run python -m pipeline.tm_predict
    pixi run python -m pipeline.tm_predict --controls
"""
from __future__ import annotations

import argparse
import sys
import time
import warnings

import numpy as np

from pipeline.common import MODELS, connect, log, manifest_write, tuning

warnings.filterwarnings("ignore")

MODEL_DIR = MODELS / "whittle_tm_v1"


def load_model():
    import joblib

    path = MODEL_DIR / "model.joblib"
    if not path.exists():
        return None
    return joblib.load(path)


def featurise(sequences: list[str], bundle: dict) -> tuple[np.ndarray, np.ndarray]:
    """(stacked features, ESM block) for a list of sequences."""
    from pipeline.tm_model import (esm_block, onehot_block, region_split,
                                   vhse_block)

    regions = [region_split(s) for s in sequences]
    oh = onehot_block(sequences, regions)
    vh = vhse_block(sequences, regions)
    es = esm_block(sequences, bundle["esm_model"], int(tuning().get("esm_batch", 8)))
    return np.hstack([oh, vh, es]), es


def predict(sequences: list[str], bundle: dict) -> tuple[np.ndarray, np.ndarray,
                                                         np.ndarray]:
    """(tm, error bar, distance to the nearest training neighbour)."""
    from pipeline.tm_model import apply_error, nearest_distance

    X, es = featurise(sequences, bundle)
    tm = bundle["model"].predict(X)
    distance = nearest_distance(es, bundle["train_esm"])
    error = apply_error(distance, bundle["calibration"])
    return tm, error, distance


def temstapro(sequences: list[str], ids: list[str]) -> dict[str, float]:
    """TemStaPro cross-check, if it is installed. Absent is fine.

    TemStaPro is a script repository with no package, so it is used only when a
    clone is present at `tools/TemStaPro` and its ProtT5 weights are cached.
    The column is left pending otherwise, and the About tab says the
    cross-check did not run rather than implying agreement.
    """
    import shutil
    import subprocess
    import tempfile
    from pathlib import Path

    from pipeline.common import ROOT

    repo = ROOT / "tools" / "TemStaPro"
    script = repo / "temstapro"
    if not script.exists():
        log("tm_predict: TemStaPro is not installed; the cross-check column "
            "stays pending")
        return {}
    weights = repo / "ProtTrans"
    if not weights.exists():
        log("tm_predict: TemStaPro's ProtT5 weights are absent; skipped")
        return {}
    with tempfile.TemporaryDirectory() as tmp:
        fasta = Path(tmp) / "in.fasta"
        fasta.write_text("".join(f">{i}\n{s}\n" for i, s in zip(ids, sequences)))
        out = Path(tmp) / "out.tsv"
        try:
            subprocess.run(
                [str(script), "-f", str(fasta), "-d", str(weights),
                 "-t", str(repo), "-e", str(Path(tmp) / "emb"),
                 "--mean-output", str(out)],
                capture_output=True, check=True, timeout=60 * 600, cwd=str(repo))
        except Exception as exc:
            log(f"tm_predict: TemStaPro failed ({type(exc).__name__}); skipped")
            return {}
        if not out.exists():
            return {}
        values: dict[str, float] = {}
        header: list[str] = []
        for line in out.read_text().splitlines():
            parts = line.split("\t")
            if not header:
                header = parts
                continue
            row = dict(zip(header, parts))
            name = row.get(header[0], "")
            # The raw probability of exceeding 65 degrees C is the single
            # number that compares best with a predicted apparent Tm.
            for key in ("raw65", "65_raw", "pred65_raw"):
                if key in row:
                    try:
                        values[name] = float(row[key])
                    except ValueError:
                        pass
                    break
        return values


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--controls", action="store_true",
                        help="predict for the control rows instead of the candidates")
    parser.add_argument("--batch", type=int, default=250)
    args = parser.parse_args()

    bundle = load_model()
    if bundle is None:
        log("tm_predict: no trained model at models/whittle_tm_v1; run tm-train first")
        return 1

    conn = connect()
    if args.controls:
        rows = conn.execute(
            "SELECT id, sequence FROM controls WHERE sequence IS NOT NULL "
            "ORDER BY id").fetchall()
    else:
        rows = conn.execute(
            "SELECT id, sequence FROM candidates ORDER BY id").fetchall()
    ids = [r["id"] for r in rows]
    sequences = [r["sequence"] for r in rows]
    log(f"tm_predict: {len(ids)} sequences, ESM-2 {bundle['esm_model']}")

    started = time.time()
    all_tm: list[float] = []
    all_err: list[float] = []
    all_dist: list[float] = []
    for i in range(0, len(sequences), args.batch):
        chunk = sequences[i:i + args.batch]
        tm, err, dist = predict(chunk, bundle)
        all_tm += [float(x) for x in tm]
        all_err += [float(x) for x in err]
        all_dist += [float(x) for x in dist]
        log(f"tm_predict: {min(i + args.batch, len(sequences))}/{len(sequences)}")

    cross = temstapro(sequences, ids)

    if args.controls:
        # The control table keeps the temperature itself, which is what the
        # observed-against-predicted scatter on Tab 2 plots.
        conn.executemany(
            "UPDATE controls SET predicted_tm = ? WHERE id = ?",
            [(round(t, 2), cid) for cid, t in zip(ids, all_tm)])
    else:
        conn.executemany(
            "INSERT OR REPLACE INTO thermostability (candidate_id, tm_pred, "
            "tm_err, train_nn_distance, temstapro, thermompnn_ddg, model_version, "
            "status) VALUES (?,?,?,?,?,?,?,?)",
            [(cid, round(t, 2), round(e, 2), round(d, 4), cross.get(cid),
              None, "whittle_tm_v1", "ok")
             for cid, t, e, d in zip(ids, all_tm, all_err, all_dist)])
    conn.commit()

    hottest = max(all_tm)
    over95 = sum(1 for t in all_tm if t > 95)
    log(f"tm_predict: {len(ids)} predictions in "
        f"{(time.time() - started) / 60:.1f} minutes. Tm {min(all_tm):.1f} to "
        f"{hottest:.1f} deg C, mean error bar "
        f"{sum(all_err) / len(all_err):.1f} deg C, {over95} above 95 deg C "
        f"(spec 9 item 3 wants none)")
    manifest_write("tm_predict", {
        "id": "controls" if args.controls else "candidates", "status": "ok",
        "n": len(ids), "min": round(min(all_tm), 2), "max": round(hottest, 2),
        "over_95": over95, "temstapro": len(cross)})
    conn.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
