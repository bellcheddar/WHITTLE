"""Measure what each co-folding knob actually costs. Spec 5.2 and 6.4.

The Tier 2 queue is the long pole, so the question "how do we make it shorter"
has to be answered with measurement rather than intuition. Boltz-2's cost
splits into a trunk (recycling passes over the pairformer, whose triangle
operations scale badly with length) and a diffusion stage (samples times
sampling steps). Which knob helps depends on that split, and the split moves
with complex size: the trunk dominates more as the antigen gets larger.

Runs one complex under several settings and reports wall-clock and the
confidence each setting produced, so a saving can be weighed against what it
costs in answer quality.

    pixi run -e cofold python tools/cofold_bench.py --target GLP1R
"""
from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline import boltz_input  # noqa: E402
from pipeline.common import RUNS, connect, log  # noqa: E402

# (label, extra boltz flags). One knob at a time first, then the combination.
SETTINGS = [
    ("baseline r3 s5 st200", ["--recycling_steps", "3", "--diffusion_samples", "5"]),
    ("recycling 1", ["--recycling_steps", "1", "--diffusion_samples", "5"]),
    ("samples 2", ["--recycling_steps", "3", "--diffusion_samples", "2"]),
    ("steps 100", ["--recycling_steps", "3", "--diffusion_samples", "5",
                   "--sampling_steps", "100"]),
    ("msa 256", ["--recycling_steps", "3", "--diffusion_samples", "5",
                 "--num_subsampled_msa", "256"]),
    ("combined r1 s3 st100", ["--recycling_steps", "1", "--diffusion_samples", "3",
                              "--sampling_steps", "100"]),
]


def confidences(run_dir: Path) -> dict:
    out = {"iptm": [], "ptm": [], "plddt": []}
    for path in sorted(run_dir.rglob("confidence_*.json")):
        try:
            payload = json.loads(path.read_text())
        except json.JSONDecodeError:
            continue
        for key, field in (("iptm", "iptm"), ("ptm", "ptm"),
                           ("plddt", "complex_plddt")):
            if isinstance(payload.get(field), (int, float)):
                out[key].append(round(float(payload[field]), 4))
    return {k: (round(sum(v) / len(v), 4) if v else None) for k, v in out.items()}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--target", default="GLP1R")
    parser.add_argument("--only", action="append", default=None)
    args = parser.parse_args()

    binary = str(ROOT / ".pixi" / "envs" / "cofold" / "bin" / "boltz")
    conn = connect()
    row = conn.execute(
        "SELECT c.id, c.sequence FROM candidates c JOIN ranks r "
        "ON r.candidate_id = c.id WHERE c.target = ? ORDER BY r.rank LIMIT 1",
        (args.target,)).fetchone()
    conn.close()
    cid, sequence = row["id"], row["sequence"]
    length = boltz_input.complex_length(args.target, sequence)
    log(f"cofold_bench: {args.target} {cid}, {length} aa complex")

    bench = RUNS / "_bench"
    input_path = boltz_input.write_input(args.target, cid, sequence, bench / "in")

    import os

    env = dict(os.environ)
    env.setdefault("PYTORCH_ENABLE_MPS_FALLBACK", "1")

    results = []
    for label, extra in SETTINGS:
        if args.only and label not in args.only:
            continue
        out_dir = bench / label.replace(" ", "_")
        shutil.rmtree(out_dir, ignore_errors=True)
        out_dir.mkdir(parents=True, exist_ok=True)
        cmd = [binary, "predict", str(input_path), "--out_dir", str(out_dir),
               "--accelerator", "mps", "--output_format", "pdb",
               "--num_workers", "0", "--seed", "0", "--override"] + extra
        started = time.time()
        proc = subprocess.run(cmd, capture_output=True, text=True,
                              timeout=60 * 180, env=env)
        seconds = time.time() - started
        conf = confidences(out_dir)
        ok = proc.returncode == 0 and any(out_dir.rglob("confidence_*.json"))
        results.append({"setting": label, "seconds": round(seconds, 1),
                        "minutes": round(seconds / 60, 2), "ok": ok, **conf})
        log(f"cofold_bench: {label:22s} {seconds / 60:6.2f} min  "
            f"ipTM {conf['iptm']}  pTM {conf['ptm']}  pLDDT {conf['plddt']}"
            + ("" if ok else "  FAILED"))

    base = next((r for r in results if r["setting"].startswith("baseline")), None)
    if base:
        log("cofold_bench: saving against baseline")
        for r in results:
            if r is base:
                continue
            saving = 1 - r["seconds"] / base["seconds"]
            d_iptm = (None if r["iptm"] is None or base["iptm"] is None
                      else round(r["iptm"] - base["iptm"], 4))
            log(f"cofold_bench:   {r['setting']:22s} {saving * 100:5.1f}% faster, "
                f"ipTM change {d_iptm}")
    (bench / "results.json").write_text(json.dumps(results, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
