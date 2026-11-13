"""Tier 2: co-folding with boltz-community. Spec 6.4.

A multi-day, resumable queue. Design follows BoltzMaker's campaign driver (see
REUSE.md): one input per complex, a manifest row per complex, and a restart
that reads the manifest and skips what is done. Assume the process is killed at
any moment; it is, repeatedly.

Order is **cheapest first: GLP1R, PRLR, TNFa, TL1A, EGFR** (spec 6.4 item 2),
so a broken setup surfaces on a 241 aa complex in minutes rather than on EGFR
overnight.

The antigen MSA is built **once per target** and reused across every candidate
for that target (spec 6.4 item 1).

Fallback ladder for an out-of-memory complex, in order, per spec 0.2:
  1. drop to a single diffusion sample
  2. add `--no_kernels`
  3. skip the complex and log it
**Never** substitute a monomer for a trimer. For TNFa and TL1A the trimer is
the point: the epitope commonly spans two subunits.

    pixi run -e cofold python -m pipeline.cofold --smoke
    pixi run -e cofold python -m pipeline.cofold
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

from pipeline import boltz_input
from pipeline.common import (MSA_CACHE, RUNS, TIER2_ORDER, connect, log,
                             manifest_read, manifest_write, targets_config,
                             tuning)

INPUTS = RUNS / "inputs"


def boltz_binary() -> str | None:
    """The cofold environment's boltz, not whatever is first on PATH."""
    from pipeline.common import ROOT

    local = ROOT / ".pixi" / "envs" / "cofold" / "bin" / "boltz"
    if local.exists():
        return str(local)
    return shutil.which("boltz")


def verify_fork() -> dict:
    """Spec 4.3: verify `boltz` is the community fork before queueing anything."""
    from pipeline.common import ROOT

    python = ROOT / ".pixi" / "envs" / "cofold" / "bin" / "python"
    if not python.exists():
        return {"ok": False, "detail": "no cofold environment"}
    script = (
        "import importlib.metadata as md, json, torch\n"
        "names = {d.metadata['Name'].lower() for d in md.distributions() "
        "if d.metadata and d.metadata.get('Name')}\n"
        "out = {'community': 'boltz-community' in names, "
        "'upstream': 'boltz' in names and 'boltz-community' not in names, "
        "'torch': torch.__version__, 'mps': torch.backends.mps.is_available()}\n"
        "out['version'] = md.version('boltz-community') if out['community'] else None\n"
        "print(json.dumps(out))\n")
    try:
        proc = subprocess.run([str(python), "-c", script], capture_output=True,
                              text=True, timeout=180)
        info = json.loads(proc.stdout.strip().splitlines()[-1])
    except Exception as exc:
        return {"ok": False, "detail": f"{type(exc).__name__}: {exc}"}
    info["ok"] = bool(info.get("community")) and not info.get("upstream")
    return info


# --------------------------------------------------------------------- MSA

def build_msa(target: str, binary: str, env: dict) -> bool:
    """Build the antigen MSA once for this target and cache it.

    Runs boltz on the antigen alone with `--use_msa_server`, which writes the
    MSA next to its processed inputs; the a3m is then copied into
    `data/msa_cache/<target>.a3m` and every candidate for this target reuses it.
    """
    cached = MSA_CACHE / f"{target}.a3m"
    if cached.exists() and cached.stat().st_size > 0:
        return True
    MSA_CACHE.mkdir(parents=True, exist_ok=True)
    work = RUNS / "msa" / target
    work.mkdir(parents=True, exist_ok=True)
    spec = targets_config()["targets"][target]
    antigen = boltz_input.antigen_sequence(target)
    # One chain only: the MSA is per sequence, and the trimer shares it.
    doc = {"version": 1, "sequences": [{"protein": {"id": "A", "sequence": antigen}}]}
    import yaml

    input_path = work / f"{target}_antigen.yaml"
    input_path.write_text(yaml.safe_dump(doc, sort_keys=False))

    log(f"cofold: building the {target} antigen MSA once "
        f"({len(antigen)} aa, reused across every candidate)")
    cmd = [binary, "predict", str(input_path), "--out_dir", str(work),
           "--accelerator", "cpu", "--use_msa_server", "--num_workers", "0",
           "--diffusion_samples", "1", "--sampling_steps", "10",
           "--recycling_steps", "0", "--output_format", "pdb"]
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True,
                              timeout=60 * 90, env=env)
    except subprocess.TimeoutExpired:
        log(f"cofold: the MSA server timed out for {target}")
        return False
    found = sorted(work.rglob("*.a3m"))
    if not found:
        log(f"cofold: no MSA produced for {target} "
            f"(rc {proc.returncode}): {(proc.stderr or '')[-300:]}")
        return False
    biggest = max(found, key=lambda p: p.stat().st_size)
    shutil.copy(biggest, cached)
    depth = sum(1 for line in cached.read_text().splitlines()
                if line.startswith(">"))
    log(f"cofold: {target} antigen MSA cached, {depth} sequences deep")
    manifest_write("msa", {"id": target, "status": "ok", "depth": depth,
                           "path": str(cached)})
    return True


# ------------------------------------------------------------------- the run

def run_one(target: str, cid: str, sequence: str, binary: str, env: dict,
            samples: int, timeout: int) -> dict:
    """One complex, with the documented OOM ladder. Never a monomer fallback."""
    out_dir = RUNS / target / cid
    out_dir.mkdir(parents=True, exist_ok=True)
    input_path = boltz_input.write_input(target, cid, sequence, INPUTS / target)
    started = time.time()
    attempts: list[str] = []

    ladder = [
        ["--diffusion_samples", str(samples)],
        ["--diffusion_samples", "1"],
        ["--diffusion_samples", "1", "--no_kernels"],
    ]
    for i, extra in enumerate(ladder, 1):
        cmd = [binary, "predict", str(input_path), "--out_dir", str(out_dir),
               "--accelerator", str(tuning().get("cofold_accelerator", "mps")),
               "--recycling_steps", str(tuning().get("cofold_recycling", 3)),
               "--output_format", "pdb", "--num_workers", "0",
               "--seed", "0", "--override"] + extra
        try:
            proc = subprocess.run(cmd, capture_output=True, text=True,
                                  timeout=timeout, env=env)
        except subprocess.TimeoutExpired:
            attempts.append(f"attempt {i}: timeout")
            continue
        if proc.returncode == 0 and list(out_dir.rglob("confidence*.json")):
            return {"id": cid, "target": target, "status": "ok",
                    "samples": int(extra[1]), "attempt": i,
                    "seconds": round(time.time() - started, 1),
                    "path": str(out_dir), "ladder": attempts}
        tail = (proc.stderr or proc.stdout or "")[-300:]
        attempts.append(f"attempt {i}: rc {proc.returncode} {tail[-160:]}")
        low_memory = any(w in tail.lower() for w in
                         ("out of memory", "mps backend out of memory", "alloc"))
        if not low_memory and i == 1:
            # Not a memory problem: the next rungs of the ladder are for OOM,
            # so do not burn two more attempts on a different fault.
            break
    return {"id": cid, "target": target, "status": "failed:cofold",
            "seconds": round(time.time() - started, 1), "ladder": attempts}


def queue(conn, targets: list[str], per_target: int) -> list[tuple[str, str, str]]:
    """Top N per target by desirability, cheapest target first (spec 6.4 item 2)."""
    done = {str(r["id"]) for r in manifest_read("cofold")
            if r.get("status") == "ok"}
    jobs: list[tuple[str, str, str]] = []
    for target in targets:
        rows = conn.execute(
            "SELECT c.id, c.sequence FROM candidates c "
            "JOIN ranks r ON r.candidate_id = c.id "
            "WHERE c.target = ? ORDER BY r.rank LIMIT ?",
            (target, per_target)).fetchall()
        for row in rows:
            if row["id"] in done:
                continue
            jobs.append((target, row["id"], row["sequence"]))
    return jobs


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--target", action="append", default=None)
    parser.add_argument("--top-n", type=int, default=None)
    parser.add_argument("--samples", type=int, default=None)
    parser.add_argument("--timeout", type=int, default=60 * 240)
    parser.add_argument("--smoke", action="store_true",
                        help="one GLP1R complex, timed, before anything is queued")
    parser.add_argument("--msa-only", action="store_true")
    args = parser.parse_args()

    cfg = targets_config()
    order = args.target or [t for t in TIER2_ORDER
                            if t in cfg["targets"]]
    per_target = args.top_n or int(cfg["run"].get("tier2_top_n_per_target", 20))
    samples = args.samples or int(cfg["tier2"].get("seeds_per_complex", 5))

    info = verify_fork()
    if not info.get("ok"):
        log(f"cofold: boltz is not the community fork or is missing: {info}. "
            f"Spec 4.2 forbids upstream boltz. Run: pixi run -e cofold install-boltz")
        return 1
    log(f"cofold: boltz-community {info.get('version')}, torch {info.get('torch')}, "
        f"MPS {info.get('mps')}")

    binary = boltz_binary()
    if binary is None:
        log("cofold: no boltz binary")
        return 1

    env = dict(os.environ)
    env.update({k: str(v) for k, v in (tuning().get("env") or {}).items()})
    env.setdefault("PYTORCH_ENABLE_MPS_FALLBACK", "1")

    conn = connect()
    for target in order:
        if not build_msa(target, binary, env):
            log(f"cofold: {target} has no cached MSA; deferring it to the end of "
                f"the queue (spec 0.2) and continuing with the others")
    if args.msa_only:
        conn.close()
        return 0

    if args.smoke:
        # Spec Phase 3: smoke-test one GLP1R complex and record the wall-clock
        # before queueing anything else.
        row = conn.execute(
            "SELECT c.id, c.sequence FROM candidates c JOIN ranks r "
            "ON r.candidate_id = c.id WHERE c.target = 'GLP1R' "
            "ORDER BY r.rank LIMIT 1").fetchone()
        log(f"cofold: smoke test on {row['id']}, the cheapest complex "
            f"({boltz_input.complex_length('GLP1R', row['sequence'])} aa)")
        record = run_one("GLP1R", row["id"], row["sequence"], binary, env,
                         samples, args.timeout)
        manifest_write("cofold", record)
        log(f"cofold: smoke test {record['status']} in "
            f"{record['seconds'] / 60:.1f} minutes")
        if record["status"] != "ok":
            log(f"cofold: ladder {record.get('ladder')}")
            conn.close()
            return 1
        # Project the whole queue from the one measurement.
        for target in order:
            length = boltz_input.complex_length(target, row["sequence"])
            base = boltz_input.complex_length("GLP1R", row["sequence"])
            # Attention is roughly quadratic in length (spec 5.2).
            estimate = record["seconds"] * (length / base) ** 2 * per_target / 3600
            log(f"cofold: {target} projected at about {estimate:.1f} hours for "
                f"{per_target} complexes ({length} aa each)")
        conn.close()
        return 0

    jobs = queue(conn, order, per_target)
    conn.close()
    if not jobs:
        log("cofold: queue empty, every complex is already done")
        return 0
    log(f"cofold: {len(jobs)} complexes queued, cheapest target first: "
        f"{', '.join(order)}")

    started = time.time()
    ok = failed = 0
    for n, (target, cid, sequence) in enumerate(jobs, 1):
        record = run_one(target, cid, sequence, binary, env, samples, args.timeout)
        manifest_write("cofold", record)
        if record["status"] == "ok":
            ok += 1
        else:
            failed += 1
            log(f"cofold: {cid} FAILED: {record.get('ladder')}")
        elapsed = time.time() - started
        rate = n / max(elapsed, 1e-9)
        log(f"cofold: {n}/{len(jobs)} ({ok} ok, {failed} failed), "
            f"{record['seconds'] / 60:.1f} min for {cid}, "
            f"about {(len(jobs) - n) / max(rate, 1e-9) / 3600:.1f} hours left")
    log(f"cofold: finished {ok} ok, {failed} failed, in "
        f"{(time.time() - started) / 3600:.1f} hours")
    return 0


if __name__ == "__main__":
    sys.exit(main())
