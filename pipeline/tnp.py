"""TNP, the Therapeutic Nanobody Profiler. Spec 6.2 item 2.

TNP is the whole non-specific-binding axis. It takes a FASTA directly, folds
each sequence with NanoBodyBuilder2 itself, and reports six metrics against
thresholds calibrated on 36 clinical-stage nanobodies:

    Total CDR Length, CDR3 Length, CDR3 Compactness, PSH, PPC, PNC

each with a green, amber or red flag (grey when a metric could not be
computed). `tnp_flags` in the database is the count of metrics outside the
clinical-stage green band, which is the single number the specificity axis
uses; the six raw values and their individual flags are kept for the drawer.

One call per target, as the spec asks, because TNP amortises its model loading
across a FASTA and a per-sequence call would pay it 500 times.

Note that TNP folds the sequences again rather than reading the models
`pipeline.fold` produced: its API is sequence in, metrics out, and its PSH, PPC
and PNC are defined against its own models. Re-folding is the honest way to
keep those numbers the ones TNP's thresholds were calibrated against.

    pixi run python -m pipeline.tnp
    pixi run python -m pipeline.tnp --target GLP1R
"""
from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import time
from pathlib import Path

from pipeline.common import (DATA, TARGETS, connect, log, manifest_read,
                             manifest_write, tuning)

TNP_DIR = DATA / "tnp"

# Thresholds from TNP's own `thresholds.py`, recorded here only so the About
# tab can state what a flag means. TNP computes the flags; nothing below
# recomputes them.
GREEN_BANDS = {
    "L": "24 to 38 total CDR residues",
    "L3": "8 to 22 CDR3 residues",
    "C": "CDR3 compactness 0.81 to 1.57",
    "PSH": "patch surface hydrophobicity 79.6 to 126.8",
    "PPC": "patch positive charge up to 0.39",
    "PNC": "patch negative charge up to 1.47",
}
METRIC_KEYS = ["Total CDR Length", "CDR3 Length", "CDR3 Compactness",
               "PSH", "PPC", "PNC"]


def tnp_version() -> str:
    try:
        import theraprofnano

        return str(getattr(theraprofnano, "__version__", "0.1.0"))
    except Exception:
        return "0.1.0"


def write_fasta(conn, target: str, path: Path) -> int:
    rows = conn.execute(
        "SELECT id, sequence FROM candidates WHERE target = ? ORDER BY id",
        (target,)).fetchall()
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as fh:
        for row in rows:
            fh.write(f">{row['id']}\n{row['sequence']}\n")
    return len(rows)


def run_target(target: str, fasta: Path, out_dir: Path, ncores: int) -> dict:
    started = time.time()
    binary = shutil.which("TNP")
    if binary is None:
        return {"id": target, "status": "failed:TNP not on PATH"}
    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = [binary, "--file", str(fasta), "--output", str(out_dir),
           "--ncores", str(ncores)]
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=60 * 300)
    except subprocess.TimeoutExpired:
        return {"id": target, "status": "failed:timeout"}
    results = out_dir / "TNP_Results_Multientry.json"
    if not results.exists():
        return {"id": target, "status": f"failed:no output (rc {proc.returncode})",
                "detail": (proc.stderr or proc.stdout)[-400:]}
    return {"id": target, "status": "ok",
            "seconds": round(time.time() - started, 1),
            "results": str(results)}


def ingest(target: str, results_path: Path) -> int:
    """One manifest row per candidate, so Tier 1 can pick them up."""
    payload = json.loads(results_path.read_text())
    written = 0
    for cid, rec in payload.items():
        flags = rec.get("Flags", {}) or {}
        # The specificity axis uses one number: how many of the six metrics sit
        # outside the clinical-stage green band. Grey means TNP could not
        # compute the metric and is not counted as a failure.
        non_green = sum(1 for v in flags.values() if v in ("amber", "red"))
        manifest_write("tnp", {
            "id": cid, "target": target, "status": "ok",
            "psh": rec.get("PSH"), "ppc": rec.get("PPC"), "pnc": rec.get("PNC"),
            "total_cdr_length": rec.get("Total CDR Length"),
            "cdr3_length": rec.get("CDR3 Length"),
            "cdr3_compactness": rec.get("CDR3 Compactness"),
            "flags": non_green,
            "flag_detail": flags,
            "n_red": sum(1 for v in flags.values() if v == "red"),
            "n_grey": sum(1 for v in flags.values() if v == "grey"),
            "version": tnp_version(),
        })
        written += 1
    return written


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--target", action="append", default=None)
    parser.add_argument("--ncores", type=int, default=None)
    args = parser.parse_args()

    targets = args.target or TARGETS
    ncores = args.ncores or max(2, min(int(tuning().get("cpu_workers", 8)), 12))
    done = {r["id"] for r in manifest_read("tnp_target") if r.get("status") == "ok"}

    conn = connect()
    for target in targets:
        if target in done:
            log(f"tnp: {target} already done, skipping")
            continue
        fasta = TNP_DIR / f"{target}.fasta"
        out_dir = TNP_DIR / target
        n = write_fasta(conn, target, fasta)
        log(f"tnp: {target}, {n} sequences, {ncores} cores. TNP folds each "
            f"sequence itself, so this is the slow step")
        record = run_target(target, fasta, out_dir, ncores)
        if record["status"] != "ok":
            manifest_write("tnp_target", record)
            log(f"tnp: {target} FAILED {record['status']} "
                f"{record.get('detail', '')[:200]}")
            continue
        written = ingest(target, Path(record["results"]))
        record["rows"] = written
        manifest_write("tnp_target", record)
        log(f"tnp: {target} complete, {written} of {n} sequences profiled in "
            f"{record['seconds'] / 60:.1f} minutes")
    conn.close()

    rows = [r for r in manifest_read("tnp") if r.get("status") == "ok"]
    log(f"tnp: {len(rows)} candidates profiled in total")
    return 0


if __name__ == "__main__":
    sys.exit(main())
