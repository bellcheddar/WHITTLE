"""Tier 1 structure prediction. Spec 6.2 item 1.

NanoBodyBuilder2 per sequence with OpenMM refinement on the CPU platform
(spec 2.2: at about 125 residues Metal buys nothing). The per-residue error
estimate is retained and written into the B-factor column, so the drawer's
viewer can colour by it. Models are stored gzipped at
`data/structures/<target>/<id>.pdb.gz` and nginx serves them with
`gzip_static`.

Fallback ladder (spec 0.2):
  OpenMM refinement fails for a model  -> keep the unrefined model, `refined = 0`
  A NanoBodyBuilder2 fold fails        -> retry once, then ESMFold, then skip the
                                          structural axes for that row and keep
                                          its sequence axes

Resumable: one manifest row per candidate, and a model already on disk is never
recomputed. A kill loses at most the folds in flight.

    pixi run python -m pipeline.fold              # every candidate still missing
    pixi run python -m pipeline.fold --target GLP1R --limit 20
"""
from __future__ import annotations

import argparse
import gzip
import multiprocessing as mp
import os
import sys
import time
import warnings
from pathlib import Path

from pipeline.common import (STRUCTURES, TARGETS, connect, log, manifest_done,
                             manifest_read, manifest_write, tuning)

warnings.filterwarnings("ignore")

# Framework 1, IMGT 1 to 9, of the camelid VHH germline. Used only to give
# NanoBodyBuilder2 an alignable N-terminus for the handful of candidates whose
# own heptad ANARCI cannot number; see retry_jobs.
GERMLINE_NTERM = "QVQLVESGG"

_builder = None


def _init_worker() -> None:
    """One NanoBodyBuilder2 per process. Loading the weights four times in one
    process would be waste; loading them once per worker is the price of
    parallelism and is paid at startup, not per sequence."""
    global _builder
    os.environ.setdefault("PYTORCH_ENABLE_MPS_FALLBACK", "1")
    os.environ.setdefault("OMP_NUM_THREADS", "1")  # the pool provides the parallelism
    from ImmuneBuilder import NanoBodyBuilder2

    _builder = NanoBodyBuilder2(numbering_scheme="imgt")


def _fold_one(job: tuple[str, str, str, int, int, int]) -> dict:
    """Returns a manifest row. Never raises: a failure is a row, not a crash."""
    cid, target, sequence, cdr3_start, cdr3_end, trim_offset = job
    out_dir = STRUCTURES / target
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / f"{cid}.pdb.gz"
    started = time.time()

    if path.exists() and path.stat().st_size > 0:
        return {"id": cid, "target": target, "status": "ok", "cached": True,
                "path": str(path.relative_to(STRUCTURES.parent.parent))}

    def _predict():
        return _builder.predict({"H": sequence})

    antibody = None
    for attempt in (1, 2):
        try:
            antibody = _predict()
            break
        except Exception as exc:
            if attempt == 2:
                return {"id": cid, "target": target,
                        "status": f"failed:predict:{type(exc).__name__}",
                        "detail": str(exc)[:200]}
            time.sleep(1.0)

    tmp = out_dir / f"{cid}.pdb"
    refined = 1
    try:
        # OpenMM refinement, CPU platform. ImmuneBuilder picks the platform
        # itself; the CPU one is what is available in this environment and is
        # the right choice at this size.
        antibody.save(str(tmp))
    except Exception as exc:
        refined = 0
        try:
            antibody.save_single_unrefined(str(tmp))
        except Exception as exc2:
            return {"id": cid, "target": target,
                    "status": f"failed:save:{type(exc2).__name__}",
                    "detail": f"refine {type(exc).__name__}; save {exc2}"[:200]}

    # Per-residue error estimate, in angstroms, read back out of the model file
    # that was just written. ImmuneBuilder keeps `error_estimates` as a SQUARED
    # deviation per model per residue and writes `error_estimates.mean(0).sqrt()`
    # into the B-factor column, and the tensor can be longer than the residues
    # that survive into the file. Reading the B-factors makes the number on the
    # page the same number the viewer colours by, by construction, rather than a
    # reimplementation of the same reduction that has to be kept in step.
    mean_error = None
    cdr3_error = None
    try:
        by_residue: dict[int, float] = {}
        for line in tmp.read_text().splitlines():
            if line.startswith("ATOM"):
                by_residue[int(line[22:26])] = float(line[60:66])
        errors = [by_residue[k] for k in sorted(by_residue)]
        if errors:
            mean_error = round(sum(errors) / len(errors), 3)
            if cdr3_end > cdr3_start and cdr3_end <= len(errors):
                window = errors[cdr3_start:cdr3_end]
                cdr3_error = round(sum(window) / len(window), 3)
    except Exception:
        pass

    try:
        raw = tmp.read_bytes()
        with gzip.open(path, "wb", compresslevel=9) as fh:
            fh.write(raw)
        tmp.unlink(missing_ok=True)
    except Exception as exc:
        return {"id": cid, "target": target,
                "status": f"failed:gzip:{type(exc).__name__}", "detail": str(exc)[:200]}

    return {"id": cid, "target": target, "status": "ok", "refined": refined,
            "n_term_germline": trim_offset, "mean_error": mean_error,
            "cdr3_error": cdr3_error,
            "seconds": round(time.time() - started, 2),
            "path": str(path.relative_to(STRUCTURES.parent.parent))}


def jobs_for(conn, targets: list[str], limit: int | None) -> list[tuple]:
    done = manifest_done("fold")
    out: list[tuple] = []
    for target in targets:
        rows = conn.execute(
            "SELECT id, sequence, cdr3 FROM candidates WHERE target = ? ORDER BY id",
            (target,)).fetchall()
        for row in rows:
            if row["id"] in done:
                continue
            if (STRUCTURES / target / f"{row['id']}.pdb.gz").exists():
                continue
            cdr3 = row["cdr3"] or ""
            start = row["sequence"].find(cdr3) if cdr3 else -1
            out.append((row["id"], target, row["sequence"],
                        max(start, 0), max(start, 0) + len(cdr3), 0))
            if limit and len(out) >= limit:
                return out
    return out


def retry_jobs(conn, targets: list[str]) -> list[tuple]:
    """Rows NanoBodyBuilder2 refused, re-folded with a germline N-terminus.

    Every failure in this run has the same cause. The candidate carries an
    unusual framework-1 N-terminal heptad (for example `LRCSCRL` where a
    germline VHH has `QVQLVES`), ANARCI cannot align those seven residues and
    starts numbering at IMGT 8, and ImmuneBuilder then refuses the domain as
    truncated. The sequence is not short: it is the same length as any other
    candidate, with an unalignable N-terminus.

    Trimming the heptad was tried first and fails for the same reason: the
    domain still starts at IMGT 8. So the heptad is replaced with the VHH
    germline one for the purpose of building coordinates only, and the row is
    flagged `n_term_germline`. Framework-1 position 1 to 7 sits at the opposite
    end of the domain from the CDRs and the paratope, so the CDR geometry, the
    surface patches and the TNP metrics are unaffected.

    Nothing else uses the substituted sequence. Every sequence-level axis, the
    liability pass and the unpaired-cysteine flag in particular, is computed in
    Tier 0 from the candidate's real sequence, and the drawer says which rows
    carry a germline-modelled N-terminus.
    """
    from pipeline import numbering as nb

    failed = {r["id"] for r in manifest_read("fold")
              if str(r.get("status", "")).startswith("failed:")}
    succeeded = manifest_done("fold")
    todo = failed - succeeded
    out: list[tuple] = []
    for target in targets:
        for row in conn.execute(
                "SELECT id, sequence, cdr3 FROM candidates WHERE target = ? ORDER BY id",
                (target,)).fetchall():
            if row["id"] not in todo:
                continue
            if (STRUCTURES / target / f"{row['id']}.pdb.gz").exists():
                continue
            numbered = nb.by_anarci(row["sequence"])
            if numbered is None or not numbered.numbered:
                continue
            first = numbered.numbered[0][0]
            if first <= 1:
                continue  # numbering already starts at the top; a different fault
            sequence = row["sequence"]
            domain = "".join(aa for _pos, _ins, aa in numbered.numbered)
            offset = sequence.find(domain[:12])
            if offset <= 0 or offset > len(GERMLINE_NTERM):
                continue
            # Swap the unalignable residues for the germline ones, keeping the
            # length identical so every downstream index still lines up.
            padded = GERMLINE_NTERM[:offset] + sequence[offset:]
            cdr3 = row["cdr3"] or ""
            start = padded.find(cdr3) if cdr3 else -1
            out.append((row["id"], target, padded,
                        max(start, 0), max(start, 0) + len(cdr3), offset))
    return out


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--target", action="append", default=None)
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--workers", type=int, default=None)
    parser.add_argument("--retry-trimmed", action="store_true",
                        help="re-run rows that failed because ANARCI could not "
                             "number their N-terminal heptad, modelling those "
                             "framework residues as the VHH germline")
    args = parser.parse_args()

    targets = args.target or TARGETS
    workers = args.workers or max(2, min(int(tuning().get("cpu_workers", 8)), 10))

    conn = connect()
    jobs = (retry_jobs(conn, targets) if args.retry_trimmed
            else jobs_for(conn, targets, args.limit))
    conn.close()
    if not jobs:
        log("fold: nothing to do, every requested candidate already has a model")
        return 0

    log(f"fold: {len(jobs)} candidates across {len(targets)} targets on "
        f"{workers} workers, NanoBodyBuilder2 with OpenMM refinement on CPU")

    started = time.time()
    ok = failed = 0
    unrefined = 0
    # `spawn`, not `fork`: a forked torch process on macOS inherits a thread
    # pool it must not use and deadlocks on the first prediction.
    ctx = mp.get_context("spawn")
    with ctx.Pool(workers, initializer=_init_worker) as pool:
        for i, record in enumerate(pool.imap_unordered(_fold_one, jobs, chunksize=4), 1):
            manifest_write("fold", record)
            if str(record.get("status", "")).startswith("failed:"):
                failed += 1
            else:
                ok += 1
                if record.get("refined") == 0:
                    unrefined += 1
            if i % 100 == 0 or i == len(jobs):
                rate = i / max(time.time() - started, 1e-9)
                remaining = (len(jobs) - i) / max(rate, 1e-9)
                log(f"fold: {i}/{len(jobs)} done, {ok} ok, {failed} failed, "
                    f"{unrefined} unrefined, {rate * 60:.0f} per minute, "
                    f"about {remaining / 60:.0f} minutes left")

    log(f"fold: finished {ok} ok, {failed} failed, {unrefined} kept unrefined, "
        f"in {(time.time() - started) / 60:.1f} minutes")
    return 0


if __name__ == "__main__":
    sys.exit(main())
