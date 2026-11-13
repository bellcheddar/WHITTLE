"""Nativeness and humanness. Spec 6.2 items 5 and 6.

**Nativeness** is the pseudo-perplexity of the sequence under AbLang2
(BSD-3, Olsen et al.), computed with the model's own
`pseudo_log_likelihood` mode: one masked copy per residue, so every position is
scored against the rest of the sequence. Odd predicted sequences score badly
for a reason, and this is the best generic sanity filter in the whole pipeline.
Reported as perplexity (exp of the negative mean log-likelihood), so lower is
better.

The caveat belongs on the About tab and is stated in `DECISIONS.md` 19: AbLang2
is trained on paired human and mouse antibodies, so for a VHH it is a general
antibody-nativeness prior rather than a nanobody-specific one. nanoBERT would
be the better-matched model but is CC-BY-NC-SA and cannot ship in a public,
permissively licensed repository.

**Humanness** is the OASis peptide identity: the fraction of the sequence's
9-mers that occur in at least 10 percent of the human repertoire subjects in
the Observed Antibody Space. Computed with `promb` (MIT, by the BioPhi
authors), which carries the same peptide set in memory rather than requiring
BioPhi's 22 GB SQLite. For a VHH this is a descriptor, not a validated
humanness call, and the app says so.

    pixi run python -m pipeline.nativeness
    pixi run python -m pipeline.nativeness --target GLP1R --limit 50
"""
from __future__ import annotations

import argparse
import math
import multiprocessing as mp
import os
import sys
import time
import warnings

from pipeline.common import (TARGETS, connect, log, manifest_done, manifest_write,
                             tuning)

warnings.filterwarnings("ignore")

_model = None
_threads = 4


def _init_worker(threads: int) -> None:
    global _model, _threads
    _threads = threads
    os.environ.setdefault("TOKENIZERS_PARALLELISM", "false")
    import torch

    torch.set_num_threads(threads)
    import ablang2

    # device='cpu' is deliberate, not a fallback. AbLang2's encoder calls
    # .numpy() on a tensor without moving it off the device first, and its
    # pseudo-log-likelihood path uses torch.isin, which has been patchy on MPS.
    # float32 on CPU is both correct and the project's numeric regime anyway.
    _model = ablang2.pretrained(model_to_use="ablang2-paired", random_init=False,
                                ncpu=threads, device="cpu")
    _model.freeze()


def _score_batch(job: tuple[list[str], list[str]]) -> list[dict]:
    import numpy as np

    ids, sequences = job
    try:
        # A VHH is a single heavy chain; the light slot is empty, which AbLang2
        # reads as an unpaired heavy chain.
        values = _model([[s, ""] for s in sequences], mode="pseudo_log_likelihood")
    except Exception as exc:
        return [{"id": cid, "status": f"failed:ablang2:{type(exc).__name__}",
                 "detail": str(exc)[:160]} for cid in ids]
    out = []
    for cid, pll in zip(ids, np.atleast_1d(values)):
        pll = float(pll)
        out.append({
            "id": cid, "status": "ok",
            "pseudo_log_likelihood": round(pll, 5),
            # AbLang2 returns the per-residue mean, so exp(-pll) is the
            # pseudo-perplexity directly.
            "pseudo_perplexity": round(math.exp(-pll), 4),
            "model": "ablang2-paired",
        })
    return out


def run_nativeness(targets: list[str], limit: int | None, workers: int,
                   threads: int, batch: int) -> None:
    conn = connect()
    done = manifest_done("nativeness")
    jobs: list[tuple[list[str], list[str]]] = []
    pending: list[tuple[str, str]] = []
    for target in targets:
        for row in conn.execute(
                "SELECT id, sequence FROM candidates WHERE target = ? ORDER BY id",
                (target,)):
            if row["id"] in done:
                continue
            pending.append((row["id"], row["sequence"]))
            if limit and len(pending) >= limit:
                break
        if limit and len(pending) >= limit:
            break
    conn.close()
    if not pending:
        log("nativeness: nothing to do")
        return
    for i in range(0, len(pending), batch):
        chunk = pending[i:i + batch]
        jobs.append(([c[0] for c in chunk], [c[1] for c in chunk]))

    log(f"nativeness: AbLang2 pseudo-perplexity for {len(pending)} sequences, "
        f"{workers} workers x {threads} threads, batches of {batch}")
    started = time.time()
    scored = failed = 0
    ctx = mp.get_context("spawn")
    with ctx.Pool(workers, initializer=_init_worker, initargs=(threads,)) as pool:
        for n, records in enumerate(pool.imap_unordered(_score_batch, jobs), 1):
            for rec in records:
                manifest_write("nativeness", rec)
                if rec["status"] == "ok":
                    scored += 1
                else:
                    failed += 1
            if n % 5 == 0 or n == len(jobs):
                elapsed = time.time() - started
                rate = (scored + failed) / max(elapsed, 1e-9)
                left = (len(pending) - scored - failed) / max(rate, 1e-9)
                log(f"nativeness: {scored + failed}/{len(pending)}, "
                    f"{rate * 60:.0f} per minute, about {left / 60:.0f} minutes left")
    log(f"nativeness: {scored} scored, {failed} failed, in "
        f"{(time.time() - started) / 60:.1f} minutes")


def run_humanness(targets: list[str]) -> None:
    """OASis peptide identity. Fast enough to be a single pass."""
    try:
        from promb import db, init_db
    except Exception as exc:
        log(f"humanness: promb unavailable ({exc}); the column stays pending")
        return
    started = time.time()
    oasis = init_db(db.HUMAN_OASIS_DB_PATH, verbose=False)
    conn = connect()
    done = manifest_done("humanness")
    n = 0
    for target in targets:
        for row in conn.execute(
                "SELECT id, sequence FROM candidates WHERE target = ? ORDER BY id",
                (target,)):
            if row["id"] in done:
                continue
            try:
                identity = float(oasis.compute_peptide_content(row["sequence"]))
                peptides = len(oasis.chop_seq_peptides(row["sequence"]))
            except Exception as exc:
                manifest_write("humanness", {
                    "id": row["id"], "status": f"failed:promb:{type(exc).__name__}"})
                continue
            manifest_write("humanness", {
                "id": row["id"], "target": target, "status": "ok",
                "oasis_identity": round(identity, 5),
                "n_peptides": peptides,
                "peptide_length": 9,
                "min_percent_subjects": 10,
                "source": "promb, OASis 9-mers at 10 percent of subjects",
            })
            n += 1
    conn.close()
    log(f"humanness: OASis peptide identity for {n} sequences in "
        f"{time.time() - started:.0f} seconds")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--target", action="append", default=None)
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--threads", type=int, default=None)
    parser.add_argument("--batch", type=int, default=16)
    parser.add_argument("--skip-nativeness", action="store_true")
    parser.add_argument("--skip-humanness", action="store_true")
    args = parser.parse_args()

    targets = args.target or TARGETS
    cores = int(tuning().get("cpu_workers", 8))
    threads = args.threads or max(2, cores // max(args.workers, 1))

    if not args.skip_humanness:
        run_humanness(targets)
    if not args.skip_nativeness:
        run_nativeness(targets, args.limit, args.workers, threads, args.batch)
    return 0


if __name__ == "__main__":
    sys.exit(main())
