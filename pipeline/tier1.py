"""Tier 1: structure level. Spec 6.2.

Folds are produced separately by `pipeline.fold`, because that is the hour-long
step and it has to be resumable on its own. This module gathers everything that
hangs off a finished model, plus the two sequence-level solubility predictors,
and writes the `structures` table:

  1. the model path, whether OpenMM refinement succeeded, and the per-residue
     error estimate (from `pipeline.fold`'s manifest)
  2. surface hydrophobic and charged patches (`pipeline.patches`)
  3. SWI and Protein-Sol (`pipeline.solubility`)
  4. TNP developability, where TNP is installed (`pipeline.tnp`)
  5. nativeness and humanness, where those are installed
  6. the nearest experimental neighbour in SAbDab-nano and the experimental
     badge (`pipeline.neighbours`)

Each of 4, 5 and 6 is optional and absent rather than faked. A column with no
value anywhere is reported as pending by the app, never as zero.

    pixi run tier1
    pixi run python -m pipeline.tier1 --target GLP1R
"""
from __future__ import annotations

import argparse
import multiprocessing as mp
import sys
import time
from pathlib import Path

from pipeline import solubility
from pipeline.common import (ROOT, STRUCTURES, TARGETS, connect, log,
                             manifest_read, manifest_write, tuning)


def _patch_one(job: tuple[str, str, str]) -> dict:
    from pipeline.patches import patches

    cid, target, path = job
    try:
        result = patches(Path(path))
    except Exception as exc:
        return {"id": cid, "target": target,
                "status": f"failed:patches:{type(exc).__name__}"}
    if not result:
        return {"id": cid, "target": target, "status": "failed:patches:empty"}
    return {"id": cid, "target": target, "status": "ok", **result}


def fold_records() -> dict[str, dict]:
    """The latest manifest row per candidate from the folding stage."""
    out: dict[str, dict] = {}
    for rec in manifest_read("fold"):
        if "id" in rec:
            out[str(rec["id"])] = rec
    return out


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--target", action="append", default=None)
    parser.add_argument("--workers", type=int, default=None)
    parser.add_argument("--skip-patches", action="store_true")
    args = parser.parse_args()

    targets = args.target or TARGETS
    workers = args.workers or max(2, min(int(tuning().get("cpu_workers", 8)), 12))
    conn = connect()
    folds = fold_records()

    # ------------------------------------------------------------- patches
    patch_results: dict[str, dict] = {}
    if not args.skip_patches:
        done = {str(r["id"]): r for r in manifest_read("patches")
                if r.get("status") == "ok"}
        jobs: list[tuple[str, str, str]] = []
        for target in targets:
            for row in conn.execute(
                    "SELECT id FROM candidates WHERE target = ? ORDER BY id",
                    (target,)):
                cid = row["id"]
                if cid in done:
                    patch_results[cid] = done[cid]
                    continue
                path = STRUCTURES / target / f"{cid}.pdb.gz"
                if path.exists():
                    jobs.append((cid, target, str(path)))
        if jobs:
            log(f"tier1: surface patches for {len(jobs)} models on {workers} workers")
            started = time.time()
            ctx = mp.get_context("spawn")
            with ctx.Pool(workers) as pool:
                for i, rec in enumerate(pool.imap_unordered(_patch_one, jobs,
                                                            chunksize=8), 1):
                    manifest_write("patches", rec)
                    if rec.get("status") == "ok":
                        patch_results[str(rec["id"])] = rec
                    if i % 250 == 0 or i == len(jobs):
                        rate = i / max(time.time() - started, 1e-9)
                        log(f"tier1: patches {i}/{len(jobs)}, {rate * 60:.0f} per minute")
        else:
            log(f"tier1: surface patches already computed for "
                f"{len(patch_results)} models")

    # ------------------------------------------- optional external predictors
    tnp = _load_side_table("tnp")
    nativeness = _load_side_table("nativeness")
    humanness = _load_side_table("humanness")
    neighbours = _load_side_table("neighbours")
    for name, table in (("TNP", tnp), ("nativeness", nativeness),
                        ("humanness", humanness), ("neighbours", neighbours)):
        log(f"tier1: {name} has values for {len(table)} candidates"
            + ("" if table else " (not run; the column stays pending)"))

    # --------------------------------------------------------------- write
    written = 0
    for target in targets:
        rows = conn.execute(
            "SELECT id, sequence FROM candidates WHERE target = ? ORDER BY id",
            (target,)).fetchall()
        records = []
        for row in rows:
            cid = row["id"]
            fold = folds.get(cid, {})
            patch = patch_results.get(cid, {})
            sol = solubility.both(row["sequence"])
            t = tnp.get(cid, {})
            n = nativeness.get(cid, {})
            h = humanness.get(cid, {})
            nb = neighbours.get(cid, {})
            path = STRUCTURES / target / f"{cid}.pdb.gz"
            status = fold.get("status", "not_run")
            # A handful of candidates carry a framework-1 N-terminal heptad
            # ANARCI cannot number, so those seven residues were modelled as
            # the VHH germline to get coordinates at all. The sequence axes
            # still use the real sequence; say so in the drawer rather than
            # leaving it in a manifest nobody reads.
            germline = int(fold.get("n_term_germline") or 0)
            note = None
            if germline:
                note = (f"The first {germline} framework residues were modelled "
                        f"as the VHH germline: this candidate's own N-terminal "
                        f"heptad could not be numbered, so no structure "
                        f"predictor would close the domain. Every sequence-level "
                        f"value on this page uses the real sequence.")
            elif str(status).startswith("failed:"):
                note = ("No predicted model: this candidate's J region carries "
                        "no framework-4 motif, so no numbering scheme can close "
                        "the domain. Its sequence axes are unaffected.")
            records.append({
                "candidate_id": cid,
                "path": str(path.relative_to(ROOT)) if path.exists() else None,
                "refined": fold.get("refined"),
                "mean_error": fold.get("mean_error"),
                "cdr3_error": fold.get("cdr3_error"),
                "psh": t.get("psh"),
                "ppc": t.get("ppc"),
                "pnc": t.get("pnc"),
                # SFvCSP is a variable-fragment charge-symmetry parameter between
                # a heavy and a light chain. A VHH has one domain, so the metric
                # is not applicable here and is stored NULL rather than zero.
                "sfvcsp_na": None,
                "tnp_flags": t.get("flags"),
                "tnp_version": t.get("version"),
                "tnp_percentile": t.get("percentile"),
                "swi": sol["swi"],
                "protein_sol": sol["protein_sol"],
                "aggrescan3d": None,
                "hydrophobic_patch_area": patch.get("hydrophobic_patch_area"),
                "charged_patch_area": patch.get("charged_patch_area"),
                "sasa_total": patch.get("sasa_total"),
                "nativeness_ppl": n.get("pseudo_perplexity"),
                "oasis_humanness": h.get("oasis_identity"),
                "nearest_pdb": nb.get("nearest_pdb"),
                "nearest_tm": nb.get("nearest_tm"),
                "nearest_identity_fw": nb.get("identity_fw"),
                "nearest_identity_cdr": nb.get("identity_cdr"),
                "has_experimental_badge": nb.get("badge", 0),
                "model_note": note,
                "status": status,
            })
        columns = list(records[0])
        conn.executemany(
            f"INSERT OR REPLACE INTO structures ({', '.join(columns)}) VALUES "
            f"({', '.join(':' + c for c in columns)})", records)
        conn.commit()
        written += len(records)
        with_model = sum(1 for r in records if r["path"])
        with_patch = sum(1 for r in records if r["hydrophobic_patch_area"] is not None)
        log(f"tier1: {target} wrote {len(records)} rows, {with_model} with a model, "
            f"{with_patch} with surface patches")

    failures = [r for r in manifest_read("fold")
                if str(r.get("status", "")).startswith("failed:")]
    trimmed = [r for r in manifest_read("fold") if r.get("trimmed")]
    conn.execute("INSERT OR REPLACE INTO qc (key, value, passed) VALUES (?,?,?)",
                 ("tier1_models",
                  f"{written} rows; {len(list(STRUCTURES.rglob('*.pdb.gz')))} models; "
                  f"{len(failures)} folds failed; {len(trimmed)} folded as the "
                  f"numbered domain after an unnumberable N-terminal extension",
                  1 if len(failures) < 0.02 * max(written, 1) else 0))
    conn.commit()
    conn.close()
    log(f"tier1: complete, {written} structure rows")
    return 0


def _load_side_table(stage: str) -> dict[str, dict]:
    """Optional predictor results, keyed by candidate id. Absent is fine."""
    out: dict[str, dict] = {}
    for rec in manifest_read(stage):
        if rec.get("id") and rec.get("status") == "ok":
            out[str(rec["id"])] = rec
    return out


if __name__ == "__main__":
    sys.exit(main())
