"""Run the control set blind through the same pipeline. Spec 8.2.

"Known nanobodies with published Kd, Tm or expression yield, run **blind**
through the same pipeline."

Blind means the published measurement is never an input. Nothing here reads
`observed_tm`, `observed_kd` or `observed_expression`; those columns exist only
so the control tab can plot observed against predicted afterwards.

The control sequences go through exactly the code paths the candidates did:
the same Tier 0 liability pass, the same TNP call, the same solubility
formulae, the same AbLang2 nativeness and the same thermostability model. The
only difference is the last step. A control's axis percentile is computed
**against its target's candidate distribution**, using scales fitted on the
candidates alone, because the question a control answers is "where would this
have ranked among the candidates", not "where does it rank among the controls".

    pixi run python -m pipeline.controls_run
    pixi run python -m pipeline.controls_run --skip-nativeness
"""
from __future__ import annotations

import argparse
import json
import subprocess
import shutil
import sys
import tempfile
import time
from pathlib import Path

from pipeline import solubility
from pipeline import tier0
from pipeline.common import (DATA, TARGETS, connect, env_binary,
                             env_with_bin, log, manifest_read, manifest_write,
                             tuning)

CONTROL_DIR = DATA / "controls"


def control_rows(conn) -> list[dict]:
    return [dict(r) for r in conn.execute(
        "SELECT id, target, sequence, kind FROM controls "
        "WHERE sequence IS NOT NULL ORDER BY id")]


# ------------------------------------------------------------------- TNP pass

def run_tnp(rows: list[dict], ncores: int) -> dict[str, dict]:
    done = {r["id"]: r for r in manifest_read("tnp_controls")
            if r.get("status") == "ok"}
    todo = [r for r in rows if r["id"] not in done]
    if not todo:
        log(f"controls_run: TNP already done for {len(done)} control sequences")
        return done
    binary = env_binary("TNP")
    if binary is None:
        log("controls_run: TNP is not on PATH; the control specificity axis "
            "stays pending")
        return done
    CONTROL_DIR.mkdir(parents=True, exist_ok=True)
    fasta = CONTROL_DIR / "controls.fasta"
    fasta.write_text("".join(f">{r['id']}\n{r['sequence']}\n" for r in todo))
    out_dir = CONTROL_DIR / "tnp"
    out_dir.mkdir(parents=True, exist_ok=True)
    log(f"controls_run: TNP over {len(todo)} control sequences on {ncores} cores")
    started = time.time()
    try:
        subprocess.run([binary, "--file", str(fasta), "--output", str(out_dir),
                        "--ncores", str(ncores)],
                       capture_output=True, text=True, timeout=60 * 240,
                       env=env_with_bin())
    except subprocess.TimeoutExpired:
        log("controls_run: TNP timed out on the controls")
        return done
    results = out_dir / "TNP_Results_Multientry.json"
    if not results.exists():
        log("controls_run: TNP produced no output for the controls")
        return done
    payload = json.loads(results.read_text())
    for cid, rec in payload.items():
        flags = rec.get("Flags", {}) or {}
        record = {
            "id": cid, "status": "ok", "psh": rec.get("PSH"),
            "ppc": rec.get("PPC"), "pnc": rec.get("PNC"),
            "flags": sum(1 for v in flags.values() if v in ("amber", "red")),
        }
        manifest_write("tnp_controls", record)
        done[cid] = record
    log(f"controls_run: TNP finished for {len(payload)} controls in "
        f"{(time.time() - started) / 60:.1f} minutes")
    # TNP returns nothing for a sequence it cannot number, and says nothing
    # about it: the entry is simply absent from the multi-entry JSON. Left
    # unlogged that is silent attrition in a control set, which is the one
    # place it must not be silent. It runs at 11 to 13 per cent on the
    # CDR-shuffled controls, evenly across targets.
    #
    # It does not bias the paired test. A shuffle that TNP drops takes its
    # whole pair out, and a shuffle preserves its parent's length exactly, so
    # the dropout acts on a property the two members of a pair share rather
    # than on the difference between them. What it costs is n, and it narrows
    # the population the result generalises to: the sequences TNP declines
    # have slightly longer CDR3 loops (16.7 against 17.4 residues).
    dropped = [r["id"] for r in todo if r["id"] not in payload]
    if dropped:
        by_target: dict[str, int] = {}
        for cid in dropped:
            parts = cid.split("-")
            key = parts[1] if len(parts) > 2 else "other"
            by_target[key] = by_target.get(key, 0) + 1
        log(f"controls_run: TNP returned nothing for {len(dropped)} of "
            f"{len(todo)} sequences ({len(dropped) / len(todo):.0%}), which it "
            f"does without reporting a failure: {by_target}. Those controls "
            f"carry no specificity axis and drop out of the paired test with "
            f"their parents. The paired comparison stays unbiased, because a "
            f"shuffle has its parent's exact length and the dropout therefore "
            f"acts on the pair rather than inside it.")
        manifest_write("tnp_controls_dropped",
                       {"id": f"run-{int(started)}", "status": "ok",
                        "n_dropped": len(dropped), "n_submitted": len(todo),
                        "per_target": by_target, "ids": sorted(dropped)})
    return done


# ------------------------------------------------------------- nativeness pass

def run_nativeness(rows: list[dict], workers: int, threads: int) -> dict[str, dict]:
    done = {r["id"]: r for r in manifest_read("nativeness_controls")
            if r.get("status") == "ok"}
    todo = [r for r in rows if r["id"] not in done]
    if not todo:
        return done
    import multiprocessing as mp

    from pipeline.nativeness import _init_worker, _score_batch

    batch = 16
    jobs = []
    for i in range(0, len(todo), batch):
        chunk = todo[i:i + batch]
        jobs.append(([c["id"] for c in chunk], [c["sequence"] for c in chunk]))
    log(f"controls_run: AbLang2 nativeness for {len(todo)} control sequences")
    ctx = mp.get_context("spawn")
    with ctx.Pool(workers, initializer=_init_worker, initargs=(threads,)) as pool:
        for records in pool.imap_unordered(_score_batch, jobs):
            for rec in records:
                manifest_write("nativeness_controls", rec)
                if rec.get("status") == "ok":
                    done[rec["id"]] = rec
    return done


# ---------------------------------------------------------------------- main

def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--skip-tnp", action="store_true")
    parser.add_argument("--skip-nativeness", action="store_true")
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()

    cores = int(tuning().get("cpu_workers", 8))
    conn = connect()
    rows = control_rows(conn)
    if not rows:
        log("controls_run: no control rows; run pipeline.controls first")
        return 1
    log(f"controls_run: {len(rows)} control sequences, run blind. No published "
        f"measurement is read by anything below.")

    tnp = {} if args.skip_tnp else run_tnp(rows, min(cores, 12))
    native = ({} if args.skip_nativeness
              else run_nativeness(rows, args.workers, max(2, cores // args.workers)))

    # The thermostability model's own prediction for each control, written by
    # `pipeline.tm_predict --controls`. Without it the thermo axis has only
    # 2.0 of its 5.0 component weight and is withheld, which is why the
    # thermo vacuous test reported "no scrambled control scored yet" for the
    # whole run: the axis was never scored for a control, so the one check
    # that could falsify it never ran.
    control_tm = {r["id"]: r["predicted_tm"] for r in conn.execute(
        "SELECT id, predicted_tm FROM controls WHERE predicted_tm IS NOT NULL")}
    if not control_tm:
        log("controls_run: no control Tm predictions; run "
            "`python -m pipeline.tm_predict --controls` first or the thermo "
            "axis stays untestable")

    # Sequence-level metrics, identical code to Tier 0 and Tier 1.
    metrics: dict[str, dict] = {}
    for row in rows:
        seq = row["sequence"]
        from pipeline import numbering as nb

        numbered = nb.number_many([(row["id"], seq)])[row["id"]]
        regions = numbered.regions or {}
        liabs = tier0.find_liabilities(seq, regions)
        charge, _pi = tier0.charge_and_pi(seq)
        sol = solubility.both(seq)
        # The remaining thermo components, by the same functions Tier 0 uses.
        _intact, hallmark_n, _detail = nb.hallmark(numbered)
        canonical, _extra = tier0.disulfides(numbered, seq)
        metrics[row["id"]] = {
            "tm_pred": control_tm.get(row["id"]),
            "hallmark_n": hallmark_n,
            "canonical_disulfide": canonical,
            "cdr3_len": len(regions.get("cdr3") or ""),
            "n_liabilities": len(liabs),
            "n_cdr_liabilities": sum(1 for l in liabs if l["region"] != "framework"),
            "abs_charge": abs(charge),
            "swi": sol["swi"],
            "protein_sol": sol["protein_sol"],
            "psh": (tnp.get(row["id"]) or {}).get("psh"),
            "tnp_flags": (tnp.get(row["id"]) or {}).get("flags"),
            "nativeness_ppl": (native.get(row["id"]) or {}).get("pseudo_perplexity"),
        }
    log(f"controls_run: sequence metrics computed for {len(metrics)} controls")

    # Percentiles against each target's candidate distribution.
    from pipeline.score import score_external

    updates: list[dict] = []
    for target in TARGETS:
        members = [r for r in rows if r["target"] == target]
        if not members:
            continue
        scored = score_external(
            conn, target,
            [{"id": r["id"], **metrics[r["id"]]} for r in members])
        for cid, axes in scored.items():
            updates.append({
                "id": cid,
                "predicted_affinity": axes.get("affinity"),
                "predicted_specificity": axes.get("specificity"),
                "predicted_thermo": axes.get("thermo"),
                "predicted_expression": axes.get("expression"),
            })
        log(f"controls_run: {len(members)} {target} controls placed against the "
            f"{target} candidate distribution")

    if updates:
        conn.executemany(
            "UPDATE controls SET predicted_affinity = :predicted_affinity, "
            "predicted_specificity = :predicted_specificity, "
            "predicted_thermo = :predicted_thermo, "
            "predicted_expression = :predicted_expression WHERE id = :id", updates)
        conn.commit()
    manifest_write("controls_run", {"id": "controls_run", "status": "ok",
                                    "scored": len(updates), "n": len(rows)})
    log(f"controls_run: {len(updates)} control rows carry axis percentiles")
    conn.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
