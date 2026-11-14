"""Bake the static JSON bundles the app loads. Spec Section 7, last paragraph.

The droplet does no request-time computation, so everything the front page needs
is precomputed here on the Studio:

    app/static/data/meta.json              run metadata, targets, the 48 columns
    app/static/data/candidates_<T>.json    one bundle per target
    app/static/data/controls.json          the control set, Tab 2
    app/static/data/about.json             references, provenance, limitations
    app/static/data/qc.json                the QC panel
    app/static/data/tm_<T>.json            the Tm distribution the drawer plots

The only thing the app reads from SQLite at request time is a single candidate
row for the drawer, which is a read of a read-only file and not a computation.

    pixi run build-db
"""
from __future__ import annotations

import datetime as dt
import json
import sys

from pipeline import columns as cols
from pipeline.common import (STATIC_DATA, TARGETS, connect, git_hash_or_blank, log,
                             run_id, targets_config, write_json)

# Which liability kinds roll up into which table counter.
LIABILITY_ROLLUP = {
    "n_glyc": ["n_glycosylation"],
    "n_deamidation": ["deamidation_ng", "deamidation_ns", "deamidation_nn"],
    "n_isomerisation": ["isomerisation_dg", "isomerisation_dp", "isomerisation_ds",
                        "fragmentation_dq"],
    "n_oxidation": ["oxidation_met", "oxidation_trp"],
    "n_unpaired_cys": ["unpaired_cysteine"],
}


def target_rows(conn) -> list[dict]:
    out = []
    cfg = targets_config()["targets"]
    for row in conn.execute("SELECT * FROM targets ORDER BY target"):
        rec = dict(row)
        spec = cfg.get(rec["target"], {})
        start, end = rec["range_start"], rec["range_end"]
        stoich = rec["stoichiometry"] or 1
        if rec["construct"] == "soluble_trimer":
            construct_label = f"soluble trimer {start}-{end} x{stoich}"
        elif rec["construct"] == "ecd_only":
            construct_label = f"ECD only {start}-{end}"
        else:
            construct_label = f"ECD {start}-{end}"
        rec["construct_label"] = construct_label
        rec["antigen_length"] = (rec["chain_length"] or 0) * stoich
        rec["complex_length"] = rec["antigen_length"] + 125
        rec["notes_list"] = spec.get("notes", [])
        out.append(rec)
    return out


def liability_counts(conn) -> dict[str, dict[str, int]]:
    counts: dict[str, dict[str, int]] = {}
    for row in conn.execute(
            "SELECT candidate_id, kind, region, COUNT(*) AS n FROM liabilities "
            "GROUP BY candidate_id, kind, region"):
        rec = counts.setdefault(row["candidate_id"], {})
        for bucket, kinds in LIABILITY_ROLLUP.items():
            if row["kind"] in kinds:
                rec[bucket] = rec.get(bucket, 0) + row["n"]
        if row["region"] != "framework":
            rec["n_cdr"] = rec.get("n_cdr", 0) + row["n"]
    return counts


def candidate_bundle(conn, target: str, liabs: dict[str, dict[str, int]]) -> dict:
    rows = conn.execute("""
        SELECT c.*, r.rank,
               s.axis_affinity, s.axis_specificity, s.axis_thermo,
               s.axis_expression, s.axis_nativeness, s.pareto_rank, s.on_front,
               s.desirability, s.tier_reached, s.axis_affinity_provisional,
               st.psh, st.ppc, st.pnc, st.tnp_flags, st.swi, st.protein_sol,
               st.nativeness_ppl, st.oasis_humanness, st.nearest_pdb,
               st.nearest_tm, st.nearest_identity_fw, st.nearest_identity_cdr,
               st.has_experimental_badge, st.aggrescan3d, st.mean_error,
               st.cdr3_error, st.refined, st.path AS structure_path,
               t.tm_pred, t.tm_err, t.temstapro, t.thermompnn_ddg,
               t.train_nn_distance,
               (SELECT AVG(f.iptm) FROM cofolds f
                 WHERE f.candidate_id = c.id AND f.status='ok') AS iptm,
               (SELECT AVG(f.ptm) FROM cofolds f
                 WHERE f.candidate_id = c.id AND f.status='ok') AS ptm,
               (SELECT AVG(f.plddt) FROM cofolds f
                 WHERE f.candidate_id = c.id AND f.status='ok') AS plddt,
               (SELECT AVG(f.prodigy_dg) FROM cofolds f
                 WHERE f.candidate_id = c.id AND f.status='ok') AS prodigy_dg,
               (SELECT AVG(f.bsa_total) FROM cofolds f
                 WHERE f.candidate_id = c.id AND f.status='ok') AS bsa_total,
               (SELECT MAX(f.n_subunits_contacted) FROM cofolds f
                 WHERE f.candidate_id = c.id AND f.status='ok') AS n_subunits_contacted,
               (SELECT MIN(f.epitope_bin) FROM cofolds f
                 WHERE f.candidate_id = c.id AND f.status='ok') AS epitope_bin,
               -- Interface quality. Averaged over the seeds like the rest, so a
               -- lucky single sample cannot carry a candidate.
               (SELECT MAX(f.pose_consensus) FROM cofolds f
                 WHERE f.candidate_id = c.id AND f.status='ok') AS pose_consensus,
               (SELECT AVG(f.sc) FROM cofolds f
                 WHERE f.candidate_id = c.id AND f.status='ok') AS sc,
               (SELECT AVG(f.cdr_interface_fraction) FROM cofolds f
                 WHERE f.candidate_id = c.id AND f.status='ok')
                 AS cdr_interface_fraction,
               (SELECT AVG(f.n_hbond) FROM cofolds f
                 WHERE f.candidate_id = c.id AND f.status='ok') AS n_hbond,
               (SELECT AVG(f.hbond_per_100a2) FROM cofolds f
                 WHERE f.candidate_id = c.id AND f.status='ok') AS hbond_per_100a2,
               (SELECT AVG(f.polar_fraction) FROM cofolds f
                 WHERE f.candidate_id = c.id AND f.status='ok') AS polar_fraction,
               (SELECT AVG(f.n_his_contacts) FROM cofolds f
                 WHERE f.candidate_id = c.id AND f.status='ok') AS n_his_contacts
        FROM candidates c
        LEFT JOIN ranks r ON r.candidate_id = c.id
        LEFT JOIN scores s ON s.candidate_id = c.id
        LEFT JOIN structures st ON st.candidate_id = c.id
        LEFT JOIN thermostability t ON t.candidate_id = c.id
        WHERE c.target = ? ORDER BY c.id""", (target,)).fetchall()

    out = []
    for row in rows:
        rec = dict(row)
        counters = liabs.get(rec["id"], {})
        flat = {
            "id": rec["id"],
            "target": rec["target"],
            "rank": rec["rank"],
            "length": rec["length"],
            "v_gene": rec["v_gene"],
            "cdr1_len": rec["cdr1_len"],
            "cdr2_len": rec["cdr2_len"],
            "cdr3_len": rec["cdr3_len"],
            "cdr3_cluster": rec["cdr3_cluster"],
            "multi_target": rec["multi_target"],
            "hallmark_n": rec["hallmark_n"],
            "net_charge_74": rec["net_charge_74"],
            "pi": rec["pi"],
            "liabilities": rec["n_liabilities"],
            "n_cdr_liabilities": counters.get("n_cdr", 0),
            "n_glyc": counters.get("n_glyc", 0),
            "n_deamidation": counters.get("n_deamidation", 0),
            "n_isomerisation": counters.get("n_isomerisation", 0),
            "n_oxidation": counters.get("n_oxidation", 0),
            "n_unpaired_cys": counters.get("n_unpaired_cys", 0),
            "on_front": rec["on_front"] or 0,
            "pareto_rank": rec["pareto_rank"],
            "desirability": rec["desirability"],
            "tier_reached": rec["tier_reached"],
            # Carried for the toolbar's CDR3 motif search (DESIGN 3.2), not a column.
            "_cdr3": rec["cdr3"],
        }
        for key in ("axis_affinity", "axis_specificity", "axis_thermo",
                    "axis_expression", "axis_nativeness", "iptm", "ptm",
                    "prodigy_dg", "bsa_total", "n_subunits_contacted",
                    "epitope_bin", "psh", "ppc", "pnc", "tnp_flags", "swi",
                    "protein_sol", "nativeness_ppl", "oasis_humanness",
                    "nearest_pdb", "tm_pred", "tm_err", "temstapro",
                    "upstream_tm_pred",
                    "pose_consensus", "sc", "cdr_interface_fraction",
                    "n_hbond", "hbond_per_100a2", "polar_fraction",
                    "n_his_contacts"):
            flat[key] = rec.get(key)
        flat["experimental_badge"] = rec.get("has_experimental_badge") or 0
        out.append(flat)

    present = {key for key in cols.BY_KEY
               if key not in ("axis_profile",)
               and any(r.get(key) is not None for r in out)}
    return {
        "target": target,
        "n": len(out),
        "rows": out,
        "columns_with_data": sorted(present),
        "columns_pending": sorted(k for k in cols.BY_KEY
                                  if k not in present and k != "axis_profile"),
    }


def tm_distribution(conn, target: str) -> dict:
    """What the drawer's Tm panel plots: the whole target's distribution."""
    rows = conn.execute(
        "SELECT t.tm_pred FROM thermostability t JOIN candidates c "
        "ON c.id = t.candidate_id WHERE c.target = ? AND t.tm_pred IS NOT NULL",
        (target,)).fetchall()
    values = sorted(float(r[0]) for r in rows)
    if not values:
        return {"target": target, "n": 0, "bins": [], "counts": [],
                "min": None, "max": None, "mean": None}
    lo, hi = values[0], values[-1]
    n_bins = 28
    width = (hi - lo) / n_bins if hi > lo else 1.0
    counts = [0] * n_bins
    for v in values:
        idx = min(int((v - lo) / width), n_bins - 1)
        counts[idx] += 1
    return {
        "target": target, "n": len(values),
        "bins": [round(lo + width * (i + 0.5), 2) for i in range(n_bins)],
        "counts": counts, "width": round(width, 3),
        "min": round(lo, 2), "max": round(hi, 2),
        "mean": round(sum(values) / len(values), 2),
    }


def meta_bundle(conn, targets: list[dict]) -> dict:
    run_cfg = targets_config().get("run", {}) or {}
    n_total = conn.execute("SELECT COUNT(*) FROM candidates").fetchone()[0]
    axes = {}
    row = conn.execute("SELECT value FROM qc WHERE key='axis_availability'").fetchone()
    if row:
        axes = json.loads(row[0])
    tier_counts = {
        "tier0": conn.execute(
            "SELECT COUNT(*) FROM candidates WHERE cdr3 IS NOT NULL").fetchone()[0],
        "tier1": conn.execute("SELECT COUNT(*) FROM structures").fetchone()[0],
        "tier1b": conn.execute("SELECT COUNT(*) FROM thermostability").fetchone()[0],
        "tier2": conn.execute(
            "SELECT COUNT(DISTINCT candidate_id) FROM cofolds").fetchone()[0],
    }
    provisional = [name for name, info in axes.items() if not info.get("available")]

    # Per-axis, per-target verdict from the vacuous tests. `vacuous_tests.py`
    # says a non-discriminating axis "must be labelled as such wherever it
    # appears", and until now it appeared only on the Control set tab: a
    # reader on the Candidates page saw an Expression percentile with nothing
    # to say it does not separate a real candidate from a shuffled one.
    verdicts: dict[str, dict[str, str]] = {}
    # The pooled block beside the per-target labels. Since DECISIONS 53 the
    # verdict is decided once per axis over all five targets, so the five
    # labels agree by construction and the number worth showing is the pooled
    # count. Without it the panel would state a verdict and hide the evidence.
    pooled: dict[str, dict] = {}
    row = conn.execute(
        "SELECT value FROM qc WHERE key='vacuous_tests'").fetchone()
    if row:
        try:
            payload = json.loads(row[0])
        except (TypeError, ValueError):
            payload = {}
        for axis, block in (payload or {}).items():
            if block.get("pooled"):
                keep = ("verdict", "strength", "pooled_wins", "pooled_n",
                        "pooled_fraction", "pooled_p", "per_target", "finding",
                        "testable_weight_fraction")
                pooled[axis] = {k: block["pooled"][k] for k in keep
                                if k in block["pooled"]}
            for target, result in (block.get("targets") or {}).items():
                # `verdict` is set by the test itself where it has something
                # more specific to say than pass or fail. "unconfirmed" and
                # "untestable" are both distinct from "pending": pending means
                # the check has not run, while these two mean it ran and could
                # not rule. Collapsing them into pending told a reader nothing
                # had been measured when something had.
                passed = result.get("passed")
                verdicts.setdefault(axis, {})[target] = result.get("verdict") or (
                    "pending" if passed is None else
                    "informative" if passed else "vacuous")
    return {
        "run_id": run_id(),
        "built_at": dt.datetime.now().strftime("%Y-%m-%d %H:%M"),
        "built_date": dt.datetime.now().strftime("%d %b %Y"),
        "git_sha": git_hash_or_blank(),
        "n_candidates": n_total,
        "n_targets": len(targets),
        "targets": targets,
        "tier_counts": tier_counts,
        "columns": cols.registry(),
        "axes": axes,
        "axes_pending": provisional,
        "axis_verdicts": verdicts,
        "axis_verdict_pooled": pooled,
        # The About page used to hard-code "250 deliberately broken
        # sequences" and "Fifty sequences per target", both of which went
        # stale the moment the control set was scaled. Statistics a reader is
        # asked to trust must come from the database that produced them.
        "n_scrambled_controls": conn.execute(
            "SELECT COUNT(*) FROM controls WHERE kind = 'scrambled'"
        ).fetchone()[0],
        "scrambled_per_target": int(targets_config()["run"].get(
            "scrambled_per_target", 50)),
        "ranking_scope": "within_target",
        "upstream_caveat": "not comparable between targets",
        # The hero draws the shortlist cut and the gates solve for it, so the
        # number has to reach the browser rather than being repeated in JS.
        "shortlist_per_target": run_cfg.get("shortlist_per_target", 24),
        "tier2_top_n_per_target": run_cfg.get("tier2_top_n_per_target"),
    }


def qc_bundle(conn) -> dict:
    rows = [dict(r) for r in conn.execute(
        "SELECT key, value, passed FROM qc ORDER BY key")]
    parsed = []
    for row in rows:
        value = row["value"]
        if value and value.startswith(("{", "[")):
            try:
                value = json.loads(value)
            except json.JSONDecodeError:
                pass
        parsed.append({"key": row["key"], "value": value,
                       "passed": None if row["passed"] is None
                       else bool(row["passed"]),
                       "state": ("pending" if row["passed"] is None
                                 else "pass" if row["passed"] else "attention")})
    return {"checks": parsed,
            "built_at": dt.datetime.now().strftime("%Y-%m-%d %H:%M")}


def controls_bundle(conn) -> dict:
    rows = [dict(r) for r in conn.execute("SELECT * FROM controls ORDER BY kind, id")]
    return {"n": len(rows), "rows": rows,
            "coverage_caveat":
                "Only EGFR and TNFa have target-matched published nanobodies, so this "
                "tab checks the general axes for all five targets and the "
                "target-aware axis for those two only. TL1A, PRLR and GLP1R have no "
                "published nanobody benchmark."}


def about_bundle(conn) -> dict:
    refs = [dict(r) for r in conn.execute(
        'SELECT * FROM "references" ORDER BY sort_order, key')]
    prov = {r["key"]: r["value"] for r in conn.execute(
        "SELECT key, value FROM provenance ORDER BY key")}
    limitations = targets_config().get("limitations_statement", [])
    upstream_defects = [
        "The upstream composite score is not comparable between targets: EGFR and "
        "TNFa blend four terms including s_his, the other three blend three. It is "
        "kept as upstream_score and the whole table is never sorted by it.",
        "The upstream Tm column spans 71 to 101 degrees C. An apparent Tm above 100 "
        "degrees C is not a value a DSF run gives a VHH in solution, so the column is "
        "treated as uncalibrated and shown only beside this build's own prediction.",
        "cdr_novelty and cdr_identity are empty in all five source files and are "
        "ingested as null. This build computes its own equivalents.",
        "target_cdr_identity exists only for EGFR and TNFa, the two targets with "
        "known nanobodies in SAbDab-nano. It is left null elsewhere rather than "
        "fabricated, and the UI says there is no target-matched reference.",
        "Ten sequences appear under more than one target. They are not deduplicated: "
        "they are flagged and used as internal promiscuity controls.",
        "The upstream hallmark column implements a three-of-four rule over the "
        "camelid hallmark tetrad, not all four. This build reports the strict "
        "all-four call and the graded 0 to 4 count side by side.",
        "The upstream CDR histidine counts use a wider CDR definition than IMGT, so "
        "they agree with the recomputed counts for about 91 percent of rows.",
    ]
    return {"references": refs, "provenance": prov, "limitations": limitations,
            "upstream_defects": upstream_defects}


def main() -> int:
    conn = connect(read_only=False)
    targets = target_rows(conn)
    liabs = liability_counts(conn)

    write_json(STATIC_DATA / "meta.json", meta_bundle(conn, targets))
    for target in TARGETS:
        bundle = candidate_bundle(conn, target, liabs)
        write_json(STATIC_DATA / f"candidates_{target}.json", bundle)
        write_json(STATIC_DATA / f"tm_{target}.json", tm_distribution(conn, target))
        log(f"build-db: {target} bundle, {bundle['n']} rows, "
            f"{len(bundle['columns_with_data'])} columns with data, "
            f"{len(bundle['columns_pending'])} pending")
    write_json(STATIC_DATA / "controls.json", controls_bundle(conn))
    write_json(STATIC_DATA / "about.json", about_bundle(conn))
    write_json(STATIC_DATA / "qc.json", qc_bundle(conn))
    conn.close()

    total = sum(path.stat().st_size for path in STATIC_DATA.glob("*.json"))
    log(f"build-db: {len(list(STATIC_DATA.glob('*.json')))} bundles, "
        f"{total / 1024:.0f} kB total")
    return 0


if __name__ == "__main__":
    sys.exit(main())
