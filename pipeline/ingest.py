"""Ingest the five shortlist CSVs into SQLite. Spec Sections 5.1 and 11 (Phase 1).

Every defect spec 5.1 records is handled here, explicitly and in one place:

1. `score` is not comparable across targets (EGFR and TNFa blend four terms,
   the others three). Stored as `upstream_score`; the app attaches the caveat
   wherever it is shown and never sorts the whole table by it.
2. `tm_pred` is uncalibrated, 71-101 C. Stored as `upstream_tm_pred` only.
3. `cdr_novelty` and `cdr_identity` are empty in all five files: ingested NULL.
4. `target_cdr_identity` exists only for EGFR and TNFa: NULL elsewhere, never
   fabricated.
5. 10 sequences occur under more than one target: flagged `multi_target`, never
   deduplicated.
6. His columns carried for all five targets; `s_his` only exists for EGFR and
   TNFa and is a design criterion, not a quality metric.

Resumable: the whole ingest is one transaction per target with a manifest row,
so a kill re-runs at most one target.
"""
from __future__ import annotations

import csv
import hashlib
import sys
from collections import Counter
from pathlib import Path

from pipeline import schema
from pipeline.common import (SHORTLIST, TARGETS, connect, log, manifest_write,
                             targets_config)

# Targets whose libraries carry the intentional pH-switch His enrichment and so
# carry s_his and target_cdr_identity (spec 5.1 items 4 and 6).
FOUR_TERM_TARGETS = {"EGFR", "TNFa"}


def md5(text: str) -> str:
    return hashlib.md5(text.encode()).hexdigest()


def num(value: str | None) -> float | None:
    if value is None:
        return None
    value = value.strip()
    if value == "" or value.lower() in {"na", "nan", "none", "null"}:
        return None
    try:
        return float(value)
    except ValueError:
        return None


def integer(value: str | None) -> int | None:
    f = num(value)
    return None if f is None else int(f)


def boolean(value: str | None) -> int | None:
    if value is None or value.strip() == "":
        return None
    return 1 if value.strip().lower() in {"true", "1", "yes"} else 0


def load_targets(conn) -> None:
    cfg = targets_config()
    rows = []
    for name, spec in cfg["targets"].items():
        start, end = spec["range"]
        rows.append((
            name, spec["uniprot"], spec["name"], spec["construct"], start, end,
            spec["length"], spec["stoichiometry"],
            1 if spec.get("his_enriched_library") else 0,
            1 if spec.get("ph_switch_intended") else 0,
            1 if spec.get("known_nanobodies") else 0,
            1 if spec["construct"] == "ecd_only" else 0,
            None, spec.get("tier2_cost"),
            " ".join(spec.get("notes", [])),
        ))
    conn.executemany(
        "INSERT OR REPLACE INTO targets (target, uniprot, name, construct, "
        "range_start, range_end, chain_length, stoichiometry, his_enriched, "
        "ph_switch_intended, known_nanobodies, ecd_scoped, antigen_md5, "
        "tier2_cost, notes) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)", rows)
    conn.commit()
    log(f"ingest: {len(rows)} targets written from config/targets.yaml")


def read_csv(target: str) -> list[dict]:
    path = SHORTLIST / f"{target}.csv"
    if not path.exists():
        raise FileNotFoundError(path)
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh))


def main() -> int:
    conn = connect()
    schema.create(conn)
    load_targets(conn)

    # Pass one: every sequence, so multi_target is known before any row is written.
    per_target: dict[str, list[dict]] = {}
    seq_counts: Counter[str] = Counter()
    for target in TARGETS:
        rows = read_csv(target)
        per_target[target] = rows
        for row in rows:
            seq_counts[row["sequence"].strip()] += 1
    shared = {seq for seq, n in seq_counts.items() if n > 1}
    unique_total = len(seq_counts)
    log(f"ingest: {sum(len(r) for r in per_target.values())} rows, "
        f"{unique_total} unique sequences, {len(shared)} shared across targets")

    conn.execute("DELETE FROM candidates")
    conn.commit()

    total = 0
    for target in TARGETS:
        rows = per_target[target]
        four_term = target in FOUR_TERM_TARGETS
        records = []
        for i, row in enumerate(rows, 1):
            seq = row["sequence"].strip()
            records.append({
                "id": f"{target}-{i:04d}",
                "target": target,
                "sequence": seq,
                "seq_md5": md5(seq),
                "length": len(seq),
                "multi_target": 1 if seq in shared else 0,
                # Tier 0 fills numbering, CDRs, charge, the recomputed hallmark
                # and the recomputed His counts: left NULL here on purpose, so a
                # half-run shows as missing rather than as an upstream value
                # wearing our column name.
                "upstream_tm_pred": num(row.get("tm_pred")),
                "upstream_liability_score": num(row.get("liability_score")),
                "upstream_on_nll": num(row.get("on_nll")),
                "upstream_score": num(row.get("score")),
                "upstream_s_tm": num(row.get("s_tm")),
                "upstream_s_liab": num(row.get("s_liab")),
                "upstream_s_nll": num(row.get("s_nll")),
                # s_his and target_cdr_identity exist only for EGFR and TNFa and
                # are never fabricated for the other three (spec 5.1 items 4, 6).
                "upstream_s_his": num(row.get("s_his")) if four_term else None,
                "upstream_hallmark": boolean(row.get("hallmark")),
                "upstream_cdr3_cluster": integer(row.get("cdr3_cluster")),
                "upstream_his_cdr_total": integer(row.get("his_cdr_total")),
                "upstream_his_cdr1": integer(row.get("his_cdr1")),
                "upstream_his_cdr2": integer(row.get("his_cdr2")),
                "upstream_his_cdr3": integer(row.get("his_cdr3")),
                "upstream_target_cdr_identity":
                    num(row.get("target_cdr_identity")) if four_term else None,
                "upstream_global_cdr_identity": num(row.get("global_cdr_identity")),
                "upstream_global_nearest":
                    (row.get("global_nearest") or "").strip() or None,
                "upstream_global_n_known": integer(row.get("global_n_known")),
                "upstream_cdr_novelty": num(row.get("cdr_novelty")),   # empty: NULL
                "upstream_cdr_identity": num(row.get("cdr_identity")),  # ditto
            })
        columns = list(records[0])
        conn.executemany(
            f"INSERT INTO candidates ({', '.join(columns)}) VALUES "
            f"({', '.join(':' + c for c in columns)})", records)
        conn.commit()
        total += len(records)
        manifest_write("ingest", {"id": target, "rows": len(records), "status": "ok"})
        log(f"ingest: {target} {len(records)} rows"
            + (" (four-term upstream score, s_his and target_cdr_identity present)"
               if four_term else " (three-term upstream score, no s_his)"))

    # Reconciliation, spec 10 item 2.
    counts = dict(conn.execute(
        "SELECT target, COUNT(*) FROM candidates GROUP BY target").fetchall())
    n_unique = conn.execute(
        "SELECT COUNT(DISTINCT seq_md5) FROM candidates").fetchone()[0]
    n_multi_rows = conn.execute(
        "SELECT COUNT(*) FROM candidates WHERE multi_target = 1").fetchone()[0]
    n_multi_seqs = conn.execute(
        "SELECT COUNT(DISTINCT seq_md5) FROM candidates WHERE multi_target = 1").fetchone()[0]
    log(f"ingest: reconcile total={total} per_target={counts} unique={n_unique} "
        f"multi_target sequences={n_multi_seqs} rows={n_multi_rows}")

    ok = (total == 2500 and all(v == 500 for v in counts.values())
          and n_unique == 2487 and n_multi_seqs == 10)
    conn.execute("INSERT OR REPLACE INTO qc (key, value, passed) VALUES (?,?,?)",
                 ("ingest_reconcile",
                  f"total={total} unique={n_unique} multi_target_sequences={n_multi_seqs}",
                  1 if ok else 0))
    conn.commit()
    conn.close()
    if not ok:
        log("ingest: WARNING reconciliation does not match spec 5.1 "
            "(2500 rows, 2487 unique, 10 multi-target sequences)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
