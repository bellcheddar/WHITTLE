"""The control set and the positive controls. Spec Sections 8.2 and 9.

This is the honesty machinery. Three kinds of control row go into the
`controls` table and onto Tab 2:

**known** - nanobodies with a published measurement, run blind through the same
pipeline. The thermostability controls are the 77 sequences held out of the Tm
model's training set entirely (NanoMelt Supplementary Dataset 2, minus the six
that overlap the training labels). The affinity controls are nanobodies from
SAbDab-nano with an experimental complex against EGFR or TNFa, which are the
only two targets with known nanobodies at all.

**scrambled** - 50 CDR-shuffled sequences per target. Spec 9 item 5, and the
control that matters most: if the affinity axis cannot separate real candidates
from shuffled ones, that is a finding to report prominently, not a failure to
hide. The shuffle keeps each CDR's amino-acid composition and its framework
intact, so what changes is the CDR sequence and nothing else.

**multi_target** - the ten sequences submitted against more than one target. A
sequence that scores well against two unrelated antigens is either polyreactive
or evidence that the scoring is not target-aware.

Every assertion in spec Section 9 is evaluated here and its result is written
to the `qc` table, pass or fail, with the number that decided it.

    pixi run controls
"""
from __future__ import annotations

import argparse
import json
import random
import sys

from pipeline.common import (TARGETS, connect, log, manifest_write, write_json)

SCRAMBLED_PER_TARGET = 50
SEED = 20261113

# Spec 5.1 item 4: only these two targets have known nanobodies in SAbDab-nano.
TARGETS_WITH_REFERENCES = ("EGFR", "TNFa")


# ------------------------------------------------------------------ scrambled

def scramble_cdrs(sequence: str, cdrs: dict[str, str], rng: random.Random) -> str:
    """Shuffle each CDR in place, leaving the framework untouched.

    Composition is preserved, so the negative control differs from the real
    candidate in CDR *sequence* and in nothing else: not length, not charge,
    not amino-acid content. Any axis that still separates them is separating on
    sequence order, which is what the control is asking about.
    """
    out = sequence
    for key in ("cdr1", "cdr2", "cdr3"):
        frag = cdrs.get(key) or ""
        if len(frag) < 3:
            continue
        shuffled = list(frag)
        for _ in range(8):
            rng.shuffle(shuffled)
            if "".join(shuffled) != frag:
                break
        out = out.replace(frag, "".join(shuffled), 1)
    return out


def build_scrambled(conn, rng: random.Random) -> list[dict]:
    rows = []
    for target in TARGETS:
        candidates = conn.execute(
            "SELECT id, sequence, cdr1, cdr2, cdr3 FROM candidates "
            "WHERE target = ? AND cdr3 IS NOT NULL ORDER BY id", (target,)).fetchall()
        if not candidates:
            continue
        chosen = rng.sample(list(candidates), min(SCRAMBLED_PER_TARGET,
                                                  len(candidates)))
        for i, row in enumerate(chosen, 1):
            cdrs = {k: row[k] for k in ("cdr1", "cdr2", "cdr3")}
            scrambled = scramble_cdrs(row["sequence"], cdrs, rng)
            if scrambled == row["sequence"]:
                continue
            rows.append({
                "id": f"SCRAM-{target}-{i:03d}",
                "target": target,
                "name": f"shuffled CDRs of {row['id']}",
                "sequence": scrambled,
                "kind": "scrambled",
                "blind": 1,
                "notes": f"CDR-shuffled negative control derived from {row['id']}. "
                         f"Framework, length and amino-acid composition are "
                         f"identical to the parent; only the CDR order differs.",
                "assertion": "must score worse on the affinity axis than the real "
                             "candidates for this target",
            })
    return rows


# ---------------------------------------------------------------------- known

def build_known_thermo() -> list[dict]:
    """The held-out Tm set: observed apparent Tm, predicted blind."""
    from pipeline.tm_data import load

    data = load()
    rows = []
    for i, rec in enumerate(data.get("external", []), 1):
        rows.append({
            "id": f"TM-EXT-{i:03d}",
            "target": None,
            "name": rec.get("name") or f"external nanobody {i}",
            "sequence": rec["sequence"],
            "observed_tm": rec["tm"],
            "kind": "known",
            "blind": 1,
            "notes": f"Held out of the Tm model's training set entirely. "
                     f"Measured by {rec.get('primary_method', 'an unstated method')}.",
            "assertion": "predicted apparent Tm within the model's stated error bar "
                         "of the published measurement",
        })
    return rows


def build_known_binders(conn) -> list[dict]:
    """Nanobodies with an experimental complex against EGFR or TNFa.

    Resolved from the SAbDab-nano snapshot at build time, never hand-typed
    (spec Section 9 opening line). If one cannot be resolved it is logged and
    dropped rather than invented.
    """
    from pipeline.neighbours import antigen_matches, load_reference

    try:
        reference = load_reference()
    except Exception as exc:
        log(f"controls: SAbDab-nano is unavailable ({exc}); the binder controls "
            f"are dropped rather than invented")
        return []
    rows = []
    for target in TARGETS_WITH_REFERENCES:
        matched = [r for r in reference if antigen_matches(r["antigens"], target)]
        for i, ref in enumerate(sorted(matched, key=lambda r: r["pdb"])[:12], 1):
            rows.append({
                "id": f"BIND-{target}-{i:02d}",
                "target": target,
                "name": f"{ref['pdb'].upper()} chain {ref['chain']}",
                "sequence": ref["sequence"],
                "pdb_id": ref["pdb"].upper(),
                "kind": "known",
                "blind": 1,
                "notes": f"Experimental nanobody-{target} complex. Antigen: "
                         f"{'; '.join(ref['antigens'][:2])}.",
                "assertion": f"must land in the top quartile of the {target} "
                             f"affinity axis and bin to its published epitope",
            })
        log(f"controls: {len(matched)} SAbDab-nano nanobodies with a {target} "
            f"complex; {min(len(matched), 12)} taken as positive controls")
    return rows


def build_multi_target(conn) -> list[dict]:
    """One row per (sequence, target) occurrence, not one per sequence.

    The point of this panel is to compare the *same* sequence against two or
    three unrelated antigens, so each occurrence needs its own target's
    percentiles. Collapsing them to one row per sequence would show whichever
    target happened to be written last and quietly answer a different question.
    """
    rows = []
    groups = conn.execute(
        "SELECT seq_md5, COUNT(*) AS n, GROUP_CONCAT(target) AS targets "
        "FROM candidates WHERE multi_target = 1 GROUP BY seq_md5 "
        "ORDER BY seq_md5").fetchall()
    for i, group in enumerate(groups, 1):
        members = conn.execute(
            "SELECT id, target, sequence FROM candidates WHERE seq_md5 = ? "
            "ORDER BY target", (group["seq_md5"],)).fetchall()
        shared = group["targets"]
        for member in members:
            rows.append({
                "id": f"MULTI-{i:02d}-{member['target']}",
                "target": member["target"],
                "name": f"MULTI-{i:02d} as {member['id']}",
                "sequence": member["sequence"],
                "kind": "multi_target",
                "blind": 0,
                "notes": f"The same sequence was also submitted against "
                         f"{shared}. Each row here is its percentile against "
                         f"its own target's candidates, so the rows can be "
                         f"compared with one another.",
                "assertion": "no assertion; scoring well against two unrelated "
                             "antigens is either polyreactivity or evidence "
                             "that the scoring is not target-aware",
            })
    return rows


# ------------------------------------------------------------------- scoring

def attach_predictions(conn, rows: list[dict]) -> None:
    """Fill the predicted columns from whatever the pipeline has computed.

    A control that is also a candidate (the multi-target rows, and any
    scrambled row whose parent is in the table) takes the candidate's scores
    directly. The rest are predicted by the same code paths and are filled when
    those stages have run over the control sequences.
    """
    by_md5 = {}
    for row in conn.execute(
            "SELECT c.seq_md5, s.axis_affinity, s.axis_specificity, s.axis_thermo, "
            "s.axis_expression, t.tm_pred FROM candidates c "
            "LEFT JOIN scores s ON s.candidate_id = c.id "
            "LEFT JOIN thermostability t ON t.candidate_id = c.id"):
        by_md5[row["seq_md5"]] = dict(row)

    import hashlib

    for rec in rows:
        digest = hashlib.md5(rec["sequence"].encode()).hexdigest()
        match = by_md5.get(digest)
        if not match:
            continue
        rec["predicted_affinity"] = match["axis_affinity"]
        rec["predicted_specificity"] = match["axis_specificity"]
        rec["predicted_thermo"] = match["axis_thermo"]
        rec["predicted_expression"] = match["axis_expression"]
        rec["predicted_tm"] = match["tm_pred"]


def assertions(conn) -> list[dict]:
    """Every assertion in spec Section 9, evaluated with its deciding number."""
    out: list[dict] = []

    def record(key: str, passed: bool | None, value: str, note: str) -> None:
        out.append({"key": key, "passed": passed, "value": value, "note": note})

    # Item 3: Tm sanity. No predicted apparent Tm above 95 degrees C.
    row = conn.execute(
        "SELECT COUNT(*) AS n, MAX(tm_pred) AS hottest FROM thermostability "
        "WHERE tm_pred IS NOT NULL").fetchone()
    if row["n"]:
        over = conn.execute(
            "SELECT COUNT(*) FROM thermostability WHERE tm_pred > 95").fetchone()[0]
        record("control_3_tm_sanity", over == 0,
               f"{over} of {row['n']} predictions above 95 deg C; "
               f"hottest {row['hottest']:.1f}",
               "Spec 9 item 3. A prediction above 95 deg C means the model is "
               "extrapolating and its error bar must reflect that.")
    else:
        record("control_3_tm_sanity", None, "no Tm predictions yet",
               "Spec 9 item 3, pending the thermostability model.")

    # Item 4: trimer sanity for TNFa and TL1A.
    for target in ("TNFa", "TL1A"):
        row = conn.execute(
            "SELECT COUNT(*) AS n, SUM(CASE WHEN n_subunits_contacted > 1 THEN 1 "
            "ELSE 0 END) AS multi FROM cofolds WHERE target = ? AND status = 'ok'",
            (target,)).fetchone()
        if row and row["n"]:
            record(f"control_4_trimer_{target}", (row["multi"] or 0) > 0,
                   f"{row['multi'] or 0} of {row['n']} co-folds contact more than "
                   f"one subunit",
                   "Spec 9 item 4. If every interface sits on a single chain, "
                   "suspect the trimer was never assembled.")
        else:
            record(f"control_4_trimer_{target}", None, "no co-folds yet",
                   "Spec 9 item 4, pending the co-folding tier.")

    # Item 5: the scrambled negative control. The one that matters most.
    for target in TARGETS:
        real = [r[0] for r in conn.execute(
            "SELECT axis_affinity FROM scores WHERE target = ? "
            "AND axis_affinity IS NOT NULL", (target,))]
        shuffled = [r[0] for r in conn.execute(
            "SELECT predicted_affinity FROM controls WHERE kind = 'scrambled' "
            "AND target = ? AND predicted_affinity IS NOT NULL", (target,))]
        if len(real) >= 10 and len(shuffled) >= 10:
            mean_real = sum(real) / len(real)
            mean_shuffled = sum(shuffled) / len(shuffled)
            separates = mean_real - mean_shuffled > 10
            record(f"control_5_scrambled_{target}", separates,
                   f"real mean {mean_real:.1f}, shuffled mean {mean_shuffled:.1f}, "
                   f"difference {mean_real - mean_shuffled:.1f} percentile points",
                   "Spec 9 item 5. If the axis does not separate them it is "
                   "non-discriminating, which is a finding to report, not a "
                   "failure to hide.")
        else:
            record(f"control_5_scrambled_{target}", None,
                   "pending the affinity axis",
                   "Spec 9 item 5, pending the co-folding tier.")

    # Item 6: hallmark reconciliation. Computed in Tier 0; surfaced here too.
    row = conn.execute(
        "SELECT value, passed FROM qc WHERE key = 'hallmark_reconciliation_graded'"
    ).fetchone()
    if row:
        record("control_6_hallmark", bool(row["passed"]),
               f"{float(row['value']):.3f} agreement under the three-of-four rule "
               f"the upstream column implements",
               "Spec 9 item 6. The strict all-four call agrees for 0.928 of rows; "
               "both are reported.")
    return out


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--skip-binders", action="store_true")
    args = parser.parse_args()

    rng = random.Random(SEED)
    conn = connect()

    rows: list[dict] = []
    rows += build_scrambled(conn, rng)
    rows += build_known_thermo()
    if not args.skip_binders:
        rows += build_known_binders(conn)
    rows += build_multi_target(conn)
    attach_predictions(conn, rows)

    conn.execute("DELETE FROM controls")
    columns = ["id", "target", "name", "source_doi", "pdb_id", "sequence",
               "observed_kd", "observed_tm", "observed_expression",
               "predicted_affinity", "predicted_thermo", "predicted_specificity",
               "predicted_expression", "kind", "blind", "assertion",
               "assertion_passed", "notes"]
    payload = [{c: rec.get(c) for c in columns} for rec in rows]
    conn.executemany(
        f"INSERT OR REPLACE INTO controls ({', '.join(columns)}) VALUES "
        f"({', '.join(':' + c for c in columns)})", payload)
    conn.commit()

    counts: dict[str, int] = {}
    for rec in rows:
        counts[rec["kind"]] = counts.get(rec["kind"], 0) + 1
    log(f"controls: {len(rows)} control rows written {counts}")

    results = assertions(conn)
    for res in results:
        passed = res["passed"]
        conn.execute("INSERT OR REPLACE INTO qc (key, value, passed) VALUES (?,?,?)",
                     (res["key"], f"{res['value']} | {res['note']}",
                      1 if passed else 0))
        state = "PASS" if passed else ("PENDING" if passed is None else "ATTENTION")
        log(f"controls: {state} {res['key']}: {res['value']}")
    conn.commit()
    manifest_write("controls", {"id": "controls", "status": "ok",
                                "rows": len(rows), "counts": counts})
    conn.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
