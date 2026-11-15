"""The shortlist stability measurement. `tools/shortlist_stability.py`.

The whole measurement is relative to the ranking the app serves, so the one
thing that must hold is that re-deriving the order under equal weights
reproduces the shipped `ranks` table exactly. If it does not, every number the
tool prints describes a shortlist nobody is being shown, and the tool raises
rather than reporting it as a finding. These tests pin that, and pin that the
ensemble is capable of moving a shortlist place at all, so a reassuring
"everything is firm" cannot come from a comparison that never varies.
"""
from __future__ import annotations

import random
import sqlite3
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tools"))

from shortlist_stability import measure, ranked, schemes  # noqa: E402


def _scored_target(tmp_path, n=60):
    """A target with a full `scores` table and a `ranks` table built from it.

    The ranks are written with `score.main`'s own sort key so the fixture
    cannot pass by agreeing with a copy of the code under test: the key is
    desirability, then the Pareto rank, then the id.
    """
    from pipeline import schema
    from pipeline.score import desirability

    conn = sqlite3.connect(tmp_path / "s.sqlite")
    conn.row_factory = sqlite3.Row
    conn.executescript(schema.SCHEMA)
    conn.execute("INSERT INTO targets (target) VALUES ('EGFR')")
    rng = random.Random(11)
    rows = []
    for i in range(1, n + 1):
        cid = f"EGFR-{i:04d}"
        conn.execute("INSERT INTO candidates (id, target, sequence, seq_md5, "
                     "length) VALUES (?, 'EGFR', 'QVQLV', ?, 5)", (cid, cid))
        spec, thermo, expr = (rng.uniform(0, 100) for _ in range(3))
        value = desirability(
            {"specificity": spec, "thermo": thermo, "expression": expr},
            {"specificity": 1.0, "thermo": 1.0, "expression": 1.0},
            ["specificity", "thermo", "expression"])
        pareto = rng.randint(1, 4)
        conn.execute(
            "INSERT INTO scores (candidate_id, target, axis_specificity, "
            "axis_thermo, axis_expression, axis_affinity, pareto_rank, "
            "on_front, desirability) VALUES (?,'EGFR',?,?,?,NULL,?,?,?)",
            (cid, spec, thermo, expr, pareto, 1 if pareto == 1 else 0, value))
        rows.append((value, pareto, cid))
    for position, (_, _, cid) in enumerate(
            sorted(rows, key=lambda r: (-r[0], r[1], r[2])), start=1):
        conn.execute("INSERT INTO ranks (candidate_id, target, rank) "
                     "VALUES (?, 'EGFR', ?)", (cid, position))
    conn.commit()
    return conn


ORDERING = ["specificity", "thermo", "expression"]


def test_equal_weights_reproduce_the_shipped_shortlist(tmp_path) -> None:
    conn = _scored_target(tmp_path)
    result = measure(conn, "EGFR", 24, ORDERING, [])
    assert result["firm"] <= 24
    stored = {r[0] for r in conn.execute(
        "SELECT candidate_id FROM ranks WHERE target='EGFR' AND rank <= 24")}
    equal = set(ranked([dict(r) for r in conn.execute(
        "SELECT candidate_id AS id, axis_specificity, axis_thermo, "
        "axis_expression, axis_affinity, pareto_rank, on_front, desirability "
        "FROM scores WHERE target='EGFR'")],
        {a: 1.0 for a in ORDERING}, ORDERING)[:24])
    assert equal == stored, (
        "the ensemble is measured relative to the shipped ranking, so equal "
        "weights must reproduce it exactly")


def test_a_shortlist_that_disagrees_with_the_ranks_table_raises(tmp_path) -> None:
    """The guard, not a finding. A silent mismatch would invalidate the lot."""
    conn = _scored_target(tmp_path)
    # Swap one shortlisted candidate's rank with one far down the table, so the
    # stored ranking no longer matches the desirability it was built from.
    conn.execute("UPDATE ranks SET rank = 400 WHERE target='EGFR' AND rank = 3")
    conn.commit()
    with pytest.raises(SystemExit, match="do not reproduce the shipped"):
        measure(conn, "EGFR", 24, ORDERING, [])


def test_the_ensemble_can_actually_move_a_place(tmp_path) -> None:
    """Otherwise "every candidate is firm" would be vacuously true."""
    conn = _scored_target(tmp_path)
    result = measure(conn, "EGFR", 24, ORDERING, ["expression"])
    assert len(schemes(ORDERING)) == 7, "equal plus two per ordering axis"
    assert result["contested_in"] or result["contested_out"], (
        "on random percentiles some place must change hands under a doubled "
        "or halved axis; if nothing moves the ensemble is not being applied")
    assert result["band_last_rank"] >= 24
    assert result["dropped"]["expression"] <= 24


def test_the_boundary_numbers_are_the_ones_the_table_stores(tmp_path) -> None:
    conn = _scored_target(tmp_path)
    result = measure(conn, "EGFR", 24, ORDERING, [])
    ds = sorted((r[0] for r in conn.execute(
        "SELECT desirability FROM scores WHERE target='EGFR'")), reverse=True)
    assert result["d_last"] == pytest.approx(ds[23])
    assert result["d_next"] == pytest.approx(ds[24])
    assert result["gap"] == pytest.approx(round(ds[23] - ds[24], 4))
    assert result["margin_over_median"] > 0


def test_a_multi_target_control_takes_its_own_target_s_percentiles(tmp_path) -> None:
    """`controls.attach_predictions` keyed its lookup on the sequence alone.

    A multi-target control is the same sequence submitted against more than
    one target, so a key of the sequence alone collapses every one of its rows
    onto whichever candidate the query returned last. Observed on 2026-11-14:
    MULTI-06 carried PRLR-0247's affinity percentile of 89.47 under GLP1R,
    where its own candidate row GLP1R-0443 scores 30.36. Every axis is a
    within-target percentile against a different distribution, so a value
    carried across targets is meaningless, and it lands on the one control
    least able to survive it: the promiscuity control asks whether a sequence
    scores alike against unrelated antigens, and identical numbers can only
    ever answer yes.
    """
    import hashlib
    import sqlite3

    from pipeline import schema
    from pipeline.controls import attach_predictions

    conn = sqlite3.connect(tmp_path / "c.sqlite")
    conn.row_factory = sqlite3.Row
    conn.executescript(schema.SCHEMA)
    sequence = "QVQLVESGGGLVQAGGSLRLSCAAS"
    digest = hashlib.md5(sequence.encode()).hexdigest()
    for target, affinity in (("GLP1R", 30.36), ("PRLR", 89.47)):
        conn.execute("INSERT INTO targets (target) VALUES (?)", (target,))
        cid = f"{target}-0001"
        conn.execute("INSERT INTO candidates (id, target, sequence, seq_md5, "
                     "length) VALUES (?,?,?,?,?)",
                     (cid, target, sequence, digest, len(sequence)))
        conn.execute("INSERT INTO scores (candidate_id, target, axis_affinity) "
                     "VALUES (?,?,?)", (cid, target, affinity))
    conn.commit()

    rows = [{"id": "MULTI-06-GLP1R", "target": "GLP1R", "sequence": sequence},
            {"id": "MULTI-06-PRLR", "target": "PRLR", "sequence": sequence}]
    attach_predictions(conn, rows)
    got = {r["id"]: r["predicted_affinity"] for r in rows}
    assert got["MULTI-06-GLP1R"] == 30.36, (
        "a multi-target control must take the percentile of the candidate row "
        f"for its own target, not another target's; got {got}")
    assert got["MULTI-06-PRLR"] == 89.47
    assert got["MULTI-06-GLP1R"] != got["MULTI-06-PRLR"], (
        "identical values across targets is the defect this test exists for")
