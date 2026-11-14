"""What the scrambled control can and cannot test.

The negative control shuffles each CDR in place and preserves amino-acid
composition, length and the framework exactly. That is deliberate: it isolates
CDR *sequence order* as the only thing that differs. The consequence is easy
to miss and it changes how a result must be read: any metric computed from
composition alone is bit-identical between a candidate and its scrambled twin,
so it contributes nothing to the comparison, and an axis built mostly from
such metrics cannot be tested this way at all.

That is not a hypothetical. The expression axis is two thirds composition-only
by weight, and it was being reported as "non-discriminating" on all five
targets on the strength of an AUC that could not have come out any other way.
"""
from __future__ import annotations

import random

import pytest

from pipeline import solubility, tier0
from pipeline.controls import scramble_cdrs
from pipeline.score import AXIS_COMPONENTS

# A real VHH framework with three identifiable CDRs.
SEQ = ("QVQLVESGGGLVQAGGSLRLSCAASGRTFSSYAMGWFRQAPGKEREFVAAISWSGGSTYYA"
       "DSVKGRFTISRDNAKNTVYLQMNSLKPEDTAVYYCAAGSYYSGSYYSRGQGTQVTVSS")
CDRS = {"cdr1": "GRTFSSYA", "cdr2": "ISWSGGST", "cdr3": "AAGSYYSGSYYSRG"}


def _scrambled(seed: int = 7) -> str:
    return scramble_cdrs(SEQ, CDRS, random.Random(seed))


def test_the_shuffle_actually_changes_the_sequence() -> None:
    out = _scrambled()
    assert out != SEQ, "a control identical to the candidate tests nothing"
    assert len(out) == len(SEQ), "the shuffle must preserve length"
    assert sorted(out) == sorted(SEQ), "the shuffle must preserve composition"


def test_the_framework_survives_the_shuffle() -> None:
    """Only the CDRs move, so framework-derived metrics cannot respond."""
    out = _scrambled()
    assert out.startswith("QVQLVESGGGLVQAGGSLRLSCAAS")
    assert out.endswith("QGTQVTVSS"), "the J-region framework must not move"


def test_composition_only_metrics_cannot_move_under_the_shuffle() -> None:
    """The reason an axis of composition metrics cannot be tested this way.

    If this ever fails, the metric has stopped being composition-only and the
    `shuffle_invariant` flags in `pipeline/score.py` need revisiting.
    """
    out = _scrambled()
    assert solubility.both(SEQ)["swi"] == solubility.both(out)["swi"]
    assert solubility.both(SEQ)["protein_sol"] == solubility.both(out)["protein_sol"]
    assert tier0.charge_and_pi(SEQ) == tier0.charge_and_pi(out)


def test_an_axis_knows_how_much_of_it_the_control_can_reach() -> None:
    """Spec 9 item 5 is only meaningful for an axis the shuffle can move."""
    reach = {}
    for axis, components in AXIS_COMPONENTS.items():
        total = sum(c.weight for c in components)
        movable = sum(c.weight for c in components if not c.shuffle_invariant)
        reach[axis] = movable / total if total else 0.0

    # Affinity is entirely structural, so the shuffle reaches all of it.
    assert reach["affinity"] == 1.0
    # Expression is mostly composition, so it reaches a minority. An AUC near
    # chance there is the control's limit, not the axis failing.
    assert reach["expression"] < 0.5
    # Specificity is driven by patch metrics that depend on sequence order.
    assert reach["specificity"] > 0.8


def test_sign_test_is_two_sided_in_both_directions():
    """A metric that runs backwards must be reported, not scored p = 1.

    Until 2026-11-14 this took only the upper tail, so `bsa_total` beating its
    own parent on 40 of 50 pairs scored p = 1 and read as carrying no signal.
    An inverted metric is worse than an uninformative one, because its weight
    moves the ranking away from the candidates it should favour.
    """
    import sys
    sys.path.insert(0, "tools")
    from vacuous_tests import sign_test

    assert sign_test(10, 50) < 1e-4, "an inverted count must be significant"
    assert sign_test(41, 50) < 1e-4, "a winning count must be significant"
    # Symmetric: the same distance either side of the middle, the same p.
    assert sign_test(10, 50) == pytest.approx(sign_test(40, 50))
    assert sign_test(0, 50) == pytest.approx(sign_test(50, 50))
    # The middle is not evidence of anything, and a p-value is a probability.
    assert sign_test(25, 50) == pytest.approx(1.0)
    assert sign_test(5, 10) == pytest.approx(1.0)
    for wins in range(0, 51):
        p = sign_test(wins, 50)
        assert 0.0 < p <= 1.0, f"{wins}/50 gave p = {p}"
    assert sign_test(0, 0) is None


def test_raising_the_control_count_does_not_move_the_original_fifty():
    """Adding controls must extend the set, never redraw it.

    250 scrambled controls have been scored and 240 pairs measured. The draw
    is `rng.sample`, and `rng.sample(pop, 50)` is not a prefix of
    `rng.sample(pop, 200)`, so widening the original call would hand
    SCRAM-EGFR-001 a different parent. `controls.main()` carries predictions
    across a rebuild by id, so the old blind prediction would then be attached
    to a new sequence with nothing failing and nothing warning.
    """
    import random
    import sqlite3

    from pipeline import common, controls, schema

    conn = sqlite3.connect(":memory:")
    conn.row_factory = sqlite3.Row
    conn.executescript(schema.SCHEMA)
    conn.execute("INSERT INTO targets (target) VALUES (?)",
                 (common.TARGETS[0],))
    rng = random.Random(1)
    for n in range(300):
        seq = "".join(rng.choice("ACDEFGHIKLMNPQRSTVWY") for _ in range(120))
        conn.execute(
            "INSERT INTO candidates (id, target, sequence, seq_md5, length, "
            "cdr1, cdr2, cdr3) VALUES (?,?,?,?,?,?,?,?)",
            (f"X-{n:04d}", common.TARGETS[0], seq, f"md5{n}", len(seq),
             seq[25:33], seq[50:58], seq[95:110]))
    conn.commit()

    cfg = common.targets_config()
    original = cfg["run"].get("scrambled_per_target")
    try:
        cfg["run"]["scrambled_per_target"] = 50
        small = {r["id"]: r["sequence"] for r in
                 controls.build_scrambled(conn, random.Random(controls.SEED))}
        cfg["run"]["scrambled_per_target"] = 200
        large = {r["id"]: r["sequence"] for r in
                 controls.build_scrambled(conn, random.Random(controls.SEED))}
    finally:
        if original is None:
            cfg["run"].pop("scrambled_per_target", None)
        else:
            cfg["run"]["scrambled_per_target"] = original

    assert len(small) > 0 and len(large) > len(small)
    for control_id, sequence in small.items():
        assert large[control_id] == sequence, (
            f"{control_id} changed sequence when the total was raised, which "
            f"invalidates every pair already measured against it")
