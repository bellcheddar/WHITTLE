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
