"""Scoring invariants. Spec 6.5 and Section 10 items 2, 3 and 6.

These are the properties that must hold whatever the data looks like. They run
against synthetic inputs, so they do not need a finished pipeline and they fail
for a reason rather than because a tier has not run yet.

The one they exist for above all: **a percentile is computed within its own
target and a ranking is never pooled across targets.** That is the single
decision the whole project rests on, and it is the kind of thing that breaks
silently when a query loses a WHERE clause.
"""
from __future__ import annotations

import random

import pytest

from pipeline.score import (desirability, pareto_fronts, pearson, percentiles,
                            zscores)


# ------------------------------------------------------------------ percentiles

def test_percentile_range_and_monotonicity() -> None:
    values = [3.0, 1.0, 2.0, 5.0, 4.0]
    out = percentiles(values)
    assert min(out) == 0.0 and max(out) == 100.0
    # The ordering of the percentiles must match the ordering of the values.
    assert [p for _v, p in sorted(zip(values, out))] == sorted(out)


def test_percentile_ties_are_averaged_not_ordered_arbitrarily() -> None:
    out = percentiles([7.0] * 10)
    assert set(out) == {50.0}, "a column of identical values must read as 50"


def test_percentile_ignores_missing_values() -> None:
    out = percentiles([1.0, None, 3.0, None, 2.0])
    assert out[1] is None and out[3] is None
    present = [p for p in out if p is not None]
    assert min(present) == 0.0 and max(present) == 100.0


def test_percentile_is_scale_invariant() -> None:
    """A percentile must not change when the underlying units do.

    This is what makes two targets' percentiles comparable *as ranks* even
    though their raw values are not comparable at all.
    """
    base = [random.Random(1).uniform(0, 1) for _ in range(50)]
    assert percentiles(base) == percentiles([v * 1000 + 42 for v in base])


def test_percentiles_are_computed_within_a_population_not_across_two() -> None:
    """Two populations on different scales must each span 0 to 100.

    If the two were pooled, the lower-scaled population would occupy the bottom
    of a single range and never reach 100. That is exactly the pooled
    cross-target table spec 6.5 item 5 forbids.
    """
    cold = [10.0, 20.0, 30.0, 40.0]
    hot = [1000.0, 2000.0, 3000.0, 4000.0]
    assert percentiles(cold) == percentiles(hot) == [0.0, 33.33, 66.67, 100.0]
    pooled = percentiles(cold + hot)
    assert max(pooled[:4]) < 50, ("pooling the two populations would be the "
                                  "cross-target league table the spec forbids")


# ---------------------------------------------------------------------- Pareto

def test_pareto_front_contains_the_dominating_point() -> None:
    points = [
        {"a": 100.0, "b": 100.0},   # dominates everything
        {"a": 50.0, "b": 50.0},
        {"a": 10.0, "b": 10.0},
    ]
    ranks = pareto_fronts(points, ["a", "b"])
    assert ranks[0] == 1
    assert ranks[1] > 1 and ranks[2] > ranks[1]


def test_pareto_front_keeps_genuine_trade_offs() -> None:
    """A trade-off is not domination. Both extremes belong on the front.

    This is the whole reason the spec asks for a Pareto front rather than a
    weighted sum: a long stabilised CDR3 traded against polyreactivity is a
    real choice, and a single score hides it.
    """
    points = [
        {"a": 100.0, "b": 0.0},
        {"a": 0.0, "b": 100.0},
        {"a": 50.0, "b": 50.0},
        {"a": 10.0, "b": 10.0},
    ]
    ranks = pareto_fronts(points, ["a", "b"])
    assert ranks[0] == ranks[1] == ranks[2] == 1
    assert ranks[3] > 1


def test_pareto_agrees_with_a_brute_force_implementation() -> None:
    """Spec 10 item 6, as a property rather than a one-off check."""
    rng = random.Random(20261113)
    for _ in range(25):
        points = [{"a": rng.uniform(0, 100), "b": rng.uniform(0, 100),
                   "c": rng.uniform(0, 100)} for _ in range(40)]
        ranks = pareto_fronts(points, ["a", "b", "c"])
        fast = {i for i, r in enumerate(ranks) if r == 1}
        brute = set()
        for i, a in enumerate(points):
            if not any(
                    all(b[ax] >= a[ax] for ax in "abc")
                    and any(b[ax] > a[ax] for ax in "abc")
                    for j, b in enumerate(points) if i != j):
                brute.add(i)
        assert fast == brute


def test_pareto_tolerates_an_axis_that_has_not_been_computed() -> None:
    points = [{"a": 90.0, "b": None}, {"a": 10.0, "b": None}]
    ranks = pareto_fronts(points, ["a", "b"])
    assert ranks[0] == 1 and ranks[1] > 1


def test_pareto_with_no_axes_at_all_puts_everything_on_one_front() -> None:
    points = [{"a": None}, {"a": None}]
    assert pareto_fronts(points, ["a"]) == [1, 1]


# ---------------------------------------------------------------- desirability

def test_desirability_is_a_geometric_mean_so_one_bad_axis_hurts() -> None:
    """A weighted sum would let a strong axis mask a failing one.

    Derringer-Suich is chosen precisely so it cannot: a candidate at the bottom
    of one axis must score below one that is merely average everywhere.
    """
    weights = {"a": 1.0, "b": 1.0}
    lopsided = desirability({"a": 100.0, "b": 0.5}, weights)
    balanced = desirability({"a": 55.0, "b": 50.0}, weights)
    assert lopsided < balanced


def test_desirability_respects_weights() -> None:
    values = {"a": 90.0, "b": 10.0}
    assert (desirability(values, {"a": 9.0, "b": 1.0})
            > desirability(values, {"a": 1.0, "b": 9.0}))


def test_desirability_is_none_when_nothing_is_computed() -> None:
    assert desirability({"a": None}, {"a": 1.0}) is None


# ----------------------------------------------------------------- statistics

def test_zscores_are_centred_and_scaled() -> None:
    out = [v for v in zscores([1.0, 2.0, 3.0, 4.0, 5.0]) if v is not None]
    assert abs(sum(out)) < 1e-9
    assert out[0] < 0 < out[-1]


def test_zscores_of_a_constant_column_do_not_divide_by_zero() -> None:
    assert zscores([4.0] * 5) == [0.0] * 5


def test_pearson_detects_a_restatement() -> None:
    xs = [float(i) for i in range(30)]
    assert pearson(xs, [2 * x + 1 for x in xs]) == pytest.approx(1.0, abs=1e-6)
    assert pearson(xs, [-x for x in xs]) == pytest.approx(-1.0, abs=1e-6)


def test_pearson_is_none_without_enough_paired_values() -> None:
    assert pearson([1.0, 2.0], [1.0, 2.0]) is None


def test_desirability_refuses_a_row_missing_a_required_axis() -> None:
    """A candidate scored on fewer axes must not out-rank a complete one.

    This is a regression test for a real ranking fault: two candidates reached
    rank 1 for their target on two axes each while carrying no specificity
    value, because a geometric mean over fewer terms is systematically higher.
    One of them was not even on the Pareto front.
    """
    weights = {"a": 1.0, "b": 1.0, "c": 1.0}
    required = ["a", "b", "c"]
    partial = {"a": 98.0, "b": 95.0, "c": None}
    complete = {"a": 80.0, "b": 80.0, "c": 80.0}
    assert desirability(partial, weights, required) is None
    assert desirability(complete, weights, required) is not None
    # Without the requirement the partial row scores higher, which is the bug.
    assert desirability(partial, weights) > desirability(complete, weights)


def test_desirability_ignores_an_axis_nobody_has() -> None:
    """An axis pending for the whole run must not disqualify the whole table."""
    weights = {"a": 1.0, "b": 1.0}
    # `required` is derived from the axes actually computed, so an all-null
    # axis is simply not in it.
    assert desirability({"a": 70.0, "b": 60.0}, weights, ["a", "b"]) is not None


def test_affinity_is_displayed_and_not_ranked_on():
    """Affinity has a percentile and no say in the ordering (DECISIONS 50).

    The guard is against a change that quietly lets Tier 2 back into the
    ranking, which would move the 24 that Marc's hard requirement fixes.
    """
    from pipeline.score import RANK_AXES, ordering_axes
    ordering = ordering_axes()
    assert "affinity" in RANK_AXES, "affinity must still be computed and shown"
    assert "affinity" not in ordering, (
        "affinity must not order the table while run.tier2_reorders_shortlist "
        "is false: no component of it beats its own shuffled control")
    for axis in ("specificity", "thermo", "expression"):
        assert axis in ordering, f"{axis} is Tier 1 and must keep ordering"


def test_an_annotate_only_axis_cannot_remove_a_candidate_from_the_front():
    """A row dominated on affinity alone stays on the Pareto front."""
    from pipeline.score import ordering_axes, pareto_fronts
    ordering = ordering_axes()
    # Identical on every ordering axis, far apart on affinity.
    points = [{"affinity": 99.0, "specificity": 50.0, "thermo": 50.0,
               "expression": 50.0},
              {"affinity": 1.0, "specificity": 50.0, "thermo": 50.0,
               "expression": 50.0}]
    ranks = pareto_fronts(points, ordering)
    assert ranks == [1, 1], (
        "affinity is annotation, so it must not dominate a candidate off the "
        f"front; got {ranks}")


def test_app_meta_fallback_covers_every_template_key():
    """The pre-build fallback must provide every `meta.X` a template reads.

    `app.meta()` substitutes a stub when the built bundle is absent, so the
    site renders a shell rather than a 500 before the first build-db run. That
    only works while the stub has every key the templates touch, and it
    drifted: `axis_verdicts` went into meta.json and into candidates.html
    without going in here, and `{{ Undefined | tojson }}` raises, so the page
    the fallback exists to protect would not have rendered at all.
    """
    import ast
    import re
    from pathlib import Path

    root = Path(__file__).resolve().parent.parent.parent
    source = (root / "app/app.py").read_text()
    # The fallback dict literal inside meta(), parsed rather than imported, so
    # the test needs neither Flask nor the built bundles.
    tree = ast.parse(source)
    provided: set[str] = set()
    for node in ast.walk(tree):
        if (isinstance(node, ast.FunctionDef) and node.name == "meta"):
            for inner in ast.walk(node):
                if isinstance(inner, ast.Dict):
                    for key in inner.keys:
                        if isinstance(key, ast.Constant) and isinstance(key.value, str):
                            provided.add(key.value)
    assert "run_id" in provided, "did not find the fallback dict in app.meta()"

    used: set[str] = set()
    for template in (root / "app/templates").glob("*.html"):
        used |= set(re.findall(r"meta\.([a-z_][a-z0-9_]*)", template.read_text()))

    missing = sorted(used - provided)
    assert not missing, (
        f"templates read meta keys the pre-build fallback does not provide: "
        f"{missing}. Add them to app.meta() or the shell 500s before the "
        f"first build-db run.")
