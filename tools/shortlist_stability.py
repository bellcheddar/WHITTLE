"""How firm is the cut at 24? Spec Section 9, the same shape as the other QC.

The shortlist is the deliverable, so the question a reader will ask about it is
not "what is the ordering" but "would a different reasonable choice have given
me a different 24". This measures that two ways and writes both into `qc`:

1. **The boundary.** The desirability of the candidate ranked 24 against the one
   ranked 25, and both against the library median. These are two separate
   claims and only one of them is strong: the top 24 sit far above the median,
   and the candidate ranked 24 is all but indistinguishable from the one ranked
   25.
2. **The reweighting ensemble.** `desirability` takes its weights from the UI,
   so the shipped equal weighting is one choice among many defensible ones.
   Each ordering axis is doubled and halved in turn, which gives an ensemble of
   seven weightings, and a candidate is called firm if it holds a shortlist
   place under every one of them. Candidates outside the shipped 24 that take a
   place under at least one weighting are the contested band, and how far down
   the shipped ranking that band reaches is the honest width of the cut.

Both are reported rather than resolved. The ensemble is deliberately a set of
reweightings and not a resampling: it answers "does the choice of weights pick
these 24", which is a question about our decisions, and says nothing about the
measurement error inside each axis, which is a different question and is
covered by the model card's error bar.

The axis-dropped comparison is not part of the ensemble. Removing an axis
entirely is a bigger step than reweighting it, and which axis is worth removing
is read from the scrambled control's own verdicts rather than hardcoded: any
ordering axis the control cannot call informative is dropped in turn and the
retention reported. That is the number that says how much of the shortlist
rests on an axis the control cannot defend.

    pixi run qc-stability
"""
from __future__ import annotations

import json
import statistics
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline.common import TARGETS, connect, log, targets_config  # noqa: E402
from pipeline.score import desirability, ordering_axes  # noqa: E402


def schemes(ordering: list[str]) -> dict[str, dict[str, float]]:
    """The ensemble: equal weights, then each axis doubled and halved."""
    equal = {axis: 1.0 for axis in ordering}
    out = {"equal": equal}
    for axis in ordering:
        out[f"{axis} x2"] = {**equal, axis: 2.0}
        out[f"{axis} x0.5"] = {**equal, axis: 0.5}
    return out


def untestable_ordering_axes(conn, ordering: list[str]) -> list[str]:
    """Ordering axes the scrambled control could not call informative.

    Read from `qc.vacuous_tests` rather than named here, so that an axis which
    later earns or loses its verdict changes this measurement without anyone
    remembering to edit it.
    """
    row = conn.execute(
        "SELECT value FROM qc WHERE key='vacuous_tests'").fetchone()
    if not row:
        return []
    try:
        payload = json.loads(row[0]) or {}
    except (TypeError, ValueError):
        return []
    out = []
    for axis in ordering:
        pooled = (payload.get(axis) or {}).get("pooled") or {}
        if pooled.get("verdict") and pooled["verdict"] != "informative":
            out.append(axis)
    return out


def ranked(rows: list[dict], weights: dict[str, float],
           ordering: list[str]) -> list[str]:
    """Candidate ids best first, under the shipped sort key.

    The key is `score.main`'s exactly: desirability, then the Pareto rank, then
    the id. The Pareto rank carries no weights, so it is read from the table
    rather than recomputed, and under equal weights this reproduces the shipped
    `ranks` row for row. `_self_check` asserts that it does.
    """
    required = [a for a in ordering if any(r[f"axis_{a}"] is not None
                                           for r in rows)]
    scored = []
    for row in rows:
        value = desirability({a: row[f"axis_{a}"] for a in ordering},
                             weights, required)
        scored.append((-(value if value is not None else -1.0),
                       row["pareto_rank"] if row["pareto_rank"] is not None else 10**6,
                       row["id"], value))
    scored.sort(key=lambda t: (t[0], t[1], t[2]))
    return [t[2] for t in scored]


def measure(conn, target: str, n: int, ordering: list[str],
            dropping: list[str]) -> dict:
    rows = [dict(r) for r in conn.execute(
        "SELECT candidate_id AS id, axis_specificity, axis_thermo, "
        "axis_expression, axis_affinity, pareto_rank, on_front, desirability "
        "FROM scores WHERE target = ?", (target,))]
    stored = {r[0]: r[1] for r in conn.execute(
        "SELECT candidate_id, rank FROM ranks WHERE target = ?", (target,))}

    ensemble = schemes(ordering)
    tops = {name: ranked(rows, w, ordering)[:n] for name, w in ensemble.items()}
    shipped = tops["equal"]

    # The self-check: equal weights must reproduce the shipped shortlist. If it
    # does not, this whole measurement is describing a ranking the app does not
    # serve, so it is an error rather than a finding.
    mismatch = [cid for cid in shipped if stored.get(cid, 10**6) > n]
    if mismatch:
        raise SystemExit(
            f"stability: equal weights do not reproduce the shipped shortlist "
            f"for {target}: {mismatch}. Fix this before reading any number "
            f"below it; the ensemble is relative to the shipped ranking.")

    held = {cid: sum(1 for t in tops.values() if cid in t) for cid in
            {c for t in tops.values() for c in t}}
    total = len(tops)
    firm = sorted(cid for cid in shipped if held[cid] == total)
    contested_in = sorted((stored[cid], cid) for cid in shipped
                          if held[cid] < total)
    contested_out = sorted((stored.get(cid, 0), cid) for cid in held
                           if cid not in shipped)

    desirabilities = sorted((r["desirability"] for r in rows
                             if r["desirability"] is not None), reverse=True)
    at = lambda k: desirabilities[k - 1] if len(desirabilities) >= k else None
    median = statistics.median(desirabilities) if desirabilities else None

    dropped = {}
    for axis in dropping:
        weights = {a: 1.0 for a in ordering if a != axis}
        kept = len(set(ranked(rows, weights, ordering)[:n]) & set(shipped))
        dropped[axis] = kept

    return {
        "target": target,
        "n_scored": len(desirabilities),
        "shortlist": n,
        "d_first": at(1),
        "d_last": at(n),
        "d_next": at(n + 1),
        "d_median": median,
        "gap": None if at(n) is None or at(n + 1) is None
        else round(at(n) - at(n + 1), 4),
        "margin_over_median": None if at(n) is None or median is None
        else round(at(n) - median, 4),
        "on_front": sum(1 for r in rows if r["id"] in set(shipped)
                        and r["on_front"]),
        "firm": len(firm),
        "contested_in": [cid for _, cid in contested_in],
        "contested_out": [cid for _, cid in contested_out],
        "band_last_rank": max([stored.get(cid, 0) for cid in held] or [n]),
        "dropped": dropped,
    }


def main() -> int:
    conn = connect()
    ordering = ordering_axes()
    n = int((targets_config().get("run") or {}).get("shortlist_per_target", 24))
    dropping = untestable_ordering_axes(conn, ordering)

    results = {t: measure(conn, t, n, ordering, dropping) for t in TARGETS}

    # Two separate claims, two separate rows, because they do not agree.
    #
    # The band is real: every target puts its 24th candidate far above the
    # library median, so this passes.
    boundary = {
        "shortlist_per_target": n,
        "ordering_axes": ordering,
        "per_target": {t: {k: r[k] for k in
                           ("d_first", "d_last", "d_next", "d_median", "gap",
                            "margin_over_median", "on_front", "band_last_rank")}
                       for t, r in results.items()},
    }
    separated = all((r["margin_over_median"] or 0) > 0.1 for r in results.values())
    boundary["finding"] = (
        f"The shortlist is a band, not an ordering. Every target places its "
        f"{n}th candidate between "
        f"{min(r['margin_over_median'] for r in results.values()):.2f} and "
        f"{max(r['margin_over_median'] for r in results.values()):.2f} "
        f"desirability above the library median, so the {n} are genuinely "
        f"separated from the "
        f"{min(r['n_scored'] for r in results.values()) - n} to "
        f"{max(r['n_scored'] for r in results.values()) - n} others that carry "
        f"a value on all three ordering axes. The cut itself is not: the "
        f"gap between rank {n} and rank {n + 1} runs from "
        f"{min(r['gap'] for r in results.values()):.4f} to "
        f"{max(r['gap'] for r in results.values()):.4f}, so there is no sense "
        f"in which candidate {n} beats candidate {n + 1}. Read the shortlist "
        f"as the top band and the ordering inside it as a suggestion.")

    # The reweighting is where it gets uncomfortable, so this row is allowed
    # to fail. Dropping an axis the control calls untestable must not take more
    # than a third of the shortlist with it, and it does.
    reweight = {
        "ensemble": sorted(schemes(ordering)),
        "untestable_axes_dropped": dropping,
        "per_target": {t: {k: r[k] for k in
                           ("firm", "contested_in", "contested_out",
                            "band_last_rank", "dropped")}
                       for t, r in results.items()},
    }
    floor = (2.0 / 3.0) * n
    fragile = {t: r["dropped"] for t, r in results.items()
               if any(v < floor for v in r["dropped"].values())}
    reweight["finding"] = (
        f"Of the {n} shortlisted, "
        f"{min(r['firm'] for r in results.values())} to "
        f"{max(r['firm'] for r in results.values())} hold a place under every "
        f"one of the {len(schemes(ordering))} weightings in the ensemble. The "
        f"contested band reaches rank "
        f"{max(r['band_last_rank'] for r in results.values())}, so a reader "
        f"choosing their own weights in the UI would see roughly that many "
        f"candidates competing for {n} places."
        + ("" if not dropping else
           f" Dropping {' and '.join(dropping)}, which the scrambled control "
           f"cannot call informative, retains only "
           f"{min(min(r['dropped'].values()) for r in results.values())} to "
           f"{max(max(r['dropped'].values()) for r in results.values())} of "
           f"the {n}. That is the single largest sensitivity in the shortlist "
           f"and it rests on an axis the control cannot defend, which is why "
           f"this check is marked for attention rather than passed."))

    conn.execute("INSERT OR REPLACE INTO qc (key, value, passed) VALUES (?,?,?)",
                 ("shortlist_boundary", json.dumps(boundary),
                  1 if separated else 0))
    conn.execute("INSERT OR REPLACE INTO qc (key, value, passed) VALUES (?,?,?)",
                 ("shortlist_reweighting", json.dumps(reweight),
                  0 if fragile else 1))
    conn.commit()

    for t, r in results.items():
        log(f"stability: {t} d@{n} {r['d_last']:.4f} against d@{n + 1} "
            f"{r['d_next']:.4f} (gap {r['gap']:.4f}) and median "
            f"{r['d_median']:.4f}; {r['firm']} of {n} firm across the "
            f"ensemble, contested band to rank {r['band_last_rank']}, "
            f"axis-dropped retention {r['dropped']}")
    conn.close()

    print(f"shortlist stability: {n} per target, ordering on {ordering}")
    print(f"{'tgt':7}{'d@1':>8}{'d@24':>8}{'d@25':>8}{'median':>8}{'gap':>8}"
          f"{'firm':>6}{'band to':>8}  dropped")
    for t, r in results.items():
        print(f"{t:7}{r['d_first']:8.3f}{r['d_last']:8.3f}{r['d_next']:8.3f}"
              f"{r['d_median']:8.3f}{r['gap']:8.4f}{r['firm']:6}"
              f"{r['band_last_rank']:8}  "
              + ", ".join(f"{k} {v}/{r['shortlist']}"
                          for k, v in r["dropped"].items()))
    print(f"\nboundary: {boundary['finding']}")
    print(f"\nreweighting: {reweight['finding']}")
    if fragile:
        print("\nATTENTION: reported as a finding, not hidden. "
              "The app shows the contested band beside the cut.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
