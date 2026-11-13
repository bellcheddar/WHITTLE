"""Ranking: within-target percentiles, Pareto fronts, desirability. Spec 6.5.

Three rules this module exists to enforce:

1. **Percentiles are computed within each target and never pooled** (spec 6.5
   item 1 and Appendix B item 6). Five independent scales. There is no global
   league table anywhere in this codebase.
2. **Pareto non-dominated sorting comes first, a single score second** (6.5
   items 2 and 3). The front is what the app shows; Derringer-Suich desirability
   is offered on top for anyone who wants one ordering, with the weights set in
   the UI.
3. **Double counting is measured, not assumed away** (6.5 item 4). Hydrophobic
   patch area feeds non-specificity, aggregation and expression, so the axis
   correlation matrix is computed per target and published on the QC panel. Any
   pair above 0.8 is stated beside the ranking.

Each axis is a weighted mean of the z-scores of whichever of its component
metrics are present, then converted to a percentile. An axis with no component
available is NULL, and `meta.json` records which tier each axis currently
reflects, so the app can say "pending" rather than imply a number exists.
"""
from __future__ import annotations

import json
import math
import sys
from dataclasses import dataclass

from pipeline.common import (AXES, STATIC_DATA, TARGETS, connect, log,
                             manifest_write, write_json)


@dataclass(frozen=True)
class Component:
    """One contributing metric. `sign` is +1 when more is better."""
    field: str
    sign: int
    tier: int
    weight: float = 1.0
    note: str = ""


# Component tables. A component is used only when the column has values; this is
# what lets Phase 1 ship two honest axes and Phase 3 ship four, with no change
# here beyond the data arriving.
AXIS_COMPONENTS: dict[str, list[Component]] = {
    "affinity": [
        Component("iptm", +1, 2, 2.0, "co-folding interface confidence, averaged over seeds"),
        Component("prodigy_dg", -1, 2, 1.0, "predicted interface binding free energy"),
        Component("bsa_total", +1, 2, 0.5, "buried surface area"),
    ],
    "specificity": [
        Component("psh", -1, 1, 2.0, "TNP patch surface hydrophobicity"),
        Component("tnp_flags", -1, 1, 1.0, "TNP metrics outside the clinical-stage range"),
        Component("n_cdr_liabilities", -1, 0, 1.0, "liability motifs inside a CDR"),
        Component("abs_charge", -1, 0, 0.5, "distance of net charge at pH 7.4 from zero"),
    ],
    "thermo": [
        Component("tm_pred", +1, 1, 3.0, "predicted apparent Tm from the WHITTLE model"),
        Component("hallmark_n", +1, 0, 1.0, "camelid hallmark residues present, 0-4"),
        Component("canonical_disulfide", +1, 0, 0.5, "canonical Cys23-Cys104 pair"),
        Component("cdr3_len", -1, 0, 0.5, "CDR3 length: long loops cost stability"),
    ],
    "expression": [
        Component("swi", +1, 1, 1.5, "Solubility-Weighted Index"),
        Component("protein_sol", +1, 1, 1.5, "Protein-Sol solubility"),
        Component("nativeness_ppl", -1, 1, 1.0, "pseudo-perplexity under an antibody LM"),
        Component("n_liabilities", -1, 0, 0.5, "liability motifs anywhere in the sequence"),
    ],
}

# What the Pareto front and the desirability function run over. Nativeness is a
# fifth explorer axis in DESIGN.md but is already a component of expression, so
# it is plotted and not ranked: including it twice is the double counting spec
# 6.5 item 4 warns about.
RANK_AXES = list(AXES)


# ----------------------------------------------------------------------- statistics

def zscores(values: list[float | None]) -> list[float | None]:
    present = [v for v in values if v is not None]
    if len(present) < 2:
        return [None] * len(values)
    mean = sum(present) / len(present)
    var = sum((v - mean) ** 2 for v in present) / (len(present) - 1)
    sd = math.sqrt(var)
    if sd == 0:
        return [0.0 if v is not None else None for v in values]
    return [None if v is None else (v - mean) / sd for v in values]


def percentiles(values: list[float | None]) -> list[float | None]:
    """Within-target percentile, 0-100, ties averaged.

    Ties are averaged rather than ranked arbitrarily, so a column of identical
    values reads as 50 everywhere instead of as a spurious ordering.
    """
    indexed = [(v, i) for i, v in enumerate(values) if v is not None]
    out: list[float | None] = [None] * len(values)
    if not indexed:
        return out
    indexed.sort()
    n = len(indexed)
    i = 0
    while i < n:
        j = i
        while j + 1 < n and indexed[j + 1][0] == indexed[i][0]:
            j += 1
        # Average rank over the tie block, mapped onto 0-100.
        mean_rank = (i + j) / 2
        pct = 100.0 * mean_rank / (n - 1) if n > 1 else 50.0
        for k in range(i, j + 1):
            out[indexed[k][1]] = round(pct, 2)
        i = j + 1
    return out


def pearson(xs: list[float | None], ys: list[float | None]) -> float | None:
    pairs = [(x, y) for x, y in zip(xs, ys) if x is not None and y is not None]
    if len(pairs) < 3:
        return None
    n = len(pairs)
    mx = sum(p[0] for p in pairs) / n
    my = sum(p[1] for p in pairs) / n
    sxy = sum((p[0] - mx) * (p[1] - my) for p in pairs)
    sxx = sum((p[0] - mx) ** 2 for p in pairs)
    syy = sum((p[1] - my) ** 2 for p in pairs)
    if sxx == 0 or syy == 0:
        return None
    return round(sxy / math.sqrt(sxx * syy), 4)


# --------------------------------------------------------------------------- Pareto

def pareto_fronts(points: list[dict[str, float | None]], axes: list[str]
                  ) -> list[int]:
    """Non-dominated sorting, maximising every axis. Returns a 1-based rank.

    A candidate missing an axis cannot dominate on it and is not dominated on it
    either: the comparison runs over the axes both candidates have. This is the
    honest treatment while axes are still filling, and once every axis is
    present it reduces to ordinary non-dominated sorting.
    """
    n = len(points)
    if n == 0:
        return []
    usable = [a for a in axes
              if any(p.get(a) is not None for p in points)]
    if not usable:
        return [1] * n

    def dominates(a: dict, b: dict) -> bool:
        shared = [ax for ax in usable
                  if a.get(ax) is not None and b.get(ax) is not None]
        if not shared:
            return False
        strictly_better = False
        for ax in shared:
            if a[ax] < b[ax]:
                return False
            if a[ax] > b[ax]:
                strictly_better = True
        return strictly_better

    dominated_by: list[int] = [0] * n
    dominates_list: list[list[int]] = [[] for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if dominates(points[i], points[j]):
                dominates_list[i].append(j)
                dominated_by[j] += 1
            elif dominates(points[j], points[i]):
                dominates_list[j].append(i)
                dominated_by[i] += 1

    ranks = [0] * n
    current = [i for i in range(n) if dominated_by[i] == 0]
    front = 1
    counts = list(dominated_by)
    while current:
        nxt: list[int] = []
        for i in current:
            ranks[i] = front
            for j in dominates_list[i]:
                counts[j] -= 1
                if counts[j] == 0:
                    nxt.append(j)
        current = nxt
        front += 1
    # Anything left in a cycle (impossible for a strict partial order, but cheap
    # to guard) lands on the last front.
    for i in range(n):
        if ranks[i] == 0:
            ranks[i] = front
    return ranks


def desirability(values: dict[str, float | None], weights: dict[str, float]
                 ) -> float | None:
    """Derringer-Suich: the weighted geometric mean of the per-axis desirabilities.

    Each axis is already a 0-100 percentile, so d_i is the percentile on 0-1 and
    the overall desirability is prod(d_i ** w_i) ** (1 / sum w). A geometric mean
    is the point of the method: one axis at zero sinks the candidate, which a
    weighted sum would hide.
    """
    terms: list[tuple[float, float]] = []
    for axis, weight in weights.items():
        value = values.get(axis)
        if value is None or weight <= 0:
            continue
        # Clamp away from exactly zero so a single bottom-ranked axis does not
        # erase the information in the others entirely.
        d = max(value / 100.0, 0.001)
        terms.append((d, weight))
    if not terms:
        return None
    total = sum(w for _, w in terms)
    log_sum = sum(w * math.log(d) for d, w in terms)
    return round(math.exp(log_sum / total), 4)


# ------------------------------------------------------------------------ the driver

def gather(conn, target: str) -> list[dict]:
    """One row per candidate with every component metric resolved."""
    rows = conn.execute("""
        SELECT c.id, c.target, c.cdr3_len, c.hallmark_n, c.canonical_disulfide,
               c.net_charge_74, c.n_liabilities,
               (SELECT COUNT(*) FROM liabilities l
                 WHERE l.candidate_id = c.id AND l.region <> 'framework')
                 AS n_cdr_liabilities,
               s.psh, s.tnp_flags, s.swi, s.protein_sol, s.nativeness_ppl,
               t.tm_pred,
               (SELECT AVG(f.iptm) FROM cofolds f
                 WHERE f.candidate_id = c.id AND f.status = 'ok') AS iptm,
               (SELECT AVG(f.prodigy_dg) FROM cofolds f
                 WHERE f.candidate_id = c.id AND f.status = 'ok') AS prodigy_dg,
               (SELECT AVG(f.bsa_total) FROM cofolds f
                 WHERE f.candidate_id = c.id AND f.status = 'ok') AS bsa_total
        FROM candidates c
        LEFT JOIN structures s ON s.candidate_id = c.id
        LEFT JOIN thermostability t ON t.candidate_id = c.id
        WHERE c.target = ? ORDER BY c.id""", (target,)).fetchall()
    out = []
    for row in rows:
        rec = dict(row)
        charge = rec.get("net_charge_74")
        rec["abs_charge"] = None if charge is None else abs(charge)
        out.append(rec)
    return out


# A row scored from only a light subset of an axis's components is not on the
# same footing as one scored from all of them: the weighted mean of fewer
# standardised terms has the same expected value but a wider spread, so those
# rows drift to both extremes and look more decisive than they are. A row that
# is missing components carrying more than half the axis weight gets no value
# for that axis at all, which the app shows as pending.
MIN_WEIGHT_FRACTION = 0.5


def axis_values(records: list[dict], axis: str) -> tuple[list[float | None], list[Component]]:
    """The axis as a weighted mean of available component z-scores.

    A row is scored only where the components it does have carry at least
    `MIN_WEIGHT_FRACTION` of the weight the axis could use: see the note above.
    """
    used: list[Component] = []
    zs: list[list[float | None]] = []
    for comp in AXIS_COMPONENTS[axis]:
        raw = [r.get(comp.field) for r in records]
        if not any(v is not None for v in raw):
            continue
        z = zscores([None if v is None else float(v) * comp.sign for v in raw])
        if all(v is None for v in z):
            continue
        used.append(comp)
        zs.append(z)
    if not used:
        return [None] * len(records), []
    available_weight = sum(c.weight for c in used)
    floor = MIN_WEIGHT_FRACTION * available_weight
    out: list[float | None] = []
    for i in range(len(records)):
        num = 0.0
        den = 0.0
        for comp, z in zip(used, zs):
            if z[i] is None:
                continue
            num += comp.weight * z[i]
            den += comp.weight
        out.append(num / den if den >= floor and den > 0 else None)
    return out, used


def score_external(conn, target: str, external: list[dict]) -> dict[str, dict]:
    """Place rows that are not candidates against a target's candidate scales.

    Used for the control set (spec 8.2). The question a control answers is
    "where would this have ranked among the candidates for this target", so:

      * every component z-score uses the **candidates'** mean and standard
        deviation, not a mean recomputed over candidates plus controls, which
        would let 50 shuffled sequences move the scale they are being measured
        against;
      * the percentile is the position among the candidate values, by counting,
        so a control can legitimately land below the 0th or above the 100th
        candidate and is clamped to the range rather than extrapolated.

    Returns {id: {axis: percentile}}.
    """
    records = gather(conn, target)
    if not records:
        return {}
    out: dict[str, dict] = {row["id"]: {} for row in external}
    for axis in RANK_AXES:
        used: list[Component] = []
        stats: list[tuple[Component, float, float]] = []
        for comp in AXIS_COMPONENTS[axis]:
            raw = [r.get(comp.field) for r in records]
            present = [float(v) * comp.sign for v in raw if v is not None]
            if len(present) < 2:
                continue
            mean = sum(present) / len(present)
            var = sum((v - mean) ** 2 for v in present) / (len(present) - 1)
            sd = math.sqrt(var)
            if sd == 0:
                continue
            # The external row must have this component too, or the two sides
            # are not being compared on the same axis.
            if not any(row.get(comp.field) is not None for row in external):
                continue
            used.append(comp)
            stats.append((comp, mean, sd))
        if not used:
            continue

        floor = MIN_WEIGHT_FRACTION * sum(c.weight for c in used)

        def combine(source: dict) -> float | None:
            num = den = 0.0
            for comp, mean, sd in stats:
                value = source.get(comp.field)
                if value is None:
                    continue
                num += comp.weight * ((float(value) * comp.sign) - mean) / sd
                den += comp.weight
            # The same rule the candidates are held to, so a control is never
            # scored on terms a candidate would have been refused.
            return num / den if den >= floor and den > 0 else None

        candidate_values = [combine(r) for r in records]
        reference = sorted(v for v in candidate_values if v is not None)
        if len(reference) < 10:
            continue
        for row in external:
            value = combine(row)
            if value is None:
                continue
            below = sum(1 for v in reference if v < value)
            equal = sum(1 for v in reference if v == value)
            rank = below + (equal - 1) / 2.0 if equal else below
            pct = 100.0 * rank / (len(reference) - 1)
            out[row["id"]][axis] = round(min(max(pct, 0.0), 100.0), 2)
    return out


def main() -> int:
    conn = connect()
    weights = {axis: 1.0 for axis in RANK_AXES}  # the UI re-weights client side
    axis_tiers: dict[str, dict] = {}
    correlations: dict[str, dict] = {}
    conn.execute("DELETE FROM scores")

    for target in TARGETS:
        records = gather(conn, target)
        if not records:
            continue
        pct: dict[str, list[float | None]] = {}
        for axis in RANK_AXES:
            raw, used = axis_values(records, axis)
            pct[axis] = percentiles(raw)
            tier = max((c.tier for c in used), default=None)
            axis_tiers.setdefault(axis, {
                "tier": tier,
                "components": [{"field": c.field, "sign": c.sign, "tier": c.tier,
                                "weight": c.weight, "note": c.note} for c in used],
                "available": bool(used),
            })
        # Nativeness is plotted as the fifth explorer axis but is not ranked.
        native_raw = [None if r.get("nativeness_ppl") is None
                      else -float(r["nativeness_ppl"]) for r in records]
        pct["nativeness"] = percentiles(native_raw)

        withheld = {axis: sum(1 for v in pct[axis] if v is None)
                    for axis in RANK_AXES
                    if axis_tiers.get(axis, {}).get("available")}
        withheld = {k: v for k, v in withheld.items() if v}
        if withheld:
            log(f"score: {target} withheld an axis value from rows missing more "
                f"than half the axis weight: {withheld} of {len(records)}")

        points = [{axis: pct[axis][i] for axis in RANK_AXES}
                  for i in range(len(records))]
        ranks = pareto_fronts(points, RANK_AXES)
        desir = [desirability(p, weights) for p in points]

        # Rank within target: by desirability where there is one, Pareto rank as
        # the tie-break, then id so the order is stable across rebuilds.
        order = sorted(range(len(records)),
                       key=lambda i: (-(desir[i] if desir[i] is not None else -1),
                                      ranks[i], records[i]["id"]))
        rank_of = {i: position + 1 for position, i in enumerate(order)}

        affinity_provisional = 0 if axis_tiers.get("affinity", {}).get("tier") == 2 else 1
        conn.executemany(
            "INSERT OR REPLACE INTO scores (candidate_id, target, axis_affinity, "
            "axis_specificity, axis_thermo, axis_expression, axis_nativeness, "
            "pareto_rank, on_front, desirability, tier_reached, "
            "axis_affinity_provisional) VALUES (?,?,?,?,?,?,?,?,?,?,?,?)",
            [(records[i]["id"], target, pct["affinity"][i], pct["specificity"][i],
              pct["thermo"][i], pct["expression"][i], pct["nativeness"][i],
              ranks[i], 1 if ranks[i] == 1 else 0, desir[i],
              max((t for t in (axis_tiers[a]["tier"] for a in RANK_AXES)
                   if t is not None), default=0),
              affinity_provisional)
             for i in range(len(records))])
        # `rank` is the display ordering position within the target, kept in its
        # own table so `scores` stays exactly the shape spec Section 7 sets out.
        conn.executemany("INSERT OR REPLACE INTO ranks (candidate_id, target, rank) "
                         "VALUES (?,?,?)",
                         [(records[i]["id"], target, rank_of[i])
                          for i in range(len(records))])

        # Axis correlation matrix, per target (spec 6.5 item 4).
        matrix: dict[str, dict[str, float | None]] = {}
        flagged: list[str] = []
        for a in RANK_AXES:
            matrix[a] = {}
            for b in RANK_AXES:
                r = 1.0 if a == b else pearson(pct[a], pct[b])
                matrix[a][b] = r
                if a < b and r is not None and abs(r) > 0.8:
                    flagged.append(f"{a} and {b} correlate at r = {r}")
        correlations[target] = {"matrix": matrix, "flagged": flagged,
                                "n": len(records)}
        n_front = sum(1 for r in ranks if r == 1)
        conn.commit()
        manifest_write("score", {"id": target, "rows": len(records),
                                 "on_front": n_front, "status": "ok"})
        log(f"score: {target} {len(records)} rows, {n_front} on the Pareto front, "
            f"axes available "
            f"{[a for a in RANK_AXES if axis_tiers[a]['available']]}"
            + (f"; CORRELATION WARNING {flagged}" if flagged else ""))

    conn.execute("INSERT OR REPLACE INTO qc (key, value, passed) VALUES (?,?,?)",
                 ("axis_correlations", json.dumps(correlations),
                  0 if any(c["flagged"] for c in correlations.values()) else 1))
    conn.execute("INSERT OR REPLACE INTO qc (key, value, passed) VALUES (?,?,?)",
                 ("axis_availability", json.dumps(axis_tiers), 1))
    conn.commit()
    write_json(STATIC_DATA / "axes.json",
               {"axes": axis_tiers, "correlations": correlations,
                "rank_axes": RANK_AXES})
    conn.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
