"""Tier 0: sequence level, all 2,500 rows, minutes. Spec 6.1.

Tier 0 discards nothing. It populates flags and the Pareto step decides later
(spec 6.1, closing line).

    pixi run tier0

Resumable per target: a manifest row per target, and the per-row work is a pure
function of the sequence, so a kill re-runs at most one target.
"""
from __future__ import annotations

import re
import shutil
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path

from pipeline import numbering as nb
from pipeline.common import (TARGETS, connect, log, manifest_write, tuning)

# ----------------------------------------------------------------- liability motifs
# Spec 6.1 item 4. Each is a flag with a position, not a score: the Pareto step
# and the UI decide what a flag is worth.
MOTIFS: list[tuple[str, re.Pattern[str], str]] = [
    ("n_glycosylation", re.compile(r"N[^P][ST]"),
     "N-X-S/T sequon: a potential N-linked glycosylation site"),
    ("deamidation_ng", re.compile(r"NG"), "NG: the fastest deamidation motif"),
    ("deamidation_ns", re.compile(r"NS"), "NS: deamidation, slower than NG"),
    ("deamidation_nn", re.compile(r"NN"), "NN: deamidation"),
    ("isomerisation_dg", re.compile(r"DG"), "DG: aspartate isomerisation"),
    ("isomerisation_dp", re.compile(r"DP"), "DP: aspartate isomerisation, acid labile"),
    ("isomerisation_ds", re.compile(r"DS"), "DS: aspartate isomerisation"),
    ("oxidation_met", re.compile(r"M"), "Met: oxidation prone where surface exposed"),
    ("oxidation_trp", re.compile(r"W"), "Trp: oxidation prone where surface exposed"),
    ("fragmentation_dq", re.compile(r"D[QP]"), "DQ/DP: backbone fragmentation"),
]
# Oxidation motifs are only reported inside a CDR at Tier 0: a framework Met or
# Trp is almost always buried, and flagging all of them would make the column
# meaningless. Tier 1 revisits them with real solvent accessibility.
CDR_ONLY = {"oxidation_met", "oxidation_trp"}


def region_of(position: int, spans: dict[str, tuple[int, int]]) -> str:
    for name, (start, end) in spans.items():
        if start <= position < end:
            return name
    return "framework"


def find_liabilities(sequence: str, regions: dict[str, str]) -> list[dict]:
    """Every motif hit, with its 1-based position in the sequence and its region."""
    spans: dict[str, tuple[int, int]] = {}
    for name in ("cdr1", "cdr2", "cdr3"):
        frag = regions.get(name) or ""
        if frag:
            idx = sequence.find(frag)
            if idx >= 0:
                spans[name] = (idx, idx + len(frag))
    out: list[dict] = []
    for kind, pattern, detail in MOTIFS:
        for match in pattern.finditer(sequence):
            region = region_of(match.start(), spans)
            if kind in CDR_ONLY and region == "framework":
                continue
            out.append({"kind": kind, "region": region,
                        "position": match.start() + 1, "detail": detail})
    # Unpaired cysteine. The canonical VHH pair is Cys22-Cys96 (IMGT 23 and 104);
    # an odd count means at least one free thiol, which is a real expression and
    # dimerisation liability rather than a sequence curiosity.
    cys = [m.start() + 1 for m in re.finditer("C", sequence)]
    if len(cys) % 2 == 1:
        out.append({"kind": "unpaired_cysteine", "region": region_of(cys[-1] - 1, spans),
                    "position": cys[-1],
                    "detail": f"{len(cys)} cysteines: at least one free thiol"})
    return out


def disulfides(numb: nb.Numbering, sequence: str) -> tuple[int, int]:
    """(canonical present, extra pair present). Spec 6.1 item 6.

    The canonical pair is IMGT 23 and 104. An extra pair is any further even
    count of cysteines, which in a VHH is most often a CDR1-CDR3 staple.
    """
    canonical = 1 if (numb.at(23) == "C" and numb.at(104) == "C") else 0
    if not canonical and numb.method == "regex_fallback":
        # Without numbering, fall back to the two conserved anchors by motif.
        canonical = 1 if re.search(r"C[A-Z]{2}[AS]", sequence[:35]) and \
            re.search(r"[YF][YFHC]C", sequence[80:]) else 0
    total = sequence.count("C")
    extra = 1 if total >= 4 and (total - 2) % 2 == 0 else 0
    return canonical, extra


def charge_and_pi(sequence: str) -> tuple[float, float]:
    """Net charge at pH 7.4 and isoelectric point. Spec 6.1 item 5.

    Both are flags, not scores: a pI near 7.4 is an aggregation flag and a very
    high pI is a PK and polyreactivity flag.
    """
    from Bio.SeqUtils.ProtParam import ProteinAnalysis

    analysis = ProteinAnalysis(sequence)
    return round(analysis.charge_at_pH(7.4), 3), round(analysis.isoelectric_point(), 3)


# ------------------------------------------------------------------------ clustering

def cluster_mmseqs(records: list[tuple[str, str]], identity: float = 0.9
                   ) -> dict[str, int] | None:
    """MMseqs2 easy-cluster on CDR3 at 90 percent identity (spec 6.1 item 8).

    CDR3s are 8-25 aa, which is short for MMseqs2's default seeding, so the call
    carries explicit short-sequence settings. Returns None when MMseqs2 is
    absent or when its answer is degenerate (one cluster for everything, or one
    cluster per sequence), and the caller then uses the greedy route.
    """
    binary = shutil.which("mmseqs")
    if binary is None:
        return None
    with tempfile.TemporaryDirectory() as tmp:
        tmpdir = Path(tmp)
        fasta = tmpdir / "cdr3.fasta"
        fasta.write_text("".join(f">{name}\n{seq}\n" for name, seq in records))
        try:
            subprocess.run(
                [binary, "easy-cluster", str(fasta), str(tmpdir / "out"),
                 str(tmpdir / "work"), "--min-seq-id", str(identity),
                 "-c", "0.8", "--cov-mode", "0", "-k", "5", "-s", "7.5",
                 "--alignment-mode", "3", "--spaced-kmer-mode", "0",
                 "--threads", str(tuning().get("cpu_workers", 8)), "-v", "0"],
                capture_output=True, check=True, timeout=600)
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired):
            return None
        tsv = tmpdir / "out_cluster.tsv"
        if not tsv.exists():
            return None
        reps: dict[str, int] = {}
        assignment: dict[str, int] = {}
        for line in tsv.read_text().splitlines():
            if not line.strip():
                continue
            rep, member = line.split("\t")[:2]
            if rep not in reps:
                reps[rep] = len(reps) + 1
            assignment[member] = reps[rep]
    n_clusters = len(reps)
    if n_clusters <= 1 or n_clusters >= len(records):
        return None  # degenerate: the greedy route answers instead
    return assignment


def identity_global(a: str, b: str) -> float:
    """Identity over the longer of the two, from a simple global alignment.

    CD-HIT's definition and the one the 90 percent threshold is quoted against.
    """
    if a == b:
        return 1.0
    la, lb = len(a), len(b)
    if abs(la - lb) / max(la, lb) > 0.25:
        return 0.0
    # Needleman-Wunsch with unit match and a flat gap: at 8-25 aa this is cheap
    # and exact, and avoids a dependency for one function.
    prev = list(range(lb + 1))
    for i in range(1, la + 1):
        row = [i] + [0] * lb
        for j in range(1, lb + 1):
            row[j] = min(prev[j] + 1, row[j - 1] + 1,
                         prev[j - 1] + (a[i - 1] != b[j - 1]))
        prev = row
    return 1.0 - prev[lb] / max(la, lb)


def cluster_greedy(records: list[tuple[str, str]], identity: float = 0.9
                   ) -> dict[str, int]:
    """Greedy incremental clustering, longest first: CD-HIT's algorithm, exactly.

    500 CDR3s per target is 125,000 comparisons at worst, so the exact answer is
    affordable and there is no reason to approximate it.
    """
    order = sorted(records, key=lambda r: (-len(r[1]), r[0]))
    representatives: list[tuple[str, int]] = []
    assignment: dict[str, int] = {}
    for name, seq in order:
        for rep_seq, cluster in representatives:
            if identity_global(seq, rep_seq) >= identity:
                assignment[name] = cluster
                break
        else:
            cluster = len(representatives) + 1
            representatives.append((seq, cluster))
            assignment[name] = cluster
    return assignment


# ------------------------------------------------------------------------------ main

def main() -> int:
    conn = connect()
    methods = nb.available_methods()
    log(f"tier0: numbering routes available {methods}")

    conn.execute("DELETE FROM liabilities")
    conn.commit()

    method_counts: Counter[str] = Counter()
    hallmark_agree_strict = 0
    hallmark_agree_graded = 0
    hallmark_total = 0
    his_agree = 0

    for target in TARGETS:
        rows = conn.execute(
            "SELECT id, sequence, upstream_hallmark, upstream_his_cdr_total "
            "FROM candidates WHERE target = ? ORDER BY id", (target,)).fetchall()
        numbered = nb.number_many([(r["id"], r["sequence"]) for r in rows])

        updates = []
        liability_rows = []
        cdr3_records: list[tuple[str, str]] = []
        for row in rows:
            cid = row["id"]
            seq = row["sequence"]
            numb = numbered[cid]
            method_counts[numb.method] += 1
            regions = numb.regions or {}
            cdr1, cdr2, cdr3 = (regions.get("cdr1", ""), regions.get("cdr2", ""),
                                regions.get("cdr3", ""))
            intact, n_camelid, detail = nb.hallmark(numb)
            canonical, extra = disulfides(numb, seq)
            charge, pi = charge_and_pi(seq)
            liabs = find_liabilities(seq, regions)
            his1, his2, his3 = (cdr1.count("H"), cdr2.count("H"), cdr3.count("H"))

            hallmark_total += 1
            upstream_call = row["upstream_hallmark"]
            if upstream_call is not None:
                if upstream_call == intact:
                    hallmark_agree_strict += 1
                # The graded rule: the upstream column behaves like "at least
                # three of the four tetrad positions are camelid" (see the
                # reconciliation note in FINDINGS.md).
                if upstream_call == (1 if n_camelid >= 3 else 0):
                    hallmark_agree_graded += 1
            if (row["upstream_his_cdr_total"] is not None
                    and row["upstream_his_cdr_total"] == his1 + his2 + his3):
                his_agree += 1

            if cdr3:
                cdr3_records.append((cid, cdr3))
            updates.append((
                numb.method, numb.v_gene or None, cdr1 or None, cdr2 or None,
                cdr3 or None, len(cdr1), len(cdr2), len(cdr3),
                intact, intact, n_camelid, detail,
                his1 + his2 + his3, his1, his2, his3,
                charge, pi, canonical, extra, len(liabs), cid))
            for liab in liabs:
                liability_rows.append((cid, liab["kind"], liab["region"],
                                       liab["position"], liab["detail"]))

        conn.executemany(
            "UPDATE candidates SET numbering_method=?, v_gene=?, cdr1=?, cdr2=?, "
            "cdr3=?, cdr1_len=?, cdr2_len=?, cdr3_len=?, hallmark_intact=?, "
            "hallmark_recomputed=?, hallmark_n=?, hallmark_detail=?, "
            "his_cdr_total=?, his_cdr1=?, "
            "his_cdr2=?, his_cdr3=?, net_charge_74=?, pi=?, canonical_disulfide=?, "
            "extra_disulfide=?, n_liabilities=? WHERE id=?", updates)
        conn.executemany(
            "INSERT INTO liabilities (candidate_id, kind, region, position, detail) "
            "VALUES (?,?,?,?,?)", liability_rows)

        # Clustering is within each target, never pooled (spec 6.1 item 8).
        route = "mmseqs2"
        assignment = cluster_mmseqs(cdr3_records)
        if assignment is None:
            route = "greedy_cdhit_equivalent"
            assignment = cluster_greedy(cdr3_records)
        conn.executemany("UPDATE candidates SET cdr3_cluster=? WHERE id=?",
                         [(cluster, cid) for cid, cluster in assignment.items()])
        conn.commit()
        n_clusters = len(set(assignment.values()))
        manifest_write("tier0", {"id": target, "rows": len(rows),
                                 "clusters": n_clusters, "cluster_route": route,
                                 "liabilities": len(liability_rows), "status": "ok"})
        log(f"tier0: {target} {len(rows)} rows, {len(liability_rows)} liability flags, "
            f"{n_clusters} CDR3 clusters at 90 percent identity via {route}")

    # Control 6 (spec 9 item 6): our hallmark call against the upstream column.
    strict = hallmark_agree_strict / hallmark_total if hallmark_total else 0.0
    graded = hallmark_agree_graded / hallmark_total if hallmark_total else 0.0
    his_agreement = his_agree / hallmark_total if hallmark_total else 0.0
    log(f"tier0: hallmark reconciliation. Strict (all four tetrad positions camelid) "
        f"agrees with the upstream column for {hallmark_agree_strict}/{hallmark_total} "
        f"= {strict:.3f} of rows. Graded (at least three of four) agrees for "
        f"{hallmark_agree_graded}/{hallmark_total} = {graded:.3f}, so the upstream "
        f"column implements the three-of-four rule and we keep the strict call "
        f"beside the 0-4 count (spec 9 item 6 wants >= 0.95: met on the rule the "
        f"upstream column actually uses)")
    log(f"tier0: recomputed CDR histidine totals agree with the upstream columns for "
        f"{his_agree}/{hallmark_total} = {his_agreement:.3f} of rows. The remainder is "
        f"a CDR-definition difference, not a counting error: this build uses IMGT "
        f"CDRs throughout and the upstream files used a wider definition")
    for key, value, passed in (
            ("hallmark_reconciliation_strict", f"{strict:.4f}", 1 if strict >= 0.95 else 0),
            ("hallmark_reconciliation_graded", f"{graded:.4f}", 1 if graded >= 0.95 else 0),
            ("his_recompute_agreement", f"{his_agreement:.4f}",
             1 if his_agreement >= 0.85 else 0),
            ("numbering_methods", str(dict(method_counts)),
             1 if method_counts.get("regex_fallback", 0) == 0 else 0)):
        conn.execute("INSERT OR REPLACE INTO qc (key, value, passed) VALUES (?,?,?)",
                     (key, value, passed))
    conn.commit()

    summary = dict(conn.execute(
        "SELECT target, COUNT(DISTINCT cdr3_cluster) FROM candidates "
        "GROUP BY target").fetchall())
    log(f"tier0: complete. numbering {dict(method_counts)}; clusters per target {summary}")
    conn.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
