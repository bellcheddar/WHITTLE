"""The 41 table columns, in the five groups DESIGN.md Section 3.4 names.

One registry, baked into app/static/data/meta.json, so the Tabulator table, the
column chooser and the CSV export cannot drift apart.

DESIGN.md gives the group order (Identity, Binding, Specificity and liabilities,
Thermostability, Expression and evidence) and a total of 41. Its per-group counts
(6, 7, 11, 5, 8) sum to 37, so they cannot all hold at once: the order and the
total of 41 are honoured here, and the chooser renders the real per-group counts
from this registry rather than any hard-coded number (see DECISIONS.md).

Each column carries:
  key        the field name in the row bundle
  label      what the header says
  group      one of the five group ids
  unit       shown in ink-3 on the right of the chooser row, and under the header
  kind       num | pct | text | chip | spark | badge  (drives the formatter)
  tier       0, 1, 2 or None: which tier fills it, so the UI can say "pending"
  default    True if visible with no column state in the URL
  help       one sentence, shown as the header tooltip and in the chooser
  decimals   for numeric formatting
  higher_is_better  True, False or None: drives the percentile colour scale
"""
from __future__ import annotations

GROUPS = [
    ("identity", "Identity"),
    ("binding", "Binding"),
    ("specificity", "Specificity and liabilities"),
    ("thermo", "Thermostability"),
    ("expression", "Expression and evidence"),
]

PCT_NOTE = "Percentile within this target. Not comparable between targets."

COLUMNS: list[dict] = [
    # ------------------------------------------------------------------ identity (9)
    dict(key="id", label="Candidate", group="identity", unit="", kind="text",
         tier=0, default=True, frozen=True,
         help="Candidate id. The target is shown beneath it."),
    dict(key="target", label="Target", group="identity", unit="", kind="text",
         tier=0, default=False,
         help="The target this candidate was designed against. One target at a time."),
    dict(key="rank", label="Rank", group="identity", unit="of 500", kind="num",
         tier=0, default=True, decimals=0, higher_is_better=False,
         help="Rank within this target by desirability. Pareto rank breaks ties."),
    dict(key="axis_profile", label="Axis profile", group="identity", unit="4 axes",
         kind="spark", tier=1, default=True,
         help="Affinity, specificity, thermostability and expression percentiles as "
              "four bars, so the shape of a candidate reads without opening it."),
    dict(key="length", label="Length", group="identity", unit="aa", kind="num",
         tier=0, default=False, decimals=0, higher_is_better=None,
         help="Sequence length."),
    dict(key="v_gene", label="V gene", group="identity", unit="IMGT", kind="text",
         tier=0, default=False,
         help="Nearest germline V gene assigned by ANARCI."),
    dict(key="cdr3_len", label="CDR3 length", group="identity", unit="aa", kind="num",
         tier=0, default=False, decimals=0, higher_is_better=None,
         help="IMGT CDR3 length. Long CDR3s buy surface and cost polyreactivity."),
    dict(key="cdr3_cluster", label="CDR3 cluster", group="identity", unit="id",
         kind="num", tier=0, default=False, decimals=0, higher_is_better=None,
         help="Cluster id at 90 percent CDR3 identity within this target. A diverse "
              "shortlist spreads across clusters."),
    dict(key="multi_target", label="Multi-target", group="identity", unit="flag",
         kind="chip", tier=0, default=False,
         help="This exact sequence also appears under another target. Treated as an "
              "internal promiscuity control, not as a better candidate."),

    # ------------------------------------------------------------------- binding (7)
    dict(key="axis_affinity", label="Affinity", group="binding", unit="percentile",
         kind="pct", tier=2, default=False, decimals=0, higher_is_better=True,
         help="Co-folding confidence turned into a percentile within this target. "
              "A weak prior for specificity, not an affinity measurement. " + PCT_NOTE),
    dict(key="iptm", label="ipTM", group="binding", unit="0-1", kind="num",
         tier=2, default=False, decimals=3, higher_is_better=True,
         help="Interface predicted TM-score, averaged over seeds. Estimated, not measured."),
    dict(key="ptm", label="pTM", group="binding", unit="0-1", kind="num",
         tier=2, default=False, decimals=3, higher_is_better=True,
         help="Predicted TM-score for the whole complex, averaged over seeds."),
    dict(key="prodigy_dg", label="PRODIGY dG", group="binding", unit="kcal/mol",
         kind="num", tier=2, default=False, decimals=2, higher_is_better=False,
         help="Predicted interface binding free energy from the co-folded model."),
    dict(key="bsa_total", label="Buried area", group="binding", unit="A2", kind="num",
         tier=2, default=False, decimals=0, higher_is_better=True,
         help="Total buried surface area across the interface."),
    dict(key="n_subunits_contacted", label="Subunits hit", group="binding",
         unit="of 1 or 3", kind="num", tier=2, default=False, decimals=0,
         higher_is_better=None,
         help="For the TNFa and TL1A trimers: how many subunits the interface touches. "
              "An epitope spanning the subunit cleft contacts more than one."),
    dict(key="epitope_bin", label="Epitope bin", group="binding", unit="id", kind="num",
         tier=2, default=False, decimals=0, higher_is_better=None,
         help="Candidates clustered by their interface residue set, so a shared face "
              "is visible."),

    # ------------------------------------------------- specificity and liabilities (11)
    dict(key="axis_specificity", label="Specificity", group="specificity",
         unit="percentile", kind="pct", tier=1, default=False, decimals=0,
         higher_is_better=True,
         help="Non-specific binding risk inverted into a percentile, driven by TNP. "
              + PCT_NOTE),
    dict(key="liabilities", label="Liabilities", group="specificity", unit="flags",
         kind="chip", tier=0, default=True,
         help="Sequence liability motifs found: glycosylation sequons, deamidation, "
              "isomerisation, unpaired cysteine, exposed Met and Trp."),
    dict(key="psh", label="PSH", group="specificity", unit="TNP", kind="num",
         tier=1, default=False, decimals=2, higher_is_better=False,
         help="Patch surface hydrophobicity from TNP. High values predict "
              "non-specific binding."),
    dict(key="ppc", label="PPC", group="specificity", unit="TNP", kind="num",
         tier=1, default=False, decimals=2, higher_is_better=None,
         help="Patch positive charge from TNP."),
    dict(key="pnc", label="PNC", group="specificity", unit="TNP", kind="num",
         tier=1, default=False, decimals=2, higher_is_better=None,
         help="Patch negative charge from TNP."),
    dict(key="tnp_flags", label="TNP flags", group="specificity", unit="count",
         kind="num", tier=1, default=False, decimals=0, higher_is_better=False,
         help="How many TNP metrics fall outside the clinical-stage antibody range."),
    dict(key="n_glyc", label="Glyc sequons", group="specificity", unit="count",
         kind="num", tier=0, default=False, decimals=0, higher_is_better=False,
         help="N-X-S/T sequons. A sequon in a CDR is a developability problem."),
    dict(key="n_deamidation", label="Deamidation", group="specificity", unit="count",
         kind="num", tier=0, default=False, decimals=0, higher_is_better=False,
         help="NG, NS and NN motifs. NG is the fastest."),
    dict(key="n_isomerisation", label="Isomerisation", group="specificity",
         unit="count", kind="num", tier=0, default=False, decimals=0,
         higher_is_better=False,
         help="DG, DP and DS motifs: aspartate isomerisation and acid-labile bonds."),
    dict(key="net_charge_74", label="Net charge", group="specificity", unit="at pH 7.4",
         kind="num", tier=0, default=False, decimals=2, higher_is_better=None,
         help="Net charge at pH 7.4. A flag, not a score: very high charge is a PK "
              "and polyreactivity flag."),
    dict(key="pi", label="pI", group="specificity", unit="pH", kind="num",
         tier=0, default=False, decimals=2, higher_is_better=None,
         help="Isoelectric point. A pI near 7.4 is an aggregation flag."),

    # -------------------------------------------------------------- thermostability (6)
    dict(key="axis_thermo", label="Thermostability", group="thermo", unit="percentile",
         kind="pct", tier=1, default=False, decimals=0, higher_is_better=True,
         help="Predicted apparent Tm turned into a percentile within this target. "
              + PCT_NOTE),
    dict(key="tm_pred", label="Tm", group="thermo", unit="deg C", kind="num",
         tier=1, default=True, decimals=1, higher_is_better=True,
         help="Predicted apparent Tm from the WHITTLE model. An estimate with an "
              "error bar, not a measurement."),
    dict(key="tm_err", label="Tm error", group="thermo", unit="+/- deg C", kind="num",
         tier=1, default=False, decimals=1, higher_is_better=False,
         help="Error estimate from the distance to the nearest training neighbour. "
              "Wide means the model is extrapolating."),
    dict(key="temstapro", label="TemStaPro", group="thermo", unit="0-1", kind="num",
         tier=1, default=False, decimals=3, higher_is_better=True,
         help="Independent thermostability cross-check. Reported, not blended in."),
    dict(key="upstream_tm_pred", label="Tm (upstream)", group="thermo", unit="deg C",
         kind="num", tier=0, default=False, decimals=1, higher_is_better=None,
         help="The Tm column from the source files. Uncalibrated: its range reaches "
              "101 deg C, which is not a value a DSF run gives a VHH in solution. "
              "Shown for comparison only."),
    dict(key="hallmark_n", label="Hallmark", group="thermo", unit="of 4", kind="num",
         tier=0, default=False, decimals=0, higher_is_better=True,
         help="How many of the four camelid hallmark positions (IMGT 42, 49, 50, 52) "
              "carry a camelid residue. A VH-type residue leaves the former VL "
              "interface hydrophobic and exposed."),

    # ------------------------------------------------- expression and evidence (8)
    dict(key="axis_expression", label="Expression", group="expression",
         unit="percentile", kind="pct", tier=1, default=False, decimals=0,
         higher_is_better=True,
         help="Expression and solubility proxy as a percentile within this target. "
              + PCT_NOTE),
    dict(key="desirability", label="Desirability", group="expression", unit="0-1",
         kind="num", tier=1, default=False, decimals=3, higher_is_better=True,
         help="Derringer-Suich desirability over the four axes with the weights set "
              "in the toolbar. One ordering for anyone who wants one."),
    dict(key="swi", label="SWI", group="expression", unit="0-1", kind="num",
         tier=1, default=False, decimals=3, higher_is_better=True,
         help="Solubility-Weighted Index from the published formula."),
    dict(key="protein_sol", label="Protein-Sol", group="expression", unit="0-1",
         kind="num", tier=1, default=False, decimals=3, higher_is_better=True,
         help="Protein-Sol solubility from the published formula. 0.45 is the "
              "average for the E. coli reference set."),
    dict(key="nativeness_ppl", label="Nativeness", group="expression",
         unit="pseudo-ppl", kind="num", tier=1, default=False, decimals=2,
         higher_is_better=False,
         help="Pseudo-perplexity under an antibody language model. Odd predicted "
              "sequences score badly for a reason."),
    dict(key="oasis_humanness", label="Humanness", group="expression", unit="OASis",
         kind="num", tier=1, default=False, decimals=3, higher_is_better=True,
         help="BioPhi OASis identity. Relevant to immunogenicity, not to expression "
              "in E. coli."),
    dict(key="nearest_pdb", label="Nearest PDB", group="expression", unit="entry",
         kind="text", tier=1, default=True,
         help="Closest experimental nanobody structure in SAbDab-nano, with its "
              "TM-score and per-region identity in the drawer."),
    dict(key="experimental_badge", label="Experimental", group="expression",
         unit="badge", kind="badge", tier=1, default=False,
         help="A near-identical sequence has a PDB entry. A positive signal only: "
              "its absence says nothing about expressibility."),
]

BY_KEY = {c["key"]: c for c in COLUMNS}
DEFAULT_KEYS = [c["key"] for c in COLUMNS if c.get("default")]


def group_counts() -> dict[str, int]:
    out: dict[str, int] = {gid: 0 for gid, _ in GROUPS}
    for col in COLUMNS:
        out[col["group"]] += 1
    return out


def registry() -> dict:
    """What goes into meta.json."""
    return {
        "groups": [{"id": gid, "label": label, "count": group_counts()[gid]}
                   for gid, label in GROUPS],
        "columns": COLUMNS,
        "total": len(COLUMNS),
        "defaults": DEFAULT_KEYS,
    }


def _self_check() -> None:
    assert len(COLUMNS) == 41, f"the registry must hold 41 columns, found {len(COLUMNS)}"
    assert len(BY_KEY) == len(COLUMNS), "duplicate column key"
    assert len(DEFAULT_KEYS) == 6, f"six default columns, found {len(DEFAULT_KEYS)}"
    assert {c["group"] for c in COLUMNS} == {gid for gid, _ in GROUPS}


_self_check()

if __name__ == "__main__":
    import json

    print(json.dumps({"total": len(COLUMNS), "groups": group_counts(),
                      "defaults": DEFAULT_KEYS}, indent=2))
