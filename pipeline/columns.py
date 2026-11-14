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
    dict(key="pose_consensus", label="Pose agreement", group="binding",
         unit="0-1", kind="num", tier=2, default=False, decimals=2,
         higher_is_better=True,
         help="How far the five diffusion samples agree about where this binds, as "
              "the mean overlap of their epitopes. A candidate can bury a lot of "
              "surface somewhere different every time, and a single good-looking "
              "sample out of five is a lucky draw rather than evidence."),
    dict(key="sc", label="Shape complementarity", group="binding", unit="0-1",
         kind="num", tier=2, default=False, decimals=3, higher_is_better=True,
         help="Lawrence-Colman shape complementarity: how well the two surfaces fit, "
              "independent of how large the interface is. Experimental "
              "nanobody-antigen interfaces measure about 0.57 to 0.65 under this "
              "implementation. A predicted pose scores lower, so treat it as a "
              "ranking signal within this target rather than as a figure to compare "
              "with the literature."),
    dict(key="cdr_interface_fraction", label="CDR paratope", group="binding",
         unit="fraction", kind="num", tier=2, default=False, decimals=2,
         higher_is_better=True,
         help="Share of the nanobody's contacting residues that sit in a CDR rather "
              "than the framework. Framework binding is either a modelling artefact "
              "or a contact unlikely to be specific, since the framework is the "
              "conserved part."),
    dict(key="n_hbond", label="H-bonds", group="binding", unit="count", kind="num",
         tier=2, default=False, decimals=0, higher_is_better=True,
         help="Hydrogen bonds across the interface, detected with PLIP on the "
              "highest-ipTM sample. A count rises with interface size, so read it "
              "beside the density."),
    dict(key="hbond_per_100a2", label="H-bond density", group="binding",
         unit="per 100 A2", kind="num", tier=2, default=False, decimals=2,
         higher_is_better=True,
         help="Hydrogen bonds per 100 square angstroms of buried surface. This is "
              "the part of the interaction inventory that does not simply restate "
              "how large the interface is."),
    dict(key="polar_fraction", label="Polar fraction", group="binding",
         unit="0-1", kind="num", tier=2, default=False, decimals=2,
         higher_is_better=None,
         help="Share of detected interactions that are polar (hydrogen bonds and "
              "salt bridges) rather than hydrophobic. A flag for the character of "
              "the interface, not a score."),
    dict(key="n_his_contacts", label="His contacts", group="binding",
         unit="count", kind="num", tier=2, default=False, decimals=0,
         higher_is_better=None,
         help="Interface contacts made by a nanobody histidine. The EGFR and TNFa "
              "libraries are histidine-enriched for an intended pH switch, so a "
              "candidate from those libraries whose histidines never reach the "
              "antigen has not done what it was designed to do."),

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


# ---------------------------------------------------------------------------
# The glossary that the About tab renders.
#
# `help` above is one sentence for a header tooltip. That is the wrong length
# for a methods page, where a reader who is not a protein engineer needs three
# different things: what the number means in plain words, how it was arrived
# at, and why it should change anyone's mind about a candidate. Those are the
# three fields here.
#
# This lives beside the registry rather than in a YAML file of its own so that
# a column cannot be added without the glossary entry being noticed: the
# self-check at the bottom of this module fails if any column lacks one.
#
# Each entry:
#   plain  what the number is, in lay terms
#   how    how it is computed, and from what
#   why    what it changes about whether a candidate is shortlisted
#   src    the reference key in config/references.yaml, or "" if it is ours

G = dict  # a short alias, so the table below stays readable

GLOSSARY: dict[str, dict] = {

    # ---------------------------------------------------------- identity
    "id": G(plain="The candidate's name. Every candidate belongs to exactly one "
                  "target and keeps the same name everywhere on this site.",
            how="Assigned when the library was read in, as the target name and a "
                "four-digit number.",
            why="It is the handle for ordering a gene. Nothing is scored from it.",
            src=""),
    "target": G(plain="Which of the five drug targets this nanobody was designed "
                      "to bind.",
                how="Taken from the source library. Each of the five libraries "
                    "holds 500 candidates.",
                why="Everything is ranked inside one target. The five libraries "
                    "differ in size, construct and design strategy, so a score "
                    "from one of them says nothing about another.",
                src=""),
    "rank": G(plain="Where this candidate sits in the queue for its own target, "
                    "1 being the most promising.",
              how="Candidates are sorted by desirability across the four axes. "
                  "Where two tie, the one on a better Pareto front wins.",
              why="It is the summary of every other column on this page, and the "
                  "order the shortlist is taken from. It is a suggestion about "
                  "ranking, not a prediction that one binder is better than "
                  "another.",
              src=""),
    "axis_profile": G(plain="Four small bars showing how this candidate scores on "
                            "binding, non-specific binding, thermostability and "
                            "expression at once.",
                      how="The four axis percentiles drawn as bar lengths, in a "
                          "fixed order so the shape is comparable between rows.",
                      why="A candidate that is excellent on one axis and poor on "
                          "another is a different proposition from one that is "
                          "adequate on all four, and the shape shows that "
                          "without opening anything.",
                      src=""),
    "length": G(plain="How many amino acids long the nanobody is.",
                how="Counted from the sequence.",
                why="VHH domains run about 110 to 135 residues. A sequence far "
                    "outside that is either mis-annotated or carries a tag, and "
                    "is worth looking at before it is ordered.",
                src=""),
    "v_gene": G(plain="The germline gene the camelid immune system started from "
                      "when it made this shape.",
                how="ANARCI assigns the closest germline V gene while numbering "
                    "the sequence under the IMGT scheme.",
                why="A sequence far from any germline is more likely to be a "
                    "design artefact. It is also the first thing anyone asks "
                    "when comparing a hit to a published nanobody.",
                src="anarci"),
    "cdr3_len": G(plain="The length of the long loop that does most of the "
                        "binding.",
                  how="The IMGT CDR3 boundaries from ANARCI, counted in residues.",
                  why="Camelid CDR3s are long by antibody standards and that "
                      "length is what lets a nanobody reach into a cleft. Length "
                      "buys surface area and costs specificity: very long loops "
                      "are flexible, more often sticky, and harder to express.",
                  src="anarci"),
    "cdr3_cluster": G(plain="A group number shared by candidates whose binding "
                            "loop is nearly the same sequence.",
                      how="CDR3 sequences within one target are clustered at 90 "
                          "per cent identity. Every member of a cluster gets the "
                          "same number.",
                      why="Twenty-four candidates drawn from one cluster are "
                          "twenty-four attempts at the same binder. Spreading "
                          "the shortlist across clusters is the cheapest "
                          "insurance available, because a whole cluster can fail "
                          "for one reason.",
                      src=""),
    "multi_target": G(plain="This exact sequence also appears in another target's "
                            "library.",
                      how="Exact sequence match across all 2,500 candidates.",
                      why="A sequence designed against two different antigens "
                          "cannot be specific to both. It is treated as an "
                          "internal promiscuity control, never as a candidate "
                          "that is twice as good.",
                      src=""),

    # ----------------------------------------------------------- binding
    "axis_affinity": G(plain="How confident the structure predictor is that this "
                             "nanobody and this antigen form a sensible complex.",
                       how="A weighted blend of the interface metrics from the "
                           "co-folded model, converted to a percentile inside "
                           "this target. The weights are set from how well each "
                           "metric separates real candidates from CDR-shuffled "
                           "controls.",
                       why="This is the weakest axis on the site and it is "
                           "labelled as such everywhere. Published benchmarks "
                           "show co-folding confidence does not reliably tell a "
                           "cognate nanobody-antigen pair from a non-cognate "
                           "one. Treat a low score as a reason to look harder "
                           "and a high score as no promise at all.",
                       src="boltz2"),
    "iptm": G(plain="The predictor's own confidence in the interface between the "
                    "two chains, from 0 to 1.",
              how="Boltz reports an interface predicted TM-score for every "
                  "diffusion sample. Five samples are run per complex and the "
                  "value shown is their mean.",
              why="It is the headline co-folding number and the one most often "
                  "quoted. Above about 0.8 the model is confident; that "
                  "confidence is about geometry being self-consistent, not "
                  "about the complex existing in solution.",
              src="boltz2"),
    "ptm": G(plain="The predictor's confidence in the whole complex, not just the "
                   "interface.",
             how="Predicted TM-score over all chains, averaged over the five "
                 "diffusion samples.",
             why="A high interface score sitting under a low whole-complex score "
                 "usually means the antigen itself modelled badly, which makes "
                 "the interface number unsafe to read.",
             src="boltz2"),
    "prodigy_dg": G(plain="An estimate of how tightly the two proteins are held "
                          "together, in kcal/mol. More negative is tighter.",
                    how="PRODIGY counts contacts by residue type across the "
                        "interface of the co-folded model and puts them through "
                        "a published linear model fitted to measured affinities.",
                    why="It is the only column expressed in energy units, which "
                        "makes it the easiest to over-read. PRODIGY was fitted "
                        "to crystal structures, and this is a predicted pose, so "
                        "it ranks within this target and does not convert to a "
                        "Kd.",
                    src="prodigy"),
    "bsa_total": G(plain="How much surface area the two proteins bury against "
                         "each other, in square angstroms.",
                   how="Solvent-accessible area of each chain alone, minus the "
                       "area in the complex, computed with FreeSASA.",
                   why="Real nanobody interfaces bury roughly 600 to 900 square "
                       "angstroms. Far below that is barely an interface; far "
                       "above it often means the model has pressed the chains "
                       "together rather than docked them.",
                   src="freesasa"),
    "n_subunits_contacted": G(plain="For the three-part targets, how many of the "
                                    "three subunits this nanobody touches.",
                              how="Counted from the contacting residues of the "
                                  "co-folded model. TNFa and TL1A are co-folded "
                                  "as homotrimers; the others are single chains.",
                              why="Many real TNF-family epitopes sit in the "
                                  "groove between two subunits. A candidate "
                                  "touching two subunits is engaging the "
                                  "biologically interesting face; one touching "
                                  "a single subunit may be binding a surface "
                                  "that is hidden in the assembled trimer.",
                              src=""),
    "epitope_bin": G(plain="A group number shared by candidates that land on the "
                           "same patch of the antigen.",
                     how="Candidates are clustered by how much their sets of "
                         "contacted antigen residues overlap.",
                     why="A shortlist that all bins into one epitope is a single "
                         "bet. Two candidates from different bins fail "
                         "independently, which is the point of carrying "
                         "twenty-four rather than one.",
                     src=""),
    "pose_consensus": G(plain="How often the five independent predictions put the "
                              "nanobody in the same place on the antigen. 1 means "
                              "all five agree.",
                        how="Mean pairwise overlap of the contacted antigen "
                            "residues across the five diffusion samples, as a "
                            "Jaccard index.",
                        why="A candidate can bury plenty of surface in a "
                            "different place every time, and one good-looking "
                            "sample out of five is a lucky draw rather than "
                            "evidence. This column is the difference between a "
                            "prediction and a coincidence, and it is the reason "
                            "five seeds are run instead of one.",
                        src=""),
    "sc": G(plain="How well the two surfaces fit together, like two halves of a "
                  "mould. 1 would be perfect.",
            how="The Lawrence and Colman shape complementarity statistic, "
                "computed over the interface dot surface of the co-folded model "
                "with a 1.7 angstrom probe.",
            why="It measures fit independently of size, so a small neat "
                "interface is not punished for being small. Experimental "
                "nanobody-antigen interfaces measure about 0.57 to 0.65 under "
                "this implementation and predicted poses sit lower, so read it "
                "as a ranking signal inside this target and not as a number to "
                "compare with a paper.",
            src="lawrence_colman"),
    "cdr_interface_fraction": G(plain="What share of the nanobody's contact with "
                                      "the antigen is made by the three binding "
                                      "loops rather than the scaffold.",
                                how="Contacting nanobody residues are assigned to "
                                    "CDR or framework using the IMGT numbering, "
                                    "and the CDR share is reported.",
                                why="The framework is the conserved part shared "
                                    "by every nanobody, so a contact there "
                                    "cannot be specific to this antigen. A low "
                                    "share is either a modelling artefact or a "
                                    "binder that would stick to a great many "
                                    "other things.",
                                src=""),
    "n_hbond": G(plain="How many hydrogen bonds the predicted model makes across "
                       "the interface.",
                 how="PLIP profiles the highest-confidence sample with "
                     "protein-protein detection enabled. Duplicate detections "
                     "are collapsed, keeping the best geometry.",
                 why="Hydrogen bonds are the directional part of binding and the "
                     "part that confers specificity. The raw count rises with "
                     "interface size, so it should always be read beside the "
                     "density below it.",
                 src="plip"),
    "hbond_per_100a2": G(plain="How densely packed those hydrogen bonds are, per "
                               "100 square angstroms of buried surface.",
                         how="Hydrogen bond count divided by the buried surface "
                             "area, scaled to 100 square angstroms.",
                         why="This is the part of the contact inventory that does "
                             "not simply restate how large the interface is. A "
                             "large interface with few bonds per unit area is a "
                             "greasy patch pressed against another greasy patch, "
                             "which is the signature of a sticky candidate "
                             "rather than a selective one.",
                         src="plip"),
    "polar_fraction": G(plain="What share of the detected interactions are polar "
                              "rather than greasy.",
                        how="Hydrogen bonds and salt bridges divided by all "
                            "detected interactions, hydrophobic contacts "
                            "included.",
                        why="It describes the character of the interface rather "
                            "than its quality, which is why it is a flag and not "
                            "a score. A strongly hydrophobic interface is the "
                            "kind that also binds the column, the membrane and "
                            "the next protein along.",
                        src="plip"),
    "n_his_contacts": G(plain="How many contacts with the antigen are made by a "
                              "histidine on the nanobody.",
                        how="Counted from the interface residues of the co-folded "
                            "model.",
                        why="The EGFR and TNFa libraries were deliberately "
                            "enriched in histidine so that binding would weaken "
                            "in the acidic endosome and release the antigen. A "
                            "candidate from those libraries whose histidines "
                            "never reach the antigen has not done the thing it "
                            "was designed to do, whatever else it scores.",
                        src=""),

    # ------------------------------------------------- specificity group
    "axis_specificity": G(plain="How likely this nanobody is to stick to things "
                                "it was not designed for. A high percentile means "
                                "a low risk.",
                          how="Driven by the Therapeutic Antibody Profiler patch "
                              "metrics on the predicted structure, with sequence "
                              "liabilities folded in, inverted so that higher is "
                              "better and converted to a percentile within the "
                              "target.",
                          why="Non-specific binding is the single most common "
                              "reason a good binder fails later: it drives fast "
                              "clearance, off-target signal in an assay and "
                              "unreadable imaging. This is the axis WHITTLE is "
                              "most confident about, because it rests on "
                              "structure-derived surface properties rather than "
                              "on predicted binding.",
                          src="tnp"),
    "liabilities": G(plain="Short stretches of sequence known to degrade, "
                           "rearrange or attract sugar.",
                     how="Pattern matching over the sequence for glycosylation "
                         "sequons, deamidation and isomerisation motifs, unpaired "
                         "cysteines and exposed methionine and tryptophan.",
                     why="Each one is a known manufacturing or shelf-life "
                         "problem, and one sitting inside a CDR is far worse than "
                         "the same motif in the framework because it degrades the "
                         "part that binds. They are cheap to avoid at this stage "
                         "and expensive to discover after a batch is made.",
                     src=""),
    "psh": G(plain="How large and how greasy the biggest hydrophobic patches on "
                   "the surface are.",
             how="Patch surface hydrophobicity from the Therapeutic Antibody "
                 "Profiler, computed on the predicted structure over the CDR "
                 "vicinity.",
             why="It is the strongest single predictor of non-specific binding in "
                 "the clinical-stage antibody literature, and it is the column "
                 "that most often decides whether a candidate survives Tier 1. "
                 "Lower is better.",
             src="tnp"),
    "ppc": G(plain="How much positive charge is concentrated in one patch on the "
                   "surface.",
             how="Patch positive charge from the Therapeutic Antibody Profiler, "
                 "over the CDR vicinity of the predicted structure.",
             why="A concentrated positive patch binds the negatively charged cell "
                 "surface and nucleic acids, which shows up as fast clearance and "
                 "as background in a binding assay. It is a flag rather than a "
                 "score, because some targets need it.",
             src="tnp"),
    "pnc": G(plain="How much negative charge is concentrated in one patch on the "
                   "surface.",
             how="Patch negative charge from the Therapeutic Antibody Profiler, "
                 "over the CDR vicinity of the predicted structure.",
             why="The counterpart to the positive patch. Extremes in either "
                 "direction are what matter, which is why neither has a "
                 "preferred direction here.",
             src="tnp"),
    "tnp_flags": G(plain="How many of the five surface measures fall outside the "
                         "range seen in antibodies that reached the clinic.",
                   how="Each Therapeutic Antibody Profiler metric is compared "
                       "with the published clinical-stage range and the "
                       "out-of-range ones are counted.",
                   why="It converts five continuous numbers into one count that "
                       "can be scanned. Nothing is rejected on this column alone, "
                       "because the reference set is conventional antibodies and "
                       "a VHH is not one.",
                   src="tnp"),
    "n_glyc": G(plain="How many places a sugar chain could be attached.",
                how="Counts of the N-X-S/T pattern, where X is any residue but "
                    "proline.",
                why="A sugar attached inside a CDR can block binding outright, "
                    "and glycosylation differs between expression systems, so the "
                    "molecule made in yeast is not the one made in mammalian "
                    "cells. A sequon in the framework is usually tolerable.",
                src=""),
    "n_deamidation": G(plain="How many places the sequence can slowly chemically "
                             "rearrange on the shelf.",
                       how="Counts of the NG, NS and NN motifs. NG is by far the "
                           "fastest.",
                       why="Deamidation changes the charge of a residue over "
                           "weeks, so a batch drifts after it is made. In a CDR "
                           "it changes binding; everywhere it complicates the "
                           "analytics.",
                       src=""),
    "n_isomerisation": G(plain="How many places an aspartate can isomerise, or a "
                               "bond can break under mild acid.",
                         how="Counts of the DG, DP and DS motifs.",
                         why="The same shelf-life argument as deamidation, with "
                             "the added problem that an isomerised aspartate is "
                             "hard to separate from the intact one.",
                         src=""),
    "net_charge_74": G(plain="The overall electrical charge of the molecule at "
                             "blood pH.",
                       how="Summed from the ionisable residues using standard pKa "
                           "values at pH 7.4.",
                       why="Strongly positive molecules are cleared quickly and "
                           "bind off-target surfaces; strongly negative ones can "
                           "be hard to purify. It is a flag rather than a score, "
                           "because the right charge depends on the target and "
                           "the format.",
                       src=""),
    "pi": G(plain="The pH at which the molecule carries no net charge.",
            how="Computed from the sequence with standard pKa values.",
            why="A protein is least soluble at its isoelectric point, so a pI "
                "near 7.4 means the molecule is at its least soluble in the "
                "buffer it will be used in. It also sets which purification "
                "chemistry will work.",
            src=""),

    # ------------------------------------------------- thermostability
    "axis_thermo": G(plain="How well the molecule is expected to hold its shape "
                           "when heated.",
                     how="The predicted apparent melting temperature, converted "
                         "to a percentile within this target and widened by its "
                         "own error estimate.",
                     why="Thermostability travels with resistance to aggregation, "
                         "to freeze-thaw and to long storage. It is the property "
                         "most often wished for after the fact, because a weakly "
                         "folded binder works in an assay on the day and fails "
                         "every time it is shipped.",
                     src="nanomelt_method"),
    "tm_pred": G(plain="The temperature at which half the molecules are expected "
                       "to have unfolded, in degrees Celsius.",
                 how="A gradient-boosted model trained here on 650 nanobodies "
                     "with measured melting temperatures, using sequence and "
                     "structure-derived features, validated by grouped nested "
                     "cross-validation so that near-identical sequences cannot "
                     "sit on both sides of a split.",
                 why="VHH domains typically melt between 55 and 80 degrees. This "
                     "is an estimate with a real error bar beside it, not a "
                     "measurement, and the error bar is the more important of "
                     "the two columns.",
                 src="nanomelt_method"),
    "tm_err": G(plain="How much to distrust the melting temperature beside it.",
                how="Derived from how far this sequence sits from the nearest "
                    "sequences in the training set. Far away means the model is "
                    "extrapolating.",
                why="A predicted 78 degrees with an error of 3 is a useful "
                    "statement. The same 78 with an error of 15 says only that "
                    "the model has not seen anything like this sequence, and "
                    "ranking on the first number alone would be reading noise.",
                src=""),
    "temstapro": G(plain="A second, independent opinion on thermostability from a "
                         "published tool.",
                   how="TemStaPro, a protein language model classifier, run on "
                       "the sequence.",
                   why="It is reported beside our own model rather than blended "
                       "into it, so that agreement is visible as agreement. Two "
                       "models trained on overlapping data agreeing is weaker "
                       "evidence than it looks, but two models disagreeing is a "
                       "clear warning.",
                   src="temstapro"),
    "upstream_tm_pred": G(plain="The melting temperature that came with the source "
                                "files, shown for comparison only.",
                          how="Copied unchanged from the supplied library "
                              "columns.",
                          why="Its range reaches 101 degrees, which is not a "
                              "value a melting experiment gives a VHH in "
                              "solution, so the column is uncalibrated. It is "
                              "kept visible rather than deleted so that anyone "
                              "comparing WHITTLE with the source can see exactly "
                              "where the two differ.",
                          src=""),
    "hallmark_n": G(plain="How many of the four signature camelid positions are "
                          "intact, out of four.",
                    how="IMGT positions 42, 49, 50 and 52 are checked against the "
                        "camelid consensus.",
                    why="Those four positions are what let a VHH be a single "
                        "domain: they replace the greasy surface that would "
                        "otherwise be covered by a light chain. A candidate "
                        "missing one has an exposed hydrophobic face, which "
                        "shows up later as aggregation and as non-specific "
                        "binding.",
                    src="anarci"),

    # --------------------------------------------- expression and evidence
    "axis_expression": G(plain="How likely the molecule is to be made in useful "
                               "amounts and stay dissolved.",
                         how="A blend of two published solubility predictors with "
                             "the language-model nativeness score, converted to a "
                             "percentile within this target.",
                         why="A candidate that cannot be made cannot be tested, "
                             "and this is the cheapest failure to avoid. It is "
                             "also the axis the scrambled control cannot test, "
                             "because shuffling a loop leaves the amino acid "
                             "composition unchanged and these measures read "
                             "composition.",
                         src="protein_sol"),
    "desirability": G(plain="One overall score from 0 to 1 combining all four "
                            "axes, for anyone who wants a single ordering.",
                      how="The Derringer-Suich desirability function: each axis "
                          "is mapped to 0 to 1, then the geometric mean is taken "
                          "with the weights set in the toolbar.",
                      why="The geometric mean is used deliberately, because it "
                          "goes to zero if any one axis does. A candidate that is "
                          "outstanding on three axes and hopeless on the fourth "
                          "scores badly here, which is the correct behaviour for "
                          "a molecule that has to do four things at once.",
                      src="derringer_suich"),
    "swi": G(plain="A solubility estimate read straight off the sequence.",
             how="The Solubility-Weighted Index, a published per-residue weighted "
                 "average.",
             why="It is fast, it needs no structure, and it was fitted to "
                 "bacterial expression data, which is how these candidates would "
                 "first be made. It reads composition only, so it cannot "
                 "distinguish two sequences built from the same residues in a "
                 "different order.",
             src="swi"),
    "protein_sol": G(plain="A second solubility estimate from a different "
                           "published formula.",
                     how="The Protein-Sol sequence model. The average for its "
                         "E. coli reference set is 0.45, so that value is the "
                         "break-even point rather than the middle of the scale.",
                     why="Carried beside SWI for the same reason TemStaPro is "
                         "carried beside the melting model: two independent "
                         "formulas agreeing is worth more than either alone. It "
                         "is also composition-only.",
                     src="protein_sol"),
    "nativeness_ppl": G(plain="How ordinary this sequence looks to a model trained "
                              "on real antibodies. Lower is more ordinary.",
                        how="Pseudo-perplexity under AbLang2, an antibody "
                            "language model.",
                        why="Designed sequences can be perfectly valid and still "
                            "look nothing like anything an immune system has "
                            "made. Odd sequences express badly and fold badly "
                            "often enough that this is worth a column, and it is "
                            "one of the few expression-side measures that "
                            "responds to residue order rather than composition.",
                        src="ablang2"),
    "oasis_humanness": G(plain="How much of this sequence would be recognised as "
                               "ordinary by a human immune system.",
                         how="BioPhi OASis identity against human antibody "
                             "repertoires.",
                         why="It speaks to immunogenicity in a human, not to "
                             "whether the molecule expresses in bacteria. It is "
                             "grouped here because it is evidence about the "
                             "sequence rather than about binding, and it matters "
                             "only if a candidate is ever taken towards the "
                             "clinic.",
                         src="biophi"),
    "nearest_pdb": G(plain="The closest nanobody that somebody has already solved "
                           "a structure of.",
                     how="Structural search against SAbDab-nano with TM-align, "
                         "reporting the best match with its TM-score and "
                         "per-region sequence identity.",
                     why="A close experimental neighbour means the predicted "
                         "model is standing on something real, and it gives a "
                         "published starting point for expression conditions and "
                         "constructs. Only derived values are stored here; the "
                         "databases themselves are not redistributed.",
                     src="sabdab_nano"),
    "experimental_badge": G(plain="A near-identical sequence already has a solved "
                                  "structure in the Protein Data Bank.",
                            how="Set when the nearest experimental neighbour "
                                "passes both a structural and a sequence identity "
                                "threshold.",
                            why="It is a positive signal only. Its presence says "
                                "this shape has been made and crystallised by "
                                "somebody. Its absence says nothing at all, "
                                "because most sequences that express perfectly "
                                "well have never been crystallised.",
                            src="sabdab_nano"),
}

del G

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
        "glossary": GLOSSARY,
    }


def _self_check() -> None:
    assert len(COLUMNS) == 48, f"the registry must hold 48 columns, found {len(COLUMNS)}"
    assert len(BY_KEY) == len(COLUMNS), "duplicate column key"
    assert len(DEFAULT_KEYS) == 6, f"six default columns, found {len(DEFAULT_KEYS)}"
    assert {c["group"] for c in COLUMNS} == {gid for gid, _ in GROUPS}
    # Every column carries a glossary entry, so the About tab's parameter table
    # cannot silently fall behind the registry when a column is added.
    missing = [c["key"] for c in COLUMNS if c["key"] not in GLOSSARY]
    assert not missing, f"no glossary entry for {missing}"
    stray = [k for k in GLOSSARY if k not in BY_KEY]
    assert not stray, f"glossary entry for a column that does not exist: {stray}"
    for key, entry in GLOSSARY.items():
        assert set(entry) == {"plain", "how", "why", "src"}, \
            f"glossary entry {key} has the wrong fields"


_self_check()

if __name__ == "__main__":
    import json

    print(json.dumps({"total": len(COLUMNS), "groups": group_counts(),
                      "defaults": DEFAULT_KEYS}, indent=2))
