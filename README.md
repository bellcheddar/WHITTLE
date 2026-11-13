# 🧪 WHITTLE

> **Two and a half thousand candidate nanobodies, five targets, cheapest checks first, and an honest account of the one axis that does not work.**

![python](https://img.shields.io/badge/python-3.12-3776AB?logo=python&logoColor=white) ![pixi](https://img.shields.io/badge/env-pixi%20%7C%20conda--forge-FF8C00) ![apple](https://img.shields.io/badge/compute-Apple%20silicon%20%7C%20float32-000000?logo=apple&logoColor=white) ![boltz](https://img.shields.io/badge/co--folding-boltz--community-3457A6) ![sklearn](https://img.shields.io/badge/Tm%20model-scikit--learn%20stacked%20GPR%2BSVR-F7931E?logo=scikitlearn&logoColor=white) ![flask](https://img.shields.io/badge/server-Flask-000000?logo=flask&logoColor=white) ![candidates](https://img.shields.io/badge/candidates-2%2C500-2F7D62) ![licence](https://img.shields.io/badge/licence-MIT-1C244B) ![author](https://img.shields.io/badge/author-Marc%20C.%20Deller,%20D.Phil.-1C244B)

<table>
<tr>
<td>🌐 <b>Live</b></td><td><a href="https://whittle.mdeller.com" target="_blank" rel="noopener noreferrer">whittle.mdeller.com</a></td>
<td>✉️ <b>Contact</b></td><td><a href="mailto:marc@marcdeller.com">marc@marcdeller.com</a></td>
<td>🐙 <b>GitHub</b></td><td><a href="https://github.com/bellcheddar/WHITTLE" target="_blank" rel="noopener noreferrer">bellcheddar/WHITTLE</a></td>
</tr>
</table>

---

<img src="app/static/img/workflow.svg" alt="The WHITTLE pipeline: 2,500 candidate sequences through Tier 0 sequence checks, Tier 1 structure and developability, the Tier 1b thermostability model and Tier 2 co-folding, into within-target percentiles and a Pareto front." width="100%">

**W**hole-repertoire **H**ierarchical **I**n-silico **T**riage: **T**hermostability, **L**iabilities, **E**xpression.

A design campaign hands you five hundred candidate nanobodies per target and no way to
choose between them. Ordering genes for all of them is not a plan, and the scores that came
with them are not comparable across targets. WHITTLE runs every candidate through the same
funnel, cheapest checks first: sequence liabilities in minutes, predicted structure and
developability in an hour, a thermostability model trained here rather than licensed, and
co-folding against the real antigen construct for the survivors. What comes out is about
twenty candidates per target on four scoring axes, each one a percentile **within its own
target**, because the five targets are different sizes, different constructs and different
libraries and a pooled league table would be meaningless.

The part most tools leave out is on the front page here: **the affinity axis is the weak
one.** Published benchmarks show co-folding confidence does not reliably tell a cognate
nanobody-antigen pair from a non-cognate one. So WHITTLE is strong at discarding and only
suggestive about ranking, it says so wherever that axis appears, and a scrambled-CDR
negative control is run through the same pipeline to show what the axis actually separates.
If it separates nothing, that is reported as a finding and the axis stops driving the
ranking.

## 🔧 How it works

| Tier | What it computes | Cost |
|---|---|---|
| **0 · sequence** | IMGT numbering and CDRs, the camelid hallmark tetrad, liability motifs, charge and pI, disulfides, CDR3 clustering | minutes, all 2,500 |
| **1 · structure** | NanoBodyBuilder2 models with OpenMM refinement, TNP developability, surface patches, SWI and Protein-Sol, AbLang2 nativeness, OASis humanness, nearest experimental neighbour | about an hour per target |
| **1b · thermostability** | A stacked Gaussian-process and support-vector model over one-hot, VHSE and ESM-2 features, trained on 650 curated apparent-Tm nanobodies | one training run |
| **2 · co-folding** | boltz-community against the real antigen construct, cheapest target first, five samples per complex, PRODIGY interface energetics and epitope binning | days, top 20 per target |
| **ranking** | Within-target percentiles, Pareto non-dominated fronts, Derringer-Suich desirability, and the axis correlation matrix published beside them | seconds |

### Things this build refuses to do

- **No pooled ranking.** Five targets, five independent scales, a target selector rather
  than an "all targets" view. The spec says it, the scoring enforces it, and a second
  independent implementation checks it.
- **No fabricated columns.** `target_cdr_identity` exists for EGFR and TNFa and is left
  null for the other three, because those are the only two targets with known nanobodies.
  The two CDR columns that are empty in every source file are ingested as null.
- **No silent deduplication.** Ten sequences were submitted against more than one target.
  They are flagged and used as internal promiscuity controls.
- **No academic-only components.** Every tool is permissively licensed and self-hosted.
  NanoMelt, CamSol and PyRosetta are not installed. nanoBERT would have been the
  better-matched nativeness model and is CC-BY-NC-SA, so AbLang2 is used instead and the
  cost of that choice is stated on the About tab.
- **No request-time compute.** The droplet gets a read-only SQLite file, some JSON and the
  gzipped models. Everything is computed on the Mac Studio.

## 📊 What the run found

Full detail in **[FINDINGS.md](FINDINGS.md)**, generated from the database rather than
written by hand. The headlines:

- The upstream `hallmark` column implements an **at-least-three-of-four** rule over the
  camelid tetrad, not the four-of-four it appears to. Both calls are reported.
- **Framework identity tells you nothing.** Against the nearest experimental nanobody,
  framework identity has a median of 93 percent and the TM-score a median of 0.92, for
  essentially every candidate. CDR identity, median 60 percent, is the number that carries
  information. Only 2 of 2,500 candidates are close enough to a solved structure to earn an
  experimental badge.
- The thermostability model's **grouped** cross-validated R squared is 0.39 against 0.70
  row-wise. That 0.31 gap is the near-duplicate leak a row-wise split hides.
- On genuinely held-out sequences the model **orders** nanobodies well (Spearman 0.67) and
  **does not place them on the temperature scale** (R squared below zero). That is why the
  axis is a percentile and the degrees carry an error bar. See
  **[MODEL_CARD.md](MODEL_CARD.md)**.

## 🚀 Running it

```bash
pixi install                       # conda-forge and bioconda, osx-arm64
pixi run probe                     # writes config/hardware.toml
pixi run ingest && pixi run tier0  # 2,500 rows, minutes
pixi run fold                      # NanoBodyBuilder2, about 20 minutes on 12 workers
pixi run tnp                       # TNP developability, one call per target
pixi run neighbours                # SAbDab-nano nearest neighbour and TM-align
pixi run tm-data && pixi run tm-train && pixi run tm-predict
pixi run score && pixi run build-db
pixi run serve                     # http://127.0.0.1:8021

pixi run -e cofold install-boltz   # the community fork, never upstream boltz
pixi run antigens
./tools/run_queue.sh --smoke       # one GLP1R complex, timed, before anything else
./tools/run_queue.sh               # the full queue; days, resumable, under caffeinate
./tools/status.sh --watch          # progress from another shell

pixi run qc && pixi run qc-pareto && pixi run qc-vacuous
bash deploy/deploy.sh
```

Every stage writes `data/manifests/<stage>.jsonl` and skips what is already done, because
Tier 2 runs for days and will be interrupted.

## 📁 Layout

| Path | What is in it |
|---|---|
| `pipeline/` | Ingest, the four tiers, scoring, the Tm model, the control set |
| `app/` | Flask app: Candidates, Control set, About |
| `tools/` | QC, the independent Pareto check, the vacuous-axis detector, the queue, the secret scan |
| `config/` | `targets.yaml` and `references.yaml` drive everything; nothing is hard-coded |
| `deploy/` | nginx, systemd and the rsync deploy to the droplet |
| `docs/` | The Mermaid workflow source |

`data/` and `models/` are gitignored in full. No reference database is mirrored into this
repository or onto the droplet: only derived columns ship, and each source is cited with its
own licence and snapshot date on the About tab.

## ⚖️ Licence and citation

Code MIT. Every component is permissively licensed; the About tab lists each with its
version, its licence and a DOI-verified citation.

> Deller, M. C. *WHITTLE: whole-repertoire hierarchical in-silico triage of nanobody
> candidates (thermostability, liabilities, expression).* https://whittle.mdeller.com

**Nothing on this site is a measurement.** Every value is predicted or estimated.
