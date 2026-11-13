# WHITTLE: Build Specification

**v1.0** · Whole-repertoire Hierarchical In-silico Triage: Thermostability, Liabilities, Expression

Project root: `/Users/dellboy/Documents/Vibe_Coding/WHITTLE`
Live at: `https://whittle.mdeller.com`
Repo: `https://github.com/bellcheddar/WHITTLE` (public, created by this build)

---

## 0. Autonomy directive

### 0.1 Build to completion

1. **Run end to end without asking.** Pick the most defensible option, write one line to `DECISIONS.md` (decision, alternative, reason, reverse-by), continue. No plan-approval pauses, no confirmations, no "shall I proceed".
2. **Ship early, then improve.** Deploy at the end of every phase, including Phase 1 with only sequence-level metrics. A live app with four columns beats a perfect one on localhost. Each deploy is `deploy/deploy.sh`, run by the build itself.
3. **Reuse before writing.** Section 3 is the reuse map. Grep the siblings first; write new code only where nothing fits. Record every lift in `REUSE.md` (source path, destination, what changed). Never import across project folders. Never modify a sibling.
4. **Fallback ladder for any failure:** retry once with a fix, then the documented alternative, then the next-best tool, then degrade (skip the item, write `status = "failed:<reason>"` to its row, continue). One sequence, one fold, one optional service never blocks the build.
5. **Resumable.** Every stage writes `data/manifests/<stage>.jsonl`; on restart read `PROGRESS.md` and skip done items. Assume the process is killed at any moment. Tier 2 is a multi-day queue: it MUST survive interruption and resume mid-target.
6. **Log as you go.** `BUILD_LOG.md` (timestamped lines) and `PROGRESS.md` (phase checklist with percentages, rewritten after every phase and at least every 30 minutes during compute).
7. **Never idle.** While folds or co-folds run, build the next phase's app surface against whatever rows are already complete; backfill when manifests close.
8. **Commit and push at the end of every phase** (`phase N: <summary>`). Pushing to the public repo is pre-approved. See 0.4 for what must never be committed.
9. **Stay inside the project root** except package caches, the droplet paths in `config/deploy.toml`, and read-only reads of sibling projects for reuse. Never `sudo` locally.
10. **Keep the Mac awake** with `caffeinate -dimsu` around any step over ten minutes. Tier 2 needs it for days.
11. **House style:** British English, no em dashes (colons or parentheses), none of groundbreaking, revolutionary, paradigm-shifting, game-changing, cutting-edge, unprecedented, seamless, leverage (verb), delve. App copy never overstates a prediction: "predicted", "estimated", "ranked", never "confirmed" or "validated".

### 0.2 No approval gates. Stop only when genuinely blocked.

There are no review gates in this build. Deployment, repo creation, pushing and overwriting your own prior output are all pre-approved. Stop and ask Marc **only** when a step has no fallback at all:

| Blocker | Opens only when | Needs from Marc |
|---|---|---|
| **B1 Credentials** | The droplet is unreachable and no usable SSH target or key is found anywhere under Section 4.4 | SSH host, user, key path |
| **B2 Disk** | The project volume is under 40 GB free and nothing under `data/` is safe to prune | What to prune |
| **B3 Science** | A positive control from Section 9 fails in a way that would invert a headline ranking, after the full fallback ladder | "Ship with the provisional banner" or "change the method" |

Everything else is a fallback, never a stop:

| Missing or broken | Fallback (automatic) |
|---|---|
| ANARCI will not build (HMMER) | `anarcii` (PyTorch, no HMMER). If both fail, derive CDRs by regex against the VHH framework consensus and mark `numbering_method = "regex_fallback"` |
| OpenMM refinement fails for a model | Keep the unrefined NanoBodyBuilder2 model, flag `refined = 0` |
| A NanoBodyBuilder2 fold fails | Retry once; then ESMFold; then skip the structural axes for that row and keep its sequence axes |
| Boltz co-fold OOMs | Drop to 1 seed, then `--no_kernels`, then skip that complex and log it. Never silently substitute a monomer for a trimer |
| ColabFold MSA server unreachable | Use the cached antigen MSA; if there is none for that target, defer that target to the end of the Tier 2 queue and continue with the others |
| Aggrescan3D unavailable | SWI plus Protein-Sol formulae alone; note the missing term in the About tab |
| Crossref unreachable for a reference | Keep the hand-written fallback title in `config/references.yaml`, flag `doi_verified = 0`, show it unverified |
| Playwright Chromium absent | Skip screenshots, log degraded, never `playwright install` |

**Blocker protocol:** write `BLOCKED.md` (which blocker, exactly what is needed, the command that unblocks it, what is already done, what happens next); notify with
`osascript -e 'display notification "<one line>" with title "WHITTLE blocked <id>" sound name "Glass"'`;
then keep working on everything the blocker does not block. The expected end state of a clean run is **no `BLOCKED.md`, the repo public, and the app live at whittle.mdeller.com**.

### 0.3 Launch recipe

```bash
cd /Users/dellboy/Documents/Vibe_Coding/WHITTLE
./launch.sh          # tmux session "whittle"; Claude Code reads this spec and executes it
```

`.claude/settings.json` and `CLAUDE.md` are already in the folder.

### 0.4 Secret hygiene (the repo is PUBLIC)

This is the one rule with no fallback and no exception.

1. **Never commit a credential.** Not in code, not in a config, not in a log, not in a manifest, not in `BUILD_LOG.md`, not in a commit message, not in a test fixture.
2. `config/deploy.toml` holds the droplet host, user and key path. It is **gitignored**. `config/deploy.toml.example` is committed with placeholder values.
3. The droplet SSH target is **read at deploy time** from the sibling AlphaFraud `.env` (Section 4.4) or from `config/deploy.toml`. Never copy `.env` into this project. Never cat a key into a file the repo tracks.
4. Before the first `git push`, and again before the Phase 4 push, run `tools/secret_scan.py`: grep the full tree and the git history for `BEGIN .*PRIVATE KEY`, `ssh-rsa`, `ssh-ed25519`, `api[_-]?key`, `secret`, `password`, `token`, `AKIA`, `ghp_`, `github_pat_`, and any literal IPv4 address. Any hit outside `.gitignore`d paths or this spec's own placeholder text **aborts the push** and writes `BLOCKED.md`.
5. If a secret is ever committed, do not rewrite history and carry on quietly: stop, write `BLOCKED.md`, and tell Marc the key needs rotating.
6. `.gitignore` is written in Phase 1 **before** `git init` is followed by any `git add`.

---

## 1. Goals

1. **Triage 2,500 candidate nanobody sequences** (500 against each of five targets) down to a defensible shortlist of about 20 per target, cheapest checks first.
2. **Four scoring axes per candidate:** binding, non-specific binding, thermostability, expression. Reported as percentiles **within each target**, never pooled across targets.
3. **Honest about the weak axis.** Published benchmarks show co-folding confidence does not reliably discriminate cognate from non-cognate nanobody-antigen pairs. WHITTLE is a de-risking tool: strong at discarding, suggestive about ranking. That statement appears on the app, not only in the docs.
4. **Keep the upstream numbers visible.** The source CSVs are already scored. Those columns ship as `upstream_*` beside recomputed values so the two can be compared, with the known defects flagged (Section 5.1).
5. **Every tool and database cited.** A DOI-verified reference table and a provenance record (tool versions, database snapshot dates) on the About tab.
6. **An honesty tab.** Known nanobodies with published measurements run blind through the same pipeline, with observed against predicted plotted per axis.
7. **Licence-clean.** Every component permissively licensed and self-hosted. No academic-only binaries, no webserver dependencies at serve time, no PyRosetta, no NanoMelt, no CamSol.
8. **Public and citable.** Public GitHub repo, public site, deep-linkable rows, URL-encoded view state.

---

## 2. Hardware and placement

| Where | What | Role |
|---|---|---|
| **Mac Studio, M1 Ultra, 128 GB** | All compute: folding, TNP, Tm model, co-folding | Everything in Sections 5 and 6 runs here |
| **MacBook Pro** | Currently holds `shortlist/` and `data/` | Source data lives here now: see 2.1 |
| **Droplet, 8 GB / 80 GB** | nginx + gunicorn + Flask, read-only | Displays finished results only. No torch, no ML stack, no request-time compute |

### 2.1 First action of Phase 1

The shortlist CSVs are on the MacBook Pro. Confirm the project root is on the Studio before any compute. If `sysctl -n machdep.cpu.brand_string` plus `system_profiler SPHardwareDataType` do not show an M1 Ultra, write one line to `BUILD_LOG.md`, continue with Phase 1 (which is cheap), and resolve the sync before Phase 2 by rsyncing the project root to the Studio. Do not run Tier 1 or Tier 2 on the laptop.

### 2.2 Apple Silicon constraints (M1 generation, no hardware bf16)

- Install everything through **miniforge**, native arm64. One Rosetta-translated dependency poisons the environment.
- Export `PYTORCH_ENABLE_MPS_FALLBACK=1` in every compute shell. Unsupported MPS ops (`aten::linalg_svd` among them) must fall back to CPU rather than crash at sequence 1,900.
- **float32 everywhere.** bf16 on Metal needs M2 or later; `bf16-mixed` is measurably *slower* on M1. Never let an upstream default re-enable it.
- OpenMM from conda-forge, CPU platform. At ~125 residues Metal buys nothing.
- `tools/hwprobe.py` writes `config/hardware.toml`; all batch sizes and worker counts come from there. Nothing hard-coded.

---

## 3. Reuse map (read before writing anything)

Siblings live under `/Users/dellboy/Documents/Vibe_Coding/`. Copy and adapt; never import across folders; never modify a sibling. Record each lift in `REUSE.md`.

| Need | Source | Destination |
|---|---|---|
| Droplet deploy by rsync then service restart, nginx vhost, systemd web service, provision script, release script | `AlphaFraud/deploy/*` | `deploy/*` |
| Droplet SSH target, user and `/opt/<app>` layout | `AlphaFraud/.env` (read only, never copied) | `config/deploy.toml` (gitignored) |
| gunicorn config for a read-only workload | `debris/deploy/gunicorn.conf.py` | `deploy/gunicorn.conf.py` |
| nginx static-asset block, long-cache snippet, gzip types | `debris/deploy/nginx.debris.conf` | `deploy/nginx.whittle.conf` |
| Deploy script with `--data-only` data swap | `debris/deploy/deploy.sh` | `deploy/deploy.sh` |
| Flask tab layout, Tabulator table, detail drawer, `/healthz`, vendored Plotly with `?v=<mtime>` cache busting | `CODSWALLOP`, `debris/app/`, `GOBSMACKED/app/` | `app/` |
| Flask app factory, `db.py`, `config.py`, `routes/`, `templates/base.html` shape | `GOBSMACKED/app/*` | `app/*` |
| Boltz campaign driver, YAML input generation, chain-id namespacing, seed handling, affinity parsing, resumable campaign queue | `BoltzMaker/landlord/campaign.py`, `build.py`, `validate.py`, `config.py` | `pipeline/cofold.py`, `pipeline/boltz_input.py` |
| Boltz runtime install and fallback ladder | `BoltzMaker/landlord/install.py`, `fallback.py` | `pipeline/boltz_env.py` |
| PLIP / interface-contact labelling from a predicted complex | `BoltzMaker/.plip_env/cif2plip`, `label_residues.py`; `GOBSMACKED/app/services` | `pipeline/interface.py` |
| Fold → analyse bundle shape, run manifest, per-run directory convention | `GOBSMACKED/bundle_template/*` | `data/runs/<target>/` |
| Hardware probe and derived tuning | `KINFOLK/tools/hwprobe.py` | `tools/hwprobe.py` |
| pixi workspace shape, task table, activation env | `KINFOLK/pixi.toml`, `BoltzMaker/pixi.toml` | `pixi.toml` |
| Manifest read/write, resumable-stage pattern, one JSONL per stage | `debris/pipeline/common.py` | `pipeline/common.py` |
| tenacity-wrapped HTTP with a process-wide rate limiter, polite-pool User-Agent | `debris/pipeline/common.py` (`http_get`, `RateLimiter`) | `pipeline/common.py` |
| Crossref DOI verification with title matching | `debris/pipeline/verify_doi.py` | `pipeline/references.py` |
| SQLite schema module and read-only serve pattern | `debris/pipeline/schema.py`, `GOBSMACKED/app/db.py` | `pipeline/schema.py`, `app/db.py` |
| QC harness and smoke test shape | `debris/tools/smoke_test.py`, `qc_ui.py` | `tools/smoke_test.py`, `tools/qc_ui.py` |
| Findings writer | `debris/tools/findings.py`, `NOBBLE/tools/findings.py` | `tools/findings.py` |
| Composite scoring cross-validation, matched stratified sampling | `NOBBLE/tools/composite_crossval.py`, `matched_arm_stratified.py` | `pipeline/validate_scores.py` |
| Vacuous-test detector (guards against a metric that cannot fail) | `NOBBLE/tools/vacuous_tests.py` | `tools/vacuous_tests.py` |
| Tier-queue runner and status reporting | `NOBBLE/tools/run_queue.sh`, `tier2.sh`, `status.sh` | `tools/run_queue.sh`, `tools/status.sh` |
| Disk headroom check before a long run | `NOBBLE/tools/space_check.py` | `tools/space_check.py` |
| Mol* viewer embed, lazy init, dispose on close | `GOBSMACKED/app/static/js`, `chatPDB-web` | `app/static/js/viewer.js` |
| Icon house style and favicon set | `marcs-vibe-icon` skill | `icons/` |
| README house standard | `marcs-vibe-coding` skill (`readme_forge`), `debris/tools/readme.py` | `README.md` from `README.source.md` |
| Reuse log and decision log conventions | `debris/REUSE.md`, `debris/DECISIONS.md` | `REUSE.md`, `DECISIONS.md` |

**Grep BoltzMaker first for anything co-folding related.** It is a campaign-scale Boltz driver and this build is a campaign. Do not rewrite its queue, its input generation or its chain-id handling.

---

## 4. Environment

### 4.1 Python

pixi workspace, conda-forge, osx-arm64. `pixi.toml` shaped after `KINFOLK/pixi.toml` with tasks: `ingest`, `tier0`, `tier1`, `tm-train`, `tier2`, `controls`, `refs`, `build-db`, `serve`, `deploy`, `qc`, `readme`.

Install order matters: **OpenMM and HMMER from conda-forge first**, then pip for the rest.

### 4.2 Components and licences

| Axis | Tool | Licence |
|---|---|---|
| Numbering, CDR definition | ANARCI, else `anarcii` | BSD-3 |
| Structure prediction | NanoBodyBuilder2 (ImmuneBuilder) | BSD-3 |
| Non-specific binding, developability | TNP (`oxpig/TNP`, `pip install .`, takes a FASTA directly) | BSD-3 |
| Thermostability | **Own model**: Ridge-stacked GPR + SVR over one-hot, VHSE and ESM-2 embeddings, trained on the public NbThermo set (Section 6.3) | ours; scikit-learn BSD, ESM-2 MIT |
| Thermostability cross-check | TemStaPro | open |
| Structural stability | ThermoMPNN **only if its repo LICENSE is permissive**: check first, drop the axis if not | check at build time |
| Solubility, aggregation | SWI and Protein-Sol formulae, Aggrescan3D | public formulae, GPL-3 |
| Nativeness | AbLang2 or nanoBERT pseudo-perplexity | BSD-3 / open weights |
| Humanness | BioPhi OASis | MIT |
| Co-folding, affinity prior | `boltz-community` (Apple Silicon MPS fork) | MIT |
| Interface energetics | PRODIGY (`prodigy-prot`) | Apache-2.0 |
| Comparison, clustering | MMseqs2, US-align, FreeSASA, Biopython | open |

**Explicitly excluded, do not install, do not substitute:** NanoMelt (Cambridge Enterprise licence), CamSol (webserver-only academic registration), PyRosetta (no third-party service use), upstream `boltz` (declares CUDA-only `cuequivariance` and will not install on macOS).

### 4.3 Boltz invocation

```
boltz predict <input.yaml> --accelerator mps --out_dir <dir>
```
with `PYTORCH_ENABLE_MPS_FALLBACK=1` and float32. Install without the `[cuda]` extra. Verify at the start of Phase 3 that `boltz` resolves to the community fork, not upstream.

### 4.4 Droplet credentials

Resolution order, stopping at the first that works:

1. `config/deploy.toml` in this project, if present.
2. The sibling `AlphaFraud/.env` (the host and user recorded there; `debris/REUSE.md` documents this as the house source). Read it, write the host, user and key path into `config/deploy.toml`, and **never copy the file or its contents into anything the repo tracks**.
3. `debris/config/deploy.toml`, `GOBSMACKED/deploy/deploy.sh` and `CODSWALLOP` for the same host under a different app name.
4. `~/.ssh/config` for an `mdeller` or droplet host alias.

If all four fail, that is blocker **B1**.

### 4.5 Port

Ports already taken on the droplet: 8000 (AlphaFraud), 8006 (CODSWALLOP), 8009 (GOBSMACKED), 8014 (debris), and others not visible from the laptop. **Probe before choosing:** `ssh <host> "ss -ltnp | grep -oE ':80[0-9][0-9]'"` and pick the lowest free port at or above 8021. Record it in `config/deploy.toml` and in `deploy/nginx.whittle.conf`. Do not assume a port is free.

---

## 5. Input data

### 5.1 The shortlist CSVs

`shortlist/{EGFR,GLP1R,PRLR,TL1A,TNFa}.csv`, 500 rows each, 2,500 total, 2,487 unique sequences, 106-133 aa. Known properties, all verified on 2026-11-13:

1. **The upstream composite `score` is not comparable across targets.** EGFR and TNFa carry two extra columns (`s_his`, `target_cdr_identity`) and blend four terms; the other three blend three. Ranges differ accordingly (EGFR 0.61-0.88, TNFa 0.61-0.92, others 0.77-0.95). Ingest it as `upstream_score`, and wherever it is displayed attach the note "not comparable between targets". Never sort the whole table by it.
2. **`upstream_tm_pred` is uncalibrated.** Range 71-101 degrees C. Apparent Tm above 100 degrees C is not a measurement a DSF run produces for a VHH in solution. Keep it as `upstream_tm_pred`, display it beside our own prediction, and say so on the About tab.
3. **`cdr_novelty` and `cdr_identity` are empty in all five files.** Ingest as NULL. We compute our own equivalents.
4. **`target_cdr_identity` exists only for EGFR and TNFa**, the only two targets with known nanobodies in SAbDab-nano. The target-aware comparison axis and the control set are therefore limited to those two. Do not fabricate the column for the other three: leave it NULL and show "no target-matched reference" in the UI.
5. **10 sequences appear under more than one target.** Do **not** deduplicate. Flag them `multi_target = 1` and surface them on the control tab as internal promiscuity controls: a candidate scoring well against two unrelated targets is either polyreactive or evidence the scoring is not target-aware.
6. **His enrichment is in the libraries, not the scoring.** EGFR and TNFa are 500/500 with at least one CDR histidine (mean 1.37 and 1.32); the other three sit at 15-18 percent, the background rate. This was an intentional pH-switch design for those two targets only (confirmed 2026-11-13). Carry `his_cdr_total`, `his_cdr1/2/3` for all five, but score `s_his` only for EGFR and TNFa, and label it in the UI as a design criterion, not a quality metric.

### 5.2 Targets

`config/targets.yaml` is authoritative and already in the folder. Boundaries are UniProt numbering, accessed 2026-11-13. Do not re-derive them, do not round them, do not substitute a PDB construct without writing a DECISIONS.md line.

| Target | UniProt | Construct | Range | Antigen | Complex | Tier 2 cost |
|---|---|---|---|---|---|---|
| GLP1R | P43220 | N-terminal ECD only | 24-139 | 116 aa | ~241 aa | lowest |
| PRLR | P16471 | ECD | 25-234 | 210 aa | ~335 aa | moderate |
| TNFa | P01375 | soluble homotrimer | 77-233 x3 | 471 aa | ~596 aa | high |
| TL1A | O95150 | secreted homotrimer | 72-251 x3 | 540 aa | ~665 aa | high |
| EGFR | P00533 | ECD | 25-645 | 621 aa | ~746 aa | highest |

Rules that follow from the table:

- **TNFa and TL1A are co-folded as trimers.** Three identical chains. The epitope commonly spans two subunits, so a monomer run would be wrong, not merely cheaper. Report interface metrics **per subunit and in total**, so a candidate binding across the subunit cleft is distinguishable from one binding a single face. If a trimer will not fit in memory, reduce seeds, then skip: never fall back to a monomer.
- **GLP1R is ECD only, by explicit decision.** The 7TM bundle and the extracellular loops at 202-227, 291-305 and 371-383 are excluded. Every GLP1R row and the GLP1R target header carry a visible "ECD-scoped" label. The app must never let a GLP1R result read as a claim about the intact receptor in a membrane.
- **EGFR is the cost driver**, roughly an order of magnitude above GLP1R per co-fold because attention scales about quadratically with length. Full ECD for this run. If Tier 2 wall-clock exceeds four days, the documented cut is to scope EGFR to domain III (approximately 311-480 in mature numbering, the usual nanobody epitope); take it, write a DECISIONS.md line, and label the EGFR rows accordingly.

### 5.3 Databases

| Source | Use | Note |
|---|---|---|
| SAbDab-nano | Nearest experimental neighbour, TM-score, per-region identity | Snapshot date recorded; derived values only in the repo |
| PLAbDab-nano | Patent and literature sequence context | Derived values only |
| NbThermo | Tm training labels (Section 6.3) | The curated set is published and freely available |
| PDB | Experimental-structure badge | CC0, so PDB-derived values ship without restriction |
| UniProt | Antigen sequences for the constructs in 5.2 | Fetch once, cache under `data/antigens/` |

**Do not mirror SAbDab-nano or PLAbDab-nano into the repo or onto the droplet.** Store only derived columns (nearest PDB id, TM-score, identity, same-target hit) and link out. Record each snapshot date in `config/about.toml` for the provenance table.

---

## 6. Methods, by tier

Tiers run cheapest first. Each tier writes its own manifest and is independently resumable.

### 6.1 Tier 0: sequence level, all 2,500, minutes

1. Numbering and CDR definition (ANARCI or `anarcii`), IMGT. Region-resolved sequence per row.
2. CDR1/2/3 lengths against the IMGT distribution; V-gene assignment.
3. **VHH hallmark tetrad** at positions 42, 49, 50, 52: intact, or VH-type hydrophobics present (a solubility and stability liability). The source files claim `hallmark = True` for 2,497 of 2,500: recompute and reconcile, do not trust it.
4. **Sequence liability regex pass:** unpaired cysteine, N-glycosylation sequon (N-X-S/T), NG and NS deamidation, DG and DP isomerisation, surface Met and Trp oxidation motifs.
5. Net charge at pH 7.4, isoelectric point. A pI near 7.4 is an aggregation flag; a very high pI is a PK and polyreactivity flag. Both are flags, not scores.
6. Canonical Cys22-Cys96 pair present; any additional CDR1-CDR3 disulfide noted.
7. His content per CDR, recomputed (Section 5.1 item 6).
8. Clustering: MMseqs2 (or CD-HIT) on CDR3 at 90 percent identity, within each target. Cluster id per row. The shortlist must be diverse across clusters, not 20 variants of one family.
9. `multi_target` flag for the 10 cross-target sequences.

Tier 0 discards nothing on its own. It populates flags; the Pareto step decides.

### 6.2 Tier 1: structure level, all survivors, about an hour per target

1. **NanoBodyBuilder2** model per sequence, OpenMM refinement on the CPU platform. Per-residue error estimate retained. Write `.pdb.gz`.
2. **TNP** over each target's FASTA in one call (`TNP --file <fasta> --output <dir>`). All five TNP metrics plus the clinical-stage comparison. This is the whole non-specific-binding axis.
3. Surface hydrophobic and charged patches from the refined model (FreeSASA plus the TNP patch definitions).
4. Solubility: SWI and Protein-Sol formulae; Aggrescan3D on the model if available.
5. Nativeness: AbLang2 or nanoBERT pseudo-perplexity. Odd predicted sequences score badly for a reason, and this is the best generic sanity filter.
6. Humanness: BioPhi OASis.
7. **Nearest experimental neighbour** in SAbDab-nano: US-align TM-score plus per-region identity. Report per region, not whole-chain: framework TM will be about 0.9 for everything and tells you nothing.
8. Experimental-structure badge where a near-identical sequence has a PDB entry. **Positive signal only.** Absence says nothing about expressibility and must not reduce a score.

### 6.3 Tier 1b: the thermostability model (ours)

Build it, do not license one.

1. Pull the public NbThermo set plus the published curated apparent-Tm dataset. Expect roughly 600 unique sequences with apparent Tm.
2. Features: one-hot, VHSE descriptors, ESM-2 embeddings (650M; drop to 150M if memory demands it).
3. Model: `sklearn.ensemble.StackingRegressor`, Ridge meta-learner over Gaussian Process and Support Vector regressors, one base model per feature block.
4. **Validation is the part that matters.** Repeated stratified nested cross-validation, **grouped by sequence similarity**: cluster the training set at 90 percent CDR3 identity and split on clusters, never on rows. A row-wise split will report a flattering R squared that evaporates on 2,500 predicted sequences.
5. Report per-row **error estimate** from the distance to the nearest training neighbour. A candidate far from the training manifold gets a wide bar, not a crisp number. The UI shows the bar.
6. Cross-check with TemStaPro; report both; if they disagree systematically, say so on the About tab rather than picking a winner.
7. Write `MODEL_CARD.md`: training set, size, split strategy, cross-validated metrics, known failure modes, and the explicit statement that this is an apparent-Tm estimate and not a measurement.

### 6.4 Tier 2: co-folding, top 20 per target, days

1. Build antigen constructs from `config/targets.yaml` (Section 5.2). Cache each antigen's MSA **once per target** and reuse it across all candidates for that target: that single decision removes most of the MSA cost. Nanobody MSAs add little since the framework is near-germline.
2. Run order **cheapest first: GLP1R, PRLR, TNFa, TL1A, EGFR.** A broken co-folding setup then surfaces on a 241 aa complex in minutes rather than on EGFR overnight.
3. Five seeds per complex. Retain ipTM, pTM, pLDDT, PAE and the per-interface confidences.
4. PRODIGY interface delta G; buried surface area; contact map. For the trimers, per-subunit and total.
5. Epitope binning across candidates per target: cluster by interface residue set, so "these 6 candidates hit the same face" is visible.
6. For EGFR and TNFa only: compare the predicted epitope against known nanobody complexes for the same target from SAbDab-nano.
7. **Queue discipline.** `tools/run_queue.sh` with per-complex manifest rows; `caffeinate -dimsu`; resume mid-target on restart; `tools/status.sh` prints progress. Roughly 100 co-folds across three large antigens: this is a multi-day job and will be interrupted.

### 6.5 Ranking

1. Convert each axis to a **percentile within its own target**. Five independent scales.
2. **Pareto non-dominated sorting** per target across the four axes. Show the front. Nanobody engineering is full of genuine trade-offs (a long stabilised CDR3 against polyreactivity; high charge against solubility and PK) and a weighted sum hides exactly the information the user needs.
3. A **Derringer-Suich desirability** function on top, with user-set weights in the UI, for anyone who wants a single ordering.
4. **Guard against double counting.** Hydrophobic patch area feeds non-specificity, aggregation and expression. Compute the axis correlation matrix per target, publish it on the QC panel, and if two axes correlate above 0.8 say so beside the ranking instead of pretending they are independent.
5. Never emit a global cross-target ranking. The app has a target selector, not a pooled table.

---

## 7. Data model

SQLite, built on the Studio, shipped read-only. One file, `data/whittle.sqlite`, plus `app/static/data/*.json` for anything the page loads directly.

```
targets(target PK, uniprot, name, construct, range_start, range_end,
        chain_length, stoichiometry, his_enriched, ph_switch_intended,
        known_nanobodies, ecd_scoped, antigen_md5, notes)

candidates(id PK, target FK, sequence, seq_md5, length, multi_target,
           numbering_method, v_gene, cdr1, cdr2, cdr3,
           cdr1_len, cdr2_len, cdr3_len, cdr3_cluster,
           hallmark_intact, hallmark_recomputed,
           his_cdr_total, his_cdr1, his_cdr2, his_cdr3,
           net_charge_74, pi, canonical_disulfide, extra_disulfide,
           -- upstream, comparison only
           upstream_tm_pred, upstream_liability_score, upstream_on_nll,
           upstream_score, upstream_s_tm, upstream_s_liab, upstream_s_nll,
           upstream_s_his, upstream_target_cdr_identity,
           upstream_global_cdr_identity, upstream_global_nearest)

liabilities(candidate_id FK, kind, region, position, detail)

structures(candidate_id PK, path, refined, mean_error, cdr3_error,
           psh, ppc, pnc, sfvcsp_na, tnp_flags, tnp_version,
           swi, protein_sol, aggrescan3d, nativeness_ppl, oasis_humanness,
           nearest_pdb, nearest_tm, nearest_identity_fw, nearest_identity_cdr,
           has_experimental_badge)

thermostability(candidate_id PK, tm_pred, tm_err, train_nn_distance,
                temstapro, thermompnn_ddg, model_version)

cofolds(candidate_id FK, target FK, seed, iptm, ptm, plddt, pae,
        prodigy_dg, bsa_total, bsa_per_subunit_json, epitope_bin,
        interface_residues_json, path, status)

scores(candidate_id PK, target FK, axis_affinity, axis_specificity,
       axis_thermo, axis_expression, pareto_rank, on_front,
       desirability, tier_reached)

controls(id PK, target, name, source_doi, pdb_id,
         observed_kd, observed_tm, observed_expression,
         predicted_affinity, predicted_thermo, predicted_specificity,
         predicted_expression, blind)

references(key PK, kind, title, authors, year, container, doi,
           doi_verified, url, repo, commit, licence, version_used,
           snapshot_date, role)

provenance(key PK, value)   -- tool versions, db snapshot dates, run id, git sha
```

Rules: the droplet gets the SQLite file and the gzipped structures, nothing else. No request-time computation. Everything the front page needs is precomputed into JSON.

---

## 8. The app

Flask, three tabs, served read-only. `DESIGN.md` in this folder is authoritative for layout and visual detail; this section is the functional contract.

### 8.1 Tab 1: Candidates

- **Target selector first.** Five targets, one at a time. There is no pooled view. The selected target's construct and any scoping caveat (GLP1R ECD-only) show in the header.
- **Four-axis explorer** pinned above the table: Plotly `parcoords`, one axis per criterion, Pareto-front candidates drawn over dominated ones, brushing wired to filter the table. Front and dominated differ in **both colour and line weight**, never colour alone.
- **Table**: Tabulator with virtual DOM rendering. Default 6 visible columns; the other 35 behind a grouped column chooser (Identity, Binding, Specificity and liabilities, Thermostability, Expression and evidence). Freeze the id column. Percentile cells carry a colour scale. 500 rows per target is small; 500 x 41 cells in a naive DOM table is not, hence Tabulator.
- **Row click opens a docked drawer**, not an overlay: Structure, Metrics, Sequence, Neighbours. One `/api/candidate/<id>.json` fetch on open. In the drawer, the **nearest-neighbour superposition** (model over nearest experimental structure, with TM-score and per-region identity) is the primary panel, and the Tm panel shows this candidate against the full distribution with its error bar.
- **Mol\* discipline.** Browsers cap WebGL contexts at about 16. Keep **one** plugin instance and reload structures into it, or call `plugin.dispose()` on drawer close. Without this the viewer dies silently part-way through a browsing session and looks like a random failure.
- **Lazy-load.** Mol\* (about 4 MB) and Plotly must not be in the initial bundle. Build a custom partial Plotly bundle with `scatter`, `bar` and `parcoords` only (`parcoords` is absent from `plotly.js-basic-dist-min`). Dynamic-import Mol\* on first drawer open.
- **Deep links.** `#/t/<target>/nb/<id>` opens that drawer directly. Column selection, brush state and the Pareto toggle persist in the URL, so a configured view is shareable and citable.

### 8.2 Tab 2: Control set

The honesty tab. A first-class page, not a footnote.

- Known nanobodies with published Kd, Tm or expression yield, run **blind** through the same pipeline.
- Observed against predicted scatter per axis, with R squared and a residual plot.
- The 10 `multi_target` sequences as internal promiscuity controls.
- **Coverage caveat, stated plainly:** only EGFR and TNFa have target-matched references, so this tab validates the general axes everywhere but the target-aware axis only for those two. TL1A, PRLR and GLP1R have no published nanobody benchmark.
- Same table component as Tab 1, different data file.

### 8.3 Tab 3: About and methods

- **Workflow diagram.** Author as Mermaid, render to **static SVG at build time** with `mermaid-cli`; ship the SVG. Do not load the Mermaid runtime on the droplet. Make each node an anchor into the matching reference row.
- **Reference table**, generated not hand-maintained: `config/references.yaml` holds DOIs and roles; `pipeline/references.py` resolves metadata via the Crossref API **at build time on the Studio** and bakes static JSON. Never call Crossref per visitor. Allow `doi` **or** `repo` plus `commit` as the citable identifier, since some tools have only a repo or a Zenodo DOI. Databases show an accession date, not a DOI alone.
- **Provenance**, which is what makes this more than a bibliography: version of every tool used in the run (TNP, ImmuneBuilder, boltz-community, our Tm model's training-set hash), database snapshot dates, run id, git sha.
- **Limitations**, verbatim from `config/targets.yaml`'s `limitations_statement`, plus the upstream-column defects from Section 5.1. This section does more for credibility than any feature: do not bury it, do not soften it.
- Licence (code) and citation line (analysis).

### 8.4 Serving

- nginx serves `.pdb.gz` with `gzip_static`; Mol\* decompresses client side.
- `/healthz` returns 200 and the run id.
- Public, no authentication, no user uploads. No submission form: all compute is offline on the Studio.
- `robots.txt` **disallow plus a `noindex` meta tag** until Phase 4 completes, then both are removed in the final deploy. Public but unlisted while the data is provisional.
- Works at phone width: the explorer collapses to a strip, the drawer stacks under the table, wide tables scroll inside a box.

---

## 9. Positive controls (asserted, never hand-typed)

Resolve each from SAbDab-nano or the PDB at build time; if a control cannot be resolved, log it and drop it rather than inventing a value.

1. **A known anti-EGFR nanobody** with an experimental complex must land in the top quartile of the EGFR affinity axis and must bin to its published epitope.
2. **A known anti-TNFa nanobody** with an experimental complex must do the same for TNFa.
3. **Tm sanity:** no predicted apparent Tm above 95 degrees C in our own model's output. If one appears, the model is extrapolating and the error bar must reflect it. (The upstream column's 101 degrees C is exactly the failure this guards against.)
4. **Trimer sanity:** for TNFa and TL1A, at least some top candidates must show contacts on more than one subunit. If every interface is confined to a single chain, suspect the trimer was not assembled and check before trusting the run.
5. **Scrambled-sequence negative control:** 50 shuffled CDR sequences per target must score worse on the affinity axis than the real candidates. If they do not, the axis is vacuous: run `tools/vacuous_tests.py` and say so on the QC panel rather than shipping a meaningless column.
6. **Hallmark reconciliation:** our recomputed hallmark call must agree with the upstream `hallmark` column for at least 95 percent of rows; investigate the rest.

Control 5 is the one that matters most. Given what the published benchmarks say about co-folding confidence and nanobody specificity, there is a real chance the affinity axis cannot separate real from scrambled. **If that happens it is a finding, not a failure:** report it prominently on the About and Control tabs, keep the axis clearly marked as non-discriminating, and do not let it drive the Pareto front. That is blocker **B3** only if it would invert a headline ranking.

---

## 10. QC (run at the end of Phase 3 and again in Phase 4)

1. `tools/smoke_test.py`: every tab returns 200; every API route returns valid JSON; the drawer loads for a random 20 candidates per target.
2. Row counts reconcile: 2,500 ingested, per-target 500, unique 2,487, `multi_target` = 10.
3. No candidate has a percentile computed against the wrong target.
4. `tools/vacuous_tests.py` on all four axes.
5. Axis correlation matrix per target, published on the QC panel.
6. Arithmetic re-run independently for the Pareto fronts: a second implementation, not the same function called twice.
7. Mol\* context count stays bounded after 30 consecutive drawer opens.
8. Lighthouse or equivalent at phone width; no horizontal page scroll.
9. `tools/secret_scan.py` clean (Section 0.4).
10. Every reference row either `doi_verified = 1` or visibly marked unverified.

---

## 11. Phase plan

Each phase ends with: `REUSE.md`, `DECISIONS.md`, `BUILD_LOG.md` and `PROGRESS.md` current; a commit; a push; and a deploy.

### Phase 1: Foundation, ingest, Tier 0, app skeleton live

- Confirm hardware (2.1). Write `config/hardware.toml` via `tools/hwprobe.py`.
- `.gitignore` **first**, then `git init`, then `gh repo create bellcheddar/WHITTLE --public --source=. --remote=origin`.
- pixi environment; OpenMM and HMMER from conda-forge; ANARCI or `anarcii` verified on one known VHH.
- Ingest all five CSVs into SQLite with the `upstream_*` mapping and every defect from Section 5.1 handled.
- Tier 0 in full (6.1).
- Flask app, three tabs, target selector, Tabulator table with sequence-level columns only, empty drawer shell, `/healthz`.
- `deploy/*` adapted from AlphaFraud and debris; port probed (4.5); nginx vhost; systemd service; TLS.
- **Deploy. The site is live at the end of Phase 1** with Tier 0 columns and a "provisional, structural axes pending" banner.
- Push.

### Phase 2: Tier 1, the Tm model, four axes, Pareto fronts

- NanoBodyBuilder2 across all 2,500; TNP per target; patches, solubility, nativeness, humanness.
- Nearest-neighbour search against SAbDab-nano; experimental badges.
- Train the Tm model (6.3) with grouped nested CV; write `MODEL_CARD.md`.
- Percentiles within target; Pareto fronts; desirability function; axis correlation matrix.
- Drawer fully populated: structure viewer, neighbour superposition, Tm distribution panel.
- Explorer `parcoords` wired to the table with brushing.
- Column chooser with all five groups.
- **Deploy.** Push.

### Phase 3: Tier 2 co-folding and the control set

- Verify `boltz` is the community fork; smoke-test one GLP1R complex and record wall-clock before queueing anything else.
- Antigen constructs and per-target MSA cache; queue in cheapest-first order; `caffeinate`; resumable.
- Interface metrics, per-subunit for the trimers; epitope binning; EGFR and TNFa epitope comparison against known complexes.
- Control set assembled and run blind; Tab 2 built; scrambled negative control (Section 9 item 5).
- Full QC pass (Section 10).
- **Deploy.** Push.

### Phase 4: About, provenance, polish, public release

- `config/references.yaml` resolved through Crossref; reference table baked.
- Mermaid workflow diagram rendered to static SVG with clickable nodes.
- Provenance table; limitations section; licence and citation.
- Icons and favicon set via the `marcs-vibe-icon` house style; `README.md` via `readme_forge` from `README.source.md`.
- `FINDINGS.md`: what the control set showed, where the axes disagree with the upstream columns, what the scrambled control did to the affinity axis.
- Second QC pass; `tools/secret_scan.py`; remove `robots.txt` disallow and the `noindex` tag.
- **Final deploy.** Push. `PROGRESS.md` at 100 percent, no `BLOCKED.md`.

---

## 12. Deliverables

| Deliverable | Where |
|---|---|
| Live app, three tabs | https://whittle.mdeller.com |
| Public repo | https://github.com/bellcheddar/WHITTLE |
| Ranked candidates, five targets | `data/whittle.sqlite`, exportable as CSV from the UI |
| Structures | `data/structures/<target>/<id>.pdb.gz`, co-folds under `data/runs/<target>/` |
| Thermostability model and card | `models/whittle_tm_v1/`, `MODEL_CARD.md` |
| Reference table, DOI-verified | About tab, `config/references.yaml` |
| Provenance record | About tab, `provenance` table |
| Control-set validation | Control tab, `FINDINGS.md` |
| Build record | `BUILD_LOG.md`, `DECISIONS.md`, `REUSE.md`, `PROGRESS.md` |

---

## Appendix A: what is already in the folder

| File | Status |
|---|---|
| `WHITTLE_BUILD_SPEC.md` | this document, authoritative |
| `DESIGN.md` | layout and visual contract for Section 8 |
| `CLAUDE.md` | short pointer, read first by Claude Code |
| `.claude/settings.json` | permissions for an unattended run |
| `config/targets.yaml` | the five targets, authoritative |
| `config/references.yaml` | seed DOI list for the About tab, extend as tools are added |
| `config/deploy.toml.example` | placeholder; real values are gitignored |
| `.gitignore` | written before `git init` |
| `launch.sh` | tmux launcher |
| `shortlist/*.csv` | the 2,500 input sequences |
| `data/` | existing working data, 168 entries, untouched by this build unless pruning is needed |

## Appendix B: decisions already taken (do not relitigate)

| # | Decision | Date |
|---|---|---|
| 1 | Upstream scores retained for comparison, prefixed `upstream_`, never used for ranking | 2026-11-13 |
| 2 | GLP1R scoped to the N-terminal ECD only (24-139) | 2026-11-13 |
| 3 | TNFa and TL1A co-folded as homotrimers, per-subunit interface reporting | 2026-11-13 |
| 4 | pH-switch His design applies to EGFR and TNFa only; TL1A libraries are correct as generated | 2026-11-13 |
| 5 | Thermostability model built in-house on NbThermo; NanoMelt not licensed | 2026-11-13 |
| 6 | Ranking is within-target percentile only; no pooled league table | 2026-11-13 |
| 7 | All compute on the Studio; the droplet is a read-only display layer | 2026-11-13 |
| 8 | Repo public; no user uploads; no request-time compute | 2026-11-13 |
