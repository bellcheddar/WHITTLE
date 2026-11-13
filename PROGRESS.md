# WHITTLE progress

Rewritten after every phase and at least every 30 minutes during compute
(spec 0.1 item 6). Live at <https://whittle.mdeller.com>.

**Overall: 25 percent** (Phase 1 of 4 complete)

---

## Phase 1: foundation, ingest, Tier 0, app skeleton live — 100 percent

- [x] Hardware confirmed and `config/hardware.toml` written (M2 Ultra Studio, 128 GB, 646 GB free)
- [x] `.gitignore` before `git init`; public repo `bellcheddar/WHITTLE` created and pushed
- [x] `tools/secret_scan.py` passing on tree and history
- [x] pixi workspace solved: `default` (python 3.12, Tier 0/1) and `cofold` (python 3.11, Tier 2)
- [x] OpenMM 8.6, HMMER, MMseqs2 from conda-forge and bioconda; torch 2.14 with MPS
- [x] ANARCI verified on a known VHH; `anarcii` and a regex route behind it as fallbacks
- [x] All five CSVs ingested: 2,500 rows, 2,487 unique, 10 multi-target sequences
- [x] Every spec 5.1 defect handled at ingest and reconciled
- [x] Tier 0 complete across 2,500 rows: numbering, CDRs, hallmark, liabilities, charge, pI, disulfides, His, CDR3 clusters
- [x] Hallmark reconciliation investigated and reported (spec 9 item 6)
- [x] Within-target percentiles, Pareto fronts, desirability, axis correlation matrix
- [x] Flask app: three tabs, target selector, parcoords explorer, Tabulator table, column chooser, docked drawer, deep links, `/healthz`
- [x] Partial Plotly bundle (scatter, bar, parcoords) at 1.2 MB against the 4.6 MB full dist
- [x] Mol* vendored and dynamic-imported, one instance, disposed on drawer close
- [x] `deploy/*` adapted from AlphaFraud and debris; port 8021 probed free
- [x] Droplet provisioned: service user, venv, systemd unit, nginx vhost, TLS
- [x] **Deployed and live** at <https://whittle.mdeller.com> with a provisional banner
- [x] Smoke test 34/34 against the live host
- [x] Committed and pushed

## Phase 2: Tier 1, the Tm model, four axes, Pareto fronts — 0 percent

- [ ] NanoBodyBuilder2 across all 2,500, OpenMM refinement on the CPU platform
- [ ] TNP per target, one call per FASTA
- [ ] Surface patches, SWI, Protein-Sol, Aggrescan3D
- [ ] Nativeness (AbLang2 or nanoBERT) and BioPhi OASis humanness
- [ ] Nearest neighbour against SAbDab-nano; experimental badges
- [ ] Tm model trained on NbThermo with grouped nested CV; `MODEL_CARD.md`
- [ ] TemStaPro cross-check
- [ ] Drawer fully populated: viewer, superposition, Tm distribution with error bars
- [ ] Deploy, push

## Phase 3: Tier 2 co-folding and the control set — 0 percent

- [ ] boltz-community verified as the community fork; one GLP1R smoke co-fold timed
- [ ] Antigen constructs and per-target MSA cache
- [ ] Queue in cheapest-first order, resumable, under `caffeinate`
- [ ] Interface metrics per subunit for the trimers; epitope binning
- [ ] Control set assembled and run blind; scrambled negative control
- [ ] Full QC pass
- [ ] Deploy, push

## Phase 4: about, provenance, polish, public release — 0 percent

- [ ] `config/references.yaml` resolved through Crossref; reference table baked
- [ ] Mermaid workflow diagram rendered to static SVG
- [ ] Provenance table, limitations, licence and citation
- [ ] Icons and favicon set; `README.md` from `README.source.md`
- [ ] `FINDINGS.md`
- [ ] Second QC pass; secret scan; robots and noindex removed
- [ ] Final deploy, push

---

## Open items carried forward

- The affinity axis is null until Tier 2. The app says so in the banner, the
  caption and the drawer footnote.
- The experimental badge, nearest-PDB and all TNP columns are null until Tier 1.
- Expression and specificity correlate at r = 0.80 on TL1A at Tier 0, because
  both currently draw on liability counts. Expected to decorrelate when Tier 1
  gives them distinct drivers; re-checked every time `pipeline/score.py` runs.
