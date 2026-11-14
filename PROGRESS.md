# WHITTLE progress

Rewritten after every phase and at least every 30 minutes during compute
(spec 0.1 item 6). Live at <https://whittle.mdeller.com>.

**Overall: about 93 percent.** Phases 1 and 2 are complete and deployed.
Phase 3's co-folding is done: 250 candidate complexes (50 per target) and 30
scrambled controls, all on MPS, 0 failed. The 40 paired control parents are
finishing now, and the interface, interaction and complementarity metrics plus
the affinity calibration run straight after them. Phase 4 waits on that
calibration and on the final QC pass.

---

## Phase 1: foundation, ingest, Tier 0, app skeleton live: 100 percent

- [x] Hardware confirmed, `config/hardware.toml` written (M2 Ultra Studio, 128 GB)
- [x] `.gitignore` before `git init`; public repo created and pushed
- [x] `tools/secret_scan.py` passing on tree and history before every push
- [x] pixi workspace: `default` (python 3.12) and `cofold` (python 3.11, boltz)
- [x] ANARCI verified on a known VHH; `anarcii` and a regex route behind it
- [x] 2,500 rows ingested, 2,487 unique, 10 multi-target; every spec 5.1 defect handled
- [x] Tier 0 complete: numbering, CDRs, hallmark, liabilities, charge, pI, disulfides, clusters
- [x] Hallmark reconciliation investigated and reported
- [x] Within-target percentiles, Pareto fronts, desirability, correlation matrix
- [x] Flask app: three tabs, target selector, parcoords explorer, Tabulator table, 41-column chooser, docked drawer, deep links
- [x] Droplet provisioned on port 8021, nginx, systemd, TLS
- [x] **Live** with a provisional banner; smoke test 34/34 against the live host

## Phase 2: Tier 1, the Tm model, four axes, Pareto fronts: complete

- [x] NanoBodyBuilder2 across all 2,500: 2,499 models in 21 minutes
- [x] The 56 refusals diagnosed and 55 recovered with a germline N-terminus
- [x] SWI and Protein-Sol reimplemented and validated against the authors' own outputs
- [x] Surface patches (SAP-style neighbourhood sum, after a cluster rule proved degenerate)
- [x] OASis humanness for all 2,500
- [x] Nearest experimental neighbour and TM-align for 2,499
- [x] Tm training set assembled: 650 training, 77 genuinely held out
- [x] Tm model trained with grouped nested CV; `MODEL_CARD.md` written
- [x] Error bars recalibrated empirically after the first pass under-covered
- [x] TNP across all five targets: 2,397 of 2,500 profiled after the padded retry
- [x] AbLang2 nativeness across all 2,500, 92.6 minutes, 0 failed
- [x] `pipeline.tier1` assembling the `structures` table, 2,500 rows
- [x] `pipeline.tm_predict` over all 2,500: 58.8 to 91.1 degrees C, none above 95
- [x] Rescored: 13 to 25 candidates on the Pareto front per target
- [x] Deployed and pushed

## Phase 3: Tier 2 co-folding and the control set: machinery done, compute pending

- [x] boltz-community 2.13 installed in its own environment, MPS available
- [x] The macOS libomp conflict resolved with the fork's own fix
- [x] Antigen constructs fetched and sliced; all five match `config/targets.yaml`
- [x] `pipeline/boltz_input.py`: trimers as one entry with three chain ids, cached antigen MSA, no nanobody MSA
- [x] `pipeline/cofold.py`: cheapest-first queue, resumable, OOM ladder that never falls back to a monomer
- [x] `pipeline/interface.py`: BSA per subunit, interface residues, PRODIGY, epitope binning
- [x] `pipeline/controls.py`: 250 scrambled, 94 known, 10 multi-target; every spec 9 assertion evaluated
- [x] 5 anti-EGFR and 13 anti-TNFa nanobodies resolved from SAbDab-nano, not hand-typed
- [x] `tools/qc_pareto.py`: an independent Pareto implementation, agrees on all five targets
- [x] `tools/vacuous_tests.py`: five ways an axis can be uninformative, checked per target
- [x] Antigen MSA cache built for all five targets
- [x] Smoke co-fold on GLP1R, timed, before the queue opened
- [x] The queue itself: **250 complexes, 50 per target, 11.7 hours, 0 failed**
- [x] Controls run blind through the same pipeline
- [x] **30 scrambled controls co-folded**, 2.0 hours, 0 failed
- [x] Paired tests added for the sequence axes; the thermo verdict changed
- [~] The 40 paired control parents: 35 of 40 done, 0 failed
- [ ] Interface, interaction and complementarity metrics over the controls
- [ ] Affinity weights decided from the paired calibration
- [ ] Full QC pass

## Phase 4: about, provenance, polish, public release: about 85 percent

- [x] `config/references.yaml` resolved through Crossref; 29 of 29 DOIs verified
- [x] Three seed DOIs found resolving to the wrong paper and corrected
- [x] Provenance table: 32 entries, tool versions, snapshot dates, run id, git sha
- [x] Mermaid workflow rendered to a static SVG with nine reference anchors
- [x] Icon and favicon set in the house style
- [x] `README.md`
- [x] `FINDINGS.md`, generated from the database
- [x] About tab rewritten in lay terms with a glossary for all 48 parameters
- [x] Workflow diagram re-laid out tall and made legible
- [x] Hero panel rebuilt as the shortlist ladder and the whittling gates
- [x] Drawer slides over the page; white viewer; CDRs coloured in all three views
- [ ] Final QC pass and secret scan
- [ ] Remove the robots disallow and the noindex tag
- [ ] Final deploy and push

---

## In flight as of 2026-11-14 12:25

- **Paired control parents**: 35 of 40 done, 0 failed, about 25 minutes left.
  Each scrambled control's own parent is co-folded so the affinity comparison
  is paired. Without it the comparison is confounded: the co-folded candidates
  are the top 50 by Tier 1 rank, which is driven by patch surface
  hydrophobicity, so they sit below the population (median PSH 93.6 against
  106.4 across all 500 GLP1R candidates) while the controls were shuffled from
  a random sample whose parents sit at 107.9 (DECISIONS 35).
- **Then**, straight after, `tools/calibrate_run.sh` runs `pipeline.interface`,
  `pipeline.interactions` and `pipeline.complementarity` over about 320
  complexes and then `tools/calibrate_affinity.py`. Roughly 1.5 to 2 hours on
  CPU, measured from GLP1R's 300 samples last night: interface a few minutes,
  interactions 14, complementarity 4.
- **Affinity weights are provisional.** The axis blends seven interface
  metrics weighted by argument. `tools/calibrate_affinity.py` proposes weights
  from how well each separates a real candidate from its co-folded shuffled
  twin. Applying them is a deliberate edit, so the reasoning lands in a commit.
  `polar_fraction` was the strongest discriminator on the unpaired read at AUC
  0.841 and is not currently a component; if it survives the paired test it has
  earned a place in the axis.

## Control-set verdicts as of 2026-11-14 12:23

Each axis is tested twice. **Unpaired** compares every scrambled control
against all 500 candidates, which is the number spec Section 9 asks for.
**Paired** compares each control with the one candidate it was shuffled from,
holding composition, length and framework fixed so CDR order is the only
difference left. Between-sequence variance drops out of the paired form and it
is far the more powerful of the two, so it decides the verdict.

| Axis | Unpaired AUC | Paired, parent beats its own shuffle | Verdict |
|---|---|---|---|
| Affinity | pending | pending | tonight's co-folded controls decide |
| Specificity | 0.535 to 0.679 | 144 of 240, p = 0.0023 | informative on PRLR, weak on GLP1R and TL1A, **non-discriminating on EGFR and TNFa** |
| Thermo | 0.589 to 0.629 | **210 of 250, p = 5e-29** | **informative on all five** |
| Expression | 0.439 to 0.522 | 136 of 250, p = 0.18 | untestable by this control |

- **The thermo axis is not vacuous.** The unpaired test called it vacuous on
  EGFR (0.591) and TL1A (0.589) and the site said so. Paired, a candidate
  beats its own shuffled twin on 84 per cent of pairs and every target clears
  p < 1e-4 on its own (PRLR 47 of 50). The unpaired AUC was low because
  comparing one shuffle against 500 different sequences is swamped by the
  variance between them, not because the axis carries nothing.
- **Specificity is weak and target-dependent.** Real on PRLR (32 of 48,
  p = 0.029), unconfirmed on GLP1R and TL1A, nothing on TNFa (24 of 49). EGFR
  is non-discriminating, which is the target where it matters most.
- The scrambled control preserves amino-acid composition exactly, so a
  composition-only metric cannot move under it. Two thirds of the expression
  axis is composition-only, so its AUC sits at chance by construction. An
  axis the control cannot reach is reported untestable rather than
  non-discriminating.
- **Open**: the expression axis needs a composition-varying negative control
  before anything can be claimed about it either way. Not built, and note that
  a synthetic one risks circularity: SWI and Protein-Sol are defined from
  composition, so a deliberately hydrophobic variant scores worse by
  construction rather than by evidence. There is no observed expression or
  solubility measurement anywhere in the control set to anchor it against.
- Four verdict states are now reported, not two: informative, weak and
  unconfirmed, non-discriminating, and untestable by this control. A hard
  threshold at AUC 0.6 had been splitting indistinguishable evidence into
  opposite labels (DECISIONS 46).

## Open items carried forward

- Mol* fetches the RCSB entry list from the viewer's own UI on every
  initialisation, which fails and logs a console error on the live site. It is
  cosmetic and predates the complex viewer, but it is an outbound call from
  each visitor's browser that the read-only display layer does not need.
- The affinity axis is null until Tier 2 completes. The banner, the caption and
  the drawer footnote all say so.
- One candidate, TNFa-0500, has no structural model: its J region is corrupt and
  no numbering scheme can close the domain. It keeps its sequence axes.
- Expression and specificity correlated at r = 0.80 on TL1A at Tier 0, because
  both drew on liability counts. Re-checked on every scoring run; expected to
  fall once Tier 1 gives them distinct drivers.
- The Tm model orders well and is not calibrated in absolute degrees on
  sequences from a different measurement regime. Stated in `MODEL_CARD.md` and
  enforced by making the axis a percentile.
