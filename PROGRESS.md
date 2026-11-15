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

## Phase 3: Tier 2 co-folding and the control set: complete

- [x] boltz-community 2.13 installed in its own environment, MPS available
- [x] The macOS libomp conflict resolved with the fork's own fix
- [x] Antigen constructs fetched and sliced; all five match `config/targets.yaml`
- [x] `pipeline/boltz_input.py`: trimers as one entry with three chain ids, cached antigen MSA, no nanobody MSA
- [x] `pipeline/cofold.py`: cheapest-first queue, resumable, OOM ladder that never falls back to a monomer
- [x] `pipeline/interface.py`: BSA per subunit, interface residues, PRODIGY, epitope binning
- [x] `pipeline/controls.py`: 1000 scrambled (200 per target), 94 known, 23 multi-target; every spec 9 assertion evaluated
- [x] 5 anti-EGFR and 13 anti-TNFa nanobodies resolved from SAbDab-nano, not hand-typed
- [x] `tools/qc_pareto.py`: an independent Pareto implementation, agrees on all five targets
- [x] `tools/vacuous_tests.py`: five ways an axis can be uninformative, checked per target
- [x] Antigen MSA cache built for all five targets
- [x] Smoke co-fold on GLP1R, timed, before the queue opened
- [x] The queue itself: **250 complexes, 50 per target, 11.7 hours, 0 failed**
- [x] Controls run blind through the same pipeline
- [x] **30 scrambled controls co-folded**, 2.0 hours, 0 failed
- [x] Paired tests added for the sequence axes; the thermo verdict changed
- [x] The 40 paired control parents co-folded, 0 failed
- [x] Interface, interaction and complementarity metrics over all 340 complexes
- [x] Per-chain pLDDT split out of the B-factor column; the gap localised to the nanobody chain (DECISIONS 49)
- [x] Affinity weights decided from the paired calibration: five components, the two inverted ones dropped (DECISIONS 51)
- [x] Tier 2 annotates the shortlist rather than reordering it (DECISIONS 50)
- [x] Control set scaled to 200 shuffles per target, additively, leaving the original 250 rows untouched (DECISIONS 55)
- [x] The specificity verdict re-measured after the artefact was found and the old one withdrawn (DECISIONS 59)
- [x] Full QC pass: 44 tests, JS parse, template render, secret scan, style

## Phase 4: about, provenance, polish, public release: complete

- [x] `config/references.yaml` resolved through Crossref; 29 of 29 DOIs verified
- [x] Three seed DOIs found resolving to the wrong paper and corrected
- [x] Provenance table: 32 entries, tool versions, snapshot dates, run id, git sha
- [x] Mermaid workflow rendered to a static SVG with nine reference anchors
- [x] Icon and favicon set in the house style
- [x] `README.md`
- [x] `FINDINGS.md`, generated from the database
- [x] About tab rewritten in lay terms with a glossary for all 50 parameters
- [x] Workflow diagram re-laid out tall and made legible
- [x] Hero panel rebuilt as the shortlist ladder and the whittling gates
- [x] Drawer slides over the page; white viewer; CDRs coloured in all three views
- [x] Final QC pass and secret scan
- [x] `tools/js_check.py` and `tools/template_check.py` gate the deploy
- [x] `tools/smoke_test.py` 34 of 34, checked against the registry not a literal
- [x] Independent Pareto check agrees on all five, asking the same question
- [x] Remove the robots disallow and the noindex tag
- [x] Final deploy and push

---

## In flight as of 2026-11-14 18:35

- **Two EGFR shortlist members are being co-folded.** Because Tier 2 now
  annotates rather than reorders (DECISIONS 50), the ordering dropped the
  affinity axis, and two candidates moved into EGFR's top 24 after the co-fold
  set had already been chosen. 22 of EGFR's 24 carry the Tier 2 annotation and
  the other four targets carry all 24. The shortlist itself is complete and
  fully scored on the three ordering axes for all five targets.
- Everything else in the pipeline has run. The affinity axis is populated for
  the 290 co-folded candidates and displayed; it does not order the table.

## Control-set verdicts as of 2026-11-14 18:25

Each axis is tested twice. **Unpaired** compares every scrambled control
against all 500 candidates, which is the number spec Section 9 asks for.
**Paired** compares each control with the one candidate it was shuffled from,
holding composition, length and framework fixed so CDR order is the only
difference left. Between-sequence variance drops out of the paired form and it
is far the more powerful of the two, so it decides the verdict. The verdict is
decided once per axis, pooled over all five targets, with the per-target counts
shown beside it (DECISIONS 53).

| Axis | Paired, parent beats its own shuffle | Verdict |
|---|---|---|
| Affinity | pending | no scrambled control is scored on this axis yet |
| Specificity | **679 of 839, 81 per cent, p = 9e-77** | **informative, and every target clears it alone** |
| Thermo | **837 of 1000, 84 per cent, p = 1e-109** | **informative, and every target clears it alone** |
| Expression | 982 of 1000, p = 3e-263 | untestable by this control |

- **The 60 per cent specificity result is withdrawn (DECISIONS 59).** It was an
  artefact, not a finding. `pipeline.controls_run` first ran before the TNP
  manifest held anything, so `score_external` found no control carrying `psh`
  (2.0 of the axis's 4.5 weight) or `tnp_flags` (1.0) and dropped both, then
  rebuilt the candidate reference on that same reduced basis. The control
  percentiles were internally consistent, so no log looked wrong; but the
  paired test compares them against `scores.axis_specificity`, built from all
  4.5 of the weight. A third of the axis was being ranked against all of it.
  Feeding the two composition components alone back through `score_external`
  reproduces the old per-target table pair for pair, which is how it was
  confirmed. `score_external` now refuses to place an external row unless it
  used exactly the components the candidates were scored on.
- **The thermo axis is not vacuous.** The unpaired test called it vacuous on
  EGFR (0.591) and TL1A (0.589) and the site said so. Paired, a candidate
  beats its own shuffled twin on 84 per cent of pairs and every target clears
  it on its own. The unpaired AUC was low because comparing one shuffle against
  500 different sequences is swamped by the variance between them, not because
  the axis carries nothing. It read 84 per cent before and after the control
  set was enlarged, which is what first showed the specificity change was in
  the measurement rather than in the biology.
- The scrambled control preserves amino-acid composition exactly, so a
  composition-only metric cannot move under it. Two thirds of the expression
  axis is composition-only, so its AUC sits at chance by construction. An
  axis the control cannot reach is reported untestable rather than
  non-discriminating.
- **Closed, 2026-11-14**: the composition-varying negative control will not
  be built, and the expression axis is reported untested rather than carrying
  a promise of validation that is not coming. A synthetic control cannot
  settle it. SWI and Protein-Sol carry 3.0 of the axis's 4.5 weight and both
  are defined from amino-acid composition, so a deliberately insoluble variant
  scores badly by arithmetic rather than by evidence, and the test would
  confirm only that a function responds to its own inputs. The anchor that
  would settle it does not exist here: `observed_expression` is null for all
  367 control rows, against 77 rows that do carry an observed Tm. Deciding
  this axis needs laboratory expression or solubility on real nanobodies
  (DECISIONS 47).
- TNP silently returns nothing for about 13 per cent of shuffled sequences.
  That is now logged (DECISIONS 58). It costs sample size and does not bias the
  paired comparison, because a dropped shuffle takes its whole pair out and
  shares its parent's exact length.

## Open items carried forward

- **`pipeline.interactions` does not resume.** `cofold`, `interface` and
  `controls_run` all skip what their manifests say is done; `interactions` has
  no manifest read, so adding two late candidates to EGFR cost a 22-minute
  recompute of the 345 samples it already had, for about 40 seconds of new
  work. Harmless since the write became an upsert, but the next late candidate
  pays the same 22 minutes. A gap against `CLAUDE.md`'s "every stage is
  resumable".
- **Spec QC item 8 has not been re-measured since the drawer changed.** The
  recorded pass is from before two columns were added to the drawer's per-seed
  table, and the Chrome extension was disconnected when the re-run was
  attempted. The static argument is that `.kv` is `width: 100%` with
  `overflow-wrap: anywhere` on its cells, so the table's min-content width is
  one character and it cannot exceed its container. That is reasoning, not a
  measurement, and the item deserves the measurement.

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
