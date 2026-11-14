# WHITTLE progress

Rewritten after every phase and at least every 30 minutes during compute
(spec 0.1 item 6). Live at <https://whittle.mdeller.com>.

**Overall: about 75 percent.** Phases 1 and 2 are complete and deployed.
Phase 3's machinery is proven end to end by a smoke co-fold and its queue is
the remaining long pole, about 3.65 days. Phase 4 is largely written and waits
only on the co-folding results.

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
- [ ] Antigen MSA cache built (needs the ColabFold server)
- [ ] Smoke co-fold on GLP1R, timed, before the queue opens
- [ ] The queue itself: about 100 complexes, days
- [ ] Controls run blind through the same pipeline
- [ ] Full QC pass

## Phase 4: about, provenance, polish, public release: about 70 percent

- [x] `config/references.yaml` resolved through Crossref; 25 of 25 DOIs verified
- [x] Three seed DOIs found resolving to the wrong paper and corrected
- [x] Provenance table: 32 entries, tool versions, snapshot dates, run id, git sha
- [x] Mermaid workflow rendered to a static SVG with nine reference anchors
- [x] Icon and favicon set in the house style
- [x] `README.md`
- [x] `FINDINGS.md`, generated from the database
- [ ] Final QC pass and secret scan
- [ ] Remove the robots disallow and the noindex tag
- [ ] Final deploy and push

---

## In flight (restarted 2026-11-13 20:40 on MPS)

- **Tier 2 co-folding queue**: 250 complexes, **50 per target**, cheapest
  target first, full quality (recycling 3, 5 samples, 200 steps), on **MPS**.
  The shortlist takes the best 24 of the 50; Tier 2 is the last tier, so
  co-folding exactly 24 would mean it selected nothing (DECISIONS 35).
  `tools/overnight.sh` runs the night as one unit: scrambled controls for the
  two cheapest targets first, then the candidates under
  `tools/supervise_tier2.sh`, then the remaining controls and the affinity
  calibration. Progress: `data/overnight.log`, `data/tier2_supervisor.log`
  and `data/manifests/cofold.jsonl`.
- **Affinity weights are provisional.** The axis blends seven interface
  metrics weighted by argument. `tools/calibrate_affinity.py` proposes weights
  from how well each separates a real candidate from a co-folded scrambled
  one. Applying them is a deliberate edit, so the reasoning lands in a commit.

## Control-set verdicts as of 2026-11-13 20:54

| Axis | Reachable by the CDR shuffle | Verdict |
|---|---|---|
| Affinity | 100 per cent | pending, tonight's co-folded controls decide |
| Specificity | 89 per cent | vacuous on EGFR (AUC 0.600) and TNFa (0.535), informative on the other three |
| Thermo | 60 per cent | vacuous on EGFR (0.590) and TL1A (0.589), informative on the other three |
| Expression | 33 per cent | **untestable by this control**, not vacuous |

- The scrambled control preserves amino-acid composition exactly, so a
  composition-only metric cannot move under it. Two thirds of the expression
  axis is composition-only, so its AUC sits at chance by construction. An
  axis the control cannot reach is now reported untestable rather than
  non-discriminating.
- **Open**: the expression axis needs a composition-varying negative control
  before anything can be claimed about it either way. Not built.
- EGFR is vacuous on both axes the control can test, which is the target
  where it matters most.
- No AUC anywhere exceeds 0.68. None of the sequence-level axes is a strong
  discriminator.

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
