# WHITTLE: Design Specification

Authoritative for layout and visual detail. Functional contract is `WHITTLE_BUILD_SPEC.md` Section 8.
Approved concept: docked drawer with a full-width explorer panel above the table.

---

## 1. Shell

Light product UI. Three tabs in a segmented control, top right: **Candidates · Control set · About**.

Masthead, one row, white on `#DDE2EB` bottom border:
`WHITTLE` wordmark (20px, 700, letter-spacing 0.1em) · the backronym expanded at 11px in `#5B6472`, max 430px · tabs · run metadata (`Run <date> · 2,500 VHH · 5 targets`) behind a left border.

**Target selector sits directly under the masthead and is the primary control.** Five targets, one at a time, with the selected target's construct shown beside it (`ECD 25-645`, `soluble trimer 77-233 x3`). GLP1R additionally carries a persistent `ECD-scoped` chip. There is no pooled view and no "all targets" option.

## 2. Palette

| Token | Value | Use |
|---|---|---|
| ground | `#F1F3F7` | page |
| surface | `#FFFFFF` | cards, table, drawer |
| border | `#DDE2EB` | card edges |
| rule | `#EBEFF5` | table rows, internal dividers |
| ink | `#111827` | primary text |
| ink-2 | `#434C5A` | secondary text |
| ink-3 | `#5B6472` | labels, captions (passes 4.5:1 on white and on ground) |
| accent | `#3457A6` | Pareto front, selection, primary action |
| accent-wash | `#F5F8FF` | selected row |
| good | `#1E4D35` on `#DCF0E4` | clean / high percentile |
| warn | `#B5781C` | mid percentile |
| poor | `#8A3A28` on `#FBE5E0` | flag / low percentile |
| viewer | `#0D1117` | Mol* canvas ground |

Inter throughout; tabular numerals on every numeric cell (`font-variant-numeric: tabular-nums`).

Percentile sparkline bars: accent at >=75, warn at 45-74, poor at <45. **Lightness differs as well as hue** so the scale survives greyscale and colour-vision deficiency. Never encode a distinction in hue alone.

## 3. Candidates tab

### 3.1 Shortlist panel (full width, above the table)

Card, 16px padding. Header row: `Shortlist` (13px, 600) · a hint stating what the
pane shows (11px, ink-3) · `Ladder` / `Trade-off` segmented control · `Reset filters`.

Two panes on a `minmax(0,1fr) 338px` grid: the shortlist on the left, the gates
that produced it on the right.

**Left, the ladder (default).** One row per candidate, the highest-ranked first,
24 or 32 or 50 of them. Columns: position · candidate id · a 15px track for each
of the five axes · predicted Tm with its error · three flag pips (liability in a
CDR, camelid hallmark, experimental neighbour). Each track carries the
percentile as a number and a tick at the population median, so a cell reads
against all 500 even though only 32 rows are drawn. An axis with no data is a
hatched track, never an empty one. A rule across the ladder marks where the
shortlist stops. Clicking or pressing Enter on a row opens the drawer. Every
column header sorts.

**Left, the trade-off plane (toggle).** Two chosen axes, one point per
candidate, hollow grey for dominated and solid accent for the Pareto front, with
marginal histograms in the top and right gutters. The frontier of the chosen
*pair* is drawn as a staircase and labelled as the smaller claim it is: it is
not the four-axis Pareto front. Dragging a box sets a brush on both axes, which
is the same brush the table filters on.

**Right, the gates.** Seven filters applied in order: no liability motif in a
CDR · camelid hallmark intact · specificity floor · thermostability floor ·
expression floor · one per CDR3 cluster · affinity floor where co-folded. Each
shows a bar of survivors over what the gate removed, and its own count. A
candidate with no co-fold is not failed by the affinity gate: Tier 2 covers the
top 50 per target, and failing the other 450 on a missing value reads an absence
as evidence. Beneath them, the surviving count against the shortlist size, in
`good` when they match, `warn` when there are too many and `poor` when the gates
are tighter than the library can satisfy. `Solve for 24` raises the three
percentile floors together and keeps the tightest setting that still leaves a
full shortlist standing.

The gates filter the table. The caption always states them in words, as it does
every other filter, so the view is self-describing however it is read.

Collapses to a single 44px strip (title · active filters in words · `Expand`)
when the column chooser is open. At phone width the ladder and the plane both go
(neither survives a 400px column) and the gates stay, because a filter the
reader cannot reach is one that silently removes rows.

**Why this replaced the parcoords panel.** The original specification here was a
Plotly `parcoords` plot: 500 polylines across five axes, per target. It drew the
wrong population, since about 24 of those 500 are ever ordered; its line
crossings carried no information, because the axes are correlated by
construction and the crossings are the same on every target; it could not be
sorted, which is the first thing anyone wants from a ranking; and an axis with
no data was an invisible gap rather than a stated absence. The replacement
answers the two questions the site exists to answer, which are *which ones do I
take* and *what removed the rest*. See DECISIONS.md.

### 3.2 Toolbar

Search (`Search ID, target, CDR3 motif`) · `Pareto front · N` as a pressed toggle · `Columns · 6 of 41` · `Export`. All controls >=38px high.

### 3.3 Table

Tabulator, virtual DOM. Default six columns: **Candidate** (id plus target beneath), **Rank**, **Axis profile**, **Tm**, **Nearest PDB**, **Liabilities**.

- Axis profile is a four-bar sparkline (7px wide, up to 24px tall, 3px gap) so the shape of a candidate reads without opening anything.
- **Numeric cells are shaded by their percentile within the target**, on one accent lightness ramp, with the number still printed. A column with no better direction (a flag, such as pI or net charge) is left unshaded rather than implying one. A number alone makes the reader do the ranking; the shading is what lets a column be scanned.
- Liabilities render as a chip: `Clean` in good, any flag in poor.
- Selected row: accent-wash with a 3px accent inset left border.
- Id column frozen. Wide tables scroll inside an `overflow-x: auto` box.
- Caption, 11px ink-3, states the active filters and marks sample or provisional data explicitly.

### 3.4 Column chooser

Popover anchored to its button, 460px, `0 14px 34px rgba(17,24,39,0.16)`, max-height 470px with internal scroll.

Header: `Columns` · find-a-column input · `Reset`.
Body: every column in five labelled groups, in this order: **Identity** · **Binding** · **Specificity and liabilities** · **Thermostability** · **Expression and evidence**. Each row is a real `<label>` wrapping a real checkbox, >=36px tall, with the unit in ink-3 on the right. The per-group counts are rendered from `pipeline/columns.py` rather than written here: the registry is the one source, and the counts in an earlier draft of this section summed to 37 against a stated total of 41, so they cannot both be authoritative. The registry currently holds 48 columns, each with a glossary entry on the About tab.
Footer: `N of <total> shown` · `Save as view` · `Done`.

Column state persists in the URL.

### 3.5 Drawer

**Docked beside the table, not overlaid.** `flex: 1 1 400px`; the table takes `flex: 999 1 560px; min-width: 0`, so at phone width the drawer stacks underneath.

Header: candidate id (17px, 700) · `ON FRONT` chip in good when on the Pareto front · 44px close button.

Four sub-tabs: **Structure · Metrics · Sequence · Neighbours**.

Structure panel, in order:
1. **Mol\* viewer**, 236px, on `#0D1117`. Caption bottom-left `model vs <PDB> · TM <score>`, bottom-right `CDR1 / 2 / 3`.
2. **Four axis tiles**, 2x2 grid: uppercase 10px label, 20px value, unit in ink-3, 4px progress rule beneath.
3. **Axis strips**, one 38px row per axis: every candidate for this target as a jittered dot, with this candidate marked by an accent needle. This is where the parcoords idea belongs, with the lines taken out. In the drawer the subject is already one candidate, which is the question "where does this one sit" actually asks. An axis with no data is a labelled empty strip, never a missing one.
4. **Tm against the full set**: histogram of all candidates for this target with this candidate marked by a 2.5px accent rule, labelled `this one`, **and its error bar drawn**. A bare number hides the uncertainty that matters most.
5. Footnote, 11px ink-3, above a top border: the error-bar and weak-affinity caveat.

The **Neighbours** sub-tab leads with the nearest-neighbour superposition (model over nearest experimental structure, TM-score and per-region identity beside it). This is the panel that justifies having a 3D viewer rather than decorating with one.

Deep link `#/t/<target>/nb/<id>` opens the drawer directly.

### 3.6 Mol\* rules (these are bugs if ignored)

- **One plugin instance.** Reload structures into it, or `plugin.dispose()` on drawer close. Browsers cap WebGL contexts near 16; without this the viewer dies silently after a dozen drawer opens and reads as a random failure.
- Dynamic-import Mol\* (~4 MB) on first drawer open, never in the initial bundle.
- Plotly: custom partial bundle, `scatter` + `bar` + `parcoords` only. `parcoords` is **not** in `plotly.js-basic-dist-min`. Since 3.1 stopped using it, the Candidates tab loads Plotly on first drawer open (for the Tm histogram) rather than on page load; the Control set tab still loads it for its scatters.
- Vendored assets cache-busted with `?v=<mtime>`.

## 4. Control set tab

Same table component, different data file. Above it, per axis, an observed-against-predicted scatter with R squared and a residual plot. Identity line dashed in ink-3; points in accent; any control that fails its assertion in poor.

Two panels that must be visible without scrolling:
- **Coverage caveat.** Only EGFR and TNFa have target-matched references. TL1A, PRLR and GLP1R have no published nanobody benchmark, so this tab validates the general axes everywhere but the target-aware axis only for those two.
- **Scrambled negative control.** Real against shuffled-CDR distributions per axis. If the affinity axis does not separate them, that is reported here as a finding, in plain words, and the axis is marked non-discriminating wherever it appears.

## 5. About tab

Reading width 760px, left-aligned, generous leading.

1. **What WHITTLE does**, four sentences.
2. **Workflow diagram**: static SVG rendered at build time from Mermaid. Tier 0 → 1 → 1b → 2 → ranking, with per-tier counts and the cheapest-first Tier 2 order. Each node is an anchor into its reference row.
3. **Reference table**: tool or database · role · version or snapshot date · citation · DOI badge. Generated from `config/references.yaml`, Crossref-resolved at build time. Unverified rows visibly marked.
4. **Provenance**: tool versions, database snapshot dates, run id, git sha.
5. **Limitations**, set apart with a 3px accent left border and never collapsed: the weak affinity axis, Tm error bars off the training manifold, the positive-only structure badge, GLP1R ECD scoping, within-target-only ranking, and the upstream column defects.
6. Licence and citation line.

## 6. Accessibility and responsive

- Real `<button>`, `<a href>`, `<input>` + `<label>` everywhere, including in static panels. Never `role` or click handlers on a div.
- `aria-label` on every icon-only control; `aria-pressed` on toggles; `aria-current="page"` on the active tab.
- Decorative SVG `aria-hidden="true"`; meaningful SVG `role="img"` with a label that says what the picture shows.
- Text contrast >=4.5:1, including 11px ink-3 labels.
- Touch targets >=44px.
- Phone width: explorer collapses to its strip, drawer stacks under the table, toolbar wraps, 16px side gutter, no horizontal page scroll.
- Prefers-reduced-motion respected: no drawer slide animation.

## 7. Copy rules

- "predicted", "estimated", "ranked". Never "confirmed", "validated", "proven".
- Percentile cells always labelled as percentile-within-target, never as an absolute measurement.
- `upstream_*` columns carry "not comparable between targets" wherever the upstream composite appears.
- `s_his` is labelled a **design criterion**, not a quality metric, and only appears for EGFR and TNFa.
- British English. No em dashes.
