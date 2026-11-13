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

### 3.1 Explorer panel (full width, above the table)

Card, 16px padding. Header row: `Four-axis explorer` (13px, 600) · `drag an axis to brush · N candidates on the Pareto front` (11px, ink-3) · legend · `Reset brush`.

Plotly `parcoords`, ~196px tall, five axes: Affinity · Specificity · Thermostability · Expression · Nativeness.

- Dominated candidates: `#C3CBD9`, 1px, opacity 0.5.
- Pareto front: accent, 2px, opacity 0.92, drawn **over** the dominated set.
- Legend distinguishes the two by **colour and line weight both**.
- An active brush draws as an accent rectangle on its axis with the threshold labelled beside it.
- Brushing filters the table beneath. The table caption always states the active filters in words, so the view is self-describing.

Collapses to a single 44px strip (title · active filters in words · `Expand`) when the column chooser is open or at narrow widths. With a popover open the table matters more than the chart.

### 3.2 Toolbar

Search (`Search ID, target, CDR3 motif`) · `Pareto front · N` as a pressed toggle · `Columns · 6 of 41` · `Export`. All controls >=38px high.

### 3.3 Table

Tabulator, virtual DOM. Default six columns: **Candidate** (id plus target beneath), **Rank**, **Axis profile**, **Tm**, **Nearest PDB**, **Liabilities**.

- Axis profile is a four-bar sparkline (7px wide, up to 24px tall, 3px gap) so the shape of a candidate reads without opening anything.
- Liabilities render as a chip: `Clean` in good, any flag in poor.
- Selected row: accent-wash with a 3px accent inset left border.
- Id column frozen. Wide tables scroll inside an `overflow-x: auto` box.
- Caption, 11px ink-3, states the active filters and marks sample or provisional data explicitly.

### 3.4 Column chooser

Popover anchored to its button, 460px, `0 14px 34px rgba(17,24,39,0.16)`, max-height 470px with internal scroll.

Header: `Columns` · find-a-column input · `Reset`.
Body: the 41 columns in five labelled groups with per-group counts, in this order: **Identity** (6) · **Binding** (7) · **Specificity and liabilities** (11) · **Thermostability** (5) · **Expression and evidence** (8). Each row is a real `<label>` wrapping a real checkbox, >=36px tall, with the unit in ink-3 on the right.
Footer: `N of 41 shown` · `Save as view` · `Done`.

Column state persists in the URL.

### 3.5 Drawer

**Docked beside the table, not overlaid.** `flex: 1 1 400px`; the table takes `flex: 999 1 560px; min-width: 0`, so at phone width the drawer stacks underneath.

Header: candidate id (17px, 700) · `ON FRONT` chip in good when on the Pareto front · 44px close button.

Four sub-tabs: **Structure · Metrics · Sequence · Neighbours**.

Structure panel, in order:
1. **Mol\* viewer**, 236px, on `#0D1117`. Caption bottom-left `model vs <PDB> · TM <score>`, bottom-right `CDR1 / 2 / 3`.
2. **Four axis tiles**, 2x2 grid: uppercase 10px label, 20px value, unit in ink-3, 4px progress rule beneath.
3. **Tm against the full set**: histogram of all candidates for this target with this candidate marked by a 2.5px accent rule, labelled `this one`, **and its error bar drawn**. A bare number hides the uncertainty that matters most.
4. Footnote, 11px ink-3, above a top border: the error-bar and weak-affinity caveat.

The **Neighbours** sub-tab leads with the nearest-neighbour superposition (model over nearest experimental structure, TM-score and per-region identity beside it). This is the panel that justifies having a 3D viewer rather than decorating with one.

Deep link `#/t/<target>/nb/<id>` opens the drawer directly.

### 3.6 Mol\* rules (these are bugs if ignored)

- **One plugin instance.** Reload structures into it, or `plugin.dispose()` on drawer close. Browsers cap WebGL contexts near 16; without this the viewer dies silently after a dozen drawer opens and reads as a random failure.
- Dynamic-import Mol\* (~4 MB) on first drawer open, never in the initial bundle.
- Plotly: custom partial bundle, `scatter` + `bar` + `parcoords` only. `parcoords` is **not** in `plotly.js-basic-dist-min`.
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
