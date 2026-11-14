/* Candidates tab: target selector, the shortlist hero, Tabulator table, column
   chooser, docked drawer. Spec 8.1 and DESIGN.md Section 3.

   The hero is three panels over one set of rows: the ladder of candidates in
   contention, the gates that produced it, and a two-axis trade-off plane. All
   three draw in the DOM or on a canvas, so the hero no longer loads Plotly;
   Plotly and Mol* both load on first drawer open (DESIGN 3.6). */
'use strict';

(function () {
  const META = window.WHITTLE_META || {};
  const REG = META.columns || { columns: [], groups: [], defaults: [], total: 0 };
  const BY_KEY = Object.fromEntries((REG.columns || []).map((c) => [c.key, c]));
  const TARGETS = META.targets || [];

  const el = (id) => document.getElementById(id);
  let bundle = null;       // the current target's rows
  let table = null;
  let tableReady = false;
  let plotlyReady = null;
  let axesMeta = null;
  let tmDist = null;
  let currentDetail = null;

  /* ------------------------------------------------------------ formatters */

  function liabilityChip(cell) {
    const row = cell.getRow().getData();
    const n = row.liabilities;
    if (n === null || n === undefined) return '';
    if (n === 0) return '<span class="chip good">Clean</span>';
    const inCdr = row.n_cdr_liabilities || 0;
    const label = inCdr ? `${n} &middot; ${inCdr} in CDR` : `${n}`;
    return `<span class="chip poor" title="${n} liability motifs, ${inCdr} inside a CDR">${label}</span>`;
  }

  function formatterFor(col) {
    if (col.key === 'id') {
      return (cell) => {
        const row = cell.getRow().getData();
        return `<span class="cid">${W.escape(row.id)}</span>`
          + `<span class="csub">${W.escape(row.target)}`
          + `${row.on_front ? ' &middot; on front' : ''}</span>`;
      };
    }
    if (col.key === 'axis_profile') return (cell) => W.sparkline(cell.getRow().getData());
    if (col.key === 'liabilities') return liabilityChip;
    if (col.key === 'multi_target') {
      return (cell) => (cell.getValue()
        ? '<span class="chip warn" title="This sequence also appears under another target and is used as an internal promiscuity control">Multi</span>'
        : '');
    }
    if (col.key === 'experimental_badge') {
      return (cell) => (cell.getValue()
        ? '<span class="chip good" title="A near-identical sequence has a PDB entry. A positive signal only.">PDB</span>'
        : '<span style="color:var(--ink-3)">&ndash;</span>');
    }
    if (col.kind === 'pct') return (cell) => W.pctCell(cell.getValue());
    if (col.kind === 'num') {
      // A number alone makes the reader do the ranking. The cell is shaded by
      // the value's percentile within this target, on one lightness ramp so it
      // survives greyscale and colour-vision deficiency, and a column with no
      // better direction (a flag, like pI) is left unshaded rather than
      // implying one.
      const directional = col.higher_is_better !== null
        && col.higher_is_better !== undefined;
      return (cell) => {
        const value = cell.getValue();
        if (value === null || value === undefined) {
          return '<span style="color:var(--ink-3)">pending</span>';
        }
        const text = W.fmt(value, col.decimals ?? 2);
        if (!directional) return text;
        const pct = pctOf(col.key)[cell.getRow().getData().id];
        if (pct === null || pct === undefined) return text;
        const shade = (0.05 + 0.5 * (pct / 100) * (pct / 100)).toFixed(3);
        const strong = pct >= 72 ? ';color:var(--ink);font-weight:600' : '';
        return `<span class="shaded" style="background:rgba(52,87,166,${shade})`
          + `${strong}" title="${Math.round(pct)}th percentile within this `
          + `target">${text}</span>`;
      };
    }
    return (cell) => {
      const value = cell.getValue();
      if (value === null || value === undefined || value === '') {
        return '<span style="color:var(--ink-3)">&ndash;</span>';
      }
      return W.escape(value);
    };
  }

  function columnDef(key) {
    const col = BY_KEY[key];
    if (!col) return null;
    const def = {
      title: col.label,
      field: col.key === 'axis_profile' ? 'axis_affinity' : col.key,
      headerTooltip: col.help + (col.unit ? ` (${col.unit})` : ''),
      formatter: formatterFor(col),
      cssClass: col.kind === 'num' || col.kind === 'pct' ? 'num' : '',
      headerSortTristate: true,
      minWidth: col.key === 'id' ? 150 : 86,
    };
    if (col.key === 'id') { def.frozen = true; def.width = 168; }
    if (col.key === 'axis_profile') { def.headerSort = false; def.width = 118; }
    if (col.key === 'liabilities') def.width = 128;
    if (col.kind === 'pct') def.width = 128;
    return def;
  }

  /* ------------------------------------------------------- the target selector */

  function targetInfo(name) {
    return TARGETS.find((t) => t.target === name) || null;
  }

  function renderTargetBar(name) {
    const info = targetInfo(name);
    el('targetset').querySelectorAll('button').forEach((btn) => {
      btn.setAttribute('aria-pressed', String(btn.dataset.target === name));
    });
    if (!info) return;
    el('construct').textContent = info.construct_label
      + ` · antigen ${W.fmt(info.antigen_length)} aa · complex about ${W.fmt(info.complex_length)} aa`;
    const chips = [];
    if (info.ecd_scoped) {
      chips.push('<span class="chip warn" title="Results for this target are scoped to the N-terminal extracellular domain only. An ECD co-fold cannot support any claim about binding the intact receptor in a membrane.">ECD-scoped</span>');
    }
    if (info.stoichiometry > 1) {
      chips.push(`<span class="chip accent" title="Co-folded as a homotrimer. Interface metrics are reported per subunit and in total.">Trimer x${info.stoichiometry}</span>`);
    }
    if (info.his_enriched) {
      chips.push('<span class="chip" title="This library was His-enriched on purpose as a pH-switch design. s_his is a design criterion, not a quality metric.">His-enriched library</span>');
    }
    if (!info.known_nanobodies) {
      chips.push('<span class="chip" title="No known nanobodies for this target in SAbDab-nano, so there is no target-matched reference and no target-aware comparison axis.">No target-matched reference</span>');
    }
    el('targetchips').innerHTML = chips.join(' ');
  }

  /* ----------------------------------------------------------- filtering */

  function activeRows() {
    if (!bundle) return [];
    let rows = applyGates(bundle.rows);
    if (W.state.pareto) rows = rows.filter((r) => r.on_front);
    const query = (W.state.search || '').trim().toUpperCase();
    if (query) {
      rows = rows.filter((r) => (r.id || '').toUpperCase().includes(query)
        || (r.target || '').toUpperCase().includes(query)
        || (r.v_gene || '').toUpperCase().includes(query)
        || (r._cdr3 || '').toUpperCase().includes(query));
    }
    for (const [axis, range] of Object.entries(W.state.brush)) {
      rows = rows.filter((r) => {
        const value = r[axis];
        if (value === null || value === undefined) return false;
        return value >= range[0] && value <= range[1];
      });
    }
    return rows;
  }

  function filterWords() {
    const parts = [];
    GATES.filter(gateActive).forEach((gate) => {
      parts.push(gate.kind === 'toggle'
        ? gate.name.toLowerCase()
        : `${gate.name.toLowerCase()} at the ${gateValue(gate)}th percentile`);
    });
    if (W.state.pareto) parts.push('Pareto front only');
    if (W.state.search) parts.push(`matching "${W.state.search}"`);
    for (const [axis, range] of Object.entries(W.state.brush)) {
      parts.push(`${W.AXIS_LABELS[axis] || axis} between the `
        + `${Math.round(range[0])}th and ${Math.round(range[1])}th percentile`);
    }
    return parts;
  }

  function renderCaption(rows) {
    const info = targetInfo(W.state.target);
    const words = filterWords();
    const pending = (META.axes_pending || []).map((a) => W.AXIS_LABELS['axis_' + a] || a);
    let text = `${W.fmt(rows.length)} of ${W.fmt(bundle ? bundle.n : 0)} `
      + `${W.state.target} candidates`;
    text += words.length ? `, filtered to ${words.join('; ')}.` : ', no filters applied.';
    text += ' Every column is predicted or estimated. Percentile columns are'
      + ' percentiles within this target and are not comparable between targets.';
    if (pending.length) {
      text += ` Provisional: ${pending.join(' and ')} `
        + `${pending.length > 1 ? 'axes have' : 'axis has'} not been computed yet.`;
    }
    if (info && info.ecd_scoped) {
      text += ' GLP1R results are scoped to the N-terminal ECD, not the intact receptor.';
    }
    el('caption').textContent = text;
    el('stripfilters').textContent = words.length ? words.join('; ') : 'no filters';
  }

  function refresh() {
    const rows = activeRows();
    // Tabulator warns and then throws if replaceData lands before `tableBuilt`.
    // The table is constructed with its data, so the first refresh has nothing
    // to replace.
    if (table && tableReady) table.replaceData(rows);
    renderCaption(rows);
    const front = bundle ? bundle.rows.filter((r) => r.on_front).length : 0;
    el('frontcount').textContent = W.fmt(front);
    const cols = W.state.cols || REG.defaults;
    el('colcount').textContent = `${cols.length} of ${REG.total}`;
    el('colfoot').textContent = `${cols.length} of ${REG.total} shown`;
    const params = new URLSearchParams({ target: W.state.target || '',
      columns: cols.join(','), ids: rows.map((r) => r.id).join(',') });
    el('exportbtn').href = `/api/export.csv?${params.toString()}`;
    drawHero(rows);
  }

  /* --------------------------------------------------------------- the table */

  function buildTable() {
    const cols = (W.state.cols || REG.defaults).map(columnDef).filter(Boolean);
    if (table) { table.destroy(); table = null; }
    tableReady = false;
    table = new Tabulator('#table', {
      data: activeRows(),
      columns: cols,
      layout: 'fitDataStretch',
      height: '560px',
      renderVertical: 'virtual',     // 500 x 41 cells is not a naive DOM table
      index: 'id',
      placeholder: 'No candidate matches the current filters.',
      rowFormatter: (row) => {
        if (row.getData().on_front) row.getElement().classList.add('on-front');
      },
    });
    table.on('tableBuilt', () => { tableReady = true; });
    table.on('rowClick', (_e, row) => openDrawer(row.getData().id));
  }

  /* ------------------------------------------------ 3.1 the shortlist hero

     Three panels share one set of rows: the ladder (which 24), the gates (what
     removed the other 476) and the trade-off plane (what one axis costs
     another). None of them draws a line that crosses another line, which is
     what the parcoords panel they replaced did 500 times over. */

  const HERO_AXES = [
    ['axis_affinity', 'Aff', 'Affinity'],
    ['axis_specificity', 'Spec', 'Specificity'],
    ['axis_thermo', 'Thermo', 'Thermostability'],
    ['axis_expression', 'Expr', 'Expression'],
    ['axis_nativeness', 'Native', 'Nativeness'],
  ];
  const SHORTLIST = META.shortlist_per_target || 24;

  // Lower is better for these, so their percentile is flipped before it is
  // shaded or sorted. Taken from the column registry where it is declared.
  function lowerIsBetter(field) {
    const col = BY_KEY[field];
    return col ? col.higher_is_better === false : false;
  }

  let pctCache = {};
  function pctOf(field) {
    if (pctCache[field]) return pctCache[field];
    const map = {};
    const values = [];
    (bundle ? bundle.rows : []).forEach((r) => {
      const v = r[field];
      if (typeof v === 'number' && isFinite(v)) values.push(v);
    });
    if (values.length > 1) {
      values.sort((a, b) => a - b);
      const flip = lowerIsBetter(field);
      (bundle.rows).forEach((r) => {
        const v = r[field];
        if (typeof v !== 'number' || !isFinite(v)) { map[r.id] = null; return; }
        let lo = 0, hi = values.length;
        while (lo < hi) { const mid = (lo + hi) >> 1; if (values[mid] < v) lo = mid + 1; else hi = mid; }
        const pct = (lo / (values.length - 1)) * 100;
        map[r.id] = flip ? 100 - pct : pct;
      });
    } else {
      (bundle ? bundle.rows : []).forEach((r) => { map[r.id] = null; });
    }
    pctCache[field] = map;
    return map;
  }

  function medianOf(field) {
    const values = [];
    (bundle ? bundle.rows : []).forEach((r) => {
      const v = r[field];
      if (typeof v === 'number' && isFinite(v)) values.push(v);
    });
    if (!values.length) return null;
    values.sort((a, b) => a - b);
    return values[Math.floor(values.length / 2)];
  }

  /* ------------------------------------------------------------- the gates */

  const GATES = [
    { id: 'cdrliab', name: 'No liability motif in a CDR',
      sub: 'glycosylation, deamidation, isomerisation', kind: 'toggle',
      test: (r) => !r.n_cdr_liabilities },
    { id: 'hallmark', name: 'Camelid hallmark intact',
      sub: 'all four of IMGT 42, 49, 50, 52', kind: 'toggle',
      test: (r) => r.hallmark_n === 4 },
    { id: 'spec', name: 'Specificity floor', sub: 'percentile within target',
      kind: 'range', test: (r, v) => r.axis_specificity !== null
        && r.axis_specificity !== undefined && r.axis_specificity >= v },
    { id: 'thermo', name: 'Thermostability floor', sub: 'percentile within target',
      kind: 'range', test: (r, v) => r.axis_thermo !== null
        && r.axis_thermo !== undefined && r.axis_thermo >= v },
    { id: 'expr', name: 'Expression floor', sub: 'percentile within target',
      kind: 'range', test: (r, v) => r.axis_expression !== null
        && r.axis_expression !== undefined && r.axis_expression >= v },
    { id: 'cluster', name: 'One per CDR3 cluster',
      sub: 'keeps the best of each 90 per cent group', kind: 'toggle',
      cluster: true },
    // A candidate with no co-fold is not failed by this gate. Tier 2 covers
    // the top 50 per target, so failing the other 450 on a missing value would
    // be reading an absence as evidence.
    { id: 'aff', name: 'Affinity floor, where co-folded',
      sub: 'skipped for a candidate with no co-fold', kind: 'range',
      test: (r, v) => r.axis_affinity === null || r.axis_affinity === undefined
        || r.axis_affinity >= v },
  ];

  function gateValue(gate) {
    const set = W.state.gates || {};
    if (gate.kind === 'toggle') return set[gate.id] === true;
    const v = set[gate.id];
    return typeof v === 'number' ? v : 0;
  }

  function gateActive(gate) {
    return gate.kind === 'toggle' ? gateValue(gate) === true : gateValue(gate) > 0;
  }

  // Applied in order, so each gate sees the survivors of the one above it.
  function applyGates(rows) {
    let pool = rows;
    for (const gate of GATES) {
      if (!gateActive(gate)) continue;
      if (gate.cluster) {
        const seen = new Set();
        const keep = [];
        pool.slice().sort((a, b) => (b.desirability || 0) - (a.desirability || 0))
          .forEach((r) => {
            const key = r.cdr3_cluster;
            if (key === null || key === undefined) { keep.push(r); return; }
            if (!seen.has(key)) { seen.add(key); keep.push(r); }
          });
        const kept = new Set(keep.map((r) => r.id));
        pool = pool.filter((r) => kept.has(r.id));
      } else {
        pool = pool.filter((r) => gate.test(r, gateValue(gate)));
      }
    }
    return pool;
  }

  function drawGates() {
    if (!bundle) return;
    const all = bundle.rows;
    let pool = all;
    // Order matters: the grid places gname and gcount on the first row and
    // lets the bar and any slider span both columns beneath them.
    let html = `<div class="gate entering"><span class="gname">Entering</span>`
      + `<span class="gcount">${W.fmt(all.length)}</span>`
      + `<span class="gbar"><i style="width:100%"></i></span></div>`;

    for (const gate of GATES) {
      const before = pool.length;
      if (gateActive(gate)) pool = applyGates(all).length >= 0 ? pool : pool;
      // Apply just this gate to the running pool, so the bar is the survivors
      // of everything above it rather than of this gate alone.
      let after = before;
      if (gateActive(gate)) {
        if (gate.cluster) {
          const seen = new Set();
          const keep = [];
          pool.slice().sort((a, b) => (b.desirability || 0) - (a.desirability || 0))
            .forEach((r) => {
              const key = r.cdr3_cluster;
              if (key === null || key === undefined) { keep.push(r); return; }
              if (!seen.has(key)) { seen.add(key); keep.push(r); }
            });
          const kept = new Set(keep.map((r) => r.id));
          pool = pool.filter((r) => kept.has(r.id));
        } else {
          pool = pool.filter((r) => gate.test(r, gateValue(gate)));
        }
        after = pool.length;
      }
      const kept = (after / all.length) * 100;
      const lost = (before / all.length) * 100;
      // A toggle is its own checkbox beside the name; a floor keeps the name
      // as a label and puts the slider on its own row under the bar.
      const name = gate.kind === 'toggle'
        ? `<label class="gname"><input type="checkbox" data-gate="${gate.id}"`
          + `${gateValue(gate) ? ' checked' : ''}>`
          + `<span>${W.escape(gate.name)}`
          + `<small>${W.escape(gate.sub)}</small></span></label>`
        : `<span class="gname"><span>${W.escape(gate.name)}`
          + `<small>${W.escape(gate.sub)}</small></span></span>`;
      const slider = gate.kind === 'range'
        ? `<span class="gctl"><input type="range" min="0" max="95" step="1"`
          + ` value="${gateValue(gate)}" data-gate="${gate.id}"`
          + ` aria-label="${W.escape(gate.name)}, percentile floor">`
          + `<output>${gateValue(gate)}</output></span>`
        : '';
      html += `<div class="gate${gateActive(gate) ? '' : ' off'}">`
        + name
        + `<span class="gcount">${gateActive(gate)
            ? W.fmt(after) : '<span class="goff">off</span>'}</span>`
        + `<span class="gbar"><i class="lost" style="width:${lost}%"></i>`
        + `<i style="width:${kept}%"></i></span>`
        + slider
        + `</div>`;
    }
    el('gates').innerHTML = html;

    el('gates').querySelectorAll('[data-gate]').forEach((input) => {
      const gate = GATES.find((g) => g.id === input.dataset.gate);
      input.addEventListener('input', () => {
        W.state.gates = W.state.gates || {};
        W.state.gates[gate.id] = gate.kind === 'toggle'
          ? input.checked : Number(input.value);
        if (gate.kind === 'toggle' && !input.checked) delete W.state.gates[gate.id];
        if (gate.kind === 'range' && Number(input.value) === 0) {
          delete W.state.gates[gate.id];
        }
        W.writeUrl(true);
        refresh();
      });
    });

    const survivors = pool.length;
    el('gatecount').textContent = W.fmt(survivors);
    const state = el('gatestate');
    const hint = el('gatehint');
    const anyGate = GATES.some(gateActive);
    if (!anyGate) {
      state.className = 'state';
      state.textContent = 'no gates applied';
      hint.textContent = `All ${W.fmt(all.length)} candidates are in the table. `
        + `The ladder shows the top ${SHORTLIST} by desirability.`;
    } else if (survivors === SHORTLIST) {
      state.className = 'state on';
      state.textContent = `exactly ${SHORTLIST}`;
      hint.textContent = 'These are the candidates that would be ordered.';
    } else if (survivors > SHORTLIST) {
      state.className = 'state over';
      state.textContent = `${W.fmt(survivors - SHORTLIST)} more than the shortlist`;
      hint.textContent = `Raise a floor, or take the top ${SHORTLIST} of these `
        + 'by desirability.';
    } else {
      state.className = 'state under';
      state.textContent = `${W.fmt(SHORTLIST - survivors)} short of ${SHORTLIST}`;
      hint.textContent = 'The gates are tighter than the library can satisfy. '
        + 'Which one gives is a scientific decision rather than a slider.';
    }
  }

  function solveForShortlist() {
    // Raise the three percentile floors together and keep the tightest setting
    // that still leaves at least a full shortlist standing. Stopping at the
    // first setting to fall BELOW the shortlist was worse: it overshot, and a
    // button called "Solve for 24" that returns 18 has not solved for 24.
    // Whatever the floors cannot separate, desirability orders.
    if (!bundle) return;
    W.state.gates = W.state.gates || {};
    const floors = ['spec', 'thermo', 'expr'];
    floors.forEach((id) => { delete W.state.gates[id]; });
    let best = null;
    for (let v = 0; v <= 95; v += 1) {
      floors.forEach((id) => { W.state.gates[id] = v; });
      if (applyGates(bundle.rows).length >= SHORTLIST) best = v;
      else break;
    }
    if (best === null || best === 0) floors.forEach((id) => { delete W.state.gates[id]; });
    else floors.forEach((id) => { W.state.gates[id] = best; });
    W.writeUrl(true);
    refresh();
  }

  /* ------------------------------------------------------------ the ladder */

  const LADDER_SORTS = [
    ['desirability', 'Desirability'],
    ['axis_affinity', 'Affinity'],
    ['axis_specificity', 'Specificity'],
    ['axis_thermo', 'Thermostability'],
    ['axis_expression', 'Expression'],
    ['axis_nativeness', 'Nativeness'],
    ['tm_pred', 'Predicted Tm'],
    ['rank', 'Tier 1 rank'],
  ];

  function ladderOrder(rows, key) {
    const flip = lowerIsBetter(key) || key === 'rank';
    return rows.slice().sort((a, b) => {
      const x = a[key], y = b[key];
      const xn = (x === null || x === undefined);
      const yn = (y === null || y === undefined);
      if (xn && yn) return (a.rank || 0) - (b.rank || 0);
      if (xn) return 1;
      if (yn) return -1;
      return flip ? x - y : y - x;
    });
  }

  function drawLadder(rows) {
    const node = el('ladder');
    if (!node || !bundle) return;
    const key = el('laddersort').value || W.state.ladder || 'desirability';
    const want = Number(el('laddern').value || 32);
    const pool = ladderOrder(rows, key).slice(0, want);

    if (!pool.length) {
      node.innerHTML = '<p class="laddempty">No candidate passes the current '
        + 'filters. Loosen a gate or reset the filters.</p>';
      el('ladderkey').textContent = '';
      return;
    }

    const medians = {};
    HERO_AXES.forEach(([field]) => { medians[field] = medianOf(field); });

    let html = '<div class="lrow lhead"><span class="pos">#</span>'
      + '<span>Candidate</span>';
    HERO_AXES.forEach(([, short, long]) => {
      html += `<span title="${long} percentile within this target">${short}</span>`;
    });
    html += '<span>Tm &deg;C</span><span>Flags</span></div>';

    pool.forEach((r, i) => {
      if (i === SHORTLIST && want > SHORTLIST) {
        html += `<div class="cutline">shortlist of ${SHORTLIST} stops here</div>`;
      }
      const inShort = i < SHORTLIST;
      html += `<div class="lrow${inShort ? ' inshort' : ''}" tabindex="0"`
        + ` role="button" data-id="${W.escape(r.id)}"`
        + ` aria-label="Open ${W.escape(r.id)}">`
        + `<span class="pos">${i + 1}</span>`
        + `<span class="cid">${W.escape(r.id)}</span>`;
      HERO_AXES.forEach(([field, , long]) => {
        const v = r[field];
        if (v === null || v === undefined) {
          html += `<span class="track pending" title="${long} has not been `
            + `computed for this candidate"></span>`;
          return;
        }
        const med = medians[field];
        html += `<span class="track" title="${long} ${v.toFixed(1)}th percentile`
          + ` within this target">`
          + `<i class="${W.pctClass(v)}" style="width:${Math.max(2, v)}%"></i>`
          + (med === null ? '' : `<u style="left:${med}%"></u>`)
          + `<b>${Math.round(v)}</b></span>`;
      });
      const tm = r.tm_pred === null || r.tm_pred === undefined
        ? '<span class="pend">&ndash;</span>'
        : `${r.tm_pred.toFixed(1)}<span class="err">&plusmn;${Math.round(r.tm_err || 0)}</span>`;
      html += `<span class="tmcell">${tm}</span><span class="flags">`
        + `<i class="pip ${r.n_cdr_liabilities ? 'bad' : 'ok'}" title="`
        + `${r.n_cdr_liabilities || 0} liability motifs inside a CDR"></i>`
        + `<i class="pip ${r.hallmark_n === 4 ? 'ok' : 'warn'}" title="`
        + `camelid hallmark ${r.hallmark_n === null ? 'unknown' : r.hallmark_n} of 4"></i>`
        + `<i class="pip ${r.experimental_badge ? 'ok' : ''}" title="`
        + `${r.experimental_badge ? 'a near-identical sequence has a PDB entry'
            : 'no near-identical experimental structure'}"></i>`
        + `</span></div>`;
    });
    node.innerHTML = html;

    node.querySelectorAll('.lrow[data-id]').forEach((row) => {
      row.addEventListener('click', () => openDrawer(row.dataset.id));
      row.addEventListener('keydown', (ev) => {
        if (ev.key === 'Enter' || ev.key === ' ') {
          ev.preventDefault();
          openDrawer(row.dataset.id);
        }
      });
    });

    const pending = HERO_AXES.filter(([field]) =>
      !bundle.rows.some((r) => r[field] !== null && r[field] !== undefined));
    el('ladderkey').innerHTML = 'tick = population median'
      + (pending.length ? ` &middot; hatched = ${pending
          .map(([, , long]) => long.toLowerCase()).join(' and ')} not computed` : '');
  }

  /* ------------------------------------------------------ the trade-off plane */

  const PLANE_FIELDS = [
    ['axis_specificity', 'Specificity percentile'],
    ['axis_thermo', 'Thermostability percentile'],
    ['axis_expression', 'Expression percentile'],
    ['axis_nativeness', 'Nativeness percentile'],
    ['axis_affinity', 'Affinity percentile'],
    ['tm_pred', 'Predicted Tm, deg C'],
    ['psh', 'Patch surface hydrophobicity'],
    ['cdr3_len', 'CDR3 length, aa'],
    ['desirability', 'Desirability'],
  ];
  let planeBox = null;

  function extentOf(field) {
    if (field.startsWith('axis_')) return [0, 100];
    let lo = Infinity, hi = -Infinity;
    (bundle ? bundle.rows : []).forEach((r) => {
      const v = r[field];
      if (typeof v !== 'number' || !isFinite(v)) return;
      if (v < lo) lo = v;
      if (v > hi) hi = v;
    });
    if (!isFinite(lo)) return null;
    if (lo === hi) { lo -= 1; hi += 1; }
    const pad = (hi - lo) * 0.05;
    return [lo - pad, hi + pad];
  }

  function drawPlane(rows) {
    const canvas = el('planecanvas');
    if (!canvas || !bundle) return;
    const xf = el('planex').value, yf = el('planey').value;
    const dpr = window.devicePixelRatio || 1;
    const width = canvas.clientWidth || 620, height = 352;
    canvas.width = width * dpr; canvas.height = height * dpr;
    const g = canvas.getContext('2d');
    g.setTransform(dpr, 0, 0, dpr, 0, 0);
    g.clearRect(0, 0, width, height);
    g.font = '10px Inter, system-ui, sans-serif';

    const ex = extentOf(xf), ey = extentOf(yf);
    if (!ex || !ey) {
      g.fillStyle = W.PALETTE.ink3;
      g.textAlign = 'center';
      g.fillText(`${!ex ? el('planex').selectedOptions[0].text
        : el('planey').selectedOptions[0].text} has not been computed for `
        + `${W.state.target} yet.`, width / 2, height / 2);
      el('planeread').textContent = '';
      return;
    }
    const MT = 34, MR = 42, MB = 38, ML = 52;
    const pw = width - ML - MR, ph = height - MT - MB;
    const sx = (v) => ML + ((v - ex[0]) / (ex[1] - ex[0])) * pw;
    const sy = (v) => MT + ph - ((v - ey[0]) / (ey[1] - ey[0])) * ph;

    g.strokeStyle = W.PALETTE.rule; g.lineWidth = 1;
    for (let i = 0; i <= 4; i++) {
      const yv = ey[0] + (ey[1] - ey[0]) * i / 4;
      const xv = ex[0] + (ex[1] - ex[0]) * i / 4;
      g.beginPath(); g.moveTo(ML, Math.round(sy(yv)) + 0.5);
      g.lineTo(ML + pw, Math.round(sy(yv)) + 0.5); g.stroke();
      g.beginPath(); g.moveTo(Math.round(sx(xv)) + 0.5, MT);
      g.lineTo(Math.round(sx(xv)) + 0.5, MT + ph); g.stroke();
      g.fillStyle = W.PALETTE.ink3;
      g.textAlign = 'right'; g.textBaseline = 'middle';
      g.fillText(yv.toFixed(ey[1] - ey[0] > 20 ? 0 : 2), ML - 7, sy(yv));
      g.textAlign = 'center'; g.textBaseline = 'top';
      g.fillText(xv.toFixed(ex[1] - ex[0] > 20 ? 0 : 2), sx(xv), MT + ph + 6);
    }

    const plotted = bundle.rows.filter((r) =>
      typeof r[xf] === 'number' && typeof r[yf] === 'number');
    const passing = new Set(rows.map((r) => r.id));

    // Marginal histograms in the top and right gutters: the distribution the
    // ladder cannot show.
    const bins = (field, n, e) => {
      const out = new Array(n).fill(0);
      plotted.forEach((r) => {
        const t = (r[field] - e[0]) / (e[1] - e[0]);
        out[Math.min(n - 1, Math.max(0, Math.floor(t * n)))] += 1;
      });
      return out;
    };
    const hx = bins(xf, 26, ex), hy = bins(yf, 20, ey);
    const mx = Math.max(...hx), my = Math.max(...hy);
    g.fillStyle = W.PALETTE.dominated;
    hx.forEach((count, i) => {
      const bw = pw / hx.length, h = (count / mx) * (MT - 12);
      g.fillRect(ML + i * bw + 0.5, MT - 4 - h, Math.max(1, bw - 1), h);
    });
    hy.forEach((count, i) => {
      const bh = ph / hy.length, w = (count / my) * (MR - 12);
      g.fillRect(ML + pw + 5, MT + i * bh + 0.5, w, Math.max(1, bh - 1));
    });

    let inBox = 0;
    const within = (r) => {
      if (!planeBox) return true;
      const x = sx(r[xf]), y = sy(r[yf]);
      return x >= Math.min(planeBox[0], planeBox[2]) && x <= Math.max(planeBox[0], planeBox[2])
        && y >= Math.min(planeBox[1], planeBox[3]) && y <= Math.max(planeBox[1], planeBox[3]);
    };
    plotted.forEach((r) => {
      const keep = passing.has(r.id) && within(r);
      if (keep) inBox += 1;
      g.beginPath();
      g.arc(sx(r[xf]), sy(r[yf]), r.on_front ? 3.1 : 2.3, 0, 6.2832);
      if (r.on_front) {
        g.fillStyle = W.PALETTE.accent; g.globalAlpha = keep ? 0.92 : 0.14; g.fill();
      } else {
        g.strokeStyle = W.PALETTE.dominated; g.lineWidth = 1;
        g.globalAlpha = keep ? 0.75 : 0.12; g.stroke();
      }
      g.globalAlpha = 1;
    });

    // The frontier of this PAIR, which is not the four-axis Pareto front. The
    // readout says so, because the two are easy to confuse and one of them is
    // a much weaker claim.
    const dirX = lowerIsBetter(xf) ? -1 : 1, dirY = lowerIsBetter(yf) ? -1 : 1;
    const sorted = plotted.slice().sort((a, b) => (b[xf] - a[xf]) * dirX);
    const front = [];
    let best = -Infinity;
    sorted.forEach((r) => {
      const yv = r[yf] * dirY;
      if (yv > best) { best = yv; front.push(r); }
    });
    g.strokeStyle = W.PALETTE.accent; g.lineWidth = 1.5; g.globalAlpha = 0.85;
    g.beginPath();
    front.forEach((r, i) => {
      const x = sx(r[xf]), y = sy(r[yf]);
      if (i === 0) g.moveTo(x, y);
      else { g.lineTo(x, sy(front[i - 1][yf])); g.lineTo(x, y); }
    });
    g.stroke(); g.globalAlpha = 1;

    if (planeBox) {
      g.strokeStyle = W.PALETTE.accent; g.lineWidth = 1; g.setLineDash([4, 3]);
      const bx = Math.min(planeBox[0], planeBox[2]), by = Math.min(planeBox[1], planeBox[3]);
      const bw = Math.abs(planeBox[2] - planeBox[0]), bh = Math.abs(planeBox[3] - planeBox[1]);
      g.strokeRect(bx, by, bw, bh); g.setLineDash([]);
      g.fillStyle = W.PALETTE.accent; g.globalAlpha = 0.06;
      g.fillRect(bx, by, bw, bh); g.globalAlpha = 1;
    }

    g.fillStyle = W.PALETTE.ink;
    g.font = '600 10.5px Inter, system-ui, sans-serif';
    g.textAlign = 'center'; g.textBaseline = 'alphabetic';
    g.fillText(el('planex').selectedOptions[0].text, ML + pw / 2, height - 5);
    g.save(); g.translate(12, MT + ph / 2); g.rotate(-Math.PI / 2);
    g.fillText(el('planey').selectedOptions[0].text, 0, 0); g.restore();

    el('planeread').innerHTML = `${W.fmt(plotted.length)} of `
      + `${W.fmt(bundle.rows.length)} ${W.state.target} candidates have both `
      + `values${planeBox ? `, ${W.fmt(inBox)} inside the box` : ''}. The `
      + `staircase holds ${W.fmt(front.length)} candidates that nothing beats on `
      + 'both of these two axes, which is a smaller claim than the four-axis '
      + 'Pareto front the table marks.';
  }

  // The box sets a range on each of the two axes, which is the same brush the
  // table already filters on, so one gesture drives both.
  function bindPlane() {
    const canvas = el('planecanvas');
    let from = null;
    canvas.addEventListener('pointerdown', (ev) => {
      const rect = canvas.getBoundingClientRect();
      from = [ev.clientX - rect.left, ev.clientY - rect.top];
      canvas.setPointerCapture(ev.pointerId);
    });
    canvas.addEventListener('pointermove', (ev) => {
      if (!from) return;
      const rect = canvas.getBoundingClientRect();
      planeBox = [from[0], from[1], ev.clientX - rect.left, ev.clientY - rect.top];
      drawPlane(activeRows());
    });
    canvas.addEventListener('pointerup', () => {
      from = null;
      if (planeBox && Math.abs(planeBox[2] - planeBox[0]) < 5
          && Math.abs(planeBox[3] - planeBox[1]) < 5) {
        planeBox = null;
      }
      commitPlaneBox();
    });
  }

  function commitPlaneBox() {
    const xf = el('planex').value, yf = el('planey').value;
    [xf, yf].forEach((f) => { if (f.startsWith('axis_')) delete W.state.brush[f]; });
    if (planeBox) {
      const canvas = el('planecanvas');
      const width = canvas.clientWidth, height = 352;
      const MT = 34, MR = 42, MB = 38, ML = 52;
      const pw = width - ML - MR, ph = height - MT - MB;
      const ex = extentOf(xf), ey = extentOf(yf);
      if (ex && xf.startsWith('axis_')) {
        const a = ex[0] + ((Math.min(planeBox[0], planeBox[2]) - ML) / pw) * (ex[1] - ex[0]);
        const b = ex[0] + ((Math.max(planeBox[0], planeBox[2]) - ML) / pw) * (ex[1] - ex[0]);
        W.state.brush[xf] = [Math.max(0, a), Math.min(100, b)];
      }
      if (ey && yf.startsWith('axis_')) {
        const a = ey[0] + ((MT + ph - Math.max(planeBox[1], planeBox[3])) / ph) * (ey[1] - ey[0]);
        const b = ey[0] + ((MT + ph - Math.min(planeBox[1], planeBox[3])) / ph) * (ey[1] - ey[0]);
        W.state.brush[yf] = [Math.max(0, a), Math.min(100, b)];
      }
    }
    W.writeUrl(true);
    refresh();
  }

  /* ------------------------------------------------------------ the switcher */

  function heroView() {
    return W.state.hero === 'plane' ? 'plane' : 'ladder';
  }

  function drawHero(rows) {
    const view = heroView();
    // `hidden` is not enough here: a stylesheet `display` beats the user
    // agent's `[hidden] { display: none }`, and .panelbar is display:flex, so
    // the ladder's controls stayed on screen over the plane. Set display.
    const show = (id, on, mode) => {
      const node = el(id);
      node.hidden = !on;
      node.style.display = on ? (mode || '') : 'none';
    };
    show('ladderbar', view === 'ladder', 'flex');
    show('ladder', view === 'ladder');
    show('planewrap', view === 'plane');
    el('heroview').querySelectorAll('button').forEach((b) => {
      b.setAttribute('aria-pressed', String(b.dataset.view === view));
    });
    drawGates();
    if (view === 'ladder') drawLadder(rows);
    else drawPlane(rows);
    const front = bundle ? bundle.rows.filter((r) => r.on_front).length : 0;
    el('explorerhint').textContent = view === 'ladder'
      ? `the ${SHORTLIST} that would be ordered, and what removed the rest`
      : `${W.fmt(front)} candidates on the four-axis Pareto front`;
  }

  let brushTimer = null;
  function refreshAfterBrush() {
    clearTimeout(brushTimer);
    brushTimer = setTimeout(() => {
      const rows = activeRows();
      if (table && tableReady) table.replaceData(rows);
      renderCaption(rows);
    }, 60);
  }

  /* ------------------------------------------- 3.3 where one candidate sits

     The strips are the parcoords idea with the lines taken out: one row per
     axis, every candidate a dot, and a needle on the one the drawer is open
     on. In the drawer the subject is already a single candidate, which is the
     question this answers. */

  function drawAxisStrips(d) {
    const host = el('axisstrips');
    if (!host || !bundle) return;
    const rows = bundle.rows;
    host.innerHTML = '<h3 class="stripshead">Where this candidate sits</h3>'
      + HERO_AXES.map(([field, , long]) =>
        `<div class="striprow"><span class="sname">${long}</span>`
        + `<span class="canvasbox"><canvas height="38" data-axis="${field}"></canvas></span>`
        + `</div>`).join('')
      + `<p class="stripsnote">Every ${W.escape(W.state.target)} candidate is a `
      + `dot. The marked one is this candidate. Percentiles are within this `
      + `target only.</p>`;

    host.querySelectorAll('canvas').forEach((canvas) => {
      const field = canvas.dataset.axis;
      const dpr = window.devicePixelRatio || 1;
      const width = canvas.clientWidth || 300, height = 38;
      canvas.width = width * dpr; canvas.height = height * dpr;
      const g = canvas.getContext('2d');
      g.setTransform(dpr, 0, 0, dpr, 0, 0);
      g.clearRect(0, 0, width, height);
      const ML = 6, pw = width - 12, mid = height / 2;

      g.strokeStyle = W.PALETTE.rule; g.lineWidth = 1;
      g.beginPath(); g.moveTo(ML, mid + 0.5); g.lineTo(ML + pw, mid + 0.5); g.stroke();

      const have = rows.filter((r) => typeof r[field] === 'number');
      if (!have.length) {
        g.fillStyle = W.PALETTE.ink3;
        g.font = '10px Inter, system-ui, sans-serif';
        g.textAlign = 'left'; g.textBaseline = 'middle';
        g.fillText('not computed for this target', ML + 2, mid);
        return;
      }
      have.forEach((r) => {
        const x = ML + (r[field] / 100) * pw;
        // Deterministic jitter, so the strip does not reshuffle on redraw.
        // The candidate number drives it rather than a hash of the whole id:
        // ids share a long target prefix, and folding that in collapsed the
        // spread into two bands instead of filling the strip.
        const n = parseInt(String(r.id).replace(/\D/g, ''), 10) || 0;
        const jitter = ((n * 2654435761) % 10007) / 10007;
        const y = mid + (jitter - 0.5) * (height - 11);
        g.beginPath(); g.arc(x, y, 1.8, 0, 6.2832);
        g.fillStyle = W.PALETTE.dominated; g.globalAlpha = 0.6; g.fill();
        g.globalAlpha = 1;
      });
      const mine = d ? d[field] : null;
      if (typeof mine === 'number') {
        const x = ML + (mine / 100) * pw;
        g.strokeStyle = W.PALETTE.accent; g.lineWidth = 1.5;
        g.beginPath(); g.moveTo(x, 2); g.lineTo(x, height - 2); g.stroke();
        g.fillStyle = '#FFFFFF';
        g.beginPath(); g.arc(x, mid, 4, 0, 6.2832); g.fill(); g.stroke();
      } else {
        g.fillStyle = W.PALETTE.ink3;
        g.font = '10px Inter, system-ui, sans-serif';
        g.textAlign = 'right'; g.textBaseline = 'middle';
        g.fillText('pending for this candidate', ML + pw - 2, mid);
      }
    });
  }

  /* --------------------------------------------------- 3.4 the column chooser */

  function buildChooser() {
    const chosen = new Set(W.state.cols || REG.defaults);
    const body = el('colbody');
    body.innerHTML = (REG.groups || []).map((group) => {
      const items = (REG.columns || []).filter((c) => c.group === group.id);
      const rows = items.map((col) => {
        const pending = (bundle && bundle.columns_pending || []).includes(col.key);
        return `<label class="${pending ? 'pending' : ''}" data-key="${col.key}"
            title="${W.escape(col.help)}">
          <input type="checkbox" value="${col.key}"
            ${chosen.has(col.key) ? 'checked' : ''}>
          <span class="cname">${W.escape(col.label)}${pending ? ' (pending)' : ''}</span>
          <span class="cunit">${W.escape(col.unit || '')}</span>
        </label>`;
      }).join('');
      return `<div class="colgroup"><h4>${W.escape(group.label)} &middot; ${items.length}</h4>${rows}</div>`;
    }).join('');

    body.querySelectorAll('input[type="checkbox"]').forEach((box) => {
      box.addEventListener('change', () => {
        const keys = [...body.querySelectorAll('input:checked')].map((b) => b.value);
        // Keep the registry order, not click order, so a saved view is stable.
        W.state.cols = (REG.columns || []).map((c) => c.key)
          .filter((k) => keys.includes(k));
        W.writeUrl(true);
        buildTable();
        refresh();
      });
    });
  }

  /* --------------------------------------------------------- 3.5 the drawer */

  /* The verdict the vacuous tests reached for an axis on the target in view.
     "A non-discriminating axis must be labelled as such wherever it appears"
     is the test's own instruction, and until this was wired the label lived
     only on the Control set tab: a reader here saw a percentile with nothing
     to say it does not separate a real candidate from a shuffled one. */
  let axisVerdicts = null;

  function axisVerdict(axis) {
    if (axisVerdicts === null) {
      const node = document.getElementById('axisverdicts');
      try {
        axisVerdicts = node ? JSON.parse(node.textContent) || {} : {};
      } catch (err) {
        axisVerdicts = {};
      }
    }
    const all = axisVerdicts;
    const byTarget = all[axis];
    if (!byTarget) return null;
    return byTarget[W.state.target] || null;
  }

  function tile(label, value, unit, pct, axis) {
    const width = pct === null || pct === undefined ? 0 : Math.max(2, pct);
    const missing = value === null || value === undefined;
    const shown = missing
      ? '<span style="font-size:13px;color:var(--ink-3)">pending</span>'
      : W.escape(value);
    // No unit beside "pending": there is no number for it to qualify.
    const unitSpan = missing ? '' : `<span class="tunit">${W.escape(unit || '')}</span>`;
    const verdict = axis ? axisVerdict(axis) : null;
    // Four states, not two. "unconfirmed" and "untestable" are things the
    // control set measured and could not rule on, and flattening either into
    // the unflagged case would present a weak axis as a working one.
    const FLAGS = {
      vacuous: ['non-discriminating',
        'This axis does not separate real candidates from CDR-shuffled '
        + 'controls, so its ordering carries no evidence.'],
      unconfirmed: ['weak',
        'This axis separates real candidates from CDR-shuffled controls on '
        + 'the unpaired test but not when each control is compared with the '
        + 'candidate it was shuffled from. Treat its ordering as suggestive.'],
      untestable: ['untested',
        'The CDR-shuffled control preserves amino-acid composition, which is '
        + 'most of what this axis reads, so the control cannot test it either '
        + 'way. This is not a failure of the axis.'],
    };
    const label = FLAGS[verdict];
    const flag = label
      ? `<div class="tflag${verdict === 'vacuous' ? '' : ' soft'}" `
        + `title="${W.escape(label[1])}">${label[0]}</div>`
      : '';
    return `<div class="tile"><div class="tlabel">${W.escape(label)}</div>
      <div class="tvalue">${shown}${unitSpan}</div>
      <div class="trule"><i style="width:${width}%"></i></div>${flag}</div>`;
  }

  function kvTable(pairs) {
    return '<table class="kv"><tbody>' + pairs.map(([key, value, title]) => {
      const pending = value === null || value === undefined || value === '';
      return `<tr><th scope="row"${title ? ` title="${W.escape(title)}"` : ''}>`
        + `${W.escape(key)}</th><td class="${pending ? 'pending' : ''}">`
        + `${pending ? 'pending' : W.escape(value)}</td></tr>`;
    }).join('') + '</tbody></table>';
  }

  /* The drawer slides over the page, so its open state is a class rather than
     `hidden`: a display:none element cannot animate, and the panel has to
     transition both ways. `inert` and aria-hidden keep the closed panel out of
     the tab order and off the accessibility tree, which `hidden` was doing. */
  function setDrawerOpen(open) {
    const node = el('drawer');
    if (!node) return;
    node.classList.toggle('open', open);
    node.setAttribute('aria-hidden', String(!open));
    if (open) node.removeAttribute('inert');
    else node.setAttribute('inert', '');
  }

  function drawerIsOpen() {
    const node = el('drawer');
    return Boolean(node && node.classList.contains('open'));
  }

  async function openDrawer(id) {
    W.state.id = id;
    W.writeUrl(true);
    setDrawerOpen(true);
    el('drawerid').textContent = id;
    document.querySelectorAll('.tabulator-row').forEach(
      (r) => r.classList.remove('whittle-selected'));
    const rowEl = table && table.getRow(id);
    if (rowEl) rowEl.getElement().classList.add('whittle-selected');

    const detail = await W.loadJson(`/api/candidate/${encodeURIComponent(id)}.json`);
    if (!detail) { el('panel-metrics').innerHTML = '<p>Could not load this candidate.</p>'; return; }
    currentDetail = detail;
    viewerMode = 'model';

    el('drawerchips').innerHTML = (detail.on_front
      ? '<span class="chip good">On front</span> ' : '')
      + (detail.multi_target ? '<span class="chip warn">Multi-target</span> ' : '')
      + (detail.target_detail && detail.target_detail.ecd_scoped
        ? '<span class="chip warn">ECD-scoped</span>' : '');

    renderStructurePanel(detail);
    renderMetricsPanel(detail);
    renderSequencePanel(detail);
    renderNeighboursPanel(detail);
    selectPanel(W.state.panel || 'structure');
  }

  // Which of the two structures the Structure panel is showing. Reset per
  // candidate, because a candidate with no co-fold has no Complex to stay on.
  let viewerMode = 'model';

  /* The Structure panel shows one of two things, and the toggle, the first
     render and the deferred render for a panel that was hidden all have to
     agree about which. One function owns the decision, so switching to
     Complex, leaving the sub-tab and coming back cannot quietly reload the
     monomer while the toggle still reads Complex. */
  async function loadStructureView(d) {
    const container = el('viewer');
    if (!container || !window.WViewer) return;
    const complex = viewerMode === 'complex' && d.complex_url;
    const url = complex ? d.complex_url : d.structure_url;
    if (!url) return;
    delete container.dataset.pending;
    if (complex) {
      const profile = await loadInteractions(d);
      // Guard against the reader switching away while the fetch was in
      // flight: drawing the complex into a viewer now showing the monomer
      // would leave the toggle and the picture disagreeing.
      if (viewerMode !== 'complex') return;
      await showComplex(d, profile);
      renderInteractions(d, profile);
    } else {
      await showModel(d);
      const panel = el('interactions');
      if (panel) panel.hidden = true;
    }
    writeViewerCaption(d);
  }

  /* The predicted model, with the three CDRs picked out.

     This goes through MolViewSpec rather than a plain load because a plain
     load cannot colour anything: the whole domain came back one colour, and
     the caption said "CDR1 / 2 / 3" over a picture that distinguished none of
     them. If the MVS extension is unavailable the plain load is still better
     than an empty panel, so it is the fallback. */
  async function showModel(d) {
    const container = el('viewer');
    if (!container || !d.structure_url) return;
    const result = await window.WViewer.showMvs(container, (builder) => {
      builder.canvas({ background_color: VIEWER_BG });
      const structure = builder
        .download({ url: d.structure_url })
        .parse({ format: 'pdb' })
        .modelStructure();
      structure.component({
        selector: frameworkRanges(MODEL_CHAIN, IMGT_CDR, IMGT_LAST) })
        .representation({ type: 'cartoon' })
        .color({ color: NANOBODY_COLOUR });
      for (const key of ['cdr1', 'cdr2', 'cdr3']) {
        structure.component({ selector: cdrSelector(MODEL_CHAIN, IMGT_CDR[key]) })
          .representation({ type: 'cartoon' })
          .color({ color: CDR_COLOURS[key] });
      }
    });
    if (!result || !result.ok) {
      await window.WViewer.show(container,
        [{ url: d.structure_url, format: 'pdb', binary: false }]);
    }
  }

  /* The complex, with the residues PLIP found in contact picked out.

     Mol* 4.9's MolViewSpec has no primitives node, so these cannot be drawn
     as dashed lines between the two atoms. Colouring the contacting residues
     is more legible at this size in any case, and the list below carries the
     atom names and distances. */
  async function showComplex(d, profile) {
    const container = el('viewer');
    const byType = highlightResidues(profile && profile.interactions);
    // Chain N is numbered 1..len(sequence) here, not IMGT, so the CDRs are
    // found by substring rather than at fixed positions.
    const spans = cdrSpansFromSequence(d);
    const last = (d.sequence || '').length;
    await window.WViewer.showMvs(container, (builder) => {
      // Loading an MVS state resets the canvas, so the background has to be
      // part of the state rather than set once on the plugin.
      builder.canvas({ background_color: VIEWER_BG });
      const structure = builder
        .download({ url: d.complex_url })
        .parse({ format: 'pdb' })
        .modelStructure();
      // Disjoint components, one per chain, rather than a `polymer` cartoon
      // with the nanobody drawn over it: two cartoons on the same geometry
      // z-fight and the upper one does not reliably win, which left the whole
      // complex a single grey. The CDRs are carved out of chain N for the
      // same reason.
      structure.component({ selector: ANTIGEN_CHAINS.map(
        (c) => ({ auth_asym_id: c })) })
        .representation({ type: 'cartoon' })
        .color({ color: ANTIGEN_COLOUR });
      const framework = last
        ? frameworkRanges(COMPLEX_NANOBODY_CHAIN, spans, last)
        : [{ auth_asym_id: COMPLEX_NANOBODY_CHAIN }];
      structure.component({ selector: framework })
        .representation({ type: 'cartoon' })
        .color({ color: NANOBODY_COLOUR });
      for (const key of ['cdr1', 'cdr2', 'cdr3']) {
        if (!spans[key]) continue;
        structure.component({
          selector: cdrSelector(COMPLEX_NANOBODY_CHAIN, spans[key]) })
          .representation({ type: 'cartoon' })
          .color({ color: CDR_COLOURS[key] });
      }
      // Weakest first, so a residue drawn twice ends on its strongest type.
      for (let i = IX_TYPES.length - 1; i >= 0; i -= 1) {
        const spec = IX_TYPES[i];
        const residues = byType[spec.key];
        if (!residues || !residues.length) continue;
        structure.component({ selector: residues })
          .representation({ type: 'ball_and_stick' })
          .color({ color: spec.colour });
      }
    });
  }

  function renderInteractions(d, profile) {
    const panel = el('interactions');
    if (!panel) return;
    panel.hidden = false;
    const records = (profile && profile.interactions) || [];
    if (!records.length) {
      panel.innerHTML = '<p class="ixempty">No interactions were detected '
        + 'across this interface.</p>';
      return;
    }
    const counts = {};
    for (const r of records) counts[r.type] = (counts[r.type] || 0) + 1;
    const chips = IX_TYPES.filter((s) => counts[s.key])
      .map((s) => `<span class="ixchip" style="color:${s.colour}">`
        + `${counts[s.key]} `
        + `${W.escape(counts[s.key] === 1 ? s.one : s.many)}</span>`)
      .join('');
    const order = Object.fromEntries(IX_TYPES.map((x, i) => [x.key, i]));
    const rows = records.slice().sort((a, b) => (order[a.type] - order[b.type])
        || ((a.distance || 0) - (b.distance || 0)))
      .map((r) => {
        const spec = IX_BY_KEY[r.type] || { colour: 'var(--ink-3)' };
        const side = (s) => `${W.escape(s.resname)}${s.resnum}`
          + `<span class="ixmono">${s.atom ? ' ' + W.escape(s.atom) : ''}</span>`;
        return `<tr><td style="width:14px"><span class="ixdot" `
          + `style="background:${spec.colour}"></span></td>`
          + `<td>${side(r.antigen)}</td>`
          + `<td style="color:${NANOBODY_COLOUR}">${side(r.nanobody)}</td>`
          + `<td class="ixmono" style="text-align:right">`
          + `${r.distance === null ? '' : r.distance.toFixed(2) + ' \u00c5'}</td></tr>`;
      }).join('');
    panel.innerHTML = `<div class="ixhead">${chips}</div>`
      + `<table><tbody>${rows}</tbody></table>`
      + `<p class="ixnote">Antigen residue, then nanobody residue. Detected `
      + `with PLIP on the highest-ipTM sample, which is the model shown above. `
      + `These are interactions of a predicted pose, not a measured structure.</p>`;
  }

  // Interaction types, their colour in the viewer and in the list, and the
  // order they are drawn in: a residue making several gets the strongest.
  const IX_TYPES = [
    { key: 'saltbridge', one: 'salt bridge', many: 'salt bridges',
      colour: '#E4698F' },
    { key: 'hbond', one: 'hydrogen bond', many: 'hydrogen bonds',
      colour: '#4DA3FF' },
    { key: 'pistacking', one: 'pi stack', many: 'pi stacks',
      colour: '#5EC2A0' },
    { key: 'pication', one: 'pi cation', many: 'pi cations',
      colour: '#5EC2A0' },
    { key: 'hydrophobic', one: 'hydrophobic contact',
      many: 'hydrophobic contacts', colour: '#E0A458' },
  ];
  const IX_BY_KEY = Object.fromEntries(IX_TYPES.map((x) => [x.key, x]));
  // A, B and C cover the monomers and the two homotrimers; naming a chain
  // that is not in the file is harmless.
  const ANTIGEN_CHAINS = ['A', 'B', 'C'];
  const ANTIGEN_COLOUR = '#9AA4B2';
  const NANOBODY_COLOUR = '#44506A';     // the framework, so the CDRs carry
  const VIEWER_BG = '#FFFFFF';
  const REFERENCE_COLOUR = '#C3CBD9';   // the experimental chain

  /* The three CDRs, in colours that step in lightness as well as hue so the
     key still reads in greyscale and under colour-vision deficiency
     (DESIGN 2). They are also the colours the legend prints, because a key in
     a different colour from the picture is worse than no key. */
  const CDR_COLOURS = { cdr1: '#2B7FB8', cdr2: '#E0A458', cdr3: '#9B2242' };
  const CDR_LABELS = { cdr1: 'CDR1', cdr2: 'CDR2', cdr3: 'CDR3' };

  /* Two files, two numbering schemes, and getting them the wrong way round
     colours the wrong residues.

     - The predicted model and its superposed copy are IMGT-numbered in chain
       H, which is the point of IMGT: the CDRs sit at fixed positions whatever
       the sequence does. Verified against the run: IMGT 27-38 and 56-65
       reproduce ANARCI's CDR1 and CDR2 exactly, and 105-117 covers the
       modelled part of CDR3 (NanoBodyBuilder2 does not always model every
       residue of a long one).
     - Chain N of the co-folded complex is numbered 1..len(sequence) and
       carries the candidate sequence exactly, so there the CDRs are found by
       substring, mirroring pipeline/interface.py:cdr_spans. */
  const IMGT_CDR = { cdr1: [27, 38], cdr2: [56, 65], cdr3: [105, 117] };
  const MODEL_CHAIN = 'H';
  const COMPLEX_NANOBODY_CHAIN = 'N';
  const IMGT_LAST = 128;

  function cdrSpansFromSequence(d) {
    const sequence = d.sequence || '';
    const spans = {};
    let cursor = 0;
    for (const key of ['cdr1', 'cdr2', 'cdr3']) {
      const motif = (d[key] || '').trim();
      if (!motif) continue;
      let at = sequence.indexOf(motif, cursor);
      if (at < 0) at = sequence.indexOf(motif);
      if (at < 0) continue;
      spans[key] = [at + 1, at + motif.length];
      cursor = at + motif.length;
    }
    return spans;
  }

  /* Everything that is not a CDR, as explicit ranges. The CDR cartoon cannot
     simply be drawn over a whole-chain one: two cartoons on the same geometry
     z-fight and the upper one does not reliably win, which is the bug that
     left the complex a single grey. Disjoint components avoid it. */
  function frameworkRanges(chain, spans, last) {
    const ordered = ['cdr1', 'cdr2', 'cdr3']
      .map((key) => spans[key]).filter(Boolean)
      .sort((a, b) => a[0] - b[0]);
    const out = [];
    let from = 1;
    for (const [begin, end] of ordered) {
      if (begin > from) {
        out.push({ auth_asym_id: chain, beg_auth_seq_id: from,
          end_auth_seq_id: begin - 1 });
      }
      from = end + 1;
    }
    if (from <= last) {
      out.push({ auth_asym_id: chain, beg_auth_seq_id: from,
        end_auth_seq_id: last });
    }
    return out;
  }

  function cdrSelector(chain, span) {
    return [{ auth_asym_id: chain, beg_auth_seq_id: span[0],
      end_auth_seq_id: span[1] }];
  }

  function cdrLegend() {
    return '<span class="cdrlegend">' + ['cdr1', 'cdr2', 'cdr3'].map((key) =>
      `<span class="cdrkey"><i style="background:${CDR_COLOURS[key]}"></i>`
      + `${CDR_LABELS[key]}</span>`).join('') + '</span>';
  }

  // One fetch per candidate, kept while the drawer stays on it.
  let interactionsFor = { id: null, data: null };

  async function loadInteractions(d) {
    if (!d.interactions_url) return null;
    if (interactionsFor.id === d.id) return interactionsFor.data;
    const data = await W.loadJson(d.interactions_url);
    interactionsFor = { id: d.id, data };
    return data;
  }

  /* Residues to highlight, keyed by interaction type, strongest wins.
     Returned per chain so the MVS selector can name them exactly. */
  function highlightResidues(records) {
    const claimed = new Map();          // "chain:resnum" -> type key
    const rank = Object.fromEntries(IX_TYPES.map((x, i) => [x.key, i]));
    for (const rec of records || []) {
      for (const side of ['antigen', 'nanobody']) {
        const s = rec[side];
        if (!s) continue;
        const id = `${s.chain}:${s.resnum}`;
        const prev = claimed.get(id);
        if (prev === undefined || rank[rec.type] < rank[prev]) {
          claimed.set(id, rec.type);
        }
      }
    }
    const byType = {};
    for (const [id, type] of claimed) {
      const [chain, resnum] = id.split(':');
      (byType[type] = byType[type] || []).push(
        { auth_asym_id: chain, auth_seq_id: Number(resnum) });
    }
    return byType;
  }

  function bestCofold(d) {
    const rows = (d.cofolds || []).filter((r) => r.iptm !== null
                                                 && r.iptm !== undefined);
    if (!rows.length) return null;
    return rows.reduce((a, b) => (b.iptm > a.iptm ? b : a));
  }

  function writeViewerCaption(d) {
    const left = el('vcapleft');
    const right = el('vcapright');
    if (!left || !right) return;
    if (viewerMode === 'complex') {
      const best = bestCofold(d);
      // ipTM belongs on the picture, not three panels away: a complex at 0.3
      // is a guess at the pose, and the reader has to see that while they are
      // looking at it. Spec 6.3 and the limitations statement.
      el('vcapleft').textContent = best && best.iptm !== null
        ? `co-folded complex · ipTM ${best.iptm.toFixed(2)}`
          + (best.iptm < 0.5 ? ' (pose unreliable)' : '')
        : 'co-folded complex';
      right.innerHTML = `<span>chain N = nanobody</span><span>${cdrLegend()}</span>`;
      return;
    }
    el('vcapleft').textContent = d.nearest_pdb
      ? `model vs ${d.nearest_pdb}${d.nearest_tm ? ` · TM ${d.nearest_tm.toFixed(2)}` : ''}`
      : '';
    // The key is printed in the colours the viewer draws the CDRs in, rather
    // than naming three colours the reader has to guess at.
    right.innerHTML = d.cdr1 ? cdrLegend() : '';
  }

  function showPlaceholder(text) {
    const node = el('viewerplaceholder');
    if (!node) return;
    node.hidden = false;
    node.style.display = 'grid';
    node.textContent = text;
  }

  function syncViewerMode(d) {
    const row = el('vmode');
    if (!row) return;
    const hasComplex = Boolean(d.complex_url);
    if (!hasComplex) viewerMode = 'model';
    row.hidden = !d.structure_url;
    row.querySelectorAll('button').forEach((btn) => {
      const mine = btn.dataset.mode;
      btn.setAttribute('aria-pressed', String(mine === viewerMode));
      btn.disabled = (mine === 'complex' && !hasComplex);
      btn.title = (mine === 'complex' && !hasComplex)
        ? 'This candidate was not taken through to co-folding.'
        : '';
    });
  }

  async function renderStructurePanel(d) {
    el('axistiles').innerHTML = [
      tile('Affinity', d.axis_affinity === null ? null : Math.round(d.axis_affinity),
        'pct', d.axis_affinity, 'affinity'),
      tile('Specificity', d.axis_specificity === null ? null : Math.round(d.axis_specificity),
        'pct', d.axis_specificity, 'specificity'),
      tile('Thermostability', d.axis_thermo === null ? null : Math.round(d.axis_thermo),
        'pct', d.axis_thermo, 'thermo'),
      tile('Expression', d.axis_expression === null ? null : Math.round(d.axis_expression),
        'pct', d.axis_expression, 'expression'),
    ].join('');

    drawAxisStrips(d);
    syncViewerMode(d);
    writeViewerCaption(d);

    const notes = [];
    if (d.model_note) notes.push(d.model_note);
    notes.push('Percentiles are within this target only and are not comparable '
      + 'between targets.');
    if (d.axis_affinity === null || d.axis_affinity_provisional) {
      notes.push('The affinity axis is a co-folding confidence prior, not an '
        + 'affinity measurement, and published benchmarks show it does not '
        + 'reliably separate cognate from non-cognate nanobody-antigen pairs.');
    }
    if (d.tm_err) {
      notes.push(`The predicted Tm carries an error bar of plus or minus `
        + `${d.tm_err.toFixed(1)} degrees C, which widens away from the training set.`);
    }
    el('structfootnote').textContent = notes.join(' ');

    await renderTmPanel(d);

    const container = el('viewer');
    // Only build a WebGL context for a panel the reader can see. Loading the
    // Structure panel's viewer while the Neighbours panel is the visible one
    // spends a context on a hidden div, and contexts are the scarce resource
    // (DESIGN 3.6).
    const structureVisible = (W.state.panel || 'structure') === 'structure';
    if (d.structure_url && window.WViewer && structureVisible) {
      await loadStructureView(d);
    } else if (d.structure_url && !structureVisible) {
      showPlaceholder('Select this tab to load the model.');
      container.dataset.pending = '1';
    } else {
      showPlaceholder('No predicted model for this candidate yet. '
        + 'The structural tier has not run.');
    }
  }

  // The only Plotly on this tab now: the Tm distribution in the drawer. The
  // hero used to pull it on page load for the parcoords panel; it is loaded on
  // first drawer open instead (DESIGN 3.6).
  async function ensurePlotly() {
    if (window.Plotly) return window.Plotly;
    if (plotlyReady) return plotlyReady;
    plotlyReady = new Promise((resolve, reject) => {
      const css = document.createElement('link');
      css.rel = 'stylesheet';
      css.href = '/static/vendor/plotly-whittle.min.css';
      document.head.appendChild(css);
      const script = document.createElement('script');
      script.src = '/static/vendor/plotly-whittle.min.js';
      script.async = true;
      script.onload = () => resolve(window.Plotly);
      script.onerror = () => reject(new Error('Plotly did not load'));
      document.head.appendChild(script);
    });
    return plotlyReady;
  }

  async function renderTmPanel(d) {
    const node = el('tmpanel');
    if (!tmDist) tmDist = await W.loadJson(`/static/data/tm_${W.state.target}.json`);
    if (!tmDist || !tmDist.n) {
      node.innerHTML = '<p style="font-size:11px;color:var(--ink-3);margin:0">'
        + 'The predicted Tm distribution appears once the thermostability model '
        + 'has run for this target.</p>';
      return;
    }
    let Plotly;
    try { Plotly = await ensurePlotly(); } catch (err) { return; }
    const traces = [{
      type: 'bar', x: tmDist.bins, y: tmDist.counts,
      marker: { color: W.PALETTE.dominated }, hoverinfo: 'x+y', name: 'all candidates',
    }];
    const shapes = [];
    const annotations = [];
    if (d.tm_pred !== null && d.tm_pred !== undefined) {
      const top = Math.max(...tmDist.counts);
      shapes.push({ type: 'line', x0: d.tm_pred, x1: d.tm_pred, y0: 0, y1: top,
        line: { color: W.PALETTE.accent, width: 2.5 } });
      if (d.tm_err) {
        // The error bar, drawn. A bare number hides the uncertainty that matters
        // most (DESIGN 3.5 item 3).
        shapes.push({ type: 'rect', x0: d.tm_pred - d.tm_err, x1: d.tm_pred + d.tm_err,
          y0: 0, y1: top, fillcolor: W.PALETTE.accent, opacity: 0.14, line: { width: 0 } });
      }
      annotations.push({ x: d.tm_pred, y: top, text: 'this one', showarrow: false,
        yanchor: 'bottom', font: { family: 'Inter', size: 10, color: W.PALETTE.accent } });
    }
    await Plotly.react(node, traces, {
      height: 150, margin: { l: 34, r: 10, t: 16, b: 30 },
      paper_bgcolor: 'rgba(0,0,0,0)', plot_bgcolor: 'rgba(0,0,0,0)',
      font: { family: 'Inter', size: 10, color: W.PALETTE.ink3 },
      xaxis: { title: { text: 'predicted apparent Tm, deg C', font: { size: 10 } },
        gridcolor: W.PALETTE.rule, zeroline: false },
      yaxis: { title: { text: 'candidates', font: { size: 10 } },
        gridcolor: W.PALETTE.rule, zeroline: false },
      shapes, annotations, showlegend: false, bargap: 0.04,
    }, { displayModeBar: false, responsive: true });
  }

  function renderMetricsPanel(d) {
    const sections = [];
    sections.push('<h3 style="font-size:11px;text-transform:uppercase;'
      + 'letter-spacing:.07em;color:var(--ink-3);margin:0 0 6px">Ranking</h3>');
    sections.push(kvTable([
      ['Rank within target', d.rank === null ? null : `${d.rank} of 500`],
      ['Pareto rank', d.pareto_rank, 'Front 1 is non-dominated across the four axes.'],
      ['Desirability', d.desirability === null ? null : d.desirability.toFixed(3),
        'Derringer-Suich geometric mean over the four axis percentiles.'],
      ['Tier reached', d.tier_reached],
    ]));
    sections.push('<h3 style="font-size:11px;text-transform:uppercase;'
      + 'letter-spacing:.07em;color:var(--ink-3);margin:16px 0 6px">Thermostability</h3>');
    sections.push(kvTable([
      ['Predicted Tm', d.tm_pred === null ? null
        : `${d.tm_pred.toFixed(1)} ± ${(d.tm_err || 0).toFixed(1)} deg C`,
        'An estimate from the WHITTLE model, not a measurement.'],
      ['Distance to nearest training neighbour', d.train_nn_distance,
        'A large distance means the model is extrapolating and the error bar widens.'],
      ['TemStaPro cross-check', d.temstapro],
      ['Upstream Tm', d.upstream_tm_pred === null ? null
        : `${d.upstream_tm_pred.toFixed(1)} deg C (uncalibrated)`,
        'From the source files. Its range reaches 101 deg C, which is not a value '
        + 'a DSF run gives a VHH in solution.'],
      ['Camelid hallmark', d.hallmark_n === null ? null
        : `${d.hallmark_n} of 4 · ${d.hallmark_detail || ''}`,
        'IMGT 42, 49, 50 and 52. A VH-type residue there leaves the former VL '
        + 'interface hydrophobic and exposed.'],
    ]));
    sections.push('<h3 style="font-size:11px;text-transform:uppercase;'
      + 'letter-spacing:.07em;color:var(--ink-3);margin:16px 0 6px">Specificity and developability</h3>');
    sections.push(kvTable([
      ['TNP patch surface hydrophobicity', d.psh],
      ['TNP patch positive charge', d.ppc],
      ['TNP patch negative charge', d.pnc],
      ['TNP metrics outside the clinical-stage range', d.tnp_flags],
      ['SWI solubility', d.swi],
      ['Protein-Sol solubility', d.protein_sol],
      ['Aggrescan3D', d.aggrescan3d],
      ['Nativeness, pseudo-perplexity', d.nativeness_ppl],
      ['OASis humanness', d.oasis_humanness],
      ['Net charge at pH 7.4', d.net_charge_74],
      ['Isoelectric point', d.pi, 'A pI near 7.4 is an aggregation flag.'],
    ]));
    if (d.cofolds && d.cofolds.length) {
      sections.push('<h3 style="font-size:11px;text-transform:uppercase;'
        + 'letter-spacing:.07em;color:var(--ink-3);margin:16px 0 6px">Co-folding, per seed</h3>');
      sections.push('<table class="kv"><thead><tr><th>Seed</th><th>ipTM</th>'
        + '<th>pTM</th><th>dG</th><th>Subunits</th></tr></thead><tbody>'
        + d.cofolds.map((f) => `<tr><td>${f.seed}</td><td>${W.fmt(f.iptm, 3)}</td>`
          + `<td>${W.fmt(f.ptm, 3)}</td><td>${W.fmt(f.prodigy_dg, 2)}</td>`
          + `<td>${f.n_subunits_contacted ?? ''}</td></tr>`).join('')
        + '</tbody></table>');
    }
    if (d.liabilities_detail && d.liabilities_detail.length) {
      sections.push('<h3 style="font-size:11px;text-transform:uppercase;'
        + 'letter-spacing:.07em;color:var(--ink-3);margin:16px 0 6px">'
        + `Liability motifs · ${d.liabilities_detail.length}</h3>`);
      sections.push('<table class="kv"><tbody>' + d.liabilities_detail.map((l) =>
        `<tr><th scope="row">${W.escape(l.kind.replace(/_/g, ' '))} at ${l.position}`
        + ` <span style="color:var(--ink-3)">(${W.escape(l.region)})</span></th>`
        + `<td style="white-space:normal">${W.escape(l.detail)}</td></tr>`).join('')
        + '</tbody></table>');
    }
    el('panel-metrics').innerHTML = sections.join('');
  }

  function numberingLabel(method) {
    if (method === 'regex_fallback') {
      return 'a regex fallback against the VHH framework consensus';
    }
    if (method === 'anarcii') return 'anarcii';
    return 'ANARCI';
  }

  function renderSequencePanel(d) {
    const seq = d.sequence || '';
    let marked = W.escape(seq);
    [['cdr1', 'c1'], ['cdr2', 'c2'], ['cdr3', 'c3']].forEach(([key, cls]) => {
      const frag = d[key];
      if (frag) marked = marked.replace(frag, `<mark class="${cls}">${frag}</mark>`);
    });
    el('panel-sequence').innerHTML =
      `<div class="seqblock">${marked}</div>`
      + `<p style="font-size:11px;color:var(--ink-3);margin:8px 0 14px">`
      + `CDR1 in blue, CDR2 in green, CDR3 in amber. IMGT definition, numbered by `
      + `${W.escape(numberingLabel(d.numbering_method))}.</p>`
      + kvTable([
        ['Length', `${d.length} aa`],
        ['V gene', d.v_gene],
        ['CDR1', d.cdr1 ? `${d.cdr1} (${d.cdr1_len} aa)` : null],
        ['CDR2', d.cdr2 ? `${d.cdr2} (${d.cdr2_len} aa)` : null],
        ['CDR3', d.cdr3 ? `${d.cdr3} (${d.cdr3_len} aa)` : null],
        ['CDR3 cluster', d.cdr3_cluster === null ? null
          : `${d.cdr3_cluster} (90 percent identity, within this target)`],
        ['Canonical Cys23-Cys104', d.canonical_disulfide ? 'present' : 'absent'],
        ['Additional disulfide', d.extra_disulfide ? 'present' : 'none detected'],
        ['CDR histidines', `${d.his_cdr_total} total · `
          + `CDR1 ${d.his_cdr1} · CDR2 ${d.his_cdr2} · CDR3 ${d.his_cdr3}`,
          'His enrichment was an intentional pH-switch design for EGFR and TNFa '
          + 'only. It is a design criterion, not a quality metric.'],
      ]);
  }

  function renderNeighboursPanel(d) {
    // The superposition leads: it is what justifies having a 3D viewer rather
    // than decorating with one (DESIGN 3.5). The model was moved onto the
    // experimental chain's frame by TM-align at build time, so this panel
    // loads two files and does no geometry of its own.
    const parts = [];
    if (d.nearest_pdb) {
      const canSuperpose = Boolean(d.superposed_url);
      // Placeholder outside #nviewer for the same reason as the Structure
      // panel: Mol* owns its mount point's children.
      parts.push(`<div class="viewerwrap">`
        + `<div id="nviewer" style="height:236px;background:var(--viewer);`
        + `border:1px solid var(--border);border-radius:4px;`
        + `position:relative;overflow:hidden"></div>`
        + `<p class="placeholder" style="position:absolute;inset:0;display:grid;`
        + `place-items:center;color:#8A93A5;font-size:11px;text-align:center;`
        + `padding:0 20px">${canSuperpose
          ? 'Loading the superposition.'
          : 'No superposed model for this candidate.'}</p>`
        + `</div><p style="font-size:10px;color:var(--ink-3);margin:4px 0 2px;`
        + `display:flex;justify-content:space-between;gap:10px"><span>model over `
        + `${W.escape(d.nearest_pdb)}`
        + `${d.nearest_chain ? ' chain ' + W.escape(d.nearest_chain) : ''}</span>`
        + `<span>${d.nearest_tm ? 'TM ' + d.nearest_tm.toFixed(2) : ''}`
        + `${d.nearest_rmsd ? ' &middot; RMSD ' + d.nearest_rmsd.toFixed(1) + ' A' : ''}`
        + `</span></p>`
        + `<p style="font-size:10px;color:var(--ink-3);margin:0 0 10px;`
        + `display:flex;gap:12px;align-items:center">`
        + cdrLegend()
        + `<span class="cdrlegend"><span class="cdrkey">`
        + `<i style="background:${REFERENCE_COLOUR}"></i>`
        + `${W.escape(d.nearest_pdb)}</span></span></p>`);
      parts.push(kvTable([
        ['Nearest experimental structure',
          d.nearest_pdb + (d.nearest_chain ? ` chain ${d.nearest_chain}` : '')],
        ['TM-score', d.nearest_tm,
          'Normalised by the model. Near 0.9 for almost every VHH, because '
          + 'every VHH has the same fold.'],
        ['RMSD over the alignment', d.nearest_rmsd === null
          || d.nearest_rmsd === undefined ? null : `${d.nearest_rmsd.toFixed(2)} A`],
        ['Framework identity', d.nearest_identity_fw === null
          || d.nearest_identity_fw === undefined ? null
          : `${(d.nearest_identity_fw * 100).toFixed(1)} percent`,
          'Framework identity is high for almost every VHH and tells you little.'],
        ['CDR identity', d.nearest_identity_cdr === null
          || d.nearest_identity_cdr === undefined ? null
          : `${(d.nearest_identity_cdr * 100).toFixed(1)} percent`,
          'This is the number that carries information.'],
        ['Experimental badge', d.has_experimental_badge ? 'yes' : 'no',
          'A positive signal only. Its absence says nothing about expressibility, '
          + 'and nothing in the scoring reads this column.'],
      ]));
      parts.push('<p style="font-size:11px;color:var(--ink-3);margin:8px 0 0">'
        + (d.reference_url
          ? 'Both structures are served from here: the neighbour is the single '
            + 'matching chain, extracted from its PDB entry, and the predicted '
            + 'model was superposed onto it at build time rather than in the '
            + 'browser.'
          : 'The experimental entry is fetched from the PDB by your browser. The '
            + 'predicted model was superposed onto it at build time.')
        + '</p>');
    } else {
      parts.push('<p style="font-size:12px;color:var(--ink-3)">The nearest '
        + 'experimental neighbour is found in the structural tier, which has not '
        + 'run for this candidate yet.</p>');
    }
    parts.push('<h3 style="font-size:11px;text-transform:uppercase;'
      + 'letter-spacing:.07em;color:var(--ink-3);margin:16px 0 6px">'
      + 'Upstream comparison columns</h3>');
    parts.push(kvTable([
      ['Upstream nearest sequence', d.upstream_global_nearest],
      ['Upstream global CDR identity', d.upstream_global_cdr_identity],
      ['Upstream target-matched CDR identity', d.upstream_target_cdr_identity,
        'Only EGFR and TNFa have known nanobodies in SAbDab-nano, so this is '
        + 'null for the other three targets rather than fabricated.'],
      ['Upstream composite score', d.upstream_score === null
        || d.upstream_score === undefined ? null
        : `${d.upstream_score.toFixed(4)} (not comparable between targets)`,
        'EGFR and TNFa blend four terms including s_his, the other three blend '
        + 'three, so the ranges are not on one scale. Kept for comparison only.'],
    ]));
    el('panel-neighbours').innerHTML = parts.join('');
  }

  function selectPanel(name) {
    const previous = W.state.panel;
    W.state.panel = name;
    document.querySelectorAll('.subtabs button').forEach((btn) => {
      btn.setAttribute('aria-selected', String(btn.dataset.panel === name));
    });
    ['structure', 'metrics', 'sequence', 'neighbours'].forEach((panel) => {
      el(`panel-${panel}`).hidden = panel !== name;
    });
    // Coming back to Structure after the viewer was skipped for a hidden
    // panel: load it now that it is visible.
    if (name === 'structure' && previous !== 'structure' && currentDetail
        && currentDetail.structure_url && window.WViewer) {
      const container = el('viewer');
      if (container && container.dataset.pending) {
        loadStructureView(currentDetail);
      }
    }
    if (name === 'neighbours' && currentDetail && currentDetail.nearest_pdb
        && window.WViewer) {
      const node = el('nviewer');
      if (node && !node.dataset.loaded && currentDetail.superposed_url) {
        node.dataset.loaded = '1';
        // The superposed copy, not the raw model: the raw one sits hundreds of
        // angstroms away in its own frame and would not be a superposition.
        // The neighbour is the single extracted chain served from here, not
        // the whole PDB entry from RCSB: a cryo-EM entry is a dozen chains and
        // the nanobody gets lost among them.
        showSuperposition(node, currentDetail);
      }
    }
  }

  /* The superposition: the predicted model over its nearest experimental
     neighbour, with the same CDR colours the Structure panel uses, so the
     reader can see at once whether the loops agree or only the framework
     does. The framework agreeing is the uninformative case, and it is what
     almost every pair of VHHs does.

     The superposed copy is IMGT-numbered in chain H like the model it came
     from. The experimental chain is drawn in one pale colour: it is the
     reference, not the subject. When the entry is not served from here the
     browser fetches it from the PDB, and a whole cryo-EM entry cannot be
     described in MVS without knowing its chains, so that case keeps the plain
     two-structure load. */
  async function showSuperposition(node, d) {
    const reference = d.reference_url;
    if (!reference) {
      const entries = [{ url: d.superposed_url, format: 'pdb' }];
      if (d.nearest_pdb) entries.push({ pdbId: d.nearest_pdb });
      await window.WViewer.show(node, entries);
      return;
    }
    const result = await window.WViewer.showMvs(node, (builder) => {
      builder.canvas({ background_color: VIEWER_BG });
      const experimental = builder
        .download({ url: reference })
        .parse({ format: 'pdb' })
        .modelStructure();
      experimental.component({ selector: 'polymer' })
        .representation({ type: 'cartoon' })
        .color({ color: REFERENCE_COLOUR });
      const model = builder
        .download({ url: d.superposed_url })
        .parse({ format: 'pdb' })
        .modelStructure();
      model.component({
        selector: frameworkRanges(MODEL_CHAIN, IMGT_CDR, IMGT_LAST) })
        .representation({ type: 'cartoon' })
        .color({ color: NANOBODY_COLOUR });
      for (const key of ['cdr1', 'cdr2', 'cdr3']) {
        model.component({ selector: cdrSelector(MODEL_CHAIN, IMGT_CDR[key]) })
          .representation({ type: 'cartoon' })
          .color({ color: CDR_COLOURS[key] });
      }
    });
    if (!result || !result.ok) {
      await window.WViewer.show(node, [
        { url: d.superposed_url, format: 'pdb' },
        { url: reference, format: 'pdb' },
      ]);
    }
  }

  function closeDrawer() {
    setDrawerOpen(false);
    W.state.id = null;
    W.writeUrl(true);
    // Dispose the Mol* instance: browsers cap WebGL contexts near 16 and the
    // viewer dies silently without this (DESIGN 3.6).
    if (window.WViewer) window.WViewer.dispose();
    const nviewer = el('nviewer');
    if (nviewer) delete nviewer.dataset.loaded;
    document.querySelectorAll('.tabulator-row').forEach(
      (r) => r.classList.remove('whittle-selected'));
  }

  /* ------------------------------------------------------------------ wiring */

  async function selectTarget(name, push) {
    W.state.target = name;
    W.state.brush = {};
    // Every percentile, median and extent on this tab is computed within the
    // target, so none of them survives a switch.
    pctCache = {};
    planeBox = null;
    tmDist = null;
    renderTargetBar(name);
    bundle = await W.loadJson(`/static/data/candidates_${name}.json`);
    if (!bundle) {
      el('caption').textContent = `No data bundle for ${name} yet.`;
      return;
    }
    el('targetcount').textContent = `${W.fmt(bundle.n)} candidates · `
      + `${bundle.columns_with_data.length} of ${REG.total} columns populated`;
    buildChooser();
    buildTable();
    refresh();
    if (push) W.writeUrl(true);
  }

  async function init() {
    W.parseUrl();
    axesMeta = await W.loadJson('/static/data/axes.json');
    const known = TARGETS.map((t) => t.target);
    if (!W.state.target || !known.includes(W.state.target)) {
      W.state.target = known[0] || null;
    }
    if (!W.state.target) {
      el('caption').textContent = 'No targets configured.';
      return;
    }

    el('targetset').addEventListener('click', (event) => {
      const btn = event.target.closest('button[data-target]');
      if (!btn) return;
      W.state.id = null;
      setDrawerOpen(false);
      if (window.WViewer) window.WViewer.dispose();
      selectTarget(btn.dataset.target, true);
    });

    el('search').value = W.state.search || '';
    let searchTimer = null;
    el('search').addEventListener('input', (event) => {
      clearTimeout(searchTimer);
      searchTimer = setTimeout(() => {
        W.state.search = event.target.value;
        W.writeUrl(true);
        refresh();
      }, 150);
    });

    const pareto = el('paretotoggle');
    pareto.setAttribute('aria-pressed', String(W.state.pareto));
    pareto.addEventListener('click', () => {
      W.state.pareto = !W.state.pareto;
      pareto.setAttribute('aria-pressed', String(W.state.pareto));
      W.writeUrl(true);
      refresh();
    });

    const colbtn = el('columnsbtn');
    const colpop = el('columnpop');
    colbtn.addEventListener('click', () => {
      const open = colpop.hidden;
      colpop.hidden = !open;
      colbtn.setAttribute('aria-expanded', String(open));
      // With a popover open the table matters more than the chart (DESIGN 3.1).
      el('explorercard').classList.toggle('collapsed', open);
    });
    el('coldone').addEventListener('click', () => {
      colpop.hidden = true;
      colbtn.setAttribute('aria-expanded', 'false');
      el('explorercard').classList.remove('collapsed');
    });
    el('colreset').addEventListener('click', () => {
      W.state.cols = null;
      W.writeUrl(true);
      buildChooser();
      buildTable();
      refresh();
    });
    el('colsave').addEventListener('click', async (event) => {
      const btn = event.currentTarget;
      try {
        await navigator.clipboard.writeText(window.location.href);
        btn.textContent = 'Link copied';
      } catch (err) {
        btn.textContent = 'Copy the address bar';
      }
      setTimeout(() => { btn.textContent = 'Save as view'; }, 2200);
    });
    el('colfind').addEventListener('input', (event) => {
      const query = event.target.value.trim().toLowerCase();
      el('colbody').querySelectorAll('label[data-key]').forEach((label) => {
        const text = label.textContent.toLowerCase();
        label.style.display = !query || text.includes(query) ? '' : 'none';
      });
      el('colbody').querySelectorAll('.colgroup').forEach((group) => {
        const visible = [...group.querySelectorAll('label')]
          .some((l) => l.style.display !== 'none');
        group.style.display = visible ? '' : 'none';
      });
    });

    el('resetbrush').addEventListener('click', () => {
      W.state.brush = {};
      W.state.gates = {};
      planeBox = null;
      W.writeUrl(true);
      refresh();
    });
    el('expand').addEventListener('click', () => {
      el('explorercard').classList.remove('collapsed');
    });

    /* ------------------------------------------------------- hero controls */

    const sortSelect = el('laddersort');
    LADDER_SORTS.forEach(([field, label]) => {
      const option = document.createElement('option');
      option.value = field;
      option.textContent = label;
      if (field === (W.state.ladder || 'desirability')) option.selected = true;
      sortSelect.appendChild(option);
    });
    sortSelect.addEventListener('change', () => {
      W.state.ladder = sortSelect.value;
      W.writeUrl(true);
      drawLadder(activeRows());
    });
    el('laddern').addEventListener('change', () => drawLadder(activeRows()));

    [['planex', 'axis_specificity'], ['planey', 'axis_thermo']].forEach(
      ([id, preset]) => {
        const select = el(id);
        PLANE_FIELDS.forEach(([field, label]) => {
          const option = document.createElement('option');
          option.value = field;
          option.textContent = label;
          if (field === preset) option.selected = true;
          select.appendChild(option);
        });
        select.addEventListener('change', () => {
          planeBox = null;
          commitPlaneBox();
        });
      });
    bindPlane();

    el('heroview').addEventListener('click', (event) => {
      const btn = event.target.closest('button[data-view]');
      if (!btn) return;
      W.state.hero = btn.dataset.view;
      W.writeUrl(true);
      drawHero(activeRows());
    });
    el('gatesolve').addEventListener('click', solveForShortlist);
    el('gatereset').addEventListener('click', () => {
      W.state.gates = {};
      W.writeUrl(true);
      refresh();
    });

    let heroResize = null;
    window.addEventListener('resize', () => {
      clearTimeout(heroResize);
      heroResize = setTimeout(() => {
        if (heroView() === 'plane') drawPlane(activeRows());
        if (currentDetail) drawAxisStrips(currentDetail);
      }, 150);
    });

    const vmode = el('vmode');
    if (vmode) {
      vmode.addEventListener('click', (event) => {
        const btn = event.target.closest('button[data-mode]');
        if (!btn || btn.disabled || !currentDetail) return;
        if (btn.dataset.mode === viewerMode) return;
        viewerMode = btn.dataset.mode;
        syncViewerMode(currentDetail);
        loadStructureView(currentDetail);
      });
    }
    el('drawerclose').addEventListener('click', closeDrawer);
    document.querySelectorAll('.subtabs button').forEach((btn) => {
      btn.addEventListener('click', () => selectPanel(btn.dataset.panel));
    });
    document.addEventListener('keydown', (event) => {
      if (event.key === 'Escape' && drawerIsOpen()) closeDrawer();
    });
    /* A shared deep link pasted into an open tab, and the back button, both
       arrive here. An earlier version switched target and stopped, so
       #/t/TNFa/nb/TNFa-0002 arriving while EGFR was selected changed the
       target and silently left the drawer shut: the half of the link that
       names a candidate was dropped. Handle both parts, in order, every time.
       hashchange is listened for as well as popstate because a hash assignment
       does not always produce a popstate. */
    const applyUrl = async () => {
      const previousTarget = W.state.target;
      const previousId = W.state.id;
      W.parseUrl();
      if (W.state.target && W.state.target !== previousTarget) {
        await selectTarget(W.state.target, false);
      }
      if (W.state.id) {
        if (W.state.id !== previousId || !drawerIsOpen()) {
          await openDrawer(W.state.id);
        }
      } else if (drawerIsOpen()) {
        closeDrawer();
      }
    };
    let applying = false;
    const onUrlChange = async () => {
      if (applying) return;   // popstate and hashchange can both fire
      applying = true;
      try { await applyUrl(); } finally { applying = false; }
    };
    window.addEventListener('popstate', onUrlChange);
    window.addEventListener('hashchange', onUrlChange);

    await selectTarget(W.state.target, true);
    if (W.state.id) openDrawer(W.state.id);
  }

  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', init);
  } else {
    init();
  }
})();
