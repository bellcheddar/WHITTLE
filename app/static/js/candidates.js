/* Candidates tab: target selector, parcoords explorer, Tabulator table, column
   chooser, docked drawer. Spec 8.1 and DESIGN.md Section 3.

   Plotly is loaded on demand for the explorer, Mol* on first drawer open. Neither
   is in the initial bundle (DESIGN 3.6). */
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
      return (cell) => {
        const value = cell.getValue();
        if (value === null || value === undefined) {
          return '<span style="color:var(--ink-3)">pending</span>';
        }
        return W.fmt(value, col.decimals ?? 2);
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
    let rows = bundle.rows;
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
    drawExplorer(rows);
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

  /* ------------------------------------------------- 3.1 the parcoords explorer */

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

  const EXPLORER_AXES = [...W.AXIS_KEYS, 'axis_nativeness'];

  async function drawExplorer(rows) {
    const node = el('explorer');
    if (!node || !bundle) return;
    const usable = EXPLORER_AXES.filter(
      (a) => bundle.rows.some((r) => r[a] !== null && r[a] !== undefined));
    if (!usable.length) {
      node.innerHTML = '<p style="font-size:12px;color:var(--ink-3);margin:0;'
        + 'padding:28px 0;text-align:center">No axis has been computed yet. The '
        + 'explorer appears once the first scoring tier completes.</p>';
      el('explorerhint').textContent = 'no axes computed yet';
      return;
    }
    let Plotly;
    try { Plotly = await ensurePlotly(); } catch (err) {
      node.innerHTML = '<p style="font-size:12px;color:var(--ink-3);margin:0">'
        + 'The explorer chart could not load. The table below is unaffected.</p>';
      return;
    }

    const front = rows.filter((r) => r.on_front);
    const dominated = rows.filter((r) => !r.on_front);
    el('explorerhint').textContent = `drag an axis to brush · `
      + `${W.fmt(front.length)} candidates on the Pareto front`;

    const dims = (subset) => usable.map((axis) => ({
      label: W.AXIS_LABELS[axis],
      range: [0, 100],
      values: subset.map((r) => (r[axis] === null || r[axis] === undefined
        ? null : r[axis])),
      constraintrange: W.state.brush[axis] || undefined,
    }));

    // Two traces, with identical domains so the front is drawn exactly OVER the
    // dominated set rather than Plotly stacking them into separate bands.
    //
    // DESIGN 3.1 asks for colour AND line weight. plotly.js parcoords has no
    // per-line width attribute (verified against Plotly.PlotSchema: line carries
    // colour attributes only), and the spec mandates parcoords, so the second
    // channel here is lightness and opacity rather than weight: a very light
    // #C3CBD9 at 0.5 under a dark #3457A6 at 0.92. That step survives greyscale
    // and colour-vision deficiency, which is what the rule is for. The legend
    // swatches below the header do carry the weight difference, and the Pareto
    // toggle shows the front on its own. See DECISIONS.md.
    const FULL_DOMAIN = { x: [0, 1], y: [0, 1] };
    const traces = [];
    if (dominated.length) {
      traces.push({
        type: 'parcoords',
        domain: FULL_DOMAIN,
        line: { color: W.PALETTE.dominated },
        opacity: 0.5,
        dimensions: dims(dominated),
        labelside: 'top',
        labelfont: { family: 'Inter', size: 11, color: W.PALETTE.ink3 },
        tickfont: { family: 'Inter', size: 10, color: W.PALETTE.ink3 },
        rangefont: { family: 'Inter', size: 9, color: W.PALETTE.ink3 },
      });
    }
    if (front.length) {
      traces.push({
        type: 'parcoords',
        domain: FULL_DOMAIN,
        line: { color: W.PALETTE.accent },
        opacity: 0.92,
        dimensions: dims(front),
        labelside: 'top',
        labelfont: { family: 'Inter', size: 11, color: W.PALETTE.ink },
        tickfont: { family: 'Inter', size: 10, color: W.PALETTE.ink3 },
        rangefont: { family: 'Inter', size: 9, color: W.PALETTE.ink3 },
      });
    }

    await Plotly.react(node, traces, {
      height: 212,
      // t: 48 so the axis labels clear the card edge. DESIGN 3.1 asks for about
      // 196px of chart; the extra is label gutter, not plot area.
      margin: { l: 58, r: 40, t: 48, b: 22 },
      paper_bgcolor: 'rgba(0,0,0,0)',
      plot_bgcolor: 'rgba(0,0,0,0)',
      font: { family: 'Inter', size: 11, color: W.PALETTE.ink3 },
    }, { displayModeBar: false, responsive: true });

    node.removeAllListeners?.('plotly_restyle');
    node.on('plotly_restyle', (event) => {
      // Brushing filters the table beneath (DESIGN 3.1). Read the constraint
      // ranges back off whichever trace Plotly reports and normalise them.
      let changed = false;
      for (const patch of event) {
        if (!patch || typeof patch !== 'object') continue;
        for (const [key, value] of Object.entries(patch)) {
          const match = key.match(/^dimensions\[(\d+)\]\.constraintrange$/);
          if (!match) continue;
          const axis = usable[Number(match[1])];
          if (!axis) continue;
          if (!value || (Array.isArray(value) && !value.length)) {
            delete W.state.brush[axis];
          } else {
            const range = Array.isArray(value[0]) ? value[0] : value;
            W.state.brush[axis] = [Math.min(...range), Math.max(...range)];
          }
          changed = true;
        }
      }
      if (changed) { W.writeUrl(true); refreshAfterBrush(); }
    });
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
    const flag = verdict === 'vacuous'
      ? '<div class="tflag" title="This axis does not separate real candidates '
        + 'from CDR-shuffled controls, so its ordering carries no evidence.">'
        + 'non-discriminating</div>'
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

  async function openDrawer(id) {
    W.state.id = id;
    W.writeUrl(true);
    el('drawer').hidden = false;
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
      await window.WViewer.show(container,
        [{ url, format: 'pdb', binary: false }]);
      const panel = el('interactions');
      if (panel) panel.hidden = true;
    }
    writeViewerCaption(d);
  }

  /* The complex, with the residues PLIP found in contact picked out.

     Mol* 4.9's MolViewSpec has no primitives node, so these cannot be drawn
     as dashed lines between the two atoms. Colouring the contacting residues
     is more legible at this size in any case, and the list below carries the
     atom names and distances. */
  async function showComplex(d, profile) {
    const container = el('viewer');
    const byType = highlightResidues(profile && profile.interactions);
    await window.WViewer.showMvs(container, (builder) => {
      // Loading an MVS state resets the canvas, so the dark background has
      // to be part of the state rather than set once on the plugin.
      builder.canvas({ background_color: '#0D1117' });
      const structure = builder
        .download({ url: d.complex_url })
        .parse({ format: 'pdb' })
        .modelStructure();
      // Disjoint components, one per chain, rather than a `polymer` cartoon
      // with the nanobody drawn over it: two cartoons on the same geometry
      // z-fight and the upper one does not reliably win, which left the whole
      // complex a single grey.
      structure.component({ selector: ANTIGEN_CHAINS.map(
        (c) => ({ auth_asym_id: c })) })
        .representation({ type: 'cartoon' })
        .color({ color: ANTIGEN_COLOUR });
      structure.component({ selector: [{ auth_asym_id: 'N' }] })
        .representation({ type: 'cartoon' })
        .color({ color: NANOBODY_COLOUR });
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
  const ANTIGEN_COLOUR = '#7E8796';
  const NANOBODY_COLOUR = '#2F8F7E';

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
      el('vcapright').textContent = 'chain N = nanobody';
      return;
    }
    el('vcapleft').textContent = d.nearest_pdb
      ? `model vs ${d.nearest_pdb}${d.nearest_tm ? ` · TM ${d.nearest_tm.toFixed(2)}` : ''}`
      : '';
    right.textContent = d.cdr1 ? 'CDR1 / 2 / 3' : '';
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
        + `border-radius:4px;position:relative;overflow:hidden"></div>`
        + `<p class="placeholder" style="position:absolute;inset:0;display:grid;`
        + `place-items:center;color:#8A93A5;font-size:11px;text-align:center;`
        + `padding:0 20px">${canSuperpose
          ? 'Loading the superposition.'
          : 'No superposed model for this candidate.'}</p>`
        + `</div><p style="font-size:10px;color:var(--ink-3);margin:4px 0 10px;`
        + `display:flex;justify-content:space-between"><span>model over `
        + `${W.escape(d.nearest_pdb)}`
        + `${d.nearest_chain ? ' chain ' + W.escape(d.nearest_chain) : ''}</span>`
        + `<span>${d.nearest_tm ? 'TM ' + d.nearest_tm.toFixed(2) : ''}`
        + `${d.nearest_rmsd ? ' &middot; RMSD ' + d.nearest_rmsd.toFixed(1) + ' A' : ''}`
        + `</span></p>`);
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
        const entries = [{ url: currentDetail.superposed_url, format: 'pdb' }];
        if (currentDetail.reference_url) {
          entries.push({ url: currentDetail.reference_url, format: 'pdb' });
        } else if (currentDetail.nearest_pdb) {
          entries.push({ pdbId: currentDetail.nearest_pdb });
        }
        window.WViewer.show(node, entries);
      }
    }
  }

  function closeDrawer() {
    el('drawer').hidden = true;
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
      el('drawer').hidden = true;
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
      W.writeUrl(true);
      refresh();
    });
    el('expand').addEventListener('click', () => {
      el('explorercard').classList.remove('collapsed');
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
      if (event.key === 'Escape' && !el('drawer').hidden) closeDrawer();
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
        if (W.state.id !== previousId || el('drawer').hidden) {
          await openDrawer(W.state.id);
        }
      } else if (!el('drawer').hidden) {
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
