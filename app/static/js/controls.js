/* Control set tab. Spec 8.2, DESIGN 4.

   Same table component as Tab 1, different data file. The scatter panels are
   drawn only when there is something to draw: an empty panel is more honest
   than one full of placeholder points. */
'use strict';

(function () {
  const el = (id) => document.getElementById(id);
  const META = window.WHITTLE_META || {};
  let plotlyReady = null;

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

  function rsquared(pairs) {
    if (pairs.length < 3) return null;
    const n = pairs.length;
    const mx = pairs.reduce((s, p) => s + p[0], 0) / n;
    const my = pairs.reduce((s, p) => s + p[1], 0) / n;
    let sxy = 0; let sxx = 0; let syy = 0;
    for (const [x, y] of pairs) {
      sxy += (x - mx) * (y - my);
      sxx += (x - mx) ** 2;
      syy += (y - my) ** 2;
    }
    if (sxx === 0 || syy === 0) return null;
    const r = sxy / Math.sqrt(sxx * syy);
    return r * r;
  }

  /* [observed, predicted, title, x label, y label, sameUnits]
     The thermostability panel compares degrees against degrees, so its
     reference line is y = x and the residual is a real temperature error. The
     other two compare a measurement against a percentile, where no identity
     line exists: those panels draw a trend only and say so. */
  const AXIS_PAIRS = [
    ['observed_tm', 'predicted_tm', 'Thermostability', 'observed Tm, deg C',
      'predicted Tm, deg C', true],
    ['observed_kd', 'predicted_affinity', 'Affinity', 'observed Kd, nM',
      'predicted affinity percentile', false],
    ['observed_expression', 'predicted_expression', 'Expression',
      'observed yield, mg/L', 'predicted expression percentile', false],
  ];

  async function drawScatters(rows) {
    const known = rows.filter((r) => r.kind === 'known');
    const panels = AXIS_PAIRS.filter(([obs, pred]) =>
      known.some((r) => r[obs] !== null && r[pred] !== null));
    if (!panels.length) return;
    const Plotly = await ensurePlotly().catch(() => null);
    if (!Plotly) return;

    const host = el('scatters');
    host.innerHTML = panels.map(([obs]) =>
      `<div><div class="plot" id="plot-${obs}"></div>`
      + `<div class="plot" id="resid-${obs}" style="height:150px"></div>`
      + `<p class="r2" id="r2-${obs}"></p></div>`).join('');

    for (const [obs, pred, label, xlabel, ylabel, sameUnits] of panels) {
      const points = known.filter((r) => r[obs] !== null && r[pred] !== null);
      const failed = points.filter((r) => r.assertion_passed === 0);
      const passed = points.filter((r) => r.assertion_passed !== 0);
      const r2 = rsquared(points.map((r) => [r[obs], r[pred]]));
      const xs = points.map((r) => r[obs]);
      const ys = points.map((r) => r[pred]);
      const lo = Math.min(...xs); const hi = Math.max(...xs);
      // y = x only where both axes are the same quantity in the same units.
      const refLo = sameUnits ? Math.min(lo, Math.min(...ys)) : lo;
      const refHi = sameUnits ? Math.max(hi, Math.max(...ys)) : hi;

      await Plotly.react(el(`plot-${obs}`), [
        { x: [refLo, refHi], y: sameUnits ? [refLo, refHi] : [0, 100],
          mode: 'lines', name: sameUnits ? 'y = x' : 'range',
          line: { color: W.PALETTE.ink3, dash: 'dash', width: 1 },
          hoverinfo: 'skip' },
        { x: passed.map((r) => r[obs]), y: passed.map((r) => r[pred]),
          mode: 'markers', name: 'control',
          text: passed.map((r) => `${r.name} (${r.target})`),
          marker: { color: W.PALETTE.accent, size: 8 } },
        { x: failed.map((r) => r[obs]), y: failed.map((r) => r[pred]),
          mode: 'markers', name: 'assertion failed',
          text: failed.map((r) => `${r.name} (${r.target})`),
          marker: { color: W.PALETTE.poorInk, size: 9, symbol: 'x' } },
      ], {
        title: { text: label, font: { family: 'Inter', size: 12 } },
        height: 260, margin: { l: 48, r: 14, t: 34, b: 42 },
        paper_bgcolor: 'rgba(0,0,0,0)', plot_bgcolor: 'rgba(0,0,0,0)',
        font: { family: 'Inter', size: 10, color: W.PALETTE.ink3 },
        xaxis: { title: { text: xlabel, font: { size: 10 } },
          gridcolor: W.PALETTE.rule, zeroline: false },
        yaxis: { title: { text: ylabel, font: { size: 10 } },
          gridcolor: W.PALETTE.rule, zeroline: false,
          range: sameUnits ? undefined : [0, 100] },
        showlegend: false,
      }, { displayModeBar: false, responsive: true });

      // The residual plot, per axis (spec 8.2). Where the two axes share units
      // the residual is predicted minus observed, in those units. Where they do
      // not, it is the departure from the rank the measurement implies.
      const resid = points.map((r) => {
        if (sameUnits) return r[pred] - r[obs];
        const expected = 100 * (r[obs] - lo) / ((hi - lo) || 1);
        return r[pred] - expected;
      });
      await Plotly.react(el(`resid-${obs}`), [
        { x: xs, y: resid, mode: 'markers', name: 'residual',
          marker: { color: W.PALETTE.accent, size: 7 } },
      ], {
        height: 150, margin: { l: 48, r: 14, t: 8, b: 36 },
        paper_bgcolor: 'rgba(0,0,0,0)', plot_bgcolor: 'rgba(0,0,0,0)',
        font: { family: 'Inter', size: 10, color: W.PALETTE.ink3 },
        xaxis: { title: { text: xlabel, font: { size: 10 } },
          gridcolor: W.PALETTE.rule, zeroline: false },
        yaxis: { title: { text: 'residual', font: { size: 10 } },
          gridcolor: W.PALETTE.rule, zerolinecolor: W.PALETTE.border },
        showlegend: false,
      }, { displayModeBar: false, responsive: true });

      let caption = r2 === null
        ? `${points.length} controls, too few for an R squared.`
        : `R squared ${r2.toFixed(3)} over ${points.length} controls`
          + (failed.length ? `, ${failed.length} failing their assertion.` : '.');
      if (sameUnits && points.length) {
        const mae = points.reduce((s, r) => s + Math.abs(r[pred] - r[obs]), 0)
          / points.length;
        caption += ` Mean absolute error ${mae.toFixed(1)} deg C.`;
      }
      caption += ' Every control was run with its published measurement withheld.';
      el(`r2-${obs}`).textContent = caption;
    }
  }

  function promiscuityTable(rows) {
    const multi = rows.filter((r) => r.kind === 'multi_target');
    if (!multi.length) {
      el('promisctable').innerHTML = '<p style="font-size:12px;color:var(--ink-3)">'
        + 'The promiscuity panel fills once the axes are scored for every target '
        + 'the shared sequences appear under.</p>';
      return;
    }
    new Tabulator('#promisctable', {
      data: multi,
      layout: 'fitDataStretch',
      height: '320px',
      groupBy: (row) => String(row.id).split('-').slice(0, 2).join('-'),
      groupHeader: (value, count) => `${W.escape(value)} `
        + `<span style="color:var(--ink-3);font-weight:400">`
        + `&middot; the same sequence against `
        + `${count} ${count === 1 ? 'target' : 'targets'}</span>`,
      columns: [
        { title: 'Row', field: 'id', frozen: true, width: 190 },
        { title: 'Target', field: 'target', cssClass: '' },
        { title: 'Affinity', field: 'predicted_affinity', cssClass: 'num',
          formatter: (c) => W.pctCell(c.getValue()) },
        { title: 'Specificity', field: 'predicted_specificity', cssClass: 'num',
          formatter: (c) => W.pctCell(c.getValue()) },
        { title: 'Thermostability', field: 'predicted_thermo', cssClass: 'num',
          formatter: (c) => W.pctCell(c.getValue()) },
        { title: 'Expression', field: 'predicted_expression', cssClass: 'num',
          formatter: (c) => W.pctCell(c.getValue()) },
        { title: 'Reading', field: 'notes', cssClass: '', minWidth: 240 },
      ],
    });
    const sequences = new Set(multi.map(
      (r) => String(r.id).split('-').slice(0, 2).join('-')));
    el('promisccaption').textContent = `${multi.length} rows across `
      + `${sequences.size} sequences. Each row is that sequence's percentile `
      + `against its own target's candidates, so two rows of one sequence can `
      + `be compared directly: a sequence high in both is high against two `
      + `unrelated antigens.`;
  }

  /* The QC panel holds structured results as well as sentences. Dumping the
     JSON makes the most important page in the build unreadable, so each known
     shape gets a summary and anything else is shown compactly. */
  function renderQcValues() {
    document.querySelectorAll('#qctable tr[data-key]').forEach((row) => {
      const holder = row.querySelector('.qcjson');
      const out = row.querySelector('.qcrendered');
      if (!holder || !out) return;
      let value;
      try { value = JSON.parse(holder.textContent); } catch (err) { return; }
      const key = row.dataset.key;

      if (key === 'axis_correlations') {
        const lines = Object.entries(value).map(([target, info]) => {
          const m = info.matrix || {};
          const pairs = [];
          const axes = Object.keys(m);
          axes.forEach((a, i) => axes.slice(i + 1).forEach((b) => {
            const r = m[a] && m[a][b];
            if (typeof r === 'number') pairs.push([a, b, r]);
          }));
          if (!pairs.length) return `${target}: no axis pair computed yet`;
          pairs.sort((x, y) => Math.abs(y[2]) - Math.abs(x[2]));
          const [a, b, r] = pairs[0];
          const flagged = (info.flagged || []).length;
          return `<strong>${W.escape(target)}</strong>: strongest pair `
            + `${W.escape(a)} and ${W.escape(b)} at r = ${r.toFixed(3)}`
            + (flagged ? `, ${flagged} pair(s) above 0.8` : ', none above 0.8');
        });
        out.innerHTML = lines.join('<br>');
        return;
      }

      if (key === 'axis_availability') {
        const lines = Object.entries(value).map(([axis, info]) => {
          if (!info.available) {
            return `<strong>${W.escape(axis)}</strong>: pending`;
          }
          const parts = (info.components || []).map((c) => c.note
            + (c.sign < 0 ? ' (lower is better)' : ''));
          return `<strong>${W.escape(axis)}</strong>, tier ${info.tier}: `
            + W.escape(parts.join('; '));
        });
        out.innerHTML = lines.join('<br>');
        return;
      }

      if (key === 'vacuous_tests') {
        // One verdict per axis, pooled over all five targets, with the
        // per-target counts beside it. Fifty pairs on one target cannot
        // resolve these effects: at a true 60 per cent rate the expected
        // p-value at n = 50 is 0.20, so the per-target verdicts this used to
        // print were a coin toss. It published specificity as informative on
        // PRLR and non-discriminating on TNFa, which is one 60 per cent effect
        // labelled two ways. The counts stay visible so a reader sees how thin
        // the evidence is per target (DECISIONS 53).
        const WORD = {
          informative: 'informative',
          weak: 'weak but measurable',
          vacuous: 'NOT DISCRIMINATING',
          inverted: 'RUNS BACKWARDS',
          untestable: 'this control cannot test it',
          pending: 'not computed yet',
        };
        const lines = Object.entries(value).map(([axis, entry]) => {
          const pooled = entry.pooled || {};
          const verdict = pooled.verdict || entry.verdict || 'pending';
          const word = WORD[verdict] || verdict;
          let head = `<strong>${W.escape(axis)}</strong>: ${W.escape(word)}`;
          if (pooled.pooled_n) {
            const pct = Math.round((pooled.pooled_fraction || 0) * 100);
            head += ` (pooled ${pooled.pooled_wins}/${pooled.pooled_n} paired,`
              + ` ${pct}%, p = ${Number(pooled.pooled_p).toPrecision(2)})`;
          }
          const per = Object.entries(pooled.per_target || {}).sort()
            .map(([t2, count]) => `${W.escape(t2)} ${W.escape(count)}`)
            .join(', ');
          const note = pooled.finding ? W.escape(pooled.finding) : '';
          // Structural faults are still per target: no variance or a modal
          // collapse is a property of one target's distribution, so those
          // stay listed against the target that has them.
          const own = Object.entries(entry.targets || {})
            .filter(([, t2]) => (t2.findings || []).length)
            .map(([t2, block]) => `${W.escape(t2)}: `
              + W.escape((block.findings || []).join('; ')));
          return head
            + (per ? `<br><span class="qcnote">per target ${per}</span>` : '')
            + (note ? `<br><span class="qcnote">${note}</span>` : '')
            + (own.length
              ? `<br><span class="qcnote">${own.join('<br>')}</span>` : '');
        });
        out.innerHTML = lines.join('<br>');
        return;
      }

      if (key === 'pareto_independent_check') {
        const lines = Object.entries(value).map(([target, info]) =>
          `<strong>${W.escape(target)}</strong>: ${info.front_stored} on the `
          + `stored front, ${info.front_recomputed} on the independently `
          + `recomputed one, ${info.agree ? 'they agree' : 'THEY DISAGREE'}`);
        out.innerHTML = lines.join('<br>');
        return;
      }

      const text = JSON.stringify(value);
      out.textContent = text.length > 300 ? text.slice(0, 300) + '...' : text;
    });
  }

  async function init() {
    renderQcValues();
    const data = await W.loadJson('/static/data/controls.json');
    const rows = (data && data.rows) || [];
    if (!rows.length) {
      el('promisctable').innerHTML = '<p style="font-size:12px;color:var(--ink-3)">'
        + 'The control set is assembled in the co-folding phase and is not yet '
        + 'populated for this run.</p>';
      return;
    }
    promiscuityTable(rows);
    await drawScatters(rows);
    const scrambled = rows.filter((r) => r.kind === 'scrambled');
    if (scrambled.length) {
      const real = rows.filter((r) => r.kind === 'known')
        .map((r) => r.predicted_affinity).filter((v) => v !== null);
      const shuffled = scrambled.map((r) => r.predicted_affinity)
        .filter((v) => v !== null);
      if (real.length && shuffled.length) {
        const mr = real.reduce((a, b) => a + b, 0) / real.length;
        const ms = shuffled.reduce((a, b) => a + b, 0) / shuffled.length;
        const separates = mr - ms > 10;
        el('scrambledsummary').innerHTML = separates
          ? `Real candidates average the ${Math.round(mr)}th percentile on the `
            + `affinity axis against the ${Math.round(ms)}th for shuffled CDRs, so `
            + `the axis does separate them. It remains a co-folding confidence `
            + `prior, not an affinity measurement.`
          : `<strong>Finding, not a failure.</strong> Real candidates average the `
            + `${Math.round(mr)}th percentile on the affinity axis against the `
            + `${Math.round(ms)}th for shuffled CDRs. The axis does not separate `
            + `them, so it is marked non-discriminating wherever it appears and it `
            + `does not drive the Pareto front.`;
      }
    }
  }

  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', init);
  } else { init(); }
})();
