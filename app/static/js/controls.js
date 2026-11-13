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

  const AXIS_PAIRS = [
    ['observed_tm', 'predicted_thermo', 'Thermostability', 'observed Tm, deg C',
      'predicted thermostability percentile'],
    ['observed_kd', 'predicted_affinity', 'Affinity', 'observed Kd, nM',
      'predicted affinity percentile'],
    ['observed_expression', 'predicted_expression', 'Expression',
      'observed yield, mg/L', 'predicted expression percentile'],
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

    for (const [obs, pred, label, xlabel, ylabel] of panels) {
      const points = known.filter((r) => r[obs] !== null && r[pred] !== null);
      const failed = points.filter((r) => r.assertion_passed === 0);
      const passed = points.filter((r) => r.assertion_passed !== 0);
      const r2 = rsquared(points.map((r) => [r[obs], r[pred]]));
      const xs = points.map((r) => r[obs]);
      const lo = Math.min(...xs); const hi = Math.max(...xs);

      await Plotly.react(el(`plot-${obs}`), [
        { x: [lo, hi], y: [0, 100], mode: 'lines', name: 'identity',
          line: { color: W.PALETTE.ink3, dash: 'dash', width: 1 }, hoverinfo: 'skip' },
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
          gridcolor: W.PALETTE.rule, zeroline: false, range: [0, 100] },
        showlegend: false,
      }, { displayModeBar: false, responsive: true });

      // The residual plot, per axis (spec 8.2).
      const resid = points.map((r) => {
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

      el(`r2-${obs}`).textContent = r2 === null
        ? `${points.length} controls, too few for an R squared.`
        : `R squared ${r2.toFixed(3)} over ${points.length} controls`
          + (failed.length ? `, ${failed.length} failing their assertion.` : '.');
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
      columns: [
        { title: 'Sequence', field: 'id', frozen: true, width: 160 },
        { title: 'Targets', field: 'target', cssClass: '' },
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
    el('promisccaption').textContent = `${multi.length} rows. Percentiles are `
      + `within each row's own target, so a sequence high in two columns is high `
      + `against two unrelated antigens.`;
  }

  async function init() {
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
