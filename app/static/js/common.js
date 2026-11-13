/* WHITTLE shared helpers: formatting, the percentile scale, URL view state.

   DESIGN 2: the percentile scale is accent at >= 75, warn at 45-74, poor below
   45, and lightness differs as well as hue so it survives greyscale and
   colour-vision deficiency. Nothing in this file encodes a distinction in hue
   alone. */
'use strict';

const W = {};

W.PALETTE = {
  accent: '#3457A6',
  dominated: '#C3CBD9',
  warn: '#B5781C',
  poorInk: '#8A3A28',
  ink: '#111827',
  ink3: '#5B6472',
  rule: '#EBEFF5',
  border: '#DDE2EB',
  goodBg: '#DCF0E4',
  goodInk: '#1E4D35',
};

/* ------------------------------------------------------------------ numbers */

W.fmt = (value, decimals) => {
  if (value === null || value === undefined || value === '') return '';
  if (typeof value !== 'number') return String(value);
  return value.toLocaleString('en-GB', {
    minimumFractionDigits: decimals ?? 0,
    maximumFractionDigits: decimals ?? 0,
  });
};

W.pctClass = (value) => {
  if (value === null || value === undefined) return 'na';
  if (value >= 75) return 'hi';
  if (value >= 45) return 'mid';
  return 'lo';
};

/* The four-bar axis sparkline (DESIGN 3.3): 7px wide, up to 24px tall, 3px gap,
   so the shape of a candidate reads without opening anything. */
W.AXIS_KEYS = ['axis_affinity', 'axis_specificity', 'axis_thermo', 'axis_expression'];
W.AXIS_LABELS = {
  axis_affinity: 'Affinity',
  axis_specificity: 'Specificity',
  axis_thermo: 'Thermostability',
  axis_expression: 'Expression',
  axis_nativeness: 'Nativeness',
};

W.sparkline = (row) => {
  const bars = W.AXIS_KEYS.map((key) => {
    const value = row[key];
    if (value === null || value === undefined) {
      return '<i class="na" title="pending"></i>';
    }
    const height = Math.max(3, Math.round((value / 100) * 24));
    const label = `${W.AXIS_LABELS[key]} ${Math.round(value)}th percentile within target`;
    return `<i class="${W.pctClass(value)}" style="height:${height}px" title="${label}"></i>`;
  }).join('');
  const described = W.AXIS_KEYS
    .map((k) => `${W.AXIS_LABELS[k]} ${row[k] === null || row[k] === undefined
      ? 'pending' : Math.round(row[k])}`)
    .join(', ');
  return `<span class="spark" role="img" aria-label="${described}">${bars}</span>`;
};

W.pctCell = (value) => {
  if (value === null || value === undefined) {
    return '<span class="pv" style="color:var(--ink-3)">pending</span>';
  }
  const cls = W.pctClass(value);
  return `<span class="pctcell"><span class="bar" aria-hidden="true">`
    + `<i class="${cls}" style="width:${Math.max(2, value)}%"></i></span>`
    + `<span class="pv">${Math.round(value)}</span></span>`;
};

/* ------------------------------------------------------------- URL view state

   Spec 8.1: the deep link is #/t/<target>/nb/<id>, and column selection, brush
   state and the Pareto toggle persist in the URL so a configured view is
   shareable and citable. The hash carries the route; the query carries the
   view. */

W.state = {
  target: null,
  id: null,
  cols: null,      // array of column keys, or null for the defaults
  pareto: false,
  brush: {},       // axis key -> [lo, hi]
  search: '',
  panel: 'structure',
};

W.parseUrl = () => {
  const hash = window.location.hash.replace(/^#/, '');
  const route = hash.match(/^\/t\/([A-Za-z0-9]+)(?:\/nb\/([A-Za-z0-9-]+))?/);
  if (route) {
    W.state.target = route[1];
    W.state.id = route[2] || null;
  }
  const params = new URLSearchParams(window.location.search);
  if (params.has('cols')) W.state.cols = params.get('cols').split(',').filter(Boolean);
  if (params.has('pareto')) W.state.pareto = params.get('pareto') === '1';
  if (params.has('q')) W.state.search = params.get('q');
  if (params.has('panel')) W.state.panel = params.get('panel');
  if (params.has('brush')) {
    try {
      W.state.brush = JSON.parse(decodeURIComponent(params.get('brush')));
    } catch (err) { W.state.brush = {}; }
  }
  return W.state;
};

W.writeUrl = (replace) => {
  const params = new URLSearchParams();
  if (W.state.cols) params.set('cols', W.state.cols.join(','));
  if (W.state.pareto) params.set('pareto', '1');
  if (W.state.search) params.set('q', W.state.search);
  if (W.state.id && W.state.panel !== 'structure') params.set('panel', W.state.panel);
  if (Object.keys(W.state.brush).length) {
    params.set('brush', encodeURIComponent(JSON.stringify(W.state.brush)));
  }
  const hash = W.state.id
    ? `#/t/${W.state.target}/nb/${W.state.id}`
    : (W.state.target ? `#/t/${W.state.target}` : '');
  const query = params.toString();
  const url = `${window.location.pathname}${query ? '?' + query : ''}${hash}`;
  if (replace) window.history.replaceState(null, '', url);
  else window.history.pushState(null, '', url);
};

/* --------------------------------------------------------------------- fetch */

W.cache = new Map();

/* The data bundles change on every rebuild and are fetched without a version,
   so stamp them with this build's timestamp. Belt and braces against the
   nginx Cache-Control rule: a proxy that ignores the header still cannot
   serve last week's table. */
W.stamp = (path) => {
  if (!path.startsWith('/static/data/')) return path;
  const built = (window.WHITTLE_META || {}).built_at || '';
  if (!built) return path;
  const sep = path.includes('?') ? '&' : '?';
  return `${path}${sep}v=${encodeURIComponent(built)}`;
};

W.loadJson = async (path) => {
  if (W.cache.has(path)) return W.cache.get(path);
  const promise = fetch(W.stamp(path), { credentials: 'same-origin' })
    .then((r) => (r.ok ? r.json() : null))
    .catch(() => null);
  W.cache.set(path, promise);
  return promise;
};

W.escape = (text) => String(text === null || text === undefined ? '' : text)
  .replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;')
  .replace(/"/g, '&quot;');

window.W = W;
