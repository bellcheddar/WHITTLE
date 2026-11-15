/* The alignment panel at the foot of the Candidates tab.

   Nothing here aligns anything. `pipeline/alignment.py` has already placed
   every residue on its IMGT key, and two residues share a column exactly when
   they share a key, so the grid is read off rather than computed. That is the
   point: a progressive aligner is free to slide a CDR against a framework to
   buy a cheaper gap, and then a column no longer means one structural position
   across the rows, which is the only thing this panel is for.

   The selection follows the table. Gates, the Pareto toggle, the search box
   and a brush on the trade-off plane all change which candidates are in
   contention, and the alignment follows that rather than a fixed 24, because
   the reader changing the selection is the reason they are looking.

   Identity, similarity and the consensus are computed here, over the rows that
   are actually shown, for the same reason. A consensus precomputed over all
   500 would describe a population the reader is not looking at. */
'use strict';

(function () {
  const el = (id) => document.getElementById(id);
  const META = window.WHITTLE_META || {};

  const BLOSUM_ALPHABET = 'ARNDCQEGHILKMFPSTWYV';
  const BLOSUM62 = [
    '  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0',  // A
    ' -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3',  // R
    ' -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3',  // N
    ' -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3',  // D
    '  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1',  // C
    ' -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2',  // Q
    ' -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2',  // E
    '  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3',  // G
    ' -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3',  // H
    ' -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3',  // I
    ' -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1',  // L
    ' -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2',  // K
    ' -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1',  // M
    ' -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1',  // F
    ' -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2',  // P
    '  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2',  // S
    '  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0',  // T
    ' -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3',  // W
    ' -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1',  // Y
    '  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4',  // V
  ].map((line) => line.trim().split(/\s+/).map(Number));
  const BLOSUM_AT = Object.fromEntries(
    [...BLOSUM_ALPHABET].map((aa, i) => [aa, i]));

  const blosum = (a, b) => {
    const i = BLOSUM_AT[a]; const j = BLOSUM_AT[b];
    return (i === undefined || j === undefined) ? null : BLOSUM62[i][j];
  };

  // Kyte and Doolittle 1982, the scale every hydropathy plot is drawn on.
  const KD = { I: 4.5, V: 4.2, L: 3.8, F: 2.8, C: 2.5, M: 1.9, A: 1.8, G: -0.4,
    T: -0.7, S: -0.8, W: -0.9, Y: -1.3, P: -1.6, H: -3.2, E: -3.5, Q: -3.5,
    D: -3.5, N: -3.5, K: -3.9, R: -4.5 };

  /* Clustal X residue classes, at the saturation the rest of the site works
     at. The classes are the standard ones, so a reader who knows the scheme
     reads it without a legend; the hues are muted so a 24 by 133 grid of them
     is still something you can look at. */
  const CLUSTAL = {
    hydrophobic: { set: 'AILMFWV', fill: '#D3E0F4', label: 'hydrophobic' },
    cysteine: { set: 'C', fill: '#F2CFE4', label: 'cysteine' },
    positive: { set: 'KR', fill: '#F7D5D0', label: 'positive, K and R' },
    negative: { set: 'ED', fill: '#EBD4EC', label: 'negative, E and D' },
    polar: { set: 'NQST', fill: '#D5E9DC', label: 'polar, N Q S T' },
    glycine: { set: 'G', fill: '#F8DEC6', label: 'glycine' },
    proline: { set: 'P', fill: '#F7EFC6', label: 'proline' },
    aromatic: { set: 'HY', fill: '#CDE6E9', label: 'aromatic, H and Y' },
  };
  const CLUSTAL_FILL = {};
  Object.values(CLUSTAL).forEach((group) => {
    [...group.set].forEach((aa) => { CLUSTAL_FILL[aa] = group.fill; });
  });

  const SSE_FILL = { H: '#EFD3E4', E: '#FAEAC4', C: '#F1F4F9', '.': '' };
  const SSE_WORD = { H: 'helix', E: 'strand', C: 'neither helix nor strand',
    '.': 'not assigned: the model carries a different residue at this position' };

  const REGION_FILL = { cdr1: '#E6EDFA', cdr2: '#DCF0E4', cdr3: '#FBF0DC',
    framework: '' };

  const SCHEMES = [
    ['clustal', 'Clustal X, residue class'],
    ['blosum', 'BLOSUM62 against the column consensus'],
    ['conserve', 'Conservation across the rows shown'],
    ['hydro', 'Hydrophobicity, Kyte and Doolittle'],
    ['charge', 'Charge at pH 7.4'],
    ['sse', 'Secondary structure, from each model'],
    ['region', 'Framework and CDR'],
    ['none', 'None, annotations only'],
  ];

  const SELECTIONS = [
    ['shortlist', 'The shortlist'],
    ['top50', 'The 50 that were co-folded'],
    ['filter', 'Whatever the table is showing'],
  ];
  /* A grid is something you read, and five hundred rows of it is not. When the
     table is showing more than this, the panel takes the best of them by rank
     and says so, rather than rendering sixty thousand cells nobody will look
     at. There is no second way to reach the Pareto front here: the table's own
     toggle narrows the selection and this follows it. */
  const MAX_ROWS = 60;
  let truncated = 0;

  let data = null;          // the alignment bundle for the active target
  let loading = null;       // in-flight fetch, so two events do not race
  let tableRows = [];       // ids the table currently shows, in its order
  let reference = null;     // the row every identity is measured against
  let scheme = 'clustal';
  let selection = 'shortlist';
  let variableOnly = false;
  let showControls = true;

  /* ------------------------------------------------------------ the grid */

  function regionOf(key) {
    const position = parseInt(key, 10);
    const spans = (data && data.cdr_spans) || {};
    for (const [name, [start, end]] of Object.entries(spans)) {
      if (position >= start && position <= end) return name;
    }
    return 'framework';
  }

  /** The rows to align, as bundle rows, in the order they are drawn. */
  function chosenRows() {
    if (!data) return [];
    const byId = Object.fromEntries(data.rows.map((r) => [r.id, r]));
    let rows;
    if (selection === 'filter') {
      rows = tableRows.map((id) => byId[id]).filter(Boolean);
    } else {
      const limit = selection === 'top50' ? (data.cofolded || 50)
        : (data.shortlist || 24);
      rows = data.rows.filter((r) => r.rank && r.rank <= limit);
    }
    truncated = Math.max(0, rows.length - MAX_ROWS);
    return truncated ? rows.slice(0, MAX_ROWS) : rows;
  }

  /** Residue letter per column index for one row, '' where the row has none. */
  function expand(row) {
    const columns = data.skeletons[row.sk] || [];
    const out = new Array(data.columns.length).fill('');
    const sse = new Array(data.columns.length).fill('');
    columns.forEach((columnIndex, i) => {
      out[columnIndex] = row.seq[i] || '';
      sse[columnIndex] = row.sse ? (row.sse[i] || '') : '';
    });
    const liabilities = new Array(data.columns.length).fill(null);
    (row.lia || []).forEach(([residue, kind]) => {
      const columnIndex = columns[residue];
      if (columnIndex !== undefined) liabilities[columnIndex] = kind;
    });
    return { row, aa: out, sse, lia: liabilities };
  }

  /** Per column: the commonest residue, how many rows carry it, and how many
      rows have any residue there at all. A gap is not a residue and never
      wins the column: a column that is mostly gap says so by its occupancy. */
  function consensusOf(expanded) {
    return data.columns.map((_key, c) => {
      const counts = {};
      let occupied = 0;
      expanded.forEach((e) => {
        const aa = e.aa[c];
        if (!aa) return;
        occupied += 1;
        counts[aa] = (counts[aa] || 0) + 1;
      });
      let best = ''; let n = 0;
      Object.entries(counts).forEach(([aa, count]) => {
        if (count > n) { best = aa; n = count; }
      });
      return { aa: best, n, occupied, rows: expanded.length,
               identity: occupied ? n / occupied : 0 };
    });
  }

  /* ------------------------------------------------------- similarity */

  /** Identity and BLOSUM62 similarity of one row against the reference, over
      the columns where both have a residue. Gaps are excluded from the
      denominator rather than counted as mismatches: a shorter CDR3 is a
      different loop length, not 6 substitutions. */
  function scoreAgainst(expanded, ref) {
    let shared = 0; let identical = 0; let similar = 0; let score = 0;
    let cdr3Shared = 0; let cdr3Identical = 0;
    data.columns.forEach((key, c) => {
      const a = ref.aa[c]; const b = expanded.aa[c];
      if (!a || !b) return;
      shared += 1;
      const s = blosum(a, b);
      if (a === b) identical += 1;
      if (s !== null) { score += s; if (s > 0) similar += 1; }
      if (regionOf(key) === 'cdr3') {
        cdr3Shared += 1;
        if (a === b) cdr3Identical += 1;
      }
    });
    return {
      shared,
      identity: shared ? (100 * identical) / shared : null,
      similarity: shared ? (100 * similar) / shared : null,
      bits: shared ? score / shared : null,
      cdr3: cdr3Shared ? (100 * cdr3Identical) / cdr3Shared : null,
    };
  }

  /* ---------------------------------------------------------- colouring */

  function mix(hex, amount) {
    // Toward white, so every ramp ends at the page background rather than at
    // a grey that reads as a value of its own.
    const n = parseInt(hex.slice(1), 16);
    const r = (n >> 16) & 255; const g = (n >> 8) & 255; const b = n & 255;
    const f = (v) => Math.round(255 - (255 - v) * Math.max(0, Math.min(1, amount)));
    return `rgb(${f(r)},${f(g)},${f(b)})`;
  }

  function fillFor(aa, c, expanded, consensus) {
    if (!aa) return '';
    if (scheme === 'none') return '';
    if (scheme === 'clustal') return CLUSTAL_FILL[aa] || '';
    if (scheme === 'region') return REGION_FILL[regionOf(data.columns[c])] || '';
    if (scheme === 'sse') return SSE_FILL[expanded.sse[c]] || '';
    if (scheme === 'conserve') {
      const col = consensus[c];
      if (!col.occupied) return '';
      return mix('#3457A6', (col.identity - 0.25) / 0.75);
    }
    if (scheme === 'charge') {
      if ('KR'.includes(aa)) return mix('#3457A6', 0.55);
      if ('DE'.includes(aa)) return mix('#8A3A28', 0.5);
      if (aa === 'H') return mix('#3457A6', 0.2);
      return '';
    }
    if (scheme === 'hydro') {
      const value = KD[aa];
      if (value === undefined) return '';
      return value >= 0 ? mix('#B5781C', value / 4.5)
                     : mix('#3457A6', -value / 4.5);
    }
    if (scheme === 'blosum') {
      const col = consensus[c];
      if (!col.aa) return '';
      const s = blosum(aa, col.aa);
      if (s === null) return '';
      return s > 0 ? mix('#1E4D35', Math.min(1, s / 6))
                   : (s < 0 ? mix('#8A3A28', Math.min(1, -s / 4)) : '');
    }
    return '';
  }


  /* ------------------------------------------------------------ rendering */

  function liabilityTitle(kindIndex, key) {
    const kind = (data.liability_kinds || [])[kindIndex];
    if (!kind) return '';
    return `${kind.replace(/_/g, ' ')} motif, IMGT ${key}`;
  }

  function headerCells(columnIndices, consensus) {
    const hallmarks = new Set((data.hallmark_positions || []).map(String));
    const disulfides = new Set((data.disulfide_positions || []).map(String));
    const ruler = []; const marks = []; const band = [];
    const cons = []; const bars = [];
    // The region band carries its own label, centred over the run of columns
    // it covers, so the reader does not have to count against the ruler.
    const runs = [];
    columnIndices.forEach((c, i) => {
      const region = regionOf(data.columns[c]);
      const last = runs[runs.length - 1];
      if (last && last.region === region) last.n += 1;
      else runs.push({ region, start: i, n: 1 });
    });
    const labelAt = {};
    runs.forEach((run) => {
      if (run.region === 'framework' || run.n < 4) return;
      labelAt[run.start + Math.floor(run.n / 2)] = run.region.toUpperCase();
    });

    columnIndices.forEach((c, i) => {
      const key = data.columns[c];
      const position = String(parseInt(key, 10));
      const region = regionOf(key);
      const bare = key === position;
      // A number every tenth position. The cell is 11px and a three-digit
      // number is not, so the text is allowed to overflow its own cell: the
      // next number is ten cells away, which is a hundred pixels of room.
      const show = bare && parseInt(key, 10) % 10 === 0;
      ruler.push(`<span class="alncell ruler">${show ? position : ''}</span>`);

      let mark = ''; let what = '';
      if (bare && hallmarks.has(position)) {
        mark = ' hallmark'; what = `IMGT ${position}, camelid hallmark position`;
      } else if (bare && disulfides.has(position)) {
        mark = ' disulf'; what = `IMGT ${position}, canonical disulfide`;
      }
      marks.push(`<span class="alncell mark${mark}"${what
        ? ` title="${W.escape(what)}"` : ''}>${mark ? '<i></i>' : ''}</span>`);

      band.push(`<span class="alncell band ${region}" title="${region === 'framework'
        ? 'framework' : region.toUpperCase()}, IMGT ${key}">`
        + `${labelAt[i] ? `<b>${labelAt[i]}</b>` : ''}</span>`);

      const col = consensus[c];
      /* Two different things share this cell, so they are drawn differently.
         The bar is agreement among the rows that have a residue here. The
         opacity is how many rows have one at all: a column where a single
         CDR3 is long enough to reach would otherwise show a full bar and
         read as perfectly conserved. */
      const occupancy = col.rows ? col.occupied / col.rows : 0;
      const faint = occupancy < 0.999
        ? ` style="opacity:${(0.25 + 0.75 * occupancy).toFixed(2)}"` : '';
      const where = `${col.occupied} of ${col.rows} rows have a residue at IMGT `
        + `${key}; ${col.n} of those carry ${col.aa || 'nothing'}`;
      cons.push(`<span class="alncell cons"${faint} title="${W.escape(where)}">`
        + `${col.aa || '-'}</span>`);
      bars.push(`<span class="alncell barcell" title="${W.escape(where)}">`
        + `<i style="height:${Math.round(col.identity * 100)}%`
        + `${occupancy < 0.999 ? `;opacity:${(0.3 + 0.7 * occupancy).toFixed(2)}` : ''}">`
        + '</i></span>');
    });
    return { ruler, marks, band, cons, bars };
  }

  function rowCells(expanded, columnIndices, consensus) {
    const cells = [];
    columnIndices.forEach((c) => {
      const aa = expanded.aa[c];
      const key = data.columns[c];
      if (!aa) {
        cells.push('<span class="alncell gap">-</span>');
        return;
      }
      const fill = fillFor(aa, c, expanded, consensus);
      const kind = expanded.lia[c];
      const classes = ['alncell'];
      if (kind !== null && kind !== undefined) classes.push('liab');
      // Dimmed only under the structural scheme: a residue the crystal did
      // not resolve is a fact about the structure, not about the sequence,
      // and dimming it under Clustal reads as a doubtful amino acid.
      if (scheme === 'sse' && expanded.sse[c] === '.') classes.push('nomodel');
      const title = `${expanded.row.id} ${aa}${key}`
        + (kind !== null && kind !== undefined ? ` | ${liabilityTitle(kind, key)}` : '')
        + (expanded.sse[c] ? ` | ${SSE_WORD[expanded.sse[c]] || ''}` : '');
      cells.push(`<span class="${classes.join(' ')}"`
        + `${fill ? ` style="background:${fill}"` : ''}`
        + ` title="${W.escape(title)}">${aa}</span>`);
    });
    return cells.join('');
  }

  /** The known controls for this target, as drawn rows. They are kept out of
      the consensus and the agreement bars on purpose: those describe the
      library, and a control is the thing being compared against it, not a
      member of it. */
  function chosenControls() {
    if (!data || !showControls) return [];
    return (data.controls || []).map(expand);
  }

  function controlChip(row) {
    if (row.role === 'binder') {
      return `<span class="alnrole binder" title="A published nanobody with an `
        + `experimental complex against ${W.escape(data.target)}. ${W.escape(row.note || '')}`
        + `">binder</span>`;
    }
    return '<span class="alnrole framework" title="No nanobody against this '
      + 'target has a published structure, so this is the deposited chain this '
      + "target's own shortlist is most often nearest to. It is a framework "
      + 'reference and does not bind this target.">framework</span>';
  }

  function render() {
    const host = el('alngrid');
    if (!host) return;
    if (!data) { host.innerHTML = '<p class="placeholder">Loading the alignment.</p>'; return; }
    const rows = chosenRows();
    if (!rows.length) {
      host.innerHTML = '<p class="placeholder">No candidate is in the current '
        + 'selection, so there is nothing to align.</p>';
      el('alnlegend').innerHTML = '';
      return;
    }
    const expanded = rows.map(expand);
    const controls = chosenControls();
    const consensus = consensusOf(expanded);

    // A control may reach a column no candidate does, so it keeps its column
    // even though it is not in the consensus the column is measured by.
    const controlled = new Set();
    controls.forEach((e) => e.aa.forEach((aa, c) => { if (aa) controlled.add(c); }));
    let columnIndices = data.columns.map((_k, i) => i)
      .filter((c) => consensus[c].occupied > 0 || controlled.has(c));
    if (variableOnly) {
      columnIndices = columnIndices.filter(
        (c) => consensus[c].identity < 1 || consensus[c].occupied < expanded.length
          || controls.some((e) => e.aa[c] && e.aa[c] !== consensus[c].aa));
    }

    const refRow = expanded.concat(controls).find((e) => e.row.id === reference)
      || expanded[0];
    reference = refRow.row.id;
    const scores = expanded.map((e) => scoreAgainst(e, refRow));
    const controlScores = controls.map((e) => scoreAgainst(e, refRow));
    const { ruler, marks, band, cons, bars } = headerCells(columnIndices, consensus);

    const label = (e, score, i) => {
      const isRef = e.row.id === reference;
      return `<button type="button" class="alnlabel${isRef ? ' isref' : ''}" `
        + `data-id="${W.escape(e.row.id)}" `
        + `title="${isRef ? 'The reference row: every identity on this panel is '
          + 'measured against it. Click another row to compare against that one '
          + 'instead.' : 'Compare every row against this one'}">`
        + `<span class="alnrank">${e.row.rank === null ? '-' : e.row.rank}</span>`
        + `<span class="alnid">${W.escape(e.row.id)}</span>`
        + `<span class="alnpct">${isRef ? 'ref'
          : `${W.fmt(score.identity, 1)}%`}</span>`
        + `<span class="alnpct sim">${isRef ? ''
          : `${W.fmt(score.similarity, 1)}%`}</span>`
        + `<span class="alnpct cdr">${score.cdr3 === null ? '-'
          : `${W.fmt(score.cdr3, 0)}%`}</span>`
        + '</button>';
    };

    const controlLabel = (e, score) => {
      const row = e.row;
      const isRef = row.id === reference;
      const bits = [];
      if (row.pdb) bits.push(`PDB ${row.pdb}`);
      if (row.observed_tm) bits.push(`observed Tm ${row.observed_tm} deg C`);
      bits.push(row.sse_source || 'no deposited chain on disk for this one, so '
        + 'it carries no secondary structure here');
      if (row.note) bits.push(row.note);
      return `<button type="button" class="alnlabel${isRef ? ' isref' : ''}" `
        + `data-id="${W.escape(row.id)}" title="${W.escape(bits.join(' | '))}">`
        + `<span class="alnrank">${controlChip(row)}</span>`
        + `<span class="alnid">${W.escape(row.pdb && row.chain
          ? `${row.pdb} ${row.chain}` : row.name)}</span>`
        + `<span class="alnpct">${isRef ? 'ref' : `${W.fmt(score.identity, 1)}%`}</span>`
        + `<span class="alnpct sim">${isRef ? '' : `${W.fmt(score.similarity, 1)}%`}</span>`
        + `<span class="alnpct cdr">${score.cdr3 === null ? '-'
          : `${W.fmt(score.cdr3, 0)}%`}</span>`
        + '</button>';
    };

    const head = '<div class="alnrow alnhead">'
      + '<span class="alnlabel alnhcell">'
      + '<span class="alnrank">#</span><span class="alnid">Candidate</span>'
      + '<span class="alnpct" title="Percent of the columns where both rows have a residue that carry the same one">ident</span>'
      + '<span class="alnpct sim" title="Percent of those columns whose BLOSUM62 score against the reference is positive">simil</span>'
      + '<span class="alnpct cdr" title="Identity over the CDR3 columns alone">CDR3</span>'
      + '</span>'
      + `<span class="alnseq">${ruler.join('')}</span></div>`
      + '<div class="alnrow alnhead alnmarks"><span class="alnlabel alnhcell">'
      + '<span class="alnlegendlabel" title="Blue: the four camelid hallmark '
      + 'positions, IMGT 42, 49, 50 and 52. Amber: the canonical disulfide, '
      + 'IMGT 23 and 104.">marked</span></span>'
      + `<span class="alnseq">${marks.join('')}</span></div>`
      + '<div class="alnrow alnhead alnband"><span class="alnlabel alnhcell">'
      + '<span class="alnlegendlabel">region</span></span>'
      + `<span class="alnseq">${band.join('')}</span></div>`;

    const bodyRows = expanded.map((e, i) =>
      `<div class="alnrow">${label(e, scores[i], i)}`
      + `<span class="alnseq">${rowCells(e, columnIndices, consensus)}</span></div>`
    ).join('');

    const controlRows = controls.length ? (
      '<div class="alnrow alnsep"><span class="alnlabel alnhcell">'
      + '<span class="alnlegendlabel" title="Published nanobodies on the same '
      + 'grid. They are deliberately not part of the consensus below: that '
      + 'describes the library, and these are what it is being compared '
      + 'against.">known controls</span></span>'
      + `<span class="alnseq">${columnIndices.map(() =>
        '<span class="alncell"></span>').join('')}</span></div>`
      + controls.map((e, i) =>
        `<div class="alnrow alnctrl">${controlLabel(e, controlScores[i])}`
        + `<span class="alnseq">${rowCells(e, columnIndices, consensus)}</span>`
        + '</div>').join('')
    ) : '';

    const foot = '<div class="alnrow alnfoot"><span class="alnlabel alnhcell">'
      + '<span class="alnlegendlabel" title="The commonest residue in each '
      + 'column, over the rows shown">consensus</span></span>'
      + `<span class="alnseq">${cons.join('')}</span></div>`
      + '<div class="alnrow alnfoot"><span class="alnlabel alnhcell">'
      + '<span class="alnlegendlabel" title="How many of the rows with a '
      + 'residue in this column carry the consensus residue">agreement</span>'
      + '</span>'
      + `<span class="alnseq">${bars.join('')}</span></div>`;

    host.innerHTML = head + bodyRows + controlRows + foot;
    host.querySelectorAll('.alnlabel[data-id]').forEach((node) => {
      node.addEventListener('click', () => {
        reference = node.dataset.id;
        render();
      });
    });
    renderLegend(rows, columnIndices, expanded);
  }

  function renderLegend(rows, columnIndices, expanded) {
    const node = el('alnlegend');
    if (!node) return;
    const parts = [];
    if (scheme === 'clustal') {
      parts.push(Object.values(CLUSTAL).map((g) =>
        `<span class="alnkey"><i style="background:${g.fill}"></i>${g.label}</span>`
      ).join(''));
    } else if (scheme === 'sse') {
      // The key names what is on screen. P-SEA resolves helix and strand, and
      // the immunoglobulin fold is a beta sandwich, so helix is absent from
      // every model in this build. Saying so is better than an unused swatch,
      // which reads as a colour the reader has failed to find.
      const seen = {};
      expanded.forEach((e) => [...(e.row.sse || '')].forEach((code) => {
        seen[code] = (seen[code] || 0) + 1;
      }));
      parts.push(['H', 'E', 'C'].filter((code) => seen[code]).map((code) =>
        `<span class="alnkey"><i style="background:${SSE_FILL[code]}"></i>`
        + `${SSE_WORD[code]}</span>`).join('')
        + (seen['.'] ? '<span class="alnkey"><i class="hatch"></i>not assigned, '
          + `${seen['.']} residues</span>` : '')
        + (seen.H ? '' : '<span class="alnkey">no helix in any of these models, '
          + 'which is what a beta sandwich looks like</span>'));
    } else if (scheme === 'region') {
      parts.push(['cdr1', 'cdr2', 'cdr3'].map((r) =>
        `<span class="alnkey"><i style="background:${REGION_FILL[r]}"></i>`
        + `${r.toUpperCase()}</span>`).join(''));
    } else if (scheme === 'blosum') {
      parts.push('<span class="alnkey"><i style="background:'
        + mix('#1E4D35', 0.8) + '"></i>favoured against the consensus</span>'
        + '<span class="alnkey"><i style="background:#FFFFFF;border:1px solid var(--border)"></i>neutral</span>'
        + '<span class="alnkey"><i style="background:' + mix('#8A3A28', 0.8)
        + '"></i>penalised</span>');
    } else if (scheme === 'hydro') {
      parts.push('<span class="alnkey"><i style="background:'
        + mix('#B5781C', 1) + '"></i>hydrophobic, +4.5</span>'
        + '<span class="alnkey"><i style="background:#FFFFFF;border:1px solid var(--border)"></i>0</span>'
        + '<span class="alnkey"><i style="background:' + mix('#3457A6', 1)
        + '"></i>hydrophilic, -4.5</span>');
    } else if (scheme === 'charge') {
      parts.push('<span class="alnkey"><i style="background:' + mix('#3457A6', 0.55)
        + '"></i>K and R</span><span class="alnkey"><i style="background:'
        + mix('#8A3A28', 0.5) + '"></i>D and E</span>'
        + '<span class="alnkey"><i style="background:' + mix('#3457A6', 0.2)
        + '"></i>H, partly charged at 7.4</span>');
    } else if (scheme === 'conserve') {
      parts.push('<span class="alnkey"><i style="background:' + mix('#3457A6', 1)
        + '"></i>every row agrees</span><span class="alnkey">'
        + '<i style="background:#FFFFFF;border:1px solid var(--border)"></i>'
        + 'a quarter or fewer agree</span>');
    }
    parts.push('<span class="alnkey"><i class="liabkey"></i>liability motif, '
      + 'from the sequence scan</span>'
      + '<span class="alnkey"><i class="hallkey"></i>camelid hallmark, '
      + 'IMGT 42, 49, 50, 52</span>'
      + '<span class="alnkey"><i class="disulfkey"></i>canonical disulfide, '
      + 'IMGT 23 and 104</span>');
    const withModel = expanded.filter((e) => e.row.sse).length;
    const shown = chosenControls();
    const binders = shown.filter((e) => e.row.role === 'binder').length;
    const controlNote = !shown.length ? ''
      : (binders
        ? `${binders} published nanobod${binders === 1 ? 'y' : 'ies'} with an `
          + `experimental complex against ${data.target} are aligned below the `
          + 'library, with secondary structure read off the deposited chain '
          + 'rather than a model.'
        : 'No nanobody against this target has a published structure, so the '
          + 'control row is the deposited chain this shortlist is most often '
          + 'nearest to. It is a framework reference and does not bind this '
          + 'target.');
    node.innerHTML = parts.join('')
      + `<span class="alnnote">${rows.length} rows`
      + (truncated ? `, the best by rank of the ${rows.length + truncated} the `
        + 'table is showing: narrow the selection with the gates, the search '
        + 'box or the Pareto toggle to see the rest' : '')
      + ', '
      + `${columnIndices.length} of ${data.columns.length} IMGT columns`
      + `${variableOnly ? ', the ones that vary' : ''}. `
      + `Columns are IMGT positions, not an alignment: no aligner is run, so a `
      + `column is one numbered position in every row. Secondary structure on `
      + `${withModel} of ${rows.length} rows, ${data.sse_source}. `
      + `${controlNote}</span>`;
  }

  /* -------------------------------------------------------------- wiring */

  async function load(target) {
    if (!target) return;
    if (data && data.target === target) return;
    data = null;
    render();
    loading = W.loadJson(`/static/data/alignment_${target}.json`);
    const loaded = await loading;
    // A target switch while a fetch was in flight: keep the newer one.
    if (!loaded || loaded.target !== W.state.target) return;
    data = loaded;
    reference = null;
    render();
  }

  function fillControls() {
    const schemeSelect = el('alnscheme');
    const selectionSelect = el('alnrows');
    if (schemeSelect && !schemeSelect.options.length) {
      schemeSelect.innerHTML = SCHEMES.map(([key, label]) =>
        `<option value="${key}">${W.escape(label)}</option>`).join('');
      schemeSelect.value = scheme;
      schemeSelect.addEventListener('change', () => {
        scheme = schemeSelect.value; render();
      });
    }
    if (selectionSelect && !selectionSelect.options.length) {
      selectionSelect.innerHTML = SELECTIONS.map(([key, label]) =>
        `<option value="${key}">${W.escape(label)
          .replace('The shortlist', `The shortlist, ${META.shortlist_per_target || 24}`)
          .replace('The 50 ', `The ${META.tier2_top_n_per_target || 50} `)}</option>`)
        .join('');
      selectionSelect.value = selection;
      selectionSelect.addEventListener('change', () => {
        selection = selectionSelect.value; render();
      });
    }
    const controlToggle = el('alncontrols');
    if (controlToggle && !controlToggle.dataset.wired) {
      controlToggle.dataset.wired = '1';
      showControls = controlToggle.checked;
      controlToggle.addEventListener('change', () => {
        showControls = controlToggle.checked;
        if (!showControls && reference
            && (data.controls || []).some((c) => c.id === reference)) {
          reference = null;           // the reference cannot be a hidden row
        }
        render();
      });
    }
    const varToggle = el('alnvariable');
    if (varToggle && !varToggle.dataset.wired) {
      varToggle.dataset.wired = '1';
      varToggle.addEventListener('change', () => {
        variableOnly = varToggle.checked; render();
      });
    }
    const copy = el('alncopy');
    if (copy && !copy.dataset.wired) {
      copy.dataset.wired = '1';
      copy.addEventListener('click', () => {
        const rows = chosenRows().map(expand);
        const text = rows.map((e) => {
          const seq = data.columns.map((_k, c) => e.aa[c] || '-').join('');
          return `>${e.row.id} rank=${e.row.rank === null ? 'none' : e.row.rank}\n${seq}`;
        }).join('\n');
        navigator.clipboard.writeText(text + '\n').then(() => {
          copy.textContent = 'Copied';
          setTimeout(() => { copy.textContent = 'Copy as FASTA'; }, 1600);
        }).catch(() => { copy.textContent = 'Could not copy'; });
      });
    }
  }

  function init() {
    if (!el('alngrid')) return;
    fillControls();
    document.addEventListener('whittle:selection', (event) => {
      const detail = event.detail || {};
      tableRows = detail.ids || [];
      if (detail.target && (!data || data.target !== detail.target)) {
        load(detail.target);
      } else if (selection === 'filter' || selection === 'front') {
        render();
      }
    });
    if (W.state.target) load(W.state.target);
  }

  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', init);
  } else {
    init();
  }
})();
