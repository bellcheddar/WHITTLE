// Partial Plotly bundle: scatter, bar and parcoords only (spec 8.1, DESIGN 3.6).
// parcoords is NOT in plotly.js-basic-dist-min, which is why this file exists.
// The full dist is 4.6 MB; this comes out around a quarter of that.
import Plotly from 'plotly.js/lib/core';
import scatter from 'plotly.js/lib/scatter';
import bar from 'plotly.js/lib/bar';
import parcoords from 'plotly.js/lib/parcoords';

Plotly.register([scatter, bar, parcoords]);
window.Plotly = Plotly;
export default Plotly;
