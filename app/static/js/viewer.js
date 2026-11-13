/* Mol* viewer. DESIGN 3.6, and these are bugs if ignored:

   - ONE plugin instance. Structures are reloaded into it; the instance is
     disposed when the drawer closes. Browsers cap WebGL contexts near 16, and
     without this the viewer dies silently after a dozen drawer opens and reads
     as a random failure.
   - Mol* (about 5 MB) is dynamic-imported on first drawer open, never in the
     initial bundle.
   - Structures are served .pdb.gz with gzip_static; the browser decompresses.

   Adapted from GOBSMACKED/app/static/js and chatPDB-web (see REUSE.md): the
   single-instance discipline and the dispose-on-close are the parts carried over.
*/
'use strict';

(function () {
  let plugin = null;         // the one instance
  let loading = null;        // the in-flight dynamic import
  let contextCount = 0;      // for the QC check: must stay bounded (spec 10 item 7)

  const CSS_ID = 'molstar-css';

  async function ensureLibrary() {
    if (window.molstar) return window.molstar;
    if (loading) return loading;
    loading = new Promise((resolve, reject) => {
      if (!document.getElementById(CSS_ID)) {
        const link = document.createElement('link');
        link.id = CSS_ID;
        link.rel = 'stylesheet';
        link.href = '/static/vendor/molstar.css';
        document.head.appendChild(link);
      }
      const script = document.createElement('script');
      script.src = '/static/vendor/molstar.js';
      script.async = true;
      script.onload = () => resolve(window.molstar);
      script.onerror = () => reject(new Error('Mol* did not load'));
      document.head.appendChild(script);
    });
    return loading;
  }

  async function ensurePlugin(container) {
    const molstar = await ensureLibrary();
    if (!molstar) throw new Error('Mol* unavailable');
    if (plugin) return plugin;
    plugin = await molstar.Viewer.create(container, {
      layoutIsExpanded: false,
      layoutShowControls: false,
      layoutShowSequence: false,
      layoutShowLog: false,
      layoutShowLeftPanel: false,
      viewportShowExpand: true,
      viewportShowSelectionMode: false,
      viewportShowAnimation: false,
      pdbProvider: 'rcsb',
      emdbProvider: 'rcsb',
    });
    contextCount += 1;
    return plugin;
  }

  /* Load one or two structures into the single instance.
     `entries`: [{url, format, label}]: the model first, then the experimental
     neighbour for the superposition panel. */
  async function show(container, entries) {
    const placeholder = container.querySelector('.placeholder');
    if (!entries || !entries.length) {
      if (placeholder) placeholder.hidden = false;
      return { ok: false, reason: 'nothing to show' };
    }
    try {
      const view = await ensurePlugin(container);
      if (placeholder) placeholder.hidden = true;
      await view.plugin.clear();
      for (const entry of entries) {
        if (entry.url) {
          await view.loadStructureFromUrl(entry.url, entry.format || 'pdb',
            Boolean(entry.binary));
        } else if (entry.pdbId) {
          await view.loadPdb(entry.pdbId);
        }
      }
      return { ok: true, contexts: contextCount };
    } catch (err) {
      if (placeholder) {
        placeholder.hidden = false;
        placeholder.textContent = 'The 3D viewer could not start in this browser.';
      }
      return { ok: false, reason: String(err && err.message) };
    }
  }

  /* Dispose on drawer close. This is the half that keeps the context count
     bounded; creating one instance is not enough on its own if it is never
     released when the drawer is dismissed for a long stretch. */
  function dispose() {
    if (!plugin) return;
    try {
      if (typeof plugin.dispose === 'function') plugin.dispose();
      else if (plugin.plugin && typeof plugin.plugin.dispose === 'function') {
        plugin.plugin.dispose();
      }
    } catch (err) { /* already gone */ }
    plugin = null;
    contextCount = Math.max(0, contextCount - 1);
  }

  window.WViewer = { show, dispose, stats: () => ({ contexts: contextCount,
    live: Boolean(plugin) }) };
})();
