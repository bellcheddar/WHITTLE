/* Mol* viewer. DESIGN 3.6, and these are bugs if ignored:

   - **Exactly one plugin instance, ever.** Browsers cap WebGL contexts near 16.
     Without this the viewer dies silently after a dozen drawer opens and reads
     as a random failure.
   - Mol* (about 5 MB) is dynamic-imported on first drawer open, never in the
     initial bundle.
   - Structures are requested as .pdb and served with Content-Encoding gzip, so
     the browser decompresses and Mol* reads plain PDB.

   The drawer has two viewer containers, one in the Structure panel and one in
   the Neighbours panel. An earlier version kept "one instance" by caching the
   first plugin it built, which went wrong twice over: the cached plugin stayed
   bound to the container it was created against, so the second panel rendered
   nothing at all, and a concurrent second call built another context anyway,
   leaving two live. The instance is now keyed to its container, moving to a
   different container disposes the previous one first, and a creation already
   in flight is awaited rather than duplicated. The count stays at one and the
   picture appears where it was asked for.

   Adapted from GOBSMACKED/app/static/js and chatPDB-web (see REUSE.md): the
   single-instance discipline and the dispose-on-close are the parts carried
   over.
*/
'use strict';

(function () {
  let plugin = null;          // the one instance
  let boundTo = null;         // the container it was created against
  let creating = null;        // in-flight creation, so two calls cannot race
  let libraryPromise = null;  // in-flight dynamic import
  let created = 0;            // how many contexts this page has ever built
  let disposed = 0;

  const CSS_ID = 'molstar-css';

  /* Mol* replaces the children of the element it mounts into, so anything the
     page put inside the container is gone once the viewer starts. The
     placeholder therefore lives beside the container, and this looks in both
     places so an older layout still works. */
  /* `hidden` alone is not enough: both placeholders set `display` (one in the
     stylesheet, one inline) and any `display` beats the user agent's
     `[hidden] { display: none }`. Without the explicit style the placeholder
     paints on top of the model it was standing in for. */
  function setPlaceholder(node, visible, text) {
    if (!node) return;
    node.hidden = !visible;
    node.style.display = visible ? 'grid' : 'none';
    if (text !== undefined) node.textContent = text;
  }

  function placeholderFor(container) {
    return container.querySelector('.placeholder')
      || (container.parentElement
          && container.parentElement.querySelector(':scope > .placeholder'))
      || null;
  }

  async function ensureLibrary() {
    if (window.molstar) return window.molstar;
    if (libraryPromise) return libraryPromise;
    libraryPromise = new Promise((resolve, reject) => {
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
    return libraryPromise;
  }

  function teardown() {
    if (!plugin) return;
    try {
      if (typeof plugin.dispose === 'function') plugin.dispose();
      else if (plugin.plugin && typeof plugin.plugin.dispose === 'function') {
        plugin.plugin.dispose();
      }
      disposed += 1;
    } catch (err) { /* already gone */ }
    plugin = null;
    boundTo = null;
  }

  function buildPlugin(container) {
    // Everything from here to assigning `plugin` runs inside one promise that
    // is stored in `creating` **synchronously** by the caller. The previous
    // version awaited ensureLibrary() between checking `creating` and setting
    // it, and two calls made in the same tick both sailed through the guard,
    // built a context each and orphaned the first: created 2, disposed 0.
    return (async () => {
      const molstar = await ensureLibrary();
      if (!molstar) throw new Error('Mol* unavailable');
      const view = await molstar.Viewer.create(container, {
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
      plugin = view;
      boundTo = container;
      created += 1;
      try {
        // DESIGN 2: the canvas sits on #0D1117, not the library's white.
        view.plugin.canvas3d?.setProps({
          renderer: { backgroundColor: 0x0D1117 },
        });
      } catch (err) { /* the viewer still works on its own background */ }
      return view;
    })();
  }

  async function ensurePlugin(container) {
    // Same container with a live plugin: reuse it and reload into it.
    if (plugin && boundTo === container) return plugin;
    if (creating) {
      // Another call is already building one; wait rather than build a second.
      try { await creating; } catch (err) { /* handled by that caller */ }
      if (plugin && boundTo === container) return plugin;
    }
    // A different container needs its own context, so release the old one
    // first rather than hold two.
    if (plugin) teardown();
    creating = buildPlugin(container);
    try {
      return await creating;
    } finally {
      creating = null;
    }
  }

  /* Load one or two structures into the single instance.
     `entries`: [{url, format, binary} | {pdbId}]: the model first, then the
     experimental neighbour for the superposition panel. */
  async function show(container, entries) {
    const placeholder = placeholderFor(container);
    if (!entries || !entries.length) {
      setPlaceholder(placeholder, true);
      return { ok: false, reason: 'nothing to show' };
    }
    try {
      const view = await ensurePlugin(container);
      await view.plugin.clear();
      let loaded = 0;
      for (const entry of entries) {
        try {
          if (entry.url) {
            await view.loadStructureFromUrl(entry.url, entry.format || 'pdb',
              Boolean(entry.binary));
          } else if (entry.pdbId) {
            await view.loadPdb(entry.pdbId);
          }
          loaded += 1;
        } catch (err) {
          // One structure failing is not the whole panel failing: a missing
          // neighbour should still leave the model on screen.
        }
      }
      setPlaceholder(placeholder, loaded === 0);
      return { ok: loaded > 0, loaded, contexts: created - disposed };
    } catch (err) {
      setPlaceholder(placeholder, true,
        'The 3D viewer could not start in this browser.');
      return { ok: false, reason: String(err && err.message) };
    }
  }

  /* Load a MolViewSpec state, which is how the complex view highlights the
     residues PLIP found in contact.

     Mol* 4.9 has no MVS `primitives()` node, so the interactions cannot be
     drawn as dashed lines between atom pairs without moving to a 5.x bundle.
     Highlighting the contacting residues reads better at this viewer's size
     anyway, and the drawer carries the atom-level detail beside it.

     `build(builder, helpers)` receives a fresh MVS builder and returns
     nothing; the caller describes the scene and this handles the plugin. */
  async function showMvs(container, build) {
    const placeholder = placeholderFor(container);
    try {
      const view = await ensurePlugin(container);
      const mvs = window.molstar
        && window.molstar.PluginExtensions
        && window.molstar.PluginExtensions.mvs;
      if (!mvs) return { ok: false, reason: 'no mvs extension' };
      const builder = mvs.MVSData.createBuilder();
      build(builder);
      const state = builder.getState();
      await mvs.loadMVS(view.plugin, state, { replaceExisting: true });
      setPlaceholder(placeholder, false);
      return { ok: true, contexts: created - disposed };
    } catch (err) {
      setPlaceholder(placeholder, true,
        'The 3D viewer could not draw this complex.');
      return { ok: false, reason: String(err && err.message) };
    }
  }

  /* Dispose on drawer close. This is the half that keeps the context count
     bounded; building one instance is not enough on its own if it is never
     released when the drawer is dismissed. */
  function dispose() {
    teardown();
  }

  window.WViewer = {
    show,
    showMvs,
    dispose,
    stats: () => ({ contexts: created - disposed, created, disposed,
      live: Boolean(plugin) }),
  };
})();
