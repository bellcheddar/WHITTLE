// Build-time asset bundler. Runs on the Studio only: the droplet gets the output.
import { build } from 'esbuild';
import { cp, mkdir, stat } from 'node:fs/promises';
import { fileURLToPath } from 'node:url';
import path from 'node:path';

const root = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '../..');
const vendor = path.join(root, 'app/static/vendor');
await mkdir(vendor, { recursive: true });

// 1. The partial Plotly bundle.
await build({
  entryPoints: [path.join(root, 'tools/bundle/plotly-whittle.js')],
  bundle: true,
  minify: true,
  format: 'iife',
  target: ['es2019'],
  // Plotly's CommonJS dependencies reference `global`, which does not exist in a
  // browser. Without this the bundle throws "global is not defined" on load and
  // the explorer never draws.
  define: { global: 'globalThis' },
  outfile: path.join(vendor, 'plotly-whittle.min.js'),
  logLevel: 'info',
});

// 2. Mol*: the prebuilt viewer bundle and its stylesheet, dynamic-imported on
//    first drawer open (DESIGN 3.6), never in the initial bundle.
const molstar = path.join(root, 'node_modules/molstar/build/viewer');
for (const [from, to] of [
  ['molstar.js', 'molstar.js'],
  ['molstar.css', 'molstar.css'],
]) {
  await cp(path.join(molstar, from), path.join(vendor, to));
}

for (const name of ['plotly-whittle.min.js', 'molstar.js', 'molstar.css']) {
  const { size } = await stat(path.join(vendor, name));
  console.log(`  ${name}  ${Math.round(size / 1024)} kB`);
}
