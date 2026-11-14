"""Render the workflow diagram to a static SVG. Spec 8.3.

"Author as Mermaid, render to static SVG at build time with mermaid-cli; ship
the SVG. Do not load the Mermaid runtime on the droplet."

Each node is an anchor into its reference row on the About tab. mermaid-cli
emits `<a>` elements for `click` directives, and this script rewrites them to
in-page anchors and strips the title tooltips that would otherwise shadow the
accessible names.

    pixi run python tools/render_diagram.py
"""
from __future__ import annotations

import re
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

SOURCE = ROOT / "docs" / "workflow.mmd"
OUTPUT = ROOT / "app" / "static" / "img" / "workflow.svg"
CONFIG = ROOT / "docs" / "mermaid.json"


def main() -> int:
    from pipeline.common import log

    if not SOURCE.exists():
        log(f"render_diagram: {SOURCE} is absent")
        return 1
    binary = shutil.which("mmdc") or str(ROOT / "node_modules" / ".bin" / "mmdc")
    if not Path(binary).exists():
        log("render_diagram: mermaid-cli (mmdc) is not installed; the About tab "
            "falls back to its tier table, which carries the same information")
        return 0

    CONFIG.write_text(
        '{\n  "theme": "base",\n'
        '  "themeVariables": {\n'
        '    "fontFamily": "Inter, system-ui, sans-serif",\n'
        '    "fontSize": "14px",\n'
        '    "primaryColor": "#FFFFFF",\n'
        '    "primaryTextColor": "#111827",\n'
        '    "primaryBorderColor": "#DDE2EB",\n'
        '    "lineColor": "#5B6472",\n'
        '    "clusterBkg": "#F1F3F7",\n'
        '    "clusterBorder": "#DDE2EB"\n'
        '  },\n'
        '  "securityLevel": "loose",\n'
        '  "flowchart": { "curve": "basis", "padding": 16, "nodeSpacing": 34,'
        ' "rankSpacing": 44, "wrappingWidth": 290 }\n}\n')

    # mermaid measures each label in the browser it renders in, then sizes the
    # node box to the measurement. Headless Chrome has no Inter installed, so it
    # measured in a fallback face and the page, which does have Inter, drew
    # wider text into boxes sized for a narrower one: every long label was
    # clipped at the node edge. Handing mmdc the same self-hosted Inter the site
    # serves makes the measurement and the rendering agree.
    # The face is inlined as a data URI rather than referenced by path: a
    # file:// subresource is not reliably fetched by the headless page, and a
    # font that fails to load silently reverts the measurement to a fallback,
    # which is the failure this exists to prevent.
    import base64

    font_css = ROOT / "docs" / "mermaid-fonts.css"
    faces = []
    for woff in sorted((ROOT / "app" / "static" / "fonts").glob("*latin.woff2")):
        b64 = base64.b64encode(woff.read_bytes()).decode()
        faces.append(
            "@font-face { font-family: 'Inter'; font-style: normal;\n"
            "  font-weight: 400 700; font-display: block;\n"
            f"  src: url(data:font/woff2;base64,{b64}) format('woff2'); }}")
    faces.append("body, .nodeLabel, .edgeLabel, .cluster-label "
                 "{ font-family: 'Inter', system-ui, sans-serif; }")
    font_css.write_text("\n".join(faces) + "\n")

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    cmd = [binary, "-i", str(SOURCE), "-o", str(OUTPUT), "-c", str(CONFIG),
           "--cssFile", str(font_css),
           "-b", "transparent", "--width", "1280"]
    # mermaid-cli renders through puppeteer, which wants its own Chrome
    # download. Point it at a browser that is already on the machine instead:
    # this is a build-time step on the Studio and nothing here ships.
    import os

    env = dict(os.environ)
    if "PUPPETEER_EXECUTABLE_PATH" not in env:
        for candidate in (
                "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome",
                "/Applications/Chromium.app/Contents/MacOS/Chromium",
                "/Applications/Microsoft Edge.app/Contents/MacOS/Microsoft Edge"):
            if Path(candidate).exists():
                env["PUPPETEER_EXECUTABLE_PATH"] = candidate
                break
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=600,
                              env=env)
    except Exception as exc:
        log(f"render_diagram: mmdc failed ({type(exc).__name__}); the About tab "
            f"falls back to its tier table")
        return 0
    if proc.returncode != 0 or not OUTPUT.exists():
        log(f"render_diagram: mmdc returned {proc.returncode}: "
            f"{(proc.stderr or '')[-300:]}")
        return 0

    svg = OUTPUT.read_text()
    # mermaid-cli writes click targets as xlink:href on an <a>; make them
    # in-page anchors and give the figure an accessible name.
    # mermaid-cli pins an inline `max-width` equal to the diagram's native
    # width, which stops the About tab scaling it up to fill its figure box.
    # The page wants it slightly larger than native for readability, so the
    # inline cap comes off and the stylesheet sets the ceiling instead.
    svg = re.sub(r'max-width:\s*[\d.]+px;\s*', "", svg, count=1)
    # mermaid sizes every label box from a measurement taken in the headless
    # browser, which has no Inter installed and therefore measures in a
    # fallback face. Inter renders about 3 per cent wider on a machine that
    # does have it, and a foreignObject clips at its own width, so long labels
    # lost their last few characters: "co-folding" became "co-foldinc".
    # mermaid-cli's --cssFile lands in the OUTPUT svg rather than in the page
    # it measures, so it cannot fix this. Each label box is widened here
    # instead, and its parent group shifted left by half the growth so the
    # text stays centred. Only the clip region changes; no node moves.
    #
    # This is a text substitution on purpose. Parsing the SVG and writing it
    # back mangles the XHTML inside each foreignObject, because those children
    # are in a second namespace and come back out prefixed, which costs every
    # <br/> in the diagram.
    def _widen(match: re.Match) -> str:
        opening, x, y, width = match[1], float(match[2]), float(match[3]), float(match[4])
        if width <= 0:                      # empty edge labels carry width 0
            return match[0]
        grown = width * 1.08 + 2
        moved = opening.replace(
            f"translate({match[2]}, {match[3]})",
            f"translate({x - (grown - width) / 2:.4f}, {y})")
        return f'{moved}<foreignObject width="{grown:.4f}"'

    # Node labels carry an empty <rect/> between the group and the
    # foreignObject; cluster and edge labels do not.
    svg, widened = re.subn(
        r'(<g [^>]*transform="translate\(([-\d.]+), ([-\d.]+)\)"[^>]*>'
        r'(?:<rect\s*/>)?)<foreignObject width="([\d.]+)"',
        _widen, svg)

    svg = svg.replace('xlink:href="#ref-', 'href="#ref-')
    svg = re.sub(r'\s(?:xlink:)?title="[^"]*"', "", svg)
    if "<title>" not in svg:
        svg = svg.replace(
            "<svg ",
            '<svg role="img" aria-label="WHITTLE workflow: 2,500 candidate '
            'sequences through Tier 0 sequence checks, Tier 1 structure and '
            'developability, the Tier 1b thermostability model and Tier 2 '
            'co-folding, into within-target percentiles and a Pareto front" ', 1)
    OUTPUT.write_text(svg)
    log(f"render_diagram: {OUTPUT.relative_to(ROOT)} written, "
        f"{len(svg) // 1024} kB, {svg.count('href=\"#ref-')} reference anchors, "
        f"{widened} label boxes widened for Inter")
    return 0


if __name__ == "__main__":
    sys.exit(main())
