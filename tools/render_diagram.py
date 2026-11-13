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
        '    "fontSize": "13px",\n'
        '    "primaryColor": "#FFFFFF",\n'
        '    "primaryTextColor": "#111827",\n'
        '    "primaryBorderColor": "#DDE2EB",\n'
        '    "lineColor": "#5B6472",\n'
        '    "clusterBkg": "#F1F3F7",\n'
        '    "clusterBorder": "#DDE2EB"\n'
        '  },\n'
        '  "securityLevel": "loose",\n'
        '  "flowchart": { "curve": "basis", "padding": 12 }\n}\n')

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    cmd = [binary, "-i", str(SOURCE), "-o", str(OUTPUT), "-c", str(CONFIG),
           "-b", "transparent", "--width", "1100"]
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
        f"{len(svg) // 1024} kB, {svg.count('href=\"#ref-')} reference anchors")
    return 0


if __name__ == "__main__":
    sys.exit(main())
