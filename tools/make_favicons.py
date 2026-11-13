"""Favicon set from the house icon. Spec Section 11, Phase 4.

The full icon carries its wordmark, which is illegible below about 64 px, so
the favicon set is built from a glyph-only variant with the artwork scaled to
fill the canvas. Everything is generated from `icons/WHITTLE_icon.svg`, so the
set cannot drift from the icon.

    pixi run python tools/make_favicons.py
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

SOURCE = ROOT / "icons" / "WHITTLE_icon.svg"
STATIC = ROOT / "app" / "static"
GLYPH_SVG = ROOT / "icons" / "WHITTLE_glyph.svg"

# The glyph band is x 30-270, y 24-224. Re-framed on its own, with a little
# padding, it fills a square favicon instead of sitting above empty space.
GLYPH_VIEWBOX = "34 20 232 212"
SIZES = {"icon-16.png": 16, "icon-32.png": 32, "icon-48.png": 48,
         "icon-180.png": 180, "apple-touch-icon.png": 180,
         "icon-192.png": 192, "icon-512.png": 512}


def build_glyph() -> str:
    svg = SOURCE.read_text()
    # Drop the wordmark paths and re-frame on the glyph band.
    svg = re.sub(r'\s*<path class="wordmark-(?:outline|fill)"[^>]*/>', "", svg)
    svg = svg.replace('viewBox="0 0 300 300"', f'viewBox="{GLYPH_VIEWBOX}"')
    svg = svg.replace('aria-label="WHITTLE: a triage funnel narrowing a '
                      'repertoire of candidate nanobodies down to a shortlist"',
                      'aria-label="WHITTLE"')
    return svg


def main() -> int:
    from pipeline.common import log

    if not SOURCE.exists():
        log(f"make_favicons: {SOURCE} is absent")
        return 1
    STATIC.mkdir(parents=True, exist_ok=True)
    glyph = build_glyph()
    GLYPH_SVG.write_text(glyph)
    (STATIC / "favicon.svg").write_text(glyph)

    import cairosvg
    from PIL import Image

    written = []
    for name, size in SIZES.items():
        cairosvg.svg2png(bytestring=glyph.encode(),
                         write_to=str(STATIC / name),
                         output_width=size, output_height=size)
        written.append(f"{name} ({size}px)")

    # A multi-resolution .ico, which is what a browser asks for by default.
    frames = [Image.open(STATIC / f"icon-{n}.png").convert("RGBA")
              for n in (16, 32, 48)]
    frames[0].save(STATIC / "favicon.ico", format="ICO",
                   sizes=[(16, 16), (32, 32), (48, 48)])
    written.append("favicon.ico (16, 32, 48)")

    (STATIC / "site.webmanifest").write_text(
        '{\n'
        '  "name": "WHITTLE",\n'
        '  "short_name": "WHITTLE",\n'
        '  "description": "Whole-repertoire hierarchical in-silico triage of '
        'nanobody candidates: thermostability, liabilities, expression.",\n'
        '  "icons": [\n'
        '    { "src": "/static/icon-192.png", "sizes": "192x192", '
        '"type": "image/png" },\n'
        '    { "src": "/static/icon-512.png", "sizes": "512x512", '
        '"type": "image/png" },\n'
        '    { "src": "/static/favicon.svg", "type": "image/svg+xml", '
        '"sizes": "any" }\n'
        '  ],\n'
        '  "theme_color": "#3457A6",\n'
        '  "background_color": "#F1F3F7",\n'
        '  "display": "browser",\n'
        '  "start_url": "/"\n'
        '}\n')
    written.append("site.webmanifest")
    log(f"make_favicons: {len(written)} files into app/static/: "
        f"{', '.join(written)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
