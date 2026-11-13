#!/usr/bin/env python3
"""Re-download the self-hosted type into app/static/fonts/.

Run from the repo root:  python tools/fetch_fonts.py

DESIGN.md Section 2 asks for Inter throughout, with tabular numerals on every
numeric cell. Inter is served from our own droplet rather than from Google's CDN:
no third-party request on any page load, nothing to leak, and no layout shift
waiting on someone else's server. Latin and latin-ext only.

Adapted from debris/deploy/fetch_fonts.py (see REUSE.md): one family instead of
three, and the output goes to app/static/fonts/.
"""

from __future__ import annotations

import pathlib
import re
import urllib.request

OUT = pathlib.Path(__file__).resolve().parent.parent / "app" / "static" / "fonts"

# A modern browser UA, because the API serves woff2 to modern browsers and much larger
# legacy formats to anything it does not recognise.
UA = ("Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 "
      "(KHTML, like Gecko) Chrome/126.0.0.0 Safari/537.36")
CSS = ("https://fonts.googleapis.com/css2"
       "?family=Inter:wght@400..700"
       "&display=swap")

HEADER = (
    "/* WHITTLE -- self-hosted type. Inter throughout (DESIGN.md Section 2).\n"
    "   Variable woff2, latin + latin-ext only, served from our own droplet: no\n"
    "   third-party font request on any page load, and no layout shift waiting on\n"
    "   one. Regenerate with tools/fetch_fonts.py. */\n\n"
)


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    req = urllib.request.Request(CSS, headers={"User-Agent": UA})
    css = urllib.request.urlopen(req, timeout=60).read().decode()

    faces, seen = [], set()
    for subset, face in re.findall(r"/\*\s*([\w\-\[\]]+)\s*\*/\s*(@font-face\s*\{[^}]*\})", css):
        if subset not in ("latin", "latin-ext"):
            continue
        family = re.search(r"font-family:\s*'([^']+)'", face).group(1)
        weight = re.search(r"font-weight:\s*([^;]+);", face).group(1).strip()
        url = re.search(r"url\((https://[^)]+\.woff2)\)", face).group(1)
        unicode_range = re.search(r"unicode-range:\s*([^;]+);", face).group(1)

        slug = weight if weight.isdigit() else "var"
        name = f"{family.replace(' ', '')}-{slug}-{subset}.woff2"
        if name not in seen:
            urllib.request.urlretrieve(url, OUT / name)
            seen.add(name)
            print(f"  {name}  {(OUT / name).stat().st_size // 1024} kB  ({family} {weight})")

        faces.append(
            f"@font-face {{\n"
            f"  font-family: '{family}';\n"
            f"  font-style: normal;\n"
            f"  font-weight: {weight};\n"
            f"  font-display: swap;\n"
            f"  src: url('{name}') format('woff2');\n"
            f"  unicode-range: {unicode_range};\n"
            f"}}\n"
        )

    (OUT / "fonts.css").write_text(HEADER + "\n".join(faces))
    print(f"wrote {OUT / 'fonts.css'} ({len(faces)} faces, {len(seen)} files)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
