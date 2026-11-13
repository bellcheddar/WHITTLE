"""House style check. Spec 0.1 item 11 and DESIGN.md Section 7.

British English, no em dashes, none of the banned marketing words, and app
copy that says "predicted" or "estimated" rather than "confirmed" or
"validated". Checks the things a reader would notice, in the files a reader
sees: templates, JavaScript strings, the generated Markdown and the config
that feeds the About tab.

    pixi run python tools/style_check.py
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

TARGETS = [
    "app/templates", "app/static/js", "app/static/css", "docs",
    "README.md", "FINDINGS.md", "MODEL_CARD.md", "PROGRESS.md",
    "DECISIONS.md", "REUSE.md", "config/references.yaml", "config/targets.yaml",
    "pipeline", "tools", "deploy",
]
SKIP_NAMES = {"vendor", "node_modules", "__pycache__", ".pixi", "data"}

BANNED_WORDS = [
    "groundbreaking", "revolutionary", "paradigm-shifting", "game-changing",
    "cutting-edge", "unprecedented", "seamless", "delve",
]
# "leverage" only as a verb; the noun is legitimate and does not appear here.
BANNED_PATTERNS = [
    (re.compile(r"—"), "em dash (use a colon or parentheses)"),
    (re.compile(r"\bleverages?\b|\bleveraging\b", re.I), "leverage as a verb"),
]
# Overstatement in anything a visitor reads.
# "cross-validated" and "cross-validation" are technical terms, not claims.
OVERSTATEMENT = re.compile(
    r"(?<!cross-)\b(confirmed|validated|proven|guarantee[sd]?)\b", re.I)
USER_FACING = ("app/templates", "app/static/js", "README.md", "FINDINGS.md")

# American spellings that have a British form this project uses.
AMERICANISMS = {
    r"\bcolor\b": "colour", r"\bcolors\b": "colours", r"\bbehavior\b": "behaviour",
    r"\banalyze\b": "analyse", r"\banalyzed\b": "analysed",
    r"\bnormalize\b": "normalise", r"\bnormalized\b": "normalised",
    r"\brecognize\b": "recognise", r"\bsummarize\b": "summarise",
    r"\bcenter\b": "centre", r"\bfavor\b": "favour", r"\blicense\b": "licence",
    r"\bmodeling\b": "modelling", r"\bmodeled\b": "modelled",
    r"\blabeled\b": "labelled", r"\bfulfill\b": "fulfil",
}
# Exempt: CSS properties, API names and library identifiers are not prose.
CODE_EXEMPT = re.compile(
    r"color\s*[:=]|colour|-color|colorscale|bgcolor|gridcolor|zerolinecolor"
    r"|linecolor|fillcolor|paper_bgcolor|plot_bgcolor|textColor|BorderColor"
    r"|clusterBkg|primaryColor|accent-color|theme_color|background_color"
    r"|License-Identifier|license\s*=|licenses|normalize\(|analyze_|\.analyze"
    r"|center\s*[:=]|text-align|place-items|justify-content|align-items"
    r"|normalized\(|center\"|\"center\"|labeled=|labels|normalise")


# Files this build was handed rather than wrote. They are inputs and are not
# rewritten to house style: config/targets.yaml in particular is authoritative
# (spec 5.2) and editing it, even a comment, would be editing the brief.
GIVEN = {"config/targets.yaml", "config/references.yaml",
         "config/deploy.toml.example", "launch.sh"}


def files() -> list[Path]:
    out: list[Path] = []
    for target in TARGETS:
        path = ROOT / target
        if path.is_file():
            out.append(path)
            continue
        if not path.is_dir():
            continue
        for item in path.rglob("*"):
            if not item.is_file():
                continue
            if any(part in SKIP_NAMES for part in item.parts):
                continue
            if item.suffix.lower() in {".py", ".js", ".html", ".css", ".md",
                                       ".yaml", ".yml", ".sh", ".mmd", ".conf"}:
                out.append(item)
    return sorted(set(out))


def main() -> int:
    issues: list[str] = []
    checked = 0
    for path in files():
        rel = str(path.relative_to(ROOT))
        if rel in GIVEN or rel == "tools/style_check.py":
            continue
        try:
            text = path.read_text(errors="ignore")
        except OSError:
            continue
        checked += 1
        user_facing = any(rel.startswith(prefix) for prefix in USER_FACING)
        for n, line in enumerate(text.splitlines(), 1):
            lowered = line.lower()
            for word in BANNED_WORDS:
                if word in lowered:
                    issues.append(f"{rel}:{n}: banned word '{word}'")
            for pattern, label in BANNED_PATTERNS:
                if pattern.search(line):
                    issues.append(f"{rel}:{n}: {label}")
            for pattern, british in AMERICANISMS.items():
                if re.search(pattern, line, re.I) and not CODE_EXEMPT.search(line):
                    issues.append(f"{rel}:{n}: American spelling, use '{british}'"
                                  f": {line.strip()[:70]}")
            if user_facing:
                match = OVERSTATEMENT.search(line)
                # "not confirmed", "never validated" and the like are the point.
                if match and not re.search(
                        r"\b(not|never|no|cannot|rather than|instead of)\b.{0,40}"
                        + re.escape(match.group(0)), line, re.I):
                    issues.append(
                        f"{rel}:{n}: overstatement '{match.group(0)}' in copy a "
                        f"visitor reads: {line.strip()[:70]}")

    print(f"style check: {checked} files")
    if issues:
        print(f"{len(issues)} issue(s):")
        for line in issues[:60]:
            print("  " + line)
        return 1
    print("house style clean: British English, no em dashes, no banned words, "
          "no overstatement in visitor-facing copy")
    return 0


if __name__ == "__main__":
    sys.exit(main())
