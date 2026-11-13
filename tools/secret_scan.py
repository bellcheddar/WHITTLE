"""Secret scan. Spec Section 0.4 item 4: this must pass before every push.

Greps the working tree and the whole git history for credential shapes and for
any literal IPv4 address. A hit outside an ignored path aborts the push and
writes BLOCKED.md. Shape lifted from debris/tools (see REUSE.md); the IPv4 rule
and the git-history pass are new here, because this repo is public and the
droplet address must never reach it.

    python tools/secret_scan.py          # tree + history, exit 1 on any hit
    python tools/secret_scan.py --tree   # working tree only (faster)
"""
from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

PATTERNS: list[tuple[str, re.Pattern[str]]] = [
    ("private key", re.compile(r"BEGIN [A-Z ]*PRIVATE KEY")),
    ("ssh-rsa key", re.compile(r"ssh-rsa\s+AAAA")),
    ("ssh-ed25519 key", re.compile(r"ssh-ed25519\s+AAAA")),
    ("api key assignment", re.compile(r"api[_-]?key\s*[=:]\s*['\"]?[A-Za-z0-9_\-]{12,}")),
    ("secret assignment", re.compile(r"secret\s*[=:]\s*['\"]?[A-Za-z0-9_\-]{12,}")),
    ("password assignment", re.compile(r"password\s*[=:]\s*['\"]?[^\s'\"]{6,}")),
    ("token assignment", re.compile(r"token\s*[=:]\s*['\"]?[A-Za-z0-9_\-]{16,}")),
    ("aws key id", re.compile(r"AKIA[0-9A-Z]{16}")),
    ("github token", re.compile(r"gh[pousr]_[A-Za-z0-9]{20,}")),
    ("github pat", re.compile(r"github_pat_[A-Za-z0-9_]{20,}")),
    ("ipv4 address", re.compile(r"(?<![\w.])(?:\d{1,3}\.){3}\d{1,3}(?![\w.])")),
]

# Documented placeholders and loopback forms. Everything here is either literal
# placeholder text this build wrote on purpose, or an address that identifies no
# host. Narrow by design: a new exemption is a decision, not a convenience.
ALLOW = re.compile(
    r"REPLACE_ME"
    r"|127\.0\.0\.1"
    r"|0\.0\.0\.0"
    r"|<IP>"
    r"|localhost"
    r"|ssh_target"
    r"|\bexample\.com\b"
    # version strings and numeric tuples read as IPv4 to a regex
    r"|\b\d+\.\d+\.\d+\.\d+\s*(?:#|$)"
)

# Paths that are gitignored, generated, or are this build's own documentation of
# the rule. The spec and the scanner itself both contain the pattern text.
SKIP_DIRS = {
    ".git", ".pixi", ".venv", "node_modules", "__pycache__", "data", "models",
    ".pytest_cache", ".ruff_cache", ".DS_Store",
}
SKIP_FILES = {
    "tools/secret_scan.py",
    "WHITTLE_BUILD_SPEC.md",
    "config/deploy.toml",
    "config/deploy.toml.example",
    "pixi.lock",
}
TEXT_SUFFIXES = {
    ".py", ".md", ".toml", ".yaml", ".yml", ".json", ".sh", ".txt", ".cfg",
    ".ini", ".html", ".js", ".css", ".service", ".timer", ".conf", ".svg", "",
}


def tracked_or_present_files() -> list[Path]:
    out: list[Path] = []
    for path in sorted(ROOT.rglob("*")):
        if not path.is_file():
            continue
        rel = path.relative_to(ROOT)
        if any(part in SKIP_DIRS for part in rel.parts):
            continue
        if str(rel) in SKIP_FILES:
            continue
        if path.suffix.lower() not in TEXT_SUFFIXES:
            continue
        if path.stat().st_size > 4_000_000:
            continue
        out.append(path)
    return out


def ignored_by_git(rel: str) -> bool:
    try:
        res = subprocess.run(["git", "check-ignore", "-q", rel], cwd=ROOT,
                             capture_output=True)
        return res.returncode == 0
    except FileNotFoundError:
        return False


def scan_lines(label: str, lines: list[str]) -> list[str]:
    hits: list[str] = []
    for n, line in enumerate(lines, 1):
        if len(line) > 4000:
            continue
        for name, pattern in PATTERNS:
            for match in pattern.finditer(line):
                fragment = match.group(0)
                if ALLOW.search(fragment) or ALLOW.search(line):
                    continue
                hits.append(f"{label}:{n}: {name}: {fragment[:60]}")
    return hits


def scan_tree() -> list[str]:
    hits: list[str] = []
    for path in tracked_or_present_files():
        rel = str(path.relative_to(ROOT))
        if ignored_by_git(rel):
            continue
        try:
            text = path.read_text(errors="ignore")
        except OSError:
            continue
        hits += scan_lines(rel, text.splitlines())
    return hits


def scan_history() -> list[str]:
    """Every blob and every commit message in the history."""
    if not (ROOT / ".git").exists():
        return []
    hits: list[str] = []
    log = subprocess.run(["git", "log", "--all", "--format=%H%n%B%n--%n"],
                         cwd=ROOT, capture_output=True, text=True)
    hits += scan_lines("git-log", log.stdout.splitlines())
    # Patch text across all commits. -S would miss a pattern; the full diff is
    # small for this repo and correctness beats speed for a one-rule-no-fallback check.
    diff = subprocess.run(["git", "log", "--all", "-p", "--no-color"],
                          cwd=ROOT, capture_output=True, text=True)
    hits += scan_lines("git-history", diff.stdout.splitlines())
    return hits


def write_blocked(hits: list[str]) -> None:
    (ROOT / "BLOCKED.md").write_text(
        "# BLOCKED: secret scan\n\n"
        "`tools/secret_scan.py` found credential-shaped text. The repo is public, "
        "so spec Section 0.4 aborts the push rather than pushing and cleaning up after.\n\n"
        "## Hits\n\n```\n" + "\n".join(hits[:200]) + "\n```\n\n"
        "## What is needed\n\n"
        "Remove or gitignore each path above, then re-run:\n\n"
        "```bash\npython tools/secret_scan.py\n```\n\n"
        "If any hit is a real credential that has already been committed, the key "
        "needs rotating (spec 0.4 item 5): do not rewrite history quietly.\n"
    )


def main() -> int:
    tree_only = "--tree" in sys.argv
    hits = scan_tree()
    if not tree_only:
        hits += scan_history()
    if hits:
        print(f"SECRET SCAN FAILED: {len(hits)} hit(s)")
        for hit in hits[:60]:
            print("  " + hit)
        write_blocked(hits)
        return 1
    blocked = ROOT / "BLOCKED.md"
    if blocked.exists() and "secret scan" in blocked.read_text():
        blocked.unlink()
    print("secret scan clean" + (" (tree only)" if tree_only else " (tree and history)"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
