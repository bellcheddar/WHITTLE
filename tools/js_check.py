#!/usr/bin/env python3
"""Parse every shipped JavaScript file, and fail if any of them will not parse.

    ./tools/js_check.py

This exists because of one defect that reached the repository. `tile()` in
`app/static/js/candidates.js` takes `label` as its first parameter, and a
four-state verdict map added on 2026-11-14 introduced `const label =
FLAGS[verdict]` in the same scope. Redeclaring a parameter with `const` is a
SyntaxError, not a shadowed variable, so it did not break that one tile: it
stopped the entire file parsing, and with it every interactive part of the
Candidates page, which is the page the site is for.

Nothing caught it. The tests are Python, `style_check.py` reads prose, and the
Flask app serves a static file without parsing it, so a dead bundle looks
exactly like a healthy one from the server side. A browser would have said so
immediately, and no step between writing the file and deploying it opened a
browser.

Parsing is deliberately all this does. It will not catch a reference to an
undefined variable or a wrong property name; it catches the class of fault that
silently removes a whole page, which is the one that got through. Run it before
every push, beside `style_check.py` and `secret_scan.py`.

A missing `node` is reported as a failure rather than a pass. A check that
quietly succeeds when it did not run is how the original defect survived, and
the same mistake in the checker would be worse than having no checker.
"""
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
TARGETS = sorted((ROOT / "app/static/js").glob("*.js"))


def main() -> int:
    node = shutil.which("node")
    if node is None:
        print("js check COULD NOT RUN: node is not on PATH, so the shipped "
              "JavaScript was not parsed. This is a failure and not a pass: "
              "a SyntaxError in one bundle takes a whole page down and the "
              "server cannot tell. Install node, or parse the files another "
              "way before deploying.", file=sys.stderr)
        return 1

    if not TARGETS:
        print("js check: no JavaScript found under app/static/js, which is "
              "itself wrong", file=sys.stderr)
        return 1

    failures: list[tuple[Path, str]] = []
    for path in TARGETS:
        done = subprocess.run([node, "--check", str(path)],
                              capture_output=True, text=True)
        if done.returncode != 0:
            detail = (done.stderr or done.stdout or "").strip().splitlines()
            # node prints the offending line, a caret, then the message. The
            # message is the useful line and the stack trace below it is not.
            message = next((ln.strip() for ln in detail
                            if "Error" in ln), detail[0] if detail else "?")
            failures.append((path, message))

    for path, message in failures:
        print(f"js check FAIL {path.relative_to(ROOT)}: {message}",
              file=sys.stderr)
    if failures:
        print(f"\n{len(failures)} of {len(TARGETS)} JavaScript files will not "
              f"parse. The pages that load them are dead in a browser, whatever "
              f"the server returns.", file=sys.stderr)
        return 1
    print(f"js check clean: {len(TARGETS)} files parse "
          f"({', '.join(p.name for p in TARGETS)})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
