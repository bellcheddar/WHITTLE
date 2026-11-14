#!/usr/bin/env python3
"""Compile every Jinja template, and render each against the states the app serves.

    ./tools/template_check.py

The companion to `tools/js_check.py`, and it exists for the same reason. A
broken bundle is served with a 200 and a broken template raises at request
time, so neither fault appears until a browser asks for the page, and nothing
between writing a template and deploying it ever asked.

Two levels, because they catch different things. Compiling finds a syntax
error: an unclosed block, a malformed macro, a stray brace. Rendering finds
the faults that compile perfectly and still fail, which on this project means
a dict method called on a Jinja Undefined. `{{ meta.axis_verdicts | tojson }}`
is the case that was live: the key went into meta.json and into
candidates.html without going into the pre-build stub in `app.meta()`, and
`tojson` on an Undefined raises, so the shell that stub exists to serve would
not have rendered at all before the first build-db run.

Both states the app really produces are rendered, and only those two. The
pre-build stub is parsed out of `app/app.py` rather than retyped, because a
copy drifts from what it mirrors. The built `meta.json` is read as served.

An earlier version of this file invented a third state, handing templates a
self-returning mock, and because that mock carried a `.get` it swallowed the
very fault described above: the check passed and the page broke. A mock
forgiving enough to never raise cannot test anything that raises.
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
TEMPLATES = ROOT / "app/templates"


def fallback_meta() -> dict:
    """The pre-build stub from `app.meta()`, parsed out of the source."""
    import ast
    tree = ast.parse((ROOT / "app/app.py").read_text())
    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef) and node.name == "meta":
            for inner in ast.walk(node):
                if isinstance(inner, ast.Dict):
                    try:
                        value = ast.literal_eval(inner)
                    except ValueError:
                        continue
                    if isinstance(value, dict) and "run_id" in value:
                        return value
    return {}


def built_meta() -> dict | None:
    """The real meta.json, which is what the app serves after a build."""
    import json
    path = ROOT / "app/static/data/meta.json"
    if not path.is_file():
        return None
    try:
        return json.loads(path.read_text())
    except ValueError:
        return None


def context(meta: dict) -> dict:
    """Mirrors `app.inject()` plus the per-route names.

    The route-level collections are empty because that is a state the app
    really serves: a route whose bundle has not been written passes empty
    lists. `meta` is a plain dict so a missing key becomes Undefined exactly
    as Flask makes it.
    """
    return {
        "tabs": [], "active": "/", "meta": meta, "backronym": "",
        "asset_version": lambda relative: relative, "indexable": False,
        "about": {}, "rows": [], "columns": [], "targets": [],
        "target": None, "candidates": [], "qc": {}, "controls": [],
        "axes": {},
        "request": {}, "url_for": lambda *a, **k: "#",
        "config": {}, "session": {},
        "get_flashed_messages": lambda **k: [],
    }


def main() -> int:
    try:
        from jinja2 import Environment, FileSystemLoader
    except ImportError:
        print("template check COULD NOT RUN: jinja2 is not importable. This "
              "is a failure and not a pass.", file=sys.stderr)
        return 1

    env = Environment(loader=FileSystemLoader(str(TEMPLATES)))
    names = sorted(p.name for p in TEMPLATES.glob("*.html"))
    if not names:
        print("template check: no templates found", file=sys.stderr)
        return 1

    states = [("pre-build fallback", fallback_meta())]
    built = built_meta()
    if built is None:
        print("template check: app/static/data/meta.json is absent, so only "
              "the pre-build state was rendered. Run pipeline.build_db.",
              file=sys.stderr)
    else:
        states.append(("built meta.json", built))

    failures: list[tuple[str, str]] = []
    for name in names:
        try:
            template = env.get_template(name)
        except Exception as exc:                      # noqa: BLE001
            failures.append((name, f"will not compile: {exc}"))
            continue
        for label, meta in states:
            try:
                template.render(**context(meta))
            except Exception as exc:                  # noqa: BLE001
                failures.append((name, f"compiles but will not render against "
                                       f"the {label}: "
                                       f"{type(exc).__name__}: {exc}"))

    for name, message in failures:
        print(f"template check FAIL {name}: {message}", file=sys.stderr)
    if failures:
        print(f"\n{len(failures)} template failure(s) of {len(names)} "
              f"templates across {len(states)} state(s). The pages that use "
              f"them raise at request time.", file=sys.stderr)
        return 1
    print(f"template check clean: {len(names)} templates compile and render "
          f"against {len(states)} state(s) ({', '.join(names)})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
