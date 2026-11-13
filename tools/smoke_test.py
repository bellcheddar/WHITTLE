"""Smoke test. Spec Section 10 item 1.

Every tab returns 200, every API route returns valid JSON, and the drawer loads
for a random 20 candidates per target. Runs against the Flask test client by
default, or against a live URL with --url.

Adapted from debris/tools/smoke_test.py (see REUSE.md): the routes and the
assertions are WHITTLE's.

    pixi run qc
    python tools/smoke_test.py --url https://whittle.mdeller.com
"""
from __future__ import annotations

import argparse
import json
import random
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

TARGETS = ["EGFR", "GLP1R", "PRLR", "TL1A", "TNFa"]


class Client:
    """One interface over the test client and a live host."""

    def __init__(self, url: str | None):
        self.url = url.rstrip("/") if url else None
        if self.url:
            import requests

            self.session = requests.Session()
            self.session.headers["User-Agent"] = "WHITTLE-smoke/1.0"
        else:
            from app.app import app

            app.config["TESTING"] = True
            self.client = app.test_client()

    def get(self, path: str) -> tuple[int, bytes]:
        if self.url:
            resp = self.session.get(self.url + path, timeout=45)
            return resp.status_code, resp.content
        resp = self.client.get(path)
        return resp.status_code, resp.data


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--url", default=None,
                        help="live host, for example https://whittle.mdeller.com")
    parser.add_argument("--per-target", type=int, default=20)
    parser.add_argument("--seed", type=int, default=20261113)
    args = parser.parse_args()

    client = Client(args.url)
    random.seed(args.seed)
    failures: list[str] = []
    checks = 0

    def check(name: str, ok: bool, detail: str = "") -> None:
        nonlocal checks
        checks += 1
        if not ok:
            failures.append(f"{name}: {detail}")
            print(f"  FAIL {name} {detail}")

    # 1. Every tab returns 200.
    for path in ("/", "/controls", "/about", "/healthz", "/robots.txt"):
        status, body = client.get(path)
        check(f"GET {path}", status == 200, f"status {status}")
        if path == "/healthz":
            try:
                payload = json.loads(body)
                check("healthz run id", bool(payload.get("run_id")),
                      "no run id in the health payload")
                print(f"  run {payload.get('run_id')} "
                      f"built {payload.get('built_at')} "
                      f"sha {payload.get('git_sha')} "
                      f"{payload.get('n_candidates')} candidates")
            except json.JSONDecodeError:
                check("healthz json", False, "not JSON")

    # 2. Every static bundle parses and reconciles.
    total_rows = 0
    for target in TARGETS:
        status, body = client.get(f"/static/data/candidates_{target}.json")
        check(f"bundle {target}", status == 200, f"status {status}")
        if status != 200:
            continue
        try:
            bundle = json.loads(body)
        except json.JSONDecodeError:
            check(f"bundle {target} json", False, "not JSON")
            continue
        check(f"{target} row count", bundle["n"] == 500, f"n = {bundle['n']}")
        total_rows += bundle["n"]
        # Spec 10 item 3: no candidate carries another target's rows.
        wrong = [r["id"] for r in bundle["rows"] if r["target"] != target]
        check(f"{target} rows all on target", not wrong,
              f"{len(wrong)} rows on the wrong target")

        # 3. The drawer loads for a random sample.
        sample = random.sample(bundle["rows"], min(args.per_target, len(bundle["rows"])))
        bad = []
        for row in sample:
            status, body = client.get(f"/api/candidate/{row['id']}.json")
            if status != 200:
                bad.append(f"{row['id']} status {status}")
                continue
            try:
                detail = json.loads(body)
            except json.JSONDecodeError:
                bad.append(f"{row['id']} not JSON")
                continue
            if detail.get("target") != target:
                bad.append(f"{row['id']} target {detail.get('target')}")
        check(f"{target} drawer x{len(sample)}", not bad, "; ".join(bad[:4]))

    check("total rows", total_rows == 2500, f"{total_rows} rows across the bundles")

    # 4. Meta and the column registry.
    status, body = client.get("/static/data/meta.json")
    check("meta bundle", status == 200, f"status {status}")
    if status == 200:
        meta = json.loads(body)
        check("41 columns", meta["columns"]["total"] == 41,
              f"{meta['columns']['total']} columns")
        check("6 default columns", len(meta["columns"]["defaults"]) == 6,
              f"{len(meta['columns']['defaults'])} defaults")
        check("5 targets", len(meta["targets"]) == 5, f"{len(meta['targets'])} targets")

    # 5. Export requires a target: there is no pooled view (spec 6.5 item 5).
    status, _ = client.get("/api/export.csv?target=EGFR")
    check("export EGFR", status == 200, f"status {status}")
    status, _ = client.get("/api/export.csv")
    check("export refuses a pooled request", status == 400, f"status {status}")

    # 6. A missing candidate is a 404, not a 500.
    status, _ = client.get("/api/candidate/NOPE-9999.json")
    check("unknown candidate 404", status == 404, f"status {status}")

    print(f"\n{checks - len(failures)}/{checks} checks passed")
    if failures:
        print("FAILURES:")
        for line in failures:
            print(f"  {line}")
        return 1
    print("smoke test clean")
    return 0


if __name__ == "__main__":
    sys.exit(main())
