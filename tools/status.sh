#!/usr/bin/env bash
# WHITTLE queue status. Spec 6.4 item 7.
# Adapted from NOBBLE/tools/status.sh (see REUSE.md).
#
#   ./tools/status.sh          # one report
#   ./tools/status.sh --watch  # refresh every 30 seconds
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

report() {
  "$ROOT/.pixi/envs/default/bin/python" - <<'PYEOF'
import json
import datetime as dt
from pathlib import Path

ROOT = Path(__file__).resolve().parents[0] if False else Path.cwd()
M = ROOT / "data" / "manifests"

def read(stage):
    path = M / f"{stage}.jsonl"
    if not path.exists():
        return []
    out = []
    for line in path.read_text().splitlines():
        if line.strip():
            try:
                out.append(json.loads(line))
            except json.JSONDecodeError:
                pass
    return out

print(f"WHITTLE status  {dt.datetime.now():%Y-%m-%d %H:%M}")
print("=" * 64)

STAGES = [
    ("ingest", "targets ingested", 5),
    ("tier0", "targets through Tier 0", 5),
    ("fold", "models built", 2500),
    ("patches", "surface patches", 2500),
    ("tnp", "TNP profiles", 2500),
    ("nativeness", "AbLang2 nativeness", 2500),
    ("humanness", "OASis humanness", 2500),
    ("neighbours", "nearest experimental neighbour", 2500),
    ("tm_predict", "Tm predictions", 1),
    ("msa", "antigen MSAs cached", 5),
    ("cofold", "co-folds", 100),
    ("interface", "interface analyses", 5),
    ("score", "targets scored", 5),
]
for stage, label, expected in STAGES:
    rows = read(stage)
    # Count distinct ids, not manifest lines: a stage that retried a row writes
    # two lines for it, and counting lines reports a retried row as a permanent
    # failure long after it succeeded.
    succeeded, failed = set(), set()
    for r in rows:
        key = str(r.get("id", len(succeeded) + len(failed)))
        if str(r.get("status", "ok")).startswith("failed"):
            failed.add(key)
        else:
            succeeded.add(key)
    failed -= succeeded
    ok, bad = len(succeeded), len(failed)
    pct = 100.0 * ok / expected if expected else 0.0
    bar = "#" * int(min(pct, 100) / 5) + "." * (20 - int(min(pct, 100) / 5))
    flag = f"  {bad} still failing" if bad else ""
    print(f"{label:34} [{bar}] {ok:>5}/{expected:<5} {pct:5.1f}%{flag}")

cofold = read("cofold")
ok = [r for r in cofold if r.get("status") == "ok" and r.get("seconds")]
if ok:
    by_target = {}
    for r in ok:
        by_target.setdefault(r.get("target", "?"), []).append(r["seconds"])
    print()
    print("co-fold wall-clock, minutes per complex")
    for target, seconds in by_target.items():
        mean = sum(seconds) / len(seconds)
        print(f"  {target:8} n={len(seconds):<4} mean {mean/60:6.1f}  "
              f"min {min(seconds)/60:5.1f}  max {max(seconds)/60:6.1f}")

blocked = ROOT / "BLOCKED.md"
print()
print("BLOCKED.md present: the build needs you" if blocked.exists()
      else "no blockers")
PYEOF
}

if [[ "${1:-}" == "--watch" ]]; then
  while true; do clear; report; sleep 30; done
else
  report
fi
