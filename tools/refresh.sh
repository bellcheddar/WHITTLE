#!/usr/bin/env bash
# Pull whatever Tier 2 has finished into the live site, without a full deploy.
#
#   ./tools/refresh.sh
#
# The co-folding queue runs for days. This takes the complexes that are done,
# recomputes the interface metrics and the ranking, rebuilds the bundles and
# swaps the data on the droplet. Safe to run while the queue is still going:
# every stage reads finished manifest rows and ignores the rest, and the data
# swap moves a complete file into place rather than writing in situ.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
PY="$ROOT/.pixi/envs/default/bin/python"

echo "==> Interface metrics for the co-folds that have finished"
"$PY" -m pipeline.interface

echo "==> Interaction profile across each finished interface (PLIP)"
"$PY" -m pipeline.interactions

echo "==> Shape complementarity across each finished interface"
"$PY" -m pipeline.complementarity

echo "==> Re-ranking"
"$PY" -m pipeline.score

echo "==> Control assertions, including the scrambled negative control"
"$PY" -m pipeline.controls

echo "==> QC"
"$PY" tools/qc_pareto.py || true
"$PY" tools/vacuous_tests.py || true

echo "==> Bundles and findings"
"$PY" -m pipeline.build_db
"$PY" tools/findings.py

echo "==> Swapping the data on the droplet"
bash deploy/deploy.sh --data-only

echo "==> Live site refreshed. Queue progress: ./tools/status.sh"
