#!/usr/bin/env bash
# Finish the blind control run and the checks that depend on it.
# Separate from finish_tier1.sh because the first pass lost TNP to a PATH
# fault and the assertions step then deleted what had been computed.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
PY="$ROOT/.pixi/envs/default/bin/python"

echo "==> TNP over the control sequences (nativeness is already cached)"
"$PY" -m pipeline.controls_run --skip-nativeness --workers 3

echo "==> Assertions, carrying the predictions through"
"$PY" -m pipeline.controls

echo "==> Vacuous-axis tests, now with a scrambled control to compare against"
"$PY" tools/vacuous_tests.py

echo "==> Bundles and findings"
"$PY" -m pipeline.build_db
"$PY" tools/findings.py
echo "==> Controls closed out."
