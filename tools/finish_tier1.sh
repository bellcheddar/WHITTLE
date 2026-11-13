#!/usr/bin/env bash
# Close out Tier 1 and get the ranking into its final state, in the one order
# that works: structures, then scoring, then the controls scored against that
# scoring, then the checks that read the controls, then the bundles.
#
# Run after the nativeness stage completes. Every step is resumable on its own,
# so a failure part-way leaves the rest re-runnable.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
PY="$ROOT/.pixi/envs/default/bin/python"

echo "==> Tier 1 structures table"
"$PY" -m pipeline.tier1 --workers 10

echo "==> Scoring: within-target percentiles, Pareto fronts, desirability"
"$PY" -m pipeline.score

echo "==> Control set rebuilt against the final ranking"
"$PY" -m pipeline.controls

echo "==> Controls run blind through the same pipeline"
"$PY" -m pipeline.controls_run --workers 3

echo "==> Tm predicted for the control sequences"
"$PY" -m pipeline.tm_predict --controls

echo "==> Spec Section 9 assertions re-evaluated"
"$PY" -m pipeline.controls

echo "==> QC"
"$PY" tools/qc_pareto.py
"$PY" tools/vacuous_tests.py

echo "==> Provenance and references refreshed"
"$PY" -m pipeline.references

echo "==> Bundles and findings"
"$PY" -m pipeline.build_db
"$PY" tools/findings.py

echo "==> Tier 1 closed out."
