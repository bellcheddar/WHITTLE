#!/usr/bin/env bash
# Co-fold the scrambled controls, then weight the affinity axis by what they
# show. Runs after the candidate queue, because it wants the same GPU.
#
#   nohup ./tools/calibrate_run.sh > data/calibrate.log 2>&1 &
#
# Why this exists: the affinity axis blends seven interface metrics whose
# weights are currently set by argument. A metric computed on a predicted pose
# can look entirely reasonable and carry no information about binding, and the
# only way to tell the difference is to run sequences that have no reason to
# bind through exactly the same code. That is what the scrambled controls are.
#
# Cost is about 10 complexes per target. At the rates this run measured that
# is a couple of hours, against a candidate queue of roughly ten.
set -uo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
PY="$ROOT/.pixi/envs/default/bin/python"
COFOLD_PY="$ROOT/.pixi/envs/cofold/bin/python"
PER_TARGET="${PER_TARGET:-10}"

stamp() { date "+%Y-%m-%d %H:%M:%S"; }

if [[ -f "$ROOT/tools/tier2.stop" ]]; then
  echo "$(stamp) stop file present, not starting the control co-folds"
  exit 0
fi

echo "$(stamp) co-folding $PER_TARGET scrambled controls per target"
# The same restart discipline as the candidate queue: this also runs for hours
# and is equally resumable, since every complex writes a manifest row.
attempt=0
while [[ $attempt -lt 5 ]]; do
  attempt=$((attempt + 1))
  caffeinate -dimsu "$COFOLD_PY" -m pipeline.cofold --controls "$PER_TARGET"
  status=$?
  [[ $status -eq 0 ]] && break
  echo "$(stamp) control queue exited $status, attempt $attempt, retrying"
  sleep $((attempt * 30))
done

echo "$(stamp) interface metrics over the controls"
"$PY" -m pipeline.interface
"$PY" -m pipeline.interactions
"$PY" -m pipeline.complementarity

echo "$(stamp) how well does each metric separate real from shuffled?"
"$PY" tools/calibrate_affinity.py

echo "$(stamp) done. The weights above are a recommendation; applying them is"
echo "$(stamp) a deliberate edit to pipeline/score.py so the reasoning is"
echo "$(stamp) recorded in the commit rather than applied silently."
