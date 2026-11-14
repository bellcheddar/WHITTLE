#!/usr/bin/env bash
# The whole Tier 2 night as one unattended unit.
#
#   nohup ./tools/overnight.sh > data/overnight.log 2>&1 &
#
# Order matters here, and not for throughput:
#
#  1. **Scrambled controls for the two cheapest targets first.** They cost
#     about half an hour and they are the only thing that can tell us whether
#     any interface metric means anything. Running them first means that if
#     the night is interrupted we still hold the negative control for GLP1R
#     and PRLR, rather than a pile of candidate co-folds no one can calibrate.
#     It also measures the real MPS rate on the production path before the
#     long run commits to it.
#  2. The candidate queue, 50 per target, cheapest target first.
#  3. The remaining scrambled controls and the calibration, which
#     supervise_tier2.sh chains itself.
#
# Everything is resumable: each complex writes a manifest row and every stage
# skips what is already done, so this can be killed and restarted at will.
set -uo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
COFOLD_PY="$ROOT/.pixi/envs/cofold/bin/python"
EARLY_TARGETS="${EARLY_TARGETS:-GLP1R PRLR}"
EARLY_N="${EARLY_N:-10}"

stamp() { date "+%Y-%m-%d %H:%M:%S"; }

if [[ -f "$ROOT/tools/tier2.stop" ]]; then
  echo "$(stamp) stop file present, standing down"
  exit 0
fi

args=()
for t in $EARLY_TARGETS; do args+=(--target "$t"); done

echo "$(stamp) stage 1: $EARLY_N scrambled controls each for $EARLY_TARGETS"
echo "$(stamp) this also measures the real MPS rate before the long run"
caffeinate -dimsu "$COFOLD_PY" -m pipeline.cofold --controls "$EARLY_N" "${args[@]}" \
  || echo "$(stamp) early controls exited non-zero; continuing to the candidates"

echo "$(stamp) stage 2: the candidate queue under its supervisor"
bash tools/supervise_tier2.sh

echo "$(stamp) overnight run complete"
