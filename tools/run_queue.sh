#!/usr/bin/env bash
# Tier 2 co-folding queue. Spec 6.4 item 7.
# Adapted from NOBBLE/tools/run_queue.sh (see REUSE.md).
#
#   ./tools/run_queue.sh --smoke      # one GLP1R complex, timed, first
#   ./tools/run_queue.sh              # the full queue, cheapest target first
#   ./tools/status.sh --watch         # progress, from another shell
#
# Resumable by construction: every complex writes a manifest row and a restart
# skips what is done. This runs for days and WILL be interrupted, so it is
# wrapped in caffeinate and it never assumes it is starting from nothing.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

LOG="${WHITTLE_QUEUE_LOG:-$ROOT/data/tier2_queue.log}"
mkdir -p "$(dirname "$LOG")"

echo "==> Disk headroom (spec blocker B2)"
"$ROOT/.pixi/envs/default/bin/python" tools/space_check.py --quiet || {
  echo "Not enough disk for Tier 2. See BLOCKED.md."; exit 1; }

echo "==> Tier 2 queue starting. Log: $LOG"
echo "    Cheapest target first: GLP1R, PRLR, TNFa, TL1A, EGFR"

# caffeinate: Tier 2 runs for days and a sleeping Mac is an interrupted queue
# (spec 0.1 item 10).
exec caffeinate -dimsu "$ROOT/.pixi/envs/cofold/bin/python" \
  -m pipeline.cofold "$@" 2>&1 | tee -a "$LOG"
