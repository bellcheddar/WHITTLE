#!/usr/bin/env bash
# Supervisor for the Tier 2 co-folding queue. Spec 0.1 item 5 and 6.4 item 7.
#
#   nohup ./tools/supervise_tier2.sh > data/tier2_supervisor.log 2>&1 &
#
# The queue is resumable by construction: every complex writes a manifest row
# and a restart skips what is done. That makes a crash recoverable, but only if
# something is there to do the restarting, and a multi-hour job left alone at
# 3am is exactly when a stray OOM or a dropped MPS context would otherwise end
# the run quietly.
#
# This loop restarts the queue whenever it exits non-zero, backing off a little
# each time so a hard failure does not spin. It stops for three reasons, all of
# them deliberate:
#   * the queue reports an empty queue, which is what completion looks like;
#   * MAX_RESTARTS consecutive failures with no progress, so a broken setup
#     does not retry for ever;
#   * tools/tier2.stop exists, which is the manual brake.
set -uo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

PY="$ROOT/.pixi/envs/cofold/bin/python"
MANIFEST="$ROOT/data/manifests/cofold.jsonl"
STOP_FILE="$ROOT/tools/tier2.stop"
MAX_RESTARTS="${MAX_RESTARTS:-20}"
EXPECTED="${EXPECTED:-250}"

stamp() { date "+%Y-%m-%d %H:%M:%S"; }
done_count() {
  [[ -f "$MANIFEST" ]] || { echo 0; return; }
  "$ROOT/.pixi/envs/default/bin/python" - <<'PYEOF'
import json
from pathlib import Path
path = Path("data/manifests/cofold.jsonl")
ok = set()
for line in path.read_text().splitlines():
    if not line.strip():
        continue
    try:
        rec = json.loads(line)
    except json.JSONDecodeError:
        continue
    if rec.get("status") == "ok" and rec.get("id"):
        ok.add(rec["id"])
print(len(ok))
PYEOF
}

echo "$(stamp) supervisor starting. Expecting $EXPECTED complexes."
echo "$(stamp) stop it cleanly with: touch tools/tier2.stop"

failures=0
attempt=0
while true; do
  if [[ -f "$STOP_FILE" ]]; then
    echo "$(stamp) stop file present, standing down with $(done_count) complexes done"
    exit 0
  fi

  before="$(done_count)"
  if [[ "$before" -ge "$EXPECTED" ]]; then
    echo "$(stamp) $before of $EXPECTED complexes done, queue complete"
    break
  fi

  attempt=$((attempt + 1))
  echo "$(stamp) attempt $attempt, $before of $EXPECTED done, starting the queue"
  caffeinate -dimsu "$PY" -m pipeline.cofold
  status=$?
  after="$(done_count)"
  echo "$(stamp) queue exited $status with $after of $EXPECTED done"

  if [[ "$status" -eq 0 ]]; then
    if [[ "$after" -ge "$EXPECTED" ]] || [[ "$after" -eq "$before" ]]; then
      # A clean exit that added nothing means the queue had nothing left.
      echo "$(stamp) clean exit and nothing left to do"
      break
    fi
    failures=0
    continue
  fi

  if [[ "$after" -gt "$before" ]]; then
    # It crashed, but it made progress first, so the fault is likely with one
    # complex rather than the setup. Reset the counter and carry on.
    echo "$(stamp) made progress before failing, resetting the failure count"
    failures=0
  else
    failures=$((failures + 1))
    echo "$(stamp) failure $failures of $MAX_RESTARTS with no progress"
  fi

  if [[ "$failures" -ge "$MAX_RESTARTS" ]]; then
    echo "$(stamp) giving up after $failures consecutive failures with no progress"
    "$ROOT/.pixi/envs/default/bin/python" - <<'PYEOF'
import sys
sys.path.insert(0, ".")
from pipeline.common import write_blocked
write_blocked(
    "Tier 2 queue",
    "The co-folding queue failed repeatedly without completing a single "
    "complex, so the supervisor stopped rather than retry for ever.",
    "tail -50 data/tier2_supervisor.log && ./tools/status.sh",
    "Everything up to and including Tier 1 is complete and deployed. The "
    "queue is resumable: whatever finished is kept.",
    "Fix the cause and restart with ./tools/supervise_tier2.sh. Nothing else "
    "in the build depends on this finishing.")
PYEOF
    exit 1
  fi

  backoff=$((failures * 30))
  echo "$(stamp) restarting in ${backoff}s"
  sleep "$backoff"
done

echo "$(stamp) queue finished with $(done_count) of $EXPECTED complexes"
echo "$(stamp) pulling the results through: interface, ranking, controls, bundles"
bash tools/refresh.sh || echo "$(stamp) refresh failed; run ./tools/refresh.sh by hand"

# The GPU is free now, so the negative control for every interface metric can
# run. Until it does, the affinity axis weights are set by argument rather
# than by evidence.
echo "$(stamp) co-folding the scrambled controls and calibrating the axis"
bash tools/calibrate_run.sh \
  || echo "$(stamp) calibration failed; run ./tools/calibrate_run.sh by hand"
echo "$(stamp) supervisor done"
