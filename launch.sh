#!/usr/bin/env bash
# WHITTLE build launcher. Spec Section 0.3.
#
# Permissions for the unattended run come from .claude/settings.json in this
# folder. Do not add a blanket permission-skip flag: the allow/deny lists are
# the mechanism, and the deny list is load-bearing (no sudo, no force-push, no
# history rewriting, no repo deletion).
set -euo pipefail
ROOT="/Users/dellboy/Documents/Vibe_Coding/WHITTLE"
cd "$ROOT"

export PYTORCH_ENABLE_MPS_FALLBACK=1

command -v tmux >/dev/null 2>&1 || { echo "tmux not found: brew install tmux" >&2; exit 1; }

if tmux has-session -t whittle 2>/dev/null; then
  echo "session 'whittle' already exists. attach with: tmux attach -t whittle"
  exit 0
fi

PROMPT='Read WHITTLE_BUILD_SPEC.md and CLAUDE.md, then execute the spec end to end under Section 0. Begin at Phase 1.'

tmux new-session -d -s whittle -c "$ROOT"
# caffeinate wraps the whole build: Tier 2 runs for days (spec 0.1 item 10)
tmux send-keys -t whittle "caffeinate -dimsu claude \"\$PROMPT\"" C-m

cat <<'MSG'
started. attach with: tmux attach -t whittle

  tmux attach -t whittle      follow the build
  ./tools/status.sh           queue progress once Phase 3 starts
  cat PROGRESS.md             phase checklist
  cat BLOCKED.md              only exists if the build needs you
MSG
