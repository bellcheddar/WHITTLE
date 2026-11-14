#!/usr/bin/env bash
# Push WHITTLE from the Studio to the droplet and restart the web service.
# Deploy to whittle.mdeller.com is pre-approved (spec 0.1 item 2) and runs at
# the end of every phase. Adapted from debris/deploy/deploy.sh (see REUSE.md).
#
#   bash deploy/deploy.sh               # code, data and static bundles, restart
#   bash deploy/deploy.sh --data-only   # SQLite, JSON bundles and structures only
#   bash deploy/deploy.sh --provision   # push, then run provision.sh on the droplet
#
# The droplet is a read-only display layer: no torch, no ML stack, no
# request-time compute (spec 2). Nothing here installs a scientific dependency.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

DATA_ONLY=0
PROVISION=0
for arg in "$@"; do
  case "$arg" in
    --data-only) DATA_ONLY=1 ;;
    --provision) PROVISION=1 ;;
    *) echo "unknown flag: $arg"; exit 2 ;;
  esac
done

# The repo is public, so this runs before anything leaves the machine
# (spec 0.4 item 4). A hit aborts rather than cleaning up afterwards.
echo "==> Secret scan"
python3 tools/secret_scan.py --tree

# A SyntaxError in one bundle takes a whole page down, and the server cannot
# tell: it serves a dead static file with a 200 exactly as it serves a healthy
# one. That shipped once, on 2026-11-14, when `const label` redeclared the
# `label` parameter of `tile()` and stopped candidates.js parsing. Deploying an
# unparseable page is worse than failing to deploy, so this gates the upload.
echo "==> JavaScript parses"
python3 tools/js_check.py

# Target from config/deploy.toml, which is gitignored and whose ssh_target was
# resolved at write time from the sibling AlphaFraud/.env (spec 4.4).
read_toml() {
  python3 -c "import tomllib;print(tomllib.load(open('config/deploy.toml','rb')).get('$1',''))"
}
SSH_TARGET="$(read_toml ssh_target)"
APP_DIR="$(read_toml app_dir)"
APP_USER="$(read_toml app_user)"
SERVICE="$(read_toml service)"
SERVER_NAME="$(read_toml server_name)"
BIND_ADDR="$(read_toml bind_addr)"

if [[ -z "$SSH_TARGET" ]]; then
  echo "config/deploy.toml has no ssh_target. See spec 4.4 (blocker B1)."
  exit 1
fi

echo "==> Target ${APP_DIR} (${SERVER_NAME} -> ${BIND_ADDR})"

sync_data() {
  ssh "$SSH_TARGET" "mkdir -p ${APP_DIR}/data/structures ${APP_DIR}/app/static/data"
  if [[ -f data/whittle.sqlite ]]; then
    # Via a temp name, moved into place on the droplet, so a reader never sees a
    # half-copied database.
    rsync -az --info=stats1 --chmod=F644 \
      data/whittle.sqlite "${SSH_TARGET}:${APP_DIR}/data/whittle.sqlite.incoming"
  fi
  rsync -az --delete --info=stats1 --chmod=D755,F644 \
    app/static/data/ "${SSH_TARGET}:${APP_DIR}/app/static/data/"
  if [[ -d data/structures ]] && [[ -n "$(ls -A data/structures 2>/dev/null)" ]]; then
    # --delete so the droplet mirrors the Mac rather than accumulating. A
    # structure withdrawn here, such as a co-fold recomputed on another
    # accelerator and discarded, must stop being served; without this the
    # drawer kept offering a complex whose metrics had been removed from the
    # database, and the viewer showed a model nothing could describe.
    rsync -az --delete --info=stats1 --chmod=D755,F644 \
      data/structures/ "${SSH_TARGET}:${APP_DIR}/data/structures/"
  fi
}

if [[ "$DATA_ONLY" -eq 1 ]]; then
  echo "==> Syncing data only"
  sync_data
  ssh "$SSH_TARGET" bash -s <<REMOTE
set -euo pipefail
cd "${APP_DIR}"
if [[ -f data/whittle.sqlite.incoming ]]; then
  mv -f data/whittle.sqlite.incoming data/whittle.sqlite
  rm -f data/whittle.sqlite-wal data/whittle.sqlite-shm
fi
chown -R ${APP_USER}:${APP_USER} data app/static/data
REMOTE
  echo "==> Data swapped. No restart needed: the app opens both per request."
  exit 0
fi

echo "==> Syncing code and static assets"
rsync -az --delete --info=stats1 --chmod=D755,F644 \
  --exclude '.venv/' --exclude '.pixi/' --exclude 'data/' --exclude '__pycache__/' \
  --exclude '*.pyc' --exclude '.git/' --exclude '.env' --exclude 'node_modules/' \
  --exclude 'models/' --exclude 'shortlist/' --exclude 'config/deploy.toml' \
  --exclude 'tests/qc_screenshots/' --exclude '.DS_Store' \
  ./ "${SSH_TARGET}:${APP_DIR}/"

echo "==> Syncing the database, the bundles and the structures"
sync_data

# The droplet's .env carries serving values only. No credential is deployed.
ssh "$SSH_TARGET" bash -s <<REMOTE
set -euo pipefail
cd "${APP_DIR}"
printf 'SERVER_NAME=%s\nBIND_ADDR=%s\nWHITTLE_DB=%s/data/whittle.sqlite\nWHITTLE_INDEXABLE=%s\n' \
  "${SERVER_NAME}" "${BIND_ADDR}" "${APP_DIR}" "${WHITTLE_INDEXABLE:-0}" > .env
if [[ -f data/whittle.sqlite.incoming ]]; then
  mv -f data/whittle.sqlite.incoming data/whittle.sqlite
  rm -f data/whittle.sqlite-wal data/whittle.sqlite-shm
fi
REMOTE

if [[ "$PROVISION" -eq 1 ]]; then
  echo "==> Provisioning on the droplet"
  ssh "$SSH_TARGET" "SERVER_NAME=${SERVER_NAME} BIND_ADDR=${BIND_ADDR} bash ${APP_DIR}/deploy/provision.sh"
else
  echo "==> Installing dependencies and restarting ${SERVICE}"
  ssh "$SSH_TARGET" bash -s <<REMOTE
set -euo pipefail
cd "${APP_DIR}"
if [[ ! -x .venv/bin/python ]]; then
  echo "No venv yet. Re-run with --provision."; exit 1
fi
sudo -u ${APP_USER} env PIP_NO_CACHE_DIR=1 ./.venv/bin/pip install --quiet -r requirements.txt
# rsync leaves new files root-owned; chown the code but prune the venv and the
# live database so an open reader is never disturbed.
find "${APP_DIR}" -path "${APP_DIR}/.venv" -prune -o -name 'whittle.sqlite*' -prune -o \
  -exec chown ${APP_USER}:${APP_USER} {} +
systemctl restart ${SERVICE}
sleep 2
systemctl --no-pager --lines=3 status ${SERVICE} || true
REMOTE
fi

# Spec Section 2: the droplet is a read-only display layer. If a scientific
# dependency ever lands there, something has quietly started computing at
# request time and the whole placement decision has been undone. Assert it
# rather than trust it.
echo "==> Asserting the droplet carries no ML stack"
ssh "$SSH_TARGET" bash -s <<'REMOTE'
set -euo pipefail
installed=$(/opt/whittle/.venv/bin/pip list --format=freeze 2>/dev/null | cut -d= -f1 | tr 'A-Z' 'a-z')
forbidden=""
for pkg in torch numpy scipy scikit-learn pandas openmm biopython fair-esm            immunebuilder anarci boltz boltz-community transformers ablang2; do
  if echo "$installed" | grep -qx "$pkg"; then forbidden="$forbidden $pkg"; fi
done
if [ -n "$forbidden" ]; then
  echo "FAIL: the droplet venv carries:$forbidden"
  echo "Spec Section 2 says the droplet does no request-time compute."
  exit 1
fi
echo "    clean: $(echo "$installed" | wc -l | tr -d ' ') packages, none of them a scientific dependency"
REMOTE

echo "==> Checking the live site"
curl -fsS --max-time 20 "https://${SERVER_NAME}/healthz" || \
  curl -fsS --max-time 20 "http://${SERVER_NAME}/healthz" || echo "healthz not answering yet"
echo
echo "==> Deployed to https://${SERVER_NAME}"
