#!/usr/bin/env bash
# One-time droplet provisioning for WHITTLE. Runs as root ON the droplet, after
# the code is present at /opt/whittle (push it first with deploy/deploy.sh).
#
#   sudo SERVER_NAME=whittle.mdeller.com bash /opt/whittle/deploy/provision.sh
#
# Idempotent. Installs the service user, the venv, the systemd web service, the
# nginx site and a Let's Encrypt certificate. Copied from
# debris/deploy/provision.sh with names, paths and port changed, and the weekly
# timer dropped: WHITTLE has no refresh job (see REUSE.md).
#
# The repo is public, so no contact address is written here. certbot reuses the
# ACME account already registered on this droplet; set CERTBOT_EMAIL in the
# environment only if a fresh account ever has to be created.
set -euo pipefail

APP_DIR=/opt/whittle
APP_USER=whittle
BIND_ADDR="${BIND_ADDR:-127.0.0.1:8021}"

if [[ -f "$APP_DIR/.env" ]]; then set -a; source "$APP_DIR/.env"; set +a; fi
SERVER_NAME="${SERVER_NAME:-whittle.mdeller.com}"

echo "==> WHITTLE provisioning for ${SERVER_NAME} on ${BIND_ADDR}"
if [[ $EUID -ne 0 ]]; then echo "Run as root."; exit 1; fi
if [[ ! -f "$APP_DIR/wsgi.py" ]]; then
  echo "No code at $APP_DIR. Push it first: bash deploy/deploy.sh (from the Mac)."
  exit 1
fi

echo "==> Service user '${APP_USER}'"
id -u "$APP_USER" &>/dev/null || \
  useradd --system --home "$APP_DIR" --shell /usr/sbin/nologin "$APP_USER"
mkdir -p "$APP_DIR/data/structures" "$APP_DIR/app/static/data"
chown -R "$APP_USER:$APP_USER" "$APP_DIR"
# useradd --home creates the directory 0700; nginx runs as www-data and serves
# app/static/ and data/structures/ straight from disk, so the path has to be
# traversable by others.
chmod 755 "$APP_DIR"

echo "==> Python venv"
if [[ ! -x "$APP_DIR/.venv/bin/python" ]]; then
  sudo -u "$APP_USER" python3 -m venv "$APP_DIR/.venv"
fi
sudo -u "$APP_USER" env PIP_NO_CACHE_DIR=1 "$APP_DIR/.venv/bin/pip" install --quiet --upgrade pip
sudo -u "$APP_USER" env PIP_NO_CACHE_DIR=1 "$APP_DIR/.venv/bin/pip" install --quiet -r "$APP_DIR/requirements.txt"

echo "==> systemd unit"
cp "$APP_DIR/deploy/whittle-web.service" /etc/systemd/system/
systemctl daemon-reload
systemctl enable --now whittle-web.service

echo "==> nginx site"
sed -e "s/__SERVER_NAME__/${SERVER_NAME}/" -e "s|__BIND_ADDR__|${BIND_ADDR}|" \
  "$APP_DIR/deploy/nginx.whittle.conf" > /etc/nginx/sites-available/whittle
ln -sf /etc/nginx/sites-available/whittle /etc/nginx/sites-enabled/whittle
nginx -t && systemctl reload nginx

echo "==> TLS certificate"
if [[ ! -d "/etc/letsencrypt/live/${SERVER_NAME}" ]]; then
  CERTBOT_ARGS=(--nginx -n --agree-tos --redirect -d "${SERVER_NAME}")
  if [[ -n "${CERTBOT_EMAIL:-}" ]]; then
    CERTBOT_ARGS+=(-m "${CERTBOT_EMAIL}")
  else
    CERTBOT_ARGS+=(--register-unsafely-without-email)
  fi
  certbot "${CERTBOT_ARGS[@]}" || \
    echo "certbot failed; the site is still served over HTTP"
fi
nginx -t && systemctl reload nginx

echo "==> Provisioned. ${SERVER_NAME} -> ${BIND_ADDR}"
systemctl --no-pager --lines=3 status whittle-web.service || true
