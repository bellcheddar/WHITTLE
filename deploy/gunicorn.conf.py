"""gunicorn config for WHITTLE. Adapted from debris/deploy/gunicorn.conf.py.

The app is read-only: it serves static JSON bundles and one read-only SQLite
query per drawer open, so sync-ish gthread workers are right and no worker holds
state. Two workers, because the droplet runs a dozen other apps and this one has
no request-time computation to parallelise (spec 2: the droplet is a display
layer).
"""
import os

bind = os.environ.get("BIND_ADDR", "127.0.0.1:8021")
workers = int(os.environ.get("WEB_WORKERS", "2"))
threads = 4
worker_class = "gthread"
timeout = 60
graceful_timeout = 30
keepalive = 5
accesslog = "-"
errorlog = "-"
loglevel = "info"
# A data swap replaces the SQLite file and the JSON bundles underneath a running
# worker. Both are opened per request, so nothing needs a restart.
preload_app = False
