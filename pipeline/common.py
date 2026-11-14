"""Shared plumbing: paths, config, resumable manifests, rate-limited HTTP, logging.

Adapted from debris/pipeline/common.py (see REUSE.md). The manifest pattern is
the load-bearing part: every stage writes data/manifests/<stage>.jsonl and reads
it back on restart, so a killed run resumes where it stopped (spec 0.1 item 5).
Tier 2 depends on this working.

DEBRIS's corpus filter, vocabulary overlay and wandb helpers are not carried
over; nothing here knows about papers.
"""
from __future__ import annotations

import datetime as dt
import json
import os
import sqlite3
import threading
import time
import tomllib
from pathlib import Path
from typing import Any, Iterable, Iterator

ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "config"
DATA = ROOT / "data"
SHORTLIST = ROOT / "shortlist"
MANIFESTS = DATA / "manifests"
ANTIGENS = DATA / "antigens"
STRUCTURES = DATA / "structures"
RUNS = DATA / "runs"
MSA_CACHE = DATA / "msa_cache"
MODELS = ROOT / "models"
DB_PATH = DATA / "whittle.sqlite"
STATIC_DATA = ROOT / "app" / "static" / "data"

for _d in (MANIFESTS, ANTIGENS, STRUCTURES, RUNS, MSA_CACHE, MODELS, STATIC_DATA):
    _d.mkdir(parents=True, exist_ok=True)

USER_AGENT = "WHITTLE/1.0 (+https://whittle.mdeller.com; mailto:marc@marcdeller.com)"

TARGETS = ["EGFR", "GLP1R", "PRLR", "TL1A", "TNFa"]
# Cheapest first, so a broken setup surfaces on a 241 aa complex (spec 6.4 item 2).
TIER2_ORDER = ["GLP1R", "PRLR", "TNFa", "TL1A", "EGFR"]

AXES = ["affinity", "specificity", "thermo", "expression"]


# --------------------------------------------------------------------------- config

def load_toml(name: str) -> dict:
    path = CONFIG / name
    if not path.exists():
        return {}
    return tomllib.loads(path.read_text())


def load_yaml(name: str) -> dict:
    import yaml

    path = CONFIG / name
    if not path.exists():
        return {}
    return yaml.safe_load(path.read_text()) or {}


_targets_cache: dict | None = None


def targets_config() -> dict:
    """config/targets.yaml, the authoritative run configuration (spec 5.2)."""
    global _targets_cache
    if _targets_cache is None:
        _targets_cache = load_yaml("targets.yaml")
    return _targets_cache


def tuning() -> dict:
    """config/hardware.toml [tuning]. Every worker count and batch size (spec 2.2)."""
    return load_toml("hardware.toml").get("tuning", {})


def hardware() -> dict:
    return load_toml("hardware.toml").get("hardware", {})


def run_id() -> str:
    return str(targets_config().get("run", {}).get("id", "whittle-dev"))


def git_sha() -> str:
    import subprocess

    try:
        return subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=ROOT,
                              capture_output=True, text=True, check=True).stdout.strip()
    except Exception:
        return ""


def env_binary(name: str, env: str = "default") -> str | None:
    """A tool from this project's pixi environment, not whatever is on PATH.

    Scripts invoked as `.pixi/envs/default/bin/python -m pipeline.x` do not get
    the environment's bin directory on PATH, so `shutil.which` returns None and
    a stage quietly degrades. That is how the control set lost its specificity
    axis: TNP was installed and working, and `which` could not see it.
    """
    import shutil

    local = ROOT / ".pixi" / "envs" / env / "bin" / name
    if local.exists():
        return str(local)
    return shutil.which(name)


def env_with_bin(env: str = "default") -> dict:
    """os.environ with this project's environment bin directory on PATH.

    Resolving the binary is not enough on its own: TNP is a console script that
    shells out to DSSP and ANARCI, and a subprocess inherits the PATH it was
    given. Running a stage as `.pixi/envs/default/bin/python -m pipeline.x`
    leaves that directory off PATH, so the tool starts and then cannot find its
    own dependencies.
    """
    import os

    environ = dict(os.environ)
    bin_dir = str(ROOT / ".pixi" / "envs" / env / "bin")
    if bin_dir not in environ.get("PATH", "").split(os.pathsep):
        environ["PATH"] = bin_dir + os.pathsep + environ.get("PATH", "")
    return environ


def git_hash_or_blank() -> str:
    """Alias kept for the bundle builder's import list."""
    return git_sha()


def now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")


def log(message: str) -> None:
    """One timestamped line to BUILD_LOG.md and to stdout (spec 0.1 item 6).

    WHITTLE_LOG repoints the file so a test run does not write the real build log.
    """
    stamp = dt.datetime.now().strftime("%Y-%m-%d %H:%M")
    line = f"- {stamp} {message}"
    print(line, flush=True)
    target = Path(os.environ.get("WHITTLE_LOG", ROOT / "BUILD_LOG.md"))
    with target.open("a") as fh:
        fh.write(line + "\n")


# ------------------------------------------------------------------------ manifests

def manifest_path(stage: str) -> Path:
    return MANIFESTS / f"{stage}.jsonl"


def manifest_write(stage: str, record: dict) -> None:
    """Append one row and flush. A kill between rows loses at most the row in flight."""
    record = {"at": now(), **record}
    with manifest_path(stage).open("a") as fh:
        fh.write(json.dumps(record, separators=(",", ":"), default=str) + "\n")
        fh.flush()
        os.fsync(fh.fileno())


def manifest_read(stage: str) -> list[dict]:
    path = manifest_path(stage)
    if not path.exists():
        return []
    out: list[dict] = []
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        try:
            out.append(json.loads(line))
        except json.JSONDecodeError:
            continue  # a row torn by a kill mid-write: skip it, the item redoes
    return out


def manifest_done(stage: str, key: str = "id", ok_only: bool = True) -> set[str]:
    """Keys already finished, so a restart skips them (spec 0.1 item 5).

    ok_only keeps a failed item out of the done set on the first restart, so it
    gets one retry; after that its row says `status` starts with "failed:" and
    the caller decides. Items marked "skipped" are never retried.
    """
    done: set[str] = set()
    for rec in manifest_read(stage):
        if key not in rec:
            continue
        status = str(rec.get("status", "ok"))
        if ok_only and status.startswith("failed:"):
            continue
        done.add(str(rec[key]))
    return done


def manifest_failed(stage: str, key: str = "id") -> dict[str, str]:
    return {str(r[key]): str(r.get("status", "")) for r in manifest_read(stage)
            if key in r and str(r.get("status", "")).startswith("failed:")}


def jsonl_write(path: Path, records: Iterable[dict], mode: str = "w") -> int:
    n = 0
    with Path(path).open(mode) as fh:
        for rec in records:
            fh.write(json.dumps(rec, separators=(",", ":"), ensure_ascii=False,
                                default=str) + "\n")
            n += 1
    return n


def jsonl_read(path: Path) -> Iterator[dict]:
    p = Path(path)
    if not p.exists():
        return
    with p.open() as fh:
        for line in fh:
            if line.strip():
                yield json.loads(line)


def write_json(path: Path, payload: Any) -> None:
    """Atomic write, so a reader never sees half a bundle."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(payload, separators=(",", ":"), ensure_ascii=False,
                              default=str))
    tmp.replace(path)


# ----------------------------------------------------------------------------- http

class RateLimiter:
    """Process-wide token spacing. Crossref asks for politeness, UniProt for less."""

    def __init__(self, per_second: float):
        self.interval = 1.0 / max(per_second, 0.05)
        self._lock = threading.Lock()
        self._next = 0.0

    def wait(self) -> None:
        with self._lock:
            wait_for = self._next - time.monotonic()
            if wait_for > 0:
                time.sleep(wait_for)
            self._next = max(time.monotonic(), self._next) + self.interval


class TransientHTTP(Exception):
    pass


def http_get(url: str, params: dict | None = None, limiter: RateLimiter | None = None,
             timeout: int = 60, headers: dict | None = None):
    import requests
    from tenacity import (retry, retry_if_exception_type, stop_after_attempt,
                          wait_exponential)

    @retry(retry=retry_if_exception_type((TransientHTTP, requests.RequestException)),
           wait=wait_exponential(multiplier=1, min=1, max=30),
           stop=stop_after_attempt(5), reraise=True)
    def _get():
        if limiter:
            limiter.wait()
        hdrs = {"User-Agent": USER_AGENT}
        if headers:
            hdrs.update(headers)
        resp = requests.get(url, params=params, timeout=timeout, headers=hdrs)
        if resp.status_code in (429, 500, 502, 503, 504):
            raise TransientHTTP(f"{resp.status_code} from {url}")
        return resp

    return _get()


# --------------------------------------------------------------------------- sqlite

def connect(path: Path | None = None, read_only: bool = False) -> sqlite3.Connection:
    target = Path(path or DB_PATH)
    if read_only:
        conn = sqlite3.connect(f"file:{target}?mode=ro", uri=True)
    else:
        target.parent.mkdir(parents=True, exist_ok=True)
        conn = sqlite3.connect(target)
        conn.execute("PRAGMA journal_mode=WAL")
        conn.execute("PRAGMA synchronous=NORMAL")
    conn.row_factory = sqlite3.Row
    return conn


# -------------------------------------------------------------------------- blocked

def write_blocked(blocker: str, needed: str, command: str, done: str,
                  next_steps: str) -> None:
    """Spec 0.2 blocker protocol. Writes BLOCKED.md and posts a notification."""
    import subprocess

    (ROOT / "BLOCKED.md").write_text(
        f"# BLOCKED: {blocker}\n\n"
        f"## What is needed\n\n{needed}\n\n"
        f"## The command that unblocks it\n\n```bash\n{command}\n```\n\n"
        f"## Already done\n\n{done}\n\n"
        f"## What happens next\n\n{next_steps}\n"
    )
    log(f"BLOCKED {blocker}: {needed}")
    try:
        subprocess.run(["osascript", "-e",
                        f'display notification "{needed[:120]}" with title '
                        f'"WHITTLE blocked {blocker}" sound name "Glass"'],
                       capture_output=True, timeout=10)
    except Exception:
        pass


def clear_blocked(blocker: str) -> None:
    path = ROOT / "BLOCKED.md"
    if path.exists() and blocker in path.read_text():
        path.unlink()
