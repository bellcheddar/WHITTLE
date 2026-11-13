# WHITTLE reuse log

Every lift from a sibling project under `/Users/dellboy/Documents/Vibe_Coding/`.
Copied and adapted, never imported; no sibling was modified.

| Source | Destination | What changed |
|---|---|---|
| `debris/deploy/deploy.sh` | `deploy/deploy.sh` | App name, paths, port, service name; the `--data-only` swap keeps its shape and also moves the gzipped structures |
| `debris/deploy/gunicorn.conf.py` | `deploy/gunicorn.conf.py` | Bind address, worker count from `config/deploy.toml` |
| `debris/deploy/nginx.debris.conf` | `deploy/nginx.whittle.conf` | Server name, port, `gzip_static` for `.pdb.gz` (spec 8.4), static alias |
| `debris/deploy/provision.sh` | `deploy/provision.sh` | Service user, app dir, no weekly timer (WHITTLE has no refresh job) |
| `debris/deploy/debris-gunicorn.service` | `deploy/whittle-web.service` | Names, paths, port |
| `AlphaFraud/.env` | `config/deploy.toml` (gitignored) | Read only at deploy time for the SSH target; never copied, never committed (spec 0.4 item 3) |
| `KINFOLK/tools/hwprobe.py` | `tools/hwprobe.py` | One output file instead of two; chip-family check from spec 2.1; float32 pin recorded as a spec decision rather than a probe result; minimal TOML writer so the probe runs before the env exists |
| `KINFOLK/pixi.toml`, `BoltzMaker/pixi.toml` | `pixi.toml` | Task table from spec 4.1; the two-environment split and the numpy<2 / python 3.11 ceiling are BoltzMaker's findings |
| `debris/pipeline/common.py` | `pipeline/common.py` | Paths, manifests, `log`, `RateLimiter`, tenacity `http_get`; DEBRIS-specific corpus and wandb helpers dropped |
| `debris/pipeline/schema.py` | `pipeline/schema.py` | Whole schema replaced with spec Section 7; the read-only connect and the build-then-swap pattern kept |
| `debris/tools/smoke_test.py` | `tools/smoke_test.py` | Routes and assertions for WHITTLE's three tabs |
| `debris/pipeline/verify_doi.py` | `pipeline/references.py` | Crossref title matching kept; `repo`+`commit` and database snapshot identifiers added (spec 8.3) |
| `NOBBLE/tools/vacuous_tests.py` | `tools/vacuous_tests.py` | Four WHITTLE axes instead of NOBBLE's composite |
| `NOBBLE/tools/run_queue.sh`, `status.sh` | `tools/run_queue.sh`, `tools/status.sh` | Tier 2 co-fold queue; per-complex manifest rows |
| `NOBBLE/tools/space_check.py` | `tools/space_check.py` | 40 GB floor from spec blocker B2 |
| `debris/app/app.py`, `debris/app/templates/base.html`, `GOBSMACKED/app/*` | `app/*` | Three tabs, target selector, read-only SQLite per request, static JSON bundles, `?v=<mtime>` cache busting; DEBRIS design system replaced by DESIGN.md |
| `debris/app/static/vendor/tabulator.min.{js,css}` | `app/static/vendor/` | Unchanged (vendored, gitignored, rsynced to the droplet) |
| `debris/deploy/fetch_fonts.py` | `tools/fetch_fonts.py` | Inter instead of Montserrat and IBM Plex Mono |
| `debris/app/static/vendor/plotly-2.35.2.min.js` | not copied | Rejected: the full dist is 4.6 MB and has no `parcoords` saving. Built a partial bundle from `plotly.js/lib/core` plus scatter, bar and parcoords instead (`tools/bundle/`), 1.2 MB |
| `debris/wsgi.py` | `wsgi.py` | App import path only |
| Protein-Sol reference package (protein-sol.manchester.ac.uk, GPL-3, Warwicker and Hebditch 2017) | `pipeline/solubility.py` | Not copied. The perl was read for its documented algorithm and its published fitting statistics; the Python is independent and is validated against the package's own worked example in `pipeline/tests/test_solubility.py` |
| SWI reference implementation (`Gardner-BinfLab/SoDoPE_paper_2020/SWI/swi.py`) | `pipeline/solubility.py` | Not copied. The 20 published weights and the two logistic constants are used as published constants |
| SAbDab2 OpenAPI (`/api/download/all-single-domain-summary`) | `pipeline/neighbours.py` | The old `sabdab/nano/summary/all/` TSV endpoint is gone: SAbDab is now a single-page app with a documented OpenAPI. Derived values only; the database is never mirrored (spec 5.3) |
| `tmtools` (TM-align bindings) | `pipeline/neighbours.py` | Used in place of a US-align binary; same TM-align algorithm, no binary to vendor |
