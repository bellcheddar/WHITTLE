# WHITTLE

- `WHITTLE_BUILD_SPEC.md` (v1.0) is authoritative. Execute it end to end under Section 0.
- `DESIGN.md` is authoritative for the app's layout and visual detail.
- **No approval gates.** Decide, log one line to `DECISIONS.md`, continue. Stop only for B1 (credentials), B2 (disk) or B3 (a control failure that would invert a headline ranking).
- Repo creation, commits, pushes and deploys to whittle.mdeller.com are all pre-approved. Run `deploy/deploy.sh` at the end of every phase.
- **The repo is PUBLIC. Never commit a credential.** Spec Section 0.4 is absolute: `.gitignore` before `git init`, `config/deploy.toml` gitignored, droplet SSH target read from the sibling `AlphaFraud/.env` at deploy time and never copied into this project. `tools/secret_scan.py` must pass before every push.
- Reuse first: copy from siblings under `/Users/dellboy/Documents/Vibe_Coding/` (BoltzMaker for the co-folding campaign, AlphaFraud and debris for droplet deploy, CODSWALLOP and GOBSMACKED for Flask tabs and Mol*, KINFOLK for hwprobe and pixi, NOBBLE for queue and cross-validation). Record every lift in `REUSE.md`. Never import across folders, never modify a sibling.
- Keep `BUILD_LOG.md`, `PROGRESS.md` and `data/manifests/` current; every stage is resumable. Tier 2 runs for days and will be interrupted.
- All compute on the Mac Studio (M1 Ultra, 128 GB). The droplet is a read-only display layer: no torch, no ML stack, no request-time compute.
- M1 has no hardware bf16: **float32 everywhere**, `PYTORCH_ENABLE_MPS_FALLBACK=1`, miniforge native arm64. Never let `bf16-mixed` back in.
- `config/targets.yaml` and `config/hardware.toml` drive everything. Nothing hard-coded.
- Rankings are **within-target percentiles**. There is never a pooled cross-target table.
- Do not install NanoMelt, CamSol, PyRosetta or upstream `boltz`. Spec 4.2 says why.
- British English. No em dashes. App copy says "predicted" or "estimated", never "confirmed" or "validated".
