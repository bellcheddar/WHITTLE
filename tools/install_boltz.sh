#!/usr/bin/env bash
# Install boltz-community into the `cofold` pixi environment. Spec 4.2 and 4.3.
#
#   pixi run -e cofold install-boltz
#
# Why this is a task rather than a `[pypi-dependencies]` entry: spec 4.2
# forbids upstream `boltz`, which declares a CUDA-only `cuequivariance` and
# will not install on macOS, and pixi's PyPI resolver reaches for packages by
# name. Installing the community fork explicitly, with no extras, is the only
# way to be certain which one landed. Phase 3 verifies the result before
# queueing anything (spec 4.3 last line).
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

PY="${PY:-$ROOT/.pixi/envs/cofold/bin/python}"
if [[ ! -x "$PY" ]]; then
  echo "No cofold environment at $PY. Run: pixi install -e cofold" >&2
  exit 1
fi

echo "==> Installing boltz-community (the Apple Silicon MPS fork), no [cuda] extra"
"$PY" -m pip install --no-input --upgrade "boltz-community"

echo "==> Verifying which boltz landed"
"$PY" - <<'PYEOF'
import importlib.metadata as md
import shutil
import sys

names = {d.metadata["Name"].lower() for d in md.distributions()
         if d.metadata and d.metadata.get("Name")}
if "boltz" in names and "boltz-community" not in names:
    sys.exit("FAIL: upstream `boltz` is installed, not the community fork. "
             "Spec 4.2 forbids it.")
if "boltz-community" not in names:
    sys.exit("FAIL: boltz-community is not installed.")
version = md.version("boltz-community")
print(f"boltz-community {version}")
try:
    import boltz
    print("import boltz ok from", boltz.__file__)
except Exception as exc:
    print(f"WARNING: boltz-community is installed but `import boltz` failed: {exc}")
binary = shutil.which("boltz")
print("boltz CLI:", binary or "not on PATH")
PYEOF

echo '==> Done. Check the accelerator with:'
echo "    pixi run -e cofold python -c 'import torch; print(torch.backends.mps.is_available())'"
