"""Hardware probe. Spec Sections 2.1 and 2.2.

Writes config/hardware.toml. Every batch size and worker count downstream is
read from it; nothing is hard-coded (spec 2.2 last bullet).

Adapted from KINFOLK/tools/hwprobe.py (see REUSE.md). Changed here: one file
rather than hardware.toml plus tuning.toml, the chip-family check that spec 2.1
asks for, and the float32 pin, which is a spec decision rather than a probe
result and so is recorded as such.

    python tools/hwprobe.py
"""
from __future__ import annotations

import datetime as dt
import platform
import re
import shutil
import subprocess
import sys
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "config"


def sh(cmd: list[str]) -> str:
    try:
        return subprocess.run(cmd, capture_output=True, text=True, check=True).stdout.strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return ""


def gpu_cores() -> int:
    out = sh(["system_profiler", "SPDisplaysDataType"])
    match = re.search(r"Total Number of Cores:\s*(\d+)", out)
    return int(match.group(1)) if match else 0


def mps_available() -> bool:
    try:
        import torch
    except Exception:
        return False
    try:
        return bool(torch.backends.mps.is_available())
    except Exception:
        return False


def torch_version() -> str:
    try:
        import torch

        return str(torch.__version__)
    except Exception:
        return ""


def probe() -> dict:
    usage = shutil.disk_usage(ROOT)
    chip = sh(["sysctl", "-n", "machdep.cpu.brand_string"])
    model = ""
    for line in sh(["system_profiler", "SPHardwareDataType"]).splitlines():
        if "Model Name" in line:
            model = line.split(":", 1)[1].strip()
    family = ""
    if match := re.search(r"Apple (M\d)", chip):
        family = match.group(1)
    return {
        "probed_at": dt.datetime.now().isoformat(timespec="seconds"),
        "hostname": platform.node(),
        "model": model,
        "chip": chip,
        "chip_family": family,
        "is_mac_studio": "Studio" in model,
        "performance_cores": int(sh(["sysctl", "-n", "hw.perflevel0.physicalcpu"]) or 0),
        "efficiency_cores": int(sh(["sysctl", "-n", "hw.perflevel1.physicalcpu"]) or 0),
        "logical_cpus": int(sh(["sysctl", "-n", "hw.ncpu"]) or 0),
        "unified_memory_gb": round(int(sh(["sysctl", "-n", "hw.memsize"]) or 0) / 2**30, 1),
        "gpu_cores": gpu_cores(),
        "free_disk_gb": round(usage.free / 2**30, 1),
        "mps_available": mps_available(),
        "torch_version": torch_version(),
        "macos": platform.mac_ver()[0],
    }


def derive(hw: dict, previous: dict | None) -> dict:
    perf = hw["performance_cores"] or 8
    cpu_workers = max(2, perf - 2)
    tuning = {
        # Leave two performance cores for the folding process and the OS.
        "cpu_workers": cpu_workers,
        "io_workers": 12,
        "fold_device": "cpu",  # ImmuneBuilder is small; OpenMM refinement is CPU-only
        "fold_batch": 1,
        "openmm_platform": "CPU",  # spec 2.2: at ~125 residues Metal buys nothing
        "esm_device": "mps" if hw["mps_available"] else "cpu",
        "esm_batch": 8 if hw["unified_memory_gb"] >= 64 else 2,
        "esm_model": "esm2_t33_650M_UR50D" if hw["unified_memory_gb"] >= 32 else "esm2_t30_150M_UR50D",
        # float32 everywhere. Spec 2.2 and CLAUDE.md: a spec decision, kept on
        # M2 as well so the Tier 1b model and the Tier 2 co-folds are numerically
        # the same run whichever Studio they land on.
        "precision": "float32",
        "precision_reason": "spec 2.2: float32 everywhere; bf16-mixed never re-enabled",
        "cofold_accelerator": "mps" if hw["mps_available"] else "cpu",
        "cofold_seeds": 5,
        "cofold_recycling": 3,
        "memory_headroom_gb": 24,
        "max_tier2_hours": 96,  # spec 5.2: beyond four days, scope EGFR to domain III
        "env": {
            "PYTORCH_ENABLE_MPS_FALLBACK": "1",
            "TOKENIZERS_PARALLELISM": "false",
            "OMP_NUM_THREADS": str(cpu_workers),
            "KMP_DUPLICATE_LIB_OK": "TRUE",
        },
    }
    if previous:
        # Anything a smoke test has already settled empirically wins over a guess.
        for key in ("fold_device", "esm_device", "esm_batch", "esm_model",
                    "cofold_accelerator", "cofold_seeds", "cofold_recycling",
                    "max_tier2_hours"):
            if key in previous and previous.get(f"{key}_status") == "measured":
                tuning[key] = previous[key]
    return tuning


def dumps(table: dict, name: str) -> str:
    """Minimal TOML writer, so the probe does not need tomli_w to run."""
    lines = [f"# Written by tools/hwprobe.py. Do not edit by hand.", f"[{name}]"]
    nested: list[tuple[str, dict]] = []
    for key, value in table.items():
        if isinstance(value, dict):
            nested.append((key, value))
            continue
        if isinstance(value, bool):
            lines.append(f"{key} = {str(value).lower()}")
        elif isinstance(value, (int, float)):
            lines.append(f"{key} = {value}")
        else:
            lines.append(f'{key} = "{value}"')
    for key, value in nested:
        lines.append("")
        lines.append(f"[{name}.{key}]")
        for k2, v2 in value.items():
            lines.append(f'{k2} = "{v2}"')
    return "\n".join(lines) + "\n"


def main() -> int:
    CONFIG.mkdir(exist_ok=True)
    hw = probe()
    path = CONFIG / "hardware.toml"
    previous = None
    if path.exists():
        try:
            previous = tomllib.loads(path.read_text()).get("tuning")
        except tomllib.TOMLDecodeError:
            previous = None
    tuning = derive(hw, previous)
    path.write_text(dumps(hw, "hardware") + "\n" + dumps(tuning, "tuning"))

    notes: list[str] = []
    if not hw["is_mac_studio"]:
        notes.append("not a Mac Studio: spec 2.1 says do not run Tier 1 or Tier 2 here")
    if hw["chip_family"] and hw["chip_family"] != "M1":
        notes.append(
            f"chip is {hw['chip']}, spec 2.2 is written for M1 Ultra. float32 kept "
            f"anyway (spec 2.2 and CLAUDE.md are authoritative)"
        )
    if hw["free_disk_gb"] < 40:
        notes.append(f"free disk {hw['free_disk_gb']} GB is under the 40 GB floor: blocker B2")
    if not hw["mps_available"]:
        notes.append("no MPS: ESM-2 and co-folding fall back to CPU")

    print(path.read_text())
    for note in notes:
        print(f"NOTE: {note}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
