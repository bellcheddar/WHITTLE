"""Disk headroom before a long run. Spec blocker B2.

Adapted from NOBBLE/tools/space_check.py (see REUSE.md); the 40 GB floor and
the B2 protocol are WHITTLE's.

    python tools/space_check.py            # report and exit non-zero if low
    python tools/space_check.py --need 80  # a different floor, in GB
"""
from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

FLOOR_GB = 40.0


def directory_size(path: Path) -> float:
    total = 0
    for item in path.rglob("*"):
        try:
            if item.is_file():
                total += item.stat().st_size
        except OSError:
            continue
    return total / 2**30


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--need", type=float, default=FLOOR_GB)
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()

    usage = shutil.disk_usage(ROOT)
    free = usage.free / 2**30
    data = ROOT / "data"
    breakdown = {}
    if data.exists():
        for child in sorted(data.iterdir()):
            if child.is_dir():
                size = directory_size(child)
                if size > 0.01:
                    breakdown[child.name] = round(size, 2)

    if not args.quiet:
        print(f"free: {free:.1f} GB (floor {args.need:.0f} GB)")
        print(f"data/: {sum(breakdown.values()):.2f} GB")
        for name, size in sorted(breakdown.items(), key=lambda kv: -kv[1]):
            print(f"  {name:20} {size:8.2f} GB")

    if free < args.need:
        from pipeline.common import write_blocked

        prunable = ", ".join(
            f"data/{n} ({s:.1f} GB)" for n, s in
            sorted(breakdown.items(), key=lambda kv: -kv[1])[:4]) or "nothing obvious"
        write_blocked(
            "B2 Disk",
            f"Only {free:.1f} GB free on the project volume; the floor is "
            f"{args.need:.0f} GB and Tier 2 needs room for the co-fold runs.",
            "python tools/space_check.py",
            f"The largest directories under data/ are: {prunable}. "
            f"data/runs/ holds the co-fold outputs and data/sabdab/pdb/ holds "
            f"cached PDB entries, both of which regenerate.",
            "Tier 2 will not start until there is headroom. Everything that does "
            "not write to disk continues.")
        print(f"BLOCKED: {free:.1f} GB free is below the {args.need:.0f} GB floor")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
