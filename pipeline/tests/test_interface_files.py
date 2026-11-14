"""Which files on disk count as a diffusion sample.

This exists for one fault that has already cost a stage. The run directories
live under iCloud-synced Documents, and a directory rewritten while a sync is
in flight leaves the superseded copy behind as `<stem>_model_4 2.pdb`. Those
are real PDBs of the same candidate from an earlier run, so a plain glob takes
them for extra samples, and parsing the index off the end of the stem then
raises `ValueError: invalid literal for int() with base 10: '4 2'`.

On 2026-11-14 that ended the first interface pass after all 250 co-folds had
finished, on the first candidate that had a duplicate. The cost was not the
crash but the hours of GPU time whose results went unread until someone looked
at the log.
"""
from __future__ import annotations

from pathlib import Path

from pipeline.interface import structures


def _sample(dir_path: Path, name: str) -> Path:
    path = dir_path / name
    path.write_text("ATOM\n")
    return path


def test_samples_are_returned_in_sample_order(tmp_path: Path) -> None:
    for n in (3, 0, 4, 1, 2):          # written out of order on purpose
        _sample(tmp_path, f"T__T-0001_model_{n}.pdb")
    got = structures(tmp_path)
    assert [p.name for p in got] == [f"T__T-0001_model_{n}.pdb" for n in range(5)]


def test_icloud_duplicates_are_skipped_not_counted(tmp_path: Path) -> None:
    for n in range(5):
        _sample(tmp_path, f"T__T-0001_model_{n}.pdb")
        _sample(tmp_path, f"T__T-0001_model_{n} 2.pdb")
    got = structures(tmp_path)
    # Five samples, not ten, and no crash parsing "4 2" as an integer.
    assert len(got) == 5
    assert all(" 2" not in p.name for p in got)


def test_pre_affinity_models_are_not_samples(tmp_path: Path) -> None:
    _sample(tmp_path, "T__T-0001_model_0.pdb")
    _sample(tmp_path, "pre_affinity_T__T-0001_model_0.pdb")
    assert [p.name for p in structures(tmp_path)] == ["T__T-0001_model_0.pdb"]


def test_nested_run_layout_is_searched(tmp_path: Path) -> None:
    nested = tmp_path / "boltz_results_T__T-0001" / "predictions" / "T__T-0001"
    nested.mkdir(parents=True)
    _sample(nested, "T__T-0001_model_0.pdb")
    assert len(structures(tmp_path)) == 1
