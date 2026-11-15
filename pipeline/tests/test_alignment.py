"""The alignment grid. The claim being tested is that a column is one IMGT
position in every row, because the whole panel is read on that basis."""
from __future__ import annotations

import sqlite3
from pathlib import Path

import pytest

from pipeline.build_db import column_key, column_order

ROOT = Path(__file__).resolve().parents[2]
DB = ROOT / "data" / "whittle.sqlite"


def test_imgt_insertions_ascend_after_the_lower_anchor() -> None:
    keys = ["110", "111", "111A", "111B"]
    assert sorted(keys, key=column_key) == keys


def test_imgt_insertions_descend_into_the_upper_anchor() -> None:
    # 112B comes before 112A, which comes before 112 itself. A plain sort puts
    # them the other way round and silently reverses half of every long CDR3.
    keys = ["112B", "112A", "112", "113"]
    assert sorted(keys, key=column_key) == keys


def test_a_long_cdr3_orders_the_way_imgt_writes_it() -> None:
    scrambled = ["113", "112", "111A", "112A", "111", "111B", "112B", "110"]
    assert sorted(scrambled, key=column_key) == [
        "110", "111", "111A", "111B", "112B", "112A", "112", "113"]


def test_column_order_rejects_an_order_no_row_could_have() -> None:
    """The guard, not the happy path: a skeleton that contradicts the column
    order means two residues of one sequence would land in one column."""
    with pytest.raises(SystemExit) as caught:
        column_order([["112", "112A"]], "TEST")
    assert "column order disagrees" in str(caught.value)


def _rows(target: str = "EGFR"):
    if not DB.exists():
        pytest.skip("no database built")
    conn = sqlite3.connect(f"file:{DB}?mode=ro", uri=True)
    conn.row_factory = sqlite3.Row
    try:
        rows = conn.execute(
            "SELECT a.row_id AS candidate_id, a.keys, a.sse, a.n_residues, "
            "       a.residues, a.offset, c.sequence "
            "FROM alignment a JOIN candidates c ON c.id = a.row_id "
            "WHERE a.target = ? AND a.kind = 'candidate'", (target,)).fetchall()
    finally:
        conn.close()
    if not rows:
        pytest.skip("the alignment stage has not run")
    return rows


def test_every_row_is_a_subsequence_of_the_shared_column_order() -> None:
    """The load-bearing property. If one row's residues do not follow the
    column order, that row is drawn in the wrong columns and every identity
    read off the panel is wrong."""
    for target in ("EGFR", "GLP1R", "PRLR", "TL1A", "TNFa"):
        rows = _rows(target)
        skeletons = [r["keys"].split(" ") for r in rows if r["keys"]]
        columns = column_order(skeletons, target)          # raises if not
        assert len(columns) == len(set(columns))


def test_no_row_repeats_an_imgt_key() -> None:
    for row in _rows():
        keys = row["keys"].split(" ") if row["keys"] else []
        assert len(keys) == len(set(keys)), row["candidate_id"]


def test_the_residues_placed_are_a_contiguous_run_of_the_sequence() -> None:
    """ANARCI declines residues at the ends, never in the middle, so the
    placed residues are a contiguous stretch and `offset` says where it
    starts. Taking them from the wrong end would slide the whole row against
    the grid, and the motifs with it."""
    for row in _rows():
        placed, offset = row["residues"], row["offset"]
        assert placed, row["candidate_id"]
        assert len(placed) == row["n_residues"], row["candidate_id"]
        assert offset >= 0, row["candidate_id"]
        assert row["sequence"][offset:offset + len(placed)] == placed, \
            row["candidate_id"]


def test_secondary_structure_is_one_letter_per_placed_residue() -> None:
    for row in _rows():
        if not row["sse"]:
            continue
        assert len(row["sse"]) == row["n_residues"], row["candidate_id"]
        assert set(row["sse"]) <= set("HEC.-"), row["candidate_id"]


def test_the_bundle_grid_reconstructs_each_candidate_s_own_sequence() -> None:
    """End to end over the shipped file: expanding a row across the columns
    and dropping the gaps has to give back the sequence the row came from,
    or the panel is drawing something that is not this candidate."""
    import json

    conn = sqlite3.connect(f"file:{DB}?mode=ro", uri=True) if DB.exists() else None
    if conn is None:
        pytest.skip("no database built")
    conn.row_factory = sqlite3.Row
    try:
        for target in ("EGFR", "TNFa"):
            path = ROOT / "app" / "static" / "data" / f"alignment_{target}.json"
            if not path.exists():
                pytest.skip("bundles not built")
            bundle = json.loads(path.read_text())
            stored = {r["id"]: r["sequence"] for r in conn.execute(
                "SELECT id, sequence FROM candidates WHERE target = ?", (target,))}
            n_columns = len(bundle["columns"])
            for row in bundle["rows"]:
                skeleton = bundle["skeletons"][row["sk"]]
                assert len(skeleton) == len(row["seq"])
                assert all(0 <= c < n_columns for c in skeleton)
                assert skeleton == sorted(skeleton)
                # A contiguous run of the sequence, from either end: the
                # numbering declines residues at the front (an unnumberable
                # N-terminal extension) and at the back (TNFa-0500's corrupt
                # J region drops 21). Slicing the tail off the sequence, as
                # this bundle first did, slid that row 21 columns.
                assert row["seq"] in stored[row["id"]], row["id"]
    finally:
        conn.close()


def test_every_liability_lands_on_the_residue_it_names() -> None:
    """The bundle gives a liability as an index into the row's own residues.
    The motif has to be on the residue the sequence scan flagged, not one
    along, or the panel underlines the wrong amino acid."""
    import json

    if not DB.exists():
        pytest.skip("no database built")
    conn = sqlite3.connect(f"file:{DB}?mode=ro", uri=True)
    conn.row_factory = sqlite3.Row
    try:
        path = ROOT / "app" / "static" / "data" / "alignment_EGFR.json"
        if not path.exists():
            pytest.skip("bundles not built")
        bundle = json.loads(path.read_text())
        kinds = bundle["liability_kinds"]
        checked = 0
        for row in bundle["rows"]:
            for residue, kind in row["lia"]:
                name = kinds[kind]
                # Indexed into the row as the panel draws it, not into the
                # full sequence: a row that dropped an unnumbered N-terminal
                # extension has no column for those residues.
                assert 0 <= residue < len(row["seq"]), (row["id"], name)
                letter = row["seq"][residue]
                if name.startswith("deamidation") or name.startswith("n_glyco"):
                    assert letter == "N", (row["id"], name, letter)
                elif name.startswith("isomerisation") or name == "fragmentation_dq":
                    assert letter == "D", (row["id"], name, letter)
                elif name == "oxidation_met":
                    assert letter == "M", (row["id"], name, letter)
                elif name == "oxidation_trp":
                    assert letter == "W", (row["id"], name, letter)
                elif name == "unpaired_cysteine":
                    assert letter == "C", (row["id"], name, letter)
                checked += 1
        assert checked > 500
    finally:
        conn.close()


def _controls(target: str):
    if not DB.exists():
        pytest.skip("no database built")
    conn = sqlite3.connect(f"file:{DB}?mode=ro", uri=True)
    conn.row_factory = sqlite3.Row
    try:
        return conn.execute(
            "SELECT row_id, label, residues, sse, status AS role "
            "FROM alignment WHERE target = ? AND kind = 'control'",
            (target,)).fetchall()
    finally:
        conn.close()


def test_every_target_gets_at_least_one_known_control_row() -> None:
    for target in ("EGFR", "GLP1R", "PRLR", "TL1A", "TNFa"):
        rows = _controls(target)
        assert rows, f"{target} has no control row in the alignment"
        assert all(r["residues"] for r in rows), target


def test_every_control_row_is_actually_a_nanobody() -> None:
    """The row says "known control nanobody", so it has to be one.

    SAbDab-nano carries conventional VH domains and TCR variable domains as
    well as VHHs, and selecting by antigen rather than by sequence picks them
    up: five of the twelve chains that matched an EGFR or TNFa complex were
    not nanobodies. ANARCI numbers all of them happily, because they are all
    V-set folds, so the hallmark tetrad is what separates them.
    """
    from pipeline import numbering as nb

    for target in ("EGFR", "GLP1R", "PRLR", "TL1A", "TNFa"):
        for row in _controls(target):
            numb = nb.number_many([("q", row["residues"])])["q"]
            _intact, n_camelid, detail = nb.hallmark(numb)
            assert n_camelid >= 3, (target, row["label"], detail)


def test_a_control_s_secondary_structure_fits_its_own_row() -> None:
    for target in ("EGFR", "GLP1R", "PRLR", "TL1A", "TNFa"):
        for row in _controls(target):
            if not row["sse"]:
                continue
            assert len(row["sse"]) == len(row["residues"]), row["label"]
            assert set(row["sse"]) <= set("HEC.-"), row["label"]


def test_a_receptor_does_not_count_as_its_cytokine() -> None:
    """"tumor necrosis factor receptor superfamily member 17" is BCMA, and it
    starts with the string that names TNF-alpha."""
    from pipeline.neighbours import antigen_matches

    assert antigen_matches(["tumor necrosis factor"], "TNFa")
    assert antigen_matches(["glycerol|sulfate ion|tumor necrosis factor"], "TNFa")
    assert not antigen_matches(
        ["tumor necrosis factor receptor superfamily member 17"], "TNFa")
    assert not antigen_matches(
        ["tumor necrosis factor receptor superfamily member 9|sulfate ion"], "TNFa")
    assert antigen_matches(["epidermal growth factor receptor"], "EGFR")
