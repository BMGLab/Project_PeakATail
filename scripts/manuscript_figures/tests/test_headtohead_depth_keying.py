#!/usr/bin/env python3
"""
Regression tests for benchmark_headtohead.depth_map (Stage 0g).

THE BUG (manuscript/10_caller_fix_plan.md, item 0g)
    depth_map() built the PeakATail name->depth mapping with
        dict(zip(names_from_pasbed_bed, rowsums_of_annotated_matrix))
    pasbed.bed is written in COORDINATE order; the matrix rows are in
    PAS-ID order (the order recorded in
    05_annotated_matrix/<strategy>/annotated_pas_ids.tsv).  The two orders
    are unrelated, so essentially every PAS was handed another PAS's depth.
    `assert len(names) == len(sums)` cannot see this: both orderings are
    permutations of the same set, so the lengths always agree.

Run:
    python3 -m pytest scripts/manuscript_figures/tests/test_headtohead_depth_keying.py -v
"""
import importlib.util
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
SCRIPT = HERE.parent / "benchmark_headtohead.py"


def _load():
    """Import benchmark_headtohead.py by path; assert it is THIS worktree's copy."""
    spec = importlib.util.spec_from_file_location("_bhh_under_test", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["_bhh_under_test"] = mod
    spec.loader.exec_module(mod)
    assert Path(mod.__file__).resolve() == SCRIPT.resolve(), mod.__file__
    return mod


# --- synthetic run directory -------------------------------------------------
# 4 PAS.  Coordinate order (pasbed.bed) and PAS-ID order (matrix rows) are
# deliberately DIFFERENT permutations, and every PAS gets a distinct depth so
# any mis-pairing is visible.
#
#   pas_id : depth   coordinate rank
#      70  :   7          2
#      12  :  30          0
#      95  :   5          3
#      41  :  10          1
TRUTH = {"70": 7, "12": 30, "95": 5, "41": 10}


@pytest.fixture
def fake_run(tmp_path):
    """Build results/benchmark_tools/<MOUSE>/peakatail/mouse1/run/… ."""
    bt = tmp_path / "benchmark_tools"
    run = bt / "ds" / "peakatail" / "mouse1" / "run"
    (run / "05_annotated_matrix" / "default").mkdir(parents=True)

    # pasbed.bed: BED6 in coordinate order, name column = pas_id
    coords = [("1", 1000, 1100), ("1", 2000, 2100), ("2", 500, 600), ("2", 900, 1000)]
    coord_order = ["41", "12", "95", "70"]      # coordinate order != id order
    with open(run / "pasbed.bed", "w") as fh:
        for (c, s, e), pid in zip(coords, coord_order):
            fh.write(f"{c}\t{s}\t{e}\t{pid}\t0\t+\n")

    # annotated_pas_ids.tsv: the matrix's true row key, numerically sorted
    row_order = sorted(TRUTH, key=int)          # ['12', '41', '70', '95']
    (run / "05_annotated_matrix" / "default" / "annotated_pas_ids.tsv").write_text(
        "pas_id\n" + "\n".join(row_order) + "\n")

    # annotated_matrix.mtx: one cell, row i holds TRUTH[row_order[i-1]]
    lines = ["%%MatrixMarket matrix coordinate integer general", "%",
             f"{len(row_order)} 1 {len(row_order)}"]
    lines += [f"{i} 1 {TRUTH[pid]}" for i, pid in enumerate(row_order, 1)]
    (run / "annotated_matrix.mtx").write_text("\n".join(lines) + "\n")
    return bt, coord_order, row_order


def test_peakatail_depth_is_keyed_on_annotated_pas_ids(fake_run, monkeypatch):
    """Each PAS must receive ITS OWN row sum, not the row at its BED position."""
    bt, coord_order, row_order = fake_run
    assert coord_order != row_order, "fixture must exercise the two orderings"
    mod = _load()
    monkeypatch.setattr(mod, "BT", bt)
    monkeypatch.setattr(mod, "MOUSE", "ds")
    assert mod.depth_map("PeakATail", "mouse1") == TRUTH


def test_length_assert_alone_cannot_detect_the_mispairing(fake_run):
    """Documents why the shipped assert passed: same set, different order."""
    bt, coord_order, row_order = fake_run
    assert len(coord_order) == len(row_order)
    assert set(coord_order) == set(row_order)
    assert coord_order != row_order


def test_depth_map_rejects_a_pas_id_set_mismatch(fake_run, monkeypatch):
    """A real guard: BED names and matrix row ids must be the SAME set."""
    bt, _, _ = fake_run
    ids = bt / "ds/peakatail/mouse1/run/05_annotated_matrix/default/annotated_pas_ids.tsv"
    ids.write_text("pas_id\n12\n41\n70\n9999\n")   # 95 -> 9999
    mod = _load()
    monkeypatch.setattr(mod, "BT", bt)
    monkeypatch.setattr(mod, "MOUSE", "ds")
    with pytest.raises(AssertionError):
        mod.depth_map("PeakATail", "mouse1")


@pytest.mark.parametrize("m", ["mouse1", "mouse2"])
def test_real_run_bed_order_is_not_row_order(m):
    """
    Guard on the SHIPPED artifacts: if these ever coincide the unit test above
    would stop discriminating.  Skips when the benchmark outputs are absent.
    """
    run = (Path("/mnt/ssd1/Projects/PeakATail_wd/results/benchmark_tools/gse104556")
           / "peakatail" / m / "run")
    bed = run / "pasbed.bed"
    ids = run / "05_annotated_matrix" / "default" / "annotated_pas_ids.tsv"
    if not (bed.exists() and ids.exists()):
        pytest.skip("benchmark outputs not present")
    names = [l.split("\t")[3] for l in bed.read_text().splitlines()]
    rows = ids.read_text().splitlines()[1:]
    assert set(names) == set(rows)
    assert names != rows, "orders coincide -- re-check the fixture assumption"
