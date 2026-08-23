"""RYA-1015 — guard the element x model availability matrix section."""
from __future__ import annotations

import csv
from pathlib import Path

import pytest

from pipeline.model_availability_matrix import (
    CANONICAL_28, MODEL_TYPES, DiskSnapshotError, build_matrix, classify_disk_file,
    load_disk_snapshot,
)

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def cells():
    return {(c["element"], c["model_type"]): c for c in build_matrix()["cells"]}


def test_matrix_covers_all_28_x_model_types():
    assert len(CANONICAL_28) == 28  # 27 (incl. Fe II) + Zn
    assert len(build_matrix()["cells"]) == 28 * len(MODEL_TYPES)


def test_oxygen_row(cells):
    """RYA-1015 pre-registered O as 1D-NLTE HAVE / <3D> REQUEST_ONLY / full-3D NONE.

    The full-3D expectation is SUPERSEDED and deliberately not asserted as NONE.
    That value came from an earlier reading in which a cell meant "can we COMPUTE full
    3D ourselves" -- under which nothing qualifies. Ryan corrected it: we DO hold
    3D-NLTE for O (the Amarsi-2019 CNO tables), and a matrix that renders real holdings
    as NONE is wrong. A cell now means "do we HOLD a runnable/applicable grid", so O
    full-3D is HAVE. The other two legs are unchanged and still pinned.
    """
    assert cells[("O", "1D_NLTE")]["state"] == "HAVE"
    assert cells[("O", "MEAN3D_NLTE")]["state"] == "REQUEST_ONLY"
    assert cells[("O", "FULL_3D_NLTE")]["state"] == "HAVE"
    assert "amarsi2019_cno" in (cells[("O", "FULL_3D_NLTE")]["code_grid"] or "")


def test_rya820_grids_are_held_and_flagged_unwired(cells):
    """Li/Mg/Na off-solar 3D-NLTE: fetched 2026-08-23, verified, NOT yet wired."""
    for element in ("Li", "Mg", "Na"):
        cell = cells[(element, "FULL_3D_NLTE")]
        assert cell["state"] == "HAVE", element
        assert "threed_offsolar" in (cell["code_grid"] or ""), element
        assert any("NOT YET WIRED" in f for f in cell["facts"]), element


def test_rya597_ti_drift_is_fixed(cells):
    """CSV cited Bergemann2011 while the code has run Mallinson2024 since RYA-545."""
    assert cells[("Ti", "1D_NLTE")]["csv_claim"] == "Ti_Mallinson2024_PySME.csv"
    assert cells[("Ti", "1D_NLTE")]["code_grid"] == "Ti_Mallinson2024_PySME.csv"
    row = next(r for r in csv.DictReader(
        (ROOT / "data" / "curation" / "nlte_grid_availability.csv").open())
        if r["element"] == "Ti" and r["subsystem"] == "registry-nlte")
    assert "Bergemann2011" not in row["grid_file"]


def test_elements_wired_outside_the_registry_are_not_false_problems(cells):
    """C/N/O run via cno-3dnlte and Fe via fe-nlte, not NLTE_CORRECTION_ELEMENTS."""
    for element in ("C", "N", "O", "Fe"):
        assert cells[(element, "1D_NLTE")]["state"] == "HAVE", element


def test_blind_scan_is_refused_rather_than_reported_as_absent(tmp_path):
    """find without -L returns nothing through the grids symlink (RYA-1013); a matrix
    built on that reports absences that are pure artifact."""
    blind = tmp_path / "blind.txt"
    blind.write_text("somedir|nlte_O_scatt_pysme.grd|1\n")
    with pytest.raises(DiskSnapshotError, match="control"):
        load_disk_snapshot(blind)
    assert "CONTROL=PASS" in (
        ROOT / "data" / "audit" / "rya1015" / "sirius_scan_raw.txt").read_text()


def test_co_molecule_is_not_parsed_as_cobalt(cells):
    """NLTEgrid4TS_CO_MARCS is carbon monoxide; mapping it to Co invents a grid."""
    assert classify_disk_file("NLTEgrid4TS_CO_MARCS_Nov-14-2024.bin") is None
    assert classify_disk_file("NLTEgrid4TS_MN_MARCS_Mar-15-2023.bin") == ("Mn", "1D_NLTE")
    assert cells[("Co", "1D_NLTE")]["state"] == "NONE"


def test_matrix_is_a_standing_935_view():
    import scripts.rya935_live_status as live
    assert hasattr(live, "collect_model_matrix")
