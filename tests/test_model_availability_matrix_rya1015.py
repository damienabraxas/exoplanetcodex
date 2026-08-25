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
    # 🔴 VOCABULARY CHANGED, and the new word is the more accurate one. RYA-821/1013
    # split the old REQUEST_ONLY into two states because they are different problems:
    # REQUEST_ONLY means neither deck nor atom, so it genuinely must be asked for;
    # T2_BUILD_OWED means no deck BUT we hold the model atom and the <3D> STAGGER
    # atmosphere, so it is reachable by the build-our-own route with nothing to fetch
    # and nobody to email. O is in `_MODEL_ATOMS`, so calling it REQUEST_ONLY told a
    # reader to go asking for something already on disk.
    #
    # 🔴 CORRECTED AGAIN, RYA-1035 — AND THE PREVIOUS WORDING HAD THE SAME DEFECT IT WAS
    # FIXING, ONE STEP FURTHER OUT. "no deck BUT we hold the model atom ... with nothing to
    # fetch" is a claim about the SOURCE derived from a scan of OUR DISK. Live listing of
    # the MPG Keeper share on 2026-08-24: 17 of 18 element folders ship a STAGGERmean3D
    # deck, and O's is one of them (NLTEgrid4TS_O_STAGGERmean3D_May-18-2021.bin.zip,
    # 44.5 MB) — on the same share Al/Cr/Eu/Y were pulled from.
    #
    # So O <3D> is an ACQUISITION, not a build, and the RYA-1013/1014 O tier-2 BUILD may
    # not be needed at all. That is a live question for those tickets; what is settled here
    # is only that the matrix must stop telling their reader there is nothing to fetch.
    #
    # 🔴 AND THEN IT WAS FETCHED (2026-08-25), so the cell moved again: T2_FETCH_OWED ->
    # T2_CONSUME_READY. That is THREE corrections to this one assertion, and each was a
    # real state change rather than churn -- REQUEST_ONLY (wrong: we held the atom) ->
    # T2_BUILD_OWED (wrong: the vendor publishes a deck) -> T2_FETCH_OWED (true for a day)
    # -> T2_CONSUME_READY (staged and verified). The value moves because the WORLD moved;
    # what must not move is that the cell is DERIVED rather than asserted.
    #
    # Verified by direct read at the solar node before this line was changed: ndep 101,
    # nlev 41 (matching atom.o41f), 21-value A(O) axis, median b -> 1.0000 at depth, and
    # tau agreement with the <3D> atmosphere of 2.45e-08. So RYA-1013's CONSUME leg works
    # -- by the Gerber deck, not the nlte.mpia.de route that leg names and that returns
    # exactly 0.000 for O I 7771.9.
    assert cells[("O", "MEAN3D_NLTE")]["state"] == "T2_CONSUME_READY"
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
