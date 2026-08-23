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
    # RYA-821: was REQUEST_ONLY. That word was wrong for O -- we hold atom.o41f AND the
    # <3D> STAGGER atmosphere, so there is nothing to request. What is missing is the
    # departure solve, i.e. the tier-2 solver. TIER2_OWED says that; REQUEST_ONLY sent a
    # reader hunting for a grid nobody publishes.
    assert cells[("O", "MEAN3D_NLTE")]["state"] == "TIER2_OWED"
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


# ── RYA-821: mean-3D reachability states say WHAT IT WOULD TAKE ───────────────

def test_mean3d_states_distinguish_the_three_situations(cells):
    """"REQUEST_ONLY" was collapsing three unrelated cases into one word:
    a deck sitting unwired, a solver we have not built, and something we must ask for.
    """
    assert cells[("Al", "MEAN3D_NLTE")]["state"] == "TIER1_WIRED"
    for el in ("Cr", "Y", "Eu"):                      # deck + atom on disk
        assert cells[(el, "MEAN3D_NLTE")]["state"] == "TIER1_WIREABLE", el
    for el in ("O", "Fe", "Ba", "Ni"):                # atom held, no deck
        assert cells[(el, "MEAN3D_NLTE")]["state"] == "TIER2_OWED", el
    for el in ("Li", "C", "K"):                       # neither
        assert cells[(el, "MEAN3D_NLTE")]["state"] == "REQUEST_ONLY", el


def test_a_solar_increment_never_masks_an_unwired_deck(cells):
    """Cr was reported HAVE on solar3d_metals -- a SOLAR-ONLY scalar -- while its full
    parameter-space <3D> deck sat unwired. A HAVE that hides a wireable deck is the
    worst kind of wrong: it looks finished."""
    cr = cells[("Cr", "MEAN3D_NLTE")]
    assert cr["state"] == "TIER1_WIREABLE"
    assert any("solar-only increment" in f for f in cr["facts"])


def test_tier2_owed_is_not_an_acquisition_task(cells):
    """These need the solver, not a download -- saying otherwise sends someone hunting
    for a grid that nobody publishes."""
    fix = cells[("O", "MEAN3D_NLTE")]["fix"]
    assert "tier-2 solver" in fix and "NOT by acquisition" in fix
