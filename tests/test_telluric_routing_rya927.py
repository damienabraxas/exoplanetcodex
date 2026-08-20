"""RYA-927 — one telluric route contract for every catalog instrument."""
from pathlib import Path

import pandas as pd

from pipeline.telluric.routing import correction_needed_at, route_for, route_for_holding

ROOT = Path(__file__).resolve().parents[1]


def test_every_catalog_instrument_has_an_explicit_route():
    df = pd.read_csv(ROOT / "data/catalog/instrument_catalog.csv")
    allowed = {"molecfit_gdas", "molecfit_gdas+clean_line_selection", "corrected_product", "clean_line_selection", "none", "audit_required"}
    for instrument in df.instrument_id:
        route = route_for(instrument)
        assert route.engine in allowed
        assert route.reason


def test_harps_and_kitt_peak_are_clean_line_routes():
    assert route_for("harps").engine == "molecfit_gdas+clean_line_selection"
    assert route_for("kpno_solar_atlas").engine == "clean_line_selection"
    assert not correction_needed_at(6631.218, "harps")
    assert not correction_needed_at(6696.185, "harps")
    assert correction_needed_at(6872.0, "harps")


def test_correction_required_instruments_use_observation_aware_engine():
    for instrument in ("crires_plus", "nirps", "spirou", "apogee", "nirspec"):
        route = route_for(instrument)
        assert route.engine == "molecfit_gdas"
        assert route.correction_required is True
        assert route.analysis_ready_required is True
    assert correction_needed_at(6872.0, "crires_plus")


def test_unresolved_basis_is_a_stop_not_a_fallback():
    route = route_for("feros")
    assert route.engine == "audit_required"
    assert route.analysis_ready_required is True


def test_crires_holding_state_overrides_instrument_capability():
    corrected = route_for_holding("solar_crires_plus_y_rya794")
    assert corrected.engine == "corrected_product"
    assert corrected.correction_required is False
    raw = route_for_holding("solar_vesta_crires_plus_idp")
    assert raw.engine == "molecfit_gdas"
    assert raw.correction_required is True
