"""RYA-925 — non-Fe Engine-A products consult their registered grid."""

import numpy as np
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))


def test_al_engine_a_uses_the_registered_per_line_resolver(monkeypatch):
    import derive_band_products as dbp
    import pipeline.nlte_corrections as nc

    seen = []

    def fake_delta(element, wave, teff, logg, feh, tol):
        seen.append((element, wave, teff, logg, feh, tol))
        return -0.025 if abs(wave - 7836.134) < 0.01 else np.nan

    monkeypatch.setattr(nc, "_mpia_element_delta", fake_delta)
    got = dbp.engine_a_delta("Al", "I", np.array([7836.134, 7361.568]))

    assert got == {7836.134: -0.025}
    assert [row[0] for row in seen] == ["Al", "Al"]
    assert all(row[-1] == 0.06 for row in seen)


def test_al_engine_a_provenance_names_amarsi_not_mpia():
    import derive_band_products as dbp

    source = dbp.engine_a_source("Al")
    assert "Al_Amarsi2020_PySME.csv" in source
    assert "Nordlander & Lind 2017" in source
    assert "MPIA" not in source
    assert dbp.engine_a_model("Al") == "amarsi"


def test_axes_allow_a_live_emitter_to_override_the_legacy_engine_a_family():
    from pipeline.treatment_axes import axes_for

    axes = axes_for("ENGINE-A", handler="ProfileFitHandler", model="amarsi")
    assert axes.model == "amarsi"
    assert axes.display == "EW · 1D-NLTE · Amarsi"


def test_al_litscan_and_validation_contract():
    from pipeline.litscan import literature_range
    from pipeline.validate_element import BandProduct, validate_element

    lit = literature_range("Al")
    assert lit is not None
    assert (lit.central, lit.min, lit.max, lit.scale) == (6.43, 6.39, 6.47, "3D-NLTE")
    verdict = validate_element("Al", [
        BandProduct("VIS", "EW · 1D-NLTE · Amarsi", 6.43,
                    sigma_statistical=0.03, sigma_systematic=0.05,
                    n_lines=2, scale="1D-NLTE"),
        BandProduct("red-optical", "EW · 1D-LTE", 6.50,
                    sigma_statistical=0.04, sigma_systematic=0.17,
                    n_lines=4, scale="1D-LTE"),
    ], asplund=6.43)
    assert verdict.validating[0].verdict == "pass"
    assert next(b for b in verdict.bands if b.band == "red-optical").verdict == "report"


def test_rya925_calibration_gate_uses_combined_uncertainty():
    from scripts.rya925_al_report import calibration_matrix, product_matrix

    got = calibration_matrix(product_matrix())
    vis = got[got.band.eq("VIS")]
    assert set(vis.calibration_verdict) == {"REPLICATED-WITHIN-1SIGMA"}
    assert vis.z_combined.max() < 1.0
    assert not got.band.eq("near-UV").any()  # attempted line is appendix-only, no product
    assert set(got.holding) == {"solar_kpno"}


def test_rya925_tracker_declares_harps_as_loader_ready():
    source = (ROOT / "scripts/rya925_al_report.py").read_text()
    assert '"id": "solar_harps"' in source
    assert '"source": "RYA-911 / PR #313"' in source
    assert '"readiness": "loader ready; Al run owed"' in source


def test_rya925_tracker_uses_rya851_band_sectioned_forest_plot():
    page = (ROOT / "data/results/rya925/tracker.html").read_text()
    assert "RYA-851 forest-plot convention" in page
    assert "forest-band" in page and "forest-row" in page
    assert "forest-lit-band" in page and "forest-reference" in page
    assert "No engines or instruments are averaged" in page


def test_rya925_3057_is_appendix_only_not_a_nearuv_product():
    from scripts.rya925_al_report import disposition, per_line_matrix, product_matrix

    assert not product_matrix().band.eq("near-UV").any()
    line = per_line_matrix().query("abs(wavelength_air_A - 3057.144) < 0.01").iloc[0]
    assert not bool(line.in_aggregate)
    assert "ABUNDANCE-INSENSITIVE" in line.excluded_reason
    disp = disposition().query("abs(lambda_A - 3057.144) < 0.01").iloc[0]
    assert disp.disposition == "rejected"
    assert disp.permitted_routes == "synth only"


def test_frontier_fit_uses_context_element_not_fe_source():
    import inspect
    from scripts import rya759_nearuv_fe_product as frontier

    source = inspect.getsource(frontier.fit_one)
    assert "ctx.get('element', 'Fe')" in source
    assert "ctx['solar_abund'], 'Fe'" not in source
