from pipeline.si_protocol_audit import ROOT, holding_specs, physical_matches
import pandas as pd


def test_identity_join_rejects_wavelength_coincidence():
    rows = [dict(species="Si II", wavelength_air_A="6371.37", ep_eV="8.121")]
    assert physical_matches(rows, "Si I", 6371.37, 8.121) == []
    assert physical_matches(rows, "Si II", 6371.37, 1.0) == []
    assert len(physical_matches(rows, "Si II", 6371.37, 8.121)) == 1
    assert len(physical_matches(rows * 2, "Si II", 6371.37, 8.121)) == 2


def test_inventory_preserves_corrected_and_raw_holdings_without_loading_flux():
    specs = holding_specs(ROOT / "scripts/measure_band_ew.py")
    assert specs["solar_kpno_molecfit_corrected"]["span"] is None
    assert specs["solar_kpno"]["reader"] != specs["solar_kpno_molecfit_corrected"]["reader"]
    assert specs["solar_iag_reiners2016"]["span"][1] == specs["solar_iag"]["span"][0]


def test_nir_pool_is_published_and_crires_measurement_stays_held():
    pool = pd.read_csv(ROOT / "data/audit/rya1218_si_protocol/si_nir_line_pool.csv")
    assert list(pool.wavelength_air_A) == [11991.57, 11984.20, 12103.54, 12031.50]
    result = ROOT / "data/results/rya1218/si_crires_nir/si_j_band_crires_coverage.csv"
    coverage = pd.read_csv(result)
    assert len(coverage) == 16
    assert set(coverage.status) == {"HOLD"}
    raw = coverage[coverage.holding == "solar_vesta_crires_plus_idp"]
    assert all("TelluricNotCorrected" in x for x in raw.reason)


def test_corrected_crires_diagnostic_covers_y_and_h_without_abundance_feed():
    import json
    p = ROOT / "data/results/rya1218/si_crires_corrected_diagnostic/si_corrected_crires_summary.json"
    d = json.loads(p.read_text())
    assert d["abundance_or_grade_pool"] is False
    assert [(x["band"], x["n_canonical_lines"], x["n_served"]) for x in d["holdings"]] == [
        ("Y", 26, 26), ("H", 236, 129)
    ]
