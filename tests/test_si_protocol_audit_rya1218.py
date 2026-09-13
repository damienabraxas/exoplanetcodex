from pipeline.si_protocol_audit import ROOT, holding_specs, physical_matches


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
