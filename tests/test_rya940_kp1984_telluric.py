"""RYA-940: correcting an atlas with no observation metadata.

The 1984 Kitt Peak atlas cannot satisfy the ratified observation-night-GDAS rule --
it records no airmass, no date and no site -- so the columns are fitted freely and
the result is refereed externally. These tests pin what makes that defensible.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
EVIDENCE = ROOT / "data" / "audit" / "rya940_kp1984_telluric"
PRODUCTS = ROOT / "data" / "processed" / "kp1984_telluric_corrected"


def _bands() -> dict:
    return {p.parent.name: json.loads(p.read_text())
            for p in EVIDENCE.glob("*/fit_manifest.json")}


def test_every_band_records_that_no_observation_night_exists():
    """The deviation must be recorded on every product, never assumed once."""
    for name, d in _bands().items():
        assert "none" in d["gdas"], name
        assert "NO observation-night profile" in d["gdas"], name


def test_o2_columns_agree_because_o2_is_well_mixed():
    """The referee that replaces RYA-931's 'column must be 1.0'.

    The airmass is unknown, but it must be the SAME unknown in every O2 band.
    Three bands 1300 A apart, fitted independently, must land together.
    """
    columns = []
    for name, d in _bands().items():
        if not d.get("accepted") or d["molecules"] != ["O2"]:
            continue
        a = d["lsf_attempts"][-1]
        columns.append(a.get("column_scale", a.get("implied_airmass")))
    assert len(columns) >= 3, "need at least three O2 bands for this to be a check"
    assert max(columns) / min(columns) < 1.15, f"O2 bands disagree on airmass: {columns}"


def test_water_columns_are_not_labelled_airmass():
    """H2O is NOT well mixed; calling its column an airmass invents an agreement."""
    for name, d in _bands().items():
        if not d.get("accepted"):
            continue
        a = d["lsf_attempts"][-1]
        if "H2O" in d["molecules"] and "column_is_airmass_proxy" in a:
            assert a["column_is_airmass_proxy"] is False, name


def test_a_band_with_no_external_reference_says_so():
    """10000-13000 A has no telluric-free product; that must be on the record."""
    ir = [d for d in _bands().values() if d["band_A"][0] >= 10000]
    assert ir, "the IR band manifest is missing"
    for d in ir:
        assert d.get("externally_validated") is False
        assert d.get("solar_reference") is None


def test_the_unconstrained_band_is_recorded_as_a_null_not_forced():
    """H2O 7160-7340 cannot constrain the LSF. Fit quality alone must not pass it."""
    d = _bands().get("h2o7160")
    assert d is not None
    assert d["accepted"] is None or d["accepted"] is False
    attempts = d["lsf_attempts"]
    assert len(attempts) >= 6, "the whole ladder must be tried before declaring a null"
    # Not `== 0.0`: one start returns 4.4e-16. The meaningful statement is that every
    # start lands BELOW the admissibility floor, not that it lands on an exact value.
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "kp", ROOT / "scripts" / "rya940_kp1984_correct.py")
    kp = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(kp)
    assert all(a["gaussfwhm"] < kp.ACCEPT_MIN_LSF_FWHM_PIX for a in attempts), (
        "the null is an LSF collapsed below the physical floor")
    assert all(not a["admissible"] for a in attempts)
    # The point of the null: the fit quality was FINE. Quality is not admissibility.
    assert min(a["rms_rel_to_mean"] for a in attempts) < 0.05


def test_the_lsf_ladder_is_derived_from_resolving_power_not_pixels():
    """Sampling runs 0.0067-0.0183 A across the atlas; a pixel ladder cannot travel."""
    source = (ROOT / "scripts" / "rya940_kp1984_correct.py").read_text()
    assert "ATLAS_RESOLVING_POWER" in source
    assert "LSF_START_LADDER_PIX" not in source, "the fixed pixel ladder must be gone"
    for d in _bands().values():
        if "lsf_start_ladder_pix" in d:
            assert d["expected_lsf_fwhm_pix"] > 0
            assert d["lsf_start_ladder_pix"][0] == pytest.approx(
                d["expected_lsf_fwhm_pix"], rel=1e-3)


def test_neither_the_continuum_nor_the_wavelength_solution_is_fitted():
    """Both ran away when free, and the atlas supplies both (RYA-911 / RYA-938)."""
    source = (ROOT / "scripts" / "rya940_kp1984_correct.py").read_text()
    assert '"--FIT_CONTINUUM=0"' in source
    assert '"--FIT_WLC=0"' in source


def test_acceptance_is_error_independent():
    """chi2 is computed against an ASSUMED error here, so it cannot be the gate."""
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "kp", ROOT / "scripts" / "rya940_kp1984_correct.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    assert not hasattr(m, "ACCEPT_MAX_REDUCED_CHI2")
    assert m.ACCEPT_MAX_RMS_REL_TO_MEAN > 0
    assert m.ACCEPT_MIN_REL_COL >= 0.5 and m.ACCEPT_MAX_REL_COL <= 5.0


def test_corrected_products_quarantine_rather_than_force():
    """Saturated cores are NaN, never divided (RYA-927 ratified policy)."""
    import numpy as np
    files = sorted(PRODUCTS.glob("kp1984_corrected_*.txt"))
    assert files, "no corrected products committed"
    for f in files:
        d = np.loadtxt(f)
        assert d.shape[1] == 3, f.name
        assert np.all(np.diff(d[:, 0]) > 0), f"{f.name}: wavelengths not increasing"
        # Assert the RELATIONSHIP, not a count. A weak band may legitimately have
        # nothing to quarantine -- O2 gamma's minimum transmission (0.470) never
        # reaches its T_min (0.087) -- and demanding NaNs there would fail a
        # perfectly correct product.
        flux, trans = d[:, 1], d[:, 2]
        quarantined = np.isnan(flux)
        if quarantined.any():
            assert trans[quarantined].max() < trans[~quarantined].min() + 1e-9, (
                f"{f.name}: quarantine is not a clean transmission threshold")
        # transmission is a transmission
        assert np.nanmin(trans) >= 0.0 and np.nanmax(trans) <= 1.0000001, f.name


def test_the_original_atlas_is_untouched():
    """A corrected product is a sibling, never an overwrite."""
    for name, d in _bands().items():
        for seg in d["segments"]:
            assert seg["segment"].startswith("lm")
            assert "md5" in seg, "source segments must be checksummed"
