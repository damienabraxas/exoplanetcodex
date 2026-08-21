"""RYA-931: the Molecfit runtime blocker was an empty SCIENCE table, not MIPAS.

RYA-927 stopped at rc=13 with molecfit reporting

    mf_config_chk_ref_atmos: Unable to open rat data file .../profiles/mipas/equ.fits!
    molecfit_model: Bad Specifications of Molecules and Reference Atmosphere

The reference atmosphere was never the defect.  The HARPS Phase-3 `ERR` column
is entirely NaN, so filtering on finite errors selected zero pixels; molecfit
failed in `molecfit_data_extract_spectrum_from_cpl_table` and only afterwards
printed the reference-atmosphere line off an already-dirty CPL error state.

These tests pin the invariants, not the example: an unusable *error* column must
never be allowed to empty the *flux* table, and a genuinely empty band must fail
by naming the flux, not by handing molecfit an empty table to misdiagnose.
"""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest
from astropy.io import fits

ROOT = Path(__file__).resolve().parents[1]
EVIDENCE = ROOT / "data" / "audit" / "rya931_molecfit_runtime"


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


runner = _load("rya931_molecfit_model")


def _harps_adp(path: Path, *, err_all_nan: bool = True, flux_finite: bool = True) -> Path:
    """A minimal stand-in for a HARPS Phase-3 S1D ADP over the O2 B band."""
    wave = np.arange(6830.0, 6915.0, 0.01)
    flux = np.full(wave.size, 1.0e5) if flux_finite else np.full(wave.size, np.nan)
    err = np.full(wave.size, np.nan) if err_all_nan else np.full(wave.size, 100.0)
    header = fits.Header()
    header["INSTRUME"] = "HARPS"
    header["SPECSYS"] = "BARYCENT"
    header["MJD-OBS"] = 60158.580736
    header["UTC"] = 50171.0
    header["HIERARCH ESO DRS BERV"] = 0.469737
    for key, value in {"ESO TEL ALT": 27.2, "ESO TEL AMBI RHUM": 8.0,
                       "ESO TEL AMBI PRES START": 773.2, "ESO TEL AMBI TEMP": 20.9,
                       "ESO TEL GEOELEV": 2400.0, "ESO TEL GEOLON": -70.7345,
                       "ESO TEL GEOLAT": -29.2584}.items():
        header[key] = value
    n = wave.size
    table = fits.BinTableHDU.from_columns([
        fits.Column(name="WAVE", format=f"{n}D", array=wave.reshape(1, -1)),
        fits.Column(name="FLUX", format=f"{n}D", array=flux.reshape(1, -1)),
        fits.Column(name="ERR", format=f"{n}D", array=err.reshape(1, -1)),
    ])
    fits.HDUList([fits.PrimaryHDU(header=header), table]).writeto(path, overwrite=True)
    return path


def test_all_nan_error_column_does_not_empty_the_science_table(tmp_path):
    """The defect: an unusable ERR column must not remove every FLUX pixel."""
    built = runner.write_inputs(_harps_adp(tmp_path / "adp.fits"), None, tmp_path / "in")
    rows = len(fits.open(tmp_path / "in" / "science.fits")[1].data)
    assert rows == built["science_rows"] > 0, "an all-NaN ERR column emptied the table again"
    assert built["source_error_column_usable"] is False
    assert built["source_error_finite_pixels"] == 0


def test_unusable_error_column_is_reported_not_silently_absorbed(tmp_path):
    """`applied` must never rest on a silently-substituted uncertainty."""
    built = runner.write_inputs(_harps_adp(tmp_path / "adp.fits"), None, tmp_path / "in")
    assert "dflux" not in [c.name for c in fits.open(tmp_path / "in" / "science.fits")[1].columns]
    usable = runner.write_inputs(
        _harps_adp(tmp_path / "ok.fits", err_all_nan=False), None, tmp_path / "in2")
    assert usable["source_error_column_usable"] is True
    assert "dflux" in [c.name for c in fits.open(tmp_path / "in2" / "science.fits")[1].columns]


def test_empty_band_fails_by_naming_the_flux(tmp_path):
    """A genuinely empty band must not be handed to molecfit to misdiagnose."""
    with pytest.raises(SystemExit) as excinfo:
        runner.write_inputs(_harps_adp(tmp_path / "adp.fits", flux_finite=False),
                            None, tmp_path / "in")
    message = str(excinfo.value)
    assert "FLUX" in message and "empty SCIENCE table" in message
    assert "mipas" not in message.lower() and "atmosphere" not in message.lower()


def test_sof_entries_survive_a_path_containing_a_space(tmp_path):
    """esorex splits SOF lines on whitespace; the store root is 'Exoplanet Codex'."""
    destination = tmp_path / "Exoplanet Codex" / "in"
    runner.write_inputs(_harps_adp(tmp_path / "adp.fits"), None, destination)
    for line in (destination / "model.sof").read_text().splitlines():
        assert len(line.split()) == 2, f"SOF line splits into >2 tokens: {line!r}"
        assert not line.startswith("/"), "absolute SOF paths break on the space in the store root"


def test_mipas_reference_atmosphere_was_never_the_defect():
    """Positive control: the blamed equ.fits is what the *successful* run loaded.

    An absence needs a positive control -- 'MIPAS is fine' is only meaningful if
    the same file demonstrably works.  The recorded run log is that control.
    """
    evidence = json.loads((EVIDENCE / "runtime_rca.json").read_text())
    assert evidence["failing_run"]["returncode"] == 13
    assert evidence["failing_run"]["science_rows"] == 0
    assert evidence["passing_run"]["returncode"] == 0
    assert evidence["passing_run"]["science_rows"] > 0
    # same reference atmosphere on both sides
    assert (evidence["failing_run"]["reference_atmosphere"]
            == evidence["passing_run"]["reference_atmosphere"])
    assert evidence["passing_run"]["reference_atmosphere_opened"] is True


def test_fitted_o2_column_is_physically_consistent():
    """O2 is well mixed, so a correct fit returns ~1.0 x the reference column.

    This is the runtime's own referee: it is not a tuned target, and the first
    RYA-931 run -- which converged, produced products and looked successful --
    failed it at 0.036 because its line-spread function had run away.
    """
    report = json.loads((EVIDENCE / "exposure_correction.json").read_text())
    assert len(report) == 10, "all ten direct-solar exposures must be corrected"
    columns = np.array([r["rel_mol_col_O2"] for r in report])
    assert np.all(np.abs(columns - 1.0) < 0.10), f"O2 column ran away: {columns}"


def test_quarantine_threshold_is_measured_not_constant():
    """T_min must track each exposure's measured transmission accuracy."""
    report = json.loads((EVIDENCE / "exposure_correction.json").read_text())
    budget = _load("rya931_correct_exposures").CORRECTED_RELATIVE_ERROR_BUDGET
    for row in report:
        assert row["t_min"] == pytest.approx(row["dt_transmission_sigma"] / budget, rel=1e-9)


def test_correction_is_confined_to_the_fitted_region():
    """`mtrans` is 0 outside the fit; dividing there would zero the blue spectrum.

    The invariant is asserted on the RAW co-add, which is where it is exactly
    true.  Normalised flux cannot carry it: continuum fitting is global, so
    recovering the O2 band moves the spline everywhere by a small, bounded
    amount.  Asserting bit-equality on the normalised product would either fail
    honestly or have to be weakened until it proved nothing.
    """
    report = json.loads((EVIDENCE / "verification.json").read_text())
    preservation = report["solar_preservation"]
    assert preservation["raw_coadd_bit_identical_outside_fitted_region"] is True
    assert (preservation["max_fractional_continuum_shift_outside"]
            <= preservation["continuum_shift_tolerance"])


def test_o2b_gate_and_solar_preservation_pass():
    report = json.loads((EVIDENCE / "verification.json").read_text())
    gate, verdict = report["o2b_gate"], report["verdict"]
    assert gate["o2b_before"]["pct_below_0.5"] > 20.0, "the before-fixture must be uncorrected"
    assert verdict["o2b_saturation_cleared"] is True
    assert verdict["control_shift_bounded"] is True
    assert verdict["solar_preserved"] is True
    assert verdict["iag_agreement_improved"] is True
    assert verdict["pass"] is True
