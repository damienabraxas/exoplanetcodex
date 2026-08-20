import datetime as dt
from pathlib import Path

import pytest
from astropy.io import fits

from pipeline.telluric import harps_preflight as hp


def _write_harps(path: Path, *, program="1102.D-0954(A)"):
    h = fits.Header()
    h["INSTRUME"] = "HARPS"
    h["MJD-OBS"] = 60158.581
    h["DATE-OBS"] = "2023-08-02T13:56:15.570"
    h["OBJECT"] = "SUN"
    h["HIERARCH ESO OBS PROG ID"] = program
    h["HIERARCH ESO TEL AIRM START"] = 2.1826
    h["HIERARCH ESO TEL AIRM END"] = 2.1826
    h["HIERARCH ESO TEL AMBI PRES START"] = 773.2
    h["HIERARCH ESO TEL AMBI TEMP"] = 20.9
    h["HIERARCH ESO TEL AMBI RHUM"] = 8.0
    fits.PrimaryHDU(header=h).writeto(path)


def test_exposure_carries_weather_and_real_profile(tmp_path, monkeypatch):
    src = tmp_path / "solar.fits"
    profile = tmp_path / "gdas.fits"
    profile.touch()
    _write_harps(src)
    monkeypatch.setattr(hp, "fetch_gdas", lambda *a, **k: profile)
    rec = hp.inspect_exposure(src, cache_dir=tmp_path)
    assert rec.program_id == "1102.D-0954(A)"
    assert rec.pressure_hpa == 773.2
    assert rec.temperature_c == 20.9
    assert rec.relative_humidity_pct == 8.0
    assert rec.gdas_slot_utc == dt.datetime(2023, 8, 2, 15).isoformat()
    assert rec.gdas_profile == str(profile)


def test_wrong_program_is_refused_before_gdas(tmp_path):
    src = tmp_path / "other.fits"
    _write_harps(src, program="not-the-solar-program")
    with pytest.raises(ValueError, match="not the RYA-927 solar sequence"):
        hp.inspect_exposure(src, cache_dir=tmp_path)


def test_sequence_requires_exact_source_set(tmp_path):
    _write_harps(tmp_path / "only-one.fits")
    with pytest.raises(ValueError, match="expected 10"):
        hp.inspect_sequence(tmp_path, cache_dir=tmp_path)
