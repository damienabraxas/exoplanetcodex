"""
tests/test_gdas_fetch_rya380.py
===============================
RYA-380 — the reusable per-night GDAS retriever. The contract that matters:
NO silent standard-atmosphere fallback. An unavailable profile RAISES GDASUnavailable
(the RYA-373 CRITICAL bug never originates here); the ESO tarball backend returns a real
per-night profile; the nearest-3-hourly slot resolution + ASCII→FITS conversion are
correct; coordinates come from constants (single source).
"""
import datetime as dt
from pathlib import Path

import numpy as np
import pytest

from pipeline.telluric import gdas_fetch as gf

_ESO_GDAS = gf._eso_gdas_dir()
_HAS_ESO = _ESO_GDAS is not None and (_ESO_GDAS / "gdas_profiles_C-70.4-24.6.tar.gz").exists()
needs_eso = pytest.mark.skipif(not _HAS_ESO, reason="ESO telluriccorr GDAS tarball not installed")

VESTA_NIGHT = dt.date(2022, 11, 23)


# ── nearest-3-hourly slot resolution ─────────────────────────────────────────
def test_nearest_3hourly_rounds_to_real_slot():
    # 01:00 UT → nearest 3-hourly slot is 00 UT (not 01, which is absent from GDAS)
    t = dt.datetime(2022, 11, 23, 1, 0)
    assert gf.nearest_3hourly(t).hour == 0
    # 02:00 UT rounds up to 03
    assert gf.nearest_3hourly(dt.datetime(2022, 11, 23, 2)).hour == 3
    # a date (no time) → 00 UT that day
    assert gf.nearest_3hourly(VESTA_NIGHT) == dt.datetime(2022, 11, 23, 0, 0)


def test_nearest_3hourly_from_mjd():
    from astropy.time import Time
    mjd = Time("2022-11-23T00:30:00", format="isot").mjd
    assert gf.nearest_3hourly(None, mjd=mjd).hour == 0


# ── the loud-fail contract (no silent standard-atm fallback) ──────────────────
def test_unknown_site_raises_not_silent():
    with pytest.raises(KeyError):
        gf.fetch_gdas("atlantis", night=VESTA_NIGHT)


def test_unavailable_profile_raises_gdas_unavailable(tmp_path):
    # a fully-specified but nonexistent site (no tarball, no staged NOAA pull) must
    # RAISE — never quietly return a standard atmosphere.
    with pytest.raises(gf.GDASUnavailable):
        gf.fetch_gdas("nowhere", night=VESTA_NIGHT, cache_dir=tmp_path,
                      lat=0.0, lon=0.0, gdas_loc="C-999.9-99.9")


def test_needs_night_or_mjd():
    with pytest.raises(ValueError):
        gf.fetch_gdas("paranal")


# ── ESO tarball backend (real per-night profile) ─────────────────────────────
@needs_eso
def test_fetch_real_vesta_profile(tmp_path):
    p = gf.fetch_gdas("paranal", night=VESTA_NIGHT, cache_dir=tmp_path)
    assert p.exists() and p.stat().st_size > 0
    assert p.suffix == ".fits" and "2022-11-23T00" in p.name
    # molecfit GDAS_PROF FITS columns
    from astropy.io import fits
    cols = {c.name for c in fits.open(p)[1].columns}
    assert {"press", "height", "temp", "relhum"} <= cols


@needs_eso
def test_cache_hit_second_call(tmp_path):
    p1 = gf.fetch_gdas("paranal", night=VESTA_NIGHT, cache_dir=tmp_path)
    mtime = p1.stat().st_mtime_ns
    p2 = gf.fetch_gdas("paranal", night=VESTA_NIGHT, cache_dir=tmp_path)  # cache hit
    assert p2 == p1 and p2.stat().st_mtime_ns == mtime


# ── ASCII → FITS conversion (m → km; CFITSIO-readable bintable) ───────────────
def test_ascii_to_fits_units_and_columns(tmp_path):
    from astropy.io import fits
    src = tmp_path / "prof.gdas"
    src.write_text("# P HGT T RELHUM\n1000.0 100.0 290.0 40.0\n900.0 1000.0 280.0 35.0\n")
    out = gf._gdas_ascii_to_fits(src, tmp_path / "prof.gdas.fits")
    d = fits.open(out)[1].data
    assert list(d["press"]) == [1000.0, 900.0]
    assert list(d["height"]) == [0.1, 1.0]          # metres → km
    assert list(d["temp"]) == [290.0, 280.0]


def test_ascii_to_fits_empty_raises(tmp_path):
    src = tmp_path / "empty.gdas"
    src.write_text("# header only\n")
    with pytest.raises(ValueError):
        gf._gdas_ascii_to_fits(src, tmp_path / "out.fits")
