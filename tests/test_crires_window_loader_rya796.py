"""
tests/test_crires_window_loader_rya796.py
=========================================
RYA-796 — the CRIRES+ Vesta arm of `measure_band_ew.load_window`.

WHAT IS WORTH PINNING
---------------------
The 18 staged IDPs are Sirius-only, so the science numbers belong on the ticket. What
these tests hold down is the set of ways this loader could return a plausible-looking
WRONG spectrum, each of which has a precedent in this repo:

  * the UNIT. `WAVE` is nanometres; read as Angstrom it yields 949-2485 and a confident
    wrong answer. RYA-794 hit the same trap from the other direction.
  * the FRAME. SPECSYS is TOPOCENT and Vesta is a moving reflector, so a wavelength
    measured off a raw frame is wrong by tens of km/s. RYA-643 was exactly this defect,
    in four copies, and it is silent.
  * the CO-ADD. Two epochs of one setting can be the same face of Vesta or the opposite
    one; averaging the second kind is not signal-to-noise.
  * the FALLBACK. An unknown instrument must still raise rather than quietly serving
    another arm's data (RYA-786).

The fixture is a synthetic IDP built to the real files' structure, so all of the above
is exercised without the 37 MB of Sirius-only data.
"""
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

fits = pytest.importorskip("astropy.io.fits")

import measure_band_ew as M  # noqa: E402


# ── a synthetic IDP, structured like the real ones ───────────────────────────

def _write_idp(path: Path, *, setting="Y1029", mjd=59904.99672, snr=100.0,
               lo_nm=1000.0, hi_nm=1010.0, n=4000, specsys="TOPOCENT",
               tunit="nm", bad_slice=None, flux_level=1.0e5):
    wave = np.linspace(lo_nm, hi_nm, n)
    flux = np.full(n, flux_level)
    qual = np.zeros(n, dtype=np.int32)
    if bad_slice is not None:
        qual[bad_slice] = 1                      # a real detector gap: QUAL != 0
    cols = fits.ColDefs([
        fits.Column(name="WAVE", format=f"{n}D", array=[wave]),
        fits.Column(name="FLUX", format=f"{n}D", array=[flux]),
        fits.Column(name="ERR", format=f"{n}D", array=[np.full(n, 10.0)]),
        fits.Column(name="QUAL", format=f"{n}J", array=[qual]),
    ])
    t = fits.BinTableHDU.from_columns(cols, name="SPECTRUM")
    t.header["TUNIT1"] = tunit
    t.header["TUNIT2"] = "adu"
    ph = fits.PrimaryHDU()
    ph.header["SPECSYS"] = specsys
    ph.header["MJD-OBS"] = mjd
    ph.header["SNR"] = snr
    ph.header["ESO INS WLEN ID"] = setting
    ph.header["ESO TEL TARG RADVEL"] = 0.0
    ph.header["INSTRUME"] = "CRIRES"
    fits.HDUList([ph, t]).writeto(path, overwrite=True)


@pytest.fixture
def staged(tmp_path, monkeypatch):
    """A staged CRIRES+ directory the loader will resolve, with the caches cleared."""
    d = tmp_path / "CRIRESPlus"
    d.mkdir()
    monkeypatch.setattr(M, "_CRIRES_CANDIDATES", (str(d),))
    monkeypatch.setattr(M, "_CRIRES_REST_CANDIDATES", (str(tmp_path / "none"),))
    M._crires_cache.clear()
    yield d
    M._crires_cache.clear()


# ── the unit ─────────────────────────────────────────────────────────────────

def test_nm_is_converted_to_angstrom(staged):
    """1000-1010 nm must present as 10000-10100 A, not 1000-1010."""
    _write_idp(staged / "a.fits", lo_nm=1000.0, hi_nm=1010.0)
    fr = M.crires_frames()[0]
    assert fr["lo"] == pytest.approx(10000.0)
    assert fr["hi"] == pytest.approx(10100.0)


def test_an_angstrom_labelled_file_is_not_scaled_twice(staged):
    _write_idp(staged / "a.fits", lo_nm=10000.0, hi_nm=10100.0, tunit="angstrom")
    fr = M.crires_frames()[0]
    assert fr["lo"] == pytest.approx(10000.0)


def test_an_unknown_wavelength_unit_stops_rather_than_guessing(staged):
    _write_idp(staged / "a.fits", tunit="furlongs")
    with pytest.raises(LookupError, match="unexpected TUNIT1"):
        M.crires_frames()


def test_a_unit_error_landing_out_of_band_is_caught(staged):
    """The belt to the unit braces: if conversion put the arm somewhere CRIRES+ cannot
    observe, that is a unit bug and must not be measured off."""
    _write_idp(staged / "a.fits", lo_nm=100.0, hi_nm=101.0, tunit="angstrom")
    with pytest.raises(LookupError, match="outside the instrument's catalogued"):
        M.crires_frames()


# ── the frame: a topocentric spectrum must not reach a measurement ───────────

def test_raw_topocentric_frames_are_refused(staged):
    _write_idp(staged / "a.fits", lo_nm=1000.0, hi_nm=1010.0)
    with pytest.raises(M.RestFrameNotConditioned) as e:
        M.load_crires_window(10050.0, 1.0)
    msg = str(e.value)
    assert "TOPOCENT" in msg and "reflected_solar_rv" in msg


def test_the_conditioning_leg_can_opt_in_explicitly(staged):
    """The leg that MEASURES the shift has to read the raw frame. It says so, and the
    provenance string carries the warning onward."""
    _write_idp(staged / "a.fits", lo_nm=1000.0, hi_nm=1010.0)
    w, f, prov = M.load_crires_window(10050.0, 1.0, allow_topocentric=True)
    assert len(w) > M.MIN_WINDOW_POINTS and len(w) == len(f)
    assert "RAW TOPOCENTRIC" in prov


def test_refusal_is_a_lookuperror_so_drivers_still_catch_it(staged):
    _write_idp(staged / "a.fits")
    assert issubclass(M.RestFrameNotConditioned, LookupError)


# ── coverage is measured, not declared ───────────────────────────────────────

def test_flagged_pixels_do_not_count_as_coverage(staged):
    """QUAL != 0 is a detector gap. A line inside one is not covered, however much the
    setting's declared span suggests otherwise (RYA-377)."""
    _write_idp(staged / "a.fits", lo_nm=1000.0, hi_nm=1010.0, n=4000,
               bad_slice=slice(1500, 2500))
    fr = M.crires_frames()[0]
    assert len(fr["wave_A"]) == 3000
    with pytest.raises(LookupError):
        M.load_crires_window(10050.0, 0.2, allow_topocentric=True)


def test_a_wavelength_no_setting_covers_says_so(staged):
    _write_idp(staged / "a.fits", lo_nm=1000.0, hi_nm=1010.0)
    with pytest.raises(LookupError, match="no CRIRES\\+ setting covers"):
        M.load_crires_window(15000.0, 1.0, allow_topocentric=True)


# ── the co-add is a rotation question ────────────────────────────────────────

def test_same_rotation_phase_epochs_are_coadded(staged):
    """28 minutes apart is 32 degrees of sub-observer longitude -- the same face."""
    _write_idp(staged / "a.fits", mjd=59904.99672, snr=86.3)
    _write_idp(staged / "b.fits", mjd=59905.01640, snr=121.1)
    w, f, prov = M.load_crires_window(10050.0, 1.0, allow_topocentric=True)
    assert "co-added 2 epochs" in prov


def test_opposite_hemisphere_epochs_are_not_coadded(staged):
    """~23.8 h apart is 166 degrees: the other side of Vesta. RYA-372 wrote its
    conditioned frames per-frame precisely because of this."""
    _write_idp(staged / "a.fits", setting="H1559", mjd=59905.01174, snr=98.8)
    _write_idp(staged / "b.fits", setting="H1559", mjd=59906.00459, snr=124.3)
    w, f, prov = M.load_crires_window(10050.0, 1.0, allow_topocentric=True)
    assert "co-added" not in prov
    assert "b.fits" in prov            # the higher-SNR epoch wins, alone


def test_settings_are_never_merged_with_each_other(staged):
    """Different settings have different blaze and continuum, and none of this is
    normalised -- merging them manufactures a continuum (RYA-794 Step 2)."""
    _write_idp(staged / "a.fits", setting="Y1028", mjd=59906.99788, snr=115.7)
    _write_idp(staged / "b.fits", setting="Y1029", mjd=59904.99672, snr=86.3)
    w, f, prov = M.load_crires_window(10050.0, 1.0, allow_topocentric=True)
    assert "not merged" in prov
    assert prov.count("+") == 0 or "co-added" not in prov


def test_rotation_clustering_is_pure_and_uses_the_real_period():
    """The clustering rule on its own, no FITS involved."""
    P = M.CRIRES_ROTATION_H
    close = [dict(mjd=0.0, snr=1), dict(mjd=(0.05 * P) / 24.0, snr=1)]
    far = [dict(mjd=0.0, snr=1), dict(mjd=(0.5 * P) / 24.0, snr=1)]
    assert len(M._coaddable(close)) == 1
    assert len(M._coaddable(far)) == 2
    assert P == pytest.approx(5.342, abs=0.01)


# ── registration ─────────────────────────────────────────────────────────────

def test_crires_idp_is_not_pre_normalised():
    """adu, not residual flux. Getting this wrong is the RYA-713 continuum defect.

    RE-KEYED BY RYA-904. This asserted `PRE_NORMALISED["crires_plus"]`, and the
    instrument-level answer was the bug: the SAME instrument also serves the RYA-794
    corrected Y arm, which IS normalised. One key could only ever be right about one of
    them. The property being pinned is unchanged — the raw IDPs are adu — but it is now
    stated about the holding that has it.
    """
    assert M.PRE_NORMALISED["solar_vesta_crires_plus_idp"] is False
    assert M.PRE_NORMALISED["solar_crires_plus_y_rya794"] is True
    assert "crires_plus" not in M.PRE_NORMALISED, (
        "PRE_NORMALISED is holding-keyed now; an instrument key here would be a "
        "silent, always-wrong answer for one of the two CRIRES+ holdings (RYA-904)")


def test_unknown_instrument_still_raises(staged):
    with pytest.raises(LookupError, match="no window loader for instrument"):
        M.load_window("jwst_nirspec_prism", 10050.0, 1.0)


def test_load_window_routes_crires_and_returns_the_common_shape(staged):
    _write_idp(staged / "a.fits", lo_nm=1000.0, hi_nm=1010.0)
    # TWO gates now stand in front of this arm, and the ORDER is the physics (RYA-806).
    # Tellurics are stationary in the topocentric frame, so the correction happens there
    # and the RV shift comes after (RYA-373) -- telluric is the earlier blocker, so it is
    # what a caller is told about first. Before RYA-806 this raised
    # RestFrameNotConditioned, because the telluric defect was not being checked at all.
    with pytest.raises(M.TelluricNotCorrected):
        M.load_window("crires_plus", 10050.0, 1.0)      # gated, by design
    # Clear the telluric gate the way the correction leg does, and the rest-frame gate
    # is still there underneath it -- neither refusal masks the other.
    with pytest.raises(M.RestFrameNotConditioned):
        M.load_window("crires_plus", 10050.0, 1.0, allow_uncorrected=True)
    w, f, prov = M.load_crires_window(10050.0, 1.0, allow_topocentric=True)
    assert isinstance(w, np.ndarray) and isinstance(f, np.ndarray)
    assert w.shape == f.shape and isinstance(prov, str)
    assert np.all(np.diff(w) >= 0)                      # ascending, like the other arms


def test_catalog_registers_the_traps_that_bite():
    from data.catalog.instruments import load_catalog
    row = next(r for r in load_catalog() if r["instrument_id"] == "crires_plus")
    assert row["telluric_basis"] == "correction_required"
    assert row["pipeline_ingest_status"] == "loader_ready"
    req = row["reduction_requirements"].lower()
    for token in ("nanometre", "topocent", "qual", "normalis"):
        assert token in req, f"reduction_requirements does not record {token!r}"
