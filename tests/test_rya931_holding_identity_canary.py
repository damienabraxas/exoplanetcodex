"""RYA-931/927 canary: corrected and uncorrected HARPS holdings cannot collapse.

Two products of the same star, from the same instrument, differing only by
whether the telluric correction was applied, are the easiest pair in the project
to confuse.  The failure modes are specific and each is checked here:

* a corrected product whose telluric state is *asserted* in the registry rather
  than *derived* from the product, so a registry typo silently relabels pixels;
* a module-level catalogue cache that answers for the holding loaded first;
* two holding ids that resolve to the same pixels.

The fixtures are deliberately distinguishable: the uncorrected one carries a
saturated absorption trough that the corrected one does not.
"""
from __future__ import annotations

import numpy as np
import pytest
from astropy.io import fits

from pipeline import telluric_intake
from pipeline.telluric_policy import TELLURIC_BANDS

O2_B = next((lo, hi) for lo, hi, name in TELLURIC_BANDS if "O2 B" in name)


def _harps_product(path, *, corrected: bool):
    wave = np.arange(6840.0, 6912.0, 0.01)
    flux = np.ones(wave.size)
    trough = (wave >= O2_B[0]) & (wave <= O2_B[1])
    transmission = np.ones(wave.size)
    transmission[trough] = 0.2 + 0.8 * np.abs(np.sin((wave[trough] - O2_B[0]) * 3.0))
    if not corrected:
        flux = flux * transmission

    header = fits.Header()
    header["INSTRUME"] = "HARPS"
    header["SPECSYS"] = "BARYCENT"
    header["ESO PRO REC1 ID"] = "harps_s1d"
    if corrected:
        header["TELLCORR"] = True
        header["TELLENG"] = "molecfit-4.4.4"
    n = wave.size
    hdus = [fits.PrimaryHDU(header=header),
            fits.BinTableHDU.from_columns([
                fits.Column(name="WAVE", format=f"{n}D", array=wave.reshape(1, -1)),
                fits.Column(name="FLUX", format=f"{n}D", array=flux.reshape(1, -1)),
            ])]
    if corrected:
        hdus.append(fits.ImageHDU(data=transmission, name="MTRANS"))
    fits.HDUList(hdus).writeto(path, overwrite=True)
    return path


@pytest.fixture
def pair(tmp_path):
    return (_harps_product(tmp_path / "uncorrected.fits", corrected=False),
            _harps_product(tmp_path / "corrected.fits", corrected=True))


def test_fixtures_are_actually_distinguishable(pair):
    """A canary that cannot tell the two apart proves nothing."""
    uncorrected, corrected = pair
    fu = fits.open(uncorrected)[1].data["FLUX"][0]
    fc = fits.open(corrected)[1].data["FLUX"][0]
    wave = fits.open(uncorrected)[1].data["WAVE"][0]
    band = (wave >= O2_B[0]) & (wave <= O2_B[1])
    assert fu[band].min() < 0.5 <= fc[band].min()
    assert not np.array_equal(fu, fc)


def test_state_is_derived_from_the_product_not_asserted(pair):
    """`telluric_applied` must come from the product's own evidence."""
    uncorrected, corrected = pair
    assert telluric_intake.from_headers(corrected).value == telluric_intake.APPLIED
    assert telluric_intake.from_headers(uncorrected).value == telluric_intake.NOT_APPLIED


def test_evidence_cites_the_transmission_not_a_keyword_alone(pair):
    """The strongest available evidence must be the persisted transmission."""
    _, corrected = pair
    evidence = telluric_intake.from_headers(corrected)
    assert "MTRANS" in evidence.transmission_hdus
    assert evidence.transmission_all_unity is False
    assert any("transmission" in r.lower() or "molecfit" in r.lower()
               for r in evidence.reasons), evidence.reasons


def test_an_all_unity_transmission_is_not_a_correction(tmp_path):
    """Computing a model and dividing by nothing is NOT-APPLIED, not APPLIED."""
    wave = np.arange(6840.0, 6912.0, 0.01)
    n = wave.size
    header = fits.Header()
    header["INSTRUME"] = "HARPS"
    header["ESO PRO REC1 ID"] = "harps_s1d"
    fits.HDUList([
        fits.PrimaryHDU(header=header),
        fits.BinTableHDU.from_columns([
            fits.Column(name="WAVE", format=f"{n}D", array=wave.reshape(1, -1)),
            fits.Column(name="FLUX", format=f"{n}D", array=np.ones(n).reshape(1, -1))]),
        fits.ImageHDU(data=np.ones(n), name="MTRANS"),
    ]).writeto(tmp_path / "unity.fits", overwrite=True)
    assert telluric_intake.from_headers(tmp_path / "unity.fits").value == \
        telluric_intake.NOT_APPLIED


def test_repeated_interleaved_reads_do_not_reuse_stale_state(pair):
    """A-B-A: the third read must reproduce the first, not echo the second."""
    uncorrected, corrected = pair
    first = telluric_intake.from_headers(uncorrected).value
    second = telluric_intake.from_headers(corrected).value
    third = telluric_intake.from_headers(uncorrected).value
    fourth = telluric_intake.from_headers(corrected).value
    assert first == third == telluric_intake.NOT_APPLIED
    assert second == fourth == telluric_intake.APPLIED
    assert telluric_intake.from_headers(uncorrected).path != \
        telluric_intake.from_headers(corrected).path


def test_from_many_refuses_to_average_a_mixed_folder(pair):
    """A folder holding both must not resolve to a single confident verdict."""
    uncorrected, corrected = pair
    verdict, evidence = telluric_intake.from_many([uncorrected, corrected])
    assert len(evidence) == 2
    assert verdict == telluric_intake.UNKNOWN, (
        "a mixed corrected/uncorrected folder must not collapse to one state")


def _normalised_csv(path, *, depth: float):
    """A minimal stand-in for a `spectra_normalize` product."""
    import pandas as pd
    wave = np.arange(6840.0, 6912.0, 0.01)
    flux = np.ones(wave.size)
    band = (wave >= O2_B[0]) & (wave <= O2_B[1])
    flux[band] = depth
    pd.DataFrame({"wavelength_air_A": wave, "flux_raw": flux * 1.0e5,
                  "continuum": np.full(wave.size, 1.0e5),
                  "flux_normalized": flux}).to_csv(path, index=False)
    return path


def test_the_harps_reader_is_keyed_by_product_not_by_instrument(tmp_path, monkeypatch):
    """RYA-904's defect shape: one instrument silently standing for one holding.

    HARPS now serves two holdings.  A single-slot cache would hand the second
    caller the first caller's pixels, and an unreachable or mis-served holding
    leaves no trace -- it reads exactly like having no data.
    """
    import sys
    sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parents[1] / "scripts"))
    monkeypatch.setenv("CODEX_KP_ATLAS", str(tmp_path))
    (tmp_path / "lm0296").touch()
    import measure_band_ew as M

    uncorrected = _normalised_csv(tmp_path / "a.csv", depth=0.02)
    corrected = _normalised_csv(tmp_path / "b.csv", depth=0.80)
    M._harps_cache.clear()

    a1 = M.harps_spectrum(uncorrected)[1]
    b = M.harps_spectrum(corrected)[1]
    a2 = M.harps_spectrum(uncorrected)[1]

    assert not np.array_equal(a1, b), "two products served identical pixels"
    assert np.array_equal(a1, a2), "the second read echoed the other product"
    assert {str(uncorrected), str(corrected)} <= set(M._harps_cache)


def test_both_harps_holdings_are_addressable_by_name():
    """An unaddressable holding reads to every caller as if it did not exist."""
    from pipeline.loader_coverage import addressable_holdings
    addressable = addressable_holdings()
    assert addressable.get("solar_harps") == "harps"
    assert addressable.get("solar_harps_molecfit_corrected") == "harps"


def test_the_corrected_holding_does_not_displace_the_uncorrected_one():
    """Wiring must not silently switch which product existing measurements use."""
    import sys
    sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parents[1] / "scripts"))
    import measure_band_ew as M
    order = [h.holding_id for h in M._INSTRUMENT_HOLDINGS["harps"]]
    assert order.index("solar_harps") < order.index("solar_harps_molecfit_corrected")
    assert M.PRE_NORMALISED["solar_harps_molecfit_corrected"] is True
