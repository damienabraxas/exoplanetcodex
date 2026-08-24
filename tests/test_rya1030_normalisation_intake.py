"""RYA-1030: normalisation is DETERMINED from the flux at intake, never asserted.

Synthetic fixtures, deliberately. The real controls that fixed every constant in
`normalization_intake` -- KP1984's two columns at identical wavelengths, KP2005's
`irradthu` vs its shipped `irradrelwl` residual, and a 55-window scan of the atlas --
live on data drives that are not mounted everywhere. A policy test that passes or fails
on whether a drive is mounted is an environment check wearing a policy check's label.
The measured numbers are recorded beside the constants; these fixtures reproduce each
class's SHAPE.
"""
import numpy as np
import pytest

from pipeline.normalization_intake import (
    CONTINUUM_LIMITED_BLUE_EDGE_A, ENVELOPE_SLOPE_MAX, NORMALISED, UNKNOWN,
    UNITY_TOLERANCE, UN_NORMALISED, NormalisationStateMismatch, cross_check, detect)

_RNG = np.random.default_rng(1030)

#: Comfortably inside the validated domain: above the blue edge, clear of every
#: registered telluric band.
CLEAN_LO, CLEAN_HI = 5000.0, 5100.0


def _spectrum(level, slope=0.0, n=20000, lo=CLEAN_LO, hi=CLEAN_HI, zeros=0):
    """An absorption spectrum whose upper envelope sits at `level` and tilts by `slope`.

    The tilt is CENTRED on the band so it does not also move the median level -- without
    that, a slope large enough to fail the flatness bound would drag the level out of
    tolerance too and the level test would silently take the credit for both.
    """
    w = np.linspace(lo, hi, n)
    x = (w - lo) / (hi - lo) - 0.5
    envelope = level * (1.0 + slope * x)
    depth = np.where(_RNG.random(n) < 0.04, _RNG.uniform(0.1, 0.8, n), 0.0)
    f = envelope * (1.0 - depth)
    if zeros:
        f[:zeros] = 0.0
    return w, f


# ── the two tests, and why neither is sufficient ─────────────────────────────

def test_detects_a_normalised_product():
    """KP1984 col1 measured 0.9800; KP2005's shipped residual 0.9883."""
    assert detect(*_spectrum(1.0)).value == NORMALISED


def test_detects_absolute_flux():
    """KP1984 col2 -- the PAIRED control -- measured 213.2875."""
    assert detect(*_spectrum(213.0, slope=0.02)).value == UN_NORMALISED


def test_the_LEVEL_test_catches_what_the_slope_test_cannot():
    """KP2005 `irradthu` as served sits at 2.14 while its slope over 100 A is only
    -0.044 -- INSIDE the flatness bound. Only the level test catches it."""
    ev = detect(*_spectrum(2.14, slope=-0.02))
    assert abs(ev.slope_across_band) <= ENVELOPE_SLOPE_MAX   # slope test would pass it
    assert ev.value == UN_NORMALISED                         # level test does not


def test_the_SLOPE_test_catches_what_the_level_test_cannot():
    """A blaze or SED riding on a product that happens to sit near unity."""
    ev = detect(*_spectrum(1.0, slope=0.30))
    assert abs(ev.median_p95 - 1.0) <= UNITY_TOLERANCE       # level test would pass it
    assert ev.value == UNKNOWN                               # not waved through
    assert "TILTED" in ev.reason


# ── the regimes where the envelope means nothing, and the module says so ─────

def test_below_the_blue_edge_is_UNKNOWN_not_an_accusation():
    """🔴 THE CONSTANT THAT WAS WRONG. Near-UV line blanketing depresses the 95th
    percentile of EVERY window whatever the normalisation: scanned across its range, the
    KNOWN-NORMALISED KP1984 atlas failed its own test at every probe below 4200 A. The
    module must decline to answer there rather than accuse a good product."""
    lo = CONTINUUM_LIMITED_BLUE_EDGE_A - 120.0
    ev = detect(*_spectrum(0.86, lo=lo, hi=lo + 100.0))
    assert ev.value == UNKNOWN
    assert "blue edge" in ev.reason


def test_a_telluric_band_is_excluded_from_the_envelope():
    """Saturated atmospheric absorption pushes the envelope down for a reason that has
    NOTHING to do with normalisation, and on a telluric-retaining product it fakes the
    un-normalised signature -- KP1984 failed at 9300 and 9500 A, both inside the
    registered H2O complex, and nowhere else between 4500 and 9900."""
    from pipeline.telluric_policy import TELLURIC_BANDS
    lo, hi = float(TELLURIC_BANDS[2][0]), float(TELLURIC_BANDS[2][1])   # O2 A-band
    ev = detect(*_spectrum(1.0, lo=lo - 2.0, hi=hi + 2.0))
    assert ev.n_telluric_windows > 0
    # what remains outside the band is too little to describe an envelope
    assert ev.value == UNKNOWN


def test_the_band_list_is_CALLED_not_re_enumerated():
    """A second copy of the authoritative telluric set is the RYA-845 defect shape."""
    import inspect

    from pipeline import normalization_intake
    src = inspect.getsource(normalization_intake)
    assert "from pipeline.telluric_policy import TELLURIC_BANDS" in src
    for edge in ("7594", "6867", "9280"):
        assert edge not in src, f"{edge} re-enumerated instead of imported"


def test_a_fill_value_is_not_flux():
    """🔴 FOUND ON REAL DATA. `solar_harps` returns EXACTLY 0.000 across eight
    consecutive 2 A windows at 5321-5337 A. Counting those zeros as flux made the
    envelope ramp 0.00 -> 0.97 and this module reported a BLAZE -- the right verdict
    with the WRONG DIAGNOSIS, which sends the reader after a defect that does not exist.
    """
    w, f = _spectrum(1.0, zeros=6000)
    ev = detect(w, f)
    assert ev.n_dropped >= 6000
    assert ev.value == NORMALISED, ev.citation()


def test_too_few_windows_is_UNKNOWN_never_a_default():
    """`unknown` is a real answer, not a failure to try -- the third value is
    load-bearing exactly as in `telluric_intake` (RYA-833)."""
    ev = detect(*_spectrum(1.0, n=400, lo=CLEAN_LO, hi=CLEAN_LO + 4.0))
    assert ev.value == UNKNOWN


# ── the cross-check rule ─────────────────────────────────────────────────────

def test_agreement_passes_and_returns_its_evidence():
    ev = cross_check("x", *_spectrum(1.0), declared=True)
    assert ev.value == NORMALISED
    assert "median rolling-P95" in ev.citation()


def test_flux_disagreeing_with_the_flag_is_a_LOUD_STOP():
    """A declared flag and a MIS-ROUTED file agree with each other perfectly and are both
    wrong. The detector informs the declared state; it must never auto-fix."""
    with pytest.raises(NormalisationStateMismatch, match="THE FLUX AND THE FLAG DISAGREE"):
        cross_check("solar_kpno_kurucz2005_corrected", *_spectrum(2.14),
                    declared=True, where="test")


def test_an_undeclared_product_must_be_scanned_before_a_continuum_stage():
    """`declared=None` is not permission to proceed -- it is the case the ticket exists
    for. The raise carries the detected state, which is what is needed to declare it."""
    with pytest.raises(NormalisationStateMismatch, match="UNDECLARED"):
        cross_check("some_new_arm", *_spectrum(1.0), declared=None)


def test_unknown_does_not_raise():
    """Turning 'the data could not speak' into a hard failure would make every short or
    blanketed window a blocker."""
    ev = cross_check("x", *_spectrum(1.0, n=400, lo=CLEAN_LO, hi=CLEAN_LO + 4.0),
                     declared=True)
    assert ev.value == UNKNOWN


def test_the_tolerance_sits_in_the_measured_gap_not_against_a_wall():
    """Derived, not chosen: the normalised population reaches 0.032 from unity and the
    nearest un-normalised case 0.696, with nothing between. A cut pressed against either
    wall misclassifies the next product."""
    assert 0.032 < UNITY_TOLERANCE < 0.696
    assert detect(*_spectrum(1.0 - 0.032)).value == NORMALISED
    assert detect(*_spectrum(1.0 + 0.696)).value == UN_NORMALISED
