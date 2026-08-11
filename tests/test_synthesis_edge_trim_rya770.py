"""
tests/test_synthesis_edge_trim_rya770.py
========================================
RYA-770 — the SynthesisHandler must not compare observed flux against iSpec's zeroed
synthesis edges.

THE DEFECT THIS PINS
--------------------
`ispec.generate_spectrum` returns exactly 0.0 for the first synthesis pixel and the
last 3-4 of every window. `_fit_synth_flux` synthesises on a fine grid and interpolates
onto the observed pixels, so an observed pixel landing in a zeroed span is compared
against ~0 instead of ~1. With the fit's sigma = 0.01 that single pixel contributes
(1/0.01)^2 = 1e4 to chi2 and dominates a few-dozen-pixel window.

Measured on the 12 optical control lines: median red_chi2 92.17 as-is vs 2.77 with only
those pixels removed, against a banked med red_chi2 of 2.976. Individual lines moved
301.59 -> 0.84, 162.71 -> 2.72, 16.58 -> 1.80.

These tests exercise the trim GEOMETRY, which is pure and needs no synthesiser. The
end-to-end proof is the optical control run on Sirius, recorded on the ticket; a unit
test cannot stand in for that, and does not pretend to.
"""
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from pipeline.measure.synthesis import (          # noqa: E402
    clean_span, edge_pixels_to_drop, EDGE_ZERO_FLUX)


def _synthetic_with_zeroed_edges(n=200, n_front=1, n_back=4):
    """A synthetic flux array shaped like iSpec's: real in the middle, 0.0 at the ends."""
    f = np.full(n, 0.99)
    f[:n_front] = 0.0
    if n_back:
        f[-n_back:] = 0.0
    return f


def test_clean_span_finds_the_real_interior():
    sw = np.linspace(500.0, 500.06, 200)
    f = _synthetic_with_zeroed_edges(200, 1, 4)
    lo, hi, n_zero = clean_span(sw, f)
    assert n_zero == 5
    assert lo == pytest.approx(sw[1])
    assert hi == pytest.approx(sw[-5])


def test_clean_span_is_the_intersection_over_all_trials():
    """A pixel real at one trial abundance and zeroed at another is NOT clean: the fit
    evaluates many trials and any of them can produce the spike."""
    sw = np.linspace(500.0, 500.06, 200)
    a = _synthetic_with_zeroed_edges(200, 1, 4)
    b = _synthetic_with_zeroed_edges(200, 3, 2)      # different edge counts
    lo, hi, _ = clean_span(sw, a, b)
    assert lo == pytest.approx(sw[3])                # widest leading zero wins
    assert hi == pytest.approx(sw[-5])               # widest trailing zero wins


def test_clean_span_reports_all_zero_rather_than_guessing():
    sw = np.linspace(500.0, 500.06, 50)
    lo, hi, n_zero = clean_span(sw, np.zeros(50))
    assert (lo, hi) == (None, None) and n_zero == 50


def test_the_sentinel_is_distinguished_from_a_deep_line_core():
    """A saturated core is small but REAL and must never be trimmed as an edge."""
    sw = np.linspace(500.0, 500.06, 200)
    f = _synthetic_with_zeroed_edges(200, 1, 4)
    f[100] = 0.02                                    # a very deep line core
    lo, hi, n_zero = clean_span(sw, f)
    assert n_zero == 5                               # the core is not counted
    assert lo < sw[100] < hi                         # and it stays inside the span


def test_only_in_window_pixels_are_dropped():
    """`_fit_synth_flux` masks out-of-window pixels itself; dropping them here would
    change nothing and would inflate the reported count."""
    obs = np.array([499.90, 499.995, 500.00, 500.03, 500.06, 500.065, 500.20])
    drop = edge_pixels_to_drop(obs, wave_base=500.0, wave_top=500.06,
                               clean_lo=500.005, clean_hi=500.055)
    assert list(drop) == [False, False, True, False, True, False, False]


def test_no_pixels_dropped_when_the_whole_window_is_clean():
    obs = np.linspace(500.0, 500.06, 25)
    drop = edge_pixels_to_drop(obs, 500.0, 500.06, clean_lo=500.0, clean_hi=500.06)
    assert not drop.any()


def test_a_pixel_exactly_on_the_clean_boundary_is_kept():
    """The boundary sample IS a real synthesised value; excluding it would throw away
    good data and make the trim depend on floating-point luck."""
    obs = np.array([500.005, 500.055])
    drop = edge_pixels_to_drop(obs, 500.0, 500.06, clean_lo=500.005, clean_hi=500.055)
    assert not drop.any()


def test_trim_removes_exactly_the_pixels_that_would_read_as_zero():
    """End-to-end on the geometry: interpolating the zeroed grid onto the surviving
    observed pixels must never return the sentinel."""
    from scipy.interpolate import interp1d
    sw = np.linspace(500.0, 500.06, 300)
    f = _synthetic_with_zeroed_edges(300, 1, 4)
    obs = np.linspace(500.0, 500.06, 40)             # samples the zeroed ends

    before = interp1d(sw, f, bounds_error=False, fill_value=1.0)(obs)
    assert (before < EDGE_ZERO_FLUX).any(), "fixture must actually hit the zeroed span"

    lo, hi, _ = clean_span(sw, f)
    keep = ~edge_pixels_to_drop(obs, 500.0, 500.06, lo, hi)
    after = interp1d(sw, f, bounds_error=False, fill_value=1.0)(obs[keep])
    assert not (after < EDGE_ZERO_FLUX).any()
    assert keep.sum() >= 5                           # still enough pixels to fit


# ── RYA-342, re-pinned: the EW saturation ceiling is an EW-PATH concept ──────────
# It came back through the data container: `LineMeasurement.__post_init__` applied the
# REW ceiling to every line whatever produced it, so a flux-space synthesis fit of a
# strong line was quarantined as "saturated" — the exact class of line flux-space
# synthesis exists to recover. Measured on the banked synth-v2 lines: the handler fitted
# 5 of 6 at red_chi2 1.7-8.8, median A(Fe I) = 7.529 against a banked 7.520, and the
# ceiling excluded all but one, so the control reported "1 line, +0.181 dex".
_STRONG = dict(element="Fe", ion="I", wavelength_air_A=5506.779,
               instrument="harps", ew_mA=151.63, ew_method="x")


def test_ew_path_still_gates_saturated_lines():
    """The ceiling must keep working where it is correct — this is not a relaxation."""
    from pipeline.band_products import LineMeasurement, REW_SATURATION_CEILING
    lm = LineMeasurement(**_STRONG, abundance=7.4)
    assert lm.rew > REW_SATURATION_CEILING
    assert lm.in_aggregate is False
    assert "saturation" in lm.excluded_reason


def test_flux_space_fit_is_not_gated_on_the_ew_saturation_ceiling():
    from pipeline.band_products import LineMeasurement
    lm = LineMeasurement(**_STRONG, abundance=7.4, ew_inversion=False)
    assert lm.in_aggregate is True
    assert lm.excluded_reason == ""


def test_rew_is_still_recorded_for_synthesis_lines():
    """Lifting the gate must not lose the diagnostic — it is reported, just not acted on."""
    from pipeline.band_products import LineMeasurement
    lm = LineMeasurement(**_STRONG, abundance=7.4, ew_inversion=False)
    assert lm.rew is not None and lm.rew == pytest.approx(-4.56, abs=0.01)


def test_the_synthesis_handler_opts_out_at_every_construction_site():
    """Both the accept and the quarantine path must opt out, or a strong line is still
    silently dropped on one of them."""
    src = (ROOT / "pipeline" / "measure" / "synthesis.py").read_text()
    assert src.count("LineMeasurement(") == src.count("ew_inversion=False"), (
        "every LineMeasurement built by the synthesis handler must set "
        "ew_inversion=False")


def test_ew_domain_callers_are_unchanged_by_default():
    """The opt-out defaults to the EW behaviour, so no other handler moves."""
    from pipeline.band_products import LineMeasurement
    import inspect
    sig = inspect.signature(LineMeasurement)
    assert sig.parameters["ew_inversion"].default is True


def test_adapter_takes_the_synthesis_step_from_the_core_not_a_local_copy():
    """The trim is only correct if it rebuilds the grid the FIT uses. A second literal
    in this file would drift silently (the RYA-768 two-sources-of-truth shape)."""
    import pipeline.measure.synthesis as S
    assert not hasattr(S, "SYNTH_WSTEP_NM"), "adapter must not define its own step"
    from pipeline.abundances_derive import SYNTH_FIT_WSTEP_NM
    assert SYNTH_FIT_WSTEP_NM == 0.0002
