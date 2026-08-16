"""
tests/test_nearuv_sigma_delta_rya846.py — RYA-846
=================================================
σ_δ is the second half of the near-UV continuum term (RYA-841 measured the first half, the
lever). What is pinned here is mostly the reasoning that keeps the number honest, because
every way this measurement could have gone wrong goes wrong QUIETLY:

  * the obvious route — the atlas's own overlapping segments — returns a spuriously perfect
    zero, because the overlaps are duplicated data;
  * a vacuum/air mix-up is 0.94 A here and reads as a normalisation difference;
  * Fortran overflow markers, if handed to a loader that drops bad rows, delete exactly the
    pixels where something went wrong;
  * and most of the apparent Kurucz-vs-Wallace disagreement is Wallace disagreeing with
    itself, which is only visible if you look.

The measurement needs both atlases and runs on Sirius; what runs here is the arithmetic,
the guards, and the reported numbers.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

RES = ROOT / "data" / "results" / "rya846"
SIGMA = RES / "rya846_sigma_delta.json"
RESPONDERS = RES / "rya846_responders.json"


@pytest.fixture(scope="module")
def sigma():
    if not SIGMA.exists():
        pytest.skip("sigma_delta summary absent (Sirius artifact)")
    return json.loads(SIGMA.read_text())


@pytest.fixture(scope="module")
def responders():
    if not RESPONDERS.exists():
        pytest.skip("responder summary absent")
    return json.loads(RESPONDERS.read_text())


# ── the controls, which are the whole reason to believe the number ────────────

def test_the_wavelength_convention_control_is_armed_and_passed(sigma):
    """Wallace ships wavenumbers => 1e8/nu is VACUUM; Kitt Peak is AIR. That is 0.94 A at
    3300 A, wider than a line, and it would masquerade as a normalisation difference. The
    control cross-correlates and the run aborts if it is off."""
    c = sigma["controls"]
    assert abs(c["residual_shift_A"]) <= 0.05, (
        "the two atlases are no longer aligned — a vacuum/air regression would look "
        "exactly like this")
    assert c["spectral_correlation"] > 0.90, (
        "the atlases no longer correlate; a ratio between them would be meaningless")


def test_the_wallace_internal_control_exists_and_is_subtracted(sigma):
    """THE CONTROL THAT HALVES THE ANSWER. Wallace's own two regions overlap and disagree
    with each other by rms 4.07%, against a Kurucz-vs-Wallace 4.94% at the lines. Without
    subtracting it the term would be quoted at 0.125 dex when the part attributable to a
    genuine reduction difference is 0.071."""
    dec = sigma["decomposition"]
    assert dec is not None
    assert dec["wallace_internal"] > 0.02, (
        "Wallace's internal scatter vanished — if that is real the net term grows, so "
        "re-derive rather than accepting a tighter number")
    assert dec["net_quadrature"] < dec["raw_at_lines"], (
        "the net figure must be smaller than the raw one; if not, the subtraction is "
        "no longer doing what it claims")
    assert dec["term_net_dex"] < dec["term_raw_dex"]


def test_the_dead_route_is_recorded(sigma):
    """The atlas's own segment overlaps are byte-identical duplicates and would have
    returned sigma_delta ~ 0. A negative result that flatters the product is exactly the
    one that gets quietly re-discovered, so it ships with the artifact."""
    assert "0.000e+00" in sigma["dead_route"] or "identical" in sigma["dead_route"].lower()


# ── the measurement ───────────────────────────────────────────────────────────

def test_all_forty_product_lines_are_covered(sigma):
    """Region 6 alone reaches 23 of 40; the composite reaches all of them. A sigma_delta
    quoted over half the lines would not be the product's sigma_delta."""
    at = sigma["sigma_delta_at_product_lines"]
    assert at["available"]
    assert at["n_lines"] == at["n_lines_total"] == 40
    assert at["n_lines_outside_wallace"] == 0
    assert sigma["band_fraction_compared"] == pytest.approx(1.0, abs=1e-6)


def test_sigma_delta_is_percent_scale_not_dex(sigma):
    """The original error was equating a fractional flux with a dex. sigma_delta is a
    FRACTION; the term is what carries dex, and only after multiplying by the lever."""
    sd = sigma["sigma_delta"]
    for v in sd.values():
        assert 0.0 < v < 0.5, "sigma_delta is a fraction of the continuum, not a dex"


def test_the_term_uses_the_linear_only_lever(sigma):
    """12 of 40 lines have no derivative. Averaging a slope over lines that have no slope
    is not a measurement, so the lever is the median over the 28 that do."""
    lv = sigma["lever"]
    assert lv["value"] == pytest.approx(2.54, abs=0.01)
    assert "+/-4%" in lv["source"] or "4%" in lv["source"], (
        "the lever's measured range must travel with it — RYA-841 forbids extrapolating "
        "it off the grid it was measured on")


def test_the_measured_term_brackets_the_assumed_one(sigma):
    """The honest result: 0.100 was never derived and its derivation was wrong, but it
    lands inside the measured range. That is worth pinning precisely because it would be
    easy to remember this ticket as having overturned the value."""
    dec = sigma["decomposition"]
    assert dec["term_net_dex"] < sigma["assumed_term_dex"] < dec["term_raw_dex"]


def test_the_revised_cells_move_only_the_systematic(sigma):
    """RYA-845 corrected these cells; this ticket only re-sizes one term inside them."""
    cells = sigma["revised_cells"]
    assert {c["treatment"] for c in cells} == {"1D-LTE", "1D-LTE-LABGF"}
    for c in cells:
        assert c["syst_with_measured_net"] < c["published_syst_dex"] \
            < c["syst_with_measured_raw"], (
            f"{c['treatment']}: the published cell no longer sits between the net and raw "
            f"revisions — re-derive before quoting either")


# ── the responders ────────────────────────────────────────────────────────────

def test_the_non_linear_responders_are_not_a_continuum_artifact(responders):
    """+0.467 dex at p=0.003 with identical line strength. A continuum misplacement cannot
    do that selectively to 12 of 40 lines."""
    c = responders["comparisons"]["a_unperturbed"]
    assert c["difference"] > 0.3
    assert c["permutation_p"] < 0.01
    depth = responders["comparisons"].get("theo_depth")
    if depth:
        assert depth["permutation_p"] > 0.05, (
            "the two populations now differ in line strength too, which would give the "
            "continuum a way to produce the split after all")


def test_the_mechanism_is_not_claimed(responders):
    """An earlier draft asserted unmodelled opacity. The blend metrics do not support it
    (core dominance p=0.26), and the write-up must not out-run the evidence."""
    v = responders["verdict"].lower()
    assert "not established" in v or "mechanism" in v
    cd = responders["comparisons"].get("depth_frac_core")
    if cd:
        assert cd["permutation_p"] > 0.05


def test_the_responders_are_carried_not_dropped(responders):
    """RYA-711/777/844. Dropping them would move the product 7.488 -> 7.413 and the
    scatter 0.413 -> 0.347 — a tighter bar bought by removing lines, which needs a stated
    per-line physical reason this ticket does not have."""
    d = responders["if_dropped"]
    assert d["median_linear_only"] < d["median_all"]
    assert d["scatter_linear_only"] < d["scatter_all"]
    assert len(responders["non_linear_lines"]) == responders["n_non_linear"] == 12


def test_every_line_keeps_its_row(responders):
    csv = RES / "rya846_responder_diagnosis.csv"
    if not csv.exists():
        pytest.skip("per-line diagnosis absent")
    assert len(pd.read_csv(csv)) == responders["n_lines"] == 40
