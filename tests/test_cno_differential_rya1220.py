"""RYA-1220 WP D -- the target/Sun CNO differential.

Pins the finding that the 2.404 dex CN_red difference is a broken SOLAR fit rather than
a stellar one, and that carbon is deliverable from data already in hand.
"""
from __future__ import annotations

import json
import pathlib

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
ART = ROOT / "data/output/rya1220/cno_target_sun_differential.json"

pytestmark = pytest.mark.skipif(not ART.exists(), reason="differential not generated")


def _doc():
    return json.loads(ART.read_text())


def test_carbon_is_measured_from_four_matched_indicators():
    c = _doc()["by_element"]["C"]
    assert c["state"] == "MEASURED"
    assert c["n_usable"] == 4
    assert sorted(c["indicators"]) == ["C2_Swan", "CH_Gband", "CI_5052", "CI_5380"]
    assert c["median_of_differences_dex"] == pytest.approx(-0.057, abs=1e-3)


def test_the_cn_red_difference_is_a_broken_solar_fit_not_a_star():
    pair = {p["key"]: p for p in _doc()["pairs"]}["CN_red"]
    assert pair["delta_dex"] == pytest.approx(2.404, abs=1e-3)
    assert pair["usable"] is False
    # the solar end is what is broken
    assert pair["solar_sigma_fit"] > 0.8
    assert pair["solar_frac_rise_weaker"] < 1e-4
    assert pair["target_sigma_fit"] < 0.1
    assert any("SOLAR end is unconstrained" in b for b in pair["blockers"])


def test_a_rejected_diagnostic_cannot_carry_a_differential():
    pair = {p["key"]: p for p in _doc()["pairs"]}["CN_red"]
    assert any("rejected as an abundance route" in b for b in pair["blockers"])


def test_oxygen_is_blocked_on_the_solar_side_too():
    o = _doc()["by_element"]["O"]
    assert o["state"] == "BLOCKED"
    assert o["blocked_indicators"] == ["OI_6300"]


def test_the_differential_is_median_of_differences_not_difference_of_medians():
    """Pairing is what makes shared terms cancel; aggregating first discards it."""
    src = (ROOT / "scripts/rya1220_cno_differential.py").read_text()
    assert "median of per-indicator DIFFERENCES" in src
    c = _doc()["by_element"]["C"]
    assert "median_of_differences_dex" in c and "difference_of_medians_dex" in c


def test_spread_is_carried_as_the_honest_uncertainty():
    c = _doc()["by_element"]["C"]
    assert c["spread_dex"] == pytest.approx(0.318, abs=1e-3)
    assert "Jacobians that would let it be decomposed are still owed" in c["note"]
