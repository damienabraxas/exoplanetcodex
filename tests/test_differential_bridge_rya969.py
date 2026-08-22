"""
RYA-969 — the differential bridge.

What these pin:

  * the chain rule is SHORTEST PATH TO THE SUN, not nearest neighbour — and the α Cen B / 55 Cnc
    pair is the case that separates the two;
  * a hop that cannot be bridged FAILS LOUDLY and names the remedy, never bridges silently;
  * matching is on wavelength AND excitation potential (RYA-780);
  * R1–R4: the bridge transports a difference and can never adjust a reference;
  * both products are always reported, and the zero-point does not shrink.
"""
from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline import differential_bridge as B  # noqa: E402
from config.constants import STAR_PARAMS  # noqa: E402


def _th(**o):
    t = B.BridgeThresholds(min_shared_lines=5, rew_min=-6.0, rew_max=-4.8)
    for k, v in o.items():
        setattr(t, k, v)
    return t


def _lines(n, *, a0=7.5, wl0=5000.0, ep=3.0, step=1.0, rew=-5.2):
    return [dict(wavelength_air_A=wl0 + step * i, ep_eV=ep, rew=rew,
                 abundance=a0 + 0.01 * ((i % 5) - 2)) for i in range(n)]


# ── the chain rule ───────────────────────────────────────────────────────────────────
def test_every_current_star_reaches_the_sun_directly():
    """Measured, not assumed: our largest gaps are far inside Jofré's failure regime."""
    p = {k: v for k, v in STAR_PARAMS.items() if k != "synthetic_no_logg"}
    for star in p:
        if star == B.SUN:
            continue
        chain = B.plan_chain(star, p)
        assert len(chain) == 1, f"{star} needed {len(chain)} hops; expected a direct hop"
        assert chain[0][1] == B.SUN


def test_the_chain_rule_is_not_nearest_neighbour():
    """🔴 α Cen B's nearest neighbour is 55 Cnc A, but its bridge reference is the Sun."""
    p = {k: v for k, v in STAR_PARAMS.items() if k != "synthetic_no_logg"}
    nearest = min((B.hop_distance(p["alpha_cen_b"], p[n]), n)
                  for n in p if n != "alpha_cen_b")[1]
    assert nearest == "55cnc_a", "the premise of this test has changed; recheck the geometry"
    assert B.plan_chain("alpha_cen_b", p)[0][1] == B.SUN, (
        "the chain followed the nearest neighbour and added a hop for no gain")


def test_a_hop_beyond_the_observed_regime_needs_an_intermediate():
    p = {"solar": STAR_PARAMS["solar"],
         "mid": dict(teff=4200.0, logg=4.7, feh_ref=0.0),
         "cool": dict(teff=3300.0, logg=4.9, feh_ref=0.0)}
    assert not B.within_known_limits(p["cool"], p["solar"])[0]
    chain = B.plan_chain("cool", p)
    assert [h[1] for h in chain] == ["mid", "solar"], "the two-step bridge did not engage"


def test_an_unreachable_star_fails_loudly_and_names_the_remedy():
    p = {"solar": STAR_PARAMS["solar"], "cool": dict(teff=3200.0, logg=5.0, feh_ref=0.0)}
    with pytest.raises(B.HopTooLarge, match="needs a new benchmark"):
        B.plan_chain("cool", p)


def test_a_star_with_no_parameters_cannot_be_planned():
    with pytest.raises(KeyError, match="no fundamental parameters"):
        B.plan_chain("tau_boo", {k: v for k, v in STAR_PARAMS.items()})


# ── the gate ─────────────────────────────────────────────────────────────────────────
def test_declared_or_nothing():
    for name in ("min_shared_lines", "rew_min", "rew_max"):
        with pytest.raises(B.ThresholdNotDeclared, match="not declared"):
            B.BridgeThresholds().require(name)


def test_a_thin_shared_set_refuses_and_forbids_lowering_the_bar():
    th = _th(min_shared_lines=20)
    s = B.line_sharing_gate(_lines(6), _lines(6), th)
    assert not s.ok
    assert "TOO LARGE" in s.reason and "do not lower the minimum" in s.reason


def test_saturated_lines_are_excluded_from_the_shared_set():
    th = _th()
    sat = _lines(10, rew=-4.0)          # outside the declared linear window
    assert B.line_sharing_gate(sat, sat, th).n_shared == 0


def test_matching_requires_excitation_potential_not_just_wavelength():
    """RYA-780: one 0.06 Å coincidence manufactured a 2.85 dex discrepancy."""
    th = _th(min_shared_lines=1)
    tgt = [dict(wavelength_air_A=5000.0, ep_eV=1.0, rew=-5.2, abundance=7.5)]
    ref = [dict(wavelength_air_A=5000.0, ep_eV=4.5, rew=-5.2, abundance=7.5)]
    assert B.line_sharing_gate(tgt, ref, th).n_shared == 0, (
        "two different transitions at the same wavelength were paired")


def test_an_ambiguous_match_is_dropped_not_guessed():
    th = _th(min_shared_lines=1)
    tgt = [dict(wavelength_air_A=5000.0, ep_eV=3.0, rew=-5.2, abundance=7.5)]
    ref = [dict(wavelength_air_A=5000.001, ep_eV=3.0, rew=-5.2, abundance=7.5),
           dict(wavelength_air_A=5000.002, ep_eV=3.0, rew=-5.2, abundance=7.4)]
    assert B.line_sharing_gate(tgt, ref, th).n_shared == 0


# ── R1–R4: the firewall ──────────────────────────────────────────────────────────────
def test_bridge_hop_is_never_given_a_reference_abundance():
    """R1 by signature: there is no parameter through which a reference could be tuned."""
    import inspect
    params = set(inspect.signature(B.bridge_hop).parameters)
    assert not (params & {"reference_abundance", "solar_absolute_dex", "published_dex"})


def test_a_hop_returns_a_difference_and_the_absolute_is_assembled_once():
    """R2."""
    th = _th()
    h = B.bridge_hop(_lines(20, a0=7.60), _lines(20, a0=7.50), th,
                     target="x", reference="solar")
    assert h.delta_dex == pytest.approx(0.10, abs=1e-9)
    out = B.chain_to_sun([h], star="x", solar_absolute_dex=7.50, zero_point_dex=0.059)
    assert out.absolute_dex == pytest.approx(7.60, abs=1e-9)
    assert out.differential_dex == pytest.approx(0.10, abs=1e-9)


def test_the_differential_is_invariant_to_the_reference_scale():
    """🔴 The cancellation, asserted: shifting BOTH stars leaves the difference untouched."""
    th = _th()
    base = B.bridge_hop(_lines(20, a0=7.60), _lines(20, a0=7.50), th, target="x", reference="s")
    for d in (-0.5, 0.3, 1.0):
        moved = B.bridge_hop(_lines(20, a0=7.60 + d), _lines(20, a0=7.50 + d), th,
                             target="x", reference="s")
        assert moved.delta_dex == pytest.approx(base.delta_dex, abs=1e-12)


def test_the_benchmark_comparison_is_a_report_not_a_term():
    """R4."""
    th = _th()
    h = B.bridge_hop(_lines(20, a0=7.60), _lines(20, a0=7.50), th, target="x", reference="solar")
    out = B.chain_to_sun([h], star="x", solar_absolute_dex=7.50, zero_point_dex=0.059)
    before = (out.differential_dex, out.absolute_dex)
    r = B.compare_to_benchmark(out, 7.20, source="GBS Jofre+2015")
    assert (out.differential_dex, out.absolute_dex) == before, "the comparison changed the product"
    assert "never fed back" in r["status"]
    assert r["difference_dex"] == pytest.approx(0.40, abs=1e-9)


def test_chain_refuses_rather_than_returning_zero_when_there_are_no_hops():
    with pytest.raises(ValueError, match="refusal rather than a zero"):
        B.chain_to_sun([], star="x")


# ── the two products ─────────────────────────────────────────────────────────────────
def test_hops_accumulate_in_quadrature():
    th = _th()
    h1 = B.bridge_hop(_lines(20, a0=7.60), _lines(20, a0=7.50), th, target="a", reference="b")
    h2 = B.bridge_hop(_lines(20, a0=7.50), _lines(20, a0=7.45), th, target="b", reference="solar")
    out = B.chain_to_sun([h1, h2], star="a", solar_absolute_dex=7.45, zero_point_dex=0.059)
    assert out.differential_dex == pytest.approx(h1.delta_dex + h2.delta_dex, abs=1e-12)
    assert out.differential_sigma_dex == pytest.approx(
        np.hypot(h1.sigma_dex, h2.sigma_dex), abs=1e-12)


def test_the_zero_point_survives_the_chain_untouched():
    """🔴 Precision from the differential, scale from the absolute — and the scale never shrinks."""
    th = _th()
    h = B.bridge_hop(_lines(200, a0=7.60), _lines(200, a0=7.50), th, target="x", reference="solar")
    out = B.chain_to_sun([h], star="x", solar_absolute_dex=7.50, zero_point_dex=0.059)
    assert out.zero_point_dex == 0.059
    assert out.differential_sigma_dex < out.zero_point_dex, (
        "this test is meant to exercise the case where the differential is far tighter than "
        "the scale it rests on")
    assert "does not shrink" in out.report()


def test_both_products_are_reported():
    th = _th()
    h = B.bridge_hop(_lines(20, a0=7.60), _lines(20, a0=7.50), th, target="x", reference="solar")
    txt = B.chain_to_sun([h], star="x", solar_absolute_dex=7.50, zero_point_dex=0.059).report()
    assert "differential" in txt and "absolute" in txt and "gf cancelled" in txt
