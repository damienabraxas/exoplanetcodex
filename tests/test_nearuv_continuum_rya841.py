"""
tests/test_nearuv_continuum_rya841.py — RYA-841
===============================================
RYA-836 made the pseudo-continuum the dominant near-UV systematic. This ticket found that
the term was never measured, and that the lever it silently assumes is wrong. Two classes
of thing are pinned here:

  * THE ARITHMETIC THAT IS STILL IN THE CODE. `error_budget.py` converts "~10% uncertain
    normalisation" into "0.10 dex" one-to-one. That conversion is still there, so the test
    that describes it must keep passing until it is fixed -- and must FAIL the moment
    somebody changes the number without changing the reasoning.

  * THE NEGATIVE RESULTS. Nothing measured predicts the near-UV line-to-line scatter, and
    the one cut that looks like it tightens the pool is contradicted by the control pool.
    Negative results are exactly what gets quietly forgotten and re-discovered, and the
    dangerous version of forgetting this one is "drop the blended lines and the near-UV
    tightens" -- which the data refute.

The fits live on Sirius; what runs here is the reasoning, the reported numbers, and the
guards that would catch the tempting mistakes.
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

RES = ROOT / "data" / "results" / "rya841"
DRIVER_JSON = RES / "rya841_scatter_driver.json"
DRIVER_CSV = RES / "rya841_scatter_driver_per_line.csv"


@pytest.fixture(scope="module")
def driver():
    if not DRIVER_JSON.exists():
        pytest.skip("scatter-driver summary absent (Sirius artifact)")
    return json.loads(DRIVER_JSON.read_text())


@pytest.fixture(scope="module")
def per_line():
    if not DRIVER_CSV.exists():
        pytest.skip("scatter-driver per-line absent (Sirius artifact)")
    return pd.read_csv(DRIVER_CSV)


# ── the term as the code actually computes it ─────────────────────────────────

def test_the_pseudo_continuum_term_is_still_a_declared_constant():
    """0.10 dex is written into the budget, not derived from anything.

    This is not a style complaint. The number is added whenever a band policy mentions a
    pseudo-continuum, with no reference to how strongly THAT band's abundances respond to
    a continuum shift -- so it cannot be right for more than one band by construction.
    """
    src = (ROOT / "pipeline" / "error_budget.py").read_text()
    assert 'Term("pseudo-continuum", 0.10' in src, (
        "the pseudo-continuum term changed — re-derive RYA-841's conclusions before "
        "editing this test")


def test_the_budget_converts_percent_flux_into_dex_one_to_one():
    """The reasoning, pinned. '~10% uncertain normalisation' -> 0.10 dex asserts that
    dA/d(delta) = 1.0 dex per unit fractional continuum error. RYA-841 measures that
    derivative and it is not 1. When this comment is fixed, this test should be updated
    in the same commit as the number."""
    src = (ROOT / "pipeline" / "error_budget.py").read_text()
    assert "uncertain at the" in src and "10% level" in src
    assert "dA/" not in src, (
        "error_budget.py now mentions a derivative — if the conversion was made explicit, "
        "update this test and the 0.10 together")


# ── the selector does not do what its name suggests ───────────────────────────

def test_min_sep_is_a_sampling_rule_not_a_blend_rejection(driver):
    """`select_lines(min_sep_A=4.0)` compares each candidate only against lines ALREADY
    SELECTED, so it spreads the sample and never clears the fit window. The measured gap
    between selected lines is >4 A while the real nearest neighbour is ~0.01 A away."""
    s = driver["selector_is_a_sampling_rule_not_a_blend_rejection"]
    assert s["min_sep_between_selected_A"] >= 4.0
    assert s["median_nn_sep_full_list_A"] < 0.05, (
        "the nearest actual neighbour is no longer sub-0.05 A — the crowding premise of "
        "this band has changed")
    assert s["median_other_lines_in_window"] > 20


def test_no_near_uv_window_is_dominated_by_its_own_line(driver):
    """Window-wide, the target line supplies ~8% of the absorption in its own fit window.
    This is why 'just use the clean lines' is not available here: there are none."""
    for pool in driver["pools"].values():
        assert pool["depth_frac_median"] < 0.20
        assert pool["depth_frac_max"] < 0.40


# ── the negative results, which are the actual findings ───────────────────────

def test_nothing_measured_predicts_the_per_line_residual(driver):
    """Seven metrics, two independent pools, every one null. This extends RYA-759's
    'correlates with nothing' from raw line properties to blend-aware ones."""
    keys = [k for k in next(iter(driver["pools"].values())) if k.startswith("spearman_")]
    assert len(keys) >= 5
    for name, pool in driver["pools"].items():
        for k in keys:
            rho = pool[k]
            if rho is None or not np.isfinite(rho):
                continue
            assert abs(rho) < 0.35, (
                f"{name}/{k} = {rho:+.3f} is no longer null — a driver for the near-UV "
                f"scatter may have been found, which changes RYA-841's conclusion")


def test_the_dominance_cut_is_contradicted_by_the_control_pool(driver):
    """THE TRAP THIS TICKET EXISTS TO CLOSE. Cutting the 832 pool on core dominance takes
    its scatter 0.412 -> 0.321, which looks like a tightening worth having. The same cut
    on the lab pool makes it WORSE. A real physical criterion cannot do both, so the
    apparent gain is sampling noise -- and RYA-161 forbids taking it."""
    p832 = driver["pools"].get("rya832")
    p836 = driver["pools"].get("rya836")
    if not (p832 and p836):
        pytest.skip("both pools required")

    def at(pool, t):
        for row in pool["dominance_sweep"]:
            if abs(row["t"] - t) < 1e-9:
                return row
        return None

    a, b = at(p832, 0.50), at(p836, 0.50)
    assert a["scatter"] < p832["scatter"], "832 no longer tightens under the cut"
    assert b["scatter"] > p836["scatter"], "836 no longer worsens under the cut"
    # and neither move is significant
    assert a["n_sigma_vs_full"] < 3.0 and b["n_sigma_vs_full"] < 3.0


def test_saturation_does_not_explain_the_pool_difference(driver):
    """The lab pool sits mostly above the 832 selector's depth ceiling and its saturated
    lines come out +0.117 dex high — which looks like the explanation and is not one.
    Pinned WITH its significance so the offset is never quoted bare."""
    sat = driver["pools"]["rya836"]["saturation"]
    assert sat["n_saturated"] > sat["n_unsaturated"]
    assert sat["permutation_p"] > 0.05
    lo, hi = sat["ci95"]
    assert lo < 0 < hi, "the saturation offset became significant — re-derive"


def test_the_pools_are_effectively_disjoint(driver):
    """RYA-836's owed item 2, answered: the labs and the clean-line criterion barely
    overlap, so the 0.238 dex difference is between two nearly disjoint samples and
    cannot be decomposed line by line."""
    assert driver["lab_and_selected_overlap"]["n"] <= 3


def test_the_pools_do_not_differ_in_wavelength_distribution(driver):
    """The obvious explanation, ruled out: 759 found the scatter worsens blueward, but
    the two pools sit at the same median wavelength with the same blue fraction."""
    pools = driver["wavelength_profile"]["pools"]
    meds = [p["median_wave_A"] for p in pools.values()]
    fracs = [p["frac_below_3300"] for p in pools.values()]
    assert max(meds) - min(meds) < 100.0
    assert max(fracs) - min(fracs) < 0.10


def test_the_pools_run_opposite_ways_with_wavelength(driver):
    """Unexplained and flagged rather than smoothed over: the depth-selected pool gets
    TIGHTER redward while the lab pool gets WORSE. Whatever separates them is not a line
    property measured here."""
    bins = driver["wavelength_profile"]["bins"]
    blue = bins[0]
    red = bins[-1]
    names = [k for k in blue if k not in ("lo_A", "hi_A")]
    if len(names) < 2:
        pytest.skip("need both pools")
    trends = {}
    for n in names:
        if blue[n]["scatter"] and red[n]["scatter"]:
            trends[n] = red[n]["scatter"] - blue[n]["scatter"]
    assert len(trends) == 2
    assert min(trends.values()) < 0 < max(trends.values()), (
        "the two pools no longer trend oppositely with wavelength — re-derive")


# ── guards on the analysis itself ─────────────────────────────────────────────

def test_every_line_is_retained(per_line):
    """RYA-711 quarantine-not-cull: the analysis reports subsets, it never drops rows."""
    counts = per_line.pool.value_counts()
    assert counts.get("rya832", 0) == 40
    assert counts.get("rya836", 0) == 60


def test_scatter_is_never_quoted_without_a_sampling_error(driver):
    """Without s/sqrt(2(n-1)) every subset that happens to scatter less looks like a
    discovery. If a sweep row ever loses its error, the guard against that is gone."""
    for pool in driver["pools"].values():
        for row in pool["dominance_sweep"] + pool["depth_sweep"]:
            if row.get("scatter") is not None:
                assert row.get("scatter_err", 0) > 0
