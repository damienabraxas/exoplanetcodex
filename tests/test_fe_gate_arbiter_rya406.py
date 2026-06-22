"""
tests/test_fe_gate_arbiter_rya406.py
====================================
RYA-406 — the solar Fe gate scores its ionization verdict on the ratified
SYNTHESIS arbiter (DECISION 2), not the blend-limited EW-path Fe II (high by
design, RYA-352). The EW ΔFe and the EW-vs-synth Fe II spread are demoted to
flagged diagnostics. Per-star: only a star with a cited arbiter constant uses it;
otherwise the EW path is scored but flagged loud (no silent verdict swap).
"""
import numpy as np
import pandas as pd
import pytest

from scripts.validate_fe_rya238 import _check_gates, GATES
from config.constants import FE_IONIZATION_SYNTH_ARBITER, FE_EW_SYNTH_SPREAD_BAND


def _abund(a_fe1, a_fe2, n_fe2=10):
    """Minimal abundances frame: Fe I/II with NLTE applied (a_nlte = A_X_nlte)."""
    return pd.DataFrame([
        {'element': 'Fe', 'ion': 'I',  'A_X': a_fe1, 'A_X_nlte': a_fe1,
         'nlte_flag': 'NLTE', 'A_X_std_nlte': 0.09, 'A_X_std': 0.09,
         'n_lines': 60, 'delta_nlte_mean': 0.009},
        {'element': 'Fe', 'ion': 'II', 'A_X': a_fe2, 'A_X_nlte': a_fe2,
         'nlte_flag': 'NLTE', 'A_X_std_nlte': 0.13, 'A_X_std': 0.13,
         'n_lines': n_fe2, 'delta_nlte_mean': 0.0},
    ])


def test_arbiter_constant_is_wired():
    """The Sun has a cited single-source arbiter (no hardcoded gate value)."""
    arb = FE_IONIZATION_SYNTH_ARBITER['solar']
    assert {'dFe', 'fe2_synth', 'provenance'} <= set(arb)
    assert abs(arb['dFe']) < GATES['solar']['dFe_max']   # the ratified Sun passes


def test_solar_ionization_scored_on_synth_arbiter():
    """ΔFe verdict = the synth arbiter (|−0.015|), NOT the EW path (|7.516−7.700|)."""
    arb = FE_IONIZATION_SYNTH_ARBITER['solar']
    res = _check_gates('solar', _abund(7.516, 7.700), pd.DataFrame(), 1.0)

    assert res['dFe'] == pytest.approx(abs(arb['dFe']))     # arbiter, not EW
    assert res['dFe'] != pytest.approx(res['dFe_ew'])
    assert res['dFe_ew'] == pytest.approx(0.184, abs=1e-3)  # EW diagnostic preserved
    assert 'synth arbiter' in res['dFe_src']
    assert res['dFe_pass'] is True                          # |−0.015| < 0.05


def test_ew_vs_synth_spread_is_a_flagged_diagnostic():
    """The EW-vs-synth Fe II spread is reported (a2 − fe2_synth), flagged if it
    exceeds the documented blend-bias band — but it is NOT the verdict."""
    arb = FE_IONIZATION_SYNTH_ARBITER['solar']
    res = _check_gates('solar', _abund(7.516, 7.700), pd.DataFrame(), 1.0)
    assert res['ew_synth_spread'] == pytest.approx(7.700 - arb['fe2_synth'])
    # +0.200 spread sits within the blend-bias band → diagnostic, not flagged
    assert (abs(res['ew_synth_spread']) > FE_EW_SYNTH_SPREAD_BAND) == res['spread_flag']
    assert res['spread_flag'] is False


def test_no_arbiter_star_falls_back_to_ew_path_flagged():
    """A star without a ratified arbiter (procyon) scores the EW ΔFe, but the source
    is flagged loud as the EW path and no synth spread is computed."""
    assert 'procyon' not in FE_IONIZATION_SYNTH_ARBITER
    res = _check_gates('procyon', _abund(7.50, 7.62, n_fe2=5), pd.DataFrame(), 1.66)
    assert res['dFe'] == pytest.approx(res['dFe_ew'])
    assert res['dFe'] == pytest.approx(abs(7.50 - 7.62), abs=1e-3)
    assert 'EW path' in res['dFe_src']
    assert not np.isfinite(res['ew_synth_spread'])
