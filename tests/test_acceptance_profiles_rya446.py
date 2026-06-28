"""
tests/test_acceptance_profiles_rya446.py
========================================
RYA-446 — guard the G-anchor Fe I scatter acceptance profile. The threshold is the
cited RYA-407 honest floor (sigma(1D)=0.1398, GES-EW 1D-LTE solar Fe I pool); the
legacy-universal 0.10 (FE_SCATTER_GATE) is superseded for the Sun but retained for
not-yet-profiled types. Threshold-only: the solar Fe science (A(Fe), ionization, REW
slope, n_lines) is unchanged — that regression is the CLI smoke test, not CI.
"""
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from config import constants as C  # noqa: E402


def test_g_profile_is_cited_rya407_floor():
    assert C.ACCEPTANCE_PROFILES['G']['fe1_scatter_max'] == 0.1398   # RYA-407 sigma(1D)
    assert 'RYA-407' in C.ACCEPTANCE_PROFILES['G']['provenance']
    assert C.fe1_scatter_threshold('G') == 0.1398


def test_legacy_retained_but_superseded_for_sun():
    # The legacy-universal value is kept (for unprofiled types) but is NOT the G value.
    assert C.FE_SCATTER_GATE == 0.10
    assert C.fe1_scatter_threshold('G') != C.FE_SCATTER_GATE


def test_threshold_flip_is_exactly_right():
    # The whole point: solar sigma 0.138 FAILS the legacy 0.10 but PASSES the cited floor.
    sigma = 0.138
    assert sigma > C.FE_SCATTER_GATE                      # FAILed legacy
    assert sigma < C.fe1_scatter_threshold('G')           # PASSes the G floor


def test_km_profiles_fail_loud_no_silent_default():
    # RYA-277: F is now POPULATED (Procyon anchor); K/M remain physics-note stubs that
    # must fail loud, never default silently.
    assert C.ACCEPTANCE_PROFILES['F'] is not None        # populated by RYA-277
    assert C.ACCEPTANCE_PROFILES['K'].get('status') == 'TODO'   # stub w/ physics note
    assert C.ACCEPTANCE_PROFILES['M'].get('status') == 'TODO'
    for t in ('K', 'M'):
        with pytest.raises(KeyError, match='No acceptance profile populated'):
            C.fe1_scatter_threshold(t)
    with pytest.raises(KeyError):
        C.fe1_scatter_threshold('Z')                      # unknown type -> loud


def test_star_spectral_type_mapping():
    assert C.STAR_SPECTRAL_TYPE['solar'] == 'G'
    assert C.STAR_SPECTRAL_TYPE['procyon'] == 'F'         # F now populated (RYA-277 anchor)


def test_no_hardcoded_scatter_threshold_in_gate_script():
    # The gate must read the threshold from the profile, not hardcode it. RYA-277
    # generalized the RYA-446 fe1_scatter_threshold('G') call into _gates_from_profile()
    # → get_acceptance_profile(), which routes EVERY per-type threshold from the profile.
    src = (ROOT / 'scripts' / 'validate_fe_rya238.py').read_text()
    assert 'get_acceptance_profile' in src               # profile-routed, not hardcoded
    # No hardcoded numeric scatter ceiling ASSIGNMENT (comments may mention values):
    assert "'scatter_max': 0." not in src                 # routed via p['fe1_scatter_max']
    assert "scatter_max = 0." not in src
