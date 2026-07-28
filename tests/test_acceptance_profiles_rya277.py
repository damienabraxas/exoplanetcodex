"""
tests/test_acceptance_profiles_rya277.py
========================================
RYA-277 — per-spectral-type acceptance profiles. Acceptance criteria are a FUNCTION of
stellar type, not universal. This guards the full RYA-277 build on top of the RYA-446
G-anchor seed:

  • G = solar anchor, current values RELABELLED (pure routing — solar result unchanged;
    the bit-for-bit regression is the CLI smoke test, not CI).
  • F = Procyon anchor, populated from real data: Fe I scatter floor 0.222 FINAL
    (RYA-281), nlte_available=False (Fe NLTE grid runs out > 6500 K), A(Fe)/dFe/n_Fe2/vmic
    from the RYA-273 F-anchor run.
  • K, M = physics-note STUBS, not populated → the accessor fails LOUD (no silent default).
  • The validate gate table is ROUTED from the profiles (single source), not hardcoded.
"""
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from config import constants as C  # noqa: E402


# ── Structure ────────────────────────────────────────────────────────────────
def test_g_and_f_anchors_populated():
    g = C.get_acceptance_profile('G')
    f = C.get_acceptance_profile('F')
    # Every gate criterion present on both populated anchors.
    for prof in (g, f):
        for key in C._ACCEPTANCE_GATE_KEYS:
            assert key in prof
        assert prof['method'] == 'EW'
        assert 'provenance' in prof


def test_f_anchor_is_the_rya281_final_floor():
    f = C.get_acceptance_profile('F')
    assert f['fe1_scatter_max'] == 0.222          # RYA-281 FINAL, not the provisional 0.282
    assert f['nlte_available'] is False           # Fe NLTE grid runs out > 6500 K (Procyon 6554 K)
    assert (f['a_fe_lo'], f['a_fe_hi']) == (7.38, 7.54)
    assert f['n_Fe2_min'] == 3
    assert f['vmic_lit'] == 1.66
    assert 'RYA-281' in f['provenance']


def test_g_anchor_is_a_pure_relabel_of_the_solar_gate_values():
    # The G profile must carry the EXACT current solar gate values (the regression proof
    # is that routing solar through them changes nothing).
    g = C.get_acceptance_profile('G')
    assert g['fe1_scatter_max'] == 0.1398
    assert g['nlte_available'] is True
    assert g['a_fe_lo'] == C.SOLAR_ASPLUND2021['Fe'] + C.FE_GATE_LOWER   # 7.41
    assert g['a_fe_hi'] == C.SOLAR_ASPLUND2021['Fe'] + C.FE_GATE_UPPER   # 7.51
    assert g['dFe_max'] == C.FE_IONISATION_GATE                          # 0.05
    assert g['n_Fe2_min'] == 8
    assert g['vmic_lit'] == 1.00


def test_scatter_floor_rises_g_to_f():
    # The whole point of a per-type gate: the floor RISES from cool G to warm F.
    assert C.get_acceptance_profile('G')['fe1_scatter_max'] < \
           C.get_acceptance_profile('F')['fe1_scatter_max']


# ── Stubs fail loud ──────────────────────────────────────────────────────────
def test_km_are_physics_note_stubs():
    for t in ('K', 'M'):
        stub = C.ACCEPTANCE_PROFILES[t]
        assert stub['status'] == 'TODO'
        assert stub['physics_note']            # carries the §2 physics rationale
    # M anticipates a METHOD change (synthesis), not just thresholds.
    assert C.ACCEPTANCE_PROFILES['M']['method'] == 'synthesis'


def test_get_acceptance_profile_fails_loud_on_stub_and_unknown():
    for t in ('K', 'M'):
        with pytest.raises(KeyError, match='deferred RYA-277 stub'):
            C.get_acceptance_profile(t)
    with pytest.raises(KeyError):
        C.get_acceptance_profile('Z')


# ── Routing ──────────────────────────────────────────────────────────────────
def test_acceptance_profile_for_star_routes_by_type():
    spec_g, prof_g = C.acceptance_profile_for_star('solar')
    spec_f, prof_f = C.acceptance_profile_for_star('procyon')
    assert spec_g == 'G' and prof_g['nlte_available'] is True
    assert spec_f == 'F' and prof_f['fe1_scatter_max'] == 0.222


def test_acceptance_profile_for_star_unknown_is_loud():
    with pytest.raises(KeyError):
        C.acceptance_profile_for_star('betelgeuse')   # no STAR_SPECTRAL_TYPE mapping


def test_validate_gate_table_is_profile_routed_not_hardcoded():
    # The gate table must be BUILT from ACCEPTANCE_PROFILES, with no hardcoded per-star
    # threshold literals left in the script (the F floor 0.15 → profile 0.222).
    src = (ROOT / 'scripts' / 'validate_fe_rya238.py').read_text()
    assert '_gates_from_profile' in src
    assert "'scatter_max': 0.15" not in src       # the old hardcoded F floor is gone
    assert 'get_acceptance_profile' in src
