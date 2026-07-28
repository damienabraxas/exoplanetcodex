"""
tests/test_frame_object_contract_rya481.py
===========================================
RYA-481 — the standing loader convention, made executable. Guards both halves of
the contract and the fail-loud principle that ties them together:

  §A  OBJECT attribution by authoritative header — the central canonical map
      resolves aliases / HD numbers / composite headers, refuses ambiguity, and
      `assert_object` raises on a wrong-star or unresolvable OBJECT (the
      Vesta-as-solar / 55-Cnc-in-Procyon-tree / α-Cen-mislabel disease).
  §B  velocity frame — `corrections_for_specsys` refuses an unknown frame,
      `VelocityFrame.validate` catches the double-BERV / under-correction trap,
      and `verify_line_position` is the BERV sign-check (a fixed offset → raise).
"""
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from pipeline import frame_object_contract as foc   # noqa: E402


# ════════════════════════════════════════════════════════════════════════════
# §A  OBJECT attribution
# ════════════════════════════════════════════════════════════════════════════

def test_normalize_is_case_and_separator_insensitive():
    for raw in ('alf Cen A', 'ALF-CEN-A', 'alf_cen_a', 'AlfCenA'):
        assert foc.normalize_object(raw) == 'alfcena'
    assert foc.normalize_object('HD 128620') == 'hd128620'


@pytest.mark.parametrize('raw, star', [
    ('Sun', 'solar'), ('SUN', 'solar'), ('Sol', 'solar'),
    ('Procyon', 'procyon'), ('HD 61421', 'procyon'), ('alpha CMi', 'procyon'),
    ('alf Cen A', 'alpha_cen_a'), ('HD128620', 'alpha_cen_a'), ('Alpha Centauri A', 'alpha_cen_a'),
    ('alf Cen B', 'alpha_cen_b'), ('HD 128621', 'alpha_cen_b'),
    ('55 Cnc', '55cnc_a'), ('rho-Cnc', '55cnc_a'), ('srho01cnc', '55cnc_a'), ('HD75732', '55cnc_a'),
])
def test_resolve_object_known_identities(raw, star):
    res = foc.resolve_object(raw)
    assert res.resolved and res.canonical == star, (raw, res)


def test_resolve_composite_helios_solar_header():
    # NIRPS HELIOS solar frames carry OBJECT='SUN,FP,G2V' (Sun on fibre A, FP on B).
    res = foc.resolve_object('SUN,FP,G2V')
    assert res.resolved and res.canonical == 'solar', res


def test_55cnc_is_not_swept_into_the_cen_suffix_fallback():
    # '55 Cancri' contains no 'cen' token — must resolve to 55cnc_a, never α Cen.
    assert foc.resolve_object('55 Cancri').canonical == '55cnc_a'


def test_resolve_refuses_missing_ambiguous_and_unknown():
    assert not foc.resolve_object(None).resolved
    assert not foc.resolve_object('').resolved
    amb = foc.resolve_object('HD128620 HD128621')   # two aliases → both A and B
    assert not amb.resolved and 'ambiguous' in amb.reason
    unk = foc.resolve_object('Betelgeuse')
    assert not unk.resolved


def test_assert_object_returns_canonical_on_match():
    assert foc.assert_object('alpha CMi', expected='procyon') == 'procyon'
    assert foc.assert_object('HD128620', expected='alpha_cen_a') == 'alpha_cen_a'


def test_assert_object_raises_on_wrong_star():
    # the 55-Cnc-in-the-Procyon-tree incident (#4): wrong star, plausible photons.
    with pytest.raises(foc.ObjectAttributionError, match='expected'):
        foc.assert_object('55 Cnc', expected='procyon')


def test_assert_object_raises_on_unresolvable():
    with pytest.raises(foc.ObjectAttributionError):
        foc.assert_object('Star S5', expected='alpha_cen_a')   # the quarantine contaminant
    with pytest.raises(foc.ObjectAttributionError):
        foc.assert_object(None, expected='solar')


def test_assert_object_honours_a_custom_exception_type():
    class MyLoaderError(RuntimeError):
        pass
    with pytest.raises(MyLoaderError):
        foc.assert_object('55 Cnc', expected='procyon', exc=MyLoaderError)


def test_expected_star_accepts_canonical_id_or_alias():
    assert foc.canonical_star_id('alpha_cen_a') == 'alpha_cen_a'
    assert foc.canonical_star_id('HD 128620') == 'alpha_cen_a'
    with pytest.raises(foc.ObjectAttributionError):
        foc.canonical_star_id('NoSuchStar')


def test_canonical_map_is_one_to_one():
    # the index refuses to map one alias to two stars (the silent-ambiguity guard)
    seen = {}
    for star, aliases in foc._STAR_ALIASES.items():
        for a in aliases:
            assert foc.normalize_object(a) == a, f"alias {a!r} is not pre-normalized"
            assert a not in seen, f"alias {a!r} shared by {seen.get(a)} and {star}"
            seen[a] = star


def test_register_aliases_extends_and_refuses_collision():
    foc.register_aliases('procyon', ['Alpha Canis Minoris'])
    try:
        assert foc.resolve_object('Alpha Canis Minoris').canonical == 'procyon'
        with pytest.raises(foc.ObjectAttributionError):
            foc.register_aliases('solar', ['procyon'])   # alias already owned by procyon
    finally:
        foc._STAR_ALIASES['procyon'].discard('alphacanisminoris')
        foc._index_aliases()


# ════════════════════════════════════════════════════════════════════════════
# §B  velocity frame
# ════════════════════════════════════════════════════════════════════════════

def test_corrections_for_specsys_table():
    assert foc.corrections_for_specsys('TOPOCENT').berv_needed is True
    assert foc.corrections_for_specsys('BARYCENT').berv_needed is False
    assert foc.corrections_for_specsys('heliocen').berv_needed is False   # case-insensitive


def test_corrections_for_specsys_refuses_unknown_frame():
    with pytest.raises(foc.VelocityFrameError, match='not a recognized'):
        foc.corrections_for_specsys('WOBBLY')


def test_velocity_frame_validate_catches_double_berv():
    # BARYCENT + BERV applied = the HARPS/α-Cen double-correction trap (RYA-479)
    with pytest.raises(foc.VelocityFrameError, match='double-correction'):
        foc.VelocityFrame('BARYCENT', berv_applied=True, systemic_rv_applied=False,
                          wave_units='air').validate()


def test_velocity_frame_validate_catches_under_correction():
    with pytest.raises(foc.VelocityFrameError, match='under-correction'):
        foc.VelocityFrame('TOPOCENT', berv_applied=False, systemic_rv_applied=False,
                          wave_units='air').validate()


def test_velocity_frame_valid_cases_pass_and_declare():
    uves = foc.VelocityFrame('TOPOCENT', berv_applied=True, systemic_rv_applied=False,
                             wave_units='air').validate()
    harps = foc.VelocityFrame('BARYCENT', berv_applied=False, systemic_rv_applied=True,
                              wave_units='air').validate()
    assert 'BERV applied' in uves.declare() and 'TOPOCENT' in uves.declare()
    assert 'stellar-rest' in harps.declare()


def test_velocity_frame_units_must_be_declared():
    with pytest.raises(foc.VelocityFrameError):
        foc.VelocityFrame('TOPOCENT', berv_applied=True, systemic_rv_applied=False,
                          wave_units='neither')


# ── the BERV sign-check (synthetic Gaussian absorption line) ──────────────────
def _synthetic_line(center_A, n=400, span=6.0, depth=0.6, width=0.15):
    w = np.linspace(center_A - span / 2, center_A + span / 2, n)
    f = 1.0 - depth * np.exp(-0.5 * ((w - center_A) / width) ** 2)
    return w, f


def test_verify_line_position_passes_at_rest():
    HA = 6562.79
    w, f = _synthetic_line(HA)
    v = foc.verify_line_position(w, f, HA, tol_kms=2.0)
    assert abs(v) < 1.0


def test_verify_line_position_flags_a_fixed_offset():
    # a frame left un-corrected by ~30 km/s lands the line off rest → must raise
    HA = 6562.79
    shift = 1.0 + 30.0 / foc._C_KMS
    w, f = _synthetic_line(HA)
    with pytest.raises(foc.VelocityFrameError, match='km/s'):
        foc.verify_line_position(w * shift, f, HA, tol_kms=8.0)


def test_line_velocity_recovers_the_shift_sign():
    HA = 6562.79
    w, f = _synthetic_line(HA)
    v_pos, _ = foc.line_velocity_kms(w * (1.0 + 20.0 / foc._C_KMS), f, HA)
    v_neg, _ = foc.line_velocity_kms(w * (1.0 - 20.0 / foc._C_KMS), f, HA)
    assert v_pos > 15 and v_neg < -15            # sign preserved
    assert abs((v_pos - v_neg) - 40.0) < 2.0     # magnitude ≈ 40 km/s apart


def test_verify_line_position_refuses_empty_or_flat_window():
    with pytest.raises(foc.VelocityFrameError):
        foc.verify_line_position(np.array([4000., 4001., 4002.]), np.array([1., 1., 1.]), 6562.79)
    flat_w = np.linspace(6560, 6566, 100)
    with pytest.raises(foc.VelocityFrameError, match='contrast'):
        foc.verify_line_position(flat_w, np.ones_like(flat_w), 6562.79)


# ════════════════════════════════════════════════════════════════════════════
# loader integration — the contract is wired, not just available
# ════════════════════════════════════════════════════════════════════════════

def test_hst_loader_target_guard_routes_through_the_central_map():
    from pipeline.loaders import hst_uv_loader as H
    H._assert_procyon('HD61421', 'f.fits')             # accept Procyon
    H._assert_procyon('alpha CMi', 'f.fits')
    for bad in ('55 Cnc', 'rho-Cnc', 'SOME OTHER STAR'):
        with pytest.raises(H.HstUvProductError):
            H._assert_procyon(bad, 'f.fits')           # reject contaminants / strays


def test_loaders_expose_an_expected_star_hook():
    # both single-file loaders accept the §A verification hook (default None = no-op)
    import inspect
    from pipeline.loaders.uves_loader import UVESLoader
    from pipeline.loaders.spirou_loader import SPIRouLoader
    assert 'expected_star' in inspect.signature(UVESLoader.__init__).parameters
    assert 'expected_star' in inspect.signature(SPIRouLoader.__init__).parameters
