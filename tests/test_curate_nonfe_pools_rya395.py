"""
tests/test_curate_nonfe_pools_rya395.py
=======================================
RYA-395 — the curation discipline (the cardinal rule), not the science numbers.

These tests pin the *contract*: the cull is abundance-blind, thresholds are fixed
constants, blend vetting is strength-weighted (never vald_proximity), and the
Cr/Mn overshoot contract holds. They do NOT assert any A(X) lands on Asplund —
that would be the cardinal sin (tuning). The science verdict is validated by the
`--verify` smoke run, which needs iSpec; these run without it.
"""
import numpy as np
import pandas as pd
import pytest

from pipeline import curate_nonfe_pools as cur


# ── Fixed-threshold / scope sanity ───────────────────────────────────────────

def test_thresholds_are_fixed_constants():
    for v in (cur.EW_MIN_MA, cur.REW_LINEAR_CEILING, cur.SAT_EW_CEILING_MA,
              cur.EW_ERR_FRAC_MAX, cur.BLEND_FRAC_MAX, cur.RESID_TOL_DEX):
        assert isinstance(v, float)
    assert cur.REW_LINEAR_CEILING < -4.0          # a real COG knee, not a tuned number
    assert cur.PHASE1 == ['Mg', 'Si', 'Ca', 'Ni']
    assert cur.PHASE2 == ['Ti', 'Cr', 'Na', 'Al', 'Mn']


# ── The blind cull ───────────────────────────────────────────────────────────

def test_cull_reasons_are_quality_only():
    weak = cur.cull_reasons({'ew_mA': 2.0, 'rew': -6.0, 'err_frac': 0.1})
    assert weak == ['WEAK']
    sat = cur.cull_reasons({'ew_mA': 60.0, 'rew': -4.0, 'err_frac': 0.1})
    assert 'SAT' in sat                            # REW above the knee
    hierr = cur.cull_reasons({'ew_mA': 40.0, 'rew': -5.5, 'err_frac': 0.9})
    assert hierr == ['HIERR']
    blend = cur.cull_reasons({'ew_mA': 40.0, 'rew': -5.5, 'err_frac': 0.1,
                              'blend_ratio': 1.0})
    assert blend == ['BLEND']
    badgf = cur.cull_reasons({'ew_mA': 40.0, 'rew': -5.5, 'err_frac': 0.1,
                              'gf_tier': 'CULL'})
    assert badgf == ['BADGF']
    clean = cur.cull_reasons({'ew_mA': 40.0, 'rew': -5.5, 'err_frac': 0.1,
                              'blend_ratio': 0.0, 'gf_tier': 'MED'})
    assert clean == []


def test_cull_refuses_abundance_columns():
    """The cardinal rule, structurally: the cull cannot even be shown A(X)."""
    with pytest.raises(AssertionError):
        cur.cull_reasons({'ew_mA': 40.0, 'A_lte': 7.0})
    with pytest.raises(AssertionError):
        cur.cull_reasons({'ew_mA': 40.0, 'residual': 0.0})


def test_cull_mask_invariant_under_abundance_shuffle():
    rng = np.random.default_rng(1)
    n = 40
    df = pd.DataFrame({
        'ew_mA': rng.uniform(3, 140, n),
        'wavelength_air_A': rng.uniform(4000, 7000, n),
        'ew_err_mA': rng.uniform(0.5, 5, n),
        'blend_flag': False,
        'nist_grade': np.nan,
        'gf_tier': 'MED',
        'A_lte': rng.uniform(5, 9, n),
    })
    df['rew'] = np.log10(df['ew_mA'] / 1000.0 / df['wavelength_air_A'])
    df['err_frac'] = df['ew_err_mA'] / df['ew_mA']
    df['blend_ratio'] = rng.uniform(0, 0.4, n)
    base = cur.apply_cull(df)['cull_reason'].values
    df2 = df.copy()
    df2['A_lte'] = rng.permutation(df2['A_lte'].values)
    shuf = cur.apply_cull(df2)['cull_reason'].values
    assert (base == shuf).all()


# ── gf provenance tiering ────────────────────────────────────────────────────

def test_gf_tiering():
    assert cur._gf_tier('A', 'whatever') == 'HIGH'
    assert cur._gf_tier('B', 'K10') == 'HIGH'      # NIST grade wins
    assert cur._gf_tier('D', 'VALD3') == 'CULL'
    assert cur._gf_tier(np.nan, 'K07') == 'LOW'    # Kurucz theoretical
    assert cur._gf_tier(np.nan, 'KP') == 'LOW'
    assert cur._gf_tier(np.nan, 'VALD3') == 'MED'
    assert cur._gf_tier(np.nan, np.nan) == 'LOW'


# ── Strength-weighted blend metric (the vald_proximity fix) ──────────────────

def test_blend_ratio_strength_weighted():
    """A weak distant neighbour barely contributes; a strong on-core neighbour of
    another species blows the ratio up. (vald_proximity flagged both identically.)"""
    target = pd.DataFrame([{'element': 'Si', 'ion': 'I',
                            'wavelength_air_A': 5000.0, 'log_gf': -1.0,
                            'excitation_potential_eV': 4.0}])

    class _LL:
        def __init__(self, recs):
            self._r = recs
        def __getitem__(self, k):
            return np.array([r[k] for r in self._r])

    # weak, far Fe neighbour
    weak = _LL([{'wave_A': 5000.0, 'element': 'Si 1', 'loggf': -1.0,
                 'lower_state_eV': 4.0, 'molecule': 'F'},
                {'wave_A': 5000.09, 'element': 'Fe 1', 'loggf': -5.0,
                 'lower_state_eV': 5.0, 'molecule': 'F'}])
    r_weak = cur.add_blend_ratio(target, weak)['blend_ratio'].iloc[0]
    # strong, on-core Fe neighbour
    strong = _LL([{'wave_A': 5000.0, 'element': 'Si 1', 'loggf': -1.0,
                   'lower_state_eV': 4.0, 'molecule': 'F'},
                  {'wave_A': 5000.01, 'element': 'Fe 1', 'loggf': 0.5,
                   'lower_state_eV': 1.0, 'molecule': 'F'}])
    r_strong = cur.add_blend_ratio(target, strong)['blend_ratio'].iloc[0]
    assert r_weak < 0.1
    assert r_strong > cur.BLEND_FRAC_MAX
    assert r_strong > 50 * max(r_weak, 1e-6)


def test_same_species_neighbour_not_a_blend():
    """An HFS / same-species neighbour in the core is part of the line, not a blend."""
    target = pd.DataFrame([{'element': 'Mn', 'ion': 'I',
                            'wavelength_air_A': 6000.0, 'log_gf': -1.0,
                            'excitation_potential_eV': 3.0}])

    class _LL:
        def __init__(self, recs):
            self._r = recs
        def __getitem__(self, k):
            return np.array([r[k] for r in self._r])

    ll = _LL([{'wave_A': 6000.0, 'element': 'Mn 1', 'loggf': -1.0,
               'lower_state_eV': 3.0, 'molecule': 'F'},
              {'wave_A': 6000.01, 'element': 'Mn 1', 'loggf': 0.0,
               'lower_state_eV': 3.0, 'molecule': 'F'}])      # HFS partner
    assert cur.add_blend_ratio(target, ll)['blend_ratio'].iloc[0] == 0.0


# ── NLTE solar Δ comes from the grids, signs as expected ─────────────────────

def test_solar_nlte_delta_from_grids():
    # RYA-396 canaries: Cr and Mn are POSITIVE (the overshoot trap).
    assert cur.solar_nlte_delta('Cr') > 0
    assert cur.solar_nlte_delta('Mn') > 0
    # Na is the real NLTE assist (negative).
    assert cur.solar_nlte_delta('Na') < 0
    # Ni / Al have no vendored grid.
    assert np.isnan(cur.solar_nlte_delta('Ni'))
    assert np.isnan(cur.solar_nlte_delta('Al'))


# ── Verdict logic (no Asplund in the cut path) ───────────────────────────────

def test_verdict_classification():
    # tight & on-target → CLEAN ; persistent offset → escalate ; thin → low-conf
    el = 'Cr'
    asp = cur.SOLAR_ASPLUND2021[el]
    d = cur.solar_nlte_delta(el)
    on = pd.DataFrame({'element': el, 'ion': 'I', 'kept': True,
                       'A_lte': np.full(10, asp - d), 'rew': np.linspace(-5.5, -5.0, 10),
                       'gf_tier': 'HIGH'})
    assert cur.element_diagnostic(el, on)['verdict'] == 'CLEAN'
    hi = on.copy(); hi['A_lte'] = asp - d + 1.0
    assert cur.element_diagnostic(el, hi)['verdict'] == 'GF_SCALE_RYA161'
    thin = on.iloc[:2].copy()
    assert cur.element_diagnostic(el, thin)['verdict'] == 'LOW_CONFIDENCE'


def test_overshoot_contract_target():
    """Cr/Mn LTE target = Asplund − Δ_NLTE so the positive NLTE lifts onto scale."""
    for el in ('Cr', 'Mn'):
        d = cur.solar_nlte_delta(el)
        df = pd.DataFrame({'element': el, 'ion': 'I', 'kept': True,
                           'A_lte': np.full(8, cur.SOLAR_ASPLUND2021[el] - d),
                           'rew': np.linspace(-5.6, -5.0, 8), 'gf_tier': 'MED'})
        rec = cur.element_diagnostic(el, df)
        assert abs((rec['target_lte'] + rec['nlte_delta_solar']) - rec['asplund']) < 1e-6
