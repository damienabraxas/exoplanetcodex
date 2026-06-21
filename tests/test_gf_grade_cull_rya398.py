"""
tests/test_gf_grade_cull_rya398.py
==================================
RYA-398 — the independent-gf (grade) cull + the validation↔survey firewall.

Pins the contract: the GRADE cull is abundance-blind and only fires under
--grade-restrict; accepted tiers are HIGH/MED; the solar-validation path raises
on RYA-161 astrophysical-gf columns; the graded verdict reads against the fixed
Asplund band (VALIDATED / RESIDUAL / LOW_CONFIDENCE) and a RESIDUAL is reported,
never tuned. No iSpec needed.
"""
import numpy as np
import pandas as pd
import pytest

from pipeline import curate_nonfe_pools as cur


# ── grade cull is opt-in, tier-based, abundance-blind ────────────────────────

def test_grade_reason_only_under_restrict():
    low = {'ew_mA': 40.0, 'rew': -5.5, 'err_frac': 0.1, 'gf_tier': 'LOW'}
    assert cur.cull_reasons(low) == []                       # 395 mode keeps it (flagged only)
    assert cur.cull_reasons(low, grade_restrict=True) == ['GRADE']


def test_grade_keeps_accepted_tiers():
    for tier in ('HIGH', 'MED'):
        row = {'ew_mA': 40.0, 'rew': -5.5, 'err_frac': 0.1, 'gf_tier': tier}
        assert cur.cull_reasons(row, grade_restrict=True) == []
    assert cur.ACCEPTED_GF_TIERS == frozenset({'HIGH', 'MED'})


def test_grade_stacks_with_quality_cuts():
    # a saturated Kurucz line is culled for BOTH SAT and GRADE under restrict
    row = {'ew_mA': 40.0, 'rew': -4.0, 'err_frac': 0.1, 'gf_tier': 'LOW'}
    r = cur.cull_reasons(row, grade_restrict=True)
    assert 'SAT' in r and 'GRADE' in r


def test_grade_cull_abundance_blind():
    rng = np.random.default_rng(3)
    n = 30
    df = pd.DataFrame({
        'ew_mA': rng.uniform(10, 90, n), 'wavelength_air_A': rng.uniform(4000, 7000, n),
        'ew_err_mA': rng.uniform(0.5, 3, n), 'blend_flag': False, 'nist_grade': np.nan,
        'gf_tier': rng.choice(['HIGH', 'MED', 'LOW'], n), 'A_lte': rng.uniform(5, 9, n),
    })
    df['rew'] = np.log10(df['ew_mA'] / 1000.0 / df['wavelength_air_A'])
    df['err_frac'] = df['ew_err_mA'] / df['ew_mA']
    df['blend_ratio'] = 0.0
    base = cur.apply_cull(df, grade_restrict=True)['cull_reason'].values
    df2 = df.copy(); df2['A_lte'] = rng.permutation(df2['A_lte'].values)
    shuf = cur.apply_cull(df2, grade_restrict=True)['cull_reason'].values
    assert (base == shuf).all()
    assert any('GRADE' in r for r in base)                   # the LOW lines were grade-culled


# ── validation↔survey firewall ───────────────────────────────────────────────

@pytest.mark.parametrize('col', ['log_gf_astro', 'loggf_astro', 'delta_log_gf',
                                 'delta_loggf', 'gf_astro'])
def test_firewall_raises_on_astrophysical_gf(col):
    df = pd.DataFrame({'wavelength_air_A': [5000.0], col: [0.1]})
    with pytest.raises(ValueError, match='firewall'):
        cur.assert_no_astrophysical_gf(df, 'test')


def test_firewall_passes_clean_pool():
    df = pd.DataFrame({'wavelength_air_A': [5000.0], 'log_gf': [-1.0]})
    cur.assert_no_astrophysical_gf(df, 'test')               # must not raise


# ── graded verdict + acceptance band ─────────────────────────────────────────

def test_validation_tol_band():
    assert cur.validation_tol('Cr') == pytest.approx(max(2 * 0.04, 0.10))
    assert cur.validation_tol('Mn') == pytest.approx(max(2 * 0.05, 0.10))


def _graded_pool(element, a_value, n):
    return pd.DataFrame({'element': element, 'ion': 'I', 'kept': True,
                         'A_lte': np.full(n, a_value),
                         'rew': np.linspace(-5.6, -5.0, n), 'gf_tier': 'MED'})


def test_graded_verdict_validated_residual_lowconf():
    el = 'Cr'; asp = cur.SOLAR_ASPLUND2021[el]; d = cur.solar_nlte_delta(el)
    # NLTE-applied lands on Asplund → VALIDATED (enough lines)
    on = _graded_pool(el, asp - d, 8)
    assert cur.element_diagnostic(el, on, graded=True)['verdict'] == 'VALIDATED'
    # gross offset gone but +0.4 residual survives → RESIDUAL (the Cr case)
    res = _graded_pool(el, asp - d + 0.40, 8)
    rec = cur.element_diagnostic(el, res, graded=True)
    assert rec['verdict'] == 'RESIDUAL'
    assert rec['residual_nlte_vs_asplund'] == pytest.approx(0.40, abs=1e-6)
    # too few graded lines → LOW_CONFIDENCE regardless of where it lands
    thin = _graded_pool(el, asp - d, 2)
    assert cur.element_diagnostic(el, thin, graded=True)['verdict'] == 'LOW_CONFIDENCE'


def test_graded_applies_nlte_to_curated():
    el = 'Cr'; d = cur.solar_nlte_delta(el)
    df = _graded_pool(el, 5.95, 7)
    rec = cur.element_diagnostic(el, df, graded=True)
    assert rec['A_nlte_curated'] == pytest.approx(5.95 + d, abs=1e-3)
