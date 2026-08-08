"""
tests/test_ew_integrity_rya458.py
=================================
RYA-458 — the EW-verification layer contract (the EW half of RYA-451).

Pins the DISCIPLINE, not the science numbers:
  * the layer FLAGS / disposes; it never mutates a measured EW and never pulls a
    value toward a literature EW (assert_no_ew_mutation is the cardinal guard),
  * thresholds are fixed module constants, uniform across elements,
  * the three charter cases land their dispositions: C I 5380 BAD_FIT + excluded,
    Li 6707 UPPER_LIMIT, Eu 6645 RECOVERED (RYA-102) / FITTER_INCOMPLETE,
  * the reference table is cited (provenance), and a missing reference is silent.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'scripts'))

import pipeline.ew_integrity as ei  # noqa: E402
from pipeline import data_namespace as ns  # noqa: E402  RYA-469 gold solar reference

PROC = ROOT / 'data' / 'processed'
# RYA-469: phase_c reads the FROZEN gold solar reference (committed), so the verdict
# wiring tests run in CI rather than skipping on the gitignored working file.
_GOLD_SOLAR = ns.reference_path(ns.current_version())


from tests.gold_scale_blocker import (  # noqa: E402  RYA-681/669, parameterised RYA-674
    verdict_gold_version, xfail_if_regeneration_blocked)

def _frame(rows):
    return pd.DataFrame(rows, columns=['element', 'ion', 'wavelength_air_A', 'ew_mA',
                                       'ew_err_mA', 'chi2', 'a_lte'])


# ── cardinal rule: never mutate / never pull toward the reference ────────────

def test_ew_integrity_never_mutates_a_measured_ew():
    df = _frame([
        ['Li', 'I', 6707.840, 1.77, 0.33, 0.0, np.nan],   # reference exists (2.0) — must NOT pull
        ['Fe', 'I', 5000.000, 45.0, 1.0, 0.01, 7.5],
    ])
    out = ei.flag_ew_integrity(df)
    # EWs byte-identical to the input (the layer only ADDS columns)
    assert list(out['ew_mA']) == list(df['ew_mA'])
    # Li still 1.77 even though the cited reference is 2.0 (no shrink toward literature)
    assert float(out[out['element'] == 'Li'].iloc[0]['ew_mA']) == 1.77


def test_assert_no_ew_mutation_raises_on_change():
    before = _frame([['Fe', 'I', 5000.0, 45.0, 1.0, 0.01, 7.5]])
    after = before.copy()
    after['ew_mA'] = [44.0]
    with pytest.raises(AssertionError):
        ei.assert_no_ew_mutation(before, after)


def test_thresholds_are_fixed_constants():
    for v in (ei.EW_FIT_CHI2_MAX, ei.SYNTH_FIT_CHI2R_MAX, ei.ABUND_OUTLIER_NSIGMA,
              ei.REW_LINEAR_CEILING, ei.SAT_EW_CEILING_MA, ei.LIT_DEVIATION_FRAC,
              ei.EU_RECOVER_LO_MA, ei.EU_RECOVER_HI_MA):
        assert isinstance(v, float)


# ── the three charter cases ──────────────────────────────────────────────────

def _charter_from_pool():
    pool = pd.read_csv(ROOT / 'data' / 'measured' / 'sol_ew_results_v1.csv')
    pool = pool[(pool['ew_mA'] > 0) & pool['ew_mA'].notna()].reset_index(drop=True)
    return ei.charter_summary(ei.flag_ew_integrity(pool))


def test_charter_c_5380_is_bad_fit_and_excluded():
    c = _charter_from_pool()['C_I_5380']
    assert c['present']
    assert 'BAD_FIT' in c['ew_integrity']
    assert c['ew_disposition'] == 'BAD_FIT'
    assert c['ew_excluded'] is True


def test_charter_li_6707_is_upper_limit_not_excluded():
    li = _charter_from_pool()['Li_6707']
    assert li['present']
    assert li['ew_disposition'] == 'UPPER_LIMIT'
    assert li['ew_excluded'] is False        # an UL is reported, not a quality strike


def test_charter_eu_6645_recovered_by_hfs():
    eu = _charter_from_pool()['Eu_6645']
    assert eu['present']
    assert ei.EU_RECOVER_LO_MA <= eu['ew_mA'] <= ei.EU_RECOVER_HI_MA
    assert eu['ew_disposition'] == 'RECOVERED'    # RYA-102 HFS-summing recovered it


def test_eu_below_noise_would_be_fitter_incomplete():
    # synthetic Eu at the pre-RYA-102 0.3 mA (below the HFS-summed window)
    df = _frame([['Eu', 'II', 6645.127, 0.3, 0.1, 0.0, np.nan]])
    out = ei.flag_ew_integrity(df)
    assert out.iloc[0]['ew_disposition'] == 'FITTER_INCOMPLETE'


# ── the per-line flags (Mechanism A) ─────────────────────────────────────────

def test_cog_flag_on_saturated_line():
    df = _frame([['Fe', 'I', 6000.0, 150.0, 5.0, 0.01, 7.5]])    # EW>100, above the COG knee
    out = ei.flag_ew_integrity(df)
    assert 'COG_FLAG' in out.iloc[0]['ew_integrity']
    assert bool(out.iloc[0]['ew_excluded'])


def test_bad_fit_on_high_profile_chi2_and_on_synth_chi2r():
    df = pd.DataFrame([
        {'element': 'Fe', 'ion': 'I', 'wavelength_air_A': 5001.0, 'ew_mA': 40.0, 'chi2': 11.0},
        {'element': 'Fe', 'ion': 'I', 'wavelength_air_A': 5002.0, 'ew_mA': 40.0, 'synth_chi2r': 120.0},
        {'element': 'Fe', 'ion': 'I', 'wavelength_air_A': 5003.0, 'ew_mA': 40.0, 'chi2': 1.0},
    ])
    out = ei.flag_ew_integrity(df)
    assert 'BAD_FIT' in out.iloc[0]['ew_integrity']     # profile chi2 11 > 10
    assert 'BAD_FIT' in out.iloc[1]['ew_integrity']     # synth chi2r 120 > 100
    assert 'BAD_FIT' not in out.iloc[2]['ew_integrity']  # clean


def test_abund_outlier_robust():
    # a tight Fe group with one gross A(X) outlier
    rows = [['Fe', 'I', 5000.0 + k, 40.0, 1.0, 0.5, 7.50 + 0.01 * k] for k in range(6)]
    rows.append(['Fe', 'I', 5100.0, 40.0, 1.0, 0.5, 9.20])   # outlier
    out = ei.flag_ew_integrity(_frame(rows))
    assert 'ABUND_OUTLIER' in out.iloc[-1]['ew_integrity']
    assert not any('ABUND_OUTLIER' in f for f in out.iloc[:-1]['ew_integrity'])


# ── literature cross-check (Mechanism B) — cited, silent when absent ─────────

def test_lit_deviation_fires_only_with_a_cited_reference():
    df = _frame([
        ['O', 'I', 6300.304, 20.0, 0.5, 0.0, np.nan],    # cited ref 5.4 → >50% dev
        ['O', 'I', 6300.304, 5.8, 0.5, 0.0, np.nan],     # within band → no flag
        ['Fe', 'I', 5500.000, 999.0, 0.5, 0.0, np.nan],  # no reference → silent (COG only)
    ])
    out = ei.flag_ew_integrity(df)
    assert 'LIT_DEVIATION' in out.iloc[0]['ew_integrity']
    assert 'LIT_DEVIATION' not in out.iloc[1]['ew_integrity']
    assert 'LIT_DEVIATION' not in out.iloc[2]['ew_integrity']   # missing ref is silent


def test_reference_table_is_cited():
    ref = ei.load_reference_table()
    for el, ion, wav in (('Li', 'I', 6707.840), ('Eu', 'II', 6645.127),
                         ('C', 'I', 5380.337), ('O', 'I', 6300.304)):
        r = ref[(ref['element'] == el) & (ref['ion'] == ion)]
        assert not r.empty, (el, wav)
        assert isinstance(r.iloc[0]['reference'], str) and len(r.iloc[0]['reference']) > 5


# ── verdict wiring (needs the run output; skips cleanly otherwise) ────────────

@pytest.mark.skipif(not _GOLD_SOLAR.exists(),
                    reason="needs a solar run (solar_abundances.csv is gitignored)")
def test_charter_c_still_pass_after_5380_exclusion():
    # RYA-674: regenerate against the gold version the COMMITTED artifact declares it
    # was built from — not a hardcoded CURRENT. Same named input, guard at full strength.
    gold = verdict_gold_version()
    xfail_if_regeneration_blocked(gold)
    import phase_c_verdict_rya371 as P
    ab, ew, phase_a, _gold = P._load(gold)
    rows = {r['element']: r for r in P.build_verdicts(ab, ew, phase_a)}
    c = rows['C']
    assert c['verdict'] == 'PASS'                       # C survives the exclusion
    assert 'BAD_FIT' in c['owed'] and '5380' in c['owed']
    assert c['sigma'] is not None and c['sigma'] < 0.149   # spread tightened vs the full set


@pytest.mark.skipif(not _GOLD_SOLAR.exists(),
                    reason="needs a solar run (gitignored output)")
def test_charter_li_reported_upper_limit_in_verdict():
    # RYA-674: regenerate against the gold version the COMMITTED artifact declares it
    # was built from — not a hardcoded CURRENT. Same named input, guard at full strength.
    gold = verdict_gold_version()
    xfail_if_regeneration_blocked(gold)
    import phase_c_verdict_rya371 as P
    ab, ew, phase_a, _gold = P._load(gold)
    rows = {r['element']: r for r in P.build_verdicts(ab, ew, phase_a)}
    assert 'UPPER LIMIT' in rows['Li']['owed'] or 'UPPER_LIMIT' in rows['Li']['owed']
