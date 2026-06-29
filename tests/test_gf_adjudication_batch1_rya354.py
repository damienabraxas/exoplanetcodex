"""
tests/test_gf_adjudication_batch1_rya354.py
===========================================
RYA-354 Batch 1 (Na/Mg/Al) — guard the per-line gf adjudication on the committed
single-source canonical_gf table: the adopted NIST-graded values + cited references,
the uncitable flags, the validate-don't-tune invariant (gf barely moves), the solar
anchor preservation (no Fe row touched), and the Al unlock through the RYA-398 firewall.
"""
import sys
import warnings
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
warnings.filterwarnings('ignore')

CANON = ROOT / 'data' / 'linelists' / 'canonical_gf.csv'
AUDIT = ROOT / 'data' / 'audit' / 'gf_adjudication' / 'batch1_namgal_rya354.csv'


def _canon():
    return pd.read_csv(CANON, low_memory=False)


def _row(df, species, wave, tol=0.02):
    sub = df[df['species'] == species]
    i = (sub['wavelength_air_A'] - wave).abs().idxmin()
    r = df.loc[i]
    assert abs(float(r['wavelength_air_A']) - wave) <= tol
    return r


def test_al6696_adjudicated_to_nist_C_plus_with_underlying_reference():
    r = _row(_canon(), 'Al I', 6696.185)
    assert r['adjudication_status'] == 'adjudicated_rya354'
    assert r['nist_grade'] == 'C+'
    assert abs(float(r['log_gf']) - (-1.569)) < 1e-3            # log10(2 * 1.35e-2)
    # the UNDERLYING reference is named (not "VALD3"/"K75"/"GES v6")
    assert 'Kelleher' in r['loggf_reference'] and 'NIST' in r['loggf_reference']
    assert 'VALD3' not in r['loggf_reference']


def test_mg6319_triplet_adjudicated_to_nist_D_and_E():
    df = _canon()
    d = _row(df, 'Mg I', 6319.237)
    e = _row(df, 'Mg I', 6319.493)
    assert d['nist_grade'] == 'D' and e['nist_grade'] == 'E'
    assert d['adjudication_status'] == 'adjudicated_rya354'
    assert 'Kelleher' in d['loggf_reference'] and '1285' in d['loggf_reference']


def test_uncitable_lines_flagged_not_guessed():
    df = _canon()
    for species, wave in (('Mg I', 6083.522), ('Mg I', 6630.833), ('Al I', 6631.218)):
        r = _row(df, species, wave)
        assert r['adjudication_status'] == 'flagged_uncitable_rya354'
        # value NOT changed (kept the existing sourced value — never guessed)
        assert pd.isna(r['nist_grade']) or str(r['nist_grade']) == 'nan'


def test_validate_dont_tune_gf_barely_moves():
    # gf cancels in the differential — every adjudicated line moved < 0.05 dex.
    a = pd.read_csv(AUDIT)
    adj = a[a['action'] == 'ADJUDICATE']
    assert len(adj) == 3
    assert (adj['delta_loggf'].abs() < 0.05).all()


def test_solar_anchor_untouched_no_fe_row_changed():
    # The adjudication may only touch Na/Mg/Al rows — Fe (the solar anchor) is sacred.
    a = pd.read_csv(AUDIT)
    species = set(a['species'])
    assert species <= {'Na I', 'Mg I', 'Al I'}
    assert not any('Fe' in s for s in species)


def test_al_unlocks_through_the_firewall_na_mg_do_not():
    # The payoff: Al 6696 now carries a graded (MED-tier) gf and survives the
    # RYA-398 independent-gf firewall; Na/Mg stay blocked by saturation/EW-error.
    from pipeline import curate_nonfe_pools as C
    df = C.apply_cull(C.load_pool(['Na', 'Mg', 'Al']), grade_restrict=True)
    al = df[(df['element'] == 'Al')]
    al6696 = al[(al['wavelength_air_A'] - 6696.185).abs() < 0.02].iloc[0]
    assert al6696['gf_tier'] == 'MED'
    assert al6696['cull_reason'] == ''                          # KEPT -> Al unlocks
    # Na/Mg: no line survives (their blocker is SAT/HIERR, not gf — the finding)
    assert (df[df['element'] == 'Na']['cull_reason'] != '').all()
    assert (df[df['element'] == 'Mg']['cull_reason'] != '').all()


def test_counts_three_adjudicated_three_flagged_three_confirmed():
    a = pd.read_csv(AUDIT)
    assert (a['action'] == 'ADJUDICATE').sum() == 3
    assert (a['action'] == 'FLAG').sum() == 3
    assert (a['action'] == 'CONFIRM').sum() == 3       # Na 5682/5688 + Mg 5711 already NIST A
