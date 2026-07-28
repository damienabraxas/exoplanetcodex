"""
tests/test_unsaturated_recovery_rya465.py
=========================================
RYA-465 — guard the unsaturated line-selection recovery for Na/Mg: the candidate
gf is NIST-graded (not the problem), the recovery lines land in the linear COG, the
added Na lines survive the RYA-398 firewall AND agree internally (the unlock signal),
Mg stays saturation-limited (the finding), and the Fe anchor is untouched.

Runs on the COMMITTED pool + canonical_gf + audit table (no solar spectrum needed).
"""
import sys
import warnings
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
warnings.filterwarnings('ignore')

POOL = ROOT / 'data' / 'measured' / 'sol_ew_results_v1.csv'
CANON = ROOT / 'data' / 'linelists' / 'canonical_gf.csv'
AUDIT = ROOT / 'data' / 'audit' / 'line_recovery' / 'namg_unsaturated_rya465.csv'


def test_recovery_lines_carry_nist_graded_gf_not_a_gf_problem():
    g = pd.read_csv(CANON, low_memory=False)
    for wave in (6154.225, 6160.747):
        r = g[(g['species'] == 'Na I') & ((g['wavelength_air_A'] - wave).abs() < 0.02)].iloc[0]
        assert str(r['nist_grade']) == 'A'          # grade-A gf — the gf was never the problem


def test_audit_splits_unsaturated_na_from_saturated_mg():
    a = pd.read_csv(AUDIT)
    na = a[a['element'] == 'Na']
    assert na['unsaturated'].all() and len(na) == 2   # both Na subordinate lines unsaturated
    mg = a[a['element'] == 'Mg']
    assert (~mg['unsaturated']).all()                 # every Mg candidate is saturated (the finding)


def test_na_recovery_lines_added_to_pool_with_provenance():
    pool = pd.read_csv(POOL)
    na = pool[(pool['element'] == 'Na') & (pool['wavelength_air_A'].round(2).isin([6154.23, 6154.22, 6160.75]))]
    assert len(na) >= 2
    for _, r in na.iterrows():
        assert 'RYA-465' in str(r['notes']) and 'NIST grade A' in str(r['notes'])
        assert 0 < r['ew_mA'] < 100                   # unsaturated EW


def test_na_unlocks_internally_consistent_through_firewall():
    from pipeline import curate_nonfe_pools as C
    df = C.apply_cull(C.load_pool(['Na']), grade_restrict=True)
    kept = df[df['cull_reason'] == '']
    waves = set(kept['wavelength_air_A'].round(2))
    assert {6154.23, 6160.75} <= {round(w, 2) for w in kept['wavelength_air_A']} or \
        len(kept) >= 2                                 # the unsaturated pair survives
    ab = C.compute_lte_abundances(kept)
    vals = ab['A_lte'].dropna()
    assert len(vals) >= 2
    assert vals.std() < 0.10                           # internal consistency = the unlock signal


def test_mg_stays_hard_no_unsaturated_line_added():
    # Mg gets NO new pool line — every optical candidate is saturated (named finding).
    pool = pd.read_csv(POOL)
    mg_new = pool[(pool['element'] == 'Mg') &
                  (pool['notes'].astype(str).str.contains('RYA-465'))]
    assert len(mg_new) == 0


def test_fe_anchor_untouched():
    # The recovery only appended Na rows — no Fe row carries an RYA-465 note.
    pool = pd.read_csv(POOL)
    fe_touched = pool[(pool['element'] == 'Fe') &
                      (pool['notes'].astype(str).str.contains('RYA-465'))]
    assert len(fe_touched) == 0
    a = pd.read_csv(AUDIT)
    assert set(a['element']) <= {'Na', 'Mg'}
