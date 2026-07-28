"""
tests/test_gf_adjudication_batch2_rya354.py
===========================================
RYA-354 Batch 2 (Mn/Cu/V — iron-peak, HFS) — guard the per-line gf adjudication on the
committed single-source canonical_gf table. The Batch-2 finding is a SPLIT, not a value
pull: 3 Mn I lines are gf-grade-limited (flagged for the verbatim manual graded pull, no
value invented) while the rest of Mn is saturation/blend-limited and Cu/V are upstream-
blocked (no measured lines reach curation). Guards: the manual-pull flags, the strict
no-value-invented invariant (log_gf + nist_grade byte-untouched), the surgical canonical
edit (only those 3 adjudication_status cells changed repo-wide), the audit per-element
split + counts, and the no-unlock consequence (Mn still culls by the RYA-398 firewall
because the grade is still owed).
"""
import sys
import warnings
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
warnings.filterwarnings('ignore')

CANON = ROOT / 'data' / 'linelists' / 'canonical_gf.csv'
AUDIT = ROOT / 'data' / 'audit' / 'gf_adjudication' / 'batch2_mncuv_rya354.csv'

# The 3 grade-only-blocked Mn I unlock candidates and their byte-frozen K07 values.
MN_UNLOCK = {5495.008: -0.952, 5720.179: 0.308, 6306.365: -0.656}


def _canon():
    return pd.read_csv(CANON, low_memory=False)


def _row(df, species, wave, tol=0.02):
    sub = df[df['species'] == species]
    i = (sub['wavelength_air_A'] - wave).abs().idxmin()
    r = df.loc[i]
    assert abs(float(r['wavelength_air_A']) - wave) <= tol
    return r


def test_three_mn_lines_resolved_to_confirmed_void_rya468():
    # RYA-468 ran the manual pull: these 3 are absent from BOTH NIST ASD and Den
    # Hartog+2011 (Kurucz-only) -> the flag resolves from 'flagged_manual_pull_rya354'
    # to the confirmed lab-gf void. The value stays K07 (asserted below); only the
    # disposition advances (the owed pull is done, the answer is "void").
    df = _canon()
    for wave in MN_UNLOCK:
        r = _row(df, 'Mn I', wave)
        assert r['adjudication_status'] == 'void_lab_gf_confirmed_rya468'


def test_no_value_invented_log_gf_and_grade_byte_untouched():
    # Batch 2 invents NOTHING: the flagged lines keep their K07 value and stay ungraded
    # (the graded value is OWED via the manual pull, not fabricated in-session).
    df = _canon()
    for wave, gf in MN_UNLOCK.items():
        r = _row(df, 'Mn I', wave)
        assert abs(float(r['log_gf']) - gf) < 1e-9            # K07 value frozen
        assert pd.isna(r['nist_grade']) or str(r['nist_grade']) == 'nan'
        assert r['loggf_reference'] == 'K07'                  # source not relabelled


def test_canonical_edit_is_surgical_only_those_three_rows():
    # RYA-468 resolved the Batch-2 manual pull: the 3 Mn I candidates are now CONFIRMED
    # lab-gf void (void_lab_gf_confirmed_rya468), so NO Mn I line carries the flag and the
    # 3 carry the confirmed-void disposition. The 'flagged_manual_pull_rya354' marker is
    # SHARED across the RYA-354 track (Batch 3 uses it for Y I), so the "none flagged"
    # invariant is scoped to Mn I.
    df = _canon()
    assert df[(df['adjudication_status'] == 'flagged_manual_pull_rya354') &
              (df['species'] == 'Mn I')].empty
    void = df[df['adjudication_status'] == 'void_lab_gf_confirmed_rya468']
    assert len(void) == 3
    assert set(void['species']) == {'Mn I'}
    assert sorted(round(float(w), 3) for w in void['wavelength_air_A']) == sorted(MN_UNLOCK)


def test_audit_split_mn_grade_limited_vs_cu_v_upstream():
    a = pd.read_csv(AUDIT)
    # Mn: 3 grade-only unlock candidates + 6 saturation/blend (not gf-fixable)
    mn = a[a['species'] == 'Mn I']
    assert (mn['action'] == 'FLAG_MANUAL_PULL').sum() == 3
    assert (mn['blocker'] == 'GRADE_ONLY -> unlock candidate').sum() == 3
    assert (mn['blocker'].str.contains('not gf-fixable')).sum() == 6
    # Cu / V: NOT gf-limited — upstream (no measured lines reach curation)
    cuv = a[a['species'].isin(['Cu I', 'V I', 'V II'])]
    assert len(cuv) == 3
    assert (cuv['action'] == 'NO_ACTION_not_gf_limited').all()
    assert (cuv['blocker'] == 'UPSTREAM (not gf)').all()


def test_flagged_lines_carry_a_named_authoritative_source():
    a = pd.read_csv(AUDIT)
    flagged = a[a['action'] == 'FLAG_MANUAL_PULL']
    assert len(flagged) == 3
    for src in flagged['manual_pull_source']:
        # underlying lab reference NAMED (never "VALD3"/"GES" as a source)
        assert 'Den Hartog' in src and '2011' in src and 'NIST' in src
        assert 'VALD' not in src and 'GES' not in src


def test_unlock_candidates_are_clean_pre_grade_but_grade_culled():
    # The whole point: these 3 pass the RYA-395 pre-grade cull (kept) yet are LOW gf_tier,
    # so the RYA-398 independent-gf firewall removes them — grade is the sole blocker.
    pre = pd.read_csv(ROOT / 'data' / 'curation' / 'nonfe_pools' / 'Mn_cull_rya395.csv',
                      comment='#')
    for wave in MN_UNLOCK:
        r = pre.loc[(pre['wavelength_air_A'] - wave).abs().idxmin()]
        assert bool(r['kept']) is True                        # clean pre-grade
        assert r['gf_tier'] == 'LOW'                          # but ungraded -> firewall culls
        assert pd.isna(r['cull_reason'])                      # not SAT/BLEND/HIERR


def test_no_unlock_mn_still_culls_through_the_firewall():
    # Consequence (honest): because the grade is still OWED, Mn does NOT unlock — every
    # Mn line is still culled by the RYA-398 firewall (unlike Al in Batch 1).
    from pipeline import curate_nonfe_pools as C
    df = C.apply_cull(C.load_pool(['Mn']), grade_restrict=True)
    assert (df[df['element'] == 'Mn']['cull_reason'] != '').all()


def test_no_fe_or_off_species_row_touched():
    a = pd.read_csv(AUDIT)
    assert set(a['species']) <= {'Mn I', 'Cu I', 'V I', 'V II'}
    assert not any('Fe' in s for s in a['species'])
