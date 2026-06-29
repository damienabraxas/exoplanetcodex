"""
tests/test_gf_adjudication_batch3_rya354.py
===========================================
RYA-354 Batch 3 (Zr/Ba/Y/Eu — neutron-capture / s-process) — guard the per-line gf
adjudication on the committed canonical_gf table. The Batch-3 finding: unlike Cu/V
(0 measured lines), these 4 elements ARE measured (12 solar lines) but gf grade is NOT
the binding blocker — saturation (7/12 lines, gf-moot) + in_regions=False (the HFS
s-process lines never reach the EW abundance path) are. So the batch ADJUDICATES gf
QUALITY without an unlock: stamp the best-available graded values already in place
(Eu II = LWHS/Lawler+2001, Zr II 5372 = LNAJ/Ljung+2006), confirm Ba II is NIST-A, flag
the ungraded Kurucz Y I lines — value/grade FROZEN, nothing invented, Fe anchor untouched.

Runs on the committed canonical_gf + the batch-3 audit CSV.
"""
import sys
import warnings
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
warnings.filterwarnings('ignore')

CANON = ROOT / 'data' / 'linelists' / 'canonical_gf.csv'
AUDIT = ROOT / 'data' / 'audit' / 'gf_adjudication' / 'batch3_zrbayeu_rya354.csv'

# The lines stamped to their already-present authoritative graded source (value frozen).
STAMPED = {
    ('Eu II', 6645.0905): ('adjudicated_lwhs2001_rya354', 'LWHS', 0.4208),
    ('Zr II', 5372.466):  ('adjudicated_lnaj2006_rya354', 'LNAJ', -0.81),
}
# The ungraded Kurucz Y I lines flagged for the verbatim manual pull (value frozen).
Y_FLAGGED = {4348.786: 0.725, 6191.718: -0.68, 6435.004: -0.521}


def _canon():
    return pd.read_csv(CANON, low_memory=False)


def _row(df, species, wave, tol=0.05):
    sub = df[df['species'] == species]
    i = (sub['wavelength_air_A'] - wave).abs().idxmin()
    r = df.loc[i]
    assert abs(float(r['wavelength_air_A']) - wave) <= tol
    return r


def test_eu_and_zr2_stamped_to_their_graded_source():
    df = _canon()
    for (species, wave), (status, ref, gf) in STAMPED.items():
        r = _row(df, species, wave)
        assert r['adjudication_status'] == status
        assert r['loggf_reference'] == ref                    # cited source already in place
        assert abs(float(r['log_gf']) - gf) < 1e-9            # value FROZEN (not invented)


def test_y_i_lines_flagged_manual_pull_value_frozen():
    df = _canon()
    for wave, gf in Y_FLAGGED.items():
        r = _row(df, 'Y I', wave)
        assert r['adjudication_status'] == 'flagged_manual_pull_rya354'
        assert r['loggf_reference'] == 'K06'                  # ungraded Kurucz, not relabelled
        assert abs(float(r['log_gf']) - gf) < 1e-9            # frozen
        assert pd.isna(r['nist_grade']) or str(r['nist_grade']) in ('nan', '')


def test_batch3_canonical_edit_is_surgical():
    # exactly the 2 batch-3 stamps + 3 Y I flags — nowhere else.
    df = _canon()
    assert (df['adjudication_status'] == 'adjudicated_lwhs2001_rya354').sum() == 1
    assert (df['adjudication_status'] == 'adjudicated_lnaj2006_rya354').sum() == 1
    yflag = df[(df['species'] == 'Y I') &
               (df['adjudication_status'] == 'flagged_manual_pull_rya354')]
    assert len(yflag) == 3
    assert sorted(round(float(w), 3) for w in yflag['wavelength_air_A']) == sorted(Y_FLAGGED)


def test_ba_ii_already_nist_grade_A_unchanged():
    # Ba II 5853 is already the best-available (NIST ASD grade A) — confirmed, not re-touched.
    r = _row(_canon(), 'Ba II', 5853.668)
    assert str(r['nist_grade']) == 'A'
    assert r['adjudication_status'] == 'nist'                 # untouched by batch 3


def test_gf_is_not_the_blocker_in_regions_false():
    # The honest core finding: gf adjudication cannot unlock these because the EW abundance
    # path never matches them — in_regions=False on every stamped/flagged line (HFS s-process,
    # the RYA-456 "regions file 0 for Eu/Ba" structure). Grade is moot until that is wired.
    df = _canon()
    for species, wave in list(STAMPED) + [('Y I', w) for w in Y_FLAGGED]:
        r = _row(df, species, wave)
        assert bool(r['in_regions']) is False


def test_audit_saturation_and_disposition_split():
    a = pd.read_csv(AUDIT)
    assert len(a) == 12
    # 7 of 12 lines are saturated (gf-moot — the Na/Mg Batch-1 precedent)
    assert (a['saturation'] == 'SATURATED').sum() == 7
    # 2 stamped graded + 3 flagged manual-pull
    assert a['new_adjudication_status'].str.startswith('adjudicated').sum() == 2
    assert a['new_adjudication_status'].str.startswith('flagged').sum() == 3
    # the cleanest line in the set is the Eu II HFS recovery line (RYA-102/458)
    eu = a[a['species'] == 'Eu II'].iloc[0]
    assert eu['saturation'] == 'clean' and 'STAMP graded' in eu['disposition']


def test_flagged_y_carries_named_source_not_a_catalog():
    a = pd.read_csv(AUDIT)
    yi = a[(a['species'] == 'Y I') & a['new_adjudication_status'].str.startswith('flagged')]
    assert len(yi) == 3
    for src in yi['cited_source']:
        assert 'Hannaford' in src and '1982' in src           # named lab reference
        assert 'NIST' in src
        assert 'VALD' not in src                              # never a catalog-as-source


def test_no_fe_or_off_species_touched():
    a = pd.read_csv(AUDIT)
    assert set(a['species']) <= {'Zr I', 'Zr II', 'Y I', 'Ba II', 'Ba I', 'Eu II'}
    assert not any('Fe' in s for s in a['species'])
