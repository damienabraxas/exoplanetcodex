"""
tests/test_mn_denhartog_rya468.py
=================================
RYA-468 — guard the Mn I gf disposition on canonical_gf:
  * 5495.008 / 5720.179 / 6306.365 recorded CONFIRMED VOID (Kurucz-only; absent from
    both NIST ASD and Den Hartog+2011) — value untouched, status flipped.
  * 6013.51 / 6016.67 / 6021.82 adjudicated to the Den Hartog+2011 (ApJS 194,35) values
    — value barely moves (validate-don't-tune), ref -> the graded lab source, so the gf
    clears the RYA-398 firewall (MED tier).
The FINDING the Linear comment records: these graded lines are SAT-culled at solar
(EW ~95 mA, REW -4.78) so the EW path does NOT unlock Mn — the HFS-synthesis route
(RYA-411) is the correct one. The gf grade is no longer Mn's blocker; saturation is.
"""
import sys
import warnings
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
warnings.filterwarnings('ignore')
from pipeline.curate_nonfe_pools import _gf_tier  # noqa: E402

CANON = ROOT / 'data' / 'linelists' / 'canonical_gf.csv'
AUDIT = ROOT / 'data' / 'audit' / 'gf_adjudication' / 'mn_denhartog_rya468.csv'


def _row(wave, tol=0.02):
    g = pd.read_csv(CANON, low_memory=False)
    mn = g[g['species'] == 'Mn I']
    r = mn.iloc[(mn['wavelength_air_A'] - wave).abs().argmin()]
    assert abs(float(r['wavelength_air_A']) - wave) <= tol
    return r


def test_den_hartog_triplet_adjudicated_and_clears_gf_firewall():
    for wave, gf in ((6013.51, -0.354), (6016.67, -0.181), (6021.82, -0.054)):
        r = _row(wave)
        assert r['adjudication_status'] == 'adjudicated_rya468'
        assert abs(float(r['log_gf']) - gf) < 1e-3
        assert 'Den Hartog' in r['loggf_reference'] and 'ApJS 194' in r['loggf_reference']
        assert 'VALD3' not in r['loggf_reference'] and 'K0' not in r['loggf_reference']
        # the gf grade is no longer the blocker — these now clear the firewall (MED)
        assert _gf_tier(r['nist_grade'], r['loggf_reference']) == 'MED'


def test_three_originals_recorded_confirmed_void_not_guessed():
    for wave in (5495.008, 5720.179, 6306.365):
        r = _row(wave)
        assert r['adjudication_status'] == 'void_lab_gf_confirmed_rya468'
        assert r['loggf_reference'] == 'K07'         # value untouched — Kurucz, never guessed
        assert _gf_tier(r['nist_grade'], r['loggf_reference']) == 'LOW'


def test_validate_dont_tune_gf_barely_moves():
    a = pd.read_csv(AUDIT)
    adj = a[a['action'] == 'ADJUDICATE']
    assert len(adj) == 3 and (adj['delta_loggf'].abs() < 0.05).all()


def test_only_mn_rows_touched_fe_anchor_safe():
    a = pd.read_csv(AUDIT)
    assert set(a['species']) == {'Mn I'}
    assert (a['action'] == 'VOID').sum() == 3 and (a['action'] == 'ADJUDICATE').sum() == 3
