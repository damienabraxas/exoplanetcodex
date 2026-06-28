"""
tests/test_procyon_outlier_disposition_rya281.py
================================================
RYA-281 — guard the Procyon Fe-outlier disposition + the bad-loggf verdict.

The four current RYA-458 ABUND_OUTLIER lines are EW blend contamination (GES log gf
cross-match disproves an oscillator-strength cause for all of them). They must stay
excluded (blend_flag True, priority 0) in linelist_procyon.csv, and the committed GES
cross-match audit must record zero LOGGF-CONFIRMED lines.
"""
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

LL    = ROOT / 'data' / 'linelists' / 'linelist_procyon.csv'
AUDIT = ROOT / 'data' / 'audit' / 'procyon_outlier_loggf_rya281.csv'

OUTLIERS = [('Fe', 'I', 5247.050), ('Fe', 'I', 5432.948),
            ('Fe', 'II', 5234.623), ('Fe', 'II', 6456.380)]


def test_four_outliers_excluded_in_linelist():
    d = pd.read_csv(LL, low_memory=False, comment='#')
    for el, ion, wl in OUTLIERS:
        m = d[(d['element'] == el) & (d['ion'] == ion) &
              (d['wavelength_air_A'].sub(wl).abs() < 0.04)]
        assert len(m) == 1, f"{el} {ion} {wl}: expected 1 row, got {len(m)}"
        r = m.iloc[0]
        assert bool(r['blend_flag']) is True, f"{wl} not blend_flagged"
        assert int(r['priority']) == 0, f"{wl} priority not 0"
        assert 'RYA-281' in str(r['notes'])


def test_audit_records_bad_loggf_disproved():
    df = pd.read_csv(AUDIT)
    assert len(df) == 4
    # every outlier: gf cannot explain the excess
    assert (df['loggf_flag'] == 'LOGGF-CONFIRMED').sum() == 0
    # delta gf is small for all (under 0.2 dex), the abundance excess is large
    assert (df['delta_log_gf'].abs() <= 0.20).all()
    assert (df['delta_A'] >= 0.5).all()


def test_no_loggf_value_was_silently_corrected():
    # The disposition is EXCLUSION, not gf editing — linelist log_gf for the four lines
    # is unchanged from the VALD value the audit recorded (we never tuned a gf).
    d = pd.read_csv(LL, low_memory=False, comment='#')
    audit = pd.read_csv(AUDIT)
    for _, a in audit.iterrows():
        wl = float(a['wavelength_air_A'])
        el, ion = a['species'].split()
        r = d[(d['element'] == el) & (d['ion'] == ion) &
              (d['wavelength_air_A'].sub(wl).abs() < 0.04)].iloc[0]
        assert abs(float(r['log_gf']) - float(a['vald_log_gf'])) < 0.001
