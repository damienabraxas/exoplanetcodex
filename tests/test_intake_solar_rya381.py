"""
tests/test_intake_solar_rya381.py
=================================
RYA-381 — guard the non-optical solar intake gate + the assembled extension.

Fast checks only: the intake gate is exercised on the smallest IR delivery, and the
assembled linelist_solar.csv is checked structurally (columns + span + preserved
optical core) by reading only the wavelength column.
"""
import sys
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parent.parent
LL = ROOT / 'data' / 'linelists'
sys.path.insert(0, str(ROOT / 'scripts'))
sys.path.insert(0, str(LL))

import intake_solar_nonoptical_rya381 as intake  # noqa: E402

SMALL_IR = LL / 'vald_solar_ir_18000_25000_hfson_raw.txt'
SOLAR = LL / 'linelist_solar.csv'


@pytest.mark.skipif(not SMALL_IR.exists(), reason="IR delivery not present")
def test_intake_gate_smallest_ir_accepts():
    row = intake.scan(SMALL_IR, 18000.0, 25000.0)
    row['hfs'] = intake.hfs_verdict(row)
    row['threshold'] = intake.threshold_flag(row)
    verdict = intake.accept_verdict(row, row['hfs'])
    # HFS must be confirmed ON from content (hfs: reference tags), not the filename
    assert row['hfs_tag_count'] > 0 and row['hfs'] == 'ON'
    # solar composition
    assert row['applied_MH'] is not None and abs(row['applied_MH']) <= intake.METAL_TOL
    # complete + not truncated → ACCEPT on the hard gates
    assert not row['truncated'] and row['span_ok'] and row['all_delivered']
    assert verdict == 'ACCEPT'
    # but the 0.05 vs 0.001 threshold mismatch must be flagged (not silently mixed)
    assert 'MISMATCH' in row['threshold']


@pytest.mark.skipif(not SOLAR.exists(), reason="assembled solar list not present")
def test_assembled_solar_spans_full_range_and_preserves_core():
    cols = pd.read_csv(SOLAR, nrows=0).columns.tolist()
    assert 'vald_proximity_flag' in cols and 'blend_flag' in cols
    w = pd.to_numeric(pd.read_csv(SOLAR, usecols=['wavelength_air_A'])['wavelength_air_A'],
                      errors='coerce').dropna()
    # extended below the optical (UV) and well into the IR (K band)
    assert w.min() < 2000.0
    assert w.max() > 24000.0
    # curated optical core (3780–6910 Å) still present in full
    n_optical = int(((w >= 3780.0383) & (w <= 6909.84)).sum())
    assert n_optical == 108971, f"optical core changed: {n_optical} rows in 3780–6910 Å"
    # NIR lines exist (proves the IR arm landed)
    assert int((w >= 18000.0).sum()) > 500


@pytest.mark.skipif(not SOLAR.exists(), reason="assembled solar list not present")
def test_nonoptical_orphans_tracked_to_rya379():
    # The non-optical extension lines have no home in the optical-only canonical_gf yet;
    # they must be VISIBLE but TRACKED to RYA-379 (not raw untracked failures), and the
    # optical core must carry zero orphans (no real break).
    import scripts.check_stewardship as cs
    vs = cs.check_all_stores_resolve()
    orphan_vs = [v for v in vs if v.value == 'absent from canonical_gf.csv']
    assert orphan_vs, "expected non-optical orphan violations to be reported"
    assert all(v.ticket == 'RYA-379' for v in orphan_vs), \
        "non-optical orphans must be tracked to RYA-379"
    assert all(v.tracked for v in cs.check_all_stores_resolve())
