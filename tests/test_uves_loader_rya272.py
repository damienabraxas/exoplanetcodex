"""
tests/test_uves_loader_rya272.py
================================
RYA-272 — guard the UVES loader + registry: the load guards reject non-topocentric
/ non-UVES products (the double-BERV trap), the committed registry honors the
RYA-271 audit (45 IDPs, one O I anchor, valid telluric verdicts), and — when the
external FITS are present — the loader applies BERV with the correct sign.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from pipeline.loaders.uves_loader import (UVESLoader, UVESProductError,    # noqa: E402
                                          _WAVE_MIN_A, _WAVE_MAX_A)

REGISTRY = ROOT / 'data' / 'spectra' / 'procyon' / 'uves_registry.csv'
UVES_DIR = ROOT.parent / 'data' / 'spectra' / 'exoplanetcodex-data' / 'Procyon' / 'Procyon UVES'
QUAR = ROOT.parent / 'data' / 'spectra' / 'exoplanetcodex-data' / 'Procyon' / 'quarantine'
ANCHOR = 'ADP.2020-06-15T10:09:57.908.fits'
GES = 'ADP.2020-12-07T15:34:23.012.fits'
_HAS_DATA = (UVES_DIR / ANCHOR).exists()


# ── committed registry (always testable) ─────────────────────────────────────
def test_registry_honors_the_audit():
    df = pd.read_csv(REGISTRY)
    assert len(df) == 45                                   # the 45 accepted IDPs
    assert int(df['oi_anchor'].sum()) == 1                 # exactly one O I anchor
    anchor = df[df['oi_anchor']].iloc[0]
    assert anchor['filename'] == ANCHOR                    # 2013-10-08 RED760 co-add
    assert anchor['oi_telluric_verdict'] == 'CLEAN'
    assert anchor['snr'] > 200 and anchor['setting'] == 'RED760'


def test_registry_telluric_verdicts_are_audit_vocabulary():
    df = pd.read_csv(REGISTRY)
    verds = set(df['oi_telluric_verdict'].dropna().astype(str).str.strip()) - {''}
    assert verds <= {'CLEAN', 'CORRECTABLE', 'EXCLUDE'}
    # every verdict sits on an O I-covering spectrum, and only those
    has_verdict = df['oi_telluric_verdict'].notna() & (df['oi_telluric_verdict'].astype(str).str.strip() != '')
    assert df.loc[has_verdict, 'covers_OI_7771'].all()
    assert int(df['covers_OI_7771'].sum()) == int(has_verdict.sum()) == 5


def test_registry_records_applied_berv():
    df = pd.read_csv(REGISTRY)
    assert df['berv_applied_kms'].abs().max() > 15         # real BERV applied, not zeros
    anchor = df[df['oi_anchor']].iloc[0]
    assert 27 < anchor['berv_applied_kms'] < 30            # +28.5 from the audit


def test_guard_constants_sane():
    assert _WAVE_MIN_A < 4000 < 9000 < _WAVE_MAX_A
    assert issubclass(UVESProductError, RuntimeError)


# ── loader behaviour (needs the external FITS) ───────────────────────────────
@pytest.mark.skipif(not _HAS_DATA, reason="external UVES FITS not present")
def test_loader_loads_anchor_barycentric_with_full_meta():
    s = UVESLoader(UVES_DIR / ANCHOR).load()
    assert s.meta['specsys_original'] == 'TOPOCENT'
    assert s.meta['frame'].startswith('barycentric')
    assert 27 < s.meta['berv_kms'] < 30                    # BERV applied
    assert not s.meta['telluric_corrected']                # UVES IDP, not telluric-corrected
    assert all(k in s.meta for k in s._REQUIRED_META)      # base contract honoured
    lo, hi = s.wave_range_A
    assert _WAVE_MIN_A <= lo and hi <= _WAVE_MAX_A


@pytest.mark.skipif(not _HAS_DATA, reason="external UVES FITS not present")
def test_guard_rejects_quarantined_ges_file():
    with pytest.raises(UVESProductError, match='TOPOCENT'):
        UVESLoader(QUAR / GES).load()                      # HELIOCEN → must reject


@pytest.mark.skipif(not _HAS_DATA, reason="external UVES FITS not present")
def test_berv_sign_is_correct():
    # Load two epochs at opposite BERV extremes; after correction the Hα core must
    # land at the SAME systemic RV (the epoch-dependent BERV signature removed).
    from config.constants import PHYSICS
    c = PHYSICS['c_kms']
    HA = 6562.79

    def halpha_v(spec):
        w, f = spec.wave_A, spec.flux
        m = (w > HA - 3) & (w < HA + 3) & np.isfinite(f)
        w, f = w[m], f[m]
        a = np.clip(np.nanpercentile(f, 90) - f, 0, None)
        return c * (float(np.sum(w * a) / np.sum(a)) - HA) / HA

    neg = UVESLoader(UVES_DIR / 'ADP.2020-07-17T12:58:31.921.fits').load()   # BERV ~ −25
    pos = UVESLoader(UVES_DIR / 'ADP.2020-08-10T12:01:00.198.fits').load()   # BERV ~ +28
    assert (neg.meta['berv_kms'] < -15) and (pos.meta['berv_kms'] > 15)      # opposite extremes
    # corrected velocities agree (systemic) — would differ by ~|BERV_pos−BERV_neg|≈53 if sign wrong
    assert abs(halpha_v(neg) - halpha_v(pos)) < 8.0
