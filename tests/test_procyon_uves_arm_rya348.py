"""
tests/test_procyon_uves_arm_rya348.py
=====================================
RYA-348 Phase 2 — Procyon UVES O I 777 arm wiring. CI-safe: the registry/guard tests
need no data; the real-spectrum resolve test is skipped when the UVES anchor frame
(external data store) is absent.
"""
from pathlib import Path

import pytest

import pipeline.cno_synthesis as cs

ROOT = Path(__file__).resolve().parents[1]
_UVES_ANCHOR = (ROOT.parent / 'data' / 'spectra' / 'exoplanetcodex-data' / 'Procyon' /
                'Procyon UVES' / 'ADP.2020-06-15T10:09:57.908.fits')


def test_uves_arm_ready_and_wired():
    uv = cs.star_arm_registry('procyon')['uves']
    assert uv.loader == 'uves_rya272'
    assert uv.ready is True
    # O I 777 primary is in the diagnostic set
    assert any(d.key == 'OI_777' and d.role == 'primary' for d in uv.diagnostics)
    ready, _ = cs.available_arms('procyon')
    assert 'uves' in ready


def test_uves_loader_is_procyon_only_never_vesta():
    uv = cs.star_arm_registry('procyon')['uves']
    forced_solar = cs.ArmWiring('uves', uv.region, uv.diagnostics, 'uves_rya272', True)
    with pytest.raises(cs.ArmNotWired):
        cs.resolve_arm_spectrum('solar', forced_solar)     # anti-silent-Vesta


def test_oi777_carries_amarsi_grid_key():
    # the cited-correction layer must give O I 777 a real Amarsi-2019 NLTE delta
    assert 'OI_777' in cs._ATOMIC_GRID_KEYS
    assert cs._ATOMIC_GRID_KEYS['OI_777'][0] == 'OI'


@pytest.mark.skipif(not _UVES_ANCHOR.exists(),
                    reason='Procyon UVES anchor frame (external data store) absent')
def test_resolve_real_procyon_uves_arm():
    import numpy as np
    uv = cs.star_arm_registry('procyon')['uves']
    w_nm, flux = cs.resolve_arm_spectrum('procyon', uv)
    assert w_nm.size > 1000
    lo, hi = float(w_nm.min()) * 10, float(w_nm.max()) * 10
    assert lo < 7771.0 < hi                              # O I 777 covered (RED760)
    assert lo < 6300.30 < hi                             # [O I] 6300 covered
    assert 0.5 < float(np.nanmedian(flux)) < 1.5         # continuum-normalized ~1.0
