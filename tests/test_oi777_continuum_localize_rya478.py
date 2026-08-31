"""
tests/test_oi777_continuum_localize_rya478.py
=============================================
RYA-478 Phase 2b — UVES stellar-RV-to-rest fix + continuum localizer. CI-safe: the
localizer test is pure-numpy; the real-spectrum RV/rest test is skipped without the
Procyon UVES anchor (external data store).
"""
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
from config.constants import SPECTRA_EXT_DIR  # noqa: E402  (RYA-1140)
_UVES_ANCHOR = (SPECTRA_EXT_DIR / 'Procyon' /
                'Procyon UVES' / 'ADP.2020-06-15T10:09:57.908.fits')

from scripts.procyon_oi777_continuum_localize_rya478 import (   # noqa: E402
    localize_continuum, OI777_TRIPLET_SPAN, OI777_CONT_WINDOW)


def test_localizer_flattens_a_sloped_continuum():
    # a sloped pseudo-continuum with an absorption triplet -> localized continuum ~1.0
    w = np.linspace(7762.0, 7784.0, 800)
    slope = 1.0 + 0.002 * (w - 7773.0)                 # ~0.4% slope across the window
    flux = slope.copy()
    for c in (7771.94, 7774.17, 7775.39):              # carve absorption cores
        flux -= 0.5 * np.exp(-0.5 * ((w - c) / 0.08) ** 2)
    wl, floc = localize_continuum(w, flux)
    sh = ((wl >= 7763) & (wl <= 7769)) | ((wl >= 7778) & (wl <= 7783))
    assert abs(float(np.nanmedian(floc[sh])) - 1.0) < 0.01   # shoulders flattened to ~1.0
    core = np.argmin(np.abs(wl - 7771.94))
    assert floc[core] < 0.8                            # the absorption core is preserved


def test_localizer_excludes_triplet_span():
    assert OI777_TRIPLET_SPAN[0] < 7771.94 and OI777_TRIPLET_SPAN[1] > 7775.39
    assert OI777_CONT_WINDOW[0] < OI777_TRIPLET_SPAN[0]
    assert OI777_CONT_WINDOW[1] > OI777_TRIPLET_SPAN[1]


@pytest.mark.skipif(not _UVES_ANCHOR.exists(),
                    reason='Procyon UVES anchor (external data store) absent')
def test_uves_arm_now_rest_frame_corrected():
    # the RYA-478 arm fix must bring the O I 777 cores to rest (prior: ~-0.12 A blueshift)
    import pipeline.cno_synthesis as cs
    arm = cs.star_arm_registry('procyon')['uves']
    w_nm, fl = cs.resolve_arm_spectrum('procyon', arm)
    wA = w_nm * 10.0
    for c in (7771.94, 7774.17, 7775.39):
        seg = (wA >= c - 0.4) & (wA <= c + 0.4)
        core = wA[seg][np.nanargmin(fl[seg])]
        assert abs(core - c) < 0.07                    # within ~0.07 A of rest (was 0.11-0.14)
