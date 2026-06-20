"""
tests/test_crires_telluric_rya373.py
====================================
RYA-373 — Vesta CRIRES+ K-band loader / CO-overtone coverage / continuum-norm and
the permanent-IR-rule guards (telluric mandatory, topocentric-before-RV, no faked
engines). The molecfit + RYA-372 RV engines are external; here we verify the
molecfit-independent core on the real data and that every engine-gated step FAILS
LOUD when its engine is absent (never a silent skip / fake correction).
"""
from pathlib import Path

import numpy as np
import pytest

from pipeline import crires_telluric as ct

_HAS_DATA = ct.VESTA_CRIRES_DIR.exists() and any(ct.VESTA_CRIRES_DIR.glob('*.fits'))
needs_data = pytest.mark.skipif(not _HAS_DATA,
                                reason="Vesta CRIRES+ data not present (outside repo)")


@pytest.fixture(scope='module')
def k_co_frame():
    """A K frame with the CO(2-0) bandhead on-chip (K2192 or K2217)."""
    co = ct.co_overtone_frames(on_chip_only=True)
    assert co, "no on-chip CO frame found"
    return max(co, key=lambda f: f.snr)


# ── Loader + nm→Å boundary ───────────────────────────────────────────────────
class TestLoader:
    @needs_data
    def test_nm_to_angstrom_at_boundary(self, k_co_frame):
        # K band ~2.0-2.5 µm → 20000-25000 Å, NOT 2000-2500 (the nm trap)
        lo = min(np.nanmin(s.wave_A) for s in k_co_frame.segments)
        hi = max(np.nanmax(s.wave_A) for s in k_co_frame.segments)
        assert 19000 < lo < 25000 and 20000 < hi < 26000

    @needs_data
    def test_raw_state_topocent_uncalibrated(self, k_co_frame):
        assert k_co_frame.specsys.upper() == 'TOPOCENT'
        assert k_co_frame.fluxcal.upper() == 'UNCALIBRATED'
        assert k_co_frame.telluric_corrected is False   # not yet conditioned

    @needs_data
    def test_segments_split_by_order_detector(self, k_co_frame):
        keys = {(s.order, s.detector) for s in k_co_frame.segments}
        assert len(keys) == len(k_co_frame.segments)    # unique chip segments
        assert len(k_co_frame.segments) > 1


# ── CO-overtone coverage (gap-aware) ─────────────────────────────────────────
class TestCoCoverage:
    @needs_data
    def test_on_chip_subset_of_header_range(self):
        hdr = ct.co_overtone_frames(on_chip_only=False)
        onchip = ct.co_overtone_frames(on_chip_only=True)
        assert 0 < len(onchip) <= len(hdr)
        # the known finding: some header-covering frames have the bandhead in a gap
        assert len(onchip) < len(hdr)

    @needs_data
    def test_segment_at_returns_onchip_segment(self, k_co_frame):
        seg = k_co_frame.segment_at(ct.CO_2_0_BANDHEAD_NM)
        assert seg is not None
        wl_A = ct.CO_2_0_BANDHEAD_NM * ct.NM_TO_A
        assert np.nanmin(seg.wave_A) <= wl_A <= np.nanmax(seg.wave_A)

    @needs_data
    def test_gap_frame_header_covers_but_not_onchip(self):
        # at least one K frame covers the bandhead in header range but has it in a gap
        hdr = ct.co_overtone_frames(on_chip_only=False)
        gap = [f for f in hdr if not f.on_chip(ct.CO_2_0_BANDHEAD_NM)]
        assert gap, "expected ≥1 header-covers-but-gap frame (K2166/K2148)"
        for f in gap:
            assert f.covers_nm(ct.CO_2_0_BANDHEAD_NM)        # header says yes
            assert f.segment_at(ct.CO_2_0_BANDHEAD_NM) is None  # but off-chip


# ── Continuum normalization (molecfit-independent) ───────────────────────────
class TestContinuum:
    @needs_data
    def test_normalizes_co_order(self, k_co_frame):
        seg = k_co_frame.segment_at(ct.CO_2_0_BANDHEAD_NM)
        norm = ct.continuum_normalize(seg.wave_A, seg.flux)
        fin = np.isfinite(norm)
        assert fin.sum() > 100
        assert 0.9 < np.nanmedian(norm[fin]) < 1.1     # continuum ~ 1
        assert np.nanmin(norm[fin]) < 0.9              # absorption cores present
        assert np.nanpercentile(norm[fin], 99) < 1.2   # no runaway

    def test_continuum_handles_empty(self):
        out = ct.continuum_normalize(np.array([1.0, 2.0]), np.array([np.nan, np.nan]))
        assert np.all(np.isnan(out))


# ── Permanent-IR-rule guards (no faked engines, no silent fallback) ──────────
class TestGuards:
    def _fake_frame(self, telluric=False, specsys='TOPOCENT'):
        return ct.CriresFrame(path=Path('x.fits'), wlen_id='K2192', band='K',
                              mjd=0.0, date_obs='', ra=0.0, dec=0.0, snr=100.0,
                              specsys=specsys, fluxcal='UNCALIBRATED',
                              wmin_nm=2000.0, wmax_nm=2400.0,
                              telluric_corrected=telluric)

    def test_assert_telluric_corrected_raises(self):
        with pytest.raises(ct.TelluricNotCorrectedError):
            ct.assert_telluric_corrected(self._fake_frame(telluric=False))

    def test_assert_telluric_corrected_passes_when_set(self):
        ct.assert_telluric_corrected(self._fake_frame(telluric=True))  # no raise

    def test_molecfit_fails_loud_when_absent(self, monkeypatch):
        monkeypatch.setattr(ct, '_esorex_available', lambda: False)
        with pytest.raises(ct.MolecfitNotAvailableError):
            ct.run_molecfit_telluric(self._fake_frame(), Path('/tmp/x'))

    def test_telluric_refuses_non_topocentric(self, monkeypatch):
        # even with the engine 'available', a non-topocentric frame must be refused
        monkeypatch.setattr(ct, '_esorex_available', lambda: True)
        with pytest.raises(RuntimeError, match='TOPOCENT'):
            ct.run_molecfit_telluric(self._fake_frame(specsys='BARYCENT'), Path('/tmp/x'))

    def test_rv_refuses_non_telluric_corrected(self):
        with pytest.raises(ct.TelluricNotCorrectedError):
            ct.rv_condition(self._fake_frame(telluric=False))

    def test_rv_fails_loud_when_372_absent(self, monkeypatch):
        # telluric flag set → passes that guard; force the 372 module absent
        monkeypatch.setattr(ct, '_load_reflected_solar_rv',
                            lambda: (_ for _ in ()).throw(
                                ct.VelocityModuleNotAvailableError("absent")))
        f = self._fake_frame(telluric=True)
        with pytest.raises(ct.VelocityModuleNotAvailableError):
            ct.rv_condition(f)


# ── molecfit integration (real esorex run; slow + env-gated) ─────────────────
import os  # noqa: E402

_RUN_MOLECFIT = (os.environ.get('RYA373_RUN_MOLECFIT') == '1'
                 and ct._esorex_available() and _HAS_DATA)


@pytest.mark.skipif(not _RUN_MOLECFIT,
                    reason="set RYA373_RUN_MOLECFIT=1 with esorex+data to run the "
                           "real molecfit telluric pass (slow, ~3 min)")
def test_real_molecfit_telluric_pass(tmp_path, k_co_frame):
    f = ct.run_molecfit_telluric(k_co_frame, tmp_path, molecules=ct.TELLURIC_MOLECULES)
    assert f.telluric_corrected is True
    assert 0.0 < f._telluric_mtrans_max <= 1.0
    seg = f.segment_at(ct.CO_2_0_BANDHEAD_NM)
    assert np.isfinite(seg.flux).sum() > 100   # corrected flux produced
    ct.assert_telluric_corrected(f)            # guard now passes
