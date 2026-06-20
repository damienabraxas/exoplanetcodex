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


# ── RV anchors from source + air→vac boundary (no molecfit) ──────────────────
class TestRvAnchors:
    def test_air_to_vac_offset(self):
        # ~+6.3 Å at 2.3 µm; positive (vac > air); within physical range
        vac = float(ct._air_to_vac([22845.94])[0])
        assert 6.0 < vac - 22845.94 < 6.6

    @needs_data
    def test_anchors_resolved_from_source_not_hardcoded(self):
        a = ct._solar_rv_anchors(22835.0, 23012.0)
        assert a['source'].endswith('linelist_solar.csv')      # canonical source
        assert len(a['air']) >= 2                              # selection found lines
        # every returned wavelength is actually present in the source file (not embedded)
        import pandas as pd
        d = pd.read_csv(a['source'], comment='#', low_memory=False)
        wl = d['wavelength_air_A'].astype(float).to_numpy()
        for w in a['air']:
            assert np.min(np.abs(wl - w)) < 1e-3

    @needs_data
    def test_no_hardcoded_anchor_wavelengths_in_module(self):
        # guardrail 1: the module must not embed CO-region solar line wavelengths
        src = (ct.__file__ and open(ct.__file__).read()) or ''
        for w in ('22845.9', '22872.5', '22934.5', '22906.5'):
            assert w not in src, f"hardcoded anchor wavelength {w} found in module"


# ── coadd checksum-dedup (no molecfit) ───────────────────────────────────────
def test_coadd_requires_rest_frame(tmp_path):
    p = tmp_path / 'x.fits'
    p.write_bytes(b'dummy')                       # real file so the md5 dedup runs
    f = ct.CriresFrame(path=p, wlen_id='K2192', band='K', mjd=0.0,
                       date_obs='', ra=0.0, dec=0.0, snr=100.0, specsys='TOPOCENT',
                       fluxcal='UNCALIBRATED', wmin_nm=2000.0, wmax_nm=2400.0,
                       telluric_corrected=True, rest_frame=False,
                       segments=[ct.CriresSegment(4, 3, np.array([22900.0, 22901.0]),
                                                  np.array([1.0, 1.0]),
                                                  np.array([0.1, 0.1]), np.array([0, 0]))])
    with pytest.raises(RuntimeError, match='blind coadd|rest'):
        ct.coadd_co([f])                          # rest_frame=False → loud refuse


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
