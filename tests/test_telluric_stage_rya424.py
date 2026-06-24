"""
tests/test_telluric_stage_rya424.py
===================================
RYA-424 — the STANDING, instrument-aware telluric data-input stage. The
molecfit-independent core is fully tested here (wavelength gate, engine routing,
GDAS loud-fail, residual metric, conditioning manifest, abundance-path guard); the
real molecfit pass on the Vesta CRIRES+ set is data-/engine-gated and skipped when
absent. Every engine-gated step FAILS LOUD (never a silent skip / fake correction).
"""
import json
from pathlib import Path

import numpy as np
import pytest

from config import constants as C
from pipeline import telluric_stage as ts


# ── 1. Wavelength gate (the rule is wavelength-gated, not instrument-gated) ────────
class TestWavelengthGate:
    def test_ir_span_required(self):
        # fully red of 6800 Å → IR forest → mandatory
        assert ts.telluric_regime(22900, 24000) == C.TELLURIC_REGIME_IR
        assert ts.requires_telluric(22900, 24000) is True

    def test_red_optical_straddle_required(self):
        # crosses 6800 Å → red_optical (has forest in its red end) → mandatory
        assert ts.telluric_regime(6500, 9500) == C.TELLURIC_REGIME_RED_OPTICAL
        assert ts.requires_telluric(6500, 9500) is True

    def test_red_optical_window_required(self):
        # the O I 7772 / K I 7699 / N I 8216 red-optical window (RYA-380 explicit)
        assert ts.requires_telluric(6910, 9500) is True

    def test_mid_optical_not_auto_required(self):
        assert ts.telluric_regime(4000, 6000) == C.TELLURIC_REGIME_MID_OPTICAL
        assert ts.requires_telluric(4000, 6000) is False

    def test_blue_not_required(self):
        assert ts.telluric_regime(3000, 3500) == C.TELLURIC_REGIME_BLUE
        assert ts.requires_telluric(3000, 3500) is False

    def test_space_uv_never_required(self):
        # HST/STIS/COS sit blueward but are above the atmosphere
        assert ts.requires_telluric(1200, 1700, space_uv=True) is False

    def test_threshold_at_6800(self):
        # exactly at the threshold the span is in the forest (≥ min)
        assert ts.requires_telluric(C.TELLURIC_LAMBDA_MIN_A, 7000) is True

    def test_unordered_span(self):
        # hi<lo is tolerated (sorted internally)
        assert ts.telluric_regime(24000, 22900) == C.TELLURIC_REGIME_IR


# ── 2. Engine routing (single source of truth; no silent default) ─────────────────
class TestEngineRouting:
    def test_molecfit_instruments(self):
        for inst in ('CRIRES', 'CRIRES+', 'UVES-red', 'ESPRESSO-red', 'FEROS-red',
                     'CHIRON', 'NIRPS'):
            assert ts.select_engine(inst) == C.TELLURIC_ENGINE_MOLECFIT

    def test_spirou_is_apero_wapiti_not_molecfit(self):
        # the SPIRou permanent rule: APERO + Wapiti, NOT molecfit
        assert ts.select_engine('SPIRou') == C.TELLURIC_ENGINE_APERO_WAPITI

    def test_unknown_instrument_loud_fails(self):
        with pytest.raises(ts.EngineNotSelectedError):
            ts.select_engine('HARPS-Frankenstein')

    def test_site_routing(self):
        assert ts.site_for_instrument('CRIRES+')['key'] == 'paranal'
        assert ts.site_for_instrument('NIRPS')['key'] == 'la_silla'
        assert ts.site_for_instrument('CHIRON')['key'] == 'ctio'
        assert ts.site_for_instrument('SPIRou')['key'] == 'cfht'

    def test_site_unknown_loud_fails(self):
        with pytest.raises(ts.EngineNotSelectedError):
            ts.site_for_instrument('Nonexistent')

    def test_every_engine_has_a_site(self):
        # routing integrity: every engine-registered instrument has a GDAS/observatory site
        for inst in C.TELLURIC_ENGINES:
            assert inst in C.TELLURIC_INSTRUMENT_SITE
            assert ts.site_for_instrument(inst)['key'] in C.SITES


# ── 3. Per-night GDAS (LOUD-FAIL on silent fallback) ──────────────────────────────
class TestGDAS:
    def test_nan_mjd_loud_fails(self, tmp_path):
        with pytest.raises(ts.GDASUnavailableError):
            ts.resolve_gdas_profile(float('nan'), 'CRIRES+', tmp_path)

    def test_missing_site_tarball_loud_fails(self):
        # La Silla GDAS tarball is not installed → loud-fail, never standard atmosphere
        with pytest.raises(ts.GDASUnavailableError):
            ts.gdas_tarball_for_site(ts.site_for_instrument('NIRPS'))

    def test_non_molecfit_site_has_no_gdas(self):
        # CFHT/SPIRou uses APERO+Wapiti, not molecfit GDAS → loud-fail if asked
        with pytest.raises(ts.GDASUnavailableError):
            ts.gdas_tarball_for_site(ts.site_for_instrument('SPIRou'))

    def test_assert_real_gdas_passes_real_profile(self):
        ts.assert_real_gdas('C-70.4-24.6D2022-11-23T00.gdas.fits', 'CRIRES+')

    @pytest.mark.parametrize('bad', ['standard-profile (no GDAS)', '', None,
                                     'standard atm fallback'])
    def test_assert_real_gdas_rejects_standard_profile(self, bad):
        with pytest.raises(ts.GDASUnavailableError):
            ts.assert_real_gdas(bad, 'CRIRES+')

    @pytest.mark.skipif(ts._telluriccorr_share() is None,
                        reason='telluriccorr GDAS data not installed')
    def test_paranal_real_gdas_resolves(self, tmp_path):
        # the Paranal tarball IS installed → a real per-night profile resolves to FITS
        g = ts.resolve_gdas_profile(59500.5, 'CRIRES+', tmp_path)
        assert Path(g).exists() and g.endswith('.fits')
        ts.assert_real_gdas(Path(g).name, 'CRIRES+')   # and it is NOT the standard sentinel


# ── 4. Verification metric (residuals in known telluric windows) ──────────────────
class TestResidualMetric:
    def _band(self):
        wave = np.linspace(22900, 23100, 400)
        mt = np.ones_like(wave)
        mt[100:200] = 0.5                      # telluric-dominated, not saturated
        return wave, mt

    def test_clean_correction_passes(self):
        wave, mt = self._band()
        flux = np.ones_like(wave)              # perfectly returns to continuum
        r = ts.telluric_residual_metric(wave, flux, mt)
        assert r['passed'] is True and r['residual'] < 1e-6 and r['n_px'] == 100

    def test_residual_above_tol_fails(self):
        wave, mt = self._band()
        flux = np.ones_like(wave)
        flux[100:200] = 0.8                     # 20% under-correction in the telluric band
        r = ts.telluric_residual_metric(wave, flux, mt)
        assert r['passed'] is False and r['residual'] > C.TELLURIC_RESIDUAL_TOL

    def test_too_few_pixels_not_silently_passed(self):
        wave = np.linspace(22900, 23100, 400)
        mt = np.ones_like(wave)                 # no telluric-dominated pixels at all
        r = ts.telluric_residual_metric(wave, np.ones_like(wave), mt)
        assert r['passed'] is False and 'too few' in r['reason']

    def test_science_mask_excludes_stellar_lines(self):
        # stellar signal over the MAJORITY of the telluric band must NOT be scored as
        # telluric misfit (the metric is a median → contamination must exceed half the
        # telluric pixels to move it; that is exactly when the mask must rescue it).
        wave, mt = self._band()
        flux = np.ones_like(wave)
        flux[100:175] = 0.7                     # stellar signal over 75 of 100 telluric px
        mask = np.zeros_like(wave, bool)
        mask[100:175] = True
        masked = ts.telluric_residual_metric(wave, flux, mt, science_mask=mask)
        unmasked = ts.telluric_residual_metric(wave, flux, mt)
        assert masked['residual'] < unmasked['residual']
        assert masked['passed'] is True and unmasked['passed'] is False


# ── Conditioning manifest + abundance-path guard ──────────────────────────────────
def _manifest(**kw):
    base = dict(dataset='vesta_crires_k', source_path='/d/v.fits', instrument='CRIRES+',
                engine='molecfit', site='paranal', wave_lo_A=22900.0, wave_hi_A=24000.0,
                regime='ir', telluric_required=True, telluric_verified=True, residual=0.03,
                tolerance=0.05, n_verify_px=120,
                gdas_profile='C-70.4-24.6D2022-11-23T00.gdas.fits', mjd=59906.0)
    base.update(kw)
    return ts.TelluricManifest(**base)


class TestManifestAndGuard:
    def test_verified_ir_is_analysis_ready(self):
        assert _manifest().analysis_ready is True

    def test_unverified_ir_not_analysis_ready(self):
        assert _manifest(telluric_verified=False, residual=None).analysis_ready is False

    def test_blue_unverified_still_analysis_ready(self):
        m = _manifest(telluric_required=False, telluric_verified=False,
                      regime='blue', wave_lo_A=3000, wave_hi_A=3500)
        assert m.analysis_ready is True

    def test_guard_passes_verified(self):
        ts.require_telluric_verified(_manifest())

    def test_guard_loud_fails_unverified(self):
        with pytest.raises(ts.TelluricNotVerifiedError):
            ts.require_telluric_verified(_manifest(telluric_verified=False))

    def test_guard_passes_not_required(self):
        ts.require_telluric_verified(_manifest(telluric_required=False,
                                               telluric_verified=False))

    def test_guard_accepts_dict(self):
        d = {'telluric_required': True, 'telluric_verified': False, 'dataset': 'd'}
        with pytest.raises(ts.TelluricNotVerifiedError):
            ts.require_telluric_verified(d)

    def test_guard_span_crosscheck_overrides_stale_flag(self):
        # manifest claims not-required, but an IR span is supplied → must still fail
        d = {'telluric_required': False, 'telluric_verified': False, 'dataset': 'd'}
        with pytest.raises(ts.TelluricNotVerifiedError):
            ts.require_telluric_verified(d, wave_lo_A=23000, wave_hi_A=24000)

    def test_manifest_roundtrip(self, tmp_path):
        m = _manifest(telluric_verified=False, residual=None)
        p = ts.write_manifest(m, tmp_path / 'x.telluric.json')
        d = ts.load_manifest(p)
        assert d['analysis_ready'] is False
        assert d['engine'] == 'molecfit' and d['site'] == 'paranal'
        assert d['gdas_profile'].startswith('C-70.4-24.6')
        # the path form of the guard reads the same record
        with pytest.raises(ts.TelluricNotVerifiedError):
            ts.require_telluric_verified(p)


# ── Engine availability (fail-loud install reporting) ─────────────────────────────
class TestAvailability:
    def test_molecfit_reported(self):
        # boolean either way; on this box molecfit is installed
        assert isinstance(ts.engine_available(C.TELLURIC_ENGINE_MOLECFIT), bool)

    def test_apero_wapiti_reported(self):
        assert isinstance(ts.engine_available(C.TELLURIC_ENGINE_APERO_WAPITI), bool)

    def test_unknown_engine_not_available(self):
        assert ts.engine_available('made-up-engine') is False


# ── Vesta CRIRES+ end-to-end (data + molecfit gated; the acceptance smoke) ─────────
try:
    from pipeline import crires_telluric as _ct
    _HAS_DATA = _ct.VESTA_CRIRES_DIR.exists() and any(_ct.VESTA_CRIRES_DIR.glob('*.fits'))
except Exception:
    _HAS_DATA = False

_HAS_MOLECFIT = ts.engine_available(C.TELLURIC_ENGINE_MOLECFIT)
e2e = pytest.mark.skipif(not (_HAS_DATA and _HAS_MOLECFIT),
                         reason='Vesta CRIRES+ data and/or molecfit engine not present')


class TestVestaEndToEnd:
    @e2e
    def test_condition_one_frame_produces_manifest(self, tmp_path):
        frames = _ct.co_overtone_frames(on_chip_only=True)
        f = max(frames, key=lambda fr: fr.snr)
        man = ts.condition_crires_frame(f, tmp_path / f.wlen_id,
                                        dataset='vesta_crires_k', instrument='CRIRES+')
        # routed, corrected, real per-night GDAS, residual scored
        assert man.engine == 'molecfit' and man.site == 'paranal' and man.regime == 'ir'
        assert man.telluric_required is True
        ts.assert_real_gdas(man.gdas_profile, 'CRIRES+')   # NOT a standard-profile fallback
        assert man.residual is not None and man.n_verify_px > 0
        # honest verdict: verified iff residual within tolerance (provisional gate)
        assert man.telluric_verified == (man.residual <= man.tolerance)
