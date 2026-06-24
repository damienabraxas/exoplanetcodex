"""
tests/test_telluric_tolerance_rya437.py
=======================================
RYA-437 — the telluric gate is CALIBRATED, not a round number. Part A: the tolerance is
derived from the 13C/CO precision requirement (binding 13CO 2-0), and the stored constant
must equal that derivation (RYA-436 discipline). Part B: the residual metric is CO-region-
local for CO science (global retained as secondary).
"""
import math

import numpy as np
import pytest

from config import constants as C
from pipeline import telluric_tolerance_rya437 as tt
from pipeline import telluric_stage as ts

_LN10 = math.log(10.0)
_HAS_ATLAS = tt._ACE_ATLAS.exists()
needs_atlas = pytest.mark.skipif(not _HAS_ATLAS, reason="ACE-FTS CO atlas not present")


# ── Part A: derivation ────────────────────────────────────────────────────────────
class TestDerivation:
    @needs_atlas
    def test_13co_is_binding(self):
        res = tt.derive_tolerance()
        assert res['binding'].name == '13CO(2-0)'         # weak feature binds
        names = {f.name: f for f in res['features']}
        # 13CO shallower than 12CO -> tighter tolerance
        assert names['13CO(2-0)'].depth < names['12CO(2-0)'].depth
        assert names['13CO(2-0)'].tolerance < names['12CO(2-0)'].tolerance

    @needs_atlas
    def test_formula_r_equals_sigma_ln10_d(self):
        res = tt.derive_tolerance(sigma_dex=0.05)
        for ft in res['features']:
            assert ft.tolerance == pytest.approx(0.05 * _LN10 * ft.depth, rel=1e-9)
            assert ft.sensitivity_dex_per_resid == pytest.approx(1.0 / (ft.depth * _LN10), rel=1e-9)

    @needs_atlas
    def test_stored_constant_matches_derivation(self):
        # RYA-436: the stored constant is the derivation, not a hand-set number.
        res = tt.derive_tolerance()
        assert C.TELLURIC_RESIDUAL_TOL == pytest.approx(res['tolerance_stored'], abs=1e-6)

    @needs_atlas
    def test_bias_at_tol_equals_target(self):
        # by construction, a residual == tol biases A(13C) by ~the target sigma
        res = tt.derive_tolerance(sigma_dex=0.05)
        d13 = res['binding'].depth
        bias = res['tolerance'] / (d13 * _LN10)
        assert bias == pytest.approx(0.05, abs=1e-9)

    @needs_atlas
    def test_calibrated_is_tighter_than_round_05(self):
        assert C.TELLURIC_RESIDUAL_TOL < 0.05         # ~5x tighter
        assert 0.005 < C.TELLURIC_RESIDUAL_TOL < 0.02

    @needs_atlas
    def test_provenance_cites_inputs(self):
        prov = tt.derive_tolerance()['provenance']
        for token in ('13CO(2-0)', 'ACE-FTS', '0.05', 'sigma'):
            assert token in prov
        assert 'TELLURIC_RESIDUAL_TOL_PROVENANCE' in dir(C) or hasattr(C, 'TELLURIC_RESIDUAL_TOL_PROVENANCE')

    @needs_atlas
    def test_bandhead_depths_physical(self):
        d12 = tt.measure_bandhead_depth(tt.CO_12_2_0_BANDHEAD_A)
        d13 = tt.measure_bandhead_depth(tt.CO_13_2_0_BANDHEAD_A)
        assert 0.0 < d13 < d12 < 1.0                  # both absorb; 12CO deeper

    def test_tighter_sigma_gives_tighter_tol(self):
        # monotonic: smaller precision target -> smaller (tighter) tolerance
        if not _HAS_ATLAS:
            pytest.skip("atlas absent")
        lo = tt.derive_tolerance(sigma_dex=0.02)['tolerance']
        hi = tt.derive_tolerance(sigma_dex=0.10)['tolerance']
        assert lo < hi


# ── Part B: CO-region-local metric ─────────────────────────────────────────────────
class TestCoLocalMetric:
    def _segment(self):
        # a segment spanning blue-clean + dirty CO window; telluric-dominated everywhere
        wave = np.linspace(22000, 24000, 4000)
        mt = np.ones_like(wave)
        mt[::3] = 0.5                                  # telluric-dominated pixels
        flux = np.ones_like(wave)                      # clean everywhere...
        co = (wave >= C.TELLURIC_CO_LOCAL_WINDOW_A[0]) & (wave <= C.TELLURIC_CO_LOCAL_WINDOW_A[1])
        flux[co & (mt < 0.9)] = 0.93                   # ...except a 7% mis-correction in the CO band
        return wave, flux, mt

    def test_window_restricts_scoring(self):
        wave, flux, mt = self._segment()
        glob = ts.telluric_residual_metric(wave, flux, mt, window_A=None)
        cow = ts.telluric_residual_metric(wave, flux, mt, window_A=C.TELLURIC_CO_LOCAL_WINDOW_A)
        # global median is diluted by the clean blue end; CO-local sees the dirt
        assert cow['residual'] > glob['residual']
        assert cow['window_A'] == C.TELLURIC_CO_LOCAL_WINDOW_A

    def test_co_local_metric_reports_both_and_verdicts_on_local(self):
        wave, flux, mt = self._segment()
        r = ts.co_local_residual_metric(wave, flux, mt)
        assert r['metric_used'] == 'co_local'
        assert r['residual'] == r['residual_co_local']
        assert r['residual_co_local'] > r['residual_global']   # CO region dirtier than whole
        assert r['passed'] is False                            # 7% >> calibrated tol
        assert r['residual_global'] < r['tol'] or r['residual_global'] >= 0  # global reported

    def test_co_local_clean_passes(self):
        # a CO region corrected to within the calibrated tol passes
        wave = np.linspace(22900, 23500, 2000)
        mt = np.ones_like(wave); mt[::2] = 0.5
        flux = np.ones_like(wave)                      # near-perfect correction
        r = ts.co_local_residual_metric(wave, flux, mt)
        assert r['passed'] is True and r['residual_co_local'] < C.TELLURIC_RESIDUAL_TOL

    def test_manifest_carries_both_residuals(self):
        m = ts.TelluricManifest(
            dataset='d', source_path='/d/x.fits', instrument='CRIRES+', engine='molecfit',
            site='paranal', wave_lo_A=22900, wave_hi_A=24000, regime='ir',
            telluric_required=True, telluric_verified=False, residual=0.032, tolerance=0.0107,
            n_verify_px=121, gdas_profile='C-70.4-24.6D2022-11-25T00.gdas.fits',
            metric_used='co_local', residual_co_local=0.032, residual_global=0.065)
        assert m.analysis_ready is False
        assert m.metric_used == 'co_local'
        assert m.residual_co_local == 0.032 and m.residual_global == 0.065
