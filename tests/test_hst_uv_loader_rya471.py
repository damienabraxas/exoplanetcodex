"""
tests/test_hst_uv_loader_rya471.py
==================================
RYA-471 — HST STIS/COS UV arm loader. CI-safe: the wiring/guard tests use a
synthetic cspec FITS; the real-Procyon-STIS test is skipped when the RYA-222
whitelist / external data store is absent.
"""
import numpy as np
import pytest
from astropy.io import fits

from pipeline.loaders import hst_uv_loader as L
from pipeline.loaders.base_loader import SpectrumData
from pipeline import uv_conditioning as uvc
from pipeline import cno_synthesis as cs


# ── synthetic cspec helper (one FUV E140M-like frame around C I 1657) ─────────
def _make_cspec(tmp_path, targname='HD61421', lo=1142.0, hi=1709.0, n=4000):
    w = np.linspace(lo, hi, n)
    f = np.ones(n)
    # a Mg II-like core to exercise masking lives in NUV; here add a He II 1640 vac core dip
    f[np.abs(w - 1640.42) < 0.3] = 5.0    # emission core (chromospheric) -> must be masked
    e = np.full(n, 0.01)
    col = fits.ColDefs([
        fits.Column(name='WAVELENGTH', format='E', array=w),
        fits.Column(name='FLUX', format='E', array=f),
        fits.Column(name='ERROR', format='E', array=e)])
    sci = fits.BinTableHDU.from_columns(col, name='SCI')
    pri = fits.PrimaryHDU()
    pri.header['TARGNAME'] = targname
    path = tmp_path / f'{targname}_e140m_cspec.fits'
    fits.HDUList([pri, sci]).writeto(path, overwrite=True)
    return path


# ── vacuum<->air keystone (RYA-426/303) ──────────────────────────────────────
def test_nuv_converts_fuv_stays_vacuum():
    # NUV (Mg II k 2796.3543 vac) converts to air; FUV C I 1657 unchanged
    assert abs(float(uvc.to_pipeline_frame(np.array([2796.3543]))[0]) - 2795.528) < 0.01
    assert abs(float(uvc.to_pipeline_frame(np.array([1657.38]))[0]) - 1657.38) < 1e-9


# ── target guard (never route 55 Cnc as Procyon) ─────────────────────────────
def test_target_guard_accepts_procyon():
    L._assert_procyon('HD61421', 'f.fits')
    L._assert_procyon('alpha CMi', 'f.fits')


def test_target_guard_rejects_55cnc():
    with pytest.raises(L.HstUvProductError):
        L._assert_procyon('55 Cnc', 'srho01cnc.fits')
    with pytest.raises(L.HstUvProductError):
        L._assert_procyon('rho-Cnc', 'lctx.fits')


def test_target_guard_rejects_unknown():
    with pytest.raises(L.HstUvProductError):
        L._assert_procyon('SOME OTHER STAR', 'f.fits')


# ── grating guards ────────────────────────────────────────────────────────────
def test_excluded_grating_raises(tmp_path):
    p = _make_cspec(tmp_path)
    with pytest.raises(L.HstUvProductError):
        L.HstUvLoader(p, grating='WFC3').load()


def test_unknown_grating_raises(tmp_path):
    p = _make_cspec(tmp_path)
    with pytest.raises(L.HstUvProductError):
        L.HstUvLoader(p, grating='NOPE').load()


# ── load() on a synthetic FUV frame ───────────────────────────────────────────
def test_load_synthetic_fuv_frame(tmp_path):
    p = _make_cspec(tmp_path)
    spec = L.HstUvLoader(p, grating='E140M').load()
    assert isinstance(spec, SpectrumData)
    assert spec.meta['provenance'] == 'measured'
    assert spec.meta['regime'] == 'FUV'
    assert spec.meta['ew_allowed'] is False           # FUV synthesis-not-EW
    assert spec.meta['telluric_corrected'] is False    # space UV
    assert spec.meta['native_frame'] == 'vacuum'
    lo, hi = spec.wave_range_A
    assert lo <= 1657.38 <= hi                         # C I 1657 covered
    # chromospheric He II 1640 core masked -> NaN present
    assert np.isnan(spec.flux).any()


def test_load_rejects_55cnc_frame(tmp_path):
    p = _make_cspec(tmp_path, targname='srho01cnc')
    with pytest.raises(L.HstUvProductError):
        L.HstUvLoader(p, grating='E140M').load()


# ── pick the merged coadd over a raveled echelle x1d ─────────────────────────
def test_pick_arm_frame_prefers_cspec():
    rows = [
        {'filename': 'a_x1d.fits', 'filepath': '/x/a_x1d.fits', 'wl_min_A': '1140', 'wl_max_A': '1730'},
        {'filename': 'b_cspec.fits', 'filepath': '/x/b_cspec.fits', 'wl_min_A': '1142', 'wl_max_A': '1709'}]
    assert L._pick_arm_frame(rows)['filename'] == 'b_cspec.fits'   # cspec even though slightly narrower


# ── UV diagnostics wired from RYA-190 ─────────────────────────────────────────
def test_uv_arm_diagnostics_from_rya190():
    diags = L.uv_arm_diagnostics()
    keys = [d.key for d in diags]
    assert 'CI_1657_fuv' in keys
    # exactly the FUV verdict=USE synthesis lines (C I 1657/1931, O I 1356, S I 1474)
    assert len(diags) == 4
    ci = next(d for d in diags if d.key == 'CI_1657_fuv')
    assert ci.role == 'primary' and ci.kind == 'atomic'
    assert 'grid_owed' in ci.nlte_flag and 'RYA-359' in ci.nlte_ref   # NLTE flagged loud, not applied


# ── RYA-464 registry + resolver wiring ────────────────────────────────────────
def test_registry_uv_arm_built_but_deferred():
    uv = cs.star_arm_registry('procyon')['uv']
    assert uv.loader == 'hst_stis'
    # honest: RYA-348 Phase 3 ran the synthesis and found it not science-grade
    assert uv.ready is False
    assert len(uv.diagnostics) == 4                   # diagnostics ARE wired
    assert 'RYA-471' in uv.defer_reason               # loader provenance still cited


def test_resolver_deferred_arm_raises():
    uv = cs.star_arm_registry('procyon')['uv']
    with pytest.raises(cs.ArmNotWired):
        cs.resolve_arm_spectrum('procyon', uv)         # ready=False -> loud, never silent


def test_resolver_hst_stis_is_procyon_only():
    uv = cs.star_arm_registry('procyon')['uv']
    forced = cs.ArmWiring(uv.name, uv.region, uv.diagnostics, uv.loader, True)
    with pytest.raises(cs.ArmNotWired):
        cs.resolve_arm_spectrum('solar', forced)       # UV loader is Procyon-only today


# ── real Procyon STIS (skipped without the RYA-222 whitelist / data store) ───
@pytest.mark.skipif(not L.SCIENCE_READY_CSV.exists(),
                    reason='RYA-222 Procyon HST whitelist / external data store absent')
def test_real_procyon_e140m_loads_and_conditions():
    rows = L.science_ready_rows(grating='E140M')
    if not rows:
        pytest.skip('no on-disk whitelisted E140M frame')
    spec = L.load_procyon_uv_arm('E140M', return_spectrum=True)
    assert spec.meta['object'].replace(' ', '').upper().startswith('HD61421')
    lo, hi = spec.wave_range_A
    assert lo <= 1657.38 <= hi                         # measured-C Procyon advantage
    assert spec.meta['analysis_ready'] is True
    assert spec.meta['ew_allowed'] is False
    # resolver returns (wave_nm, flux) when the arm is forced ready
    uv = cs.star_arm_registry('procyon')['uv']
    forced = cs.ArmWiring(uv.name, uv.region, uv.diagnostics, uv.loader, True)
    w_nm, flux = cs.resolve_arm_spectrum('procyon', forced)
    assert 114.0 < float(np.nanmin(w_nm)) < 171.0
