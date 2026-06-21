"""
tests/test_reflected_solar_rv.py
================================
Tests for pipeline/reflected_solar_rv.py (RYA-372). These cover the pure-numeric
core WITHOUT the network (JPL Horizons) or the external Vesta FITS data: line-core
velocity recovery (incl. the large-shift two-pass), robust median sigma-clipping,
the Horizons-cache path, and the status logic. The full per-frame conditioning is
exercised by the smoke test (`--set vesta_espresso/--verify`, output in the Linear
comment), which needs the out-of-repo data + Horizons.
"""
import sys

import numpy as np
import pytest

from pipeline import reflected_solar_rv as rsv

C = rsv.C_KMS


def _synth_spectrum(lines, v_kms=0.0, depth=0.4, fwhm_A=0.12, step=0.01, pad=2.0):
    """Continuum-1.0 spectrum with Gaussian absorption at `lines` (air Å), the whole
    set Doppler-shifted by v_kms. Returns (wave_air, flux)."""
    lo, hi = min(lines) - pad, max(lines) + pad
    w = np.arange(lo, hi, step)
    f = np.ones_like(w)
    sig = fwhm_A / 2.3548
    for x in lines:
        xs = x * (1.0 + v_kms / C)
        f -= depth * np.exp(-0.5 * ((w - xs) / sig) ** 2)
    return w, f


class TestCoreVelocity:
    def test_recovers_small_shift(self):
        w, f = _synth_spectrum([5269.537], v_kms=-4.0)
        v = rsv._core_velocity(w, f, 5269.537, win=0.22)
        assert abs(v - (-4.0)) < 0.3

    def test_flat_continuum_returns_nan(self):
        w = np.arange(5200, 5300, 0.01)
        f = np.ones_like(w)
        assert np.isnan(rsv._core_velocity(w, f, 5269.537))

    def test_center_offset_finds_shifted_line(self):
        # A +30 km/s line is ~0.53 Å red of rest at 5270 — outside a ±0.22 window
        # unless `center` is moved to the predicted position.
        w, f = _synth_spectrum([5269.537], v_kms=30.0)
        assert np.isnan(rsv._core_velocity(w, f, 5269.537, win=0.22))      # missed
        c = 5269.537 * (1 + 30.0 / C)
        v = rsv._core_velocity(w, f, 5269.537, win=0.22, center=c)
        assert abs(v - 30.0) < 0.3                                          # found


class TestBulkVelocity:
    def test_recovers_small_shift(self):
        w, f = _synth_spectrum(rsv.VEL_LINES_AIR, v_kms=-4.5)
        r = rsv.measure_bulk_velocity(w, f)
        assert r['n_used'] >= rsv.MIN_LINES
        assert abs(r['v_med'] - (-4.5)) < 0.2

    def test_two_pass_recovers_large_shift(self):
        # +30 km/s would defeat a single narrow-window search; the coarse→fine
        # two-pass must still recover it.
        w, f = _synth_spectrum(rsv.VEL_LINES_AIR, v_kms=30.0)
        r = rsv.measure_bulk_velocity(w, f)
        assert abs(r['v_med'] - 30.0) < 0.2

    def test_insufficient_when_no_lines(self):
        w = np.arange(7000, 7100, 0.01)        # no VEL_LINES here
        f = np.ones_like(w)
        r = rsv.measure_bulk_velocity(w, f)
        assert np.isnan(r['v_med']) and r['n_used'] < rsv.MIN_LINES

    def test_train_test_disjoint(self):
        assert not (set(rsv.TRAIN_LINES) & set(rsv.TEST_LINES))


class TestRobustMedian:
    def test_clips_outlier(self):
        v = np.array([10.0, 10.1, 9.9, 10.05, 40.0])    # one gross outlier
        med, std, n = rsv._robust_median(v)
        assert abs(med - 10.0) < 0.2 and n == 4

    def test_empty(self):
        med, std, n = rsv._robust_median(np.array([np.nan, np.nan]))
        assert np.isnan(med) and n == 0


class TestEphemerisPath:
    def test_uses_cache_no_network(self):
        # Seed the cache so reflected_solar_rv resolves without hitting Horizons.
        # RYA-394: key is (mjd, body_key, obs_code); value carries the targetname.
        mjd = 12345.678
        rsv._HZ_CACHE[(round(mjd, 6), 'vesta', rsv.PARANAL)] = (1.5, 12.0, '4 Vesta (A807 FA)')
        out = rsv.reflected_solar_rv(mjd, 'vesta')
        assert out['v_helio'] == 1.5 and out['v_obs'] == 12.0
        assert out['v_total'] == pytest.approx(13.5)    # two-leg sum
        assert 'Vesta' in out['targetname']

    def test_paranal_obscode_and_default_body_key(self):
        # the default is a registry KEY, never a raw Horizons id (no bare '4'=Mars)
        assert rsv.PARANAL == '309' and rsv.DEFAULT_BODY == 'vesta'
        assert not hasattr(rsv, 'VESTA_ID')


class TestBodyIDGuard:
    """RYA-394 — the body-ID bug must be structurally impossible."""

    def test_unknown_body_key_raises(self):
        with pytest.raises(rsv.BodyIDError):
            rsv.reflected_solar_rv(60000.0, 'nibiru')

    def test_raw_horizons_id_rejected(self):
        # passing a raw id (the old bug surface) is an unknown KEY → BodyIDError,
        # never a silent Horizons lookup that could return Mars.
        with pytest.raises(rsv.BodyIDError):
            rsv.reflected_solar_rv(60000.0, '4')

    def test_registry_matches_intended_bodies(self):
        from config.constants import REFLECTED_SOLAR_BODIES as REG
        assert REG['vesta']['id_type'] == 'smallbody' and REG['vesta']['match'] == 'Vesta'
        assert REG['iris']['id_type'] == 'smallbody'
        assert REG['ganymede']['id'] == '503' and REG['ganymede']['id_type'] is None

    def test_guard_raises_on_target_mismatch(self, monkeypatch):
        # simulate Horizons resolving 'vesta' to Mars Barycenter → guard must loud-fail
        import types

        class _Eph:
            def __getitem__(self, k):
                return {'targetname': ['Mars Barycenter (4)'],
                        'r_rate': [2.3], 'delta_rate': [-3.3]}[k]

        class _H:
            def __init__(self, *a, **k):
                pass

            def ephemerides(self):
                return _Eph()

        monkeypatch.setitem(sys.modules, 'astroquery.jplhorizons',
                            types.SimpleNamespace(Horizons=_H))
        rsv._HZ_CACHE.clear()
        with pytest.raises(rsv.BodyIDError):
            rsv.reflected_solar_rv(60001.234, 'vesta')


class TestRestFrameAssert:
    """RYA-394 — the closed loop that was left report-only."""

    def test_assert_passes_within_tol(self):
        rsv.assert_rest_frame(0.2, 'test line')          # no raise

    def test_assert_loud_fails_off_rest(self):
        with pytest.raises(rsv.RestFrameError):
            rsv.assert_rest_frame(25.0, 'Vesta-as-Mars line')   # the Mars-swap signature

    def test_nonfinite_is_not_a_rest_failure(self):
        rsv.assert_rest_frame(float('nan'), 'unmeasurable')      # insufficient ≠ off-rest


def _write_fits(tmp_path, pro_catg, instrument='ESPRESSO', wave_air=None):
    from astropy.io import fits
    if wave_air is None:
        wave_air = np.linspace(5880.0, 5900.0, 50)
    ph = fits.PrimaryHDU()
    ph.header['INSTRUME'] = instrument
    ph.header['HIERARCH ESO PRO CATG'] = pro_catg
    ph.header['SPECSYS'] = 'BARYCENT'
    ph.header['MJD-OBS'] = 58425.0
    ph.header['EXPTIME'] = 600.0
    ph.header['SNR'] = 250.0
    cols = [fits.Column(name='WAVE', format='D', array=wave_air * 1.000403),  # fake vacuum
            fits.Column(name='WAVE_AIR', format='D', array=wave_air),
            fits.Column(name='FLUX', format='D', array=np.ones_like(wave_air))]
    tb = fits.BinTableHDU.from_columns(cols)
    path = tmp_path / f'{pro_catg}.fits'
    fits.HDUList([ph, tb]).writeto(path, overwrite=True)
    return path


class TestFrameLoading:
    def test_stack_a_rejected_with_reason(self, tmp_path):
        p = _write_fits(tmp_path, 'S1D_STACK_A')
        fr = rsv._load_frame(p)
        assert 'reject' in fr and 'STACK' in fr['reject']

    def test_final_a_loads_wave_air(self, tmp_path):
        p = _write_fits(tmp_path, 'S1D_FINAL_A')
        fr = rsv._load_frame(p)
        assert 'reject' not in fr
        assert 'wave_air' in fr and fr['wave_air'][0] == pytest.approx(5880.0, abs=1e-3)
        # WAVE_AIR (air) chosen over WAVE (vacuum), so the first sample is the air value
        assert fr['instrument'] == 'ESPRESSO'

    def test_units_trap_reads_angstrom_column(self, tmp_path):
        # The data column is Å (~5880), not the nm header unit — loader must return Å.
        p = _write_fits(tmp_path, 'S1D_FINAL_A')
        fr = rsv._load_frame(p)
        assert fr['wave_air'].min() > 1000   # Å, not nm


if __name__ == '__main__':
    raise SystemExit(pytest.main([__file__, '-v']))
