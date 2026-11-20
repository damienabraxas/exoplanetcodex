"""RYA-1282: the Bride interpolates collision tables exactly as MULTI 2.3 gencol does.

pipeline.multi_tautsp ports TAUTSP/PPVALU (and SPLIN for N < 4). The pinned native Fortran
oracle (bride_native_spline.NativeTautsp, mul23_opacu.f 4d6aa3c) is the referee where it
is held; oracle-free invariants run everywhere. Measured agreement: 1.4e-14 over all 1106
O41 rows, 4.5e-15 over 3000 random splines. Bound 1e-12 = measured x ~70: rounding only.

Negative control: the pre-RYA-1282 default (Lightweaver weno4 on the LINEAR rate) differs
from MULTI by up to 42% (CI) on the held atom; the oracle test must reject it.
"""
import numpy as np
import pytest

from pipeline.multi_tautsp import gencol_coefficient, ppvalu, splin_small, tautsp


def _oracle():
    try:
        from config.constants import require_codex_path
        from pipeline.bride_native_spline import NativeTautsp
        return NativeTautsp(require_codex_path('experiments.bride_multi_spline_oracle'))
    except Exception as e:                                        # noqa: BLE001
        pytest.skip(f'native spline oracle not held here: {e}')


def _o41():
    try:
        from config.constants import require_codex_path
        from pipeline.multi_atom import read_multi_atom
        return read_multi_atom(require_codex_path('experiments.bride_oxygen_atom'))
    except Exception as e:                                        # noqa: BLE001
        pytest.skip(f'O41 atom not held here: {e}')


def _native(orc, kind, g, v, t):
    out = orc.interpolate(g, v, t, logarithmic=kind in ('CE', 'CI', 'CH0'), allow_negative=kind == 'CH')
    return np.maximum(out, 0) if kind == 'CH' else out


def test_every_held_O41_row_matches_the_native_fortran():
    orc, atom = _oracle(), _o41()
    worst = 0.
    for raw in atom.collisions:
        g = np.asarray(raw.temperature, float); v = np.asarray(raw.rates, float)
        t = np.concatenate([np.linspace(g[0], g[-1], 499), g])
        ref = _native(orc, raw.kind, g, v, t); got = gencol_coefficient(raw.kind, g, v, t)
        nz = ref != 0
        assert np.all(got[~nz] == 0)
        worst = max(worst, float(np.max(abs(got[nz] / ref[nz] - 1))))
    assert worst < 1e-12


def test_random_splines_through_every_branch_match_the_native_fortran():
    orc = _oracle(); rng = np.random.default_rng(1282); worst = 0.
    for _ in range(1500):
        k = int(rng.integers(4, 14)); g = np.cumsum(rng.uniform(.1, 5, k))
        v = rng.normal(0, 1, k) * rng.choice([1, 1e-3, 1e3])
        if rng.random() < .3:
            v[rng.integers(0, k)] = v[rng.integers(0, k)]           # flat pieces: z = 0 / 1 branches
        t = np.linspace(g[0], g[-1], 61)
        ref = orc.interpolate(g, v, t, logarithmic=False, allow_negative=True)
        worst = max(worst, float(np.max(abs(ppvalu(*tautsp(g, v), t) - ref)) / max(np.max(abs(ref)), 1e-300)))
    assert worst < 1e-12


def test_negative_control_the_old_weno4_default_is_rejected():
    orc, atom = _oracle(), _o41()
    from lightweaver.utils import weno4
    worst = 0.
    for raw in atom.collisions:
        if raw.kind != 'CI':
            continue
        g = np.asarray(raw.temperature, float); v = np.asarray(raw.rates, float)
        t = np.linspace(g[0], g[-1], 499)
        worst = max(worst, float(np.max(abs(weno4(t, g, v) / _native(orc, 'CI', g, v, t) - 1))))
    assert worst > 0.1                                             # measured 0.42


# ---- oracle-free invariants (run on every host) ----------------------------------------
def test_spline_passes_through_its_nodes_and_reproduces_a_straight_line():
    g = np.array([1000., 2000., 3500., 5000., 8000., 10000.])
    v = np.array([3., -1., 4., 1., 5., 9.])
    np.testing.assert_allclose(ppvalu(*tautsp(g, v), g), v, rtol=0, atol=1e-12)
    line = 2.5 * g - 7.
    t = np.linspace(g[0], g[-1], 101)
    np.testing.assert_allclose(ppvalu(*tautsp(g, line), t), 2.5 * t - 7., rtol=1e-13)


def test_gencol_clamps_temperature_and_interpolates_log_for_ce_ci_ch0():
    g = np.array([1000., 2000., 4000., 6000., 8000., 10000.])
    v = 1e-9 * (g / 1000.) ** 0.7
    for kind in ('CE', 'CI', 'CH0'):
        out = gencol_coefficient(kind, g, v, np.array([500., 1000., 10000., 12000.]))
        np.testing.assert_allclose(out[[0, 3]], v[[0, -1]], rtol=1e-12)        # held, not extrapolated
        np.testing.assert_allclose(gencol_coefficient(kind, g, v, g), v, rtol=1e-12)
    expo = 1e-9 * np.exp(-g / 3000.)                                            # log-linear in T
    t = np.linspace(g[0], g[-1], 77)
    np.testing.assert_allclose(gencol_coefficient('CE', g, expo, t), 1e-9 * np.exp(-t / 3000.), rtol=1e-12)


def test_ch_is_linear_and_clipped_at_zero():
    g = np.array([1000., 2000., 3000., 4000., 5000.])
    v = np.array([1e-15, 0., 0., 0., 5e-14])                     # spline undershoots between
    out = gencol_coefficient('CH', g, v, np.linspace(1000., 5000., 401))
    assert np.all(out >= 0)
    with pytest.raises(ValueError):
        gencol_coefficient('CE', g, v, g)                          # log kind refuses zeros


def test_splin_small_tables_follow_multi():
    x2, y2 = np.array([1000., 3000.]), np.array([2., 6.])
    np.testing.assert_allclose(splin_small(np.array([0., 2000., 5000.]), x2, y2), [0., 4., 10.])  # extrapolates
    np.testing.assert_array_equal(splin_small(np.array([1., 9e9]), x2[:1], y2[:1]), [2., 2.])
    from scipy.interpolate import CubicSpline
    x3, y3 = np.array([1000., 2500., 7000.]), np.array([1., 4., 2.])
    t = np.linspace(1000., 7000., 121)
    np.testing.assert_allclose(splin_small(t, x3, y3), CubicSpline(x3, y3, bc_type='natural')(t), rtol=0, atol=1e-7)
    np.testing.assert_array_equal(splin_small(np.array([10., 9000.]), x3, y3), [1., 2.])           # held outside
