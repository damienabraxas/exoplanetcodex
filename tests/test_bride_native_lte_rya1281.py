"""Actual native LTEPOP checked against independent two-stage equilibrium."""
import numpy as np
import pytest

from pipeline.bride_native_lte import NativeLte


@pytest.fixture
def oracle():
    from config.constants import require_codex_path
    path = require_codex_path('experiments.bride_multi_lte_oracle')
    if not (path / 'manifest.json').exists():
        pytest.skip('Pinned native LTE oracle not installed')
    return NativeLte(path)


def test_native_saha_boltzmann_and_debye_lowering(oracle):
    # Non-square layout catches transposition; excited neutral and ion test both
    # Boltzmann excitation and electron phase-space normalization independently.
    T = np.array([4000., 6000., 9000., 12000.])
    ne = np.array([1e16, 1e17, 1e18, 1e19])
    E = np.array([0., 10000., 100000.])
    g = np.array([5., 3., 4.])
    actual, c = oracle.fractions(T, ne, E, g, [1, 1, 2])
    k, h, me = c['k_erg_K'], c['h_erg_s'], c['electron_g']
    # These coefficients are native REAL literals before DP promotion.
    # The constant-only product is evaluated in default REAL before division
    # by the DP temperature, as in the unchanged Fortran expression.
    coefficient = float(np.float32(4.98e-4) * np.float32(5040.))
    lowering = coefficient / T * np.sqrt(ne * 1e-6 * k * T)
    excitation = c['h_erg_s'] * c['c_cm_s'] * E / c['eV_erg']
    expected = np.ones((3, len(T)))
    expected[1] = g[1] / g[0] * np.exp(-c['eV_over_k_K'] * excitation[1] / T)
    expected[2] = (2 * g[2] / g[0] / (ne * 1e-6)
        * (2 * c['pi'] * me * k * T / h**2)**1.5
        * np.exp(-c['eV_over_k_K'] * (excitation[2] - lowering) / T))
    expected /= expected.sum(axis=0)
    np.testing.assert_allclose(actual, expected, rtol=3e-14)
    np.testing.assert_allclose(actual.sum(axis=0), 1., rtol=2e-15)
    # Repeated calls must reset common-block populations and dimensions.
    one, _ = oracle.fractions([6000.], [1e17], [0.], [5.], [1])
    np.testing.assert_array_equal(one, [[1.]])
    again, _ = oracle.fractions(T, ne, E, g, [1, 1, 2])
    np.testing.assert_array_equal(again, actual)


@pytest.mark.parametrize('change', [
    {'temperature': [0.]}, {'ne_m3': [-1.]}, {'energies_cm': [1.]},
    {'weights': [0.]}, {'stages': [1.5]}, {'stages': [3]},
    {'temperature': [6000.] * 401, 'ne_m3': [1e17] * 401},
    {'energies_cm': [0.] * 65, 'weights': [1.] * 65, 'stages': [1] * 65},
])
def test_reject_before_entering_native_code(oracle, change):
    args = dict(temperature=[6000.], ne_m3=[1e17], energies_cm=[0.], weights=[5.], stages=[1])
    args.update(change)
    with pytest.raises(ValueError):
        oracle.fractions(**args)
