import numpy as np
import pytest
from pipeline.bride_input_sensitivity import field_difference, diagnostic_depth_masks


def test_native_depth_scope_includes_exact_boundaries():
    native = np.array([-4., -3., -1., 0., 1.])
    scopes = dict(diagnostic_depth_masks(native))
    np.testing.assert_array_equal(scopes['logtau_-3_0'], [False, True, True, True, False])
    # Same coordinate interpolated with floating arithmetic may fall outside.
    rounded = native.copy()
    rounded[1] = np.nextafter(-3., -np.inf)
    assert np.count_nonzero((rounded >= -3) & (rounded <= 0)) == 2
    assert np.count_nonzero(scopes['logtau_-3_0']) == 3


def test_depth_selection_and_independent_electron_density_ratio():
    baseline = dict(temperature_K=np.array([10., 20.]), electron_pressure_Pa=np.array([2., 4.]), rho_kg_m3=np.array([1., 2.]))
    candidate = dict(temperature_K=np.array([100., 40.]), electron_pressure_Pa=np.array([100., 8.]), rho_kg_m3=np.array([100., 3.]))
    result = field_difference(candidate, baseline, np.array([False, True]))
    assert result == dict(max_abs_temperature_difference_K=20., max_abs_fractional_Pe_difference=1.,
                          max_abs_fractional_rho_difference=.5, max_abs_fractional_ne_difference=0.)
    assert field_difference(candidate, baseline, np.array([True, True]))['max_abs_fractional_ne_difference'] == pytest.approx(4.)


@pytest.mark.parametrize('mask', [[False, False], [1, 0], [True]])
def test_invalid_depth_scope_refused(mask):
    s = {k: np.array([1., 2.]) for k in ('temperature_K', 'electron_pressure_Pa', 'rho_kg_m3')}
    with pytest.raises(ValueError, match='mask'):
        field_difference(s, s, mask)


def test_nonphysical_source_not_hidden_by_depth_selection():
    s = {k: np.array([1., 2.]) for k in ('temperature_K', 'electron_pressure_Pa', 'rho_kg_m3')}
    bad = dict(s, rho_kg_m3=np.array([np.nan, 2.]))
    with pytest.raises(ValueError, match='field'):
        field_difference(bad, s, np.array([False, True]))
