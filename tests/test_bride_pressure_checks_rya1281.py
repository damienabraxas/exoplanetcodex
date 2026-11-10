import numpy as np
import pytest
from pipeline.bride_pressure_checks import pressure_column_metrics, pressure_column_candidate, read_pressure_table


def test_hydrostatic_column_with_independent_nonzero_pressure_intercept():
    m = np.array([2., 4., 7., 11.]); g = 3.; p = 99.+g*m
    metrics = pressure_column_metrics(p, m, g)
    assert max(metrics.values()) == 0.
    np.testing.assert_array_equal(pressure_column_candidate(p, g, m[0]), m)
    # Absolute P/g would be wrong because of the unknown intercept.
    assert not np.allclose(p/g, m)
    assert pressure_column_metrics(p, m, g*2)['max_abs_adjacent_hydrostatic_difference'] == .5


@pytest.mark.parametrize('m,g', [([1, 1, 2], 1), ([1, 2, 3], 0), ([1, np.nan, 3], 1)])
def test_bad_mass_or_gravity_refused(m, g):
    with pytest.raises(ValueError):
        pressure_column_metrics([1, 2, 3], m, g)


def test_native_pressure_units_and_row_identity(tmp_path):
    path = tmp_path/'fixture.mod'
    text = '0.0 No mass for plane-parallel models\nk lgTauR lgTau5 Depth T Pe Pg Prad Pturb\n1 -1 -1 0 5000 10 100 20 0\n2 0 0 1 6000 20 200 30 0\n'
    path.write_text(text)
    s = dict(ndep=2, logtau_ross=np.array([-1., 0.]), temperature_K=np.array([5000., 6000.]), electron_pressure_Pa=np.array([1., 2.]))
    p = read_pressure_table(path, s)
    np.testing.assert_array_equal(p['gas_pressure_Pa'], [10., 20.])
    np.testing.assert_array_equal(p['radiation_pressure_Pa'], [2., 3.])
    with pytest.raises(ValueError, match='identities'):
        read_pressure_table(path, dict(s, temperature_K=np.array([5000., 6100.])))
