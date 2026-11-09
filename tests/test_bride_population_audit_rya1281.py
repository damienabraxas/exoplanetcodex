import numpy as np
import pytest
from pipeline.multi_atom import MultiCollision
from pipeline.bride_parameter_checks import diagnostic_collision_scaling
from scripts.audit_bride_oxygen_populations import conserve_departure_replay


def test_external_replay_preserves_totals_and_level_ratios():
    lte = np.array([[100., 30.], [2., 4.], [1., 1.]])
    b = np.array([[1.1, .9], [.7, .8], [.4, .6]])
    pop, factors = conserve_departure_replay(b, lte)
    np.testing.assert_allclose(pop.sum(axis=0), lte.sum(axis=0), rtol=4*np.finfo(float).eps)
    np.testing.assert_allclose((pop/lte)[1]/(pop/lte)[0], b[1]/b[0])
    assert not np.allclose(factors, 1)
    with pytest.raises(ValueError):
        conserve_departure_replay(b[:, :1], lte)
    with pytest.raises(ValueError):
        conserve_departure_replay(b*0, lte)


def test_rate_control_scales_only_named_class_without_mutating_source():
    from types import SimpleNamespace
    source = SimpleNamespace(collisions=[MultiCollision('CE', 2, 1, np.array([2., 3.])),
                                        MultiCollision('CH', 2, 1, np.array([5., 7.]))])
    changed = diagnostic_collision_scaling(source, {'CH': .1})
    np.testing.assert_array_equal(source.collisions[1].rates, [5., 7.])
    np.testing.assert_array_equal(changed.collisions[0].rates, [2., 3.])
    np.testing.assert_allclose(changed.collisions[1].rates, [.5, .7])
    for scales in ({'CI': 2}, {'CE': -1}, {'CE': float('nan')}):
        with pytest.raises(ValueError):
            diagnostic_collision_scaling(source, scales)
