import numpy as np
import pytest

from pipeline.bride_rectangular_checks import rectangular_weights, rectangular_structure


def test_bilinear_exact_cross_term_and_input_permutation():
    points = np.array([[0., 0.], [0., 2.], [4., 0.], [4., 2.]])
    field = lambda x, y: 3 + 2*x - y + 5*x*y
    for order in ([0, 1, 2, 3], [3, 1, 0, 2]):
        p = points[order]
        w = rectangular_weights(p, [1., 1.5])
        assert np.all(w >= 0)
        assert w.sum() == pytest.approx(1)
        assert w @ field(p[:, 0], p[:, 1]) == pytest.approx(field(1., 1.5))
    for index, node in enumerate(points):
        np.testing.assert_array_equal(rectangular_weights(points, node), np.eye(4)[index])


@pytest.mark.parametrize('points,target', [
    ([[0, 0], [0, 1], [1, 0], [1, 0]], [0, 0]),
    ([[0, 0], [0, 1], [1, 0], [1, 1]], [2, 0]),
    ([[0, 0], [0, 1], [1, 0], [1, 1]], [np.nan, 0]),
])
def test_invalid_cell_or_extrapolation_refused(points, target):
    with pytest.raises(ValueError):
        rectangular_weights(points, target)


def test_log_layer_structure_has_analytic_geometry_and_no_column():
    nodes = []
    for x in (0., 1.):
        for y in (0., 1.):
            # Analytic positive fields whose logarithms are bilinear.
            scale = np.exp(x + 2*y + x*y)
            nodes.append(((x, y), dict(ndep=3, logtau_ross=np.array([-1., 0., 1.]),
                logtau_500=np.array([-2., -1., 0.]), abundances_dex={1:12.}, hydrogen_mass_fraction=.7,
                height_m=np.array([2., 0., -3.])*scale, temperature_K=np.full(3, 5000.+x*100+y*200),
                electron_pressure_Pa=np.full(3, scale), rho_kg_m3=np.full(3, 2*scale),
                column_mass_kg_m2=np.array([1., 2., 3.]))))
    s, receipt = rectangular_structure(nodes, [.5, .5], thermodynamic_scale='log_positive_layers')
    scale = np.exp(1.75)
    np.testing.assert_allclose(s['height_m'], np.array([2., 0., -3.])*scale)
    np.testing.assert_allclose(s['rho_kg_m3'], 2*scale)
    np.testing.assert_allclose(s['electron_pressure_Pa'], scale)
    np.testing.assert_allclose(s['temperature_K'], 5150.)
    assert 'column_mass_kg_m2' not in s
    assert receipt['weights'] == [.25]*4
    nodes[0][1]['hydrogen_mass_fraction'] = .6
    with pytest.raises(ValueError, match='composition'):
        rectangular_structure(nodes, [.5, .5], thermodynamic_scale='linear')
