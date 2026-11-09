import numpy as np
import pytest


def test_unit_geometry_preserves_sign_and_distinguishes_rank_loss():
    from scripts.analyze_bride_profile_geometry import direction_geometry
    matrix = np.array([[2., 0., -5.], [0., 3., -5.], [0., 0., 0.]])
    unit, result = direction_geometry(matrix)
    np.testing.assert_allclose(unit.T@unit, [[1, 0, -1/np.sqrt(2)], [0, 1, -1/np.sqrt(2)], [-1/np.sqrt(2), -1/np.sqrt(2), 1]], atol=1e-15)
    assert result['smallest_to_largest_singular_ratio'] == 0
    _, independent = direction_geometry(np.diag([2., 3., 5.]))
    assert independent['smallest_to_largest_singular_ratio'] == 1
    with pytest.raises(ValueError, match='Zero'): direction_geometry(np.zeros((4, 3)))
    with pytest.raises(ValueError, match='finite'): direction_geometry(np.full((4, 3), np.nan))
