import numpy as np
import pytest
from pipeline.bride_ts_points import point_input, point_driver_source


def test_point_units_order_and_native_real_identity():
    t = np.array([4045.9, 6000.1])
    pe = np.array([.01234, 10.01])
    text = point_input(t, pe)
    assert text.splitlines()[0] == '2'
    parsed = np.loadtxt(text.splitlines()[1:], dtype=np.float32)
    np.testing.assert_array_equal(parsed[:, 0], t.astype(np.float32))
    np.testing.assert_array_equal(parsed[:, 1], (pe/.1).astype(np.float32))


@pytest.mark.parametrize('t,pe', [([1, 2], [1]), ([1, np.nan], [1, 2]), ([0, 2], [1, 2]), ([1, 2], [1, np.inf]), ([1, 2], [1, 1e300])])
def test_bad_points_refused(t, pe):
    with pytest.raises(ValueError):
        point_input(t, pe)


def test_source_insertion_is_bounded_and_unknown_layout_refused():
    source = '      PROGRAM BABSMA\nC placeholder\n      NEWT=2\nC opacity calls\n'
    changed = point_driver_source(source)
    assert changed.endswith('      NEWT=2\nC opacity calls\n')
    assert changed.index('READ(97') < changed.index('      NEWT=2')
    assert 'BRIDE_N.NE.NTAU' in changed
    with pytest.raises(ValueError):
        point_driver_source(changed)
