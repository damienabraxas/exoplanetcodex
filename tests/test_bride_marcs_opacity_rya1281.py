"""Source-format dimensions, interleaving and opacity units, independent of RT."""
import gzip

import numpy as np
import pytest

from pipeline.bride_marcs_opacity import read_marcs_opacity


def example():
    return (' MRXF 2 5000.0\n3\n1000 5000 10000\n'
        '1 0.1 4000 2 3 4e-7 0 5\n1 2 3 4 5 6\n'
        '1 0.2 5000 3 4 5e-7 1 10\n6 5 4 3 2 1\n'
        + ' '.join(['12'] + ['-1']*91))


def test_relative_coefficients_interleaving_and_si_units(tmp_path):
    path = tmp_path/'opacity.gz'
    path.write_bytes(gzip.compress(example().encode()))
    result = read_marcs_opacity(path)
    np.testing.assert_array_equal(result['absorption_m2_kg'], [[.5, 6], [1.5, 4], [2.5, 2]])
    np.testing.assert_array_equal(result['scattering_m2_kg'], [[1, 5], [2, 3], [3, 1]])
    np.testing.assert_allclose(result['rho_kg_m3'], [4e-4, 5e-4])
    np.testing.assert_allclose(result['electron_pressure_Pa'], [.2, .3])
    np.testing.assert_array_equal(result['vturb_m_s'], [0, 1000])
    assert result['wavelength_medium'] == 'UNVERIFIED'


@pytest.mark.parametrize('transform', [
    lambda x: x.replace('MRXF', 'MXXX'),
    lambda x: x+' 1', lambda x: ' '.join(x.split()[:-1]),
    lambda x: x.replace('1000 5000 10000', '1000 1000 10000'),
    lambda x: x.replace('4000', 'nan'),
    lambda x: x.replace('1 2 3 4 5 6', '-1 2 3 4 5 6'),
])
def test_bad_source_is_rejected(tmp_path, transform):
    path = tmp_path/'opacity'
    path.write_text(transform(example()))
    with pytest.raises(ValueError):
        read_marcs_opacity(path)
