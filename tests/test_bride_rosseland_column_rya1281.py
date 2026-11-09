import numpy as np
import pytest
from pipeline.bride_rosseland_column import column_input, read_rosseland_opacity, log_opacity_integral


def test_log_cgs_opacity_and_tau_input_units():
    text = column_input([-4., -3., -2., -1.], [.1, 1., 10., 100.])
    assert text.splitlines()[0] == '4'
    values = np.loadtxt(text.splitlines()[1:])
    np.testing.assert_array_equal(values[:, 0], [-4, -3, -2, -1])
    np.testing.assert_array_equal(values[:, 1], [0, 1, 2, 3])


@pytest.mark.parametrize('tau,k', [([-3., -2., -1.], [1.]*3), ([-4., -3., -2., -1.], [1., 0., 1., 1.]),
    ([-4., -3., -2., np.nan], [1.]*4), ([-4., -3., -2., 100.], [1.]*4),
    ([-4., -3., -2., -2.+1e-10], [1.]*4)])
def test_invalid_or_native_real_collapsed_inputs_refused(tau, k):
    with pytest.raises(ValueError):
        column_input(tau, k)


def test_opacity_table_units_and_row_identity(tmp_path):
    p = tmp_path/'table.mod'
    p.write_text('k lgTauR KappaRoss Density Mu Vconv Fconv RHOX\n1 -1 2 1 1 0 0 1\n2 0 4 1 1 0 0 2\n')
    s = dict(ndep=2, logtau_ross=np.array([-1., 0.]))
    np.testing.assert_array_equal(read_rosseland_opacity(p, s), [.2, .4])
    with pytest.raises(ValueError):
        read_rosseland_opacity(p, dict(s, logtau_ross=np.array([-2., 0.])))


@pytest.mark.parametrize('power', [0., 1., 2.])
def test_declared_powerlaw_opacity_integral_against_closed_form(power):
    tau = np.array([.01, .02, .08, .3, 1.]); top = 7.
    actual = log_opacity_integral(np.log10(tau), tau**power, top)
    if power == 0:
        expected = top+tau-tau[0]
    elif power == 1:
        expected = top+np.log(tau/tau[0])
    else:
        expected = top+1/tau[0]-1/tau
    np.testing.assert_allclose(actual, expected, rtol=2e-15, atol=1e-14)


def test_actual_cap_uses_top_column_density_not_deeper_mass_profile():
    from pipeline.bride_top_cap import cap_from_structure
    first = dict(column_mass_kg_m2=np.array([1., 2., 3.]), rho_kg_m3=np.array([.1, .2, .3]))
    other = dict(first, column_mass_kg_m2=np.array([1., 4., 20.]))
    assert cap_from_structure(first).length_m == cap_from_structure(other).length_m == 10.
