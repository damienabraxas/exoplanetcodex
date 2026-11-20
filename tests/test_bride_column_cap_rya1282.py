"""RYA-1282: the Bride column path's default top boundary (ColumnCap).

A uniform-source slab of length L above the top layer replaces zero incoming radiation.
L = native overlying column mass / top density (MARCS) or the top layer's pressure scale
height P/(rho g) (STAGGER). Tests: the slab's analytic limits, per-direction coefficients
for moving columns, the refresh/convergence contract, and L's two definitions agreeing.
"""
from types import SimpleNamespace as S
import numpy as np
import pytest

pytest.importorskip('lightweaver')
from pipeline.bride_top_cap import ColumnCap, cap_intensity, scale_height_cap_length  # noqa: E402


def _ctx(chi, eta, sca, J_top, mus, wave):
    """chi/eta: (wave, mu) for the DOWNWARD (toObs=0) direction; upward filled with junk."""
    nw, nm = chi.shape
    full_chi = np.zeros((nw, nm, 2, 3)); full_eta = np.zeros_like(full_chi)
    full_chi[:, :, 0, 0] = chi; full_eta[:, :, 0, 0] = eta
    full_chi[:, :, 1, 0] = 99.; full_eta[:, :, 1, 0] = 99.           # must never be read
    return S(spect=S(wavelength=wave, J=np.tile(J_top[:, None], (1, 3))),
             depthData=S(chi=full_chi, eta=full_eta),
             background=S(sca=np.tile(sca[:, None], (1, 3))),
             atmos=S(pyAtmos=S(muz=mus)))


def test_zero_length_cap_is_vacuum_and_unrefreshed_cap_is_zero():
    wave = np.array([500., 600.]); mus = np.array([0.3, 1.0])
    cap = ColumnCap(0.0); cap.set_required_angles(None, None, mus, None)
    assert np.array_equal(cap.compute_bc(None, S(wavelength=wave)), np.zeros((2, 2, 1)))
    cap.update_from_context(_ctx(np.ones((2, 2)), np.ones((2, 2)), np.zeros(2), np.zeros(2), mus, wave), np.zeros(2))
    assert np.all(cap.compute_bc(None, S(wavelength=wave)) == 0)


def test_optically_thick_cap_returns_the_source_function_and_thin_cap_the_column_emission():
    wave = np.array([500., 600.]); mus = np.array([0.5, 1.0])
    chi = np.full((2, 2), 1e-3); eta = np.full((2, 2), 2e-3)                 # S = 2
    thick = ColumnCap(1e7); thick.set_required_angles(None, None, mus, None)
    thick.update_from_context(_ctx(chi, eta, np.zeros(2), np.zeros(2), mus, wave), np.zeros(2))
    np.testing.assert_allclose(thick.compute_bc(None, S(wavelength=wave)), 2.0, rtol=1e-12)
    thin = ColumnCap(1.0); thin.set_required_angles(None, None, mus, None)
    thin.update_from_context(_ctx(chi, eta, np.zeros(2), np.zeros(2), mus, wave), np.zeros(2))
    got = thin.compute_bc(None, S(wavelength=wave))[:, :, 0]
    exact = 2.0 * (1 - np.exp(-1e-3 * 1.0 / mus))[None, :]                     # S (1 - e^-tau)
    np.testing.assert_allclose(got, np.broadcast_to(exact, got.shape), rtol=1e-12)
    np.testing.assert_allclose(got, np.broadcast_to(2e-3 / mus[None, :], got.shape), rtol=2e-3)  # thin: eta L / mu


def test_each_incoming_ray_uses_its_own_downward_coefficients():
    wave = np.array([500., 600.]); mus = np.array([0.2, 1.0])
    chi = np.array([[1e-3, 5e-3], [1e-3, 5e-3]]); eta = chi * np.array([1., 3.])  # S differs per ray
    cap = ColumnCap(1e7); cap.set_required_angles(None, None, mus, None)
    cap.update_from_context(_ctx(chi, eta, np.zeros(2), np.zeros(2), mus, wave), np.zeros(2))
    bc = cap.compute_bc(None, S(wavelength=wave))[:, :, 0]
    np.testing.assert_allclose(bc[:, 0], 1.0, rtol=1e-12)
    np.testing.assert_allclose(bc[:, 1], 3.0, rtol=1e-12)                      # upward 99s never used


def test_scattering_emissivity_uses_the_previous_top_J():
    wave = np.array([500.]); mus = np.array([1.0])
    cap = ColumnCap(1e7); cap.set_required_angles(None, None, mus, None)
    cap.update_from_context(_ctx(np.array([[1e-3]]), np.array([[0.]]), np.array([5e-4]), np.array([4.]), mus, wave),
                            np.array([4.]))
    np.testing.assert_allclose(cap.compute_bc(None, S(wavelength=wave)), 4. * 5e-4 / 1e-3, rtol=1e-12)


def test_refresh_reports_the_relative_boundary_change_and_converges_to_zero():
    wave = np.array([500., 600.]); mus = np.array([1.0])
    c = _ctx(np.full((2, 1), 1e-3), np.full((2, 1), 1e-3), np.zeros(2), np.zeros(2), mus, wave)
    cap = ColumnCap(1e5); cap.set_required_angles(None, None, mus, None)
    assert cap.update_from_context(c, np.zeros(2)) == pytest.approx(1.0)
    assert cap.update_from_context(c, np.zeros(2)) == 0.0


def test_cap_refuses_wavelength_extrapolation_and_invalid_length():
    with pytest.raises(ValueError):
        ColumnCap(-1.)
    wave = np.array([500., 600.]); mus = np.array([1.0])
    cap = ColumnCap(1.); cap.set_required_angles(None, None, mus, None)
    cap.update_from_context(_ctx(np.ones((2, 1)), np.ones((2, 1)), np.zeros(2), np.zeros(2), mus, wave), np.zeros(2))
    with pytest.raises(ValueError):
        cap.compute_bc(None, S(wavelength=np.array([400., 500.])))


def test_scale_height_matches_the_native_marcs_overlying_column():
    """Measured on the 14 Gerber+2023 benchmark models: max |Hp/(m0/rho0) - 1| = 0.0201
    (4500 K, where the rule neglects molecules). Bound 0.025 = that measurement, rounded up."""
    from config.constants import require_codex_path
    from pipeline.bride_parameter_checks import marcs_filename
    from pipeline.bride_reference import read_marcs_structure
    from pipeline.bride_columns import marcs_column
    worst = 0.
    for teff in (4500, 5750, 6500):
        path = require_codex_path('grids.marcs_standard') / marcs_filename(teff, 4.5)
        if not path.exists():
            pytest.skip('MARCS grid not on this host')
        s = read_marcs_structure(path); c = marcs_column(s)
        hp = scale_height_cap_length(c['temperature_K'][0], c['ne_m3'][0], c['nHtot_m3'][0],
                                     s['abundances_dex'], s['logg_cgs'])
        worst = max(worst, abs(hp / c['cap_length_m'] - 1))
    assert worst < 0.025


def test_estimated_error_is_exact_on_a_geometric_sequence_and_refuses_non_contraction():
    from pipeline.bride_columns import estimated_error
    rho = 0.97
    d = [rho ** k for k in range(40)]
    err, r = estimated_error(d)
    assert r == pytest.approx(rho, rel=1e-12)
    assert err == pytest.approx(sum(rho ** k for k in range(40, 4000)), rel=1e-9)   # the true remaining distance
    assert estimated_error([1.0] * 8) == (float('inf'), 1.0)                          # stalled: never "converged"
    assert estimated_error([1.0, 2.0])[0] == float('inf')                             # too short to know
    # a small step on a slow contraction is NOT convergence: d = 1e-6 at rho 0.97 -> 3.2e-5 away
    slow = [1e-6 / rho ** k for k in range(10)][::-1]
    assert estimated_error(slow)[0] == pytest.approx(1e-6 * rho / (1 - rho), rel=1e-9)
