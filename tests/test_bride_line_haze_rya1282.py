"""RYA-1282: MULTI-style background line opacity (pipeline.bride_line_haze)."""
import numpy as np
import pytest


def test_psemi_reproduces_multis_table_and_limits():
    from pipeline.bride_line_haze import psemi
    betref = np.array([-2.0, -1.699, -1.398, -1.0, -0.699, -0.398, 0.0, 0.301, 0.602, 1.0])
    pn = [1.160, 0.956, 0.758, 0.493, 0.331, 0.209, 0.100, 0.063, 0.040, 0.023]
    pc = [1.160, 0.977, 0.788, 0.554, 0.403, 0.290, 0.214, 0.201, 0.200, 0.200]
    inner = slice(1, -1)                     # the end nodes sit on MULTI's open limits
    np.testing.assert_allclose(psemi(1, 10 ** betref[inner]), np.array(pn)[inner], rtol=1e-9)
    np.testing.assert_allclose(psemi(2, 10 ** betref[inner]), np.array(pc)[inner], rtol=1e-9)
    assert psemi(1, np.array([0.005]))[0] == pytest.approx(0.27566 * (0.577 + np.log(0.005)))
    assert psemi(1, np.array([40.]))[0] == pytest.approx(0.066 / np.sqrt(40.))
    assert psemi(2, np.array([40.]))[0] == pytest.approx(0.200)


def test_epsilon_is_multis_van_regemorter_expression():
    from pipeline.bride_line_haze import multi_line_epsilon, psemi
    w, T, ne = 3000., 5000., 1e13
    beta = 6.62607015e-27 * 2.99792458e10 * 1e8 / w / 1.380649e-16 / T
    C = np.sqrt(T) / (20.6 * (w * 1e-8) ** 3) * (0.5 / psemi(2, np.array([beta]))[0]
                                                 + 0.5 / psemi(1, np.array([beta]))[0])
    eps = multi_line_epsilon([w], [T], [ne])
    assert eps.shape == (1, 1)
    assert eps[0, 0] == pytest.approx(ne / (C + ne), rel=1e-12)
    assert 0 < eps[0, 0] < 0.01                         # UV background lines: nearly pure scattering


def test_vald_long_format_parser_keeps_atoms_and_counts_molecules(tmp_path):
    from pipeline.bride_line_haze import read_vald_stellar
    p = tmp_path / 'v.txt'
    p.write_text(
        " 6910.00000, 9500.00000, 2, 9, 1.0 Wavelength region, lines selected, lines processed, Vmicro\n"
        "header\nheader\n"
        "'Mn 1',       6910.15261, -1.380,  5.2320,  0.5,  7.0257,  0.5,  3.320,  3.310,  3.310, 8.070,-4.190,-7.270, 0.002,\n"
        "'  LS                                                                      3d6.(5D).4p z6D*'\n"
        "'  LS                                                                       3d6.(5D).5s i6D'\n"
        "'_          Kurucz MnI 2007    1 wl:K07   1 gf:K07   1 K07   Mn'\n"
        "'CO 1',       6911.00000, -9.000,  1.0000,  1.0,  2.7937,  2.0,  0.000,  0.000,  0.000, 0.000, 0.000, 0.000, 0.001,\n"
        "'Fe 2',       6912.00000, -2.000,  3.0000,  2.5,  4.7937,  2.5,  1.000,  1.000,  1.000, 8.000,-6.000, 243.25, 0.010,\n")
    L = read_vald_stellar([p])
    assert list(L['Z']) == [25, 26] and list(L['ion']) == [0, 1]
    assert L['molecular_rows_skipped'] == 1
    assert L['loggf'][0] == -1.380 and L['rad'][0] == 8.070 and L['waals'][1] == 243.25


def test_one_isolated_line_integrates_to_its_strength():
    """The Voigt profile is normalised: integral of kappa d(nu) = S (pi e^2/mc gf n/U e^-E/kT stim)."""
    from pipeline.bride_line_haze import line_haze, _stage_fractions
    from pipeline.bride_ts_background import ts_data_dir, TsMetalBackground
    irw = ts_data_dir() / TsMetalBackground.IRWIN
    if not irw.exists():
        pytest.skip('TS data not held here')
    L = dict(Z=np.array([26]), ion=np.array([0]), wave_A=np.array([5000.]), loggf=np.array([-1.]),
             elow_eV=np.array([2.]), eup_eV=np.array([4.48]), rad=np.array([8.]), stark=np.array([-6.]),
             waals=np.array([-7.5]))
    T = np.array([5000.]); ne = np.array([1e13]); nH = np.array([1e16]); ab = np.full(92, -20.); ab[0] = 12.; ab[25] = 7.5
    w = np.linspace(4995., 5005., 40001)
    k = line_haze(w, T, ne, nH, nH, 1e5, ab, L, irw, window_A=5., strong_window_A=5.)[:, 0]
    nu = 2.99792458e10 / (w * 1e-8)
    integral = -np.trapezoid(k, nu)
    f, U = _stage_fractions(irw, 26, T, ne)
    n = 10 ** (7.5 - 12) * nH * f[0] / U[0]
    S = 0.026540081 * 0.1 * n * np.exp(-2. * 1.602176634e-12 / (1.380649e-16 * T)) * -np.expm1(
        -6.62607015e-27 * 2.99792458e10 / 5000e-8 / (1.380649e-16 * T))
    assert integral == pytest.approx(S[0], rel=0.03)    # the remainder is the far Lorentz wings


def test_provider_masks_the_active_atoms_lines_and_keeps_chi_above_sca():
    pytest.importorskip('lightweaver')
    import json
    from dataclasses import replace
    import pipeline.bride_columns as BC
    from config.constants import require_codex_path
    from pipeline.multi_atom import read_multi_atom
    from pipeline.bride_reference import read_marcs_structure
    from pipeline.bride_line_haze import read_vald_stellar, default_vald_paths
    paths = default_vald_paths()
    try:
        atom = replace(read_multi_atom(require_codex_path('experiments.bride_oxygen_atom')), abundance=8.7)
        S = read_marcs_structure(require_codex_path('grids.marcs_standard')
                                 / 'p5750_g+4.5_m0.0_t01_st_z+0.00_a+0.00_c+0.00_n+0.00_o+0.00_r+0.00_s+0.00.mod')
    except Exception as e:                                         # noqa: BLE001
        pytest.skip(f'inputs not held here: {e}')
    if not paths:
        pytest.skip('VALD lists not held here')
    lines = read_vald_stellar([p for p in paths if 'redopt' in p.name])
    kw = dict(lte=True, abundances_dex=S['abundances_dex'], logg_cgs=S['logg_cgs'], vturb_m_s=S['vturb_m_s'], rays=3)
    c0 = BC.column_context(BC.marcs_column(S), atom, **kw)
    c1 = BC.column_context(BC.marcs_column(S), atom, line_haze=dict(lines=lines, abundances_dex=S['abundances_dex']), **kw)
    chi0, chi1 = np.array(c0.background.chi), np.array(c1.background.chi)
    sca0, sca1 = np.array(c0.background.sca), np.array(c1.background.sca)
    wl = np.array(c1.spect.wavelength)
    added = chi1 - chi0
    triplet = (wl > 777.3) & (wl < 777.9)                          # 777 triplet (vacuum nm) line points
    assert triplet.sum() > 10
    assert np.all(added[triplet] == 0)                             # MULTI: no background lines in the atom's lines
    assert np.all(added >= -1e-12 * chi0)
    assert np.all(chi1 >= sca1)
    assert added.max() > 0                                         # the red-optical list does something
    assert np.all((sca1 - sca0) <= added * (1 + 1e-9) + 1e-300)    # scattering part <= the added opacity
