"""RYA-1282: Bride passive metal bound-free from TS/MARCS's own OP tables.

Measured on MARCS 5750/4.5 against native TS (iSpec build, RYA-1281 runs):
  * Fe I kappa per neutral atom, port vs native TS packet, 56 depths x 1601 lambda:
    flat in lambda to 1.33e-4 (residual only at the 2-A node pairs: vac->air rounding)
  * Irwin U(Fe I) vs U implied by native TS: 2.6e-6
  * total background absorption vs native TS, 200-360 nm x 56 depths, |log10 ratio|:
    Lightweaver passive atoms median 0.226 / max 1.278 dex; with TS metals 0.025 / 0.098 dex
    (remaining: CH continuum, absent from the Bride, and the Bride-vs-TS EOS).
Bounds are those measurements rounded up; the negative control must FAIL the parity bound.
"""
from pathlib import Path
import numpy as np
import pytest

pytest.importorskip('lightweaver')

PACKET = Path.home() / 'codex/bride_inputs/rya1281/fe-nuv-packet-v1/packet.npz'
TS_RUN = Path.home() / 'codex/bride_mac_checks/rya1281-ts-nuv-v1/NUV_source_air_baseline/continuum.opa'


def _data():
    try:
        from pipeline.bride_ts_background import ts_data_dir
        d = ts_data_dir()
        if not (d / 'jonabs_vac_v19.2.dat').exists():
            raise FileNotFoundError(d)
        return d
    except Exception as e:                                         # noqa: BLE001
        pytest.skip(f'TS DATA not held here: {e}')


def test_fe_block_is_read_as_ts_reads_it():
    from pipeline.ts_jonabs import JonabsComponent, read_component
    d = _data()
    b = read_component(d / 'jonabs_vac_v19.2.dat', 'FeI ')
    assert b['xkap'].shape == (38, 18)
    assert b['lflags'] == dict(ILOGL=0, KVADL=0, MINEX=0, MAXEX=0, NLATB=38)
    assert b['tflags'] == dict(ILOGT=1, KVADT=1, MINET=1, MAXET=0, NTETB=18, ITETA=0)
    fe = JonabsComponent(d / 'jonabs_vac_v19.2.dat', 'FeI ')
    assert (fe.tbot, fe.delt, fe.ntm) == ([1000., 8000.], [700., 3000.], [10, 7])


def test_table_nodes_are_reproduced_and_outside_is_zero():
    from pipeline.ts_jonabs import JonabsComponent
    from pipeline.wavelength_util import air_to_vac
    d = _data()
    fe = JonabsComponent(d / 'jonabs_vac_v19.2.dat', 'FeI ')
    k = 24                                                     # an interior lambda node
    w_vac_nm = air_to_vac(np.array([fe.lam_air_A[k]])) / 10.
    for j in (2, 5, 12):                                       # T nodes in both segments
        got = fe(w_vac_nm, fe.tet[j])[0, 0]
        want = fe.xkap[k, j]
        assert got == pytest.approx(want, rel=1e-5)            # vac<->air round trip only
    assert np.all(fe(np.array([50., 9000.]), 5000.) == 0)      # MINEX = MAXEX = 0
    assert np.all(fe(np.array([300.]), 35000.) == 0)           # MAXET = 0: no hot extrapolation


def test_fe_port_matches_native_ts_per_atom_opacity_and_partition_function():
    if not PACKET.exists():
        pytest.skip('native TS Fe packet not held here')
    from pipeline.ts_jonabs import JonabsComponent, irwin_partition
    from pipeline.wavelength_util import vac_to_air
    d = _data(); p = np.load(PACKET)
    w = p['wavelength_vac_nm']; T = p['native_temperature_K']; lam = vac_to_air(w * 10)
    stim = 1 - np.exp(-28556. * (5040. / T)[None, :] / lam[:, None])     # TS detabs factor
    tab = JonabsComponent(d / 'jonabs_vac_v19.2.dat', 'FeI ')(w, T) * 1e-4 * stim
    r = p['opacity_net_per_neutral_Fe_m2'] / tab                         # = 1/U(T)
    assert np.max(r.max(0) / r.min(0) - 1) < 2e-4
    U = irwin_partition(d / 'IRWIN_atoms_v07.3.dat', 26, 0, T)
    np.testing.assert_allclose(U, 1 / np.median(r, 0), rtol=1e-5)


def _solar_context():
    from dataclasses import replace
    from config.constants import require_codex_path
    from pipeline.multi_atom import read_multi_atom
    from pipeline.bride_parameter_checks import marcs_filename
    from pipeline.bride_reference import read_marcs_structure
    from pipeline.bride_columns import column_context, marcs_column
    try:
        atom = replace(read_multi_atom(require_codex_path('experiments.bride_oxygen_atom')), abundance=8.70)
        S = read_marcs_structure(require_codex_path('grids.marcs_standard') / marcs_filename(5750, 4.5))
    except Exception as e:                                        # noqa: BLE001
        pytest.skip(f'O atom / MARCS not held here: {e}')
    ctx = column_context(marcs_column(S), atom, lte=True, abundances_dex=S['abundances_dex'],
                         logg_cgs=S['logg_cgs'], vturb_m_s=S['vturb_m_s'], rays=3, top='vacuum')
    return ctx


def _background(ctx, factory, wave_nm):
    eq = ctx.kwargs['eqPops']; radset = ctx.kwargs['spect'].radSet
    p = factory(eq, radset, wave_nm)
    sh = (len(wave_nm), ctx.atmos.pyAtmos.Nspace)
    chi, eta, sca = np.zeros(sh), np.zeros(sh), np.zeros(sh)
    p.compute_background(ctx.atmos, chi, eta, sca)
    return chi, eta, sca, p


def test_replacement_is_exact_lte_and_never_double_counts():
    from lightweaver import constants as c
    from lightweaver.LwCompiled import FastBackground
    from pipeline.bride_ts_background import TsMetalBackground
    _data(); ctx = _solar_context()
    w = np.linspace(200., 420., 157)
    base = _background(ctx, lambda e, r, x: FastBackground(e, r, x, Nthreads=1), w)
    none = _background(ctx, lambda e, r, x: TsMetalBackground(e, r, x, species=()), w)
    for a, b in zip(base[:3], none[:3]):
        np.testing.assert_array_equal(a, b)                      # no species -> Lightweaver exactly
    fe = _background(ctx, lambda e, r, x: TsMetalBackground(e, r, x, species=[('Fe', 0)]), w)
    T = np.asarray(ctx.atmos.pyAtmos.temperature)
    new, old = fe[3]._ts_and_old(T)
    np.testing.assert_allclose(fe[0] - base[0], new - old, rtol=1e-12, atol=1e-30)
    B = 2 * c.HPlanck * c.CLight / (w[:, None] * 1e-9) ** 3 / np.expm1(c.HC / (w[:, None] * 1e-9 * c.KBoltzmann * T[None, :]))
    np.testing.assert_allclose(fe[1] - base[1], (new - old) * B, rtol=1e-10, atol=1e-30)
    np.testing.assert_array_equal(fe[2], base[2])                 # scattering untouched


def _parity(ctx, factory):
    from pipeline.bride_ts_budget import read_ts_opacity
    from pipeline.wavelength_util import air_to_vac
    ts = read_ts_opacity(TS_RUN)
    T = ctx.atmos.pyAtmos.temperature
    assert np.max(abs(T / ts['temperature_K'] - 1)) < 1e-6
    w = air_to_vac(ts['wavelength_air_A']) / 10.
    chi, _, sca, _ = _background(ctx, factory, w)
    return np.abs(np.log10((chi - sca) / (ts['absorption_m2_kg'] * ts['rho_kg_m3'][None, :])))


def test_total_background_matches_native_ts_and_old_background_does_not():
    if not TS_RUN.exists():
        pytest.skip('native TS near-UV run not held here')
    from lightweaver.LwCompiled import FastBackground
    from pipeline.bride_ts_background import TsMetalBackground
    _data(); ctx = _solar_context()
    new = _parity(ctx, lambda e, r, x: TsMetalBackground(e, r, x))
    assert np.median(new) < 0.03 and new.max() < 0.11            # measured 0.025 / 0.098
    old = _parity(ctx, lambda e, r, x: FastBackground(e, r, x, Nthreads=1))
    assert np.median(old) > 0.15 and old.max() > 1.0              # measured 0.226 / 1.278
