"""RYA-1282: the n(H)-from-(T, Pe) diagnostic behind the Gerber-grid H-collision finding."""
import numpy as np
import pytest


def _marcs(z, a):
    from config.constants import require_codex_path
    from pipeline.bride_reference import read_marcs_structure
    return read_marcs_structure(require_codex_path('grids.marcs_standard')
                                / f'p5750_g+4.5_m0.0_t01_st_z{z}_a{a}_c+0.00_n+0.00_o{a}_r+0.00_s+0.00.mod')


def test_eos_reproduces_marcs_hydrogen_density_and_the_solar_ratio_is_metallicity_ordered():
    from pipeline.bride_reference import marcs_number_densities
    from pipeline.bride_reference_eos import nh_from_pe, solar_mix_nh_ratio
    try:
        sol = _marcs('+0.00', '+0.00')
        models = {z: _marcs(z, a) for z, a in (('+1.00', '+0.00'), ('-1.00', '+0.40'), ('-3.00', '+0.40'))}
    except Exception as e:                                       # noqa: BLE001
        pytest.skip(f'MARCS / TS data not held here: {e}')
    ratios = {}
    for z, S in [('+0.00', sol)] + list(models.items()):
        ne, nh = marcs_number_densities(S)
        T = np.asarray(S['temperature_K'])
        got = nh_from_pe(T, ne * 1e-6, S['abundances_dex'])
        assert np.max(np.abs(got / (nh * 1e-6) - 1)) < 0.02         # measured 1.1-1.5%
        lt = np.asarray(S['logtau_500']); k = int(np.argmin(abs(lt + 1)))
        ratios[z] = solar_mix_nh_ratio(T, ne * 1e-6, S['abundances_dex'], sol['abundances_dex'])[k]
    assert ratios['+0.00'] == pytest.approx(1.0, abs=1e-12)
    assert ratios['+1.00'] > 5 and 0.15 < ratios['-1.00'] < 0.35 and ratios['-3.00'] < 0.12   # 9.5, 0.26, 0.08
