"""RYA-1282: O departure grids on the Gerber nodes (pipeline.o_departure_grid)."""
import glob
import json
import os
from pathlib import Path

import numpy as np
import pytest


def test_aux_parse_and_a_axis_thinning(tmp_path):
    from pipeline.o_departure_grid import gerber_nodes, select_nodes
    rows = ["#atmos ID, Teff, logg, [Fe/H], [a/Fe], mass, vt, A(X), pointer"]
    for k, A in enumerate(np.round(np.arange(7.6, 9.61, 0.1) - 2.0, 2)):
        rows.append(f"'p5750_g+4.5_m0.0_t01_st_z-2.00_a+0.40'\t5750.0\t4.5\t-2.0\t0.4\t0.00\t1.00\t{A:.4f}\t{1001 + k}")
    rows.append("'p4000_g+1.0_m1.0_t02_st_z+0.00_a+0.00'\t4000.0\t1.0\t0.0\t0.0\t1.00\t2.00\t8.7000\t9")
    (tmp_path / 'aux').write_text('\n'.join(rows) + '\n')
    n = gerber_nodes(tmp_path / 'aux')
    assert len(n) == 2 and len(n['p5750_g+4.5_m0.0_t01_st_z-2.00_a+0.40']['abundances']) == 21
    s = select_nodes(n, teff=(4500, 7000), logg_min=3.0, a_step=0.2)
    assert list(s) == ['p5750_g+4.5_m0.0_t01_st_z-2.00_a+0.40']
    ab = np.array(s['p5750_g+4.5_m0.0_t01_st_z-2.00_a+0.40']['abundances'])
    assert len(ab) == 11 and ab[0] == pytest.approx(5.6) and ab[-1] == pytest.approx(7.6)
    assert np.allclose(np.diff(ab), 0.2)


def test_ts_writer_round_trips_in_the_gerber_layout(tmp_path):
    from pipeline.o_departure_grid import assemble, read_record, HEADER_BYTES, ID_BYTES
    lt = np.linspace(-5, 2, 56)
    b = np.random.default_rng(1).uniform(0.5, 1.5, (2, 41, 56))
    np.savez(tmp_path / 'm.npz', model_id='p5750_g+4.5_m0.0_t01_st_z+0.00_a+0.00', engine='x', teff=5750., logg=4.5,
             feh=0., alpha=0., vturb=1., abundances=np.array([8.6, 8.8]), logtau=lt, b=b,
             ok=np.array([True, False]), steps=np.array([3, -1]), seconds=np.array([1., 1.]), error=np.array(['', 'x']))
    c = assemble(tmp_path, tmp_path / 'g.bin', tmp_path / 'g.aux', header_text='test grid')
    assert c == dict(records=1, failed_nodes=1, models=1)
    aux = [l for l in (tmp_path / 'g.aux').read_text().splitlines() if not l.startswith('#')]
    ptr = int(aux[0].split()[-1])
    assert ptr == HEADER_BYTES + 1
    mid, ltr, dep = read_record(tmp_path / 'g.bin', ptr)
    assert mid.startswith('p5750') and np.allclose(ltr, lt) and np.array_equal(dep, b[0])
    assert os.path.getsize(tmp_path / 'g.bin') == HEADER_BYTES + ID_BYTES + 8 + 8 * 56 + 8 * 41 * 56   # Gerber O stride 19324


def test_multi23_in_gerbers_setup_reproduces_the_gerber_record():
    """MULTI 2.3 with a solar ABUND and unscaled background lines = the Gerber grid's setup."""
    from config.constants import codex_path
    home = codex_path('engines.mul23')
    if not (home / 'source_dp' / 'mul23lus.x').exists():
        pytest.skip('MULTI 2.3 build not held here')
    from pipeline.o_departure_grid import solve_multi23, _load_model
    rec = sorted(glob.glob(os.path.expanduser('~/codex/bride_inputs/rya1282/bench_gerber2023_O_MARCS_zladder/p5750_g+4.5_-2.00_*')))
    if not rec:
        pytest.skip('zladder record not held here')
    r = json.loads(Path(rec[0]).read_text())
    S = _load_model(r['atmosphere_file']['name'].removesuffix('.mod')); S['feh'] = -2.0
    bref = np.asarray(r['reference']['departures'], float); bref = bref if bref.shape[0] == 41 else bref.T
    lt = np.log10(np.asarray(r['reference']['tau'], float)); phot = (lt >= -3) & (lt <= 0)
    TL = [5, 7, 8, 9]
    b, _ = solve_multi23(S, float(r['reference']['abundance']), solar_abund=True, scale_lines=False)
    assert np.abs(np.log10(b[TL][:, phot] / bref[TL][:, phot])).max() < 0.004      # measured 0.0012
    bc, _ = solve_multi23(S, float(r['reference']['abundance']))                     # artefacts removed
    assert np.abs(np.log10(bc[TL][:, phot] / bref[TL][:, phot])).max() > 0.02       # measured 0.037
