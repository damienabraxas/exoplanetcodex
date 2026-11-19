"""RYA-1282 Bride G5 speed test: fixed log tau500 grid and warm start, on one column.

Cases (all loglin collision extension, fixed A(O)):
  grid_lte   - target column resampled to --npts points on log tau500 [--top, --bottom], LTE start
  donor_lte  - neighbour column on the same grid, LTE start (supplies departures)
  grid_warm  - target column on the grid, started from the donor's departure coefficients
The native-grid solve of the same column (run_bride_stagger_columns) is the accuracy reference.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from dataclasses import replace
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--column', required=True)
    p.add_argument('--donor', required=True)
    p.add_argument('--npts', type=int, default=80)
    p.add_argument('--top', type=float, default=-6.0)
    p.add_argument('--bottom', type=float, default=2.0)
    p.add_argument('--maxiter', type=int, default=6000)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    from config.constants import require_codex_path
    from multi3d_read import read_atmos_multi3d
    from pipeline.bride_columns import (continuum_logtau, extend_collisions, mu1_spectrum,
                                        resample_logtau, stagger_column)
    from pipeline.bride_parameter_checks import marcs_filename
    from pipeline.bride_reference import read_marcs_structure
    from pipeline.multi_atom import read_multi_atom
    from pipeline.wavelength_util import air_to_vac
    sys.path.insert(0, str(ROOT / 'scripts'))
    from run_bride_stagger_columns import SNAPSHOT, SNAPSHOT_LOGG, STEP_A, WINDOW_AIR_A

    base = Path(require_codex_path('engines.bride_eval'))
    cube = read_atmos_multi3d(str(base / f'mesh.{SNAPSHOT}'), str(base / f'atm3d.{SNAPSHOT}'))
    atom = replace(read_multi_atom(require_codex_path('experiments.bride_oxygen_atom')), abundance=8.70)
    ab = read_marcs_structure(require_codex_path('grids.marcs_standard') / marcs_filename(5750, 4.5))['abundances_dex']
    wave = np.arange(WINDOW_AIR_A[0], WINDOW_AIR_A[1] + STEP_A / 2, STEP_A)
    vac = air_to_vac(wave) * .1
    grid = np.linspace(a.top, a.bottom, a.npts)
    a.output.mkdir(parents=True, exist_ok=False)
    results = dict(grid=dict(npts=a.npts, top=a.top, bottom=a.bottom), cases=[])

    def run(name, xy, departures=None):
        col = stagger_column(cube, *xy)
        lt = continuum_logtau(col, atom, abundances_dex=ab, logg_cgs=SNAPSHOT_LOGG)
        rc = resample_logtau(col, lt, grid)
        use = extend_collisions(atom, rc['temperature_K'].max() * 1.0001, 'loglin') \
            if rc['temperature_K'].max() > 1e4 else atom
        t0 = time.perf_counter()
        r = mu1_spectrum(rc, use, abundances_dex=ab, logg_cgs=SNAPSHOT_LOGG, wave_vac_nm=vac,
                         maxiter=a.maxiter, seconds=3600., departures=departures)
        np.savez_compressed(a.output / f'{name}.npz', wavelength_air_A=wave, line=r['line'],
                            continuum=r['continuum'], populations=r['populations'], nstar=r['nstar'])
        rec = dict(name=name, column=list(xy), iterations=r['iterations'], seconds=time.perf_counter() - t0,
                   seconds_per_iteration=r['seconds_solve'] / r['iterations'], final=r['history'][-1])
        results['cases'].append(rec)
        (a.output / 'results.json').write_text(json.dumps(results, indent=1) + '\n')
        print(name, rec['iterations'], round(rec['seconds'], 1), flush=True)
        return r

    xy = tuple(int(v) for v in a.column.split(','))
    dxy = tuple(int(v) for v in a.donor.split(','))
    run('grid_lte', xy)
    donor = run('donor_lte', dxy)
    run('grid_warm', xy, departures=donor['populations'] / donor['nstar'])


if __name__ == '__main__':
    main()
