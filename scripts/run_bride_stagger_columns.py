"""RYA-1282 Bride G5: 1.5D O NLTE on STAGGER solar columns (Mac, serial, one thread).

Each listed column is solved as an independent plane-parallel atmosphere with the held
O41f atom (fixed A(O); no fitting) and the disk-centre (mu = 1) O I 777 spectrum is
written with its continuum, populations, convergence history and timings. Collision
tables above 10000 K are extended only by a named rule (--extension), run per rule.

  PYTHONPATH=.:~/bride_rya1118 $LW scripts/run_bride_stagger_columns.py \
      --columns 12,40 3,3 --extension clamp loglin --output <fresh dir>
"""
from __future__ import annotations

import argparse
import json
import platform
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

SNAPSHOT = 't5777g44m0005_00091_80x80x240'
SNAPSHOT_LOGG = 4.44   # STAGGER t5777g44 solar node
WINDOW_AIR_A = (7770.5, 7777.0)
STEP_A = 0.005


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--columns', nargs='*', default=[], help='ix,iy pairs')
    p.add_argument('--stride', type=int, help='every Nth column in x and y, offset N//2')
    p.add_argument('--hot-cool', action='store_true', help='add the hottest and coolest column at z=0')
    p.add_argument('--extension', nargs='+', choices=['clamp', 'loglin'], required=True)
    p.add_argument('--abundance', type=float, default=8.70)
    p.add_argument('--rays', type=int, default=3)
    p.add_argument('--seconds', type=float, default=1800.)
    p.add_argument('--maxiter', type=int, default=1500)
    args = p.parse_args()
    if platform.system() != 'Darwin':
        raise RuntimeError('Mac-only experiment')
    from dataclasses import replace
    from config.constants import require_codex_path
    from multi3d_read import read_atmos_multi3d
    from pipeline.bride_columns import extend_collisions, mu1_spectrum, stagger_column
    from pipeline.bride_parameter_checks import assert_collision_coverage, marcs_filename
    from pipeline.bride_reference import fingerprint, read_marcs_structure
    from pipeline.multi_atom import read_multi_atom
    from pipeline.wavelength_util import air_to_vac

    base = Path(require_codex_path('engines.bride_eval'))
    mesh, atm = base / f'mesh.{SNAPSHOT}', base / f'atm3d.{SNAPSHOT}'
    cube = read_atmos_multi3d(str(mesh), str(atm))
    atom_path = require_codex_path('experiments.bride_oxygen_atom')
    atom = replace(read_multi_atom(atom_path), abundance=args.abundance)
    marcs = require_codex_path('grids.marcs_standard') / marcs_filename(5750, 4.5)
    abundances = read_marcs_structure(marcs)['abundances_dex']   # solar mix for the passive atoms

    cols = [tuple(int(v) for v in c.split(',')) for c in args.columns]
    ny, nx = cube['temperature'].shape[1:]
    if args.stride:
        s = args.stride
        cols += [(ix, iy) for iy in range(s // 2, ny, s) for ix in range(s // 2, nx, s)]
    if args.hot_cool:
        k = int(np.argmin(abs(cube['z'])))
        t = cube['temperature'][k]
        for idx in (np.argmax(t), np.argmin(t)):
            iy, ix = np.unravel_index(idx, t.shape)
            cols.append((int(ix), int(iy)))
    cols = list(dict.fromkeys(cols))
    if not cols:
        raise ValueError('No columns selected')

    wave = np.arange(WINDOW_AIR_A[0], WINDOW_AIR_A[1] + STEP_A / 2, STEP_A)
    vac = air_to_vac(wave) * .1
    args.output.mkdir(parents=True, exist_ok=False)
    summary = dict(schema='codex.bride_stagger_columns/1', ticket='RYA-1282',
                   scientific_admission='HOLD: G5 1.5D experiment; fixed A(O), no fit, no full-3D claim',
                   code_commit=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
                   generator=fingerprint(__file__), adapter=fingerprint(ROOT / 'pipeline/bride_columns.py'),
                   atom=fingerprint(atom_path), abundance_O=args.abundance, snapshot=SNAPSHOT,
                   mesh=fingerprint(mesh), atmosphere=fingerprint(atm), passive_abundances=fingerprint(marcs),
                   logg_cgs=SNAPSHOT_LOGG, rays=args.rays, vturb_m_s=0., maxiter=args.maxiter,
                   horizontal_velocity='dropped (1.5D, mu = 1 only)', window_air_A=list(WINDOW_AIR_A),
                   step_A=STEP_A, columns=[])
    save = lambda: (args.output / 'summary.json').write_text(json.dumps(summary, indent=1) + '\n')
    save()
    for ix, iy in cols:
        col = stagger_column(cube, ix, iy)
        t_hi = float(col['temperature_K'].max())
        for mode in args.extension:
            name = f'x{ix:02d}_y{iy:02d}_{mode}'
            rec = dict(name=name, ix=ix, iy=iy, extension=mode, T_max_K=t_hi,
                       T_surface_K=float(col['temperature_K'][int(np.argmin(abs(col['height_m'])))]),
                       layers_above_table=int((col['temperature_K'] > 10000.).sum()))
            t0 = time.perf_counter()
            try:
                use = extend_collisions(atom, t_hi * 1.0001, mode) if t_hi > 10000. else atom
                assert_collision_coverage(use, col['temperature_K'])
                r = mu1_spectrum(col, use, abundances_dex=abundances, logg_cgs=SNAPSHOT_LOGG,
                                 wave_vac_nm=vac, rays=args.rays, seconds=args.seconds,
                                 maxiter=args.maxiter)
                np.savez_compressed(args.output / f'{name}.npz', wavelength_air_A=wave, line=r['line'],
                                    continuum=r['continuum'], populations=r['populations'], nstar=r['nstar'],
                                    temperature_K=col['temperature_K'], height_m=col['height_m'])
                (args.output / f'{name}_history.json').write_text(json.dumps(r['history']) + '\n')
                rec.update(status='CONVERGED', iterations=r['iterations'], continuum_iterations=r['continuum_iterations'],
                           seconds_solve=r['seconds_solve'], seconds_total=r['seconds_total'],
                           final=r['history'][-1])
            except Exception as e:  # keep every failure, with its iteration history
                if getattr(e, 'history', None):
                    (args.output / f'{name}_history.json').write_text(json.dumps(e.history) + '\n')
                    rec['final'] = e.history[-1]
                rec.update(status='FAILED', error=f'{type(e).__name__}: {e}', seconds_total=time.perf_counter() - t0)
            summary['columns'].append(rec)
            save()
            print(name, rec['status'], round(rec['seconds_total'], 1), rec.get('iterations', ''), rec.get('error', ''), flush=True)
    if any(c['status'] != 'CONVERGED' for c in summary['columns']):
        raise SystemExit(1)


if __name__ == '__main__':
    main()
