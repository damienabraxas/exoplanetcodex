"""RYA-1282 Bride G5: measure the observable stopping rule on many STAGGER columns.

Each column is resampled to a fixed log tau500 grid, solved in NLTE (O41f, loglin
collision extension, fixed A(O) = 8.70) to the STRICT 1e-6 joint rule, with the disk-centre
777 triplet equivalent widths recorded every --every iterations. Any observable tolerance
can then be scored post hoc against the strictly converged widths of the same column.
Columns run in parallel worker processes (one Lightweaver thread each).
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'scripts'))

EW_WINDOWS_AIR_A = ((7771.2, 7773.0), (7773.0, 7774.8), (7774.8, 7776.5))


def _setup():
    from dataclasses import replace
    from config.constants import require_codex_path
    from multi3d_read import read_atmos_multi3d
    from pipeline.bride_parameter_checks import marcs_filename
    from pipeline.bride_reference import read_marcs_structure
    from pipeline.multi_atom import read_multi_atom
    from run_bride_stagger_columns import SNAPSHOT
    base = Path(require_codex_path('engines.bride_eval'))
    cube = read_atmos_multi3d(str(base / f'mesh.{SNAPSHOT}'), str(base / f'atm3d.{SNAPSHOT}'))
    atom = replace(read_multi_atom(require_codex_path('experiments.bride_oxygen_atom')), abundance=8.70)
    ab = read_marcs_structure(require_codex_path('grids.marcs_standard') / marcs_filename(5750, 4.5))['abundances_dex']
    return cube, atom, ab


def solve_one(xy, grid, every, maxiter, out):
    from pipeline.bride_columns import (column_context, continuum_logtau, extend_collisions,
                                        resample_logtau, solve_column, stagger_column)
    from pipeline.wavelength_util import air_to_vac
    from run_bride_stagger_columns import SNAPSHOT_LOGG, STEP_A, WINDOW_AIR_A
    cube, atom, ab = _setup()
    wave = np.arange(WINDOW_AIR_A[0], WINDOW_AIR_A[1] + STEP_A / 2, STEP_A)
    vac = air_to_vac(wave) * .1

    def ew(ctx):
        intensity = np.asarray(ctx.compute_rays(wavelengths=vac, mus=1.0)).ravel()
        r = 1 - intensity / np.interp(wave, [wave[0], wave[-1]], [intensity[0], intensity[-1]])
        return [float(np.trapezoid(r[(wave >= lo) & (wave < hi)], wave[(wave >= lo) & (wave < hi)]) * 1e3)
                for lo, hi in EW_WINDOWS_AIR_A]

    rec = dict(column=list(xy))
    t0 = time.perf_counter()
    try:
        col = stagger_column(cube, *xy)
        lt = continuum_logtau(col, atom, abundances_dex=ab, logg_cgs=SNAPSHOT_LOGG)
        rc = resample_logtau(col, lt, grid)
        t_hi = float(rc['temperature_K'].max())
        use = extend_collisions(atom, t_hi * 1.0001, 'loglin') if t_hi > 1e4 else atom
        ctx = column_context(rc, use, ng=True, abundances_dex=ab, logg_cgs=SNAPSHOT_LOGG, vturb_m_s=0., rays=3)
        n, hist, t_solve, _ = solve_column(ctx, maxiter=maxiter, seconds=3600., monitor=ew, every=every)
        final = ew(ctx)
        marks = [(e['iteration'], e['monitor']) for e in hist if 'monitor' in e]
        rec.update(status='CONVERGED', iterations=len(hist), seconds_solve=t_solve,
                   seconds_per_iteration=t_solve / len(hist), T_max_K=t_hi,
                   surface_T_K=float(np.interp(0., lt, col['temperature_K'])),
                   ew_final_mA=final, ew_marks=marks)
        np.savez_compressed(out / f'x{xy[0]:02d}_y{xy[1]:02d}.npz', populations=n,
                            nstar=ctx.activeAtoms[0].nStar.copy(), logtau500=grid)
    except Exception as e:
        rec.update(status='FAILED', error=f'{type(e).__name__}: {e}')
        hist = getattr(e, 'history', None)
        if hist:   # keep the failure's trajectory: strict-rule failure is not observable failure
            rec['iterations'] = len(hist)
            rec['final'] = {k: v for k, v in hist[-1].items() if k != 'monitor'}
            rec['ew_marks'] = [(h['iteration'], h['monitor']) for h in hist if 'monitor' in h]
            rec['tail'] = [{k: v for k, v in h.items() if k != 'monitor'} for h in hist[-60:]]
    rec['seconds_total'] = time.perf_counter() - t0
    return rec


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--ncols', type=int, default=12)
    p.add_argument('--seed', type=int, default=1282)
    p.add_argument('--workers', type=int, default=4)
    p.add_argument('--npts', type=int, default=80)
    p.add_argument('--every', type=int, default=25)
    p.add_argument('--maxiter', type=int, default=6000)
    p.add_argument('--columns', nargs='*', default=None, help='explicit ix,iy list instead of the seeded draw')
    a = p.parse_args()
    rng = np.random.default_rng(a.seed)
    cols = [(int(x), int(y)) for x, y in zip(rng.integers(0, 80, a.ncols), rng.integers(0, 80, a.ncols))]
    cols += [(26, 5), (71, 20)]          # hottest / coolest at the surface
    if a.columns is not None:
        cols = [tuple(int(v) for v in c.split(',')) for c in a.columns]
    cols = list(dict.fromkeys(cols))
    grid = np.linspace(-6., 2., a.npts)
    a.output.mkdir(parents=True, exist_ok=False)
    summary = dict(schema='codex.bride_column_rule_validation/1', ticket='RYA-1282', seed=a.seed,
                   grid=dict(npts=a.npts, top=-6., bottom=2.), every=a.every, strict_rtol=1e-6,
                   collision_extension='loglin', abundance_O=8.70, ew_windows_air_A=EW_WINDOWS_AIR_A,
                   scientific_admission='HOLD: numerical rule measurement only', columns=[])
    save = lambda: (a.output / 'summary.json').write_text(json.dumps(summary, indent=1) + '\n')
    save()
    with ProcessPoolExecutor(a.workers) as pool:
        futs = [pool.submit(solve_one, xy, grid, a.every, a.maxiter, a.output) for xy in cols]
        for f in as_completed(futs):
            rec = f.result()
            summary['columns'].append(rec)
            save()
            print(rec['column'], rec['status'], rec.get('iterations', ''), round(rec['seconds_total'], 1),
                  rec.get('error', ''), flush=True)


if __name__ == '__main__':
    main()
