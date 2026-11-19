"""RYA-1282 Bride G5: 1.5D O NLTE over a regular subset of STAGGER columns (Mac, parallel).

Native snapshot depth grid (no resampling: an 80-point log tau grid biased the 777 EWs
by -0.3%, 160 points matched native but saved nothing). Each column stops on the
observable rule (777 triplet EWs change < --observable-tol over --every iterations;
measured on 12 columns: median 150 iterations, max true EW error 7.5e-4 at 1e-4) or on the
strict 1e-6 joint rule, whichever comes first. Collisions above 10000 K: loglin
(clamp and loglin gave identical 777 EWs on the hottest and coolest columns).
Writes each column's absolute mu = 1 line and continuum intensity for averaging.
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


def solve_one(xy, args, out):
    from pipeline.bride_columns import ew_monitor, extend_collisions, mu1_spectrum, stagger_column
    from pipeline.wavelength_util import air_to_vac
    from run_bride_column_rule_validation import _setup
    from run_bride_stagger_columns import SNAPSHOT_LOGG, STEP_A, WINDOW_AIR_A
    cube, atom, ab = _setup()
    wave = np.arange(WINDOW_AIR_A[0], WINDOW_AIR_A[1] + STEP_A / 2, STEP_A)
    vac = air_to_vac(wave) * .1
    name = f'x{xy[0]:02d}_y{xy[1]:02d}'
    rec = dict(column=list(xy))
    t0 = time.perf_counter()
    try:
        col = stagger_column(cube, *xy)
        t_hi = float(col['temperature_K'].max())
        use = extend_collisions(atom, t_hi * 1.0001, 'loglin') if t_hi > 1e4 else atom
        r = mu1_spectrum(col, use, abundances_dex=ab, logg_cgs=SNAPSHOT_LOGG, wave_vac_nm=vac,
                         maxiter=args.maxiter, seconds=args.seconds, monitor=ew_monitor(wave, vac),
                         every=args.every, observable_tol=args.observable_tol)
        np.savez_compressed(out / f'{name}.npz', wavelength_air_A=wave, line=r['line'], continuum=r['continuum'])
        rec.update(status='CONVERGED', stopped_by=r['stopped_by'], iterations=r['iterations'],
                   seconds_solve=r['seconds_solve'], ew_mA=r['history'][-1].get('monitor'),
                   T_max_K=t_hi)
    except Exception as e:
        rec.update(status='FAILED', error=f'{type(e).__name__}: {e}')
        hist = getattr(e, 'history', None)
        if hist:
            rec['iterations'] = len(hist)
            rec['ew_marks'] = [(h['iteration'], h['monitor']) for h in hist if 'monitor' in h]
    rec['seconds_total'] = time.perf_counter() - t0
    return rec


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--stride', type=int, default=5)
    p.add_argument('--offset', type=int, default=None, help='grid start in x and y (default stride//2)')
    p.add_argument('--workers', type=int, default=4)
    p.add_argument('--every', type=int, default=25)
    p.add_argument('--observable-tol', type=float, default=1e-4)
    p.add_argument('--maxiter', type=int, default=6000)
    p.add_argument('--seconds', type=float, default=3600.)
    args = p.parse_args()
    s = args.stride
    o = s // 2 if args.offset is None else args.offset
    cols = [(ix, iy) for iy in range(o, 80, s) for ix in range(o, 80, s)]
    args.output.mkdir(parents=True, exist_ok=False)
    from pipeline.bride_reference import fingerprint
    import subprocess
    summary = dict(schema='codex.bride_subset/1', ticket='RYA-1282', stride=s, offset=o, ncols=len(cols),
                   code_commit=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
                   adapter=fingerprint(ROOT / 'pipeline/bride_columns.py'), generator=fingerprint(__file__),
                   observable_tol=args.observable_tol, every=args.every, collision_extension='loglin',
                   abundance_O=8.70, depth_grid='native snapshot (240 layers)',
                   horizontal_velocity='dropped (1.5D, mu = 1)',
                   scientific_admission='HOLD: G5 representativeness/cost experiment; fixed A(O); no fit', columns=[])
    save = lambda: (args.output / 'summary.json').write_text(json.dumps(summary, indent=1) + '\n')
    save()
    t0 = time.perf_counter()
    with ProcessPoolExecutor(args.workers) as pool:
        futs = [pool.submit(solve_one, xy, args, args.output) for xy in cols]
        for k, f in enumerate(as_completed(futs), 1):
            rec = f.result()
            summary['columns'].append(rec)
            summary['wall_seconds'] = time.perf_counter() - t0
            save()
            print(k, rec['column'], rec['status'], rec.get('stopped_by', ''), rec.get('iterations', ''),
                  round(rec['seconds_total'], 1), rec.get('error', ''), flush=True)


if __name__ == '__main__':
    main()
