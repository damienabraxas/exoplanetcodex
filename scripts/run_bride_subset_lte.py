"""RYA-1282 Bride G5: the LTE partner of run_bride_subset.py on the same STAGGER columns.

Same column builder, atom, collision extension, window and mu = 1 output; oxygen held in
LTE (no SE). Background radiation is iterated to a declared J tolerance, never a fixed
count. Gives the per-column and averaged NLTE/LTE ratio when paired with the NLTE subset.
"""
import argparse, json, subprocess, sys, time
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / 'scripts')]


def lte_spectrum(col, atom, *, abundances_dex, logg_cgs, wave_vac_nm, rtol=1e-6, maxiter=500):
    from pipeline.bride_columns import column_context, converge_continuum
    kw = dict(abundances_dex=abundances_dex, logg_cgs=logg_cgs, vturb_m_s=0., rays=3)
    ctx = column_context(col, atom, lte=True, **kw)
    for it in range(1, maxiter + 1):
        dj = ctx.formal_sol_gamma_matrices(lambdaIterate=True).dJMax
        if dj < rtol:
            break
    else:
        raise RuntimeError(f'LTE background J not converged: dJ {dj:.2e}')
    line = np.asarray(ctx.compute_rays(wavelengths=wave_vac_nm, mus=1.0)).ravel()
    nstar = ctx.activeAtoms[0].nStar.copy() if ctx.activeAtoms else ctx.kwargs['eqPops'].atomicPops['O'].nStar.copy()
    cont_ctx = column_context(col, atom, ng=False, fixed_populations=nstar, continuum_only=True, **kw)
    converge_continuum(cont_ctx)
    cont = np.asarray(cont_ctx.compute_rays(wavelengths=wave_vac_nm, mus=1.0)).ravel()
    return line, cont, it


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--stride', type=int, default=5)
    p.add_argument('--offset', type=int, default=None)
    a = p.parse_args()
    from pipeline.bride_columns import extend_collisions, stagger_column
    from pipeline.wavelength_util import air_to_vac
    from run_bride_column_rule_validation import _setup
    from run_bride_stagger_columns import SNAPSHOT_LOGG, STEP_A, WINDOW_AIR_A
    o = a.stride // 2 if a.offset is None else a.offset
    cols = [(ix, iy) for iy in range(o, 80, a.stride) for ix in range(o, 80, a.stride)]
    a.output.mkdir(parents=True, exist_ok=False)
    cube, atom, ab = _setup()
    wave = np.arange(WINDOW_AIR_A[0], WINDOW_AIR_A[1] + STEP_A / 2, STEP_A); vac = air_to_vac(wave) * .1
    summary = dict(schema='codex.bride_subset_lte/1', ticket='RYA-1282', stride=a.stride, offset=o, ncols=len(cols),
                   code_commit=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
                   collision_extension='loglin', abundance_O=8.70, columns=[])
    t0 = time.perf_counter()
    for xy in cols:
        col = stagger_column(cube, *xy); t_hi = float(col['temperature_K'].max())
        use = extend_collisions(atom, t_hi * 1.0001, 'loglin') if t_hi > 1e4 else atom
        line, cont, it = lte_spectrum(col, use, abundances_dex=ab, logg_cgs=SNAPSHOT_LOGG, wave_vac_nm=vac)
        np.savez_compressed(a.output / f'x{xy[0]:02d}_y{xy[1]:02d}.npz', wavelength_air_A=wave, line=line, continuum=cont)
        summary['columns'].append(dict(column=list(xy), status='LTE', j_iterations=it))
        summary['wall_seconds'] = time.perf_counter() - t0
        (a.output / 'summary.json').write_text(json.dumps(summary, indent=1) + '\n')
        print(len(summary['columns']), xy, it, flush=True)


if __name__ == '__main__':
    main()
