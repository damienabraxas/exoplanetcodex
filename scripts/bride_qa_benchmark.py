"""Bride QA benchmark: solve every benchmark model and grade the model's PHYSICS (RYA-1282).

Each benchmark record is one exact node of the Gerber, Bergemann et al. 2023 (A&A 669, A43)
MARCS O grid (extracted on Sirius by scripts/prepare_bride_oxygen_reference.py). The Bride
is run on the same MARCS model, atom and A(O) through its DEFAULT path, and graded on:

  Q1 convergence      strict joint rule reached; iterations; Ng rejections
  Q2 LTE limit        max |b - 1| over all levels at the 3 deepest layers
  Q3 conservation     max |sum(n) / sum(n*) - 1| over depth (O I + O II particles)
  Q4 stat. equil.     physical rate residual after a FRESH formal solution at fixed n
  Q5 detailed balance collisional C_ij n*_j vs C_ji n*_i, max relative mismatch
  Q6 departures       max |dlog10 b| vs the record: triplet levels and all levels on
                      logtau500 [-3, 0], and all levels/all depths
  Q7 777 triplet      flux EWs (5-point Gauss-Legendre in mu), NLTE and LTE

The benchmark record is NEVER passed to the solve. Inputs are hash-checked against the
record (atom and MARCS file), so a grade can never be against a different model.

  source config/hosts/mac.env
  OMP_NUM_THREADS=1 PYTHONPATH=.:$CODEX_BRIDE_EVAL $CODEX_BRIDE_EVAL/.venv/bin/python \
      scripts/bride_qa_benchmark.py --bench <dir of records> --output <fresh dir> [--workers 2]
"""
import argparse, json, subprocess, sys, time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import replace
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT)]
TRIPLET_LEVELS = (5, 7, 8, 9)            # 0-based: 3s 5S2, 3p 5P1,2,3 (o41f)
EW_WINDOWS_AIR_A = ((7768.5, 7773.055), (7773.055, 7774.777), (7774.777, 7778.5))


def grade(record_path):
    from config.constants import require_codex_path
    from pipeline.bride_reference import fingerprint, read_marcs_structure
    from pipeline.bride_columns import column_context, marcs_column, solve_column
    from pipeline.bride_numerical_validation import physical_rate_matrix
    from pipeline.multi_atom import read_multi_atom
    from pipeline.wavelength_util import air_to_vac
    rec = json.loads(Path(record_path).read_text())
    ref = rec['reference']
    atom_path = require_codex_path('experiments.bride_oxygen_atom')
    model_path = require_codex_path('grids.marcs_standard') / rec['atmosphere_file']['name']
    if fingerprint(atom_path)['sha256'] != rec['atom']['sha256']:
        raise ValueError('Atom differs from the benchmark record')
    if fingerprint(model_path)['sha256'] != rec['atmosphere_file']['sha256']:
        raise ValueError('MARCS model differs from the benchmark record')
    atom = replace(read_multi_atom(atom_path), abundance=float(ref['abundance']))
    S = read_marcs_structure(model_path)
    col = marcs_column(S)
    bref = np.asarray(ref['departures'], float)
    bref = bref if bref.shape[0] == len(atom.levels) else bref.T
    lt = np.log10(np.asarray(ref['tau'], float))
    if len(lt) != S['ndep']:
        raise ValueError('Record depth grid differs from the model')
    wave = np.arange(7766.0, 7781.0 + 0.005, 0.01); vac = air_to_vac(wave) * .1
    xg, wg = np.polynomial.legendre.leggauss(5); mu = 0.5 * (xg + 1); w = 0.5 * wg

    def flux_ews(ctx):
        I = np.asarray(ctx.compute_rays(wavelengths=vac, mus=list(mu)))
        I = I.reshape(len(wave), len(mu)) if I.shape[0] == len(wave) else I.T
        F = (I * (mu * w)[None, :]).sum(1)
        r = 1 - F / np.interp(wave, [wave[0], wave[-1]], [F[0], F[-1]])
        return [float(np.trapezoid(r[(wave >= a) & (wave < b)], wave[(wave >= a) & (wave < b)]) * 1e3)
                for a, b in EW_WINDOWS_AIR_A]

    kw = dict(abundances_dex=S['abundances_dex'], logg_cgs=S['logg_cgs'], vturb_m_s=S['vturb_m_s'], rays=3)
    out = dict(record=Path(record_path).name, model=rec['atmosphere_file']['name'], A_O=float(ref['abundance']),
               teff=float(S.get('teff', np.nan)) if 'teff' in S else None, logg=float(S['logg_cgs']))
    t0 = time.perf_counter()
    lte = column_context(col, atom, lte=True, **kw)
    for _ in range(500):
        if lte.formal_sol_gamma_matrices(lambdaIterate=True).dJMax < 1e-6:
            break
    out['ew_lte_flux_mA'] = flux_ews(lte)
    ctx = column_context(col, atom, ng=True, **kw)
    try:
        n, hist, t_solve, _ = solve_column(ctx, maxiter=6000, seconds=3000.)
    except Exception as e:                                    # a failure is a grade, not a skip
        out.update(Q1_converged=False, error=f'{type(e).__name__}: {e}', seconds=time.perf_counter() - t0)
        return out
    a = ctx.activeAtoms[0]
    b = n / a.nStar
    out.update(Q1_converged=True, Q1_iterations=len(hist),
               Q1_ng_rejections=int(sum(bool(h.get('rejected_invalid_ng')) for h in hist)))
    out['Q2_lte_limit_max_abs_b_minus_1'] = float(np.max(np.abs(b[:, -3:] - 1)))
    out['Q3_conservation_max'] = float(np.max(np.abs(n.sum(0) / a.nStar.sum(0) - 1)))
    ctx.formal_sol_gamma_matrices(lambdaIterate=True)          # fresh radiation, n unchanged
    np.testing.assert_array_equal(a.n, n)
    R = physical_rate_matrix(a)
    bal = np.einsum('ijd,jd->id', R, n); flow = np.einsum('ijd,jd->id', abs(R), n)
    out['Q4_rate_residual_fresh'] = float(np.max(abs(bal) / np.maximum(flow, np.finfo(float).tiny)))
    a.compute_collisions(fillDiagonal=False)
    C = np.array(a.C); ns = np.array(a.nStar)
    fwd = C * ns[None, :, :]; rev = np.transpose(C, (1, 0, 2)) * ns[:, None, :]
    m = (fwd > 0) | (rev > 0); np.fill_diagonal(m[:, :, 0], False)
    out['Q5_detailed_balance_max_rel'] = float(np.max(abs(fwd - rev)[m] / np.maximum(fwd, rev)[m]))
    d = np.abs(np.log10(np.maximum(b, 1e-300)) - np.log10(np.maximum(bref, 1e-300)))
    phot = (lt >= -3) & (lt <= 0)
    out['Q6_triplet_phot_max_dlogb'] = float(d[list(TRIPLET_LEVELS)][:, phot].max())
    out['Q6_all_phot_max_dlogb'] = float(d[:, phot].max())
    out['Q6_all_phot_worst_level_1based'] = int(np.unravel_index(np.argmax(np.where(phot[None, :], d, -1)), d.shape)[0] + 1)
    out['Q6_all_depth_max_dlogb'] = float(d.max())
    out['Q6_all_depth_worst_logtau'] = float(lt[np.unravel_index(np.argmax(d), d.shape)[1]])
    out['ew_nlte_flux_mA'] = flux_ews(ctx)
    out['nlte_over_lte'] = [x / y for x, y in zip(out['ew_nlte_flux_mA'], out['ew_lte_flux_mA'])]
    out['seconds'] = time.perf_counter() - t0
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--bench', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--workers', type=int, default=2)
    p.add_argument('--only', nargs='*', default=None, help='record filenames to run')
    a = p.parse_args()
    from pipeline.bride_reference import fingerprint
    records = sorted(a.bench.glob('*.json'))
    if a.only:
        records = [r for r in records if r.name in set(a.only)]
    if not records:
        raise ValueError('No benchmark records')
    a.output.mkdir(parents=True, exist_ok=False)
    summary = dict(schema='codex.bride_qa_benchmark/1', ticket='RYA-1282',
                   code_commit=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
                   dirty=bool(subprocess.check_output(['git', 'status', '--porcelain', '--untracked-files=no'], cwd=ROOT, text=True).strip()),
                   generator=fingerprint(__file__), benchmark='Gerber+2023 A&A 669 A43 O MARCS grid, exact nodes',
                   records=[fingerprint(r) for r in records], results=[])
    with ProcessPoolExecutor(a.workers) as pool:
        futs = {pool.submit(grade, str(r)): r for r in records}
        for f in as_completed(futs):
            try:
                res = f.result()
            except Exception as e:
                res = dict(record=futs[f].name, Q1_converged=False, error=f'{type(e).__name__}: {e}')
            summary['results'].append(res)
            (a.output / 'summary.json').write_text(json.dumps(summary, indent=1) + '\n')
            print(json.dumps({k: (round(v, 6) if isinstance(v, float) else v) for k, v in res.items()
                              if k.startswith(('record', 'Q', 'error', 'nlte_over'))}), flush=True)


if __name__ == '__main__':
    main()
