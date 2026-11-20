"""RYA-1282: the column path's population operator is the RYA-1281 receipted MALI plugin.

Measured 2026-11-20 (cap top, TS metals, 41-level O): solar 5750/4.5 converges in 10
steps (pure Lambda iteration: 426-571, contraction ~0.998); 5750/3.5 in 15 (Lambda+Ng
stalls at 3e-4 for 4000 steps); 5750/4.5/-2 in 34 (Lambda 958, LTE limit 7.7e-5).
Solar MALI vs the Lambda answer converged to 1e-7: 5.5e-5 dex (= the Lambda answer's own
remaining error at rho ~0.998). The operator changes the update, never the physical rates.
"""
import numpy as np
import pytest

pytest.importorskip('lightweaver')


def _solar(operator):
    from dataclasses import replace
    import pipeline.bride_columns as BC
    from config.constants import require_codex_path
    from pipeline.multi_atom import read_multi_atom
    from pipeline.bride_parameter_checks import marcs_filename
    from pipeline.bride_reference import read_marcs_structure
    try:
        atom = replace(read_multi_atom(require_codex_path('experiments.bride_oxygen_atom')), abundance=8.70)
        S = read_marcs_structure(require_codex_path('grids.marcs_standard') / marcs_filename(5750, 4.5))
        require_codex_path('experiments.bride_mali_cap_plugin')
    except Exception as e:                                         # noqa: BLE001
        pytest.skip(f'O atom / MARCS / MALI plugin not held here: {e}')
    return BC.column_context(BC.marcs_column(S), atom, abundances_dex=S['abundances_dex'],
                             logg_cgs=S['logg_cgs'], vturb_m_s=S['vturb_m_s'], rays=3, operator=operator)


def _fresh_residual(ctx, n):
    import pipeline.bride_columns as BC
    from pipeline.bride_numerical_validation import physical_rate_matrix
    BC.converge_radiation(ctx, tol=1e-8)
    R = physical_rate_matrix(ctx.activeAtoms[0])
    return float(np.max(abs(np.einsum('ijd,jd->id', R, n)) / np.einsum('ijd,jd->id', abs(R), n)))


def test_mali_converges_the_solar_column_and_lambda_iteration_does_not_in_the_same_steps():
    import pipeline.bride_columns as BC
    ctx = _solar('mali')
    n, h, *_ = BC.solve_column(ctx, rtol=1e-7, maxiter=20, seconds=600.)
    assert len(h) <= 20
    assert _fresh_residual(ctx, n) < 1e-7
    assert np.max(abs(n.sum(0) / ctx.activeAtoms[0].nTotal - 1)) < 1e-12
    lam = _solar('lambda')
    with pytest.raises(RuntimeError):                             # negative control
        BC.solve_column(lam, rtol=1e-7, maxiter=len(h), seconds=600.)


def test_the_mali_refresh_leaves_the_physical_rate_matrix_unchanged():
    import pipeline.bride_columns as BC
    from pipeline.bride_numerical_validation import physical_rate_matrix
    ctx = _solar('mali')
    BC.converge_radiation(ctx, tol=1e-10)
    R0 = physical_rate_matrix(ctx.activeAtoms[0]).copy()
    ctx.formal_sol_gamma_matrices(lambdaIterate=False, extraParams={'BrideCapLengthM': ctx.atmos.pyAtmos.zUpperBc.length_m})
    R1 = physical_rate_matrix(ctx.activeAtoms[0])
    scale = np.max(abs(R0), axis=1, keepdims=True)
    assert np.max(abs(R1 - R0) / scale) < 1e-9                    # J moved < 1e-10: rates follow


def test_unknown_operator_is_refused():
    with pytest.raises(ValueError):
        _solar('newton')
