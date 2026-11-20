"""RYA-1282: an Ng-accelerated SE step that breaks particle conservation is rejected and the
same equation re-solved unaccelerated; never renormalised. STAGGER column x12_y77 failed G5
(2.7e-10 drift at iteration 740 with Ng; 1.8e-15 over 800 iterations without Ng)."""
from types import SimpleNamespace as S
import numpy as np
import pytest
from pipeline.bride_parameter_checks import positive_stat_equil

GOOD = np.array([[0.6, 0.3], [0.4, 0.7]])
DRIFT = GOOD * (1 + 3e-10)                       # positive, finite, but sum != nTotal


def _ctx(sequence):
    active = S(n=GOOD.copy(), nTotal=np.ones(2))
    calls = iter(sequence)
    ctx = S(cleared=0)
    ctx.stat_equil = lambda: active.n.__setitem__(slice(None), next(calls))
    ctx.clear_ng = lambda: setattr(ctx, 'cleared', ctx.cleared + 1)
    return ctx, active


def test_conservation_breaking_ng_step_is_rejected_and_resolved_unaccelerated():
    ctx, active = _ctx([DRIFT, GOOD])
    assert positive_stat_equil(ctx, active, GOOD, reject_invalid_ng=True, conservation_tol=1e-10)
    assert ctx.cleared == 1
    np.testing.assert_array_equal(active.n, GOOD)  # the retry's result, not a rescaled DRIFT


def test_without_a_tolerance_the_old_contract_is_unchanged():
    ctx, active = _ctx([DRIFT])
    assert not positive_stat_equil(ctx, active, GOOD, reject_invalid_ng=True)
    assert ctx.cleared == 0


def test_an_unaccelerated_retry_that_still_drifts_fails():
    ctx, active = _ctx([DRIFT, DRIFT])
    with pytest.raises(ValueError):
        positive_stat_equil(ctx, active, GOOD, reject_invalid_ng=True, conservation_tol=1e-10)


def test_ng_on_J_is_exact_for_a_two_mode_linear_iteration_and_refuses_bad_steps():
    """J at fixed populations follows J <- A + M J (scattering). With two decaying modes,
    an order-2 Ng step from four iterates lands on the fixed point to rounding."""
    from pipeline.bride_columns import ng_extrapolate
    rng = np.random.default_rng(1282)
    Q, _ = np.linalg.qr(rng.normal(size=(6, 6)))
    M = Q @ np.diag([0.995, 0.9, 0, 0, 0, 0]) @ Q.T           # albedo-like slow mode
    fixed = 1 + rng.random(6)
    A = fixed - M @ fixed
    x = np.ones(6) * 3; its = []
    for _ in range(4):
        x = A + M @ x; its.append(x)
    np.testing.assert_allclose(ng_extrapolate(its), fixed, rtol=1e-9)
    assert np.max(abs(its[-1] / fixed - 1)) > 1e-2               # plain iterate is far off
    flat = [np.ones(6)] * 4
    assert ng_extrapolate(flat) is None                           # singular -> plain iterate stands


def test_ali_converges_the_scattering_dominated_metal_poor_column_plain_lambda_does_not():
    """MARCS 5750/3.5/-2, LTE populations, cap top. With the full TS background (hbop
    dissolved Lyman continuum + H+H, 789fd36c) the most scattering-dominated model held:
    albedo 0.982 at 122.5 nm (5750/4.5/-2 fell from 0.998 to 0.949 once those were added).
    Measured 2026-11-20: ALI+Ng 63 formal solutions to 1e-10; plain Lambda after 63 passes
    still dJ 1.8e-2, J off by up to 89%."""
    pytest.importorskip('lightweaver')
    from dataclasses import replace
    import pipeline.bride_columns as BC
    from config.constants import require_codex_path
    from pipeline.multi_atom import read_multi_atom
    from pipeline.bride_parameter_checks import marcs_filename
    from pipeline.bride_reference import read_marcs_structure
    try:
        atom = replace(read_multi_atom(require_codex_path('experiments.bride_oxygen_atom')), abundance=7.1)
        S = read_marcs_structure(require_codex_path('grids.marcs_standard') / marcs_filename(5750, 3.5, -2.0, 0.4))
    except Exception as e:                                         # noqa: BLE001
        pytest.skip(f'O atom / MARCS not held here: {e}')
    def ctx():
        return BC.column_context(BC.marcs_column(S), atom, lte=True, abundances_dex=S['abundances_dex'],
                                 logg_cgs=S['logg_cgs'], vturb_m_s=S['vturb_m_s'], rays=3)
    c = ctx()
    n = BC.converge_radiation(c, tol=1e-10, maxiter=200)
    assert n <= 200
    le = BC.scattering_lambda_star(c)
    assert 0.97 < le.max() < 1                                    # the scattering-dominated regime
    plain = ctx(); plain.depthData.fill = True
    for _ in range(n):
        top = np.array(plain.spect.J)[:, 0].copy()
        dj = plain.formal_sol_gamma_matrices(lambdaIterate=True).dJMax
        BC.refresh_cap(plain, top)
    assert dj > 1e-6                                              # negative control: plain is far from done
    J_ref = np.array(c.spect.J)
    assert np.max(abs(np.array(plain.spect.J) / J_ref - 1)) > 1e-4
