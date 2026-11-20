"""RYA-1282 Bride G5: 1.5D NLTE = independent per-column solves on a 3D snapshot.

A 1.5D solve treats each vertical column of a 3D snapshot as its own plane-parallel
atmosphere: statistical equilibrium and radiative transfer are solved column by column
with no horizontal radiation exchange. It is NOT full 3D (G6); horizontal velocities
enter only through inclined rays, which a 1D column cannot carry, so this module
declares them dropped and is meant for disk-centre (mu = 1) observables first.

Two pieces:
- `stagger_column` turns one (ix, iy) column of a Multi3D-format snapshot (read by
  the RYA-1118 reader) into the SI arrays a Lightweaver 1D atmosphere takes.
- `column_context` builds the same physics as `bride_reference.make_reference_context`
  (same passive atoms, Ng, formal solver, boundaries) from explicit ne / nHtot / vz.
  It is kept separate so that the pinned reference builder (and the receipts that hash
  it) stays unchanged; `tests/test_bride_columns_rya1282.py` proves the two give
  bit-identical solves on a MARCS column.
"""
from __future__ import annotations

import time

import numpy as np

from pipeline.bride_reference import PASSIVE_ATOM_FACTORIES


def stagger_column(cube: dict, ix: int, iy: int) -> dict:
    """One vertical column, top first, in Lightweaver's geometric convention.

    The held Multi3D snapshot already uses HEIGHT for z (z[0] = +871 km at the top,
    decreasing inward) and vz positive along +z, i.e. upward: near z = 0 the hot gas
    rises (corr(T, vz) = +0.54 at z = 0, +0.70 at -50 km on t5777g44m0005_00091).
    That is Lightweaver's convention, so nothing is flipped. The guards refuse a cube
    in any other orientation rather than guess.
    """
    z = np.asarray(cube['z'], dtype=float)
    if not np.all(np.diff(z) < 0):
        raise ValueError('Snapshot z must be height, top first (decreasing); refusing an unknown orientation')
    sl = (slice(None), int(iy), int(ix))
    temperature = cube['temperature'][sl].copy()
    if not temperature[-1] > temperature[0]:
        raise ValueError('Column is not hotter at depth; the axis order is wrong')
    n_htot = cube['hydrogen_density'][sl] + cube['proton_density'][sl]
    col = dict(source='stagger', ix=int(ix), iy=int(iy), ndep=z.size,
               height_m=z.copy(), temperature_K=temperature,
               ne_m3=cube['electron_density'][sl].copy(), nHtot_m3=n_htot.copy(),
               vlos_m_s=cube['vz'][sl].copy())
    for key in ('temperature_K', 'ne_m3', 'nHtot_m3'):
        if not (np.isfinite(col[key]).all() and np.all(col[key] > 0)):
            raise ValueError(f'Nonpositive or nonfinite {key} in column ({ix},{iy})')
    return col


def marcs_column(structure: dict) -> dict:
    """The same column record from a MARCS structure (for the 1D equality check)."""
    from pipeline.bride_reference import marcs_number_densities
    ne, nh = marcs_number_densities(structure)
    return dict(source='marcs', ndep=structure['ndep'], height_m=structure['height_m'].copy(),
                temperature_K=structure['temperature_K'].copy(), ne_m3=ne, nHtot_m3=nh,
                vlos_m_s=np.zeros(structure['ndep']))


def column_context(col: dict, source_atom, *, abundances_dex, logg_cgs: float, vturb_m_s,
                   rays: int, qscale: int = 1, ng: bool = True, lte: bool = False,
                   fixed_populations=None, continuum_only: bool = False):
    """Lightweaver context for one column; O active (or fixed), passive atoms in LTE."""
    import lightweaver as lw
    from lightweaver import rh_atoms
    from lightweaver.atomic_table import PeriodicTable
    from pipeline.multi_atom import to_lightweaver
    n = col['ndep']
    vturb = np.broadcast_to(np.asarray(vturb_m_s, dtype=float), (n,)).copy()
    atmosphere = lw.Atmosphere.make_1d(scale=lw.ScaleType.Geometric,
        depthScale=col['height_m'].copy(), temperature=col['temperature_K'].copy(),
        ne=col['ne_m3'].copy(), nHTot=col['nHtot_m3'].copy(), vturb=vturb,
        vlos=col['vlos_m_s'].copy(), logG=logg_cgs - np.log10(100.),
        lowerBc=lw.ThermalisedRadiation(), upperBc=lw.ZeroRadiation())
    atmosphere.quadrature(rays)
    model = to_lightweaver(source_atom)
    if continuum_only:
        model.lines = []
    if qscale != 1:
        for line in model.lines:
            line.quadrature.Nlambda = qscale * (line.quadrature.Nlambda - 1) + 1
    abundance = {PeriodicTable[i + 1]: a for i, a in enumerate(abundances_dex)}
    abundance[PeriodicTable['O']] = source_atom.abundance
    radset = lw.RadiativeSet([getattr(rh_atoms, name)() for name in PASSIVE_ATOM_FACTORIES] + [model],
                             abundance=lw.AtomicAbundance(abundance))
    if lte or fixed_populations is not None:
        radset.set_detailed_static('O')
    else:
        radset.set_active('O')
    eq = radset.compute_eq_pops(atmosphere)
    if fixed_populations is not None:
        actual = eq.atomicPops[model.element].n
        if fixed_populations.shape != actual.shape or not np.all(fixed_populations > 0):
            raise ValueError('Invalid fixed population replay')
        actual[:] = fixed_populations
    spect = radset.compute_wavelength_grid()
    options = (lw.NgOptions(Norder=3, Nperiod=5, Ndelay=10, threshold=.05, lowerThreshold=0.)
               if ng and not lte else None)
    return lw.Context(atmosphere, spect, eq, Nthreads=1, formalSolver='piecewise_besser_1d',
                      conserveCharge=False, ngOptions=options)


def solve_column(ctx, *, rtol: float = 1e-6, maxiter: int = 1500, seconds: float = 900., monitor=None, every: int = 50,
                 observable_tol: float | None = None, validate_observable_rule: bool = False):
    """Joint RT + SE iteration, the G4 convergence rule: population change, normalised
    rate residual and J change all below rtol (after >= 3 iterations). Returns the
    populations and the full history. Raises on non-convergence; nothing is clipped.

    OBSERVABLE RULE (G5 speed lever, experimental): with `monitor` (a function of the
    context returning the observables, e.g. the 777 EWs) and `observable_tol`, the solve
    may stop when the observables change by less than `observable_tol` (relative) between
    two checkpoints `every` iterations apart. With `validate_observable_rule` it records
    where that rule would have stopped and carries on to the strict rule, so the rule's
    true error is measured, never assumed. Returns a 4th item: the stop record or None."""
    from pipeline.bride_numerical_validation import physical_rate_matrix
    from pipeline.bride_parameter_checks import positive_stat_equil
    active = ctx.activeAtoms[0]
    start = time.perf_counter()
    history = []
    observable_stop = None
    for iteration in range(maxiter):
        if time.perf_counter() - start > seconds:
            raise TimeoutError('Column budget exceeded')
        old = active.n.copy()
        update = ctx.formal_sol_gamma_matrices(lambdaIterate=True)
        matrix = active.Gamma
        scale = np.max(abs(matrix), axis=1, keepdims=True)
        if np.any(scale == 0):
            raise ValueError('Disconnected SE equation')
        matrix[:] /= scale
        # Safeguarded Ng (RYA-1282): with MULTI-semantics collisions an Ng extrapolation can
        # propose negative populations; reject it, clear Ng, re-solve the same equation once.
        rejected = positive_stat_equil(ctx, active, old, reject_invalid_ng=True)
        n = active.n.copy()
        rates = physical_rate_matrix(active)
        balance = np.einsum('ijd,jd->id', rates, n)
        flow = np.einsum('ijd,jd->id', abs(rates), n)
        rate = float(np.max(abs(balance) / np.maximum(flow, np.finfo(float).tiny)))
        delta = float(np.max(abs(n - old) / n))
        conservation = float(np.max(abs(n.sum(axis=0) / active.nTotal - 1)))
        history.append(dict(iteration=iteration + 1, population_change=delta, J_change=float(update.dJMax),
                            physical_rate_residual=rate, particle_conservation_error=conservation,
                            rejected_invalid_ng=bool(rejected)))
        if conservation > 1e-10:
            raise ValueError('Particle conservation failed')
        if monitor is not None and (iteration + 1) % every == 0:
            history[-1]['monitor'] = monitor(ctx)
            if observable_tol is not None:
                marks = [e['monitor'] for e in history if 'monitor' in e]
                if len(marks) >= 2:
                    prev, now = np.asarray(marks[-2]), np.asarray(marks[-1])
                    change = float(np.max(abs(now - prev) / abs(now)))
                    history[-1]['observable_change'] = change
                    if change < observable_tol and observable_stop is None:
                        observable_stop = dict(iteration=iteration + 1, observable=now.tolist(),
                                               populations=active.n.copy())
                        if not validate_observable_rule:
                            return n, history, time.perf_counter() - start, observable_stop
        if iteration >= 2 and max(delta, rate, float(update.dJMax)) < rtol:
            return n, history, time.perf_counter() - start, observable_stop
    err = RuntimeError('Joint convergence not reached')
    err.history = history
    raise err


def converge_continuum(ctx, *, tol: float = 1e-8, maxiter: int = 1000) -> int:
    for j in range(maxiter):
        if ctx.formal_sol_gamma_matrices().dJMax < tol:
            return j + 1
    raise RuntimeError('Continuum J convergence failed')


def extend_collisions(atom, t_top: float, mode: str):
    """A copy of the atom whose collision tables reach `t_top`, by a DECLARED rule.

    The held O atom's tables stop at 10000 K; STAGGER granules pass that just below
    the surface. Rates there are not in the atomic data, so any solve on those layers
    extends it. Two bracketing rules, never a silent default:
      'clamp'  - hold every rate at its last tabulated value;
      'loglin' - continue log(rate) linearly in log(T) from the last two nodes
                 (a zero rate stays zero).
    Measure the observable under both; adopt neither without that measurement.
    """
    import copy
    if mode not in ('clamp', 'loglin'):
        raise ValueError(f'Unknown collision extension {mode!r}')
    out = copy.deepcopy(atom)
    for c in out.collisions:
        grid = out.temperature if c.temperature is None else c.temperature
        if t_top <= grid[-1]:
            continue
        r = np.asarray(c.rates, dtype=float)
        if mode == 'clamp' or len(grid) < 2 or r[-1] <= 0 or r[-2] <= 0:
            new = r[-1]
        else:
            slope = np.log(r[-1] / r[-2]) / np.log(grid[-1] / grid[-2])
            new = r[-1] * (t_top / grid[-1]) ** slope
        c.temperature = np.append(grid, t_top)
        c.rates = np.append(r, new)
    if out.temperature[-1] < t_top:
        out.temperature = np.append(out.temperature, t_top)
    return out


def mu1_spectrum(col, atom, *, abundances_dex, logg_cgs, wave_vac_nm, rays=3, vturb_m_s=0.,
                 rtol=1e-6, maxiter=1500, seconds=900., departures=None,
                 monitor=None, every=25, observable_tol=None):
    """Solve one column in NLTE and return disk-centre (mu = 1) line and continuum
    intensity on `wave_vac_nm`, the populations, the history and timings."""
    t0 = time.perf_counter()
    kw = dict(abundances_dex=abundances_dex, logg_cgs=logg_cgs, vturb_m_s=vturb_m_s, rays=rays)
    ctx = column_context(col, atom, ng=True, **kw)
    if departures is not None:
        set_departures(ctx, departures)
    n, history, t_solve, stop = solve_column(ctx, rtol=rtol, maxiter=maxiter, seconds=seconds,
                                             monitor=monitor, every=every, observable_tol=observable_tol)
    line = np.asarray(ctx.compute_rays(wavelengths=wave_vac_nm, mus=1.0)).ravel()
    cont_ctx = column_context(col, atom, ng=False, fixed_populations=n, continuum_only=True, **kw)
    j = converge_continuum(cont_ctx)
    cont = np.asarray(cont_ctx.compute_rays(wavelengths=wave_vac_nm, mus=1.0)).ravel()
    return dict(line=line, continuum=cont, populations=n, nstar=ctx.activeAtoms[0].nStar.copy(),
                history=history, iterations=len(history), continuum_iterations=j,
                stopped_by='observable' if stop is not None else 'strict',
                seconds_solve=t_solve, seconds_total=time.perf_counter() - t0)


def continuum_logtau(col, atom, *, abundances_dex, logg_cgs, wavelength_nm: float = 500.):
    """log10 of the LTE background continuum optical depth at `wavelength_nm`, per layer,
    integrated down from the top (first layer: half its own thickness). Diagnostic for
    choosing a depth cut, computed from the column's own background opacity."""
    ctx = column_context(col, atom, abundances_dex=abundances_dex, logg_cgs=logg_cgs,
                         vturb_m_s=0., rays=3, lte=True, ng=False)
    wave = np.asarray(ctx.spect.wavelength)
    k = int(np.argmin(abs(wave - wavelength_nm)))
    chi = np.asarray(ctx.background.chi)[k]
    h = col['height_m']
    dtau = 0.5 * (chi[1:] + chi[:-1]) * np.abs(np.diff(h))
    tau = np.concatenate(([chi[0] * abs(h[1] - h[0]) * 0.5], chi[0] * abs(h[1] - h[0]) * 0.5 + np.cumsum(dtau)))
    return np.log10(tau)


def trim_column(col, keep) -> dict:
    """Column restricted to a contiguous top-first slice `keep` (a boolean mask)."""
    idx = np.flatnonzero(keep)
    if idx.size < 10 or not np.all(np.diff(idx) == 1) or idx[0] != 0:
        raise ValueError('Depth cut must keep a contiguous top slice of >= 10 layers')
    out = dict(col)
    for key in ('height_m', 'temperature_K', 'ne_m3', 'nHtot_m3', 'vlos_m_s'):
        out[key] = col[key][idx].copy()
    out['ndep'] = idx.size
    return out


def set_departures(ctx, b) -> None:
    """Warm start: initial populations n = b * nStar (b per level and depth), then
    renormalised to the column's own total so particle conservation holds."""
    active = ctx.activeAtoms[0]
    n = np.asarray(b) * active.nStar
    n *= active.nTotal / n.sum(axis=0)
    if not (np.isfinite(n).all() and np.all(n > 0)):
        raise ValueError('Invalid warm-start populations')
    active.n[:] = n


def resample_logtau(col, logtau, grid) -> dict:
    """The column on a fixed log tau500 grid (top first): height and T linear,
    ne and nHtot log-linear, vlos linear in log tau. An interpolation of the snapshot,
    declared as such; its effect on the observable is measured against the native grid."""
    logtau = np.asarray(logtau)
    grid = np.asarray(grid, dtype=float)
    if not np.all(np.diff(logtau) > 0):
        raise ValueError('Column log tau500 is not monotonic; cannot resample')
    if grid[0] < logtau[0] or grid[-1] > logtau[-1] or not np.all(np.diff(grid) > 0):
        raise ValueError('Resampling grid outside the column or not increasing')
    f = lambda y: np.interp(grid, logtau, y)
    out = dict(col, ndep=grid.size, logtau500=grid.copy(),
               height_m=f(col['height_m']), temperature_K=f(col['temperature_K']),
               ne_m3=np.exp(f(np.log(col['ne_m3']))), nHtot_m3=np.exp(f(np.log(col['nHtot_m3']))),
               vlos_m_s=f(col['vlos_m_s']))
    if not np.all(np.diff(out['height_m']) < 0):
        raise ValueError('Resampled heights not top-first decreasing')
    return out


#: Disk-centre O I 777 triplet equivalent-width windows (air A), one per line.
EW_WINDOWS_AIR_A = ((7771.2, 7773.0), (7773.0, 7774.8), (7774.8, 7776.5))


def ew_monitor(wave_air_A, wave_vac_nm, windows=EW_WINDOWS_AIR_A):
    """Monitor for `solve_column`: the mu = 1 equivalent widths (mA) in `windows`,
    against a straight continuum through the window edges."""
    wave_air_A = np.asarray(wave_air_A)

    def ew(ctx):
        intensity = np.asarray(ctx.compute_rays(wavelengths=wave_vac_nm, mus=1.0)).ravel()
        r = 1 - intensity / np.interp(wave_air_A, [wave_air_A[0], wave_air_A[-1]], [intensity[0], intensity[-1]])
        return [float(np.trapezoid(r[(wave_air_A >= lo) & (wave_air_A < hi)],
                                   wave_air_A[(wave_air_A >= lo) & (wave_air_A < hi)]) * 1e3)
                for lo, hi in windows]
    return ew
