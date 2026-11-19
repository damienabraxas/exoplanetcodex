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

    Multi3D z is DEPTH (increases inward) and vz is positive along +z, i.e. downward.
    Lightweaver's 1D geometric scale is HEIGHT (decreasing with index, top first) and
    its vlos is positive along +height (upward). Both flips happen here and nowhere else.
    """
    z = np.asarray(cube['z'], dtype=float)
    if not np.all(np.diff(z) > 0):
        raise ValueError('Multi3D z must increase inward (depth); refusing an unknown orientation')
    sl = (slice(None), int(iy), int(ix))
    temperature = cube['temperature'][sl].copy()
    if not temperature[-1] > temperature[0]:
        raise ValueError('Column is not hotter at depth; the axis order is wrong')
    n_htot = cube['hydrogen_density'][sl] + cube['proton_density'][sl]
    col = dict(source='stagger', ix=int(ix), iy=int(iy), ndep=z.size,
               height_m=-z.copy(), temperature_K=temperature,
               ne_m3=cube['electron_density'][sl].copy(), nHtot_m3=n_htot.copy(),
               vlos_m_s=-cube['vz'][sl].copy())
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


def solve_column(ctx, *, rtol: float = 1e-6, maxiter: int = 1500, seconds: float = 900.):
    """Joint RT + SE iteration, the G4 convergence rule: population change, normalised
    rate residual and J change all below rtol (after >= 3 iterations). Returns the
    populations and the full history. Raises on non-convergence; nothing is clipped."""
    from pipeline.bride_numerical_validation import physical_rate_matrix
    from pipeline.bride_parameter_checks import positive_stat_equil
    active = ctx.activeAtoms[0]
    start = time.perf_counter()
    history = []
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
        rejected = positive_stat_equil(ctx, active, old, reject_invalid_ng=False)  # G4 Mac baseline
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
        if iteration >= 2 and max(delta, rate, float(update.dJMax)) < rtol:
            return n, history, time.perf_counter() - start
    raise RuntimeError('Joint convergence not reached')


def converge_continuum(ctx, *, tol: float = 1e-8, maxiter: int = 1000) -> int:
    for j in range(maxiter):
        if ctx.formal_sol_gamma_matrices().dJMax < tol:
            return j + 1
    raise RuntimeError('Continuum J convergence failed')
