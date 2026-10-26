"""Analytic controls through the pinned LW formal solver, synthetic data only.

The continuum-only atmosphere/spectrum control has no active atomic models.
A BackgroundProvider supplies independently specified coefficients directly to the
actual compiled transfer backend. Separate synthetic CE/CI and line controls
exercise SE and thermal equilibrium, without claiming an admitted oxygen model.
"""
from __future__ import annotations
from types import SimpleNamespace
import numpy as np


def constant_slab(*, tau_vertical, source, lower_intensity, upper_intensity,
                  mu, depths, solver):
    """Return actual emergent intensity and analytic homogeneous-slab prediction.

    Unit-length SI slab; extinction is tau_vertical m^-1, emission chi*S.
    No scattering. Nonzero opacity is required because the backend constructs
    source eta/chi; exact vacuum has a separate undefined-source limitation.
    Numerical test tolerance is specified by the caller before evaluation.
    """
    import lightweaver as lw
    from lightweaver.LwCompiled import BackgroundProvider
    if not (np.isfinite(tau_vertical) and tau_vertical > 0 and 0 < mu <= 1
            and depths >= 3):
        raise ValueError('Require positive optical depth, valid ray and >=3 depths')

    class Fixed(lw.BoundaryCondition):
        def __init__(self, intensity): self.intensity = intensity
        def compute_bc(self, atmos, spect):
            return np.full((len(spect.wavelength), len(self.muz), 1), self.intensity)

    class Coefficients(BackgroundProvider):
        def compute_background(self, atmos, chi, eta, sca):
            chi[:] = tau_vertical
            eta[:] = tau_vertical * source
            sca[:] = 0.

    n = depths
    atmosphere = lw.Atmosphere.make_1d(scale=lw.ScaleType.Geometric,
        depthScale=np.linspace(1., 0., n), temperature=np.full(n, 6000.),
        ne=np.full(n, 1e16), nHTot=np.full(n, 1e20), vturb=np.zeros(n),
        vlos=np.zeros(n), logG=2., lowerBc=Fixed(lower_intensity),
        upperBc=Fixed(upper_intensity))
    atmosphere.rays(muz=np.array([mu]))
    spect = SimpleNamespace(wavelength=np.array([500., 600.]),
        radSet=SimpleNamespace(activeAtoms=[], detailedAtoms=[]))
    ctx = lw.Context(atmosphere, spect, {}, backgroundProvider=Coefficients,
        Nthreads=1, formalSolver=solver)
    ctx.formal_sol(upOnly=False)
    actual = np.array(ctx.spect.I, copy=True)
    t = tau_vertical / mu
    expected = lower_intensity * np.exp(-t) + source * (-np.expm1(-t))
    return actual, float(expected)


def homogeneous_error_bound(*, tau_vertical, mu, depths, source,
                            lower_intensity, solver):
    """Conservative forward bound for the inspected 0.16.1 coefficient branches.

    Uniform slab only. This is not a physical-error/dex acceptance threshold.
    LwInternal.hpp w2 truncates w0 after t^2 below 5e-4; Besser edt is
    degree 7 below .14; Bezier edt and summed weights are degree 3 below .05.
    Alternating exponential-series remainder <= t^(p+1)/(p+1)! per cell.
    Recurrence attenuation <=1 permits summing per-cell bounds. Away from the
    Taylor branches rational cancellation amplifies roundoff by t^-p. 128 eps
    per cell bounds the inspected bounded arithmetic and coefficient sums.
    """
    import math
    t = tau_vertical / mu / (depths - 1)
    n = depths - 1
    scale = max(abs(source), abs(lower_intensity), 1.)
    eps = np.finfo(np.float64).eps
    if solver == 'piecewise_linear_1d':
        p, cutoff = 2, 5e-4
    elif solver == 'piecewise_besser_1d':
        p, cutoff = 7, .14
    elif solver == 'piecewise_bezier3_1d':
        p, cutoff = 3, .05
    else:
        raise ValueError('No inspected error model for solver')
    if t < cutoff:
        truncation = 2 * scale * t**(p+1) / math.factorial(p+1)
        rounding = 128 * eps * scale
    else:
        # Rational coefficient expressions have denominators t^2 or t^3.
        denominator_power = {'piecewise_linear_1d': 0,
                             'piecewise_besser_1d': 2,
                             'piecewise_bezier3_1d': 3}[solver]
        rounding = 128 * eps * scale * (1 + t + t*t + t**3) / min(1., t**denominator_power)
        truncation = 2*scale*np.exp(-t) if t > 30 else 0.
    return n * (truncation + rounding)


def collision_lte_control(*, temperature, initial_weights, rate_scale=1., radiative=False):
    """Actual CE/CI and compiled SE kernel vs independent Saha-Boltzmann law.

    Three artificial levels/rates form a numerical fixture, NOT an oxygen atom.
    With radiative=False there is no radiation; radiative=True adds a synthetic
    CRD line inside an isothermal blackbody enclosure. Both have exact LTE.
    Passive hydrogen supplies LW's required EOS storage; it is not fitted.
    """
    import lightweaver as lw
    from lightweaver.atomic_model import AtomicModel, AtomicLevel
    from lightweaver.atomic_table import Element
    from lightweaver.collisional_rates import CE, CI
    from lightweaver.rh_atoms import H_6_atom
    # Match the pinned engine's declared constant vintage; modern CODATA is a
    # different physical setup. Only constants are shared, never its LTE routine.
    import lightweaver.constants as constants
    h, k, c, m_e = (constants.HPlanck, constants.KBoltzmann,
                    constants.CLight, constants.MElectron)

    n = 5
    t = float(temperature)
    atmosphere = lw.Atmosphere.make_1d(scale=lw.ScaleType.Geometric,
        depthScale=np.linspace(4., 0., n), temperature=np.full(n, t),
        ne=np.full(n, 1e16), nHTot=np.full(n, 1e20), vlos=np.zeros(n),
        vturb=np.zeros(n), logG=2., lowerBc=lw.ThermalisedRadiation() if radiative else lw.ZeroRadiation(),
        upperBc=lw.ThermalisedRadiation() if radiative else lw.ZeroRadiation())
    atmosphere.quadrature(2)
    levels = [AtomicLevel(E=0., g=1., label='synthetic0', stage=0),
        AtomicLevel(E=1e4, g=3., label='synthetic1', stage=0),
        AtomicLevel(E=1.1e5, g=1., label='syntheticIon', stage=1)]
    lines = []
    if radiative:
        from lightweaver.atomic_model import VoigtLine, LineType, LinearCoreExpWings
        from lightweaver.broadening import LineBroadening, RadiativeBroadening
        lines = [VoigtLine(j=1, i=0, f=.1, type=LineType.CRD,
            quadrature=LinearCoreExpWings(qCore=3.,qWing=10.,Nlambda=21),
            broadening=LineBroadening(natural=[RadiativeBroadening(1e7)],elastic=[]))]
    model = AtomicModel(element=Element(8), levels=levels, lines=lines, continua=[],
        collisions=[CE(j=1, i=0, temperature=[3000.,10000.],
                       rates=[1e-15*rate_scale]*2),
                    CI(j=2, i=0, temperature=[3000.,10000.],
                       rates=[1e-15*rate_scale]*2)])
    radset = lw.RadiativeSet([H_6_atom(), model])
    radset.set_active('O')
    eq = radset.compute_eq_pops(atmosphere)
    spect = radset.compute_wavelength_grid(extraWavelengths=np.array([500.,600.]))
    ctx = lw.Context(atmosphere, spect, eq, Nthreads=1,formalSolver='piecewise_linear_1d')
    atom = ctx.activeAtoms[0]
    atom.compute_collisions(fillDiagonal=True)
    matrix = atom.C.copy()
    # Equilibrate equations (rows), preserving the same homogeneous SE system.
    # Without this, eliminating the initially dominant population can amplify
    # conditioning errors. Diagnostics always use the unscaled physical matrix.
    row_norm = np.max(np.abs(matrix),axis=1,keepdims=True)
    if np.any(row_norm == 0):
        raise ValueError('Disconnected synthetic SE equation')
    atom.Gamma[:] = matrix / row_norm
    start = np.asarray(initial_weights, dtype=float)
    if start.shape != (3,) or np.any(start <= 0):
        raise ValueError('Three positive initial weights required')
    atom.n[:] = atom.nTotal[None, :] * (start/start.sum())[:, None]
    history = []
    for iteration in range(80 if radiative else 20):
        previous = atom.n.copy()
        if radiative:
            radiation_update = ctx.formal_sol_gamma_matrices()
            matrix = atom.Gamma.copy()  # MALI iteration operator, not physical R+C
            row_norm = np.max(np.abs(matrix),axis=1,keepdims=True)
            atom.Gamma[:] = matrix/row_norm
        ctx.stat_equil()
        actual = atom.n.copy()
        delta = np.max(np.abs(actual-previous)/np.maximum(np.abs(actual),np.finfo(float).tiny))
        if radiative:
            # Gamma is a MALI iteration operator, NOT the physical rate matrix.
            # Reconstruct actual R+C from separately accumulated transition rates.
            physical = atom.C.copy()
            for level in range(physical.shape[0]): physical[level,level,:] = 0.
            for transition in atom.trans:
                physical[transition.j,transition.i,:] += transition.Rij
                physical[transition.i,transition.j,:] += transition.Rji
            for level in range(physical.shape[0]):
                physical[level,level,:] = -np.sum(physical[:,level,:],axis=0)
        else:
            physical = matrix
        balance = np.einsum('ijd,jd->id',physical,actual)
        throughput = np.einsum('ijd,jd->id',np.abs(physical),np.abs(actual))
        relative_residual = np.max(np.abs(balance)/np.maximum(throughput,np.finfo(float).tiny))
        norm_scale = np.max(np.sum(np.abs(physical),axis=1),axis=0)*np.max(np.abs(actual),axis=0)
        backward_error = np.max(np.max(np.abs(balance),axis=0)/norm_scale)
        history.append(dict(iteration=iteration+1,population_change=float(delta),
                            normalized_rate_residual=float(relative_residual),
                            normwise_backward_error=float(backward_error)))
        if delta < 1e-11 and backward_error < (1e-11 if radiative else 4096*np.finfo(float).eps) and (not radiative or radiation_update.dJMax < 1e-11):
            break
    else:
        raise RuntimeError(f'Synthetic SE failed convergence: {history[-3:]}')
    # Independent analytic partition: energies are cm^-1; factor 2 counts free
    # electron spin. Neutral lower/upper degeneracies are 1/3, ionic g=1.
    debye_length = np.sqrt(constants.Epsilon0*k*t /
                           (2.*1e16*constants.QElectron**2))
    ion_lowering = constants.QElectron**2 / (4.*np.pi*constants.Epsilon0*debye_length)
    weights = np.array([1., 3.*np.exp(-h*c*100*1e4/(k*t)),
        2.*(2.*np.pi*m_e*k*t/h**2)**1.5/1e16
          *np.exp(-(h*c*100*1.1e5-ion_lowering)/(k*t))])
    fraction = weights/weights.sum()
    residual = np.einsum('ijd,jd->id',physical,actual)
    throughput = np.einsum('ijd,jd->id',np.abs(physical),actual)
    augmented = (matrix/row_norm)[:,:,0].copy()
    augmented[np.argmax(actual[:,0]),:] = 1.
    condition = np.linalg.cond(augmented,np.inf)
    return dict(radiative=radiative,condition_inf=float(condition),
        normwise_backward_error=float(backward_error),actual_fraction=actual/atom.nTotal[None,:],
        expected_fraction=fraction, conservation=actual.sum(axis=0)/atom.nTotal-1,
        normalized_rate_residual=np.abs(residual)/np.maximum(throughput,np.finfo(float).tiny),
        collision_matrix=matrix, physical_rate_matrix=physical, history=history)
