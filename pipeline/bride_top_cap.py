"""Experimental uniform-source cap above a finite native column-mass boundary.

This is a static 1D source-family diagnostic, not an atmosphere extension or an
adopted reconstruction of the unknown grid producer. No free thickness factor.
"""
import numpy as np
import time
from lightweaver import BoundaryCondition


def cap_intensity(chi_m_inv, eta_per_m, length_m, mus):
    """Exact slab solution with zero external illumination and constant chi/eta."""
    chi=np.asarray(chi_m_inv,dtype=float);eta=np.asarray(eta_per_m,dtype=float)
    mu=np.asarray(mus,dtype=float)
    if (chi.ndim!=1 or eta.shape!=chi.shape or mu.ndim!=1 or not mu.size
            or not np.isfinite(chi).all() or not np.isfinite(eta).all()
            or np.any(chi<0) or np.any(eta<0) or not np.isfinite(length_m)
            or length_m<0 or not np.isfinite(mu).all() or np.any(mu<=0) or np.any(mu>1)):
        raise ValueError('Cap requires finite nonnegative coefficients/length and 0<mu<=1')
    distance=length_m/mu
    attenuation=-np.expm1(-chi[:,None]*distance[None,:])
    factor=np.broadcast_to(distance,attenuation.shape).copy()
    np.divide(attenuation,chi[:,None],out=factor,where=chi[:,None]>0)
    return eta[:,None]*factor


class FiniteColumnCap(BoundaryCondition):
    """Frozen coefficient snapshot, explicitly refreshed by the outer iteration.

    Contains no Context pointer: compute_rays can safely deepcopy it and sample
    new angles/wavelengths. Snapshot interpolation is linear and in-domain only.
    """
    def __init__(self, length_m):
        if not np.isfinite(length_m) or length_m<0:
            raise ValueError('Invalid native-column cap length')
        self.length_m=float(length_m)
        self.wave=None
        self.chi=None
        self.eta=None

    def compute_bc(self, atmos, spect):
        wave=np.asarray(spect.wavelength)
        if self.wave is None:
            return np.zeros((len(wave),len(self.muz),1))
        if wave.min()<self.wave[0] or wave.max()>self.wave[-1]:
            raise ValueError('Cap snapshot cannot extrapolate wavelengths')
        chi=np.interp(wave,self.wave,self.chi)
        eta=np.interp(wave,self.wave,self.eta)
        return np.ascontiguousarray(cap_intensity(chi,eta,self.length_m,abs(self.muz))[:,:,None])

    def update_from_context(self, ctx, previous_J_top):
        """Store actual top chi and eta+sca*Jdag used by this formal solution."""
        wave=np.asarray(ctx.spect.wavelength)
        chi=np.asarray(ctx.depthData.chi)[:,:,:,0]
        eta=np.asarray(ctx.depthData.eta)[:,:,:,0]
        if chi.ndim!=3 or chi.shape[0]!=len(wave):
            raise ValueError('Wrong depth coefficient shape')
        if not np.any(chi>0):
            raise ValueError('No top-depth coefficients captured; use the full formal-solution path')
        # Static CRD diagnostic: all angles/directions must give the same chi/eta.
        np.testing.assert_allclose(chi,np.broadcast_to(chi[:,0,0,None,None],chi.shape),rtol=1e-12,atol=0)
        np.testing.assert_allclose(eta,np.broadcast_to(eta[:,0,0,None,None],eta.shape),rtol=1e-12,atol=0)
        new_chi=chi[:,0,0].copy()
        new_eta=eta[:,0,0]+np.asarray(ctx.background.sca)[:,0]*previous_J_top
        mu=abs(np.asarray(ctx.atmos.pyAtmos.muz))
        new=cap_intensity(new_chi,new_eta,self.length_m,mu)
        old=np.zeros_like(new) if self.wave is None else cap_intensity(self.chi,self.eta,self.length_m,mu)
        delta=float(np.max(abs(new-old)/np.maximum(abs(new),np.finfo(float).tiny)))
        self.wave=wave.copy();self.chi=new_chi;self.eta=new_eta.copy()
        return delta


def cap_from_structure(structure):
    """Native m/rho sets an opacity-per-mass cap; refuse invented columns."""
    if 'column_mass_kg_m2' not in structure:
        raise ValueError('Native column mass is required for the finite-column cap')
    column=np.asarray(structure['column_mass_kg_m2'])
    rho=np.asarray(structure['rho_kg_m3'])
    if (column.shape!=rho.shape or not np.isfinite(column).all()
            or np.any(column<=0) or np.any(np.diff(column)<=0)
            or not np.isfinite(rho).all() or np.any(rho<=0)):
        raise ValueError('Invalid native column-mass/density metadata')
    return FiniteColumnCap(column[0]/rho[0])


def emergent_rays(ctx, wavelengths, mus, weights):
    """Retain grid-point rays and transfer through the same cap to its outer face.

    Capture coefficients on the actual diagnostic wavelength/ray grid, avoiding
    interpolation of sharp line opacity for the emergent cap correction.
    """
    _, rayctx=ctx.compute_rays(wavelengths=wavelengths,mus=mus,returnCtx=True,upOnly=False)
    weights=np.asarray(weights,dtype=float)
    if (weights.shape!=np.asarray(mus).shape or not np.isfinite(weights).all()
            or np.any(weights<=0) or not np.isclose(weights.sum(),1.,rtol=1e-14,atol=0)):
        raise ValueError('Explicit normalized angular quadrature weights required')
    # Atmosphere.rays is for ray tracing and leaves J-quadrature weights unset.
    # Full coefficient capture also updates J, so supply the actual flux quadrature.
    rayctx.atmos.pyAtmos.wmu[:]=weights
    rayctx.depthData.fill=True
    previous_J_top=np.array(rayctx.spect.J)[:,0].copy()
    # simple formal_sol intentionally disables StoreDepthData in LW0.16.1.
    # This full formal solve updates rates/J, but never updates populations.
    rayctx.formal_sol_gamma_matrices(lambdaIterate=True)
    if not np.isfinite(np.asarray(rayctx.spect.J)).all():
        raise ValueError('Nonfinite diagnostic-ray radiation field')
    grid=np.array(rayctx.spect.I)
    cap=rayctx.atmos.pyAtmos.zUpperBc
    chi=np.asarray(rayctx.depthData.chi)[:,0,1,0]
    eta=np.asarray(rayctx.depthData.eta)[:,0,1,0]+np.asarray(rayctx.background.sca)[:,0]*previous_J_top
    outgoing=grid*np.exp(-chi[:,None]*cap.length_m/np.asarray(mus)[None,:])
    outgoing+=cap_intensity(chi,eta,cap.length_m,mus)
    return grid,outgoing


def converge_cap_radiation(ctx, *, rtol, maxiter, deadline):
    """Solve J and the cap at fixed populations before a population-map step.

    The one-refresh-per-SE-step control is retained separately. This nested
    iteration changes no physical coefficient or boundary equation.
    """
    cap=ctx.atmos.pyAtmos.zUpperBc
    original_J=np.array(ctx.spect.J).copy()
    mu=abs(np.asarray(ctx.atmos.pyAtmos.muz))
    original_boundary=(np.zeros((len(ctx.spect.wavelength),len(mu))) if cap.wave is None
                       else cap_intensity(cap.chi,cap.eta,cap.length_m,mu))
    for index in range(maxiter):
        if time.perf_counter()>deadline:
            raise TimeoutError('Fixed-population cap radiation budget exceeded')
        previous_J_top=np.array(ctx.spect.J)[:,0].copy()
        update=ctx.formal_sol_gamma_matrices(lambdaIterate=True)
        if not np.isfinite(update.dJMax):
            raise ValueError('Nonfinite cap radiation update; verify angular quadrature')
        boundary_change=cap.update_from_context(ctx,previous_J_top)
        if max(update.dJMax,boundary_change)<rtol:
            break
    else:
        raise RuntimeError(f'Fixed-population cap radiation did not converge: J={update.dJMax}, cap={boundary_change}')
    final_J=np.array(ctx.spect.J)
    final_boundary=cap_intensity(cap.chi,cap.eta,cap.length_m,mu)
    relative=lambda new,old:float(np.max(abs(new-old)/np.maximum(abs(new),np.finfo(float).tiny)))
    return dict(iterations=index+1,final_J_change=float(update.dJMax),
                final_boundary_change=boundary_change,
                outer_J_change=relative(final_J,original_J),
                outer_boundary_change=relative(final_boundary,original_boundary))
