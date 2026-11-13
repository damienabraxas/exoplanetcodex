"""Source-defined linear zero-ghost boundary for fixed-population radiation controls."""
import time
import numpy as np
from pipeline.bride_top_cap import FiniteColumnCap

def linear_ghost_intensity(chi,eta,step_m,mu):
    """Native linear solver: vacuum ghost (chi,S,I=0), full cell spacing to node.

    dtau=chi*step/(2*mu), I_node=S_node*[1-(1-exp(-dtau))/dtau].
    This is the source's discrete ghost interval, not a fitted atmosphere extension.
    """
    chi=np.asarray(chi);eta=np.asarray(eta);mu=np.asarray(mu)
    if chi.ndim!=1 or eta.shape!=chi.shape or mu.ndim!=1 or not np.isfinite(step_m)or step_m<=0 or any(not np.isfinite(x).all()for x in [chi,eta,mu])or np.any(chi<0)or np.any(eta<0)or np.any(mu<=0)or np.any(mu>1)or np.any((chi==0)&(eta!=0)):raise ValueError('Invalid native ghost interval')
    tau=chi[:,None]*step_m/(2*mu[None,:]);f=np.empty_like(tau);small=tau<1e-3;t=tau[small]
    f[small]=t*(.5+t*(-1/6+t*(1/24+t*(-1/120+t/720))))
    f[~small]=1+np.expm1(-tau[~small])/tau[~small]
    S=np.divide(eta,chi,out=np.zeros_like(eta,dtype=float),where=chi>0)
    return S[:,None]*f

class NativeLinearGhost(FiniteColumnCap):
    """Read-only geometry control; no SE preconditioner/adoption qualification."""
    def __init__(self,step_m):
        super().__init__(0.)
        if not np.isfinite(step_m)or step_m<=0:raise ValueError('Invalid source ghost spacing')
        self.step_m=float(step_m)
    def incoming(self,mu):return linear_ghost_intensity(self.chi,self.eta,self.step_m,abs(np.asarray(mu)))
    def compute_bc(self,atmos,spect):
        if self.wave is None:return np.zeros((len(spect.wavelength),len(self.muz),1))
        if not np.array_equal(spect.wavelength,self.wave):raise ValueError('Ghost control refuses changed wavelength grid')
        return np.ascontiguousarray(self.incoming(self.muz)[:,:,None])
    def update_from_context(self,ctx,previous_J_top):
        old=None if self.wave is None else self.incoming(ctx.atmos.pyAtmos.muz)
        super().update_from_context(ctx,previous_J_top) # capture/check chi, eta+sca*J only
        new=self.incoming(ctx.atmos.pyAtmos.muz)
        if old is None:old=np.zeros_like(new)
        return float(np.max(abs(new-old)/np.maximum(abs(new),np.finfo(float).tiny)))

def converge_ghost_radiation(ctx,*,rtol,maxiter,deadline):
    boundary=ctx.atmos.pyAtmos.zUpperBc
    if not isinstance(boundary,NativeLinearGhost):raise ValueError('Require explicit native ghost control')
    for i in range(maxiter):
        if time.perf_counter()>deadline:raise TimeoutError('Original fixed-radiation budget exhausted')
        top=np.array(ctx.spect.J)[:,0].copy();update=ctx.formal_sol_gamma_matrices(lambdaIterate=True);change=boundary.update_from_context(ctx,top)
        if not np.isfinite(update.dJMax):raise ValueError('Nonfinite ghost radiation')
        if max(update.dJMax,change)<rtol:return dict(iterations=i+1,final_J_change=float(update.dJMax),final_boundary_change=change,ghost_step_m=boundary.step_m)
    raise RuntimeError('Original fixed-radiation maxiter exhausted')
