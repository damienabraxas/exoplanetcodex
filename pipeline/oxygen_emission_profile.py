"""Explicit fixed-state emission/bottom-boundary axis; no active-rate changes.

Callable boundary reproduces FormalScalar.cpp's piecewise_linear_1d diffusion
condition using captured total extinction. It never substitutes I=B after the
first coefficient-capture iteration. No compiled constants are patched.
"""
import numpy as np
from lightweaver import BoundaryCondition


def profile_planck(wave, temperature, profile, native_source):
    wave=np.asarray(wave,float);temperature=np.asarray(temperature,float)
    if wave.ndim!=1 or temperature.ndim!=1 or not wave.size or temperature.size<2 or any(not np.isfinite(x).all()or np.any(x<=0)for x in [wave,temperature]):
        raise ValueError('Invalid emission profile grid/state')
    if profile=='source-native':
        from pipeline.hop_packet import native_planck
        B,_=native_planck(wave,temperature,native_source)
    elif profile=='held-LW':
        import lightweaver.constants as c
        # Match the held compiled planck_nu arithmetic and cutoff, including
        # exp(x)-1 rather than replacing it with an expm1 approximation.
        x=(c.HPlanck*c.CLight/(c.KBoltzmann*1e-9)/wave[:,None])/temperature[None,:]
        B=np.zeros_like(x);present=x<=150
        prefactor=np.broadcast_to((2*c.HPlanck*c.CLight/1e-27)/wave[:,None]**3,x.shape)
        B[present]=prefactor[present]/(np.exp(x[present])-1)
    else:raise ValueError('Unknown explicit emission profile')
    if not np.isfinite(B).all()or np.any(B<=0):raise ValueError('Emission profile requires positive finite Planck values')
    return B


def diffusion_intensity(Bbottom, chi_bottom, spacing, mus):
    B=np.asarray(Bbottom,float);chi=np.asarray(chi_bottom,float);mu=np.asarray(mus,float)
    if B.ndim!=2 or B.shape[1]!=2 or chi.shape!=B.shape or mu.ndim!=1 or not mu.size or not np.isfinite(spacing)or spacing<=0 or any(not np.isfinite(x).all()or np.any(x<=0)for x in [B,chi,mu])or np.any(mu>1):
        raise ValueError('Invalid lower diffusion boundary inputs')
    dtau=(.5/mu[None,:])*(chi[:,1]+chi[:,0])[:,None]*spacing
    intensity=B[:,1,None]-(B[:,0]-B[:,1])[:,None]/dtau
    if not np.isfinite(intensity).all()or np.any(intensity<0):raise ValueError('Negative/nonfinite diffusion boundary; no clipping')
    return intensity


class EmissionProfileBoundary(BoundaryCondition):
    """Stationary, fixed-population, exact-grid bottom diffusion snapshot."""
    def __init__(self,profile,native_source):
        if profile not in ['held-LW','source-native']:raise ValueError('Unknown explicit emission profile')
        self.profile=profile;self.native_source=native_source;self.wave=None;self.chi=None

    def compute_bc(self,atmos,spect):
        wave=np.asarray(spect.wavelength)
        B=profile_planck(wave,np.asarray(atmos.temperature)[-2:],self.profile,self.native_source)
        if self.wave is None:return np.ascontiguousarray(np.broadcast_to(B[:,1,None,None],(len(wave),len(self.muz),1)))
        if not np.array_equal(wave,self.wave)or not np.array_equal(atmos.temperature,self.temperature)or not np.array_equal(atmos.height,self.height):raise ValueError('Changed lower boundary grid/state')
        return np.ascontiguousarray(diffusion_intensity(B,self.chi,abs(self.height[-1]-self.height[-2]),abs(self.muz))[:,:,None])

    def update_from_context(self,ctx):
        atmosphere=ctx.atmos.pyAtmos;wave=np.asarray(ctx.spect.wavelength)
        chi=np.asarray(ctx.depthData.chi)[:,:,:,-2:]
        if chi.shape!=(len(wave),len(atmosphere.muz),2,2):raise ValueError('Wrong total bottom opacity capture')
        np.testing.assert_array_equal(chi,np.broadcast_to(chi[:,0,0,None,None,:],chi.shape))
        new=chi[:,0,0,:].copy()
        if self.chi is not None:np.testing.assert_array_equal(new,self.chi)
        B=profile_planck(wave,atmosphere.temperature[-2:],self.profile,self.native_source)
        intensity=diffusion_intensity(B,new,abs(atmosphere.height[-1]-atmosphere.height[-2]),abs(atmosphere.muz))
        old=B[:,1,None]if self.chi is None else intensity
        change=float(np.max(abs(intensity-old)/np.maximum(abs(intensity),np.finfo(float).tiny)))
        self.wave=wave.copy();self.chi=new;self.temperature=atmosphere.temperature.copy();self.height=atmosphere.height.copy()
        return change
