"""Observable transformations only: fixed-population experiments, no abundance fit."""
import numpy as np


def flux_quadrature(rays):
    """Gauss-Legendre nodes/weights on outward mu in [0,1]."""
    if rays<2:raise ValueError('At least two flux rays required')
    mu,weight=np.polynomial.legendre.leggauss(rays)
    return (mu+1)/2,weight/2


def integrate_flux(intensity,mu,weight):
    """F_nu = 2 pi integral I_nu(mu) mu dmu; last axis is outward mu."""
    intensity=np.asarray(intensity)
    if intensity.shape[-1]!=len(mu) or len(mu)!=len(weight):
        raise ValueError('Intensity/angle shape mismatch')
    if not np.isfinite(intensity).all() or np.any(intensity<=0):
        raise ValueError('Invalid intensities')
    return 2*np.pi*np.sum(intensity*(mu*weight),axis=-1)


def window_integral(wave,normalized,lo,hi,*,max_gap=None):
    """Bounded residual-flux area, in wavelength units; no extrapolated wings."""
    wave=np.asarray(wave);normalized=np.asarray(normalized)
    if wave.shape!=normalized.shape or not np.all(np.diff(wave)>0):
        raise ValueError('Invalid spectral grid')
    if not np.isfinite(normalized).all() or lo< wave[0] or hi>wave[-1] or not lo<hi:
        raise ValueError('Invalid or uncovered integration window')
    if max_gap is not None:
        if not np.isfinite(max_gap) or max_gap<=0:raise ValueError('Invalid maximum spectral gap')
        # Include the samples bracketing both endpoints: inserted endpoints
        # cannot disguise interpolation across an absent output window.
        first=max(0,np.searchsorted(wave,lo,side='right')-1)
        last=min(len(wave)-1,np.searchsorted(wave,hi,side='left'))
        if np.any(np.diff(wave[first:last+1])>max_gap):raise ValueError('Uncovered spectral gap in integration window')
    inside=(wave>lo)&(wave<hi)
    x=np.r_[lo,wave[inside],hi]
    return float(np.trapezoid(1-np.interp(x,wave,normalized),x))
