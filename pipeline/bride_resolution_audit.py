"""Explicit Gaussian sensitivity controls, not calibrated FTS kernels."""
import numpy as np


def extra_gaussian_sigma_pixels(wave,centre,target_R,*,native_R=None,no_extra=False):
    w=np.asarray(wave,float)
    if w.ndim!=1 or len(w)<3 or not np.isfinite(w).all() or np.any(np.diff(w)<=0):
        raise ValueError('Invalid increasing reference grid')
    if not np.isfinite(centre) or centre<=0 or not np.isfinite(target_R) or target_R<=0:
        raise ValueError('Invalid reference resolution')
    if no_extra:
        if native_R is not None:raise ValueError('No-extra control cannot specify native R')
        return 0.
    inverse_square=1/target_R**2
    if native_R is not None:
        if not np.isfinite(native_R) or native_R<=target_R:
            raise ValueError('Native R must exceed target R; no deconvolution')
        inverse_square-=1/native_R**2
    # Same centre, median AIR step, FWHM convention as the shipped atlas helper.
    fwhm=centre/target_R if native_R is None else centre*np.sqrt(inverse_square)
    return fwhm/2.3548/float(np.median(np.diff(w)))


def broaden_reference(wave,flux,centre,target_R,*,native_R=None,no_extra=False):
    from scipy.ndimage import gaussian_filter1d
    w,f=np.asarray(wave,float),np.asarray(flux,float)
    if f.shape!=w.shape or not np.isfinite(f).all():raise ValueError('Invalid reference flux')
    sigma=extra_gaussian_sigma_pixels(w,centre,target_R,native_R=native_R,no_extra=no_extra)
    return (f.copy() if sigma==0 else gaussian_filter1d(f,sigma)),dict(
        sigma_pixels=sigma,truncate=4.,boundary_mode='reflect',native_R_assumed=native_R,
        target_R=target_R,no_extra=no_extra,
        scope='Gaussian approximation only; native and target FWHM assumed Gaussian. No calibrated FTS LSF or resolution fit.')
