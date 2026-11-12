"""Explicit supported wavelength coordinates; qualified binary readers stay immutable."""
import numpy as np

def supported_air_axis(wavelength_vac_A):
    """AIR coordinates only above200nm; keep vacuum primary for UV/IR output.

    NaN denotes the explicitly unsupported AIR convention, never missing intensity.
    """
    from pipeline.wavelength_util import vac_to_air
    wave=np.asarray(wavelength_vac_A,dtype=float)
    if wave.ndim!=1 or wave.size==0 or not np.isfinite(wave).all()or np.any(wave<=0)or np.any(np.diff(wave)<=0):raise ValueError('Invalid explicit vacuum wavelength axis')
    supported=wave>2000.;air=np.full(wave.shape,np.nan)
    if np.any(supported):air[supported]=vac_to_air(wave[supported])
    return air,supported
