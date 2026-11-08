"""Diagnostic line extinction from externally supplied CN state.

No chemistry, partition, isotope, wavelength-frame or background-provider inference.
Inputs must already carry one consistent population and statistical convention.
"""
import numpy as np
from scipy.special import wofz


def external_state_extinction(offset_doppler, damping_a, doppler_width_Hz,
                              lower_population_m3, gf, stimulated_emission,
                              integrated_strength_m2_Hz):
    """Extinction m^-1: explicit dimensionless offsets x depth external state.

    integrated_strength is the frequency-integrated cross section for gf=1;
    the profile H(a,v)/(sqrt(pi)*Dnu) has unit frequency area. No default
    constants or population conventions are substituted for the caller.
    """
    v=np.asarray(offset_doppler,dtype=float)
    a,width,n,g,stim,strength=[np.asarray(x,dtype=float)for x in
        (damping_a,doppler_width_Hz,lower_population_m3,gf,stimulated_emission,integrated_strength_m2_Hz)]
    arrays=[v,a,width,n,g,stim,strength]
    if (v.ndim!=1 or not len(v)or a.ndim!=1 or not len(a)
            or any(x.shape!=a.shape for x in arrays[2:])
            or not all(np.isfinite(x).all()for x in arrays)
            or np.any(a<0)or np.any(width<=0)or np.any(n<0)
            or np.any(g<0)or np.any(stim<0)or np.any(stim>1)
            or np.any(strength<=0)):
        raise ValueError('Invalid or inconsistent external CN state')
    profile=wofz(v[:,None]+1j*a[None,:]).real/(np.sqrt(np.pi)*width[None,:])
    chi=profile*(n*g*stim*strength)[None,:]
    if not np.isfinite(chi).all()or np.any(chi<0):raise ValueError('Invalid external-state extinction')
    return chi
