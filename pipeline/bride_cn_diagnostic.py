"""Diagnostic line extinction from externally supplied CN state.

No chemistry, partition, isotope, wavelength-frame or background-provider inference.
Inputs must already carry one consistent population and statistical convention.
"""
import numpy as np
from scipy.special import wofz


def external_state_extinction(offset_doppler, damping_a, doppler_width_Hz,
                              lower_population_m3, gf, stimulated_emission,
                              integrated_strength_m2_Hz):
    """Extinction m^-1: shared 1D or explicit 2D offsets x depth state.

    integrated_strength is the frequency-integrated cross section for gf=1;
    the profile H(a,v)/(sqrt(pi)*Dnu) has unit frequency area. No default
    constants or population conventions are substituted for the caller.
    """
    v=np.asarray(offset_doppler,dtype=float)
    a,width,n,g,stim,strength=[np.asarray(x,dtype=float)for x in
        (damping_a,doppler_width_Hz,lower_population_m3,gf,stimulated_emission,integrated_strength_m2_Hz)]
    arrays=[v,a,width,n,g,stim,strength]
    if (v.ndim not in (1,2)or not len(v)or a.ndim!=1 or not len(a)
            or (v.ndim==2 and v.shape[1]!=len(a))
            or any(x.shape!=a.shape for x in arrays[2:])
            or not all(np.isfinite(x).all()for x in arrays)
            or np.any(a<0)or np.any(width<=0)or np.any(n<0)
            or np.any(g<0)or np.any(stim<0)or np.any(stim>1)
            or np.any(strength<=0)):
        raise ValueError('Invalid or inconsistent external CN state')
    profile=wofz((v[:,None]if v.ndim==1 else v)+1j*a[None,:]).real/(np.sqrt(np.pi)*width[None,:])
    chi=profile*(n*g*stim*strength)[None,:]
    if not np.isfinite(chi).all()or np.any(chi<0):raise ValueError('Invalid external-state extinction')
    return chi


def frequency_offsets(wavelength_A, lambda0_air_A, doppler_width_Hz, c_cm_s,
                      native_lambda_cm, *, input_medium, mode):
    """Explicit legacy-native versus exact coordinate map; no inferred frame.

    Native mode retains the exported rounded XL for its linear slope. Exact
    vacuum mode converts only the declared AIR line centres using the canonical
    owner; the caller supplies an already-vacuum wavelength grid. All return
    red-positive offsets, equivalent for a symmetric Voigt profile.
    """
    from pipeline.wavelength_util import air_to_vac,AIR_VACUUM_BOUNDARY_A
    wave=np.asarray(wavelength_A,dtype=float)
    centre,width,c,legacy=[np.asarray(x,dtype=float)for x in
        (lambda0_air_A,doppler_width_Hz,c_cm_s,native_lambda_cm)]
    if (wave.ndim!=1 or not len(wave)or centre.ndim!=1 or not len(centre)
            or any(x.shape!=centre.shape for x in (width,c,legacy))
            or not all(np.isfinite(x).all()for x in (wave,centre,width,c,legacy))
            or any(np.any(x<=0)for x in (wave,width,c,legacy))
            or np.any(centre<AIR_VACUUM_BOUNDARY_A)):
        raise ValueError('Invalid explicit CN wavelength/state contract')
    if input_medium=='air'and np.any(wave<AIR_VACUUM_BOUNDARY_A):
        raise ValueError('AIR wavelengths undefined below canonical boundary')
    if mode=='native_linear_air':
        if input_medium!='air':raise ValueError('Native mapping requires AIR grid')
        delta=(wave[:,None]-centre[None,:])*1e-8*c[None,:]/legacy[None,:]**2
    elif mode in ('exact_air_frequency','exact_vacuum_frequency'):
        expected='air'if mode=='exact_air_frequency'else 'vacuum'
        if input_medium!=expected:raise ValueError('Declared frame incompatible with mapping')
        rest=centre if expected=='air'else air_to_vac(centre)
        # Stable reciprocal-frequency difference, preserving zero at line centre.
        delta=c[None,:]*1e8*(wave[:,None]-rest[None,:])/(rest[None,:]*wave[:,None])
    else:raise ValueError('Unknown explicit CN frequency mapping')
    offsets=delta/width[None,:]
    if not np.isfinite(offsets).all():raise ValueError('Invalid CN frequency offsets')
    return offsets
