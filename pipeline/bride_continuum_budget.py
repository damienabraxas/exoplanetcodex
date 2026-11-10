"""Diagnostic reconstruction of passive bound-free and H/He resonance scattering.

Matches the inspected LW FastBackground source at fixed populations; it never
updates the atmosphere, populations, radiation or backend coefficients.
"""
import numpy as np


def continuum_budget(ctx):
    from lightweaver import constants as c
    from lightweaver.atomic_table import PeriodicTable
    wave=np.asarray(ctx.spect.wavelength)
    radset=ctx.kwargs['spect'].radSet;eq=ctx.kwargs['eqPops']
    atmos=ctx.atmos.pyAtmos
    thermal=-np.expm1(-c.HC/(c.KBoltzmann*1e-9*wave[:,None]*atmos.temperature[None,:]))
    absorption={}
    for model in radset.passiveAtoms:
        pop=eq.atomicPops[model.element].n
        terms=sum((transition.alpha(wave)[:,None]*thermal*pop[transition.i][None,:]
                   for transition in model.continua),np.zeros_like(thermal))
        absorption[str(model.element)]=terms
    sigma=8*np.pi/3*(c.QElectron/(np.sqrt(4*np.pi*c.Epsilon0)*
                        np.sqrt(c.MElectron)*c.CLight))**4
    scattering={'electron_Thomson':np.broadcast_to(atmos.ne[None,:]*sigma,thermal.shape).copy()}
    coefficient=2*np.pi*(c.QElectron/c.Epsilon0)*c.QElectron/c.MElectron/c.CLight
    for element in (PeriodicTable[1],PeriodicTable[2]):
        model=radset.atoms[element];terms=np.zeros_like(thermal)
        ground=eq.atomicPops[element].n[0]
        for line in model.lines:
            if line.i!=0:continue
            mask=wave>line.wavelength()[-1]
            f=line.Aji*(line.jLevel.g/line.iLevel.g)*(line.lambda0*1e-9)**2/coefficient
            terms[mask]+=sigma*f/((wave[mask,None]/line.lambda0)**2-1)**2*ground[None,:]
        scattering[str(element)+'_resonant_Rayleigh']=terms
    return absorption,scattering
