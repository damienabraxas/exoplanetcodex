"""Experimental trace LTE line opacity added to native LW continuum.

Natural Voigt broadening only; pressure broadening, isotopic splitting and
molecular opacity remain outside this prototype. No pseudo atom/partition fit.
"""
import numpy as np


def lte_line_coefficients(wavelength_nm, temperature, n_ion_over_partition,
                          vturb_m_s, *, lambda0_nm, mass_amu, excitation_eV,
                          loggf, gamma_s_inv):
    """Return extinction m^-1 and thermal emissivity per Hz/sr in SI."""
    from scipy.special import wofz
    import lightweaver.constants as c
    w=np.asarray(wavelength_nm);T=np.asarray(temperature)
    if (np.any(w<=0) or np.any(T<=0) or np.any(n_ion_over_partition<0)
            or gamma_s_inv<0 or mass_amu<=0 or not np.isfinite(loggf)):
        raise ValueError('Invalid LTE line inputs')
    nu=c.CLight/(w[:,None]*1e-9);nu0=c.CLight/(lambda0_nm*1e-9)
    doppler=nu0/c.CLight*np.sqrt(2*c.KBoltzmann*T/(mass_amu*c.Amu)+np.asarray(vturb_m_s)**2)
    profile=wofz((nu-nu0)/doppler+1j*gamma_s_inv/(4*np.pi*doppler)).real/(np.sqrt(np.pi)*doppler)
    lower_times_f=n_ion_over_partition*10**loggf*np.exp(-excitation_eV*c.EV/(c.KBoltzmann*T))
    cross_section=c.QElectron**2/(4*c.Epsilon0*c.MElectron*c.CLight)
    exponent=c.HPlanck*nu/(c.KBoltzmann*T)
    chi=cross_section*lower_times_f*profile*(-np.expm1(-exponent))
    planck=2*c.HPlanck*nu**3/c.CLight**2/np.expm1(exponent)
    return chi,chi*planck


def ni_provider(line, abundance_dex):
    """Build a provider retaining native background and adding trace Ni opacity.

    LTE Saha ionization and U(T) use the installed Kurucz table; fixed atmosphere
    ne/nH and background are retained. This does not solve a new atmosphere/EOS.
    """
    from dataclasses import replace
    from lightweaver.LwCompiled import BackgroundProvider,FastBackground
    from lightweaver.atomic_table import KuruczPfTable,PeriodicTable

    class NiBackground(BackgroundProvider):
        def __init__(self,eqpops,radset,wavelength):
            self.native=FastBackground(eqpops,radset,wavelength,Nthreads=1)
            self.wavelength=np.asarray(wavelength)
            partition=KuruczPfTable(atomicAbundance=radset.abundance)[PeriodicTable['Ni']]
            self.partition=replace(partition,abundance=10**(abundance_dex-12))

        def compute_background(self,atmos,chi,eta,sca):
            self.native.compute_background(atmos,chi,eta,sca)
            a=atmos.pyAtmos
            if a.temperature.min()<self.partition.Tpf[0] or a.temperature.max()>self.partition.Tpf[-1]:
                raise ValueError('Outside Ni partition-function temperature domain')
            ions=self.partition.lte_ionisation(a)
            U=np.exp(np.interp(a.temperature,self.partition.Tpf,self.partition.pf[0]))
            extinction,emission=lte_line_coefficients(self.wavelength,a.temperature,ions[0]/U,a.vturb,
                lambda0_nm=line['lambda0_vac_nm'],mass_amu=PeriodicTable['Ni'].mass,
                excitation_eV=line['excitation_eV'],loggf=line['loggf'],gamma_s_inv=line['gamma_rad_s_inv'])
            chi[:]+=extinction;eta[:]+=emission
    return NiBackground
