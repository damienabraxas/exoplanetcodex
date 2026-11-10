"""Opt-in CN provider for disposable H-only synthetic array tests, no RT admission."""
from dataclasses import dataclass
from functools import partial
import numpy as np
from lightweaver.LwCompiled import BackgroundProvider,FastBackground

BIND_FIELDS=('height','temperature','ne','nHTot','vturb','vlos')


def binding(atmosphere,abundance):
    return tuple((key,tuple(np.asarray(getattr(atmosphere,key),dtype=float)))for key in BIND_FIELDS),tuple(float(abundance[z])for z in range(1,93))


@dataclass(frozen=True)
class SyntheticCnPacket:
    purpose:str
    state_binding:tuple
    abundance_binding:tuple
    line_id:str
    lambda0_air_A:float
    lower_population_m3:tuple
    width_Hz:tuple
    damping_a:tuple
    stimulated_emission:tuple
    integrated_strength_m2_Hz:tuple
    gf:float

    def __post_init__(self):
        if self.purpose!='synthetic_array_test':raise ValueError('Only disposable synthetic CN packets admitted by this adapter')
        if not self.line_id:raise ValueError('Synthetic CN line identity required')
        for name in ('lower_population_m3','width_Hz','damping_a','stimulated_emission','integrated_strength_m2_Hz'):
            object.__setattr__(self,name,tuple(float(x)for x in getattr(self,name)))
        from pipeline.bride_cn_diagnostic import external_state_extinction
        external_state_extinction([0.],self.damping_a,self.width_Hz,self.lower_population_m3,np.full(len(self.width_Hz),self.gf),self.stimulated_emission,self.integrated_strength_m2_Hz)
        from pipeline.wavelength_util import AIR_VACUUM_BOUNDARY_A
        if not np.isfinite(self.lambda0_air_A)or self.lambda0_air_A<AIR_VACUUM_BOUNDARY_A:raise ValueError('Invalid AIR line centre')


class EqLookupAudit:
    def __init__(self,wrapped):self.wrapped=wrapped;self.requests=[]
    def __getattr__(self,name):return getattr(object.__getattribute__(self,'wrapped'),name)
    def __contains__(self,key):
        self.requests.append(('contains',str(key)))
        if key=='CN':raise ValueError('Baseline requested CN; duplicate-opacity inventory must be revisited')
        return key in self.wrapped
    def __getitem__(self,key):
        self.requests.append(('getitem',str(key)))
        if key=='CN':raise ValueError('Baseline consumed CN; duplicate-opacity inventory must be revisited')
        return self.wrapped[key]


def synthetic_cn_factory(packet,*,enabled,cn_already_present,support_vac_nm,mapping):
    if type(enabled)is not bool or cn_already_present is not False:raise ValueError('Explicit enable and absent-CN inventory required')
    if not isinstance(packet,SyntheticCnPacket):raise ValueError('Synthetic packet required')
    if mapping!='exact_vacuum_frequency':raise ValueError('Array adapter requires explicit exact vacuum mapping')
    support=tuple(support_vac_nm)
    if len(support)!=2 or not np.isfinite(support).all()or not 0<support[0]<support[1]:raise ValueError('Invalid CN support')
    from pipeline.wavelength_util import air_to_vac
    centre=float(air_to_vac([packet.lambda0_air_A])[0])/10
    if not support[0]<centre<support[1]:raise ValueError('CN centre outside explicit support')
    return partial(SyntheticCnBackground,packet=packet,enabled=enabled,support=support)


class SyntheticCnBackground(BackgroundProvider):
    def __init__(self,eqpops,radset,wavelength,*,packet,enabled,support):
        import lightweaver.constants as c
        from lightweaver.atomic_table import PeriodicTable
        super().__init__(eqpops,radset,wavelength)
        if len(radset.atoms)!=1 or PeriodicTable['H']not in radset.atoms:raise ValueError('Adapter restricted to H-only disposable array tests')
        self.eq=eqpops;self.radset=radset;self.packet=packet;self.enabled=enabled
        self.wave=np.array(wavelength,dtype=float,copy=True);self.original_wave=wavelength
        if self.wave.ndim!=1 or not np.isfinite(self.wave).all()or np.any(self.wave<=0)or np.any(np.diff(self.wave)<=0):raise ValueError('Invalid explicit vacuum-nm grid')
        self.mask=(self.wave>=support[0])&(self.wave<=support[1]);self.audit=EqLookupAudit(eqpops)
        self.check_state(eqpops.atmosphere)
        self.base=FastBackground(self.audit,radset,wavelength,Nthreads=1)
        from pipeline.bride_cn_diagnostic import frequency_offsets,external_state_extinction
        n=len(packet.width_Hz);width=np.asarray(packet.width_Hz)
        v=frequency_offsets(self.wave[self.mask]*10,np.full(n,packet.lambda0_air_A),width,np.full(n,c.CLight*100),np.full(n,packet.lambda0_air_A*1e-8),input_medium='vacuum',mode='exact_vacuum_frequency')
        self.added_chi=external_state_extinction(v,packet.damping_a,width,packet.lower_population_m3,np.full(n,packet.gf),packet.stimulated_emission,packet.integrated_strength_m2_Hz)
        nu=c.CLight/(self.wave[self.mask,None]*1e-9);T=np.asarray(eqpops.atmosphere.temperature)[None,:]
        planck=2*c.HPlanck*nu**3/c.CLight**2/np.expm1(c.HPlanck*nu/(c.KBoltzmann*T))
        self.added_eta=self.added_chi*planck
        if not np.isfinite(self.added_eta).all():raise ValueError('Invalid synthetic LTE emission')

    def check_state(self,atmosphere):
        actual,abundance=binding(atmosphere,self.radset.abundance)
        if actual!=self.packet.state_binding or abundance!=self.packet.abundance_binding:raise ValueError('External synthetic state/abundance binding changed')
        if not np.array_equal(self.wave,self.original_wave):raise ValueError('Bound CN wavelength grid changed')
        if len(atmosphere.temperature)!=len(self.packet.width_Hz):raise ValueError('CN depth mismatch')

    def compute_background(self,atmos,chi,eta,sca):
        self.check_state(atmos.pyAtmos)
        self.base.compute_background(atmos,chi,eta,sca)
        if self.enabled and np.any(self.added_chi):
            chi[self.mask]+=self.added_chi;eta[self.mask]+=self.added_eta
        if not all(np.isfinite(x).all()for x in (chi,eta,sca)):raise ValueError('Invalid synthetic CN background')
