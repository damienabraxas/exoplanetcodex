"""Explicit SI collider interface to held TS neutral-collision width equations."""
import math
import numpy as np


def pressure_gamma(temperature,nH0_m3,nHe0_m3,nH2_m3,recipe):
    from lightweaver.atomic_table import PeriodicTable
    import lightweaver.constants as c
    T,H,He,H2=np.broadcast_arrays(temperature,nH0_m3,nHe0_m3,nH2_m3)
    if (not all(np.isfinite(x).all()for x in (T,H,He,H2)) or np.any(T<=0)
            or any(np.any(x<0)for x in (H,He,H2))):raise ValueError('Invalid neutral collider inputs')
    density=(H+.41336*He+.85*H2)*1e-6 # m^-3 -> cm^-3
    if recipe['recipe']=='W':
        coefficient=10**recipe['log_gamma6_10000']*(T/10000.)**.3
    elif recipe['recipe']in ('A','B','C','S'):
        alpha=recipe['alpha'];sigma=recipe['sigma_au']
        if not np.isfinite([alpha,sigma]).all()or not 0<=alpha<1 or sigma<=0:raise ValueError('Invalid ABO data')
        mH=PeriodicTable['H'].mass*c.Amu*1000;mFe=PeriodicTable['Fe'].mass*c.Amu*1000
        reduced=1/(1/mH+1/mFe)
        velocity=np.sqrt(8*c.KBoltzmann*1e7*T/(np.pi*reduced))
        coefficient=2*(4/np.pi)**(alpha/2)*math.gamma(2-alpha/2)*velocity*(sigma*2.8002e-17)*(velocity/1e6)**(-alpha)
    else:raise ValueError('Unimplemented actual pressure recipe')
    result=coefficient*density
    if not np.isfinite(result).all()or np.any(result<0):raise ValueError('Invalid pressure width')
    return result


from functools import partial
from pipeline.bride_trace_fe_lines import TraceFeIBackground,validate_lines,trace_fe_factory


def pressure_fe_factory(base_factory,lines,support_vac_nm,recipes):
    trace_fe_factory(base_factory,lines,support_vac_nm) # same support/data refusal contract
    mapping={r['line_id']:dict(r) for r in recipes}
    if len(mapping)!=len(recipes)or set(mapping)!={r['line_id']for r in lines}:raise ValueError('Pressure recipe identities differ')
    if any(r['recipe']not in ('W','A','B','C','S')for r in recipes):raise ValueError('Unknown native recipe')
    return partial(PressureFeIBackground,base_factory=base_factory,lines=tuple(dict(r)for r in lines),
        support=tuple(support_vac_nm),recipes=mapping)


class PressureFeIBackground(TraceFeIBackground):
    """Explicit H/He-only diagnostic; H2 and Stark are deliberately absent."""
    def __init__(self,*args,recipes,**kwargs):
        self.recipes=recipes;self.colliders=None
        super().__init__(*args,**kwargs)

    def current_colliders(self):
        from lightweaver.atomic_table import PeriodicTable
        values=[]
        for name in ('H','He'):
            state=self.eq.atomicPops[PeriodicTable[name]]
            values.append(state.n[[l.stage==0 for l in state.model.levels]].sum(axis=0))
        return values

    def compute_background(self,atmos,chi,eta,sca):
        current=self.current_colliders()
        if self.colliders is None:self.colliders=[v.copy()for v in current]
        if any(not np.array_equal(a,b)for a,b in zip(current,self.colliders)):
            raise ValueError('Cached pressure widths refuse changed neutral colliders')
        super().compute_background(atmos,chi,eta,sca)

    def line_width(self,line,atmosphere):
        return line['gamma_rad_s_inv']+pressure_gamma(atmosphere.temperature,*self.colliders,
            np.zeros_like(atmosphere.temperature),self.recipes[line['line_id']])
