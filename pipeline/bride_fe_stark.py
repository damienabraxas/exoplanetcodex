"""Explicit alternative neutral-Fe Stark prescriptions from pinned TS DEPTH."""
import numpy as np
from functools import partial
from pipeline.bride_fe_pressure import PressureFeIBackground


def stark_gamma(temperature,ne_m3,row,mode):
    T,ne=np.broadcast_arrays(temperature,ne_m3)
    if not np.isfinite(T).all()or not np.isfinite(ne).all()or np.any(T<=0)or np.any(ne<0):raise ValueError('Invalid Stark inputs')
    if mode=='none':return np.zeros_like(T)
    if mode=='provided':
        coefficient=row['log_gamma_stark_10000']
        if not np.isfinite(coefficient)or coefficient>=0:raise ValueError('Supplied Stark log coefficient must be negative/nonzero')
        gamma=10**coefficient*(T/10000.)**(1/6)*(ne*1e-6)
    elif mode=='cowley':
        gap=row['ionization_eV']-(row['excitation_eV']+3.40*3647./row['wavelength_air_A'])
        if not np.isfinite(gap)or gap<=0:raise ValueError('Unqualified Cowley upper-level gap')
        gamma=2.26e-7*13.598/gap**2*(ne*1e-6)
    else:raise ValueError('Unknown explicit Stark mode')
    if not np.isfinite(gamma).all()or np.any(gamma<0):raise ValueError('Invalid Stark width')
    return gamma


def stark_fe_factory(base_factory,lines,support_vac_nm,recipes,stark_rows,mode):
    from pipeline.bride_fe_pressure import pressure_fe_factory
    pressure_fe_factory(base_factory,lines,support_vac_nm,recipes)
    if mode not in ('none','provided','cowley'):raise ValueError('Unknown Stark prescription')
    rows={r['line_id']:dict(r)for r in stark_rows}
    if len(rows)!=len(stark_rows)or set(rows)!={r['line_id']for r in lines}:raise ValueError('Stark line identities differ')
    return partial(StarkFeIBackground,base_factory=base_factory,lines=tuple(dict(r)for r in lines),support=tuple(support_vac_nm),
        recipes={r['line_id']:dict(r)for r in recipes},stark_rows=rows,stark_mode=mode)


class StarkFeIBackground(PressureFeIBackground):
    def __init__(self,*args,stark_rows,stark_mode,**kwargs):
        self.stark_rows=stark_rows;self.stark_mode=stark_mode
        super().__init__(*args,**kwargs)

    def line_width(self,line,atmosphere):
        return super().line_width(line,atmosphere)+stark_gamma(atmosphere.temperature,atmosphere.ne,
            self.stark_rows[line['line_id']],self.stark_mode)
