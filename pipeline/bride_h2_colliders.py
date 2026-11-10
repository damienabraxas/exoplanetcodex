"""Experimental FeI collider-only H2 intervention; native EOS/RT remain atomic.

This replaces pressure perturber densities only. It never mutates the native
H populations, electron density, continuum or O collision rates.
"""
from functools import partial
import numpy as np
from pipeline.bride_fe_stark import StarkFeIBackground,stark_fe_factory,stark_gamma
from pipeline.bride_fe_pressure import pressure_gamma


def h2_collider_factory(base_factory,lines,support_vac_nm,recipes,stark_rows,stark_mode,tables,mode):
    if mode not in ('off','molecular'):raise ValueError('Unknown H2 collider intervention')
    factory=stark_fe_factory(base_factory,lines,support_vac_nm,recipes,stark_rows,stark_mode)
    return partial(H2ColliderFeIBackground,**factory.keywords,h2_tables=tables,h2_mode=mode)


class H2ColliderFeIBackground(StarkFeIBackground):
    def __init__(self,*args,h2_tables,h2_mode,**kwargs):
        self.h2_tables=h2_tables;self.h2_mode=h2_mode;self.hydrogen_inputs=None
        self.hydrogen_solution=None;self.collider_atmosphere=None
        super().__init__(*args,**kwargs)

    def compute_background(self,atmos,chi,eta,sca):
        from lightweaver.atomic_table import PeriodicTable
        a=atmos.pyAtmos;state=self.eq.atomicPops[PeriodicTable['H']]
        # Exact dependency checks supplement the inherited complete atmosphere guard.
        inputs=[a.temperature,a.nHTot,state.n,self.eq.HminPops]
        if self.hydrogen_inputs is None:self.hydrogen_inputs=[np.array(v,copy=True) for v in inputs]
        if any(not np.array_equal(v,old) for v,old in zip(inputs,self.hydrogen_inputs)):
            raise ValueError('Cached H2 collider intervention refuses changed hydrogen inputs')
        self.collider_atmosphere=a
        super().compute_background(atmos,chi,eta,sca)

    def current_colliders(self):
        from lightweaver.atomic_table import PeriodicTable
        H,He=super().current_colliders()
        if self.h2_mode=='off':return [H,He,np.zeros_like(H)]
        state=self.eq.atomicPops[PeriodicTable['H']];a=self.collider_atmosphere
        total_atomic=state.n.sum(axis=0)
        self.hydrogen_solution=self.h2_tables.hydrogen_network(a.temperature,a.nHTot,H/total_atomic,self.eq.HminPops/H)
        return [self.hydrogen_solution['HI_m3'],He,self.hydrogen_solution['H2_m3']]

    def line_width(self,line,atmosphere):
        return (line['gamma_rad_s_inv']+pressure_gamma(atmosphere.temperature,*self.colliders,self.recipes[line['line_id']])
                +stark_gamma(atmosphere.temperature,atmosphere.ne,self.stark_rows[line['line_id']],self.stark_mode))
