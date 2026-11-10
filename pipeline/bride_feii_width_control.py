"""Opt-in fixed-column FeII widths; native W pressure and explicit ionic Stark.

H and He stage-zero populations come from the live backend EOS. H2 remains zero
as in the frozen control. This adds opacity/emission, not population feedback.
"""
from functools import partial
import numpy as np
from pipeline.bride_feii_control import FeIIControl,feii_factory
from pipeline.bride_fe_pressure import pressure_gamma
from pipeline.bride_feii_profile import feii_stark_gamma


def feii_width_factory(base_factory,line,support,recipe,stark_row,mode):
    feii_factory(base_factory,line,support)
    if line is None or mode not in ('cowley','provided'):raise ValueError('Explicit populated FeII width mode required')
    if recipe['recipe']!='W':raise ValueError('Only verified FeII W branch admitted to this control')
    if recipe['line_id']!=line['line_id'] or stark_row['line_id']!=line['line_id']:raise ValueError('FeII width identity mismatch')
    pressure_gamma(5772.,0.,0.,0.,recipe)
    feii_stark_gamma(5772.,0.,stark_row,mode)
    return partial(FeIIWidthControl,base_factory=base_factory,line=dict(line),support=tuple(support),
        recipe=dict(recipe),stark_row=dict(stark_row),mode=mode)


class FeIIWidthControl(FeIIControl):
    def __init__(self,eqpops,*args,recipe,stark_row,mode,**kwargs):
        self.eq=eqpops;self.recipe=recipe;self.stark_row=stark_row;self.mode=mode;self.colliders=None
        super().__init__(eqpops,*args,**kwargs)

    def current_colliders(self):
        from lightweaver.atomic_table import PeriodicTable
        values=[]
        for name in ('H','He'):
            state=self.eq.atomicPops[PeriodicTable[name]]
            values.append(state.n[[level.stage==0 for level in state.model.levels]].sum(axis=0))
        return values

    def compute_background(self,atmos,chi,eta,sca):
        current=self.current_colliders()
        if self.colliders is None:self.colliders=[value.copy() for value in current]
        if any(not np.array_equal(a,b) for a,b in zip(current,self.colliders)):
            raise ValueError('FeII cached neutral colliders changed')
        super().compute_background(atmos,chi,eta,sca)

    def line_width(self,atmosphere):
        return (self.line['gamma_rad_s_inv']+pressure_gamma(atmosphere.temperature,*self.colliders,
            np.zeros_like(atmosphere.temperature),self.recipe)
            +feii_stark_gamma(atmosphere.temperature,atmosphere.ne,self.stark_row,self.mode))
