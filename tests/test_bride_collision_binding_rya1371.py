"""Whole-matrix binding refuses stale physical state and double-counted rates."""
from types import SimpleNamespace as NS
import numpy as np
import pytest
pytest.importorskip('lightweaver')
from lightweaver.atomic_table import PeriodicTable
from pipeline.bride_collision_binding import FixedNativeCollisionBinding

def setup():
    state=dict(temperature_K=np.array([5000.,6000.]),ne_m3=np.array([1e16,2e16]),height_m=np.array([100.,0.]),vturb_m_s=np.array([1000.,1000.]))
    star=np.ones((41,2));matrix=np.zeros((41,41,2));matrix[1,0]=[3.,4.];nH=np.array([1e20,2e20])
    binding=FixedNativeCollisionBinding(matrix,state,star,nH,{})
    model=NS(element=PeriodicTable['O'],levels=[NS(E=float(i),g=1.,stage=0)for i in range(41)],collisions=[binding]);binding.setup(model)
    oxygen=NS(model=model,nStar=star.copy(),nTotal=star.sum(axis=0),n=star.copy())
    eq=NS(atomicPops={model.element:oxygen});atm=NS(temperature=state['temperature_K'].copy(),ne=state['ne_m3'].copy(),height=state['height_m'].copy(),vturb=state['vturb_m_s'].copy(),nHTot=nH.copy())
    return binding,atm,eq,oxygen,model

def test_matrix_copy_and_population_updates_do_not_lookup_fake_H():
    b,a,e,o,m=setup();C=np.zeros((41,41,2));b.compute_rates(a,e,C);np.testing.assert_array_equal(C,b.matrix)
    o.n[0]*=2;C.fill(0);b.compute_rates(a,e,C);assert b.calls==2
    with pytest.raises(ValueError):b.matrix[1,0,0]=99
    with pytest.raises(ValueError,match='stacking'):b.compute_rates(a,e,C)

def test_changed_thermodynamics_nstar_levels_totals_and_extra_rates_refuse():
    for name in ['temperature','ne','nHTot','height','vturb']:
        b,a,e,o,m=setup();getattr(a,name)[0]*=1.01
        with pytest.raises(ValueError,match='state'):b.compute_rates(a,e,np.zeros((41,41,2)))
    for mutate in [lambda o,m:o.nStar.__setitem__((0,0),2),lambda o,m:o.nTotal.__setitem__(0,42),lambda o,m:setattr(m.levels[0],'g',2),lambda o,m:m.collisions.append(object())]:
        b,a,e,o,m=setup();mutate(o,m)
        with pytest.raises(ValueError):b.compute_rates(a,e,np.zeros((41,41,2)))
