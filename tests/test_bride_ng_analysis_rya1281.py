import copy
import numpy as np
import pytest


def fixtures():
    history=[dict(iteration=1,rejected_invalid_ng=True,population_change=.1,physical_rate_residual=.01,J_change=.2,cap_boundary_change=.1),
        dict(iteration=2,rejected_invalid_ng=False,population_change=.02,physical_rate_residual=.003,J_change=.01,cap_boundary_change=.01)]
    updates=[]
    for i,outer,attempt,ng,valid in [(1,1,1,True,False),(2,1,2,False,True),(3,2,1,True,True)]:
        updates.append(dict(schema='codex.bride_stat_equil_update/1',call_index=i,outer_iteration=outer,
            attempt_in_outer=attempt,active_species=['O'],updatedPops=True,ngAccelerated=[ng],
            populations_finite=[True],populations_positive=[valid]))
    return history,updates


def test_join_keeps_rejected_acceleration_separate_from_accepted_steps_and_off_refuses():
    from scripts.analyze_bride_ng_controls import join_updates
    h,u=fixtures();rows=join_updates(h,u,'record')
    assert [r['accepted_ng'] for r in rows]==[False,True]
    assert [r['native_calls'] for r in rows]==[2,1]
    with pytest.raises(ValueError,match='Disabled'):join_updates(h,u,'off')
    bad=copy.deepcopy(u);bad[1]['attempt_in_outer']=1
    with pytest.raises(ValueError,match='retry'):join_updates(h,bad,'record')
    bad=copy.deepcopy(u);bad[-1]['populations_positive']=[False]
    with pytest.raises(ValueError,match='acceptance'):join_updates(h,bad,'record')


def test_replay_prefix_requires_exact_numbers_instead_of_a_relaxed_tolerance():
    from scripts.analyze_bride_ng_controls import exact_prefix
    h,_=fixtures();assert exact_prefix(h,h[:1])==1
    bad=copy.deepcopy(h);bad[0]['population_change']=float(np.nextafter(.1,1.))
    with pytest.raises(ValueError,match='iteration1'):exact_prefix(h,bad)
    with pytest.raises(ValueError,match='Empty'):exact_prefix([],h)
