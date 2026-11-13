import copy
import pytest
from pipeline.native_reference_controls import validate_observer_plan


def base():return dict(scope='qualified',controls=dict(nq_fac=2.,rtol=1e-6,threads=1),atom=dict(sha256='held'),output_grid=dict(frame='air_A',windows=[[6154.,6162.]]))


def test_only_frequency_refinement_changes_and_metadata_are_permitted():
    original=base();refined=copy.deepcopy(original);refined['controls']['nq_fac']=4.;refined['scope']='explicit numerical test'
    assert validate_observer_plan(refined,original)==4.
    assert original==base()


@pytest.mark.parametrize('change',['rtol','atom','frame','nq'])
def test_other_physical_or_control_changes_refuse(change):
    original=base();p=copy.deepcopy(original);p['controls']['nq_fac']=4.
    if change=='rtol':p['controls']['rtol']=1e-5
    elif change=='atom':p['atom']['sha256']='other'
    elif change=='frame':p['output_grid']['frame']='vacuum_nm'
    else:p['controls']['nq_fac']=8.
    with pytest.raises(ValueError):validate_observer_plan(p,original)


def test_shared_radiation_is_an_explicit_separate_one_factor_control():
    from pipeline.native_reference_controls import validate_shared_radiation_plan
    qualified=base();qualified['controls']['nq_fac']=4.
    p=copy.deepcopy(qualified);p['scope']='shared field';p['controls']['decouple_continuum']=False
    validate_shared_radiation_plan(p,qualified)
    p['controls']['rtol']=1e-5
    with pytest.raises(ValueError):validate_shared_radiation_plan(p,qualified)


def test_shared_frequency8_preserves_every_other_physical_axis():
    from pipeline.native_reference_controls import validate_shared_refinement_plan
    qualified=base();qualified['controls'].update(nq_fac=4.,decouple_continuum=False)
    p=copy.deepcopy(qualified);p['controls']['nq_fac']=8.;p['scope']='predeclared refinement'
    validate_shared_refinement_plan(p,qualified)
    p['controls']['decouple_continuum']=True
    with pytest.raises(ValueError):validate_shared_refinement_plan(p,qualified)
    p['controls']['decouple_continuum']=False;p['controls']['threads']=2
    with pytest.raises(ValueError):validate_shared_refinement_plan(p,qualified)


def test_shared16_requires_the_actual_preceding_shared8_case():
    from pipeline.native_reference_controls import validate_shared_refinement_plan
    qualified=base();qualified['controls'].update(nq_fac=8.,decouple_continuum=False)
    plan=copy.deepcopy(qualified);plan['controls']['nq_fac']=16.
    validate_shared_refinement_plan(plan,qualified,target_factor=16.)
    bad=copy.deepcopy(qualified);bad['controls']['nq_fac']=4.
    with pytest.raises(ValueError):validate_shared_refinement_plan(plan,bad,target_factor=16.)
    plan['atom']['sha256']='other'
    with pytest.raises(ValueError):validate_shared_refinement_plan(plan,qualified,target_factor=16.)
