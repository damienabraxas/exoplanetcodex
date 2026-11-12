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
