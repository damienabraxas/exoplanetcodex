from types import SimpleNamespace
import numpy as np
import pytest
from pipeline.bride_parameter_checks import controlled_structure,assert_collision_coverage


def test_microturbulence_and_depth_selection_do_not_mutate_original_structure():
    original={'ndep':6,'vturb_m_s':1000.}
    names=('height_m','temperature_K','electron_pressure_Pa','rho_kg_m3','logtau_ross','logtau_500')
    original.update({k:np.arange(6.)for k in names})
    changed=controlled_structure(original,xi=.9,drop_bottom=1)
    assert changed['vturb_m_s']==900 and changed['ndep']==5
    changed['temperature_K'][0]=999
    assert original['temperature_K'][0]==0 and original['ndep']==6
    with pytest.raises(ValueError):controlled_structure(original,xi=1,drop_bottom=4)


def test_collision_domain_checks_each_record_without_clamping():
    rates=[SimpleNamespace(temperature=np.array([1000.,10000.]),kind='CH',i=1,j=2),
           SimpleNamespace(temperature=np.array([2000.,200000.]),kind='CE',i=1,j=2)]
    atom=SimpleNamespace(collisions=rates)
    assert_collision_coverage(atom,[4000,9900])
    with pytest.raises(ValueError,match='exceeds 1 collision'):assert_collision_coverage(atom,[4000,10161.4])
