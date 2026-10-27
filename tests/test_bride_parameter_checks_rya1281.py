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


def test_full_structure_interpolation_recovers_affine_fields_endpoints_and_rejects_extrapolation():
    from pipeline.bride_parameter_checks import interpolate_marcs_structure
    nodes=[]
    for teff in (5750,6000):
        for logg in (4.,4.5):
            s=dict(ndep=4,abundances_dex=[12.,10.93],hydrogen_mass_fraction=.73826,model_id='fixture',
                vturb_m_s=1000,logg_cgs=logg,height_m=np.array([3.,2.,1.,0.]),
                temperature_K=np.arange(4.)+teff+10*logg,
                electron_pressure_Pa=np.full(4,teff/100+logg),rho_kg_m3=np.full(4,teff/10000+logg/10),
                logtau_ross=np.arange(4.),logtau_500=np.arange(4.)+.01*logg)
            nodes.append(((teff,logg),s))
    result,metadata=interpolate_marcs_structure(nodes,(5772,4.438))
    np.testing.assert_allclose(result['temperature_K'],np.arange(4.)+5772+44.38)
    np.testing.assert_allclose(result['rho_kg_m3'],5772/10000+4.438/10)
    assert sum(metadata['weights'])==pytest.approx(1)
    endpoint,_=interpolate_marcs_structure(nodes,(5750,4.))
    np.testing.assert_array_equal(endpoint['temperature_K'],nodes[0][1]['temperature_K'])
    with pytest.raises(ValueError,match='outside'):interpolate_marcs_structure(nodes,(5700,4.438))
    nodes[0][1]['logtau_ross'][0]=-.2
    with pytest.raises(ValueError,match='layer grid'):interpolate_marcs_structure(nodes,(5772,4.438))
