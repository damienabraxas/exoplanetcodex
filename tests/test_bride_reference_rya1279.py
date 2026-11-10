"""Exact held-record identity/layout and MARCS unit boundaries, synthetic bytes."""
import numpy as np
import pytest
from pipeline.bride_reference import read_exact_departure_record,read_marcs_structure


def test_non_square_level_depth_payload_and_exact_abundance(tmp_path):
    tau=np.array([.01,.1,1.],dtype='<f8')
    dep=np.array([[100*l+d for d in range(3)]for l in range(4)],dtype='<f8')
    raw=b'fixture'.ljust(500,b' ')+np.array([3,4],dtype='<i4').tobytes()+tau.tobytes()+dep.tobytes()
    p=tmp_path/'grid.bin';p.write_bytes(b' '*1000+raw)
    aux="'fixture' 5750 4.5 0 0 0 1 8.7 1001\n"
    d=read_exact_departure_record(p,aux,model_id='fixture',abundance=8.7,nlevels=4,ndep=3)
    np.testing.assert_array_equal(d['departures'],dep)
    np.testing.assert_array_equal(d['tau'],tau)
    with pytest.raises(ValueError,match='found 0'):
        read_exact_departure_record(p,aux,model_id='fixture',abundance=8.69,nlevels=4,ndep=3)
    with pytest.raises(ValueError,match='found 2'):
        read_exact_departure_record(p,aux+aux,model_id='fixture',abundance=8.7,nlevels=4,ndep=3)
    p.write_bytes(b' '*1000+raw.replace(b'fixture',b'WRONGID',1))
    with pytest.raises(ValueError,match='identity'):
        read_exact_departure_record(p,aux,model_id='fixture',abundance=8.7,nlevels=4,ndep=3)


MARCS='''fixture
10000 Surface gravity [cm/s2]
1.5 Microturbulence parameter [km/s]
0.7 0.28 0.02 are X, Y and Z
3 Number of depth points
Model structure
k lgTauR lgTau5 Depth T Pe Pg Prad Pturb
1 -2 -2.1 -100 4000 10 100 0 0
2 -1 -1.1 0 5000 20 200 0 0
3 0 -0.1 100 6000 30 300 0 0
k lgTauR KappaRoss Density Mu Vconv Fconv/F RHOX
1 -2 .1 1e-9 1 0 0 .1
2 -1 .2 2e-9 1 0 0 .2
3 0 .3 3e-9 1 0 0 .3
'''

MARCS=MARCS.replace('3 Number of depth points',
    'Logarithmic chemical number abundances, H always 12.00\n'+
    ' '.join(['12.00']+['-99.00']*91)+'\n3 Number of depth points')


def test_marcs_rows_units_and_source_identity(tmp_path):
    p=tmp_path/'fixture.mod';p.write_text(MARCS)
    d=read_marcs_structure(p)
    assert d['model_id']=='fixture'
    np.testing.assert_array_equal(d['height_m'],[1.,0.,-1.])
    np.testing.assert_array_equal(d['temperature_K'],[4000.,5000.,6000.])
    np.testing.assert_allclose(d['rho_kg_m3'],[1e-6,2e-6,3e-6],rtol=1e-15,atol=0)
    np.testing.assert_array_equal(d['electron_pressure_Pa'],[1.,2.,3.])
    np.testing.assert_allclose(d['column_mass_kg_m2'],[1.,2.,3.],rtol=1e-15,atol=0)
    assert d['hydrogen_mass_fraction']==.7
    assert d['logg_cgs']==4. and d['vturb_m_s']==1500.


@pytest.mark.parametrize('old,new',[('2 -1 .2','2 -1.5 .2'),
    ('2 -1 .2','3 -1 .2'),('100 6000','-100 6000'),('5000 20','5000 -20')])
def test_inconsistent_marcs_geometry_identity_or_physics_is_refused(tmp_path,old,new):
    p=tmp_path/'fixture.mod';p.write_text(MARCS.replace(old,new))
    with pytest.raises(ValueError):read_marcs_structure(p)
