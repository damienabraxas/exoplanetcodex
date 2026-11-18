import numpy as np
import pytest
from pipeline.native_gross_rates import physical_matrix,read_state


def test_thermal_line_and_bound_free_have_distinct_native_normalization():
    state=dict(nk=3,nz=1,nr=2,nline=1,lower=np.array([0,1]),upper=np.array([1,2]),nstar_cm3=np.array([[8.],[2.],[1.]]),norm=np.array([[.5],[1.]]))
    # BB g=2, beta=.5, th=4*.5; thermal LTE flow 8*2 == 2*8.
    R=physical_matrix(state,(np.array([[4.],[3.]]),np.array([[0.],[5.]]),{}),[2.,1.,1.])
    assert R[1,0,0]==2.;assert R[0,1,0]==8.
    assert R[2,1,0]==3.;assert R[1,2,0]==10.


def test_truncated_observer_is_refused(tmp_path):
    p=tmp_path/'bad.bin';p.write_bytes(np.array([1371,1],dtype='<i4').tobytes())
    with pytest.raises(ValueError,match='Truncated'):read_state(p)


def test_sweep_stage_and_exact_population_identity_are_required(tmp_path):
    from pipeline.native_gross_rates import read_integrals
    p=tmp_path/'integral.bin'
    state=dict(iteration=13,zs=1,nz=1,nk=2,nr=1,populations_cm3=np.array([[8.],[2.]]),nstar_cm3=np.array([[8.],[2.]]))
    def write(pop):
        with p.open('wb')as f:
            np.array([1371,1,1,12,1,1,1,1,1,2,1,1],dtype='<i4').tofile(f)
            np.array([1],dtype='<i4').tofile(f)
            np.array([4.,0.],dtype='<f4').tofile(f)
            np.asarray(pop,dtype='<f8').tofile(f)
            np.array([8.,2.],dtype='<f8').tofile(f)
    write([8.,2.]);a,b,audit=read_integrals(p,state,radiation_iteration=12)
    assert a[0,0]==4. and b[0,0]==0. and audit['records']==1
    with pytest.raises(ValueError,match='stage mismatch'):read_integrals(p,state,radiation_iteration=13)
    write([8.,2.01])
    with pytest.raises(ValueError,match='different population stage'):read_integrals(p,state,radiation_iteration=12)
