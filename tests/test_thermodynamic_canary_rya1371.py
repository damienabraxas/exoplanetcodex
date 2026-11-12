import numpy as np
import pytest
from pipeline.thermodynamic_canary import occupation,bb_thermal_defect,bf_thermal_defect
from pipeline.thermodynamic_canary import collision_balance


def test_consistent_bb_and_bf_blackbody_identities_close_at_machine_precision():
    x=np.array([[.1,1.,10.,40.]]);beta=np.exp(-x);n=occupation(x);bare=np.full_like(x,3.)
    assert np.max(abs(bb_thermal_defect(beta,n)))<5e-15
    assert np.max(abs(bf_thermal_defect(bare*n,bare,beta)))<5e-15


def test_inconsistent_star_or_background_is_discriminated():
    x=np.array([[1.,10.]]);beta=np.exp(-x);n=occupation(x)
    assert np.max(abs(bb_thermal_defect(beta*1.001,n)))>4e-4
    assert np.max(abs(bf_thermal_defect(n*.999,np.ones_like(x),beta)))>3e-4
    with pytest.raises(ValueError):occupation([0.])
    with pytest.raises(ValueError):bb_thermal_defect([1.1],[.1])


def test_collision_balance_distinguishes_reciprocity_missing_and_one_way():
    n=np.array([[2.],[1.],[3.]])
    C=np.zeros((3,3,1));C[1,0]=3;C[0,1]=6
    defect,present=collision_balance(C,n)
    assert defect[1,0,0]==0 and not present[2,0,0]
    assert np.isnan(defect[2,0,0])
    changed=n.copy();changed[1]*=2
    assert collision_balance(C,changed)[0][1,0,0]==pytest.approx(-1/3)
    C[2,0]=1
    assert collision_balance(C,n)[0][2,0,0]==1
    C[0,0]=1
    with pytest.raises(ValueError):collision_balance(C,n)
