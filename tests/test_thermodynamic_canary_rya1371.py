import numpy as np
import pytest
from pipeline.thermodynamic_canary import occupation,bb_thermal_defect,bf_thermal_defect


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
