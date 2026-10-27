import numpy as np
import pytest


def test_directed_charge_exchange_preserves_source_order_and_log_interpolation():
    pytest.importorskip('lightweaver')
    from types import SimpleNamespace as S
    from pipeline.multi_atom import MultiCollision
    from pipeline.bride_gencol_controls import install_gencol_control
    source=S(collisions=[MultiCollision('CH0',1,2,np.array([1e-12,1e-8]),np.array([4000.,6000.]))])
    model=S(collisions=[None]);ctx=S(activeAtoms=[S(atomicModel=model)])
    install_gencol_control(ctx,source,change_ch=False)
    matrix=np.zeros((2,2,1));rate=model.collisions[0]
    rate.compute_rates(S(temperature=np.array([5000.])),{'H':np.array([[1e20]])},matrix)
    assert matrix[0,1,0]==0 # opposite direction must stay absent
    np.testing.assert_allclose(matrix[1,0,0],1e4,rtol=1e-14)
    with pytest.raises(ValueError):
        rate.compute_rates(S(temperature=np.array([7000.])),{'H':np.array([[1e20]])},matrix)
