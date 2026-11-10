import numpy as np
import pytest


def test_invalid_acceleration_retries_same_equation_once_without_clipping():
    from types import SimpleNamespace as S
    from pipeline.bride_parameter_checks import positive_stat_equil
    active=S(n=np.array([[1.],[2.]]));before=active.n.copy()
    class Context:
        count=0
        cleared=False
        def stat_equil(self):
            self.count+=1
            active.n[:]=np.array([[-1.],[4.]]) if not self.cleared else np.array([[.5],[2.5]])
        def clear_ng(self):
            np.testing.assert_array_equal(active.n,before)
            self.cleared=True
    ctx=Context()
    assert positive_stat_equil(ctx,active,before,reject_invalid_ng=True)
    assert ctx.count==2
    np.testing.assert_array_equal(active.n,[[.5],[2.5]])
    ctx=Context();ctx.clear_ng=lambda:None
    with pytest.raises(ValueError):positive_stat_equil(ctx,active,before,reject_invalid_ng=True)


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


def test_legacy_ch_coefficient_is_downward_and_balances_lte():
    pytest.importorskip('lightweaver')
    from types import SimpleNamespace as S
    from pipeline.multi_atom import MultiCollision
    from pipeline.bride_gencol_controls import install_gencol_control
    grid=np.array([4000.,6000.]);rates=np.array([1e-16,1e-16])
    source=S(collisions=[MultiCollision('CH',2,1,rates*1e6,grid)])
    model=S(collisions=[S(i=0,j=1,temperature=grid,rates=rates)],levels=[S(),S()],element='O')
    ctx=S(activeAtoms=[S(atomicModel=model)])
    install_gencol_control(ctx,source,change_ch=True)
    matrix=np.zeros((2,2,1));lte=np.array([[10.],[2.]])
    eq=S(atomicPops={'O':S(nStar=lte)})
    class Populations:
        atomicPops=eq.atomicPops
        def __getitem__(self,key):return np.array([[1e20]])
    model.collisions[0].compute_rates(S(temperature=np.array([5000.])),Populations(),matrix)
    np.testing.assert_allclose(matrix[0,1],1e4)
    np.testing.assert_allclose(matrix[1,0],2e3)
    np.testing.assert_allclose(matrix[0,1]*lte[1],matrix[1,0]*lte[0])
