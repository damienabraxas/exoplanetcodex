import numpy as np
import pytest


@pytest.fixture
def oracle():
    from config.constants import codex_path
    from pipeline.bride_native_spline import NativeTautsp
    directory=codex_path('experiments.bride_multi_spline_oracle')
    if not (directory/'manifest.json').exists():pytest.skip('Held compiled native oracle unavailable on this host')
    return NativeTautsp(directory)


def test_actual_native_spline_affine_endpoint_log_laws_and_domain(oracle):
    x=np.arange(1000.,6000.,1000.)
    t=np.array([1000.,1500.,3500.,5000.])
    np.testing.assert_allclose(oracle.interpolate(x,2*x+1,t,logarithmic=False),2*t+1,rtol=1e-14)
    np.testing.assert_allclose(oracle.interpolate(x,np.exp(x/1000),t,logarithmic=True),np.exp(t/1000),rtol=1e-14)
    for query in (np.array([999.]),np.array([5001.]),np.array([np.nan])):
        with pytest.raises(ValueError):oracle.interpolate(x,x,query,logarithmic=False)


def test_native_collision_units_direction_and_lte_balance(oracle):
    pytest.importorskip('lightweaver')
    from types import SimpleNamespace as S
    from pipeline.multi_atom import MultiCollision
    from pipeline.bride_native_spline import install_native_collision_control
    import lightweaver.constants as c
    T=np.array([4000.,5000.]);ne=np.array([1e16,2e16]);H=np.array([1e20,2e20])
    grid=np.array([3000.,4000.,5000.,6000.])
    raw=[MultiCollision(kind,2,1,np.full(4,q),grid) for kind,q in [('CE',1e-9),('CI',1e-10),('CH',1e-12)]]
    raw.append(MultiCollision('CH0',1,2,np.full(4,1e-12),grid))
    model=S(levels=[S(E=0.,g=2.),S(E=1e4,g=4.)],element='O')
    atmosphere=S(temperature=T,ne=ne);ctx=S(activeAtoms=[S(atomicModel=model)],atmos=S(pyAtmos=atmosphere))
    install_native_collision_control(ctx,S(collisions=raw),oracle,source_directions=True)
    lte=np.array([[10.,20.],[2.,4.]])
    class Populations:
        atomicPops={'O':S(nStar=lte)}
        def __getitem__(self,key):return H[None,:]
    matrices=[]
    for rate in model.collisions:
        m=np.zeros((2,2,2));rate.compute_rates(atmosphere,Populations(),m);matrices.append(m)
    np.testing.assert_allclose(matrices[0][0,1],1e-15*ne*np.sqrt(T)*.5)
    np.testing.assert_allclose(matrices[1][1,0],1e-16*ne*np.sqrt(T)*np.exp(-1e4*c.HC_CM/(c.KBoltzmann*T)))
    np.testing.assert_allclose(matrices[2][0,1],1e-18*H)
    np.testing.assert_allclose(matrices[3][1,0],1e-18*H)
    np.testing.assert_array_equal(matrices[3][0,1],0.)
    for m in matrices[:3]:np.testing.assert_allclose(m[0,1]*lte[1],m[1,0]*lte[0])
