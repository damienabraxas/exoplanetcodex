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
    model=S(levels=[S(E=0.,g=2.,stage=0),S(E=1e4,g=4.,stage=1)],element='O')
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
    install_native_collision_control(ctx,S(collisions=raw),oracle,source_directions=True,reversible_ch0=True)
    paired=np.zeros((2,2,2));model.collisions[-1].compute_rates(atmosphere,Populations(),paired)
    np.testing.assert_allclose(paired[1,0],matrices[-1][1,0])
    np.testing.assert_allclose(paired[0,1]*lte[1],paired[1,0]*lte[0])


def test_h_background_control_keeps_populations_lte_and_oxygen_identity():
    pytest.importorskip('lightweaver')
    from pipeline.multi_atom import MultiAtom,MultiLevel,MultiLine
    from pipeline.bride_reference import make_reference_context
    from lightweaver.atomic_table import PeriodicTable
    atom=MultiAtom(symbol='O',abundance=8.7,weight=16.,
        levels=[MultiLevel(0.,5.,'synthetic lower',1,1),MultiLevel(12000.,7.,'synthetic upper',1,2),MultiLevel(110000.,4.,'synthetic ion',2,3)],
        lines=[MultiLine(2,1,.4,9,5.,3.89e7,0.,0.,'VOIGT')],continua=[],temperature=np.array([3000.,10000.]),collisions=[])
    abundances=[12.]+[-99.]*91;abundances[1]=10.93
    structure=dict(ndep=5,height_m=np.linspace(1000.,0.,5),temperature_K=np.full(5,6000.),
        electron_pressure_Pa=np.full(5,.01),rho_kg_m3=np.full(5,1e-5),hydrogen_mass_fraction=.7,
        abundances_dex=abundances,logg_cgs=4.5,vturb_m_s=1000.)
    plain=make_reference_context(structure,atom,rays=2,qscale=1)
    with_h=make_reference_context(structure,atom,rays=2,qscale=1,lte_h_background_lines=True)
    assert len(plain.detailedAtoms)==0 and len(with_h.detailedAtoms)==1
    np.testing.assert_array_equal(plain.activeAtoms[0].nStar,with_h.activeAtoms[0].nStar)
    hydrogen=with_h.kwargs['eqPops'].atomicPops[PeriodicTable['H']]
    np.testing.assert_array_equal(hydrogen.n,hydrogen.nStar)
    assert len(with_h.spect.wavelength)>len(plain.spect.wavelength)
