import copy
from functools import partial
import numpy as np
import pytest

pytest.importorskip("lightweaver", reason="Compiled width controls require the optional Bride backend")
from tests.test_bride_feii_control_rya1281 import context
from pipeline.bride_feii_width_control import feii_width_factory


LINE=dict(line_id='controlled',lambda0_vac_nm=778.,excitation_eV=3.199,loggf=-3.,gamma_rad_s_inv=2.95e8)
RECIPE=dict(line_id='controlled',recipe='W',log_gamma6_10000=-7.91)
STARK=dict(line_id='controlled',ionization_eV=16.188,excitation_eV=3.199,wavelength_air_A=7777.099,log_gamma_stark_10000=-6.53)


@pytest.mark.parametrize('mode',['cowley','provided'])
def test_depth_width_coefficients_emission_and_cache_with_real_backend(mode):
    from lightweaver.atomic_table import PeriodicTable,KuruczPfTable
    from pipeline.bride_lte_blend import lte_line_coefficients
    bare=context('bare')
    mixed=context(LINE,partial(feii_width_factory,recipe=RECIPE,stark_row=STARK,mode=mode))
    a=mixed.atmos.pyAtmos;eq=mixed.kwargs['eqPops'];w=np.asarray(mixed.spect.wavelength);inside=(w>=777)&(w<=779)
    neutral=[]
    for name in ('H','He'):
        state=eq.atomicPops[PeriodicTable[name]]
        neutral.append(state.n[[level.stage==0 for level in state.model.levels]].sum(axis=0))
    T=a.temperature;ne=a.ne
    pressure=10**-7.91*(T/10000.)**.3*(neutral[0]+.41336*neutral[1])*1e-6
    if mode=='provided':stark=10**-6.53*(T/10000.)**(1/6)*(ne*1e-6)
    else:
        gap=16.188-(3.199+3.40*3647./7777.099)
        stark=5.42e-7*13.598*(2/gap)**2*(ne*1e-6)/9
    gamma=2.95e8+pressure+stark
    assert np.ptp(gamma)>0 and np.all(gamma>2.95e8)
    p=KuruczPfTable(atomicAbundance=eq.abundance)[PeriodicTable['Fe']]
    U=np.exp(np.interp(T,p.Tpf,p.pf[1]))
    x,e=lte_line_coefficients(w[inside],T,p.lte_ionisation(a)[1]/U,a.vturb,
        lambda0_nm=778.,mass_amu=PeriodicTable['Fe'].mass,excitation_eV=3.199,loggf=-3.,gamma_s_inv=gamma)
    for key,added in [('chi',x),('eta',e)]:
        total=np.asarray(getattr(mixed.background,key));base=np.asarray(getattr(bare.background,key))
        assert np.all(abs((total-base)[inside]-added)<=16*np.finfo(float).eps*(abs(total[inside])+abs(base[inside])+abs(added)))
        np.testing.assert_array_equal(total[~inside],base[~inside])
    np.testing.assert_array_equal(mixed.background.sca,bare.background.sca)
    clone=copy.deepcopy(mixed);np.testing.assert_array_equal(clone.background.chi,mixed.background.chi)
    eq.atomicPops[PeriodicTable['H']].n[0,0]*=1.01
    with pytest.raises(ValueError,match='neutral colliders'):mixed.background.update_background(mixed.atmos)


def test_invalid_width_identity_or_mode_refuses():
    for recipe,stark,mode in [(dict(RECIPE,line_id='wrong'),STARK,'cowley'),(RECIPE,STARK,'neutral'),(dict(RECIPE,recipe='S'),STARK,'provided')]:
        with pytest.raises(ValueError):feii_width_factory(None,LINE,(777.,779.),recipe,stark,mode)
