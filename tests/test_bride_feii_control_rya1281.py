import copy
import numpy as np
import pytest


def context(line,wrapper=None):
    from tests.test_bride_solar_control_rya1281 import inputs
    from pipeline.bride_solar_control import bound_fe_factory
    from pipeline.bride_reference import make_reference_context
    from pipeline.bride_feii_control import feii_factory
    from pipeline.bride_top_cap import FiniteColumnCap
    folder,base,s,atom,contract=inputs('fine_T5772')
    provider,_=bound_fe_factory(folder,base,s,atom.abundance,column_contract=contract)
    if line!='bare':provider=(wrapper or feii_factory)(provider,line,(777.,779.))
    return make_reference_context(s,atom,rays=2,qscale=1,fe_background=provider,
        extra_wavelengths_vac_nm=np.array([777.5,778.,778.5]),experimental_column_cap=FiniteColumnCap(s['column_mass_kg_m2'][0]/s['rho_kg_m3'][0]))


def test_empty_identity_feii_stage_and_thermal_emission_with_real_backend():
    import lightweaver.constants as c
    from lightweaver.atomic_table import PeriodicTable,KuruczPfTable
    from pipeline.bride_lte_blend import lte_line_coefficients
    line=dict(line_id='controlled',lambda0_vac_nm=778.,excitation_eV=4.,loggf=-2.,gamma_rad_s_inv=1e8)
    bare,empty,mixed=context('bare'),context(None),context(line)
    for k in ('chi','eta','sca'):np.testing.assert_array_equal(getattr(bare.background,k),getattr(empty.background,k))
    np.testing.assert_array_equal(mixed.background.sca,bare.background.sca)
    w=np.asarray(mixed.spect.wavelength);a=mixed.atmos.pyAtmos;inside=(w>=777)&(w<=779)
    p=KuruczPfTable(atomicAbundance=mixed.kwargs['eqPops'].abundance)[PeriodicTable['Fe']]
    U=np.exp(np.interp(a.temperature,p.Tpf,p.pf[1]))
    expected,emission=lte_line_coefficients(w[inside],a.temperature,p.lte_ionisation(a)[1]/U,a.vturb,
        lambda0_nm=778.,mass_amu=PeriodicTable['Fe'].mass,excitation_eV=4.,loggf=-2.,gamma_s_inv=1e8)
    for key,added in [('chi',expected),('eta',emission)]:
        total=np.asarray(getattr(mixed.background,key));base=np.asarray(getattr(bare.background,key))
        assert np.all(abs((total-base)[inside]-added)<=16*np.finfo(float).eps*(abs(total[inside])+abs(base[inside])+abs(added)))
        np.testing.assert_array_equal((total-base)[~inside],0)
    clone=copy.deepcopy(mixed);np.testing.assert_array_equal(clone.background.chi,mixed.background.chi)
    a.vturb[0]=np.nextafter(a.vturb[0],np.inf)
    with pytest.raises(ValueError):mixed.background.update_background(mixed.atmos)


def test_invalid_or_outside_line_refuses():
    from pipeline.bride_feii_control import feii_factory
    line=dict(line_id='x',lambda0_vac_nm=778.,excitation_eV=4.,loggf=-2.,gamma_rad_s_inv=1e8)
    with pytest.raises(ValueError):feii_factory(None,dict(line,gamma_rad_s_inv=0),(777.,779.))
    with pytest.raises(ValueError):feii_factory(None,line,(779.,780.))
