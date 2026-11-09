from tests.bride_test_requirements import requires_lightweaver
import copy
import numpy as np
import pytest


def context(lines):
    from tests.test_bride_solar_control_rya1281 import inputs
    from pipeline.bride_solar_control import bound_fe_factory
    from pipeline.bride_reference import make_reference_context
    from pipeline.bride_top_cap import FiniteColumnCap
    from pipeline.bride_trace_fe_lines import trace_fe_factory
    folder,base,s,atom,contract=inputs('fine_T5772')
    provider,_=bound_fe_factory(folder,base,s,atom.abundance,column_contract=contract)
    if lines is not None:provider=trace_fe_factory(provider,lines,(777.,779.))
    return make_reference_context(s,atom,rays=2,qscale=1,fe_background=provider,
        extra_wavelengths_vac_nm=np.array([777.5,778.0,778.5]),
        experimental_column_cap=FiniteColumnCap(s['column_mass_kg_m2'][0]/s['rho_kg_m3'][0]))


def line():
    return dict(line_id='controlled',lambda0_vac_nm=778.,excitation_eV=4.,loggf=-2.,gamma_rad_s_inv=1e8)


@requires_lightweaver
def test_empty_delegate_identity_and_real_mixed_kirchhoff_scattering():
    import lightweaver.constants as c
    base=context(None);empty=context([]);mixed=context([line()])
    for k in ('chi','eta','sca'):np.testing.assert_array_equal(getattr(base.background,k),getattr(empty.background,k))
    dchi=np.array(mixed.background.chi)-np.array(base.background.chi)
    deta=np.array(mixed.background.eta)-np.array(base.background.eta)
    np.testing.assert_array_equal(mixed.background.sca,base.background.sca)
    w=np.array(base.spect.wavelength);inside=(w>=777)&(w<=779)
    assert np.all(dchi[inside]>0)
    np.testing.assert_array_equal(dchi[~inside],0)
    nu=c.CLight/(w[inside,None]*1e-9);T=base.atmos.pyAtmos.temperature
    B=2*c.HPlanck*nu**3/c.CLight**2/np.expm1(c.HPlanck*nu/(c.KBoltzmann*T))
    # Subtracting two large continuum coefficients loses relative precision in
    # weak wings. Bound the roundoff using the original operands, not the delta.
    scale=(abs(np.asarray(mixed.background.eta)[inside])+abs(np.asarray(base.background.eta)[inside])
        +(abs(np.asarray(mixed.background.chi)[inside])+abs(np.asarray(base.background.chi)[inside]))*B)
    assert np.all(abs(deta[inside]-dchi[inside]*B) <= 16*np.finfo(float).eps*scale)
    clone=copy.deepcopy(mixed)
    np.testing.assert_array_equal(clone.background.chi,mixed.background.chi)
    # The underlying complete-input refusal remains active through the wrapper.
    a=mixed.atmos.pyAtmos.temperature;old=a[0];a[0]=np.nextafter(old,np.inf)
    with pytest.raises(ValueError):mixed.background.update_background(mixed.atmos)


@requires_lightweaver
def test_bad_line_data_and_duplicate_or_outside_support_refuse():
    from pipeline.bride_trace_fe_lines import trace_fe_factory
    for key,value in [('gamma_rad_s_inv',0),('loggf',np.nan),('excitation_eV',-1)]:
        bad=line();bad[key]=value
        with pytest.raises(ValueError):trace_fe_factory(None,[bad],(777.,779.))
    with pytest.raises(ValueError):trace_fe_factory(None,[line(),line()],(777.,779.))
    with pytest.raises(ValueError):trace_fe_factory(None,[line()],(779.,780.))


@requires_lightweaver
def test_packet_revalidated_against_real_source_not_just_declared_hashes(tmp_path):
    import json
    from scripts.prepare_bride_fe_lines import source_packet,read_packet
    packet=source_packet();path=tmp_path/'lines.json';path.write_text(json.dumps(packet))
    assert read_packet(path)==packet
    packet['lines'][0]['loggf']+=.01
    path.write_text(json.dumps(packet))
    with pytest.raises(ValueError,match='fresh source/canonical'):read_packet(path)
