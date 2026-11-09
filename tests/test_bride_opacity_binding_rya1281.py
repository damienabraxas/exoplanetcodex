from tests.bride_test_requirements import requires_lightweaver
import copy
import json
import numpy as np
import pytest
from pipeline.bride_opacity_binding import local_binding, validate_binding, binding_digest, load_bound_coefficients


def structure():
    return dict(model_id='fixture', ndep=3, temperature_K=np.array([5000., 6000., 7000.]),
        electron_pressure_Pa=np.array([1., 2., 3.]), rho_kg_m3=np.array([1., 2., 3.]),
        height_m=np.array([2., 1., 0.]), logtau_ross=np.array([-2., -1., 0.]),
        logtau_500=np.array([-2., -1., 0.]), abundances_dex=np.arange(92, dtype=float).tolist(),
        hydrogen_mass_fraction=.7, vturb_m_s=1000., logg_cgs=4.438)


@requires_lightweaver
def test_exact_binding_and_explicit_missing_column():
    s = structure(); binding = local_binding(s, 8.7)
    validate_binding(binding, copy.deepcopy(s), 8.7)
    assert binding['column_mass_kg_m2'] is None
    assert len(binding['ne_m3']) == len(binding['nH_m3']) == 3
    with pytest.raises(ValueError, match='changed local'):
        validate_binding(binding, s, 8.71)


@requires_lightweaver
@pytest.mark.parametrize('key', ['temperature_K', 'electron_pressure_Pa', 'rho_kg_m3', 'height_m', 'logtau_ross', 'logtau_500', 'abundances_dex'])
def test_one_ulp_input_change_refused(key):
    s = structure(); binding = local_binding(s, 8.7); changed = copy.deepcopy(s)
    changed[key][0] = np.nextafter(changed[key][0], np.inf)
    with pytest.raises(ValueError, match='changed local'):
        validate_binding(binding, changed, 8.7)


@requires_lightweaver
@pytest.mark.parametrize('key', ['hydrogen_mass_fraction', 'vturb_m_s', 'logg_cgs'])
def test_scalar_identity_change_refused(key):
    s = structure(); binding = local_binding(s, 8.7); changed = copy.deepcopy(s)
    changed[key] = np.nextafter(changed[key], np.inf)
    with pytest.raises(ValueError, match='changed local'):
        validate_binding(binding, changed, 8.7)


@requires_lightweaver
def test_added_column_and_depth_cut_refused():
    s = structure(); binding = local_binding(s, 8.7)
    changed = copy.deepcopy(s); changed['column_mass_kg_m2'] = np.array([1., 2., 3.])
    with pytest.raises(ValueError, match='changed local'):
        validate_binding(binding, changed, 8.7)
    with pytest.raises(ValueError, match='field'):
        validate_binding(binding, dict(s, ndep=2), 8.7)


def test_unmodeled_velocity_is_not_silently_ignored():
    with pytest.raises(ValueError, match='Unmodeled'):
        local_binding(dict(structure(), vlos_m_s=np.zeros(3)), 8.7)


@requires_lightweaver
@pytest.mark.parametrize('tamper', ['array', 'source_report', 'binding_digest', 'provider_admission'])
def test_certificate_load_and_integrity_refusals(tmp_path, tamper):
    from pipeline.bride_reference import fingerprint
    s = structure(); binding = local_binding(s, 8.7)
    np.savez(tmp_path/'points.npz', wavelength_vac_nm=np.array([200., 300.]),
        requested_temperature_K=s['temperature_K'], requested_electron_pressure_Pa=s['electron_pressure_Pa'],
        local_rho_kg_m3=s['rho_kg_m3'], FeI_opacity_per_neutral_m2=np.ones((2, 3)))
    array = fingerprint(tmp_path/'points.npz')
    (tmp_path/'source.json').write_text(json.dumps(dict(status='DIAGNOSTIC_COMPLETE', cases=[dict(name='fixture', array=array)])))
    manifest = dict(schema='codex.bride_bound_point_coefficients/1', permitted_normalization='local_neutral_only',
        local_binding=binding, binding_sha256=binding_digest(binding), array=array,
        source_report=fingerprint(tmp_path/'source.json'), source_case='fixture',
        producer_native_abundances_dex=s['abundances_dex'],
        provider_admission=False, interpolated_cap_admission=False)
    (tmp_path/'manifest.json').write_text(json.dumps(manifest))
    load_bound_coefficients(tmp_path, s, 8.7)
    if tamper in ('array', 'source_report'):
        with (tmp_path/manifest[tamper]['name']).open('ab') as stream:
            stream.write(b' ')
    elif tamper == 'binding_digest':
        manifest['binding_sha256'] = 'wrong'
    else:
        manifest['provider_admission'] = True
    (tmp_path/'manifest.json').write_text(json.dumps(manifest))
    with pytest.raises(ValueError):
        load_bound_coefficients(tmp_path, s, 8.7)
