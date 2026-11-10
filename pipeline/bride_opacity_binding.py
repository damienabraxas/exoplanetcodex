"""Exact local-input certificates for experimental TS coefficients, not providers."""
import hashlib
import json
from pathlib import Path

import numpy as np


def binding_digest(binding):
    return hashlib.sha256(json.dumps(binding, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()


def local_binding(structure, oxygen_abundance_dex):
    """Bind complete static local inputs, derived number densities and O abundance.

    No tolerance, source-column construction, depth interpolation or EOS repair.
    Missing column metadata is explicit and cannot later become cap permission.
    """
    from pipeline.bride_reference import marcs_number_densities
    names = ('temperature_K', 'electron_pressure_Pa', 'rho_kg_m3', 'height_m', 'logtau_ross', 'logtau_500')
    allowed = set(names) | {'model_id', 'ndep', 'hydrogen_mass_fraction', 'vturb_m_s', 'logg_cgs', 'abundances_dex', 'column_mass_kg_m2'}
    if set(structure)-allowed:
        raise ValueError('Unmodeled fields in static local-input binding')
    n = structure['ndep']
    if not isinstance(n, int) or isinstance(n, bool) or n < 2:
        raise ValueError('Invalid bound depth count')
    fields = {}
    for name in names:
        values = np.asarray(structure[name], dtype=float)
        if values.shape != (n,) or not np.isfinite(values).all():
            raise ValueError('Invalid bound field: '+name)
        fields[name] = values.tolist()
    if any(np.any(np.asarray(fields[k]) <= 0) for k in names[:3]) or np.any(np.diff(fields['height_m']) >= 0):
        raise ValueError('Invalid bound thermodynamic/geometry fields')
    if any(np.any(np.diff(fields[k]) <= 0) for k in ('logtau_ross', 'logtau_500')):
        raise ValueError('Invalid native row ordering')
    scalars = {k: float(structure[k]) for k in ('hydrogen_mass_fraction', 'vturb_m_s', 'logg_cgs')}
    scalars['active_oxygen_abundance_dex'] = float(oxygen_abundance_dex)
    if not all(np.isfinite(v) for v in scalars.values()) or not 0 < scalars['hydrogen_mass_fraction'] < 1 or scalars['vturb_m_s'] < 0:
        raise ValueError('Invalid bound scalar')
    abundances = np.asarray(structure['abundances_dex'], dtype=float)
    if abundances.shape != (92,) or not np.isfinite(abundances).all():
        raise ValueError('Require complete native abundance vector')
    normalized = dict(structure, **{k: np.asarray(fields[k]) for k in names[:3]})
    ne, nH = marcs_number_densities(normalized)
    if not np.isfinite(ne).all() or not np.isfinite(nH).all() or np.any(ne <= 0) or np.any(nH <= 0):
        raise ValueError('Invalid derived number densities')
    column = structure.get('column_mass_kg_m2')
    if column is not None:
        column = np.asarray(column, dtype=float)
        if column.shape != (n,) or not np.isfinite(column).all() or np.any(column <= 0) or np.any(np.diff(column) <= 0):
            raise ValueError('Invalid bound column metadata')
    return dict(schema='codex.bride_local_input_binding/1', model_id=str(structure['model_id']), ndep=n,
        fields=fields, scalars=scalars, native_abundances_dex=abundances.tolist(),
        ne_m3=ne.tolist(), nH_m3=nH.tolist(), static_vlos_m_s=[0.]*n,
        column_mass_kg_m2=None if column is None else column.tolist())


def validate_binding(expected, structure, oxygen_abundance_dex):
    if binding_digest(expected) != binding_digest(local_binding(structure, oxygen_abundance_dex)):
        raise ValueError('Bound opacity refuses changed local inputs, composition, abundance or depth/column identity')


def load_bound_coefficients(directory, structure, oxygen_abundance_dex):
    """Load a checked certificate; deliberately not a BackgroundProvider factory."""
    from pipeline.bride_reference import fingerprint
    directory = Path(directory)
    manifest = json.loads((directory/'manifest.json').read_text())
    if manifest['schema'] != 'codex.bride_bound_point_coefficients/1' or manifest['permitted_normalization'] != 'local_neutral_only':
        raise ValueError('Wrong bound coefficient schema or normalization')
    if manifest['provider_admission'] is not False or manifest['interpolated_cap_admission'] is not False:
        raise ValueError('Certificate cannot enable a provider or interpolated cap')
    if binding_digest(manifest['local_binding']) != manifest['binding_sha256']:
        raise ValueError('Changed binding digest')
    validate_binding(manifest['local_binding'], structure, oxygen_abundance_dex)
    if not np.array_equal(manifest['producer_native_abundances_dex'], structure['abundances_dex']):
        raise ValueError('Producer/native composition differs from local structure')
    for key in ('array', 'source_report'):
        path = directory/manifest[key]['name']
        if fingerprint(path) != manifest[key]:
            raise ValueError('Changed bound coefficient source: '+key)
    report = json.loads((directory/manifest['source_report']['name']).read_text())
    source_case = next(c for c in report['cases'] if c['name'] == manifest['source_case'])
    if report['status'] != 'DIAGNOSTIC_COMPLETE' or source_case['array'] != manifest['array']:
        raise ValueError('Bound array does not match completed producer case')
    with np.load(directory/manifest['array']['name'], allow_pickle=False) as loaded:
        data = {k: loaded[k] for k in loaded.files}
    if any(not np.isfinite(v).all() for v in data.values()):
        raise ValueError('Nonfinite producer coefficient evidence')
    for key, local in (('requested_temperature_K', 'temperature_K'), ('requested_electron_pressure_Pa', 'electron_pressure_Pa'), ('local_rho_kg_m3', 'rho_kg_m3')):
        if not np.array_equal(data[key], structure[local]):
            raise ValueError('Producer requests do not match bound local fields')
    wave = data['wavelength_vac_nm']
    if wave.ndim != 1 or len(wave) < 2 or not np.isfinite(wave).all() or np.any(wave <= 0) or np.any(np.diff(wave) <= 0):
        raise ValueError('Invalid bound wavelength coordinate')
    coefficient = data['FeI_opacity_per_neutral_m2']
    if coefficient.shape != (len(wave), structure['ndep']) or not np.isfinite(coefficient).all() or np.any(coefficient < 0):
        raise ValueError('Invalid bound neutral Fe coefficient')
    # Producer rho is separate evidence, not an implicit local-density replacement.
    return data, manifest
