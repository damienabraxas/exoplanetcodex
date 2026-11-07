"""Create immutable exact local-input certificates, without enabling a provider."""
import argparse
import copy
import csv
import json
from pathlib import Path
import shutil
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    from config.constants import require_codex_path, get_star_params
    from pipeline.bride_reference import fingerprint, read_marcs_structure
    from pipeline.bride_parameter_checks import marcs_filename, interpolate_marcs_structure
    from pipeline.bride_opacity_binding import local_binding, binding_digest, load_bound_coefficients
    from pipeline.bride_fe_background import load_packet

    run = json.loads(args.run.read_text())
    if run['status'] != 'DIAGNOSTIC_COMPLETE':
        raise ValueError('Incomplete point coefficient run')
    plan_path = ROOT/'config/experiments/rya1281_mac_fine_temperature.json'
    if fingerprint(plan_path) != run['plan']:
        raise ValueError('Changed input plan')
    plan = json.loads(plan_path.read_text())
    store = require_codex_path('grids.marcs_standard')
    for receipt in run['source_models']:
        if fingerprint(store/receipt['name']) != receipt:
            raise ValueError('Changed native source model')
    native = read_marcs_structure(store/run['source_models'][0]['name'])
    solar_case = next(c for c in plan['cases'] if c['name'].startswith('fine_T') and c['teff'] == get_star_params('solar')['teff'])
    # Same declared fixed diagnostic O as the existing solar plan; the producer
    # composition is recorded separately and is not silently changed here.
    cases = {'native_identity': (native, solar_case['abundance'])}
    for case in [c for c in plan['cases'] if c['name'].startswith('fine_T')]:
        nodes = [((t, g), read_marcs_structure(store/marcs_filename(t, g))) for t in case['bracket_teff'] for g in case['bracket_logg']]
        s, receipt = interpolate_marcs_structure(nodes, (case['teff'], case['logg']), thermodynamic_scale='log_positive_layers')
        cases[case['name']] = (s, case['abundance'])
    if set(cases) != {c['name'] for c in run['cases']}:
        raise ValueError('Producer/local case identities differ')
    args.output.mkdir(parents=True, exist_ok=False)
    result = dict(schema='codex.bride_bound_coefficients_audit/1', scientific_admission='HOLD',
        generator=fingerprint(__file__), binding_code=fingerprint(ROOT/'pipeline/bride_opacity_binding.py'),
        source_run=fingerprint(args.run), plan=run['plan'], cases=[], refusals=[])
    checks = []
    for source_case in run['cases']:
        name = source_case['name']
        structure, oxygen = cases[name]
        source_array = args.run.parent/source_case['array']['name']
        if fingerprint(source_array) != source_case['array']:
            raise ValueError('Changed producer point array')
        directory = args.output/name
        directory.mkdir()
        shutil.copyfile(source_array, directory/source_array.name)
        shutil.copyfile(args.run, directory/'source_run.json')
        binding = local_binding(structure, oxygen)
        manifest = dict(schema='codex.bride_bound_point_coefficients/1', scientific_admission='HOLD',
            generator=fingerprint(__file__), binding_code=result['binding_code'],
            local_binding=binding, binding_sha256=binding_digest(binding), array=fingerprint(directory/source_array.name),
            source_report=fingerprint(directory/'source_run.json'), source_case=name,
            producer_native_abundances_dex=native['abundances_dex'], permitted_normalization='local_neutral_only',
            provider_admission=False, interpolated_cap_admission=False,
            interpretation='Exact static local-input certificate, not an enabled background provider. Producer raw mass coefficients remain diagnostic; do not replace local rho. Missing column is explicit, and adding a constructed column invalidates this binding.')
        (directory/'manifest.json').write_text(json.dumps(manifest, indent=2, allow_nan=False)+'\n')
        loaded, checked = load_bound_coefficients(directory, structure, oxygen)
        result['cases'].append(dict(name=name, manifest=fingerprint(directory/'manifest.json'),
            binding_sha256=manifest['binding_sha256'], ndep=structure['ndep'], nwave=len(loaded['wavelength_vac_nm']),
            local_column_present=binding['column_mass_kg_m2'] is not None, exact_identity_pass=True))

        def refuse(label, candidate, abundance=oxygen):
            try:
                load_bound_coefficients(directory, candidate, abundance)
            except ValueError as exc:
                checks.append(dict(case=name, probe=label, refused=True, reason=str(exc)))
            else:
                raise ValueError('Invalid binding accepted: '+label)

        for key in ('temperature_K', 'electron_pressure_Pa', 'rho_kg_m3', 'height_m', 'logtau_ross', 'logtau_500'):
            altered = copy.deepcopy(structure)
            altered[key][0] = np.nextafter(altered[key][0], np.inf)
            refuse('one_ulp_'+key, altered)
        for key in ('vturb_m_s', 'hydrogen_mass_fraction', 'logg_cgs'):
            altered = copy.deepcopy(structure)
            altered[key] = np.nextafter(altered[key], np.inf)
            refuse('one_ulp_'+key, altered)
        altered = copy.deepcopy(structure); altered['abundances_dex'][25] += .01
        refuse('changed_Fe_composition', altered)
        refuse('changed_active_O_abundance', structure, oxygen+.01)
        altered = copy.deepcopy(structure)
        if 'column_mass_kg_m2' in altered:
            altered['column_mass_kg_m2'][0] *= 1.001
        else:
            altered['column_mass_kg_m2'] = np.arange(1, structure['ndep']+1, dtype=float)
        refuse('changed_or_added_column', altered)
        altered = copy.deepcopy(structure); altered['ndep'] -= 1
        refuse('depth_cut', altered)
        altered = copy.deepcopy(structure); altered['vlos_m_s'] = np.zeros(structure['ndep'])
        refuse('unmodeled_velocity', altered)
        try:
            load_packet(directory)
        except ValueError as exc:
            checks.append(dict(case=name, probe='legacy_provider_schema_refusal', refused=True, reason=str(exc)))
        else:
            raise ValueError('Legacy Fe packet loader unexpectedly accepts certificate')
    # Valid inputs from another 10K point must not use a neighboring certificate.
    for left, right in zip(run['cases'][1:-1], run['cases'][2:]):
        try:
            load_bound_coefficients(args.output/left['name'], *cases[right['name']])
        except ValueError as exc:
            checks.append(dict(case=left['name'], probe='neighbor_'+right['name'], refused=True, reason=str(exc)))
        else:
            raise ValueError('Neighbor input accepted')
    with (args.output/'checks.csv').open('x') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(checks[0]), lineterminator='\n')
        writer.writeheader(); writer.writerows(checks)
    result['refusals'] = checks
    result['checks'] = fingerprint(args.output/'checks.csv')
    result['provider_admission'] = False
    (args.output/'report.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps(dict(exact_identities=len(result['cases']), refused=len(checks), provider_admission=False)))


if __name__ == '__main__':
    main()
