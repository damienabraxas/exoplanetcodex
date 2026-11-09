"""Joint native/held-out/solar pressure, EOS, geometry and column diagnostics."""
import argparse
import csv
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bindings', type=Path, required=True)
    parser.add_argument('--point-run', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    from config.constants import require_codex_path
    from pipeline.bride_reference import read_marcs_structure, fingerprint
    from pipeline.bride_pressure_checks import read_pressure_table, pressure_column_metrics, pressure_column_candidate
    from pipeline.bride_parameter_checks import interpolate_marcs_structure, marcs_filename
    from pipeline.bride_column_checks import integrated_column, candidate_column
    from pipeline.bride_opacity_binding import load_bound_coefficients
    from pipeline.bride_ts_budget import read_ts_opacity

    models, pressures, sources, rows, profiles, interpolation = {}, {}, [], [], [], []
    store = require_codex_path('grids.marcs_standard')
    for t in (5500, 5750, 6000):
        for g in (4., 4.5, 5.):
            path = store/marcs_filename(t, g)
            s = read_marcs_structure(path)
            models[t, g] = s
            pressures[t, g] = read_pressure_table(path, s)
            sources.append(fingerprint(path))

    def add(label, group, s, p, m, column, *, exact=None, producer=None):
        total = sum(p.values())
        gravity = 10**(s['logg_cgs']-2)
        integral = integrated_column(s['height_m'], s['rho_kg_m3'], m[0])
        xi_control = pressure_column_metrics(total+s['rho_kg_m3']*s['vturb_m_s']**2, m, gravity)
        row = dict(label=label, group=group, column=column,
            **pressure_column_metrics(total, m, gravity),
            rho_xi2_control_max_abs_adjacent_difference=xi_control['max_abs_adjacent_hydrostatic_difference'],
            rho_xi2_control_max_abs_top_anchored_difference=xi_control['max_abs_top_anchored_hydrostatic_difference'],
            max_fractional_mass_closure=float(abs(integral/m-1).max()),
            max_fractional_native_column_error=None,
            max_fractional_native_gas_pressure_error=None,
            max_fractional_producer_gas_pressure_difference=None,
            max_fractional_producer_rho_difference=None)
        if exact is not None:
            es, ep = exact
            row['max_fractional_native_column_error'] = float(abs(m/es['column_mass_kg_m2']-1).max())
            row['max_fractional_native_gas_pressure_error'] = float(abs(p['gas_pressure_Pa']/ep['gas_pressure_Pa']-1).max())
        if producer is not None:
            row['max_fractional_producer_gas_pressure_difference'] = float(abs(producer['gas_pressure_Pa']/p['gas_pressure_Pa']-1).max())
            row['max_fractional_producer_rho_difference'] = float(abs(producer['rho_kg_m3']/s['rho_kg_m3']-1).max())
        rows.append(row)
        for d in range(s['ndep']):
            profiles.append(dict(label=label, group=group, column=column, depth_1based=d+1,
                logtauRoss=s['logtau_ross'][d], gas_pressure_Pa=p['gas_pressure_Pa'][d],
                radiation_pressure_Pa=p['radiation_pressure_Pa'][d], turbulent_pressure_Pa=p['turbulent_pressure_Pa'][d],
                column_kg_m2=m[d], density_integrated_column_kg_m2=integral[d],
                producer_gas_pressure_Pa=None if producer is None else producer['gas_pressure_Pa'][d],
                local_rho_kg_m3=s['rho_kg_m3'][d], producer_rho_kg_m3=None if producer is None else producer['rho_kg_m3'][d]))

    for key, s in models.items():
        add(f'native_{key[0]}_{key[1]}', 'native', s, pressures[key], s['column_mass_kg_m2'], 'native_RHOX')

    def candidates(label, group, nodes, target, *, exact=None, producer=None, certificate=None):
        s, receipt = interpolate_marcs_structure(nodes, target, thermodynamic_scale='log_positive_layers')
        if 'column_mass_kg_m2' in s:
            raise ValueError('Unexpected column admission')
        indices, weights = receipt['simplex_source_indices'], receipt['weights']
        pp = [pressures[tuple(nodes[i][0])] for i in indices]
        p = dict(gas_pressure_Pa=np.exp(np.asarray(weights) @ np.log(np.stack([x['gas_pressure_Pa'] for x in pp]))),
            radiation_pressure_Pa=np.asarray(weights) @ np.stack([x['radiation_pressure_Pa'] for x in pp]),
            turbulent_pressure_Pa=np.asarray(weights) @ np.stack([x['turbulent_pressure_Pa'] for x in pp]))
        m = candidate_column([nodes[i][1]['column_mass_kg_m2'] for i in indices], weights, scale='log')
        if certificate is not None:
            load_bound_coefficients(args.bindings/label, s, certificate['local_binding']['scalars']['active_oxygen_abundance_dex'])
        interpolation.append(dict(label=label, target=list(target), receipt=receipt,
            pressure_convention='log gas pressure; linear radiation/turbulent pressure; same simplex weights, no EOS solve'))
        constructed = [('log_native_column', m), ('density_integrated_log_top', integrated_column(s['height_m'], s['rho_kg_m3'], m[0])),
                       ('pressure_difference_log_top', pressure_column_candidate(sum(p.values()), 10**(s['logg_cgs']-2), m[0]))]
        for column, value in constructed:
            add(label, group, s, p, value, column, exact=exact, producer=producer)

    outer = [((t, g), models[t, g]) for t in (5500, 6000) for g in (4., 5.)]
    for key in ((5500, 4.5), (5750, 4.), (5750, 4.5), (5750, 5.), (6000, 4.5)):
        candidates(f'heldout_{key[0]}_{key[1]}', 'heldout', outer, key, exact=(models[key], pressures[key]))

    point = json.loads(args.point_run.read_text())
    if point['status'] != 'DIAGNOSTIC_COMPLETE':
        raise ValueError('Incomplete producer point run')
    for case in point['cases']:
        name = case['name']
        certificate = json.loads((args.bindings/name/'manifest.json').read_text())
        if certificate['source_report']['sha256'] != fingerprint(args.point_run)['sha256']:
            raise ValueError('Binding and point report identities differ')
        entry = next(e for e in case['runs'] if e['kind'] == 'reporting')
        path = args.point_run.parent/(name+'_reporting')/'continuum.opa'
        if fingerprint(path) != entry['output']:
            raise ValueError('Changed producer pressure/EOS output')
        producer = read_ts_opacity(path)
        if name == 'native_identity':
            s = models[5750, 4.5]
            load_bound_coefficients(args.bindings/name, s, certificate['local_binding']['scalars']['active_oxygen_abundance_dex'])
            add(name, 'native_producer', s, pressures[5750, 4.5], s['column_mass_kg_m2'], 'native_RHOX', producer=producer)
        else:
            fields = certificate['local_binding']
            # T label is from the existing case identifier, checked against its
            # source binding after rebuilding the same interpolated fields.
            teff = int(name.removeprefix('fine_T'))
            nodes = [((t, g), models[t, g]) for t in (5750, 6000) for g in (4., 4.5)]
            candidates(name, 'solar', nodes, (teff, fields['scalars']['logg_cgs']), producer=producer, certificate=certificate)

    args.output.mkdir(parents=True, exist_ok=False)
    result = dict(schema='codex.bride_pressure_column_audit/1', scientific_admission='HOLD', provider_admission=False,
        generator=fingerprint(__file__), pressure_checks=fingerprint(ROOT/'pipeline/bride_pressure_checks.py'),
        model_sources=sources, source_point_run=fingerprint(args.point_run), source_binding_report=fingerprint(args.bindings/'report.json'),
        primary_reference=fingerprint(require_codex_path('refs.library')/'MARCS Stellar Model.pdf'),
        method='Gustafsson2008 p953 Eqs1-4, plane-parallel pressure differences; native printed Pg/Prad/Pturb. Interpolated Pg logarithmic, Prad/Pturb linear. Density-integrated and pressure-derived columns enforce separate identities by construction; top column always an explicit log-native anchor. Additional beta1 rho*xi^2 is a declared nonfitted control only, NOT identification of spectral microturbulence with support velocity (paper p954 sets grid support vt=0). No native truth at actual solar targets, threshold, EOS repair, force/flux solve or boundary admission.',
        interpolation=interpolation, cases=rows)
    for name, data in [('summary.csv', rows), ('profiles.csv', profiles)]:
        with (args.output/name).open('x') as stream:
            writer = csv.DictWriter(stream, fieldnames=list(data[0]), lineterminator='\n')
            writer.writeheader(); writer.writerows(data)
        result[name] = fingerprint(args.output/name)
    (args.output/'report.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps(dict(cases=len(rows), native_models=len(models), source_fields_unchanged=True)))


if __name__ == '__main__':
    main()
