"""Compare simplex/bilinear atmosphere diagnostics on three overlapping native cells."""
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
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    from config.constants import require_codex_path
    from pipeline.bride_reference import read_marcs_structure, fingerprint
    from pipeline.bride_parameter_checks import marcs_filename, interpolate_marcs_structure
    from pipeline.bride_rectangular_checks import rectangular_structure
    from pipeline.bride_column_checks import integrated_column, candidate_column

    store = require_codex_path('grids.marcs_standard')
    models, sources, native, rows, profiles = {}, [], [], [], []
    for teff in (5250, 5500, 5750, 6000, 6250):
        for logg in (4., 4.5, 5.):
            path = store/marcs_filename(teff, logg)
            s = read_marcs_structure(path)
            models[teff, logg] = s
            sources.append(dict(teff=teff, logg=logg, source=fingerprint(path)))
            m = s['column_mass_kg_m2']
            integral = integrated_column(s['height_m'], s['rho_kg_m3'], m[0])
            native.append(dict(teff=teff, logg=logg, max_fractional_mass_closure=float(abs(integral/m-1).max())))
    for low, high in ((5250, 5750), (5500, 6000), (5750, 6250)):
        middle = (low+high)//2
        nodes = [((t, g), models[t, g]) for t in (low, high) for g in (4., 5.)]
        targets = [(low, 4.5), (middle, 4.), (middle, 4.5), (middle, 5.), (high, 4.5)]
        for target in targets:
            exact = models[target]
            for method in ('simplex', 'bilinear'):
                for thermo in ('linear', 'log_positive', 'log_positive_layers'):
                    interpolate = interpolate_marcs_structure if method == 'simplex' else rectangular_structure
                    s, receipt = interpolate(nodes, target, thermodynamic_scale=thermo)
                    if 'column_mass_kg_m2' in s:
                        raise ValueError('Diagnostic unexpectedly admits column metadata')
                    indices = receipt['simplex_source_indices'] if method == 'simplex' else receipt['source_indices']
                    weights = receipt['weights']
                    identity = dict(cell_low_K=low, cell_high_K=high, teff=target[0], logg=target[1],
                                    method=method, thermodynamic_scale=thermo)
                    native_ne = exact['electron_pressure_Pa']/exact['temperature_K']
                    ne = s['electron_pressure_Pa']/s['temperature_K']
                    for d in range(s['ndep']):
                        profiles.append(dict(**identity, depth_1based=d+1, logtauRoss=exact['logtau_ross'][d],
                            native_temperature_K=exact['temperature_K'][d], candidate_temperature_K=s['temperature_K'][d],
                            native_electron_pressure_Pa=exact['electron_pressure_Pa'][d], candidate_electron_pressure_Pa=s['electron_pressure_Pa'][d],
                            native_rho_kg_m3=exact['rho_kg_m3'][d], candidate_rho_kg_m3=s['rho_kg_m3'][d],
                            native_height_m=exact['height_m'][d], candidate_height_m=s['height_m'][d]))
                    for column in ('linear', 'log', 'density_integrated_log_top'):
                        m = candidate_column([nodes[i][1]['column_mass_kg_m2'] for i in indices], weights,
                                             scale='linear' if column == 'linear' else 'log')
                        if column == 'density_integrated_log_top':
                            m = integrated_column(s['height_m'], s['rho_kg_m3'], m[0])
                        integral = integrated_column(s['height_m'], s['rho_kg_m3'], m[0])
                        rows.append(dict(**identity, column_scale=column,
                            max_temperature_error_K=float(abs(s['temperature_K']-exact['temperature_K']).max()),
                            max_fractional_Pe_error=float(abs(s['electron_pressure_Pa']/exact['electron_pressure_Pa']-1).max()),
                            max_fractional_ne_error=float(abs(ne/native_ne-1).max()),
                            max_fractional_rho_error=float(abs(s['rho_kg_m3']/exact['rho_kg_m3']-1).max()),
                            max_fractional_column_error=float(abs(m/exact['column_mass_kg_m2']-1).max()),
                            max_fractional_mass_closure=float(abs(integral/m-1).max()),
                            weights=json.dumps(weights), source_indices=json.dumps(indices)))
    args.output.mkdir(parents=True, exist_ok=False)
    result = dict(schema='codex.bride_rectangular_interpolation_audit/1', scientific_admission='HOLD',
        production_column_admission=False, sources=sources, native_closure=native,
        method='Three overlapping 500K by 1dex cells, five held-out targets per cell; not fifteen independent targets. Both operators, all three thermodynamic and column conventions retained. No RT/SE, EOS/equilibrium solve or threshold derived from these comparisons. ne comparison uses Pe/T; kB cancels.',
        code=[fingerprint(ROOT/p) for p in ('scripts/audit_bride_rectangular_interpolation.py',
            'pipeline/bride_rectangular_checks.py', 'pipeline/bride_column_checks.py', 'pipeline/bride_parameter_checks.py', 'pipeline/bride_reference.py')], cases=rows)
    for name, data in [('summary.csv', rows), ('profiles.csv', profiles)]:
        with (args.output/name).open('x') as stream:
            writer = csv.DictWriter(stream, fieldnames=list(data[0]), lineterminator='\n')
            writer.writeheader()
            writer.writerows(data)
        result[name] = fingerprint(args.output/name)
    (args.output/'report.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps(dict(cases=len(rows), native_models=len(sources), held_out_cell_targets=15)))


if __name__ == '__main__':
    main()
