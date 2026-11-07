"""Solar 10K input sensitivity to interpolation operator/cell; no accuracy claim."""
import argparse
import csv
import json
from pathlib import Path
import subprocess
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, default=ROOT/'config/experiments/rya1281_mac_fine_temperature.json')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    from config.constants import require_codex_path, get_star_params
    from pipeline.bride_reference import read_marcs_structure, fingerprint
    from pipeline.bride_parameter_checks import marcs_filename, interpolate_marcs_structure
    from pipeline.bride_rectangular_checks import rectangular_structure
    from pipeline.bride_column_checks import integrated_column, candidate_column
    from pipeline.bride_input_sensitivity import field_difference, diagnostic_depth_masks

    plan = json.loads(args.plan.read_text())
    cases = sorted([c for c in plan['cases'] if c['name'].startswith('fine_T')], key=lambda c: c['teff'])
    solar = get_star_params('solar')
    temperatures = [c['teff'] for c in cases]
    if len(cases) != 5 or not np.all(np.diff(temperatures) == 10) or temperatures[2] != solar['teff']:
        raise ValueError('Expected existing five-point 10K sequence centered on canonical Sun')
    if any(c['logg'] != solar['logg'] or c['bracket_teff'] != [5750, 6000] or c['bracket_logg'] != [4., 4.5] for c in cases):
        raise ValueError('Changed canonical gravity or existing nearest-cell experiment contract')
    cells = dict(nearest=([5750, 6000], [4., 4.5]), wider_temperature=([5500, 6000], [4., 4.5]),
                 wider_both=([5500, 6000], [4., 5.]))
    store = require_codex_path('grids.marcs_standard')
    models, sources = {}, []
    for t, g in sorted({(t, g) for ts, gs in cells.values() for t in ts for g in gs}):
        path = store/marcs_filename(t, g)
        models[t, g] = read_marcs_structure(path)
        sources.append(dict(teff=t, logg=g, source=fingerprint(path)))
    native_tau = next(iter(models.values()))['logtau_ross']
    if any(not np.array_equal(s['logtau_ross'], native_tau) for s in models.values()):
        raise ValueError('Native source row coordinates differ')
    scopes = diagnostic_depth_masks(native_tau)
    profiles, columns, comparisons, receipts, structures = [], [], [], [], {}
    modes = ('linear', 'log_positive', 'log_positive_layers')
    for cell, (ts, gs) in cells.items():
        nodes = [((t, g), models[t, g]) for t in ts for g in gs]
        for method in ('simplex', 'bilinear'):
            for mode in modes:
                for t in temperatures:
                    interpolate = interpolate_marcs_structure if method == 'simplex' else rectangular_structure
                    s, receipt = interpolate(nodes, (t, solar['logg']), thermodynamic_scale=mode)
                    if 'column_mass_kg_m2' in s:
                        raise ValueError('Diagnostic must withhold production column metadata')
                    if not np.allclose(s['logtau_ross'], native_tau, rtol=0, atol=1e-12):
                        raise ValueError('Candidate changes native depth coordinates')
                    structures[cell, method, mode, t] = s
                    identity = dict(cell=cell, method=method, thermodynamic_scale=mode, teff=t, logg=solar['logg'])
                    receipts.append(dict(**identity, receipt=receipt, source_nodes=[list(p) for p, s in nodes]))
                    for d in range(s['ndep']):
                        profiles.append(dict(**identity, depth_1based=d+1, logtauRoss=native_tau[d],
                            temperature_K=s['temperature_K'][d], electron_pressure_Pa=s['electron_pressure_Pa'][d],
                            rho_kg_m3=s['rho_kg_m3'][d], height_m=s['height_m'][d]))
                    indices = receipt['simplex_source_indices'] if method == 'simplex' else receipt['source_indices']
                    for column in ('linear', 'log', 'density_integrated_log_top'):
                        m = candidate_column([nodes[i][1]['column_mass_kg_m2'] for i in indices], receipt['weights'],
                                             scale='linear' if column == 'linear' else 'log')
                        if column == 'density_integrated_log_top':
                            m = integrated_column(s['height_m'], s['rho_kg_m3'], m[0])
                        integral = integrated_column(s['height_m'], s['rho_kg_m3'], m[0])
                        for scope, mask in scopes:
                            columns.append(dict(**identity, column_scale=column, scope=scope,
                                max_fractional_mass_closure=float(abs(integral[mask]/m[mask]-1).max())))

    def compare(kind, candidate_key, baseline_key):
        candidate, baseline = structures[candidate_key], structures[baseline_key]
        for scope, mask in scopes:
            comparisons.append(dict(kind=kind, cell=candidate_key[0], method=candidate_key[1],
                thermodynamic_scale=candidate_key[2], teff=candidate_key[3], baseline_cell=baseline_key[0],
                baseline_method=baseline_key[1], baseline_teff=baseline_key[3], scope=scope,
                **field_difference(candidate, baseline, mask)))

    for cell in cells:
        for mode in modes:
            for t in temperatures:
                compare('method_spread', (cell, 'bilinear', mode, t), (cell, 'simplex', mode, t))
            for method in ('simplex', 'bilinear'):
                for low, high in zip(temperatures[:-1], temperatures[1:]):
                    compare('10K_response', (cell, method, mode, high), (cell, method, mode, low))
                if cell != 'nearest':
                    for t in temperatures:
                        compare('cell_spread', (cell, method, mode, t), ('nearest', method, mode, t))
    args.output.mkdir(parents=True, exist_ok=False)
    result = dict(schema='codex.bride_solar_input_sensitivity/1', scientific_admission='HOLD', production_column_admission=False,
        code_commit=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
        plan=fingerprint(args.plan), solar_parameters=solar, stars_config=fingerprint(ROOT/'config/stars.yaml'), sources=sources,
        code=[fingerprint(ROOT/p) for p in ('scripts/audit_bride_solar_input_sensitivity.py', 'pipeline/bride_input_sensitivity.py',
            'pipeline/bride_rectangular_checks.py', 'pipeline/bride_column_checks.py', 'pipeline/bride_parameter_checks.py', 'pipeline/bride_reference.py')],
        depth_scopes=[dict(name=name, source_rows_1based=(np.flatnonzero(mask)+1).tolist()) for name, mask in scopes],
        interpretation='90 input field candidates, full native depth extent (previous solves dropped one bottom row). Method/cell spread and internal column closure are NOT independent accuracy, EOS/equilibrium validation, line contribution or abundance uncertainty. No native reference at solar targets. logtau[-3,0] uses fixed native source rows, not rounded interpolated tau, and is a declared diagnostic subset only. All three field conventions retained. No RT/SE or new opacity packet.', receipts=receipts)
    for name, data in [('profiles.csv', profiles), ('column_checks.csv', columns), ('comparisons.csv', comparisons)]:
        with (args.output/name).open('x') as stream:
            writer = csv.DictWriter(stream, fieldnames=list(data[0]), lineterminator='\n')
            writer.writeheader()
            writer.writerows(data)
        result[name] = fingerprint(args.output/name)
    (args.output/'report.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps(dict(candidates=len(receipts), comparisons=len(comparisons), sources=len(sources))))


if __name__ == '__main__':
    main()
