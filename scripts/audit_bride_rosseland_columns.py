"""Isolate official calcrhox; analytic, native, held-out and solar column checks."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import platform
import shutil
import subprocess
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bindings', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if platform.system() != 'Darwin':
        raise RuntimeError('Mac-only column diagnostic')
    from config.constants import require_codex_path
    from pipeline.bride_reference import fingerprint, read_marcs_structure
    from pipeline.bride_parameter_checks import marcs_filename, interpolate_marcs_structure
    from pipeline.bride_rosseland_column import read_rosseland_opacity, run_column_kernel, log_opacity_integral
    from pipeline.bride_column_checks import candidate_column, integrated_column
    from pipeline.bride_pressure_checks import read_pressure_table, pressure_column_metrics
    from pipeline.bride_opacity_binding import load_bound_coefficients

    source = require_codex_path('experiments.bride_marcs_interpolator')
    receipt = json.loads((source/'manifest.json').read_text())
    entry = next(r for r in receipt['files'] if r['relative_path'] == 'interpol_marcs/interpol_modeles.f')
    path = source/entry['relative_path']
    if fingerprint(path) != {k: v for k, v in entry.items() if k != 'relative_path'}:
        raise ValueError('Changed official interpolation source')
    text = path.read_text()
    start, stop = '      subroutine calcrhox(', '      subroutine calc_error('
    if text.count(start) != 1 or text.count(stop) != 1:
        raise ValueError('Changed native kernel extraction sites')
    fragment = text[text.index(start):text.index(stop)]
    if fragment.lower().count('subroutine ') != 2 or 'real function rinteg(' not in fragment:
        raise ValueError('Unknown official kernel fragment')
    args.output.mkdir(parents=True, exist_ok=False)
    build = args.output/'build'; build.mkdir()
    (build/'official_kernel.f').write_text(fragment)
    driver = """      PROGRAM BRIDE_COLUMN
      INTEGER N,I,IOS
      REAL TAU(400),KAPPA(400),MASS(400)
      READ(*,*,IOSTAT=IOS) N
      IF(IOS.NE.0.OR.N.LT.4.OR.N.GT.400) STOP 61
      DO I=1,N
        READ(*,*,IOSTAT=IOS) TAU(I),KAPPA(I)
        IF(IOS.NE.0) STOP 62
      ENDDO
      CALL CALCRHOX(TAU,KAPPA,N,MASS)
      DO I=1,N
        WRITE(*,'(ES24.16E3)') MASS(I)
      ENDDO
      END
"""
    (build/'driver.f').write_text(driver)
    compiler = shutil.which('gfortran')
    if compiler is None:
        raise RuntimeError('Existing gfortran is required')
    command = [compiler, '-O0', '-g', '-std=legacy', '-fcheck=all', '-fbacktrace', 'driver.f', 'official_kernel.f', '-o', 'column_kernel']
    compiled = subprocess.run(command, cwd=build, capture_output=True, text=True, timeout=60)
    (build/'compiler.log').write_text(compiled.stdout+compiled.stderr)
    if compiled.returncode:
        raise RuntimeError('Isolated official column kernel failed to compile')
    binary = build/'column_kernel'
    analytic = []
    for n in (4, 9, 56):
        for law in ('constant', 'linear_inverse_opacity'):
            tau = np.linspace(.01, 1., n)
            inverse = np.full(n, 2.) if law == 'constant' else 2.+3*tau
            opacity = .1/inverse  # f is cgs reciprocal opacity
            actual = run_column_kernel(binary, np.log10(tau), opacity)
            expected = 10*(2*tau if law == 'constant' else tau[0]*inverse[0]+2*(tau-tau[0])+1.5*(tau**2-tau[0]**2))
            error = float(abs(actual/expected-1).max())
            # Prospective conservative REAL operation-count budget, not fitted
            # from output or a stellar-structure admission tolerance.
            operations = 128*n
            eps = np.finfo(np.float32).eps
            budget = float(operations*eps/(1-operations*eps))
            if error > budget:
                raise ValueError('Analytic column kernel exceeds declared REAL budget')
            analytic.append(dict(law=law, ndepth=n, max_fractional_error=error, declared_REAL_budget=budget))

    store = require_codex_path('grids.marcs_standard')
    models, opacities, pressures, sources, rows, profiles, weights = {}, {}, {}, [], [], [], []
    for t in (5500, 5750, 6000):
        for g in (4., 4.5, 5.):
            path = store/marcs_filename(t, g)
            s = read_marcs_structure(path)
            models[t, g], opacities[t, g], pressures[t, g] = s, read_rosseland_opacity(path, s), read_pressure_table(path, s)
            sources.append(fingerprint(path))

    def add(label, group, s, opacity, pressure, native_column=None, log_column=None, native_structure=None):
        original = run_column_kernel(binary, s['logtau_ross'], opacity)
        anchor = log_column if log_column is not None else native_column
        methods = [('official_native_mesh', original)]
        for factor in (2, 4):
            tau = s['logtau_ross']
            refined_tau = np.r_[np.concatenate([np.linspace(tau[i], tau[i+1], factor, endpoint=False) for i in range(len(tau)-1)]), tau[-1]]
            refined_opacity = np.exp(np.interp(refined_tau, tau, np.log(opacity)))
            refined_opacity[::factor] = opacity
            mass = run_column_kernel(binary, refined_tau, refined_opacity)
            methods.append((f'official_refined_{factor}x_same_log_opacity', mass[::factor]))
        methods.append(('analytic_piecewise_log_opacity', log_opacity_integral(s['logtau_ross'], opacity, 10**s['logtau_ross'][0]/opacity[0])))
        for method, raw in methods:
            for convention, m in [('official_tau0_over_kappa0', raw), ('external_native_or_log_top', raw+(anchor[0]-raw[0]))]:
                integral = integrated_column(s['height_m'], s['rho_kg_m3'], m[0])
                row = dict(label=label, group=group, integrator=method, anchor=convention,
                    effective_top_cap_length_m=float(m[0]/s['rho_kg_m3'][0]),
                    top_fractional_native_cap_length_error=None if native_structure is None else float((m[0]/s['rho_kg_m3'][0])/(native_structure['column_mass_kg_m2'][0]/native_structure['rho_kg_m3'][0])-1),
                    max_fractional_native_column_error=None if native_column is None else float(abs(m/native_column-1).max()),
                    top_fractional_native_column_error=None if native_column is None else float(m[0]/native_column[0]-1),
                    max_fractional_log_column_difference=None if log_column is None else float(abs(m/log_column-1).max()),
                    max_fractional_mass_closure=float(abs(integral/m-1).max()),
                    **pressure_column_metrics(sum(pressure.values()), m, 10**(s['logg_cgs']-2)))
                rows.append(row)
                for d in range(s['ndep']):
                    profiles.append(dict(label=label, group=group, integrator=method, anchor=convention, depth_1based=d+1,
                        logtauRoss=s['logtau_ross'][d], opacity_m2_kg=opacity[d], column_kg_m2=m[d], density_integral_kg_m2=integral[d],
                        native_column_kg_m2=None if native_column is None else native_column[d],
                        interpolated_log_column_kg_m2=None if log_column is None else log_column[d]))

    for key, s in models.items():
        add(f'native_{key[0]}_{key[1]}', 'native', s, opacities[key], pressures[key], native_column=s['column_mass_kg_m2'], native_structure=s)

    def interpolate(label, group, nodes, target, exact=None, oxygen=None):
        s, receipt = interpolate_marcs_structure(nodes, target, thermodynamic_scale='log_positive_layers')
        indices, w = receipt['simplex_source_indices'], np.asarray(receipt['weights'])
        keys = [tuple(nodes[i][0]) for i in indices]
        opacity = np.exp(w @ np.log(np.stack([opacities[k] for k in keys])))
        pressure = {field: (np.exp(w @ np.log(np.stack([pressures[k][field] for k in keys]))) if field == 'gas_pressure_Pa' else w @ np.stack([pressures[k][field] for k in keys])) for field in pressures[keys[0]]}
        column = candidate_column([models[k]['column_mass_kg_m2'] for k in keys], w, scale='log')
        if oxygen is not None:
            load_bound_coefficients(args.bindings/label, s, oxygen)
        weights.append(dict(label=label, receipt=receipt, opacity_convention='log Rosseland opacity on common native rows; no new opacity/equilibrium solve'))
        add(label, group, s, opacity, pressure, native_column=None if exact is None else exact['column_mass_kg_m2'], log_column=column, native_structure=exact)

    nodes = [((t, g), models[t, g]) for t in (5500, 6000) for g in (4., 5.)]
    for key in ((5500, 4.5), (5750, 4.), (5750, 4.5), (5750, 5.), (6000, 4.5)):
        interpolate(f'heldout_{key[0]}_{key[1]}', 'heldout', nodes, key, exact=models[key])
    plan_path = ROOT/'config/experiments/rya1281_mac_fine_temperature.json'
    plan = json.loads(plan_path.read_text())
    for case in [c for c in plan['cases'] if c['name'].startswith('fine_T')]:
        nodes = [((t, g), models[t, g]) for t in case['bracket_teff'] for g in case['bracket_logg']]
        interpolate(case['name'], 'solar', nodes, (case['teff'], case['logg']), oxygen=case['abundance'])
    result = dict(schema='codex.bride_rosseland_column_audit/1', scientific_admission='HOLD', provider_admission=False,
        generator=fingerprint(__file__), boundary_code=fingerprint(ROOT/'pipeline/bride_rosseland_column.py'),
        native_format_history=dict(url='https://marcs.oreme.org/news/', update='April2018 RHOX recomputed by direct density integral; pressure/gravity proxy replaced. Alternative kernel reconstruction is not native producer parity.'),
        official_intake=fingerprint(source/'manifest.json'), official_source={k: v for k, v in entry.items() if k != 'relative_path'},
        extracted_fragment_sha256=hashlib.sha256(fragment.encode()).hexdigest(), build=dict(command=command,
            compiler_version=subprocess.check_output([compiler, '--version'], text=True).splitlines()[0],
            driver=fingerprint(build/'driver.f'), fragment=fingerprint(build/'official_kernel.f'), binary=fingerprint(binary), compiler_log=fingerprint(build/'compiler.log')),
        analytic=analytic, source_models=sources, plan=fingerprint(plan_path), binding_report=fingerprint(args.bindings/'report.json'),
        interpolation=weights, cases=rows,
        interpretation='Unchanged extracted calcrhox/rinteg/parcoe only, not full official interpolation or empirical calibration. Analytic REAL budget is an engineering guard, not stellar-structure accuracy. Mesh refinements retain one declared log-opacity interpolant, not new physical opacity points. Independent piecewise-log integral checks numerical quadrature, not model truth. Native and held-out RHOX are independent comparisons; actual solar has no native truth. Reanchoring is a declared constant offset, not a fitted scale. No source/certificate/production column changed or cap/provider enabled.')
    for name, data in [('summary.csv', rows), ('profiles.csv', profiles)]:
        with (args.output/name).open('x') as stream:
            writer = csv.DictWriter(stream, fieldnames=list(data[0]), lineterminator='\n')
            writer.writeheader(); writer.writerows(data)
        result[name] = fingerprint(args.output/name)
    (args.output/'report.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps(dict(analytic=len(analytic), cases=len(rows), official_kernel_executed=True, provider_admission=False)))


if __name__ == '__main__':
    main()
