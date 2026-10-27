"""Mac-only fixed-population oxygen audit; no SE, fitting or admission."""
from __future__ import annotations

import argparse
import csv
from dataclasses import replace
from importlib.metadata import version
import json
from pathlib import Path
import platform
import subprocess
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


def conserve_departure_replay(departures, nstar):
    """Apply external b to local LTE, explicitly renormalizing local O totals.

    This is a population-shape diagnostic, not identical-EOS MULTI replication.
    """
    if (departures.shape != nstar.shape or not np.isfinite(departures).all() or np.any(departures <= 0)
            or not np.isfinite(nstar).all() or np.any(nstar <= 0)):
        raise ValueError('Positive finite departures with matching layout required')
    trial = departures * nstar
    factors = nstar.sum(axis=0) / trial.sum(axis=0)
    return trial * factors, factors


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference', required=True, type=Path)
    parser.add_argument('--collision-runs', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    if platform.system() != 'Darwin' or version('lightweaver') != '0.16.1':
        raise RuntimeError('Authorized Mac / pinned LW0.16.1 required')
    from config.constants import require_codex_path
    from pipeline.bride_reference import fingerprint, read_marcs_structure, make_reference_context
    from pipeline.bride_parameter_checks import controlled_structure, marcs_filename, interpolate_marcs_structure
    from pipeline.multi_atom import read_multi_atom, to_lightweaver
    from pipeline.bride_observables import flux_quadrature, integrate_flux, window_integral
    from pipeline.wavelength_util import air_to_vac, vac_to_air
    from pipeline.kp_atlas_integrity import read_kurucz2005_residual
    from lightweaver.atomic_table import PeriodicTable
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    if fingerprint(args.reference)['sha256'] != 'c3113efdb21e3f40d79c1836b0ddb0ebb41dbceecaf04d622af6b7a8d8c080ec':
        raise ValueError('Changed independent reference')
    ref = json.loads(args.reference.read_text())
    atom_path = require_codex_path('experiments.bride_oxygen_atom')
    if fingerprint(atom_path) != ref['atom']:
        raise ValueError('Changed atom')
    source = replace(read_multi_atom(atom_path), abundance=8.70)
    store = require_codex_path('grids.marcs_standard')
    grid = require_codex_path('experiments.bride_mac_grid')
    observed = require_codex_path('data.kurucz2005_residual')
    if fingerprint(observed)['sha256'] != '75dd93479b5d3fe63c2e7828f5dbb3ce29f4c672f4f1147e0d30ce41bc8f5adb':
        raise ValueError('Changed corrected Kitt Peak holding')
    args.output.mkdir(parents=True, exist_ok=False)
    report = dict(schema='codex.bride_population_audit/1', scientific_admission='HOLD',
        status='RUNNING', host=dict(system=platform.system(), machine=platform.machine()),
        environment={k: version(k) for k in ('lightweaver', 'numpy', 'scipy', 'numba')},
        code_commit=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
        generator=fingerprint(__file__), context_builder=fingerprint(ROOT/'pipeline/bride_reference.py'),
        reference=fingerprint(args.reference), atom=fingerprint(atom_path), observation=fingerprint(observed),
        geometry='1D disk-integrated surface Fnu; five flux rays',
        controls=dict(threads=1, fixed_populations=True, SE_updates=False, J_rtol=1e-8,
                      max_iterations=1000, seconds_per_context=180),
        caveats=['No full blend/molecular spectrum; no fitted continuum/RV/width/abundance',
                 'Independent b replay renormalizes local LTE totals; not EOS-equivalent MULTI synthesis',
                 'Source b normalization and atom lineage remain open; no full3D claim'], cases=[])
    def save():
        (args.output/'report.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    save()
    model = to_lightweaver(source)
    report['oxygen_transitions'] = [dict(lower_0based=l.i, upper_0based=l.j,
        wavelength_air_A=float(vac_to_air(np.array([l.lambda0*10]))[0]),
        oscillator_f=l.f, loggf=float(np.log10(l.f*model.levels[l.i].g))) for l in model.lines]
    windows = [('triplet', 7768., 7780., [7771.944, 7774.166, 7775.388], .5),
               ('forbidden6300', 6299., 6302., [6300.304], .3),
               ('permitted6158', 6155., 6161., [6158.19], .5)]
    canonical_path = ROOT/'data/linelists/canonical_gf.csv'
    with canonical_path.open() as stream:
        canonical = list(csv.DictReader(stream))
    report['canonical_atomic_inventory_source'] = fingerprint(canonical_path)
    report['canonical_atomic_inventory_scope'] = 'Candidate rows by wavelength only; not strengths, blend completeness or admitted transitions. No molecular inventory.'
    report['canonical_atomic_inventory'] = {
        name: [row for row in canonical if any(abs(float(row['wavelength_air_A'])-c) <= half for c in centres)]
        for name, _, _, centres, half in windows}
    for line in report['oxygen_transitions']:
        candidates = [r for r in canonical if r['species'] == 'O I' and abs(float(r['wavelength_air_A'])-line['wavelength_air_A']) < .1]
        line['nearby_canonical_O_rows'] = [dict(line_id=r['line_id'], wavelength_air_A=float(r['wavelength_air_A']),
            loggf=float(r['log_gf']), atom_minus_canonical_loggf=line['loggf']-float(r['log_gf'])) for r in candidates]
    from lightweaver import rh_atoms
    report['native_passive_atom_source'] = fingerprint(rh_atoms.__file__)
    report['blend_handling'] = 'Only O is detailed static. Passive species supply native LW background; neighbouring atomic bound-bound blends and molecules are absent. Never multiply independent normalized spectra to simulate mixed opacity.'
    for name, relative in [('native5750', 'rya1281-v1/baseline'),
                           ('solar5772', 'rya1281-mesh-v1/meshfine_T5772')]:
        origin = grid/relative
        manifest = json.loads((origin/'manifest.json').read_text())
        if manifest['status'] != 'CONVERGED_EXPERIMENT':
            raise ValueError('Unconverged input')
        committed = ROOT/'data/results/rya1281/mac_parameter_checks_v1'/('coarse/baseline' if name=='native5750' else 'mesh10K/meshfine_T5772')/'manifest.json'
        if fingerprint(committed) != fingerprint(origin/'manifest.json'):
            raise ValueError('Source manifest differs from reviewed archive')
        if fingerprint(origin/'flux.csv') != manifest['flux']:
            raise ValueError('Changed source flux')
        case = manifest['case']
        if name == 'native5750':
            path = store/marcs_filename(5750, 4.5, 0., 0.)
            if fingerprint(path) != ref['atmosphere_file']:
                raise ValueError('Changed native atmosphere')
            structure = read_marcs_structure(path)
        else:
            nodes = []
            for teff in case['bracket_teff']:
                for logg in case['bracket_logg']:
                    path = store/marcs_filename(teff, logg, 0., 0.)
                    if fingerprint(path) not in manifest['atmosphere']['sources']:
                        raise ValueError('Changed interpolation source')
                    nodes.append(((teff, logg), read_marcs_structure(path)))
            structure, _ = interpolate_marcs_structure(nodes, (case['teff'], case['logg']),
                thermodynamic_scale=case['thermodynamic_scale'])
        structure = controlled_structure(structure, xi=case['xi'], drop_bottom=case['drop_bottom'])
        saved = np.load(origin/'populations.npz')
        np.testing.assert_array_equal(saved['temperature_K'], structure['temperature_K'])
        populations, nstar = saved['populations'], saved['nstar']
        modes = [('local_nlte', populations), ('local_lte', nstar)]
        result = dict(name=name, source_manifest=fingerprint(origin/'manifest.json'),
                      source_populations=fingerprint(origin/'populations.npz'), modes=[])
        if name == 'native5750':
            np.testing.assert_allclose(np.log10(ref['reference']['tau']), structure['logtau_500'], atol=2e-8, rtol=0)
            external, factors = conserve_departure_replay(np.array(ref['reference']['departures']), nstar)
            modes.append(('external_b_local_lte', external))
            result['external_population_renormalization_range'] = [float(factors.min()), float(factors.max())]
        fig, axes = plt.subplots(1, 3, figsize=(14, 4))
        for mode, pop in modes:
            mode_start = time.perf_counter()
            flux_pair = []
            mode_result = dict(name=mode, radiation_iterations=[], windows=[])
            for continuum_only in (False, True):
                start = time.perf_counter()
                ctx = make_reference_context(structure, source, rays=3, qscale=1,
                    fixed_populations=pop, continuum_only=continuum_only)
                actual = ctx.kwargs['eqPops'].atomicPops[PeriodicTable['O']]
                np.testing.assert_allclose(actual.nStar, nstar, rtol=64*np.finfo(float).eps, atol=0)
                for iteration in range(1000):
                    if time.perf_counter()-start > 180:
                        raise TimeoutError('Fixed-population context exceeded 180 seconds')
                    u = ctx.formal_sol_gamma_matrices()
                    if u.dJMax < 1e-8:
                        break
                else:
                    raise RuntimeError('Radiation failed convergence')
                np.testing.assert_array_equal(actual.n, pop)
                mode_result['radiation_iterations'].append(iteration+1)
                mu, weights = flux_quadrature(5)
                curves = []
                for _, lo, hi, _, _ in windows:
                    wave = np.arange(lo, hi+.0025, .005)
                    curves.append(integrate_flux(ctx.compute_rays(wavelengths=air_to_vac(wave)*.1, mus=mu), mu, weights))
                flux_pair.append(curves)
            for wi, (window, lo, hi, centres, halfwidth) in enumerate(windows):
                wave = np.arange(lo, hi+.0025, .005)
                norm = flux_pair[0][wi]/flux_pair[1][wi]
                if not np.isfinite(norm).all() or np.any(flux_pair[1][wi] <= 0):
                    raise ValueError('Invalid emergent flux')
                wo, fo = read_kurucz2005_residual(observed, lo, hi)
                if wo[0] > lo or wo[-1] < hi:
                    raise ValueError('Observed window not bracketed')
                rows = [dict(centre_air_A=c, model_area_mA=1000*window_integral(wave, norm, c-halfwidth, c+halfwidth),
                    observed_area_mA=1000*window_integral(wo, fo, c-halfwidth, c+halfwidth)) for c in centres]
                mode_result['windows'].append(dict(name=window, halfwidth_A=halfwidth, lines=rows))
                np.savetxt(args.output/f'{name}_{mode}_{window}.csv',
                    np.column_stack((wave, flux_pair[0][wi], flux_pair[1][wi], norm, np.interp(wave, wo, fo))),
                    delimiter=',', header='air_A,Fnu_W_m2_Hz,continuum_Fnu_W_m2_Hz,normalized_flux,kitt_peak', comments='')
                axes[wi].plot(wave, norm, label=mode)
                if mode == 'local_nlte':
                    axes[wi].plot(wo, fo, color='black', lw=.8, label='Kitt Peak corrected')
                if mode == 'local_nlte' and wi == 0:
                    old = np.genfromtxt(origin/'flux.csv', delimiter=',', names=True)
                    np.testing.assert_allclose(wave, old['wavelength_air_A'], atol=1e-8, rtol=0)
                    mode_result['triplet_replay_max_absolute_difference'] = float(np.max(abs(norm-old['normalized_flux'])))
                axes[wi].set(title=window, xlabel='Air wavelength (Angstrom)', ylabel='Normalized flux')
            mode_result['elapsed_seconds'] = time.perf_counter()-mode_start
            result['modes'].append(mode_result)
            print(name, mode, mode_result['windows'], flush=True)
        axes[0].legend(fontsize=7)
        fig.suptitle(f'{name}, fixed O8.70; oxygen-only diagnostic, no blend fit')
        fig.tight_layout(); fig.savefig(args.output/f'{name}_profiles.png', dpi=150); plt.close(fig)
        triplet = [l for l in report['oxygen_transitions'] if 7770 < l['wavelength_air_A'] < 7777]
        depth = [structure['logtau_500'], structure['temperature_K']]
        headers = ['logtau500', 'temperature_K']
        fig, ax = plt.subplots(figsize=(7, 4))
        for mode, pop in modes:
            b = pop/nstar
            for li, line in enumerate(triplet):
                lower, upper = line['lower_0based'], line['upper_0based']
                depth.extend([b[lower], b[upper], b[upper]/b[lower]])
                headers.extend([f'{mode}_{li}_b_lower', f'{mode}_{li}_b_upper', f'{mode}_{li}_bu_bl'])
                ax.plot(structure['logtau_500'], b[upper]/b[lower], label=f'{mode} {line["wavelength_air_A"]:.2f}')
        np.savetxt(args.output/f'{name}_departures.csv', np.column_stack(depth), delimiter=',', header=','.join(headers), comments='')
        ax.set(xlabel='Native log tau500', ylabel='b_upper / b_lower', title=f'{name}: population ratios; not a contribution function')
        ax.legend(fontsize=6); fig.tight_layout(); fig.savefig(args.output/f'{name}_departures.png', dpi=150); plt.close(fig)
        report['cases'].append(result); save()
    report['status'] = 'DIAGNOSTIC_COMPLETE'
    collision_summary = json.loads((args.collision_runs/'summary.json').read_text())
    if collision_summary['status'] != 'EXPERIMENTS_FINISHED' or collision_summary['atom'] != ref['atom']:
        raise ValueError('Collision control set incomplete or atom changed')
    if collision_summary['context_builder'] != report['context_builder']:
        raise ValueError('Collision controls use different context physics')
    report['collision_summary'] = fingerprint(args.collision_runs/'summary.json')
    report['collision_cases'] = []
    fig, ax = plt.subplots(figsize=(10, 4))
    for mode in ('local_lte', 'local_nlte', 'external_b_local_lte'):
        curve = np.genfromtxt(args.output/f'native5750_{mode}_triplet.csv', delimiter=',', names=True)
        ax.plot(curve['air_A'], curve['normalized_flux'], label=mode)
    for entry in collision_summary['cases']:
        directory = args.collision_runs/entry['name']
        if fingerprint(directory/'manifest.json') != entry['manifest']:
            raise ValueError('Changed collision manifest')
        manifest = json.loads((directory/'manifest.json').read_text())
        if manifest['status'] != 'CONVERGED_EXPERIMENT' or manifest['case']['teff'] != 5750 or manifest['case']['logg'] != 4.5:
            raise ValueError('Collision experiment failed or wrong reference node')
        if fingerprint(directory/'flux.csv') != manifest['flux']:
            raise ValueError('Changed collision spectrum')
        curve = np.genfromtxt(directory/'flux.csv', delimiter=',', names=True)
        rows = [dict(centre_air_A=c, model_area_mA=1000*window_integral(curve['wavelength_air_A'], curve['normalized_flux'], c-.5, c+.5))
                for c in (7771.944, 7774.166, 7775.388)]
        report['collision_cases'].append(dict(name=entry['name'], manifest=entry['manifest'],
            source_directory=directory.name, source_run=args.collision_runs.name,
            flux=manifest['flux'], populations=fingerprint(directory/'populations.npz'),
            diagnostics=manifest['final_diagnostics'], scales=manifest['diagnostic_collision_scales'], lines=rows))
        np.savetxt(args.output/f'{entry["name"]}_triplet.csv',
            np.column_stack((curve['wavelength_air_A'], curve['normalized_flux'])), delimiter=',',
            header='air_A,normalized_flux', comments='')
        ax.plot(curve['wavelength_air_A'], curve['normalized_flux'], label=entry['name'], alpha=.8)
    ax.set(xlabel='Air wavelength (Angstrom)', ylabel='Intrinsic normalized flux',
           title='Native 5750/4.5/O8.70: causal collision controls, no rate adoption')
    ax.legend(fontsize=7); fig.tight_layout(); fig.savefig(args.output/'collision_controls.png', dpi=150); plt.close(fig)
    save()


if __name__ == '__main__':
    main()
