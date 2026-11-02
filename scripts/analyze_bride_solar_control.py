"""Review one bound solar O control alongside the native Fe-sampled diagnostic."""
import argparse
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--solar-run', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    from config.constants import require_codex_path
    from pipeline.bride_reference import fingerprint, read_marcs_structure
    from pipeline.bride_opacity_binding import binding_digest
    from pipeline.bride_observables import window_integral
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    summary = json.loads((args.solar_run/'summary.json').read_text())
    if summary['status'] != 'EXPERIMENTS_FINISHED' or len(summary['cases']) != 1:
        raise ValueError('Requires the completed one-case solar experiment')
    old_receipt = ROOT/'data/results/rya1281/fe_background_v1/analysis/summary.json'
    native = next(x for x in json.loads(old_receipt.read_text())['cases'] if x['name'] == 'author_MN_cap_Fe_neutral_sampled')
    store = require_codex_path('experiments.bride_mac_grid')
    items = [(store/native['source_run']/native['name'], native),
             (args.solar_run/summary['cases'][0]['name'], None)]
    args.output.mkdir(parents=True, exist_ok=False)
    result = dict(schema='codex.bride_solar_control_review/1', scientific_admission='HOLD',
        generator=fingerprint(__file__), source_solar_summary=fingerprint(args.solar_run/'summary.json'),
        native_receipt=fingerprint(old_receipt), cases=[],
        scope='Native5750/4.5 vs experimental solar5772/4.438, fixed O8.70. T/gravity/structure/cap and exact-input coefficients change together; not a pure temperature derivative or external-reference population test. Each b uses its own local nstar.')
    fig, axes = plt.subplots(1, 3, figsize=(13, 4))
    spectrum, ax = plt.subplots(figsize=(11, 4))
    for directory, prior in items:
        manifest_path = directory/'manifest.json'
        m = json.loads(manifest_path.read_text())
        if m['status'] != 'CONVERGED_EXPERIMENT':
            raise ValueError('Unconverged control')
        if prior is not None:
            for filename, key in (('manifest.json', 'source_manifest'), ('populations.npz', 'populations'), ('flux.csv', 'flux')):
                if fingerprint(directory/filename) != prior[key]:
                    raise ValueError('Changed native control')
            s = read_marcs_structure(require_codex_path('grids.marcs_standard')/m['atmosphere']['name'])
            tau = s['logtau_500']
        else:
            if fingerprint(manifest_path) != summary['cases'][0]['manifest']:
                raise ValueError('Changed solar manifest')
            record_path = directory/'experimental_input.json'
            if fingerprint(record_path) != m['experimental_bound_input']['receipt']:
                raise ValueError('Changed fresh experimental binding')
            bound = json.loads(record_path.read_text())
            if binding_digest(bound['local_binding']) != bound['binding_sha256']:
                raise ValueError('Changed solar binding fields')
            s = {k:np.asarray(v) for k, v in bound['local_binding']['fields'].items()}
            tau = s['logtau_500']
        if m['case']['abundance'] != 8.7 or m['experimental_collision_convention'] != 'barklem2018_MN_LTE_pair':
            raise ValueError('Changed fixed O/collision contract')
        with np.load(directory/'populations.npz') as data:
            n, nstar = data['populations'].copy(), data['nstar'].copy()
            np.testing.assert_array_equal(data['temperature_K'], s['temperature_K'])
        if n.shape != (41, 56) or nstar.shape != n.shape or not np.isfinite(n).all() or not np.isfinite(nstar).all() or np.any(n <= 0) or np.any(nstar <= 0):
            raise ValueError('Invalid full-atom population arrays')
        np.testing.assert_allclose(n.sum(axis=0), nstar.sum(axis=0), rtol=1e-10, atol=0)
        history = [json.loads(x) for x in (directory/'history.jsonl').read_text().splitlines()]
        if history[-1] != m['final_diagnostics'] or len(history) != m['iterations']:
            raise ValueError('Changed convergence history')
        if max(history[-1][k] for k in ('population_change','J_change','physical_rate_residual','cap_boundary_change')) >= m['controls']['rtol']:
            raise ValueError('Incomplete joint convergence')
        flux_path = directory/'flux.csv'
        if fingerprint(flux_path) != m['flux']:
            raise ValueError('Changed emergent flux')
        flux = np.genfromtxt(flux_path, delimiter=',', names=True)
        b = n/nstar
        curve_path = args.output/(directory.name+'_departures.csv')
        np.savetxt(curve_path, np.column_stack((tau,b[5],b[9],b[39])), delimiter=',',
            header='logtau500,lower_triplet_b,upper_triplet_b,OII_ground_b', comments='')
        for panel, level, label in zip(axes, (5,9,39), ('Lower triplet', 'Upper triplet', 'O II ground')):
            panel.plot(tau, b[level], label=directory.name)
            panel.set(xlabel='log tau500', ylabel='Departure b', title=label)
        axes[-1].set_yscale('log')
        ax.plot(flux['wavelength_air_A'], flux['normalized_flux'], label=directory.name)
        metrics = [dict(centre_air_A=c, intrinsic_window_area_mA=1000*window_integral(flux['wavelength_air_A'],flux['normalized_flux'],c-.5,c+.5)) for c in (7771.944,7774.166,7775.388)]
        result['cases'].append(dict(name=directory.name, source_run=directory.parent.name,
            source_manifest=fingerprint(manifest_path), flux=fingerprint(flux_path), populations=fingerprint(directory/'populations.npz'),
            history=fingerprint(directory/'history.jsonl'), departures=fingerprint(curve_path),
            parameters={k:m['case'][k] for k in ('teff','logg','xi','abundance','drop_bottom')},
            iterations=m['iterations'], elapsed_seconds=m['elapsed_seconds'], final_diagnostics=m['final_diagnostics'],
            top_departures=dict(lower_triplet=float(b[5,0]),upper_triplet=float(b[9,0]),OII_ground=float(b[39,0])),
            effective_cap_length_m=m['finite_column_top_cap']['effective_length_m'], metrics=metrics,
            experimental_bound_input=m.get('experimental_bound_input')))
    for panel in axes:
        panel.legend(fontsize=6); panel.grid(alpha=.2)
    fig.suptitle('Fixed O8.70: native reference versus experimental solar inputs; not an error test')
    fig.tight_layout(); figure=args.output/'departures.png'; fig.savefig(figure,dpi=150); plt.close(fig)
    result['departure_figure']=fingerprint(figure)
    ax.set(xlim=(7770.5,7776.8),xlabel='AIR wavelength (Angstrom)',ylabel='Normalized surface flux')
    ax.legend(fontsize=7); ax.grid(alpha=.2); spectrum.tight_layout()
    figure=args.output/'intrinsic_profiles.png'; spectrum.savefig(figure,dpi=150); plt.close(spectrum)
    result['intrinsic_figure']=fingerprint(figure)
    (args.output/'report.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps([dict(name=c['name'],iterations=c['iterations'],top_departures=c['top_departures'],metrics=c['metrics']) for c in result['cases']],indent=2))


if __name__ == '__main__':
    main()
