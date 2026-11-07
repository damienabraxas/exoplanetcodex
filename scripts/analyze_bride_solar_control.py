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
    parser.add_argument('--solar-run', type=Path, nargs='+', required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    from config.constants import require_codex_path
    from pipeline.bride_reference import fingerprint, read_marcs_structure
    from pipeline.bride_opacity_binding import binding_digest
    from pipeline.bride_observables import window_integral
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    summaries = [(path, json.loads((path/'summary.json').read_text())) for path in args.solar_run]
    if any(s['status'] != 'EXPERIMENTS_FINISHED' for path, s in summaries):
        raise ValueError('Requires completed explicitly supplied solar experiments')
    old_receipt = ROOT/'data/results/rya1281/fe_background_v1/analysis/summary.json'
    native = next(x for x in json.loads(old_receipt.read_text())['cases'] if x['name'] == 'author_MN_cap_Fe_neutral_sampled')
    store = require_codex_path('experiments.bride_mac_grid')
    native_root=store/native['source_run']
    if fingerprint(native_root/'summary.json')!=native['source_summary']:
        raise ValueError('Changed native source summary')
    native_summary=json.loads((native_root/'summary.json').read_text())
    if any(s['atom']!=native_summary['atom'] for path,s in summaries):
        raise ValueError('Mixed atomic models in temperature comparison')
    native_controls=json.loads((native_root/native['name']/'manifest.json').read_text())['controls']
    items = [(store/native['source_run']/native['name'], native, None)]
    for path, summary in summaries:
        items.extend((path/entry['name'], None, entry['manifest']) for entry in summary['cases'])
    if len({p.name for p, prior, expected in items}) != len(items):
        raise ValueError('Duplicate case identities in supplied runs')
    args.output.mkdir(parents=True, exist_ok=False)
    result = dict(schema='codex.bride_solar_control_review/1', scientific_admission='HOLD',
        generator=fingerprint(__file__), source_solar_summaries=[dict(source_run=p.name, receipt=fingerprint(p/'summary.json')) for p, s in summaries],
        native_receipt=fingerprint(old_receipt), cases=[],
        scope='Native5750/4.5 vs experimental solar-temperature inputs at4.438, fixed O8.70. Native-to-solar changes T/gravity/structure/cap/coefficients together. Within solar series only declared Teff and its constructed fields/coefficients vary. Not external-reference population errors or an atmosphere-equilibrium certification. Each b uses its own local nstar.')
    fig, axes = plt.subplots(1, 3, figsize=(13, 4))
    spectrum, ax = plt.subplots(figsize=(11, 4))
    for directory, prior, expected_manifest in items:
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
            if fingerprint(manifest_path) != expected_manifest:
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
        if m['controls']!=native_controls or m.get('diagnostic_collision_scales') or m['case'].get('lte_h_background_lines') or m['case'].get('lte_uv_background_lines'):
            raise ValueError('Mixed numerics/collision scales/background-line physics')
        if prior is None and (m['case']['logg'],m['case']['xi'],m['case']['drop_bottom']) != (4.438,1.,0):
            raise ValueError('Changed solar temperature-control parameters')
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
        display = ('Native ' if prior is not None else 'Experimental ')+f"T{m['case']['teff']} g{m['case']['logg']}"
        curve_path = args.output/(directory.name+'_departures.csv')
        np.savetxt(curve_path, np.column_stack((tau,b[5],b[9],b[39])), delimiter=',',
            header='logtau500,lower_triplet_b,upper_triplet_b,OII_ground_b', comments='')
        for panel, level, label in zip(axes, (5,9,39), ('Lower triplet', 'Upper triplet', 'O II ground')):
            panel.plot(tau, b[level], label=display)
            panel.set(xlabel='log tau500', ylabel='Departure b', title=label)
        axes[-1].set_yscale('log')
        ax.plot(flux['wavelength_air_A'], flux['normalized_flux'], label=display)
        metrics = [dict(centre_air_A=c, intrinsic_window_area_mA=1000*window_integral(flux['wavelength_air_A'],flux['normalized_flux'],c-.5,c+.5)) for c in (7771.944,7774.166,7775.388)]
        result['cases'].append(dict(name=directory.name, display_label=display, source_run=directory.parent.name,
            source_manifest=fingerprint(manifest_path), flux=fingerprint(flux_path), populations=fingerprint(directory/'populations.npz'),
            history=fingerprint(directory/'history.jsonl'), departures=fingerprint(curve_path),
            parameters={k:m['case'][k] for k in ('teff','logg','xi','abundance','drop_bottom')},
            iterations=m['iterations'], elapsed_seconds=m['elapsed_seconds'], final_diagnostics=m['final_diagnostics'],
            top_departures=dict(lower_triplet=float(b[5,0]),upper_triplet=float(b[9,0]),OII_ground=float(b[39,0])),
            effective_cap_length_m=m['finite_column_top_cap']['effective_length_m'], metrics=metrics,
            experimental_bound_input=m.get('experimental_bound_input')))
    solar_cases=sorted([c for c in result['cases'] if c['experimental_bound_input'] is not None],key=lambda c:c['parameters']['teff'])
    result['temperature_steps']=[]
    for low,high in zip(solar_cases[:-1],solar_cases[1:]):
        step=high['parameters']['teff']-low['parameters']['teff']
        if step!=10:raise ValueError('Supplied solar series is not an adjacent10K sequence')
        result['temperature_steps'].append(dict(low_teff=low['parameters']['teff'],high_teff=high['parameters']['teff'],
            window_area_changes_mA=[b['intrinsic_window_area_mA']-a['intrinsic_window_area_mA'] for a,b in zip(low['metrics'],high['metrics'])]))
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
