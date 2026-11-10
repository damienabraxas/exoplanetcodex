"""Verify matched numerical-only endpoint replays; no observation or parameter fit."""
import argparse
import json
from pathlib import Path
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
ENDPOINTS = ('thermal_T5572', 'thermal_T5822', 'thermal_g4338', 'thermal_g4538')
CENTRES = (7771.944, 7774.166, 7775.388)  # AIR A, same held O777 apertures
RESOURCE = 'experiments.bride_thermal_extreme_refinements'


def verify_pair(base, refined, plan, *, ng_control=None, extended_budget=False):
    from pipeline.bride_thermal_sweep import effective_controls
    for r in (base, refined):
        if (r['schema'] != 'codex.bride_thermal_flux/1' or r['mode'] != 'cowley'
                or r['population_mode'] != 'coupled' or r['scientific_admission'] != 'HOLD'):
            raise ValueError('Wrong numerical replay physics')
    if base['numeric_refinement'] or not refined['numeric_refinement']:
        raise ValueError('Require baseline/refinement pair')
    for key in ('case', 'plan', 'sweep_adapter', 'certificate',
                'atmosphere_sources', 'local_binding', 'stark_source', 'stark_kernel',
                'pressure_source', 'pressure_provider', 'provider', 'line_coefficients',
                'lines', 'partition', 'source_manifest', 'source_populations',
                'source_flux', 'source_binding', 'wide_baseline'):
        if base[key] != refined[key]:
            raise ValueError('Refinement changes physical/source identity: ' + key)
    initial = refined.get('initialization', {})
    if initial.get('mode') == 'matched_physical_baseline':
        from pipeline.bride_reference import fingerprint
        if (initial['populations'] != base['coupled']['populations']
                or initial['helper'] != fingerprint(ROOT/'pipeline/bride_thermal_initial.py')
                or refined['generator'] != fingerprint(ROOT/'scripts/run_bride_thermal_sweep.py')):
            raise ValueError('Changed matched initialization/runner provenance')
    elif base['generator'] != refined['generator']:
        raise ValueError('Changed legacy generator')
    expected=effective_controls(plan['controls'], True)
    if ng_control is not None:
        if ng_control not in ('record','off') or refined.get('ng_control')!=ng_control:raise ValueError('Wrong declared Ng comparison')
        expected['ng']=ng_control=='record'
    elif refined.get('ng_control') is not None:raise ValueError('Declare the Ng comparison explicitly')
    if extended_budget:
        from pipeline.bride_compute_budget import extended_controls
        from pipeline.bride_reference import fingerprint
        declaration=ROOT/'config/experiments/rya1281_ng_extended.json'
        expected=extended_controls(expected,json.loads(declaration.read_text()),refined['case']['name'],ng_control)
        if (refined.get('compute_extension')!=fingerprint(declaration)
                or refined.get('budget_adapter')!=fingerprint(ROOT/'pipeline/bride_compute_budget.py')):raise ValueError('Changed declared budget provenance')
    elif refined.get('compute_extension') is not None:raise ValueError('Declare compute extension explicitly')
    if (base['coupled']['controls'] != plan['controls'] or
            refined['coupled']['controls'] != expected):
        raise ValueError('Changed numerical refinement controls/budgets')
    for key in ('source_spline', 'source_mn', 'helper'):
        if base['coupled'][key] != refined['coupled'][key]:
            raise ValueError('Changed coupled source: ' + key)


def independent_balance(pop, radiative, collision, stored_total, target_total, rtol):
    """Rebuild the rate operator from offdiagonal arrays and check target nuclei."""
    if (pop.ndim != 2 or radiative.shape != (pop.shape[0], pop.shape[0], pop.shape[1])
            or collision.shape != radiative.shape or stored_total.shape != radiative.shape
            or target_total.shape != (pop.shape[1],)):
        raise ValueError('Wrong rate/population shape')
    if any(not np.isfinite(v).all() for v in (pop, radiative, collision, stored_total, target_total)):
        raise ValueError('Nonfinite rate/population evidence')
    if np.any(pop <= 0) or np.any(target_total <= 0) or np.any(radiative < 0) or np.any(collision < 0):
        raise ValueError('Invalid populations or offdiagonal rates')
    total = radiative + collision
    indices = np.arange(len(total))
    if np.any(total[indices, indices] != 0):
        raise ValueError('Offdiagonal arrays have nonzero diagonal')
    outgoing = total.sum(axis=0)
    if np.any(outgoing <= 0):
        raise ValueError('Disconnected rate equation')
    total[indices, indices] = -outgoing
    np.testing.assert_allclose(total, stored_total, rtol=128*np.finfo(float).eps, atol=0)
    net = np.einsum('ijd,jd->id', total, pop)
    gross = np.einsum('ijd,jd->id', abs(total), pop)
    residual = float(np.max(abs(net)/np.maximum(gross, np.finfo(float).tiny)))
    conservation = float(np.max(abs(pop.sum(axis=0)/target_total-1)))
    if residual >= rtol or conservation > 1e-10:
        raise ValueError('Independent rate balance/particle conservation failed')
    return dict(residual=residual, particle_conservation_error=conservation)


def failed_attempt(directory, base, plan, log, matched):
    """Explicitly retain a verified budget refusal, with no inferred spectrum."""
    from pipeline.bride_reference import fingerprint
    from pipeline.bride_thermal_sweep import effective_controls
    inputs = json.loads((directory/'input_manifest.json').read_text())
    history = [json.loads(line) for line in (directory/'history.jsonl').read_text().splitlines()]
    if ((directory/'report.json').exists() or not history or not inputs['numeric_refinement']
            or inputs['schema'] != 'codex.bride_thermal_input/1'
            or inputs['mode'] != 'cowley' or inputs['population_mode'] != 'coupled'
            or inputs['controls'] != effective_controls(plan['controls'], True)
            or log.read_text().rstrip().splitlines()[-1] not in (
                'TimeoutError: Fixed-population cap radiation budget exceeded',
                'TimeoutError: Original coupled budget exceeded')):
        raise ValueError('Wrong timeout evidence')
    for key in ('case', 'plan', 'sweep_adapter', 'certificate', 'atmosphere_sources',
                'local_binding', 'source_manifest', 'source_binding', 'lines', 'pressure_source', 'stark_source'):
        if inputs[key] != base[key]: raise ValueError('Failed attempt changes physical/source input: '+key)
    if matched:
        initial = inputs['initialization']
        if (initial['mode'] != 'matched_physical_baseline' or initial['populations'] != base['coupled']['populations']
                or initial['helper'] != fingerprint(ROOT/'pipeline/bride_thermal_initial.py')
                or inputs['generator'] != fingerprint(ROOT/'scripts/run_bride_thermal_sweep.py')):
            raise ValueError('Wrong failed matched initialization')
    elif inputs['generator'] != base['generator']:
        raise ValueError('Wrong failed legacy generator')
    return dict(status='BUDGET_EXHAUSTED_NOT_ACCEPTED', case=base['case']['name'],
        input_manifest=fingerprint(directory/'input_manifest.json'), history=fingerprint(directory/'history.jsonl'),
        log=fingerprint(log), completed_iterations=len(history), last_diagnostics=history[-1])


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--plan', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--failed-attempt-log', type=Path, required=True)
    p.add_argument('--failed-refinements', nargs='*', default=[])
    p.add_argument('--failure-log-dir', type=Path)
    args = p.parse_args()
    from config.constants import require_codex_path
    from pipeline.bride_reference import fingerprint
    from pipeline.bride_thermal_sweep import validate_plan
    from scripts.analyze_bride_thermal_sweep import verify_case_binding
    from pipeline.bride_observables import window_integral
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    plan = json.loads(args.plan.read_text())
    validate_plan(plan, json.loads((ROOT/'config/experiments/rya1281_mac_solar_bound_control.json').read_text()))
    store = require_codex_path('experiments.bride_mac_grid')
    refined_root = require_codex_path(RESOURCE)
    coefficients = require_codex_path('experiments.bride_thermal_coefficients')
    result = dict(schema='codex.bride_thermal_refinement_review/1', scientific_admission='HOLD',
                  generator=fingerprint(__file__), plan=fingerprint(args.plan), cases=[],
                  scope='Combined angular/wavelength/flux-ray sampling at four predeclared endpoints plus existing nominal. No fit, new physical accuracy threshold, or individual error attribution. Observations and continuum controls unchanged.')
    failed = require_codex_path('experiments.bride_thermal_extreme_legacy')/'thermal_T5572'
    failed_base = json.loads((store/'rya1281-thermal-coupled-v1/thermal_T5572/report.json').read_text())
    result['failed_original_solar_start'] = dict(failed_attempt(failed, failed_base, plan, args.failed_attempt_log, False),
        resource_key='experiments.bride_thermal_extreme_legacy')
    if len(set(args.failed_refinements)) != len(args.failed_refinements) or not set(args.failed_refinements) <= set(ENDPOINTS):
        raise ValueError('Unknown/duplicate failed refinement')
    fig, axes = plt.subplots(2, 1, figsize=(12, 8))
    profiles = []; labels = ['wavelength_air_A']; wave = None
    for name in ENDPOINTS + ('thermal_nominal',):
        baseline = store/'rya1281-thermal-coupled-v1'/name
        refined = refined_root/name if name in ENDPOINTS else store/'rya1281-thermal-refined-v1'
        if name in args.failed_refinements:
            if args.failure_log_dir is None: raise ValueError('Explicit failed-refinement logs required')
            base = json.loads((baseline/'report.json').read_text())
            if fingerprint(baseline/'report.json') != fingerprint(ROOT/'data/results/rya1281/thermal_slice_v1/runs'/name/'report.json'):
                raise ValueError('Changed failed-case baseline')
            entry = failed_attempt(refined, base, plan, args.failure_log_dir/('rya1281-matched-'+name+'.log'), True)
            entry['resource_key'] = RESOURCE
            result['cases'].append(entry)
            continue
        receipts = [json.loads((d/'report.json').read_text()) for d in (baseline, refined)]
        base, ref = receipts
        if fingerprint(baseline/'report.json') != fingerprint(ROOT/'data/results/rya1281/thermal_slice_v1/runs'/name/'report.json'):
            raise ValueError('Baseline differs from archived physical control')
        verify_pair(base, ref, plan)
        if name in ENDPOINTS and (ref.get('initialization', {}).get('mode') != 'matched_physical_baseline'
                or ref['initialization']['report'] != fingerprint(baseline/'report.json')):
            raise ValueError('Wrong matched baseline initialization')
        if name in ENDPOINTS:
            with np.load(baseline/'populations.npz') as old, np.load(refined/'populations.npz') as new:
                np.testing.assert_array_equal(old['populations'], new['initial_populations'])
        declared = next(c for c in plan['cases'] if c['name'] == name)
        if base['case'] != declared or base['plan'] != fingerprint(args.plan):
            raise ValueError('Relabelled case/plan')
        cert = coefficients/name/'manifest.json'
        for r in receipts:
            if r['certificate'] != fingerprint(cert): raise ValueError('Changed exact certificate')
            verify_case_binding(r, json.loads(cert.read_text()))
        balance = []; arrays = []
        for d, r in zip((baseline, refined), receipts):
            for filename, key in (('input_manifest.json', 'input_manifest'), ('flux.csv', 'flux'),
                                  ('neutral_colliders.csv', 'neutral_colliders'), ('neutral_reservoir.csv', 'neutral_reservoir'), ('widths.csv', 'widths')):
                if fingerprint(d/filename) != r[key]: raise ValueError('Changed replay evidence: '+filename)
            for filename, key in (('populations.npz', 'populations'), ('rates.npz', 'rates'), ('history.jsonl', 'history')):
                if fingerprint(d/filename) != r['coupled'][key]: raise ValueError('Changed coupled evidence: '+filename)
            inputs = json.loads((d/'input_manifest.json').read_text())
            if (inputs['case'] != r['case'] or inputs['controls'] != r['coupled']['controls']
                    or inputs['local_binding'] != r['local_binding']
                    or inputs.get('initialization') != r.get('initialization')
                    or inputs['numeric_refinement'] != r['numeric_refinement']):
                raise ValueError('Input/report identity disagreement')
            history = [json.loads(line) for line in (d/'history.jsonl').read_text().splitlines()]
            if (len(history) != r['coupled']['iterations'] or len(history) < 3
                    or history[-1] != r['coupled']['final_diagnostics']):
                raise ValueError('Wrong joint convergence history')
            last = history[-1]
            if (max(last[k] for k in ('population_change', 'physical_rate_residual', 'J_change', 'cap_boundary_change')) >= plan['controls']['rtol']
                    or last['particle_conservation_error'] > 1e-10):
                raise ValueError('Joint convergence failed')
            with np.load(d/'populations.npz') as saved: pop = saved['populations'].copy()
            binding = r['local_binding']['local_binding']
            target = np.array(binding['nH_m3'])*10**(binding['scalars']['active_oxygen_abundance_dex']-12)
            with np.load(d/'rates.npz') as rates:
                balance.append(independent_balance(pop, rates['radiative_rates_offdiag_s_inv'],
                    rates['collision_rates_offdiag_s_inv'], rates['total_rates_s_inv'], target, plan['controls']['rtol']))
            arrays.append(np.genfromtxt(d/'flux.csv', delimiter=',', names=True))
        a, b = arrays
        if any(not all(np.isfinite(x[k]).all() for k in x.dtype.names) or np.any(x['continuum_Fnu_W_m2_Hz'] <= 0) for x in arrays):
            raise ValueError('Invalid flux/continuum arrays')
        for x in arrays:
            np.testing.assert_array_equal(x['normalized_flux'], x['Fnu_W_m2_Hz']/x['continuum_Fnu_W_m2_Hz'])
        np.testing.assert_array_equal(a['wavelength_air_A'], b['wavelength_air_A'])
        if wave is None: wave = a['wavelength_air_A']; profiles.append(wave)
        np.testing.assert_array_equal(wave, a['wavelength_air_A'])
        delta = b['normalized_flux'] - a['normalized_flux']
        area = [[1000*window_integral(wave, x['normalized_flux'], c-.5, c+.5) for c in CENTRES] for x in arrays]
        profile_max = [float(np.max(abs(delta[(wave >= c-.5) & (wave <= c+.5)]))) for c in CENTRES]
        entry = dict(case=name, status='CONVERGED_EXPERIMENT', baseline_report=fingerprint(baseline/'report.json'), refined_report=fingerprint(refined/'report.json'),
            independent_balance=balance, baseline_area_mA=area[0], refined_area_mA=area[1],
            area_delta_mA=(np.array(area[1])-area[0]).tolist(), max_aperture_profile_delta=profile_max,
            max_full_region_profile_delta=float(np.max(abs(delta))), iterations=ref['coupled']['iterations'], telemetry=ref['telemetry'])
        if name in ENDPOINTS:
            entry['matched_initial_populations_array_identity'] = True
            entry['external_rates'] = dict(schema='codex.bride_external_derived_artifact/1', resource_key=RESOURCE,
                relative_path=name+'/rates.npz', artifact=fingerprint(refined/'rates.npz'), source_report=fingerprint(refined/'report.json'),
                source_writer=ref['coupled']['helper'], scientific_admission='HOLD')
        result['cases'].append(entry)
        profiles.extend((a['normalized_flux'], b['normalized_flux'], delta))
        labels.extend((name+'_baseline', name+'_refined', name+'_delta'))
        axes[0].plot(wave, delta, label=name)
        axes[1].plot(CENTRES, entry['area_delta_mA'], 'o-', label=name)
    axes[0].set(xlim=(7770.5, 7776.5), xlabel='AIR wavelength (Angstrom)', ylabel='Refined minus baseline flux')
    axes[1].set(xlabel='O aperture centre, AIR (Angstrom)', ylabel='Refined minus baseline area (mA)')
    for ax in axes: ax.grid(alpha=.2); ax.legend(fontsize=8)
    fig.suptitle('Completed numerical refinement pairs; failures retained in report\nFixed source physics, no abundance/parameter fit; scientific HOLD')
    fig.tight_layout()
    args.output.mkdir(parents=True, exist_ok=False)
    fig.savefig(args.output/'endpoint_refinement.png', dpi=150); plt.close(fig)
    np.savetxt(args.output/'profiles.csv', np.column_stack(profiles), delimiter=',', header=','.join(labels), comments='')
    result['figure'] = fingerprint(args.output/'endpoint_refinement.png')
    result['profiles'] = fingerprint(args.output/'profiles.csv')
    (args.output/'report.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    for entry in result['cases']:
        if 'external_rates' in entry:
            (args.output/(entry['case']+'_rates_external.json')).write_text(json.dumps(entry['external_rates'], indent=2, allow_nan=False)+'\n')
    print(json.dumps(result['cases']))


if __name__ == '__main__': main()
