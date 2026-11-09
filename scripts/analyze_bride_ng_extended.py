"""Verify a longer cost-only replay without weakening scientific acceptance."""
import argparse,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--attempt-log',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    from config.constants import require_codex_path
    from pipeline.bride_reference import fingerprint
    from pipeline.bride_thermal_sweep import effective_controls
    from pipeline.bride_compute_budget import extended_controls
    from scripts.analyze_bride_ng_controls import exact_prefix,join_updates
    from scripts.analyze_bride_thermal_refinement import verify_pair,independent_balance,CENTRES
    from scripts.analyze_bride_thermal_sweep import verify_case_binding
    from pipeline.bride_observables import window_integral
    import matplotlib
    matplotlib.use('Agg');import matplotlib.pyplot as plt

    d=require_codex_path('experiments.bride_ng_extended')
    declaration=ROOT/'config/experiments/rya1281_ng_extended.json';budget=json.loads(declaration.read_text())
    previous=ROOT/budget['reference_run'];old_inputs=json.loads((previous/'input_manifest.json').read_text())
    plan_path=ROOT/budget['physical_plan'];plan=json.loads(plan_path.read_text())
    inputs=json.loads((d/'input_manifest.json').read_text())
    for key in ('case','plan','sweep_adapter','certificate','atmosphere_sources','local_binding','source_manifest',
                'source_binding','initial_populations','initialization','lines','pressure_source','stark_source','ng_control','ng_experiment','iteration_recorder'):
        if inputs[key]!=old_inputs[key]:raise ValueError('Compute extension changes physical/source input: '+key)
    expected=extended_controls(old_inputs['controls'],budget,inputs['case']['name'],inputs['ng_control'])
    if (inputs['controls']!=expected or inputs['compute_extension']!=fingerprint(declaration)
            or inputs['budget_adapter']!=fingerprint(ROOT/'pipeline/bride_compute_budget.py')
            or inputs['generator']!=fingerprint(ROOT/'scripts/run_bride_thermal_sweep.py')):raise ValueError('Changed explicit cost-only controls/provenance')
    if {k for k in expected if expected[k]!=old_inputs['controls'][k]}!={'seconds_per_case'}:raise ValueError('More than compute cost changed')
    cert=require_codex_path('experiments.bride_thermal_coefficients')/budget['case']/'manifest.json'
    if inputs['certificate']!=fingerprint(cert):raise ValueError('Changed coefficient certificate')
    verify_case_binding(inputs,json.loads(cert.read_text()))
    history=[json.loads(x) for x in (d/'history.jsonl').read_text().splitlines()]
    updates=[json.loads(x) for x in (d/'iteration_updates.jsonl').read_text().splitlines()]
    old_h=[json.loads(x) for x in (previous/'history.jsonl').read_text().splitlines()]
    old_u=[json.loads(x) for x in (previous/'iteration_updates.jsonl').read_text().splitlines()]
    joined=join_updates(history,updates,'record')
    report=dict(schema='codex.bride_ng_extended_review/1',scientific_admission='HOLD',generator=fingerprint(__file__),
        declaration=fingerprint(declaration),input_manifest=fingerprint(d/'input_manifest.json'),attempt_log=fingerprint(args.attempt_log),
        history=fingerprint(d/'history.jsonl'),updates=fingerprint(d/'iteration_updates.jsonl'),
        reference_input=fingerprint(previous/'input_manifest.json'),reference_history=fingerprint(previous/'history.jsonl'),
        reference_updates=fingerprint(previous/'iteration_updates.jsonl'),
        exact_history_prefix_rows=exact_prefix(old_h,history),exact_update_prefix_calls=exact_prefix(old_u,updates),
        iterations=len(history),native_ng_calls=sum(any(r['ngAccelerated']) for r in updates),controls=expected,
        scope='One cost-only extension. Previous600s failures remain failures. Same physical equations, sampling, initial state and accuracy gates; no universal convergence or physical/abundance admission.')
    if report['exact_history_prefix_rows']!=len(old_h) or report['exact_update_prefix_calls']!=len(old_u):raise ValueError('Replay did not cover complete prior history')
    if not (d/'report.json').exists():
        last=args.attempt_log.read_text().rstrip().splitlines()[-1]
        if last not in ('TimeoutError: Original coupled budget exceeded','TimeoutError: Fixed-population cap radiation budget exceeded'):raise ValueError('Unexpected failure; not a qualified timeout')
        report.update(status='BUDGET_EXHAUSTED_NOT_ACCEPTED',last_joint_residual=joined[-1]['joint_residual'])
    else:
        r=json.loads((d/'report.json').read_text())
        baseline=ROOT/'data/results/rya1281/thermal_slice_v1/runs'/budget['case'];base=json.loads((baseline/'report.json').read_text())
        verify_pair(base,r,plan,ng_control='record',extended_budget=True)
        if (r['initialization']!=inputs['initialization'] or r['coupled']['controls']!=inputs['controls']
                or r['compute_extension']!=inputs['compute_extension'] or r['budget_adapter']!=inputs['budget_adapter']):raise ValueError('Input/report provenance differs')
        for name,key in (('input_manifest.json','input_manifest'),('flux.csv','flux'),('widths.csv','widths'),('neutral_colliders.csv','neutral_colliders'),('neutral_reservoir.csv','neutral_reservoir')):
            if fingerprint(d/name)!=r[key]:raise ValueError('Changed accepted evidence: '+name)
        for name,key in (('history.jsonl','history'),('iteration_updates.jsonl','iteration_updates'),('populations.npz','populations'),('rates.npz','rates')):
            if fingerprint(d/name)!=r['coupled'][key]:raise ValueError('Changed coupled evidence: '+name)
        if (r['coupled']['final_diagnostics']!=history[-1] or r['coupled']['iterations']!=len(history)
                or joined[-1]['joint_residual']>=expected['rtol'] or history[-1]['particle_conservation_error']>1e-10):raise ValueError('Joint convergence not accepted')
        with np.load(d/'populations.npz') as saved,np.load(baseline/'populations.npz') as original:
            np.testing.assert_array_equal(saved['initial_populations'],original['populations']);pop=saved['populations'].copy()
        b=r['local_binding']['local_binding'];total=np.array(b['nH_m3'])*10**(b['scalars']['active_oxygen_abundance_dex']-12)
        with np.load(d/'rates.npz') as rates:
            balance=independent_balance(pop,rates['radiative_rates_offdiag_s_inv'],rates['collision_rates_offdiag_s_inv'],rates['total_rates_s_inv'],total,expected['rtol'])
        data=np.genfromtxt(d/'flux.csv',delimiter=',',names=True);before=np.genfromtxt(baseline/'flux.csv',delimiter=',',names=True)
        np.testing.assert_array_equal(data['wavelength_air_A'],before['wavelength_air_A'])
        if not all(np.isfinite(data[k]).all() for k in data.dtype.names) or np.any(data['continuum_Fnu_W_m2_Hz']<=0):raise ValueError('Invalid flux')
        np.testing.assert_array_equal(data['normalized_flux'],data['Fnu_W_m2_Hz']/data['continuum_Fnu_W_m2_Hz'])
        wave=data['wavelength_air_A'];delta=data['normalized_flux']-before['normalized_flux']
        areas=[[1000*window_integral(wave,x['normalized_flux'],c-.5,c+.5) for c in CENTRES] for x in (before,data)]
        report.update(status='CONVERGED_EXPERIMENT',source_report=fingerprint(d/'report.json'),balance=balance,
            initial_array_identity=True,baseline_area_mA=areas[0],refined_area_mA=areas[1],
            area_delta_mA=(np.array(areas[1])-areas[0]).tolist(),max_normalized_profile_delta=float(np.max(abs(delta))),
            max_aperture_profile_delta=[float(np.max(abs(delta[(wave>=c-.5)&(wave<=c+.5)]))) for c in CENTRES],
            telemetry=r['telemetry'],recorder_overhead_seconds=r['coupled']['recorder_overhead_seconds'],
            external_rates=dict(schema='codex.bride_external_derived_artifact/1',scientific_admission='HOLD',
                resource_key='experiments.bride_ng_extended',relative_path='rates.npz',artifact=fingerprint(d/'rates.npz'),
                source_report=fingerprint(d/'report.json'),source_writer=r['coupled']['helper'],generator=fingerprint(__file__)))
    args.output.mkdir(parents=True,exist_ok=False)
    fig,ax=plt.subplots(figsize=(11,5));ax.semilogy([r['iteration'] for r in joined],[r['joint_residual'] for r in joined],label='1200s allowance')
    flagged=[r for r in joined if r['accepted_ng']]
    ax.scatter([r['iteration'] for r in flagged],[r['joint_residual'] for r in flagged],s=12,label='Actual native Ng')
    ax.axvline(len(old_h),color='grey',ls=':',label='Earlier600s stop, exact prefix')
    ax.axhline(expected['rtol'],color='black',ls='--',label='Unchanged joint criterion')
    ax.set(xlabel='Outer iteration',ylabel='Maximum joint residual',title='Lower-gravity cost-only replay: '+report['status']+'\nFixed accuracy/physical inputs; scientific HOLD')
    ax.legend();ax.grid(alpha=.2);fig.tight_layout();figure=args.output/'extended_history.png';fig.savefig(figure,dpi=150);plt.close(fig)
    report['figure']=fingerprint(figure)
    if report['status']=='CONVERGED_EXPERIMENT':
        fig,axes=plt.subplots(2,1,figsize=(11,7),sharex=True)
        axes[0].plot(wave,before['normalized_flux'],label='Unrefined baseline')
        axes[0].plot(wave,data['normalized_flux'],ls='--',label='Refined, longer cost allowance')
        axes[1].plot(wave,delta)
        for ax in axes:ax.set_xlim(7770.5,7776.5);ax.grid(alpha=.2)
        axes[0].legend();axes[0].set_ylabel('Intrinsic normalized model flux')
        axes[1].set(xlabel='AIR wavelength (Angstrom)',ylabel='Refined minus baseline')
        fig.suptitle('Same lower-gravity physical state; numerical comparison only\nNo observed spectrum/continuum fit; scientific HOLD')
        fig.tight_layout();path=args.output/'extended_spectrum.png';fig.savefig(path,dpi=150);plt.close(fig)
        report['spectrum_figure']=fingerprint(path)
    import csv
    path=args.output/'iterations.csv'
    with path.open('x') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(joined[0]),lineterminator='\n');writer.writeheader();writer.writerows(joined)
    report['iteration_table']=fingerprint(path)
    (args.output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    if 'external_rates' in report:(args.output/'rates_external.json').write_text(json.dumps(report['external_rates'],indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:report[k] for k in ('status','iterations','native_ng_calls','exact_history_prefix_rows','balance','area_delta_mA','telemetry') if k in report}))


if __name__=='__main__':main()
