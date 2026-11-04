"""Audit native acceleration calls and independently verify bounded O controls."""
import argparse
import json
from pathlib import Path
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def join_updates(history,updates,mode):
    groups={i:[] for i in range(1,len(history)+1)}
    for call,r in enumerate(updates,1):
        if (r['schema']!='codex.bride_stat_equil_update/1' or r['call_index']!=call
                or r['outer_iteration'] not in groups or r['updatedPops'] is not True
                or len(r['active_species'])!=1):raise ValueError('Wrong native update accounting')
        for key in ('ngAccelerated','populations_finite','populations_positive'):
            if len(r[key])!=1 or not all(type(v) is bool for v in r[key]):raise ValueError('Missing native species flags')
        if updates and r['active_species']!=updates[0]['active_species']:raise ValueError('Changed active species')
        if mode=='off' and any(r['ngAccelerated']):raise ValueError('Disabled Ng actually accelerated')
        group=groups[r['outer_iteration']]
        if r['attempt_in_outer']!=len(group)+1:raise ValueError('Wrong retry index')
        group.append(r)
    joined=[]
    for i,h in enumerate(history,1):
        group=groups[i];expected=2 if h['rejected_invalid_ng'] else 1
        if h['iteration']!=i or len(group)!=expected:raise ValueError('Missing/unexplained native calls')
        valid=lambda r:all(r['populations_finite']) and all(r['populations_positive'])
        if not valid(group[-1]) or h['rejected_invalid_ng']==valid(group[0]):raise ValueError('Native rejection/acceptance mismatch')
        joint=max(h[k] for k in ('population_change','physical_rate_residual','J_change','cap_boundary_change'))
        joined.append(dict(iteration=i,accepted_ng=any(group[-1]['ngAccelerated']),native_calls=len(group),
            rejected_invalid_ng=h['rejected_invalid_ng'],joint_residual=joint,
            population_change=h['population_change'],physical_rate_residual=h['physical_rate_residual'],
            J_change=h['J_change'],cap_boundary_change=h['cap_boundary_change']))
    return joined


def exact_prefix(left,right):
    length=min(len(left),len(right))
    if not length:raise ValueError('Empty replay prefix')
    for i,(a,b) in enumerate(zip(left,right),1):
        if a!=b:raise ValueError(f'Recorded numerical prefix differs at iteration{i}')
    return length


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    from config.constants import require_codex_path
    from pipeline.bride_reference import fingerprint
    from pipeline.bride_thermal_sweep import effective_controls
    from scripts.analyze_bride_thermal_sweep import verify_case_binding
    from scripts.analyze_bride_thermal_refinement import verify_pair,independent_balance,CENTRES
    from pipeline.bride_observables import window_integral
    import matplotlib
    matplotlib.use('Agg');import matplotlib.pyplot as plt
    root=require_codex_path('experiments.bride_ng_controls')
    declaration=ROOT/'config/experiments/rya1281_ng_control.json';declared=json.loads(declaration.read_text())
    plan_path=ROOT/declared['physical_plan'];plan=json.loads(plan_path.read_text());case=declared['case']
    baseline=ROOT/'data/results/rya1281/thermal_slice_v1/runs'/case
    base=json.loads((baseline/'report.json').read_text())
    cert=require_codex_path('experiments.bride_thermal_coefficients')/case/'manifest.json'
    if base['coupled']['helper']!=fingerprint(ROOT/'pipeline/bride_coupled_control.py'):raise ValueError('Changed inherited solver')
    attempts=[json.loads(x) for x in (root/'attempts.jsonl').read_text().splitlines()]
    if [r['mode'] for r in attempts]!=declared['modes']:raise ValueError('Require both declared attempts, no favorable filter')
    report=dict(schema='codex.bride_ng_review/1',scientific_admission='HOLD',generator=fingerprint(__file__),
        declaration=fingerprint(declaration),baseline=fingerprint(baseline/'report.json'),attempts=fingerprint(root/'attempts.jsonl'),
        cases=[],scope='Numerical Ng policy intervention at one fixed lower-gravity atmosphere. Same physical operator, starting populations, sampling and acceptance limits. Not new atomic physics, an abundance fit or universal solver policy.')
    histories={};joined={};spectra={};csvs=[]
    for attempt in attempts:
        mode=attempt['mode'];d=root/mode
        if attempt['generator']!=fingerprint(ROOT/'scripts/run_bride_ng_controls.py') or attempt['plan']!=fingerprint(declaration) or attempt['log']!=fingerprint(root/(mode+'.log')):raise ValueError('Changed attempt receipt')
        inputs=json.loads((d/'input_manifest.json').read_text())
        controls=effective_controls(plan['controls'],True);controls['ng']=mode=='record'
        if (inputs['ng_control']!=mode or inputs['ng_experiment']!=fingerprint(declaration)
                or inputs['iteration_recorder']!=fingerprint(ROOT/'pipeline/bride_iteration_trace.py')
                or inputs['controls']!=controls or not inputs['numeric_refinement']
                or inputs['case']!=base['case'] or inputs['generator']!=fingerprint(ROOT/'scripts/run_bride_thermal_sweep.py')):raise ValueError('Changed requested Ng policy/input')
        for key in ('plan','sweep_adapter','certificate','atmosphere_sources','local_binding','source_manifest','source_binding','lines','pressure_source','stark_source'):
            if inputs[key]!=base[key]:raise ValueError('Ng comparison changes physical input: '+key)
        if inputs['plan']!=fingerprint(plan_path) or inputs['certificate']!=fingerprint(cert):raise ValueError('Changed case plan/certificate')
        verify_case_binding(inputs,json.loads(cert.read_text()))
        initial=inputs['initialization']
        if (initial['mode']!='matched_physical_baseline' or initial['populations']!=base['coupled']['populations']
                or initial['report']!=fingerprint(baseline/'report.json') or initial['helper']!=fingerprint(ROOT/'pipeline/bride_thermal_initial.py')):raise ValueError('Changed initial populations')
        h=[json.loads(x) for x in (d/'history.jsonl').read_text().splitlines()]
        u=[json.loads(x) for x in (d/'iteration_updates.jsonl').read_text().splitlines()]
        histories[mode]=h;joined[mode]=join_updates(h,u,mode)
        entry=dict(mode=mode,history=fingerprint(d/'history.jsonl'),updates=fingerprint(d/'iteration_updates.jsonl'),
            native_calls=len(u),native_ng_calls=sum(any(r['ngAccelerated']) for r in u),
            accepted_ng_steps=sum(r['accepted_ng'] for r in joined[mode]),iterations=len(h),
            initial=initial,wall_seconds=attempt['wall_seconds'],returncode=attempt['returncode'])
        if attempt['returncode']:
            last=(root/(mode+'.log')).read_text().rstrip().splitlines()[-1]
            if (d/'report.json').exists() or last not in ('TimeoutError: Fixed-population cap radiation budget exceeded','TimeoutError: Original coupled budget exceeded'):raise ValueError('Unexpected failed attempt; not a qualified budget refusal')
            entry.update(status='BUDGET_EXHAUSTED_NOT_ACCEPTED',last_joint_residual=joined[mode][-1]['joint_residual'])
        else:
            r=json.loads((d/'report.json').read_text());verify_pair(base,r,plan,ng_control=mode)
            if r['initialization']!=inputs['initialization'] or r['ng_experiment']!=inputs['ng_experiment']:
                raise ValueError('Input/report control provenance disagrees')
            for filename,key in (('input_manifest.json','input_manifest'),('flux.csv','flux')):
                if fingerprint(d/filename)!=r[key]:raise ValueError('Changed accepted input/flux')
            for filename,key in (('history.jsonl','history'),('populations.npz','populations'),('rates.npz','rates'),('iteration_updates.jsonl','iteration_updates')):
                if fingerprint(d/filename)!=r['coupled'][key]:raise ValueError('Changed coupled evidence')
            if (r['coupled']['recorder']!=fingerprint(ROOT/'pipeline/bride_iteration_trace.py') or r['coupled']['iterations']!=len(h)
                    or r['coupled']['final_diagnostics']!=h[-1] or joined[mode][-1]['joint_residual']>=controls['rtol']
                    or h[-1]['particle_conservation_error']>1e-10):raise ValueError('Joint acceptance check failed')
            with np.load(d/'populations.npz') as pop, np.load(baseline/'populations.npz') as old:
                np.testing.assert_array_equal(pop['initial_populations'],old['populations']);n=pop['populations'].copy()
            b=r['local_binding']['local_binding'];target=np.array(b['nH_m3'])*10**(b['scalars']['active_oxygen_abundance_dex']-12)
            with np.load(d/'rates.npz') as rates:
                balance=independent_balance(n,rates['radiative_rates_offdiag_s_inv'],rates['collision_rates_offdiag_s_inv'],rates['total_rates_s_inv'],target,controls['rtol'])
            data=np.genfromtxt(d/'flux.csv',delimiter=',',names=True)
            if not all(np.isfinite(data[k]).all() for k in data.dtype.names) or np.any(data['continuum_Fnu_W_m2_Hz']<=0):raise ValueError('Invalid accepted flux')
            np.testing.assert_array_equal(data['normalized_flux'],data['Fnu_W_m2_Hz']/data['continuum_Fnu_W_m2_Hz'])
            spectra[mode]=data
            entry.update(status='CONVERGED_EXPERIMENT',balance=balance,telemetry=r['telemetry'],
                recorder_overhead_seconds=r['coupled']['recorder_overhead_seconds'],source_report=fingerprint(d/'report.json'),
                intrinsic_area_mA=[1000*window_integral(data['wavelength_air_A'],data['normalized_flux'],c-.5,c+.5) for c in CENTRES],
                external_rates=dict(schema='codex.bride_external_derived_artifact/1',scientific_admission='HOLD',
                    generator=fingerprint(__file__),resource_key='experiments.bride_ng_controls',relative_path=mode+'/rates.npz',
                    artifact=fingerprint(d/'rates.npz'),source_report=fingerprint(d/'report.json'),source_writer=r['coupled']['helper']))
        report['cases'].append(entry)
    old_root=ROOT/'data/results/rya1281/thermal_extremes_v1'
    previous=json.loads((old_root/'analysis/report.json').read_text())
    previous_case=next(c for c in previous['cases'] if c['case']==case)
    old_path=old_root/'runs'/case/'history.jsonl'
    if fingerprint(old_path)!=previous_case['history']:raise ValueError('Changed previous timeout history')
    old=[json.loads(x) for x in old_path.read_text().splitlines()]
    report['legacy_prefix_reference']=dict(report=fingerprint(old_root/'analysis/report.json'),history=fingerprint(old_path))
    report['recorded_on_exact_legacy_prefix_rows']=exact_prefix(old,histories['record'])
    accelerated=[r['iteration'] for r in joined['record'] if r['accepted_ng']]
    report['first_actual_ng_iteration']=accelerated[0] if accelerated else None
    before=accelerated[0]-1 if accelerated else min(len(histories['record']),len(histories['off']))
    report['on_off_exact_prefix_before_first_ng_rows']=exact_prefix(histories['record'][:before],histories['off'][:before])
    report['first_on_off_different_iteration']=next((i for i,(a,b) in enumerate(zip(histories['record'],histories['off']),1) if a!=b),None)
    base_data=np.genfromtxt(baseline/'flux.csv',delimiter=',',names=True)
    for entry in report['cases']:
        if entry['mode'] in spectra:
            data=spectra[entry['mode']];np.testing.assert_array_equal(data['wavelength_air_A'],base_data['wavelength_air_A'])
            base_areas=[1000*window_integral(base_data['wavelength_air_A'],base_data['normalized_flux'],c-.5,c+.5) for c in CENTRES]
            entry['area_delta_vs_unrefined_mA']=(np.array(entry['intrinsic_area_mA'])-base_areas).tolist()
    if set(spectra)=={'record','off'}:
        a=spectra['record'];b=spectra['off'];np.testing.assert_array_equal(a['wavelength_air_A'],b['wavelength_air_A'])
        report['off_minus_on']=dict(max_normalized_flux_delta=float(np.max(abs(a['normalized_flux']-b['normalized_flux']))),
            area_delta_mA=(np.array(report['cases'][1]['intrinsic_area_mA'])-report['cases'][0]['intrinsic_area_mA']).tolist())
    args.output.mkdir(parents=True,exist_ok=False)
    fig,ax=plt.subplots(figsize=(11,5))
    for mode,rows in joined.items():
        ax.semilogy([r['iteration'] for r in rows],[r['joint_residual'] for r in rows],label=mode)
        marked=[r for r in rows if r['accepted_ng']]
        if marked:ax.scatter([r['iteration'] for r in marked],[r['joint_residual'] for r in marked],s=13,label=mode+' actual Ng')
        fields=list(rows[0]);path=args.output/(mode+'_iterations.csv')
        import csv
        with path.open('x') as stream:
            writer=csv.DictWriter(stream,fieldnames=fields);writer.writeheader();writer.writerows(rows)
        csvs.append(fingerprint(path))
    ax.axhline(plan['controls']['rtol'],color='black',ls='--',label='unchanged joint tolerance')
    ax.set(xlabel='Outer iteration',ylabel='Maximum joint residual',title='Same lower-gravity state: native Ng on/off, actual flags recorded\nTime-limited attempts remain explicit; scientific HOLD')
    ax.legend();ax.grid(alpha=.2);fig.tight_layout();figure=args.output/'ng_history.png';fig.savefig(figure,dpi=150);plt.close(fig)
    report.update(figure=fingerprint(figure),iteration_tables=csvs,
        caveats='One fixed atmosphere and start. Timing depends on host load. No universal acceleration/default policy, new physical model accuracy, abundance, continuum/telluric fit or 3D/all-band admission. Full unscaled rates and original particle/joint tolerances govern acceptance.')
    (args.output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    for entry in report['cases']:
        if 'external_rates' in entry:
            (args.output/(entry['mode']+'_rates_external.json')).write_text(json.dumps(entry['external_rates'],indent=2,allow_nan=False)+'\n')
    print(json.dumps(report['cases']))


if __name__=='__main__':main()
