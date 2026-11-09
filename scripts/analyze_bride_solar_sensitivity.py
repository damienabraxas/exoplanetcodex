"""Fresh rate replay and fixed-observable response of independent solar T/xi probes."""
import argparse,csv,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    from config.constants import require_codex_path
    from pipeline.bride_reference import fingerprint
    from pipeline.bride_rate_audit import process_flow
    from scripts.analyze_bride_O_shape_milestone import shape_metrics
    from scripts.analyze_bride_H_kitt_orch import validate_processed_flux
    plan_path=ROOT/'config/experiments/rya1281_solar_sensitivity.json';plan=json.loads(plan_path.read_text());checks=[]
    def check(path,receipt):
        if fingerprint(path)!=receipt:raise ValueError('Changed sensitivity source '+str(path))
        checks.append((path,receipt))
    check(plan_path,fingerprint(plan_path));check(Path(__file__),fingerprint(__file__))
    shape_path=ROOT/'config/experiments/rya1281_O_shape_milestone.json';shape=json.loads(shape_path.read_text());check(shape_path,fingerprint(shape_path))
    parent_path=ROOT/'data/results/rya1281/H_kitt_orch_v1/report.json';parent=json.loads(parent_path.read_text());check(parent_path,fingerprint(parent_path))
    sys.path.insert(0,str(require_codex_path('engines.ispec')));import ispec,inspect
    check(Path(inspect.getfile(ispec.apply_post_fundamental_effects)),parent['ispec_effects'])
    op=ROOT/'data/results/rya1281/O_shape_milestone_v1/run/metrics.csv';check(op,fingerprint(op))
    observed={(r['rv'],float(r['rest_air_A']),r['scheme'],r['label']):r for r in csv.DictReader(op.open())if r['kind']=='observation'}
    cases=[c for c in plan['cases']if c['name']!='baseline'];rates=[];metrics=[];comparisons=[];profiles={};reports=[];fractions=None
    for case in cases:
        name=case['name'];folder=a.source/name;r=json.loads((folder/'report.json').read_text());initial=json.loads((folder/'input.json').read_text())
        if r['case_parameters']!=case or r['plan']!=fingerprint(plan_path)or any(r.get(k)!=v for k,v in initial.items())or r.get('dry_run'):raise ValueError('Wrong/mutated source case')
        if r['controls']!=plan['controls']or r['iteration_plugin']['receipt']!=plan['iteration_plugin']['build_receipt']:raise ValueError('Changed solver controls/plugin')
        for f,receipt in [('input.json',fingerprint(folder/'input.json')),('report.json',fingerprint(folder/'report.json')),('flux.csv',r['flux']),('populations.npz',r['coupled']['populations']),('independent_rates.npz',r['independent']['rates'])]:check(folder/f,receipt)
        with np.load(folder/'populations.npz')as data:
            pop=data['populations'];nstar=data['nstar'];seed=data['initial_populations'];current=seed/seed.sum(axis=0)
            if any(not np.isfinite(x).all()or np.any(x<=0)for x in [pop,nstar,seed]):raise ValueError('Invalid positive population state')
        if fractions is None:fractions=current.copy()
        if np.max(abs(current/fractions-1))>16*41*np.finfo(float).eps:raise ValueError('Different common initial fractions')
        with np.load(folder/'independent_rates.npz')as d:net,gross=process_flow(d['collision_rates_s_inv']+d['radiative_rates_s_inv'],pop)
        rate=float(np.max(abs(net)/np.maximum(gross,np.finfo(float).tiny)));conservation=float(np.max(abs(pop.sum(axis=0)/nstar.sum(axis=0)-1)))
        if rate>=r['controls']['rtol']or conservation>1e-10:raise ValueError('Independent rates/conservation failed')
        if reports:
            for key in ['atom','source_nodes','cap_length_m','background','backend','LW_constants','context_helper','state_helper','interpolator','density_helper','warm_report','warm_populations','source_collision_contract','iteration_plugin']:
                if r[key]!=reports[0][key]:raise ValueError('Changed common physics '+key)
        reports.append(r)
        if np.any(np.asarray(r['state']['vturb_m_s'])!=1000*case.get('xi_km_s',plan['xi_km_s'])):raise ValueError('Wrong microturbulence binding')
        if name.startswith('xi_'):
            for key,value in reports[0]['state'].items():
                if key!='vturb_m_s'and value!=r['state'][key]:raise ValueError('Xi probe changed another physical field '+key)
        rates.append(dict(case=name,teff_K=case['teff'],logg=case['logg'],xi_km_s=case.get('xi_km_s',plan['xi_km_s']),rate_residual=rate,conservation=conservation,iterations=r['coupled']['iterations'],seconds=r['telemetry']['seconds']))
        f=np.genfromtxt(folder/'flux.csv',delimiter=',',names=True);w=f['wavelength_air_A'];raw=f['normalized_flux'];np.testing.assert_array_equal(raw,f['Fnu_W_m2_Hz']/f['continuum_Fnu_W_m2_Hz'])
        if not np.isfinite(raw).all()or np.any(raw<=0):raise ValueError('Invalid source spectrum')
        if name=='solar':
            prior_folder=ROOT/'data/results/rya1281/solar_full_cap_v1/rya1281-solar-full-cap-v1/solar';prior=json.loads((prior_folder/'report.json').read_text());check(prior_folder/'flux.csv',prior['flux']);old=np.genfromtxt(prior_folder/'flux.csv',delimiter=',',names=True);np.testing.assert_array_equal(w,old['wavelength_air_A']);replay=float(np.max(abs(raw-old['normalized_flux'])))
            if replay>1e-6:raise ValueError('Solar endpoint replay screen failed')
        seg=np.array([(w[0]/10,w[-1]/10)],dtype=[('wave_base',float),('wave_top',float)])
        for rv in ['no_gravity','plus_gravity']:
            b=parent['broadening'];model=ispec.apply_post_fundamental_effects(w/10,raw.copy(),seg,macroturbulence=b['vmac_km_s'],vsini=b['vsini_km_s'],limb_darkening_coeff=b['limb_darkening'],R=b['R'],vrad=(parent['RV_km_s'][rv],));validate_processed_flux(w,model,shape['centres_air_A']);profiles[name+'_'+rv]=model
            beta=parent['RV_km_s'][rv]*1000/299792458.
            for rest in shape['centres_air_A']:
                centre=rest*np.sqrt((1+beta)/(1-beta))
                for scheme,shoulders in shape['shoulders'].items():
                    m=shape_metrics(w,model,centre,core_half_A=shape['core_half_A'],feature_half_A=shape['feature_half_A'],**shoulders)
                    row=dict(case=name,teff_K=case['teff'],xi_km_s=case.get('xi_km_s',plan['xi_km_s']),rv=rv,rest_air_A=rest,scheme=scheme,contrast=m['contrast'],area_mA=m['relative_feature_area_mA']);metrics.append(row)
                    for label in shape['observation_legs']:
                        o=observed[rv,rest,scheme,label];comparisons.append(dict(**row,observation=label,contrast_residual=m['contrast']-float(o['contrast']),area_residual_mA=m['relative_feature_area_mA']-float(o['relative_feature_area_mA'])))
    responses=[]
    for row in metrics:
        solar=next(r for r in metrics if r['case']=='solar'and all(row[k]==r[k]for k in ['rv','rest_air_A','scheme']))
        responses.append(dict(**row,area_change_mA=row['area_mA']-solar['area_mA'],contrast_change=row['contrast']-solar['contrast']))
    derivatives=[]
    for row in metrics:
        if row['case']!='solar':continue
        low,high=[next(r for r in metrics if r['case']==c and all(row[k]==r[k]for k in ['rv','rest_air_A','scheme']))for c in ['T_5762','T_5782']]
        derivatives.append(dict(rv=row['rv'],rest_air_A=row['rest_air_A'],scheme=row['scheme'],central_area_mA_per_K=(high['area_mA']-low['area_mA'])/20,lower_10K_area_change_mA=row['area_mA']-low['area_mA'],upper_10K_area_change_mA=high['area_mA']-row['area_mA'],second_difference_mA=high['area_mA']+low['area_mA']-2*row['area_mA']))
    for path,receipt in checks:
        if fingerprint(path)!=receipt:raise ValueError('Source changed during analysis')
    a.output.mkdir(parents=True,exist_ok=False)
    for name,rows in [('rates.csv',rates),('metrics.csv',responses),('comparisons.csv',comparisons),('temperature_derivatives.csv',derivatives)]:
        with(a.output/name).open('w',newline='')as f:writer=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');writer.writeheader();writer.writerows(rows)
    np.savez_compressed(a.output/'profiles.npz',wave_air_A=w,**profiles)
    import matplotlib
    matplotlib.use('Agg');import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,3,figsize=(14,4))
    for ax,rest in zip(axes,shape['centres_air_A']):
        rows=sorted([r for r in responses if r['rv']=='plus_gravity'and r['scheme']=='inner'and r['rest_air_A']==rest and r['teff_K']==5772],key=lambda x:x['xi_km_s'])
        ax.plot([r['xi_km_s']for r in rows],[r['area_change_mA']for r in rows],'o-');ax.axhline(0,color='gray',lw=.5);ax.set_title(str(rest)+' AIR A');ax.set_xlabel('RT microturbulence (km/s)');ax.set_ylabel('Relative feature-area change (mA)')
    fig.suptitle('Fixed atmosphere/abundance; independent NLTE solves, no best-fit selection');fig.tight_layout();fig.savefig(a.output/'microturbulence_response.png',dpi=150);plt.close(fig)
    report=dict(schema='codex.bride_solar_sensitivity/1',scientific_admission='HOLD',generator=fingerprint(__file__),plan=fingerprint(plan_path),scope=plan['scope'],cases=len(cases),solar_replay_max_flux_difference=replay,rates=rates,protected_inputs=[dict(path=str(p),receipt=r)for p,r in checks],metric_rows=len(responses),comparison_rows=len(comparisons),total_seconds=sum(r['seconds']for r in rates),limitations=['10K sampling is not10K physical accuracy.','Microturbulence changes RT Doppler velocities only, on a fixed xi1 source atmosphere.','Shoulder-relative diagnostics use unchanged shipped observations; separate ORCH continuum parity findings remain unresolved.','No inferred abundance, preferred input, physical3D/all-band or provider admission.'])
    (a.output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n');print(json.dumps({k:report[k]for k in ['cases','solar_replay_max_flux_difference','metric_rows','comparison_rows','total_seconds']}))

if __name__=='__main__':main()
