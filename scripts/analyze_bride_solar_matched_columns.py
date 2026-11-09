"""Fresh four-column O rate admission, solar shape and finite T/g accounting."""
import argparse,csv,inspect,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--run-plan',type=Path,default=ROOT/'config/experiments/rya1281_solar_matched_columns.json');a=p.parse_args()
 from config.constants import require_codex_path
 from pipeline.bride_reference import fingerprint
 from pipeline.bride_rate_audit import process_flow
 from pipeline.bride_observables import window_integral
 from pipeline.bride_crossed_continuum import finite_accounting
 from scripts.analyze_bride_O_shape_milestone import shape_metrics
 from scripts.analyze_bride_H_kitt_orch import validate_processed_flux
 plan_path=ROOT/'config/experiments/rya1281_solar_matched_analysis.json';plan=json.loads(plan_path.read_text());run_plan=a.run_plan;rp=json.loads(run_plan.read_text());checks=[]
 def check(path,receipt):
  if fingerprint(path)!=receipt:raise ValueError('Changed matched analysis input '+str(path))
  checks.append((path,receipt))
 for name,receipt in plan['source_files'].items():check(ROOT/name,receipt)
 reports=[];populations=[];fluxes=[];rates=[]
 for case in plan['case_order']:
  folder=a.source/case;r=json.loads((folder/'report.json').read_text());initial=json.loads((folder/'input.json').read_text())
  if r['case']!=case or r.get('dry_run',False)or any(r.get(k)!=v for k,v in initial.items())or r['case_parameters']!=next(c for c in rp['cases']if c['name']==case)or r['plan']!=fingerprint(run_plan):raise ValueError('Wrong/mutated matched solve')
  if r['controls']!=rp['controls']or r.get('numerical_algorithm')!=rp.get('numerical_algorithm'):raise ValueError('Wrong numerical control')
  if rp.get('iteration_plugin')and r.get('iteration_plugin',{}).get('receipt')!=rp['iteration_plugin']['build_receipt']:raise ValueError('Wrong iteration plugin build')
  for name,receipt in [('report.json',fingerprint(folder/'report.json')),('flux.csv',r['flux']),('populations.npz',r['coupled']['populations']),('history.jsonl',r['coupled']['history']),('independent_rates.npz',r['independent']['rates'])]:check(folder/name,receipt)
  with np.load(folder/'populations.npz')as data:pop={k:data[k].copy()for k in data.files}
  if pop['populations'].shape!=(41,56)or any(not np.isfinite(v).all()or np.any(v<=0)for v in pop.values()):raise ValueError('Invalid final populations')
  with np.load(folder/'independent_rates.npz')as data:net,gross=process_flow(data['collision_rates_s_inv']+data['radiative_rates_s_inv'],pop['populations'])
  residual=float(np.max(abs(net)/np.maximum(gross,np.finfo(float).tiny)));conservation=float(np.max(abs(pop['populations'].sum(axis=0)/pop['nstar'].sum(axis=0)-1)))
  if residual>=r['controls']['rtol']or conservation>1e-10:raise ValueError('Fresh physical rate/conservation replay refused')
  f=np.genfromtxt(folder/'flux.csv',delimiter=',',names=True)
  if not all(np.isfinite(f[k]).all()for k in f.dtype.names)or any(np.any(f[k]<=0)for k in ['Fnu_W_m2_Hz','continuum_Fnu_W_m2_Hz']):raise ValueError('Invalid matched flux')
  np.testing.assert_array_equal(f['normalized_flux'],f['Fnu_W_m2_Hz']/f['continuum_Fnu_W_m2_Hz']);reports.append(r);populations.append(pop);fluxes.append(f);rates.append(dict(case=case,rate_residual=residual,conservation=conservation,iterations=r['coupled']['iterations'],seconds=r['telemetry']['seconds'],peak_rss_bytes=r['telemetry']['process_peak_rss_bytes']))
 base=reports[0];mode=rp.get('initialization_mode','common_fractions');seed0=populations[0]['initial_populations']/(populations[0]['nstar']if mode=='common_departure_pattern'else 1);base_fraction=seed0/seed0.sum(axis=0);seed_deltas=[]
 for r,pop,f in zip(reports,populations,fluxes):
  if r.get('iteration_plugin')!=base.get('iteration_plugin')or r['coupled']['helper']!=base['coupled']['helper']:raise ValueError('Changed numerical helper/plugin across matrix')
  for key in ['atom','source_nodes','column_preparation','field_ownership','controls','cap_length_m','background','backend','LW_constants','context_helper','state_helper','interpolator','density_helper','warm_report','warm_populations','source_collision_contract']:
   if r[key]!=base[key]:raise ValueError('Changed common matched input '+key)
  if r['state']['case']!='model_mass_H'or r['state']['ne_owner']!=base['state']['ne_owner']or r['state']['nH_owner']!=base['state']['nH_owner']or r['state']['vturb_m_s']!=base['state']['vturb_m_s']:raise ValueError('Changed field ownership/microturbulence')
  np.testing.assert_array_equal(f['wavelength_air_A'],fluxes[0]['wavelength_air_A']);seed=pop['initial_populations']/(pop['nstar']if mode=='common_departure_pattern'else 1);fraction=seed/seed.sum(axis=0);delta=float(np.max(abs(fraction/base_fraction-1)));seed_deltas.append(delta)
  if delta>16*41*np.finfo(float).eps:raise ValueError('Different common starting fractions')
 old=np.genfromtxt(ROOT/'data/results/rya1281/H_nlte_pair_v2/model_mass_H/flux.csv',delimiter=',',names=True);np.testing.assert_array_equal(old['wavelength_air_A'],fluxes[0]['wavelength_air_A']);replay=float(np.max(abs(old['normalized_flux']-fluxes[0]['normalized_flux'])));creplay=float(np.max(abs(fluxes[0]['continuum_Fnu_W_m2_Hz']/old['continuum_Fnu_W_m2_Hz']-1)))
 if replay>plan['baseline_max_normalized_flux_difference']or creplay>plan['baseline_max_relative_continuum_difference']:raise ValueError('Source precision baseline replay screen failed')
 endpoint_replays=[]
 for endpoint in plan['endpoint_controls']:
  directory=ROOT/endpoint['path'];legacy=json.loads((directory/'report.json').read_text());i=plan['case_order'].index(endpoint['case']);current=reports[i]
  for key in ['state','case_parameters','source_nodes','cap_length_m','atom','context_helper','backend','LW_constants','source_collision_contract']:
   if legacy[key]!=current[key]:raise ValueError('Endpoint physical inputs differ '+key)
  expected_controls=dict(legacy['controls'])
  if rp.get('numerical_algorithm',{}).get('method')in ['anderson_physical_direct','transition_local_safeguarded_anderson']:expected_controls['ng']=False
  if rp.get('numerical_algorithm',{}).get('method')=='transition_local_mali'and rp['numerical_algorithm'].get('native_ng')is False:expected_controls['ng']=False
  if current['controls']!=expected_controls:raise ValueError('Endpoint controls differ beyond declared acceleration')
  if legacy['independent']['physical_rate_residual']>=legacy['controls']['rtol']:raise ValueError('Unadmitted legacy endpoint')
  check(directory/'flux.csv',legacy['flux']);prior=np.genfromtxt(directory/'flux.csv',delimiter=',',names=True);np.testing.assert_array_equal(prior['wavelength_air_A'],fluxes[i]['wavelength_air_A']);df=float(np.max(abs(prior['normalized_flux']-fluxes[i]['normalized_flux'])));dc=float(np.max(abs(fluxes[i]['continuum_Fnu_W_m2_Hz']/prior['continuum_Fnu_W_m2_Hz']-1)))
  if df>plan['endpoint_max_normalized_flux_difference']or dc>plan['endpoint_max_relative_continuum_difference']:raise ValueError('Legacy endpoint numerical replay screen failed')
  endpoint_replays.append(dict(case=endpoint['case'],legacy_report=fingerprint(directory/'report.json'),max_normalized_flux_difference=df,max_relative_continuum_difference=dc,scope='Same physical inputs under separate initialization/operator; numerical screen only, not independent physical accuracy'))
 shape_plan=json.loads((ROOT/'config/experiments/rya1281_O_shape_milestone.json').read_text());parent=json.loads((ROOT/'data/results/rya1281/H_kitt_orch_v1/report.json').read_text());check(require_codex_path('data.kurucz2005_residual'),parent['observation']);check(require_codex_path('data.solar_iag_baker2020'),parent['iag'])
 sys.path.insert(0,str(require_codex_path('engines.ispec')));import ispec
 check(Path(inspect.getfile(ispec.apply_post_fundamental_effects)),parent['ispec_effects'])
 with np.load(ROOT/'data/results/rya1281/O_shape_milestone_v1/run/profiles.npz')as data:observed_profiles={k:data[k].copy()for k in data.files if k.startswith('observed_')or k.startswith('wave_')}
 with(ROOT/'data/results/rya1281/O_shape_milestone_v1/run/metrics.csv').open()as stream:observed={(r['rv'],float(r['rest_air_A']),r['scheme'],r['label']):r for r in csv.DictReader(stream)if r['kind']=='observation'}
 mw=fluxes[0]['wavelength_air_A'];segment=np.array([(mw[0]/10,mw[-1]/10)],dtype=[('wave_base',float),('wave_top',float)]);b=parent['broadening'];models={};profiles={};metrics=[];comparisons=[];intrinsic=[]
 for rv in plan['RVs']:
  beta=parent['RV_km_s'][rv]*1000/299792458.
  for case,raw in zip(plan['case_order'],fluxes):
   model=ispec.apply_post_fundamental_effects(mw/10,raw['normalized_flux'].copy(),segment,macroturbulence=b['vmac_km_s'],vsini=b['vsini_km_s'],limb_darkening_coeff=b['limb_darkening'],R=b['R'],vrad=(parent['RV_km_s'][rv],));validate_processed_flux(mw,model,plan['centres_air_A']);models[rv,case]=model
   for ci,rest in enumerate(plan['centres_air_A']):
    centre=float(rest*np.sqrt((1+beta)/(1-beta)));grid=observed_profiles[f'wave_{rv}_{ci}'];offset=grid-centre;core=abs(offset)<=shape_plan['core_half_A']+1e-12;mf=np.interp(grid,mw,model);profiles[f'model_{rv}_{ci}_{case}']=mf;profiles[f'wave_{rv}_{ci}']=grid
    for scheme in plan['shoulders']:
     m=shape_metrics(mw,model,centre,core_half_A=shape_plan['core_half_A'],feature_half_A=shape_plan['feature_half_A'],**shape_plan['shoulders'][scheme]);metrics.append(dict(case=case,rv=rv,rest_air_A=rest,scheme=scheme,contrast=m['contrast'],relative_core_area_mA=m['relative_core_area_mA'],relative_feature_area_mA=m['relative_feature_area_mA'],shoulder_level=m['shoulder_level']))
     for observation in shape_plan['observation_legs']:
      om=observed[rv,rest,scheme,observation];of=observed_profiles[f'observed_{rv}_{ci}_{observation}'];profiles[f'observed_{rv}_{ci}_{observation}']=of;comparisons.append(dict(case=case,rv=rv,rest_air_A=rest,scheme=scheme,observation=observation,contrast_residual=m['contrast']-float(om['contrast']),relative_feature_area_residual_mA=m['relative_feature_area_mA']-float(om['relative_feature_area_mA']),raw_core_rms=float(np.sqrt(np.mean((mf[core]-of[core])**2))),shoulder_ratio_core_rms=float(np.sqrt(np.mean((mf[core]/m['shoulder_level']-of[core]/float(om['shoulder_level']))**2)))))
 for case,raw in zip(plan['case_order'],fluxes):
  for rest in plan['centres_air_A']:intrinsic.append(dict(case=case,rest_air_A=rest,area_mA=1000*window_integral(mw,raw['normalized_flux'],rest-.5,rest+.5)))
 factorial=[]
 for rv in plan['RVs']:
  for rest in plan['centres_air_A']:
   for scheme in plan['shoulders']:
    by={r['case']:r for r in metrics if r['rv']==rv and r['rest_air_A']==rest and r['scheme']==scheme}
    for statistic in ['contrast','relative_feature_area_mA']:
     corners={key:dict(value=by[case][statistic])for key,case in [('before_before','baseline'),('after_before','temperature'),('before_after','gravity'),('after_after','solar')]};change=finite_accounting(corners,'value')
     if not change['available']or abs(change['closure'])>1e-10:raise ValueError('Finite T/g accounting failed')
     factorial.append(dict(rv=rv,rest_air_A=rest,scheme=scheme,statistic=statistic,total=change['total'],temperature_path_average=change['support_path_average'],gravity_path_average=change['denominator_path_average'],interaction=change['interaction'],closure=change['closure']))
 if(len(metrics),len(comparisons),len(factorial))!=(48,144,24):raise ValueError('Incomplete matched analysis')
 for path,receipt in checks:
  if fingerprint(path)!=receipt:raise ValueError('Protected source changed during matched analysis')
 a.output.mkdir(parents=True,exist_ok=False)
 for name,rows in [('rates.csv',rates),('metrics.csv',metrics),('comparisons.csv',comparisons),('factorial.csv',factorial),('intrinsic.csv',intrinsic)]:
  with(a.output/name).open('w',newline='')as stream:writer=csv.DictWriter(stream,fieldnames=list(rows[0]),lineterminator='\n');writer.writeheader();writer.writerows(rows)
 np.savez(a.output/'profiles.npz',**profiles)
 import matplotlib
 matplotlib.use('Agg');import matplotlib.pyplot as plt
 fig,axes=plt.subplots(2,3,figsize=(15,8));rv='plus_gravity';scheme='inner'
 for ci,rest in enumerate(plan['centres_air_A']):
  grid=profiles[f'wave_{rv}_{ci}'];beta=parent['RV_km_s'][rv]*1000/299792458.;offset=grid-rest*np.sqrt((1+beta)/(1-beta))
  for observation in shape_plan['observation_legs']:
   f=profiles[f'observed_{rv}_{ci}_{observation}'];om=observed[rv,rest,scheme,observation];axes[0,ci].plot(offset,f/float(om['shoulder_level']),label=observation)
  for case in plan['case_order']:
   f=profiles[f'model_{rv}_{ci}_{case}'];m=next(x for x in metrics if x['case']==case and x['rv']==rv and x['rest_air_A']==rest and x['scheme']==scheme);axes[0,ci].plot(offset,f/m['shoulder_level'],label=case)
   baseline=next(x for x in metrics if x['case']=='baseline'and x['rv']==rv and x['rest_air_A']==rest and x['scheme']==scheme);axes[1,ci].plot(offset,f/m['shoulder_level']-profiles[f'model_{rv}_{ci}_baseline']/baseline['shoulder_level'],label=case)
  axes[0,ci].set_title(f'{rest:.3f} AIR A');axes[0,ci].set_ylabel('Fixed-shoulder ratio');axes[1,ci].set_ylabel('Case minus matched baseline ratio')
  for ax in axes[:,ci]:ax.set_xlabel('AIR offset from fixed +gravity centre (A)');ax.legend(fontsize=7)
 fig.suptitle('Matched native/T-only/g-only/solar O columns: source abundance unchanged\nFixed cap/physics; no optimized model, broadening, RV or continuum');fig.tight_layout();fig.savefig(a.output/'matched_columns.png',dpi=140);plt.close(fig)
 result=dict(schema='codex.bride_solar_matched_analysis/1',scientific_admission='HOLD',scope=plan['scope'],generator=fingerprint(__file__),plan=fingerprint(plan_path),run_plan=fingerprint(run_plan),source_reports=[fingerprint(a.source/case/'report.json')for case in plan['case_order']],rates=rates,source_files=plan['source_files'],baseline_replay=dict(max_absolute_normalized_flux_difference=replay,max_relative_continuum_difference=creplay,source_field_changes=base['binding_change_from_prior']),endpoint_replays=endpoint_replays,initialization_mode=mode,max_common_seed_pattern_difference=max(seed_deltas),rows=dict(metrics=len(metrics),comparisons=len(comparisons),factorial=len(factorial)),total_solved_seconds=sum(r['seconds']for r in rates),limitations=['Exact canonical labels are interpolated input controls, not independent atmosphere/EOS/interpolation accuracy or physical model validation.','Fixed original geometric cap is an experimental control, not an interpolated/validated column mass.','Only this model-H source ownership; exported native-network packet not transplanted.','Finite T/g effects/interactions are correlated controlled differences, not independent uncertainty components or inverse parameter fits.','Current source atom/collision/cap/MULTI/opacity/reference/3D/all-band and CNO gates remain open; no adopted abundance/model/continuum/production.'])
 (a.output/'report.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(json.dumps(dict(baseline_replay=result['baseline_replay'],rates=rates,total_seconds=result['total_solved_seconds'])))

if __name__=='__main__':main()
