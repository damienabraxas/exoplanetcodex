"""Matched four-column solar O NLTE experiment with unchanged numerical controls."""
import argparse,json,os,platform,resource,sys,time
from dataclasses import asdict
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--case',choices=['baseline','temperature','gravity','solar'],required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--dry-run',action='store_true');p.add_argument('--plan',type=Path,default=ROOT/'config/experiments/rya1281_solar_matched_defect.json');a=p.parse_args()
 import lightweaver as lw
 import lightweaver.constants as c
 import lightweaver.LwCompiled as compiled
 from lightweaver.atomic_table import PeriodicTable
 from importlib.metadata import version
 from config.constants import require_codex_path,get_star_params
 from pipeline.bride_reference import fingerprint,read_marcs_structure
 from pipeline.multi_atom import read_multi_atom
 from pipeline.bride_solar_matched import prepare_columns,column_state,departure_seed
 from pipeline.bride_H_nlte_control import make_context,scale_initial_fractions
 from pipeline.bride_coupled_defect_control import solve_coupled_defect as solve_coupled
 from pipeline.bride_top_cap import converge_cap_radiation,emergent_rays
 from pipeline.bride_observables import flux_quadrature,integrate_flux
 from pipeline.bride_rate_audit import collision_classes,process_flow
 from pipeline.wavelength_util import air_to_vac
 if platform.system()!='Darwin'or version('lightweaver')!='0.16.1':raise RuntimeError('Mac/LW0.16.1 required')
 plan_path=a.plan;plan=json.loads(plan_path.read_text());warm_report=ROOT/plan['warm_report'];warm=json.loads(warm_report.read_text());controls=plan['controls']
 if controls!=warm['controls']:raise ValueError('Changed inherited joint solver controls')
 checks=[]
 def check(path,receipt):
  if fingerprint(path)!=receipt:raise ValueError('Changed predeclared input '+str(path))
  checks.append((path,receipt))
 check(warm_report,fingerprint(warm_report));check(plan_path,fingerprint(plan_path));check(Path(__file__),fingerprint(__file__))
 warm_path=ROOT/plan['warm_populations'];check(warm_path,warm['coupled']['populations']);check(ROOT/'pipeline/bride_coupled_control.py',warm['coupled']['helper'])
 with np.load(warm_path)as data:pop=data['populations'].copy();seed_nstar=data['nstar'].copy()
 check(Path(compiled.__file__),warm['backend']);check(Path(c.__file__),warm['LW_constants'])
 solar=get_star_params('solar');target=next(x for x in plan['cases']if x['name']=='solar')
 if(target['teff'],target['logg'],plan['xi_km_s'])!=(solar['teff'],solar['logg'],solar['xi']):raise ValueError('Changed canonical solar experiment target')
 columns,preparation=prepare_columns(plan,require_codex_path('grids.marcs_standard'));structure=columns[a.case];state=column_state(structure);case=next(x for x in plan['cases']if x['name']==a.case)
 from pipeline.bride_parameter_checks import marcs_filename
 for node in plan['nodes']:check(require_codex_path('grids.marcs_standard')/marcs_filename(node['teff'],node['logg']),node['source'])
 for helper in ['pipeline/bride_defect_correction.py','pipeline/bride_coupled_defect_control.py','pipeline/bride_solar_matched.py','pipeline/bride_parameter_checks.py','pipeline/bride_reference.py','pipeline/bride_H_nlte_control.py']:check(ROOT/helper,fingerprint(ROOT/helper))
 binding_change={key:dict(max_absolute=float(np.max(abs(np.array(value)-np.array(warm['state'][key])))),max_relative=None if key=='height_m'else float(np.max(abs(np.array(value)/np.array(warm['state'][key])-1))))for key,value in asdict(state).items()if isinstance(value,tuple)}
 atom_path=require_codex_path('experiments.bride_oxygen_atom');check(atom_path,warm['atom']);atom=read_multi_atom(atom_path)
 if atom.abundance!=8.70 or len(atom.levels)!=41 or pop.shape!=(41,56):raise ValueError('Changed source atom/warm layout')
 lo,hi=plan['range_air_A'];wave=np.linspace(lo,hi,int(round((hi-lo)/plan['step_air_A']))+1);vac=air_to_vac(wave)/10;extra=np.array([vac[0],vac[-1]])
 cap_length=plan['cap_length_m']
 if cap_length!=warm['cap_length_m']or not np.isfinite(cap_length)or cap_length<=0:raise ValueError('Changed fixed geometric cap policy')
 a.output.mkdir(parents=True,exist_ok=False)
 result=dict(schema='codex.bride_solar_matched_column/1',scientific_admission='HOLD',case=a.case,scope=plan['scope'],generator=fingerprint(__file__),plan=fingerprint(plan_path),state=asdict(state),state_helper=fingerprint(ROOT/'pipeline/bride_solar_matched.py'),context_helper=fingerprint(ROOT/'pipeline/bride_H_nlte_control.py'),case_parameters=case,source_nodes=preparation['sources'],column_preparation=preparation,field_ownership=plan['field_ownership'],binding_change_from_prior= binding_change,density_helper=fingerprint(ROOT/'pipeline/bride_reference.py'),interpolator=fingerprint(ROOT/'pipeline/bride_parameter_checks.py'),atom=fingerprint(atom_path),warm_report=fingerprint(warm_report),warm_populations=fingerprint(warm_path),initialization=plan['initialization'],initialization_mode=plan.get('initialization_mode','common_fractions'),solver_operator=plan['solver_operator'],controls=controls,cap_length_m=cap_length,background=plan['background'],backend=fingerprint(compiled.__file__),LW_constants=fingerprint(c.__file__))
 (a.output/'input.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');started=time.perf_counter()
 def construct(*,fixed=None,continuum_only=False):return make_context(state,structure,read_multi_atom(atom_path),rays=controls['rays'],qscale=controls['qscale'],extra_vac_nm=extra,cap_length_m=cap_length,ng=controls['ng'],fixed=fixed,continuum_only=continuum_only)
 ctx,audits=construct();active=ctx.activeAtoms[0];nstar=active.nStar.copy()
 mode=plan.get('initialization_mode','common_fractions')
 if mode=='common_fractions':active.n[:]=scale_initial_fractions(pop,active.nTotal)
 elif mode=='common_departure_pattern':active.n[:]=departure_seed(pop,seed_nstar,nstar,active.nTotal)
 else:raise ValueError('Unknown predeclared initialization mode')
 result['initial_total_conservation_error']=float(np.max(abs(active.n.sum(axis=0)/active.nTotal-1)));result['SE_wavelength_rows']=len(ctx.spect.wavelength)
 if result['SE_wavelength_rows']!=warm['SE_wavelength_rows']:raise ValueError('Changed spectral layout')
 result['source_collision_contract']=dict(spline=warm['coupled']['source_spline'],MN=warm['coupled']['source_mn']);print(json.dumps(dict(case=a.case,dry_run=a.dry_run,SE_rows=result['SE_wavelength_rows'],load=os.getloadavg())),flush=True)
 if a.dry_run:
  result.update(dry_run=True,telemetry=dict(seconds=time.perf_counter()-started,process_peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss));(a.output/'report.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');return
 original=dict(controls=controls,experimental_collision_convention='barklem2018_MN_LTE_pair',native_spline_control=warm['coupled']['source_spline'],author_charge_control=warm['coupled']['source_mn'])
 try:final,coupled=solve_coupled(ctx,atom,original,a.output,started)
 except Exception as error:
  failed=dict(result,completed=False,SE_balance_claimed=False,initialization_only=True,error_type=type(error).__name__,error=str(error),elapsed_seconds=time.perf_counter()-started)
  current=active.n.copy();conservation=float(np.max(abs(current.sum(axis=0)/active.nTotal-1)))
  if np.isfinite(current).all()and np.all(current>0)and conservation<1e-10:
   np.savez(a.output/'failure_seed.npz',populations=current,nstar=active.nStar.copy());failed['failure_seed']=fingerprint(a.output/'failure_seed.npz')
  failed['conservation']=conservation if np.isfinite(conservation)else None
  (a.output/'failure.json').write_text(json.dumps(failed,indent=2,allow_nan=False)+'\n')
  raise
 np.testing.assert_array_equal(active.nStar,nstar)
 groups=collision_classes(active.atomicModel,ctx.atmos.pyAtmos,ctx.kwargs['eqPops']);collision=sum(groups.values());active.compute_collisions(fillDiagonal=False);np.testing.assert_allclose(collision,active.C,rtol=256*np.finfo(float).eps,atol=0)
 radiation=np.zeros_like(collision)
 for line in active.trans:radiation[line.j,line.i]+=line.Rij;radiation[line.i,line.j]+=line.Rji
 net,gross=process_flow(collision+radiation,final);residual=float(np.max(abs(net)/np.maximum(gross,np.finfo(float).tiny)));conservation=float(np.max(abs(final.sum(axis=0)/active.nTotal-1)))
 if residual>=controls['rtol']or conservation>1e-10:raise ValueError('Independent final physical balance failed')
 np.savez_compressed(a.output/'independent_rates.npz',collision_rates_s_inv=collision,radiative_rates_s_inv=radiation,net_flow_m3_s=net,gross_flow_m3_s=gross)
 continuum,continuum_audits=construct(fixed=final,continuum_only=True);cont_rad=converge_cap_radiation(continuum,rtol=1e-8,maxiter=1000,deadline=started+controls['seconds_per_case']);mu,weights=flux_quadrature(controls['flux_rays']);_,I=emergent_rays(ctx,vac,mu,weights);flux=integrate_flux(I,mu,weights);_,Ic=emergent_rays(continuum,vac,mu,weights);cont=integrate_flux(Ic,mu,weights)
 for context in(ctx,continuum):np.testing.assert_array_equal(context.kwargs['eqPops'].atomicPops[PeriodicTable['O']].n,final)
 if time.perf_counter()>started+controls['seconds_per_case']:raise TimeoutError('Matched-column whole-case budget exceeded')
 if any(not np.isfinite(x).all()or np.any(x<=0)for x in [flux,cont]):raise ValueError('Invalid emergent/continuum flux')
 np.savetxt(a.output/'flux.csv',np.column_stack((wave,flux,cont,flux/cont)),delimiter=',',header='wavelength_air_A,Fnu_W_m2_Hz,continuum_Fnu_W_m2_Hz,normalized_flux',comments='')
 for path,receipt in checks:
  if fingerprint(path)!=receipt:raise ValueError('Protected source changed during paired solve')
 inventory=sorted(set(key for audit in audits+continuum_audits for _,key in audit.requests))
 if 'CN'in inventory:raise ValueError('Unexpected CN baseline')
 result.update(coupled=coupled,independent=dict(physical_rate_residual=residual,particle_conservation_error=conservation,collision_classes=list(groups),rates=fingerprint(a.output/'independent_rates.npz'),kernel=fingerprint(ROOT/'pipeline/bride_rate_audit.py')),continuum_radiation=cont_rad,flux=fingerprint(a.output/'flux.csv'),baseline_lookup_inventory=inventory,telemetry=dict(seconds=time.perf_counter()-started,process_peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,threads=1,load=os.getloadavg()),limitations=['Explicit model-H/source-owned column experiment; no exported native-network packet transplanted; same initial fractions across four fresh cases.','Same experimental collision spline/MN and finite-cap hypotheses, not independent physical model/reference admission.','No added Fe/CN line opacity, validated interpolation/column mass, abundance fit, EOS/charge/thermal feedback or preferred source/model admission.'])
 (a.output/'report.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(json.dumps(dict(case=a.case,iterations=coupled['iterations'],physical_rate_residual=residual,conservation=conservation,telemetry=result['telemetry'])),flush=True)

if __name__=='__main__':main()
