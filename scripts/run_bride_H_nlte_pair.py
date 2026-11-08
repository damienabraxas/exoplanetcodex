"""Bounded paired-H O NLTE column diagnostic with inherited numerical controls."""
import argparse,json,os,platform,resource,sys,time
from dataclasses import asdict
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--case',choices=['model_mass_H','native_network_H'],required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--dry-run',action='store_true');a=p.parse_args()
 import lightweaver as lw
 import lightweaver.constants as c
 import lightweaver.LwCompiled as compiled
 from lightweaver.atomic_table import PeriodicTable
 from importlib.metadata import version
 from config.constants import require_codex_path
 from pipeline.bride_reference import fingerprint,read_marcs_structure
 from pipeline.multi_atom import read_multi_atom
 from pipeline.bride_binding_controls import control_state
 from pipeline.bride_H_nlte_control import make_context,scale_initial_fractions
 from pipeline.bride_coupled_control import solve_coupled
 from pipeline.bride_top_cap import converge_cap_radiation,emergent_rays
 from pipeline.bride_observables import flux_quadrature,integrate_flux
 from pipeline.bride_rate_audit import collision_classes,process_flow
 from pipeline.wavelength_util import air_to_vac
 from scripts.analyze_bride_cn_binding_closure import read_binding_table
 if platform.system()!='Darwin'or version('lightweaver')!='0.16.1':raise RuntimeError('Mac/LW0.16.1 required')
 plan_path=ROOT/'config/experiments/rya1281_H_nlte_pair.json';plan=json.loads(plan_path.read_text());warm_report=ROOT/plan['warm_report'];warm=json.loads(warm_report.read_text());controls=plan['controls']
 if controls!=warm['controls']:raise ValueError('Changed inherited joint solver controls')
 checks=[]
 def check(path,receipt):
  if fingerprint(path)!=receipt:raise ValueError('Changed predeclared input '+str(path))
  checks.append((path,receipt))
 warm_path=ROOT/plan['warm_populations'];check(warm_path,warm['coupled']['populations']);check(ROOT/'pipeline/bride_coupled_control.py',warm['coupled']['helper'])
 with np.load(warm_path)as data:pop=data['populations'].copy()
 paired=ROOT/plan['paired_array_archive'];paired_report=json.loads((paired/'report.json').read_text());check(Path(compiled.__file__),paired_report['backend']);check(Path(c.__file__),paired_report['LW_constants'])
 source=require_codex_path('experiments.bride_cn_native_binding');archive=ROOT/'data/results/rya1281/cn_native_binding_v1/native'
 for name in ['report.json','state_binding.csv','composition.json']:check(source/name,fingerprint(archive/name))
 native=json.loads((source/'report.json').read_text());records=read_binding_table(source/'state_binding.csv');model_path=require_codex_path('grids.marcs_standard')/native['model']['name'];check(model_path,native['model']);structure=read_marcs_structure(model_path);state=control_state(structure,records,a.case,kBoltzmann_SI=c.KBoltzmann)
 if asdict(state)!={k:tuple(v)if k not in ('case','ne_owner','nH_owner')else v for k,v in json.loads((paired/a.case/'input.json').read_text()).items()if k in asdict(state)}:raise ValueError('Changed paired binding identity')
 atom_path=require_codex_path('experiments.bride_oxygen_atom');check(atom_path,warm['atom']);atom=read_multi_atom(atom_path)
 if atom.abundance!=8.70 or len(atom.levels)!=41 or pop.shape!=(41,56):raise ValueError('Changed source atom/warm layout')
 lo,hi=plan['range_air_A'];wave=np.linspace(lo,hi,int(round((hi-lo)/plan['step_air_A']))+1);vac=air_to_vac(wave)/10;extra=np.array([vac[0],vac[-1]])
 cap_length=float(structure['column_mass_kg_m2'][0]/structure['rho_kg_m3'][0]);a.output.mkdir(parents=True,exist_ok=False)
 result=dict(schema='codex.bride_H_nlte_pair/1',scientific_admission='HOLD',case=a.case,scope=plan['scope'],generator=fingerprint(__file__),plan=fingerprint(plan_path),state=asdict(state),state_helper=fingerprint(ROOT/'pipeline/bride_binding_controls.py'),context_helper=fingerprint(ROOT/'pipeline/bride_H_nlte_control.py'),native_inventory=fingerprint(source/'report.json'),paired_report=fingerprint(paired/'report.json'),source_model=native['model'],atom=fingerprint(atom_path),warm_report=fingerprint(warm_report),warm_populations=fingerprint(warm_path),initialization=plan['initialization'],controls=controls,cap_length_m=cap_length,background=plan['background'],backend=fingerprint(compiled.__file__),LW_constants=fingerprint(c.__file__))
 (a.output/'input.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');started=time.perf_counter()
 def construct(*,fixed=None,continuum_only=False):return make_context(state,structure,read_multi_atom(atom_path),rays=controls['rays'],qscale=controls['qscale'],extra_vac_nm=extra,cap_length_m=cap_length,ng=controls['ng'],fixed=fixed,continuum_only=continuum_only)
 ctx,audits=construct();active=ctx.activeAtoms[0];nstar=active.nStar.copy();active.n[:]=scale_initial_fractions(pop,active.nTotal)
 result['initial_total_conservation_error']=float(np.max(abs(active.n.sum(axis=0)/active.nTotal-1)));result['SE_wavelength_rows']=len(ctx.spect.wavelength);result['source_collision_contract']=dict(spline=warm['coupled']['source_spline'],MN=warm['coupled']['source_mn']);print(json.dumps(dict(case=a.case,dry_run=a.dry_run,SE_rows=result['SE_wavelength_rows'],load=os.getloadavg())),flush=True)
 if a.dry_run:
  result.update(dry_run=True,telemetry=dict(seconds=time.perf_counter()-started,process_peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss));(a.output/'report.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');return
 original=dict(controls=controls,experimental_collision_convention='barklem2018_MN_LTE_pair',native_spline_control=warm['coupled']['source_spline'],author_charge_control=warm['coupled']['source_mn'])
 final,coupled=solve_coupled(ctx,atom,original,a.output,started);np.testing.assert_array_equal(active.nStar,nstar)
 groups=collision_classes(active.atomicModel,ctx.atmos.pyAtmos,ctx.kwargs['eqPops']);collision=sum(groups.values());active.compute_collisions(fillDiagonal=False);np.testing.assert_allclose(collision,active.C,rtol=256*np.finfo(float).eps,atol=0)
 radiation=np.zeros_like(collision)
 for line in active.trans:radiation[line.j,line.i]+=line.Rij;radiation[line.i,line.j]+=line.Rji
 net,gross=process_flow(collision+radiation,final);residual=float(np.max(abs(net)/np.maximum(gross,np.finfo(float).tiny)));conservation=float(np.max(abs(final.sum(axis=0)/active.nTotal-1)))
 if residual>=controls['rtol']or conservation>1e-10:raise ValueError('Independent final physical balance failed')
 np.savez_compressed(a.output/'independent_rates.npz',collision_rates_s_inv=collision,radiative_rates_s_inv=radiation,net_flow_m3_s=net,gross_flow_m3_s=gross)
 continuum,continuum_audits=construct(fixed=final,continuum_only=True);cont_rad=converge_cap_radiation(continuum,rtol=1e-8,maxiter=1000,deadline=started+controls['seconds_per_case']);mu,weights=flux_quadrature(controls['flux_rays']);_,I=emergent_rays(ctx,vac,mu,weights);flux=integrate_flux(I,mu,weights);_,Ic=emergent_rays(continuum,vac,mu,weights);cont=integrate_flux(Ic,mu,weights)
 for context in(ctx,continuum):np.testing.assert_array_equal(context.kwargs['eqPops'].atomicPops[PeriodicTable['O']].n,final)
 np.savetxt(a.output/'flux.csv',np.column_stack((wave,flux,cont,flux/cont)),delimiter=',',header='wavelength_air_A,Fnu_W_m2_Hz,continuum_Fnu_W_m2_Hz,normalized_flux',comments='')
 for path,receipt in checks:
  if fingerprint(path)!=receipt:raise ValueError('Protected source changed during paired solve')
 inventory=sorted(set(key for audit in audits+continuum_audits for _,key in audit.requests))
 if 'CN'in inventory:raise ValueError('Unexpected CN baseline')
 result.update(coupled=coupled,independent=dict(physical_rate_residual=residual,particle_conservation_error=conservation,collision_classes=list(groups),rates=fingerprint(a.output/'independent_rates.npz'),kernel=fingerprint(ROOT/'pipeline/bride_rate_audit.py')),continuum_radiation=cont_rad,flux=fingerprint(a.output/'flux.csv'),baseline_lookup_inventory=inventory,telemetry=dict(seconds=time.perf_counter()-started,process_peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,threads=1,load=os.getloadavg()),limitations=['Bare FastBackground diagnostic, different from prior Fe/bound-blanketed spectra; previous fractions used only as initialization.','Same experimental collision spline/MN and finite-cap hypotheses, not independent physical model/reference admission.','No added CN line opacity, abundance fit, EOS/charge/thermal feedback or preferred H normalization.'])
 (a.output/'report.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(json.dumps(dict(case=a.case,iterations=coupled['iterations'],physical_rate_residual=residual,conservation=conservation,telemetry=result['telemetry'])),flush=True)

if __name__=='__main__':main()
