"""Fixed-fraction transfer intervention against completed paired O NLTE controls."""
import argparse,json,platform,resource,sys,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--plan',type=Path,default=ROOT/'config/experiments/rya1281_H_fixed_transfer.json');a=p.parse_args()
 from importlib.metadata import version
 import lightweaver.constants as constants
 import lightweaver.LwCompiled as compiled
 from lightweaver.atomic_table import PeriodicTable
 from config.constants import require_codex_path
 from pipeline.bride_binding_controls import BindingControl
 from pipeline.bride_H_nlte_control import make_context,scale_initial_fractions
 from pipeline.bride_reference import fingerprint,read_marcs_structure
 from pipeline.multi_atom import read_multi_atom
 from pipeline.bride_top_cap import converge_cap_radiation,emergent_rays
 from pipeline.bride_observables import flux_quadrature,integrate_flux,window_integral
 from pipeline.wavelength_util import air_to_vac
 if platform.system()!='Darwin'or version('lightweaver')!='0.16.1':raise RuntimeError('Mac/LW0.16.1 required')
 plan_path=a.plan;plan=json.loads(plan_path.read_text());source=ROOT/plan['paired_archive'];reports=[];parents=[];pops=[];checks=[]
 def check(path,receipt):
  if fingerprint(path)!=receipt:raise ValueError('Changed source '+str(path))
  checks.append((path,receipt))
 for case in plan['cases']:
  folder=source/case;r=json.loads((folder/'report.json').read_text())
  if r['case']!=case or r.get('dry_run',False)or r['independent']['physical_rate_residual']>=r['controls']['rtol']or r['independent']['particle_conservation_error']>1e-10:raise ValueError('Unadmitted parent')
  check(folder/'report.json',fingerprint(folder/'report.json'));check(folder/'populations.npz',r['coupled']['populations']);check(folder/'flux.csv',r['flux'])
  with np.load(folder/'populations.npz')as data:n=data['populations'].copy()
  if n.shape!=(41,56)or not np.isfinite(n).all()or np.any(n<=0):raise ValueError('Invalid parent population')
  reports.append(r);pops.append(n);parents.append(np.genfromtxt(folder/'flux.csv',delimiter=',',names=True))
 A,B=reports
 for key in ['atom','source_model','controls','cap_length_m','background','backend','LW_constants','context_helper','state_helper','source_collision_contract']:
  if A[key]!=B[key]:raise ValueError('Changed common parent input '+key)
 for key in ['height_m','temperature_K','electron_pressure_Pa','ne_m3','vturb_m_s','ne_owner']:
  if A['state'][key]!=B['state'][key]:raise ValueError('Changed common parent state '+key)
 check(ROOT/'pipeline/bride_H_nlte_control.py',A['context_helper']);check(ROOT/'pipeline/bride_binding_controls.py',A['state_helper']);check(Path(compiled.__file__),A['backend']);check(Path(constants.__file__),A['LW_constants'])
 model_path=require_codex_path('grids.marcs_standard')/A['source_model']['name'];atom_path=require_codex_path('experiments.bride_oxygen_atom');check(model_path,A['source_model']);check(atom_path,A['atom']);structure=read_marcs_structure(model_path)
 np.testing.assert_array_equal(parents[0]['wavelength_air_A'],parents[1]['wavelength_air_A']);wave=parents[0]['wavelength_air_A'];vac=air_to_vac(wave)/10;extra=np.array([vac[0],vac[-1]]);mu,weights=flux_quadrature(A['controls']['flux_rays'])
 np.testing.assert_allclose(pops[1].sum(axis=0)/pops[0].sum(axis=0),np.array(B['state']['nH_m3'])/np.array(A['state']['nH_m3']),rtol=256*np.finfo(float).eps,atol=0)
 fraction_index=plan['cases'].index(plan['fraction_source'])
 a.output.mkdir(parents=True,exist_ok=False);spectra=[];fixed_fields=[]
 for case,r,target in zip(plan['cases'],reports,pops):
  folder=a.output/case;folder.mkdir();started=time.perf_counter();deadline=started+plan['seconds_per_case'];state=BindingControl(**{k:tuple(v)if isinstance(v,list)else v for k,v in r['state'].items()});fixed=scale_initial_fractions(pops[fraction_index],target.sum(axis=0));fixed_fields.append(fixed)
  manifest=dict(schema='codex.bride_H_fixed_transfer/1',scientific_admission='HOLD',case=case,scope=plan['scope'],generator=fingerprint(__file__),plan=fingerprint(plan_path),parent_report=fingerprint(source/case/'report.json'),fraction_source=fingerprint(source/plan['fraction_source']/'populations.npz'),state=r['state'],controls=r['controls'],cap_length_m=r['cap_length_m'],source_model=r['source_model'],atom=r['atom'],context_helper=r['context_helper'],backend=r['backend'],LW_constants=r['LW_constants'])
  (folder/'input.json').write_text(json.dumps(manifest,indent=2,allow_nan=False)+'\n')
  def construct(continuum_only):return make_context(state,structure,read_multi_atom(atom_path),rays=r['controls']['rays'],qscale=r['controls']['qscale'],extra_vac_nm=extra,cap_length_m=r['cap_length_m'],ng=False,fixed=fixed,continuum_only=continuum_only)
  ctx,audits=construct(False)
  if len(ctx.activeAtoms)!=0 or len(ctx.spect.wavelength)!=r['SE_wavelength_rows']:raise ValueError('Fixed-transfer atom/grid identity failed')
  rad=converge_cap_radiation(ctx,rtol=plan['radiation_rtol'],maxiter=plan['radiation_maxiter'],deadline=deadline);continuum,cont_audits=construct(True);crad=converge_cap_radiation(continuum,rtol=plan['radiation_rtol'],maxiter=plan['radiation_maxiter'],deadline=deadline)
  _,I=emergent_rays(ctx,vac,mu,weights);_,Ic=emergent_rays(continuum,vac,mu,weights);flux=integrate_flux(I,mu,weights);cont=integrate_flux(Ic,mu,weights)
  if time.perf_counter()>deadline:raise TimeoutError('Fixed-transfer case budget exceeded')
  if any(not np.isfinite(x).all()or np.any(x<=0)for x in [flux,cont]):raise ValueError('Invalid fixed-transfer flux')
  for context in [ctx,continuum]:np.testing.assert_array_equal(context.kwargs['eqPops'].atomicPops[PeriodicTable['O']].n,fixed)
  inventory=sorted(set(key for audit in audits+cont_audits for _,key in audit.requests))
  if 'CN'in inventory:raise ValueError('Unexpected CN lookup')
  conservation=float(np.max(abs(fixed.sum(axis=0)/target.sum(axis=0)-1)))
  if conservation>1e-10:raise ValueError('Fixed-total conservation failed')
  np.savez(folder/'populations.npz',fixed_populations=fixed);np.savetxt(folder/'flux.csv',np.column_stack((wave,flux,cont,flux/cont)),delimiter=',',header='wavelength_air_A,Fnu_W_m2_Hz,continuum_Fnu_W_m2_Hz,normalized_flux',comments='')
  for path,receipt in checks:
   if fingerprint(path)!=receipt:raise ValueError('Source changed during intervention')
  manifest.update(radiation=rad,continuum_radiation=crad,SE_updated=False,SE_balance_claimed=False,population_identity='exact before/after full and continuum ray solves',conservation=conservation,baseline_lookup_inventory=inventory,flux=fingerprint(folder/'flux.csv'),populations=fingerprint(folder/'populations.npz'),telemetry=dict(seconds=time.perf_counter()-started,process_peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,threads=1))
  (folder/'report.json').write_text(json.dumps(manifest,indent=2,allow_nan=False)+'\n');spectra.append(dict(flux=flux,continuum=cont,normalized=flux/cont));print(json.dumps(dict(case=case,telemetry=manifest['telemetry'],radiation=rad,conservation=conservation)),flush=True)
 fractions=[n/n.sum(axis=0)for n in fixed_fields];fraction_delta=float(np.max(abs(fractions[1]/fractions[0]-1)))
 if fraction_delta>16*41*np.finfo(float).eps:raise ValueError('Different frozen fractions')
 replay=float(np.max(abs(spectra[fraction_index]['normalized']-parents[fraction_index]['normalized_flux'])));creplay=float(np.max(abs(spectra[fraction_index]['continuum']/parents[fraction_index]['continuum_Fnu_W_m2_Hz']-1)))
 if replay>plan['replay_max_absolute_normalized_flux']or creplay>plan['replay_max_relative_continuum']:raise ValueError('Same-H fixed/coupled replay screen failed')
 coupled=parents[1]['normalized_flux']-parents[0]['normalized_flux'];frozen=spectra[1]['normalized']-spectra[0]['normalized'];feedback=coupled-frozen;rows=[]
 for centre in plan['centres_air_A']:
  area=lambda n:1000*window_integral(wave,n,centre-plan['half_width_A'],centre+plan['half_width_A'])
  c=area(parents[1]['normalized_flux'])-area(parents[0]['normalized_flux']);f=area(spectra[1]['normalized'])-area(spectra[0]['normalized']);rows.append(dict(centre_air_A=centre,coupled_delta_mA=c,fixed_delta_mA=f,conditional_remainder_mA=c-f))
 import csv
 with(a.output/'windows.csv').open('w',newline='')as stream:w=csv.DictWriter(stream,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
 np.savez(a.output/'contrasts.npz',wavelength_air_A=wave,coupled=coupled,fixed=frozen,conditional_remainder=feedback)
 import matplotlib
 matplotlib.use('Agg');import matplotlib.pyplot as plt
 fig,axes=plt.subplots(3,1,figsize=(10,9))
 for ax,centre in zip(axes,plan['centres_air_A']):
  mask=abs(wave-centre)<=.65
  for values,label in [(coupled,'Coupled pair'),(frozen,'Fixed common fractions'),(feedback,'Conditional remainder')]:ax.plot(wave[mask],values[mask],label=label)
  ax.set_xlabel('AIR A');ax.set_ylabel('Native-H minus model-H normalized flux');ax.set_title(f'{centre:.3f} A');ax.legend()
 fig.suptitle('Fixed-fraction transfer intervention; no unique physical decomposition');fig.tight_layout();fig.savefig(a.output/'contrasts.png',dpi=150);plt.close(fig)
 result=dict(schema='codex.bride_H_fixed_transfer_analysis/1',scientific_admission='HOLD',scope=plan['scope'],generator=fingerprint(__file__),plan=fingerprint(plan_path),source_reports=[fingerprint(source/case/'report.json')for case in plan['cases']],case_reports=[fingerprint(a.output/case/'report.json')for case in plan['cases']],fraction_max_relative_difference=fraction_delta,same_H_replay_max_normalized_flux_difference=replay,same_H_replay_max_relative_continuum_difference=creplay,max_absolute_coupled_contrast=float(np.max(abs(coupled))),max_absolute_fixed_contrast=float(np.max(abs(frozen))),max_absolute_conditional_remainder=float(np.max(abs(feedback))),max_relative_fixed_continuum_contrast=float(np.max(abs(spectra[1]['continuum']/spectra[0]['continuum']-1))),windows=rows,limitations=['Frozen O fractions are deliberately not recomputed in SE; no physical rate-balance claim.','Transfer intervention includes H background and rescaled O total at fixed fractions; conditional remainder is not a unique collision contribution.','Replay screen and inherited radiation/joint controls are numerical diagnostics, not independent physical accuracy.','No chemistry/EOS/model/abundance/reference-data or all-band admission.'])
 (a.output/'report.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(json.dumps(result),flush=True)

if __name__=='__main__':main()
