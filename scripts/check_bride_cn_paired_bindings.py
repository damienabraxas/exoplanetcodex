"""Paired native-column H controls: background/LTE arrays, no RT or provider."""
import argparse,csv,json,os,platform,resource,sys,time
from dataclasses import asdict
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();before=time.perf_counter()
 import lightweaver as lw
 import lightweaver.constants as c
 import lightweaver.LwCompiled as compiled
 from lightweaver import rh_atoms
 from lightweaver.atomic_table import PeriodicTable
 from config.constants import require_codex_path
 from pipeline.bride_reference import fingerprint,read_marcs_structure,PASSIVE_ATOM_FACTORIES
 from pipeline.multi_atom import read_multi_atom,to_lightweaver
 from pipeline.bride_binding_controls import control_state
 from pipeline.bride_cn_provider import EqLookupAudit
 from pipeline.bride_cn_diagnostic import frequency_offsets,external_state_extinction
 from pipeline.wavelength_util import air_to_vac
 from scripts.analyze_bride_cn_binding_closure import read_binding_table
 if platform.system()!='Darwin':raise ValueError('Mac paired diagnostic required')
 plan_path=ROOT/'config/experiments/rya1281_cn_paired_bindings.json';plan=json.loads(plan_path.read_text());source=require_codex_path('experiments.bride_cn_native_binding');archive=ROOT/'data/results/rya1281/cn_native_binding_v1/native';checks=[]
 def check(path,receipt):
  if fingerprint(path)!=receipt:raise ValueError('Changed bound input '+str(path))
  checks.append((path,receipt))
 for name in ['report.json','state_binding.csv','composition.json']:
  if(source/name).read_bytes()!=(archive/name).read_bytes():raise ValueError('Changed native inventory '+name)
  check(source/name,fingerprint(archive/name))
 native=json.loads((source/'report.json').read_text());records=read_binding_table(source/'state_binding.csv');model_path=require_codex_path('grids.marcs_standard')/native['model']['name'];check(model_path,native['model']);structure=read_marcs_structure(model_path)
 atom_path=require_codex_path('experiments.bride_oxygen_atom');atom=read_multi_atom(atom_path);check(atom_path,fingerprint(atom_path))
 if atom.abundance!=plan['oxygen_abundance_dex']or len(atom.levels)!=41:raise ValueError('Changed O benchmark input')
 cn_source=require_codex_path('experiments.bride_cn_coefficients')/'CN/states.csv';cn_archive=ROOT/'data/results/rya1281/cn_coefficients_v1/CN/states.csv';check(cn_source,fingerprint(cn_archive));cn=np.genfromtxt(cn_source,delimiter=',',names=True)
 if len(cn)!=168:raise ValueError('Changed sampled CN state count')
 for key in ['temperature_K','electron_pressure_dyn_cm2','gas_pressure_dyn_cm2','rho_g_cm3']:
  if not np.array_equal(cn[key].reshape(3,56),np.tile(records[key],(3,1))):raise ValueError('CN/native binding mismatch '+key)
 lo,hi=plan['range_air_A'];air=np.linspace(lo,hi,int(round((hi-lo)/plan['step_air_A']))+1);wave=air_to_vac(air)/10;sample=np.array([int(np.argmin(abs(air-x)))for x in plan['sample_air_A']])
 cn_chi=np.zeros((len(wave),56))
 for block in cn.reshape(3,56):
  centre=air_to_vac(block['wavelength_air_A'])*1e-8;ratio=block['wavelength_cm']/centre;width=block['doppler_width_Hz']*ratio;damping=block['voigt_a']/ratio;stim=-np.expm1(-block['h_erg_s']*block['c_cm_s']/(block['boltz_erg_K']*centre*block['temperature_K']))
  offsets=frequency_offsets(wave*10,block['wavelength_air_A'],width,block['c_cm_s'],block['wavelength_cm'],input_medium='vacuum',mode='exact_vacuum_frequency')
  population=block['lower_population_per_g']*block['rho_g_cm3']*1e6;strength=np.pi*block['e_esu']**2/(block['electron_mass_g']*block['c_cm_s'])*1e-4
  cn_chi+=external_state_extinction(offsets,damping,width,population,block['gf'],stim,strength)
 a.output.mkdir(parents=True,exist_ok=False);cases={};report=dict(schema='codex.bride_cn_paired_bindings/1',scientific_admission='HOLD',scope=plan['scope'],generator=fingerprint(__file__),plan=fingerprint(plan_path),state_helper=fingerprint(ROOT/'pipeline/bride_binding_controls.py'),CN_kernel=fingerprint(ROOT/'pipeline/bride_cn_diagnostic.py'),native_inventory=fingerprint(source/'report.json'),model=native['model'],atom=fingerprint(atom_path),CN_states=fingerprint(cn_source),backend=fingerprint(compiled.__file__),LW_constants=fingerprint(c.__file__),cases=[])
 for name in plan['cases']:
  state=control_state(structure,records,name,kBoltzmann_SI=c.KBoltzmann);directory=a.output/name;directory.mkdir();packet=asdict(state);packet.update(binding_kind='explicit_normalization_control_not_provider_admission',native_inventory=report['native_inventory'],model=report['model'],CN_population_owner=report['CN_states'],rho_metadata_owner='retained native atmosphere and separate native EOS; neither rescaled to force nH')
  (directory/'input.json').write_text(json.dumps(packet,indent=2,allow_nan=False)+'\n');started=time.perf_counter()
  atmos=lw.Atmosphere.make_1d(scale=lw.ScaleType.Geometric,depthScale=np.array(state.height_m),temperature=np.array(state.temperature_K),ne=np.array(state.ne_m3),nHTot=np.array(state.nH_m3),vturb=np.array(state.vturb_m_s),vlos=np.zeros(56),logG=structure['logg_cgs']-np.log10(100.));atmos.quadrature(2)
  abundance={PeriodicTable[z+1]:value for z,value in enumerate(structure['abundances_dex'])};abundance[PeriodicTable['O']]=atom.abundance
  rad=lw.RadiativeSet([getattr(rh_atoms,x)()for x in PASSIVE_ATOM_FACTORIES]+[to_lightweaver(read_multi_atom(atom_path))],abundance=lw.AtomicAbundance(abundance));rad.set_active('O');eq=rad.compute_eq_pops(atmos,mols=lw.MolecularTable([]));audit=[]
  def factory(e,r,w):
   proxy=EqLookupAudit(e);audit.append(proxy);return lw.LwCompiled.FastBackground(proxy,r,w,Nthreads=1)
  spect=rad.compute_wavelength_grid(extraWavelengths=wave).subset_configuration(wave)
  ctx=lw.Context(atmos,spect,eq,Nthreads=1,conserveCharge=False,backgroundProvider=factory)
  if not np.array_equal(ctx.spect.wavelength,wave):raise ValueError('Paired wavelength identity failed')
  if not np.array_equal(atmos.temperature,state.temperature_K)or not np.array_equal(atmos.ne,state.ne_m3)or not np.array_equal(atmos.nHTot,state.nH_m3):raise ValueError('LW changed explicit binding')
  arrays={key:np.array(getattr(ctx.background,key))for key in ['chi','eta','sca']}
  if any(x.shape!=cn_chi.shape or not np.isfinite(x).all()or np.any(x<=0)for x in arrays.values()):raise ValueError('Unexpected invalid/zero VIS background coefficient')
  pops=eq.atomicPops[PeriodicTable['O']].nStar.copy();total=pops.sum(axis=0);expected=np.array(state.nH_m3)*rad.abundance[PeriodicTable['O']];error=float(np.max(abs(total/expected-1)));bound=16*len(atom.levels)*np.finfo(float).eps
  if not np.isfinite(pops).all()or np.any(pops<0)or error>bound:raise ValueError('O LTE particle budget failed')
  np.savez(directory/'background_window.npz',wavelength_vac_nm=wave,**arrays)
  np.savez(directory/'oxygen_LTE.npz',populations_m3=pops,nH_m3=state.nH_m3,ne_m3=state.ne_m3,temperature_K=state.temperature_K)
  queries=sorted(set(key for proxy in audit for _,key in proxy.requests))
  if 'CN'in queries or len(eq.molecularTable.indices):raise ValueError('Unexpected CN/molecule background inventory')
  cases[name]=dict(state=state,arrays=arrays,pops=pops,total=total,abundance=np.array([rad.abundance[z]for z in range(1,93)]))
  report['cases'].append(dict(case=name,seconds=time.perf_counter()-started,input=fingerprint(directory/'input.json'),background_window=fingerprint(directory/'background_window.npz'),oxygen_LTE=fingerprint(directory/'oxygen_LTE.npz'),particle_budget_max_relative_difference=error,particle_budget_prospective_bound=float(bound),base_lookup_inventory=queries,LW_molecular_equilibrium_species=list(eq.molecularTable.indices)))
  print(json.dumps(dict(case=name,shape=list(arrays['chi'].shape),particle_error=error)),flush=True)
 first,second=(cases[x]for x in plan['cases']);Hratio=np.array(second['state'].nH_m3)/np.array(first['state'].nH_m3);Oratio=second['total']/first['total']
 for field in ['height_m','temperature_K','electron_pressure_Pa','ne_m3','vturb_m_s']:
  if getattr(first['state'],field)!=getattr(second['state'],field):raise ValueError('More than H changed in paired state')
 if not np.array_equal(first['abundance'],second['abundance']):raise ValueError('Paired composition changed')
 fractions=[case['pops']/case['total'][None,:]for case in [first,second]];fraction_delta=float(np.max(abs(fractions[1]-fractions[0])));O_ratio_error=float(np.max(abs(Oratio/Hratio-1)))
 if O_ratio_error>16*len(atom.levels)*np.finfo(float).eps:raise ValueError('O/H budget ratio identity failed')
 B=lw.planck(np.array(first['state'].temperature_K)[None,:],wave[:,None]);cn_eta=cn_chi*B;samples=dict(requested_air_A=np.array(plan['sample_air_A']),actual_air_A=air[sample],wavelength_vac_nm=wave[sample],common_CN_chi=cn_chi[sample],common_CN_eta=cn_eta[sample]);rows=[]
 for label,case in cases.items():
  for key,value in case['arrays'].items():samples[label+'_'+key]=value[sample]
  samples[label+'_postprocessed_chi']=(case['arrays']['chi']+cn_chi)[sample]
 for j in range(56):
  row=dict(depth_index=j+1,logtau500=float(structure['logtau_500'][j]),nH_ratio=float(Hratio[j]),oxygen_LTE_total_ratio=float(Oratio[j]),O_fraction_max_absolute_change=float(np.max(abs(fractions[1][:,j]-fractions[0][:,j]))))
  for key in ['chi','eta','sca']:
   ratio=second['arrays'][key][:,j]/first['arrays'][key][:,j];row[key+'_ratio_min']=float(ratio.min());row[key+'_ratio_max']=float(ratio.max());row[key+'_max_relative_departure_from_H_scaling']=float(np.max(abs(ratio/Hratio[j]-1)))
  mixed=(second['arrays']['chi'][:,j]+cn_chi[:,j])/(first['arrays']['chi'][:,j]+cn_chi[:,j]);row.update(postprocessed_chi_ratio_min=float(mixed.min()),postprocessed_chi_ratio_max=float(mixed.max()),max_common_CN_over_model_background=float(np.max(cn_chi[:,j]/first['arrays']['chi'][:,j])))
  rows.append(row)
 np.savez(a.output/'sampled_backgrounds.npz',**samples)
 with(a.output/'summary.csv').open('w',newline='')as f:w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
 import matplotlib
 matplotlib.use('Agg');import matplotlib.pyplot as plt
 x=np.array([row['logtau500']for row in rows]);fig,axes=plt.subplots(2,1,figsize=(10,7),sharex=True)
 axes[0].plot(x,100*(Hratio-1),label='native-network / model-mass H');axes[0].plot(x,100*(Oratio-1),'--',label='O LTE total ratio');axes[0].fill_between(x,100*(np.array([r['sca_ratio_min']for r in rows])-1),100*(np.array([r['sca_ratio_max']for r in rows])-1),alpha=.25,label='background scattering ratio range');axes[0].axhline(0,color='gray',lw=.7);axes[0].set_ylabel('Change (%)');axes[0].legend()
 axes[1].semilogy(x,[r['max_common_CN_over_model_background']for r in rows]);axes[1].set_ylabel('Max sampled CN / baseline extinction');axes[1].set_xlabel('Source log tau500 (not contribution function)');fig.suptitle('Fixed T/Pe, paired H normalization; background/LTE arrays only, no RT');fig.tight_layout();fig.savefig(a.output/'normalization.png',dpi=150);plt.close(fig)
 for path,receipt in checks:
  if fingerprint(path)!=receipt:raise ValueError('Protected input changed during paired controls')
 report.update(common_fields_exact=True,relative_abundance_exact=True,common_CN_population_fixed=True,sampled_CN_lines=3,oxygen_total_ratio_identity_max_relative_difference=O_ratio_error,oxygen_LTE_fraction_max_absolute_difference=fraction_delta,summary={key:[float(min(r[key]for r in rows)),float(max(r[key]for r in rows))]for key in rows[0]if key not in ['depth_index','logtau500']},seconds_main=time.perf_counter()-before,process_peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,limitations=['LTE initialization/background arrays on two deliberately different H budgets, not an O NLTE solve or spectrum.','LW molecular equilibrium explicitly empty; fixed native CN opacity from only three sampled transitions is postprocessed, not installed as provider or chemistry.','Native CN constants/coordinate rescaling retained; LTE emission diagnostic uses explicit LW Planck convention.','No physical normalization chosen, density correction applied, primary provenance or EOS/reference admission.'])
 (a.output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n');print(json.dumps({k:report[k]for k in ['oxygen_total_ratio_identity_max_relative_difference','oxygen_LTE_fraction_max_absolute_difference','summary','seconds_main','process_peak_rss_bytes']}))

if __name__=='__main__':main()
