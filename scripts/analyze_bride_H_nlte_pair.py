"""Paired H-normalization O NLTE diagnostics; no fit or preferred normalization."""
import argparse,csv,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
 from config.constants import require_codex_path
 from pipeline.bride_reference import fingerprint,read_marcs_structure
 from pipeline.bride_observables import window_integral
 from pipeline.bride_rate_audit import process_flow
 plan_path=ROOT/'config/experiments/rya1281_H_nlte_pair_analysis.json';plan=json.loads(plan_path.read_text());reports=[];spectra=[];populations=[];sources=[];residuals=[]
 for case in plan['cases']:
  folder=a.source/case;report=json.loads((folder/'report.json').read_text());initial=json.loads((folder/'input.json').read_text())
  if report['case']!=case or report.get('dry_run',False)or any(report.get(k)!=v for k,v in initial.items()):raise ValueError('Invalid/mutated solved manifest')
  for path,receipt in [(folder/'flux.csv',report['flux']),(folder/'populations.npz',report['coupled']['populations']),(folder/'history.jsonl',report['coupled']['history']),(folder/'independent_rates.npz',report['independent']['rates'])]:
   if fingerprint(path)!=receipt:raise ValueError('Changed solved artifact '+str(path))
  flux=np.genfromtxt(folder/'flux.csv',delimiter=',',names=True)
  if not all(np.isfinite(flux[k]).all()for k in flux.dtype.names)or np.any(flux['Fnu_W_m2_Hz']<=0)or np.any(flux['continuum_Fnu_W_m2_Hz']<=0):raise ValueError('Invalid solved flux')
  np.testing.assert_array_equal(flux['normalized_flux'],flux['Fnu_W_m2_Hz']/flux['continuum_Fnu_W_m2_Hz'])
  with np.load(folder/'populations.npz')as raw:pop={k:raw[k].copy()for k in raw.files}
  if pop['populations'].shape!=(41,56)or any(not np.isfinite(v).all()or np.any(v<=0)for v in pop.values()):raise ValueError('Invalid O populations')
  with np.load(folder/'independent_rates.npz')as raw:net,gross=process_flow(raw['collision_rates_s_inv']+raw['radiative_rates_s_inv'],pop['populations'])
  residual=float(np.max(abs(net)/np.maximum(gross,np.finfo(float).tiny)))
  if residual>=report['controls']['rtol']or report['independent']['particle_conservation_error']>1e-10:raise ValueError('Independent physical admission failed')
  residuals.append(residual);reports.append(report);spectra.append(flux);populations.append(pop);sources.append(fingerprint(folder/'report.json'))
 A,B=reports;pa,pb=populations;fa,fb=spectra
 for key in ['atom','source_model','warm_report','warm_populations','controls','cap_length_m','background','backend','LW_constants','source_collision_contract']:
  if A[key]!=B[key]:raise ValueError('Changed common NLTE input '+key)
 for key in ['height_m','temperature_K','electron_pressure_Pa','ne_m3','vturb_m_s','ne_owner']:
  if A['state'][key]!=B['state'][key]:raise ValueError('Changed common state '+key)
 np.testing.assert_array_equal(fa['wavelength_air_A'],fb['wavelength_air_A']);Hratio=np.array(B['state']['nH_m3'])/np.array(A['state']['nH_m3'])
 fractions=[p['populations']/p['populations'].sum(axis=0)[None,:]for p in populations];departures=[p['populations']/p['nstar']for p in populations];relative=departures[1]/departures[0]-1
 seed_delta=float(np.max(abs((pb['initial_populations']/pb['initial_populations'].sum(axis=0))/(pa['initial_populations']/pa['initial_populations'].sum(axis=0))-1)))
 if seed_delta>16*41*np.finfo(float).eps:raise ValueError('Different warm fraction seed')
 wave=fa['wavelength_air_A'];delta=fb['normalized_flux']-fa['normalized_flux'];windows=[]
 for centre in plan['centres_air_A']:
  areas=[1000*window_integral(wave,f['normalized_flux'],centre-plan['half_width_A'],centre+plan['half_width_A'])for f in spectra];windows.append(dict(centre_air_A=centre,model_mass_H_mA=areas[0],native_network_H_mA=areas[1],delta_mA=areas[1]-areas[0]))
 model_path=require_codex_path('grids.marcs_standard')/A['source_model']['name']
 if fingerprint(model_path)!=A['source_model']:raise ValueError('Changed source model')
 structure=read_marcs_structure(model_path);a.output.mkdir(parents=True,exist_ok=False)
 with(a.output/'windows.csv').open('w',newline='')as f:w=csv.DictWriter(f,fieldnames=list(windows[0]),lineterminator='\n');w.writeheader();w.writerows(windows)
 level_rows=[dict(level_index=j+1,max_relative_departure_change=float(np.max(abs(relative[j]))),max_absolute_departure_change=float(np.max(abs(departures[1][j]-departures[0][j]))),max_absolute_fraction_change=float(np.max(abs(fractions[1][j]-fractions[0][j]))))for j in range(41)]
 with(a.output/'levels.csv').open('w',newline='')as f:w=csv.DictWriter(f,fieldnames=list(level_rows[0]),lineterminator='\n');w.writeheader();w.writerows(level_rows)
 np.savez(a.output/'departure_comparison.npz',model_mass_H=departures[0],native_network_H=departures[1],relative_change=relative,logtau500=structure['logtau_500'])
 import matplotlib
 matplotlib.use('Agg');import matplotlib.pyplot as plt
 fig,axes=plt.subplots(2,2,figsize=(12,8))
 for ax,centre in zip(axes.flat[:3],plan['centres_air_A']):
  mask=abs(wave-centre)<=.65;ax.plot(wave[mask],fa['normalized_flux'][mask],label='model-mass H');ax.plot(wave[mask],fb['normalized_flux'][mask],'--',label='native-network H');ax.set_title(f'O window {centre:.3f} AIR A');ax.set_ylabel('Intrinsic Fnu / own true continuum');ax.set_xlabel('AIR A');ax.legend()
 axes.flat[3].plot(wave,delta);axes.flat[3].set_xlabel('AIR A');axes.flat[3].set_ylabel('Native-network minus model-mass normalized flux');fig.suptitle('Paired O NLTE normalization diagnostic; no observation/fit');fig.tight_layout();fig.savefig(a.output/'profiles.png',dpi=150);plt.close(fig)
 fig,ax=plt.subplots(figsize=(11,6));bound=float(np.max(abs(relative*100)));image=ax.pcolormesh(structure['logtau_500'],np.arange(1,42),relative*100,cmap='RdBu_r',vmin=-bound,vmax=bound,shading='nearest');fig.colorbar(image,ax=ax,label='Relative departure change (%)');ax.set_xlabel('Source log tau500 (not contribution function)');ax.set_ylabel('Source O level index');ax.set_title('All41levels: native-network / model-mass departure coefficients');fig.tight_layout();fig.savefig(a.output/'departures.png',dpi=150);plt.close(fig)
 result=dict(schema='codex.bride_H_nlte_pair_analysis/1',scientific_admission='HOLD',scope=plan['scope'],generator=fingerprint(__file__),plan=fingerprint(plan_path),source_reports=sources,independent_rate_residuals=residuals,warm_fraction_max_relative_difference=seed_delta,Hratio_range=[float(Hratio.min()),float(Hratio.max())],max_relative_departure_change=float(np.max(abs(relative))),max_absolute_departure_change=float(np.max(abs(departures[1]-departures[0]))),max_absolute_population_fraction_change=float(np.max(abs(fractions[1]-fractions[0]))),max_normalized_flux_difference=float(np.max(abs(delta))),max_relative_continuum_difference=float(np.max(abs(fb['continuum_Fnu_W_m2_Hz']/fa['continuum_Fnu_W_m2_Hz']-1))),windows=windows,limitations=['Both inherit joint1e-6accuracy; paired difference is not a validated abundance error or physical accuracy bound.','Bare native backgrounds and experimental cap/collision hypotheses; no direct earlier Fe-blanketed/profile/observation comparison.','No CN added, chemistry/EOS update, source parameter fit, preferred normalization or solar/model admission.'])
 (a.output/'report.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(json.dumps({k:result[k]for k in ['independent_rate_residuals','max_relative_departure_change','max_normalized_flux_difference','max_relative_continuum_difference','windows']}))

if __name__=='__main__':main()
