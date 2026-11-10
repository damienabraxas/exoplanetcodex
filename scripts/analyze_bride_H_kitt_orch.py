"""Pinned ORCH continuum replay and Kitt Peak residuals for paired-H O controls."""
import argparse,csv,importlib.util,json,os,subprocess,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def validate_processed_flux(wave,flux,centres):
 """Retain iSpec zero padding, but prohibit it in any analysed +/-5A window."""
 wave=np.asarray(wave);flux=np.asarray(flux);support=np.zeros(len(wave),bool)
 for centre in centres:support|=abs(wave-centre)<=5
 if flux.shape!=wave.shape or not np.isfinite(flux).all()or np.any(flux<0)or np.any(flux[support]<=0):raise ValueError('Invalid broadened flux in comparison support')
 return np.flatnonzero(flux==0).tolist()


def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
 from config.constants import require_codex_path
 from pipeline.bride_reference import fingerprint
 from pipeline.kp_atlas_integrity import read_kurucz2005_residual
 from pipeline.reflected_solar_rv import SOLAR_GRAVITATIONAL_REDSHIFT
 from pipeline.telluric_policy import gate_holding
 from pipeline.prenormalised_guard import assert_not_renormalising
 from pipeline.bride_observables import window_integral
 from pipeline.bride_continuum_audit import selection,trace_fit
 plan_path=ROOT/'config/experiments/rya1281_H_kitt_orch.json';plan=json.loads(plan_path.read_text());orch_plan=ROOT/plan['orch_plan'];pin=json.loads(orch_plan.read_text())['orch_commit'];orch=require_codex_path('experiments.bride_current_orch_checkout');checks=[]
 def check(path,receipt):
  if fingerprint(path)!=receipt:raise ValueError('Changed source '+str(path))
  checks.append((path,receipt))
 if subprocess.check_output(['git','rev-parse','HEAD'],cwd=orch,text=True).strip()!=pin:raise ValueError('Changed ORCH pin')
 native=orch/'pipeline/local_continuum.py'
 if native.read_bytes()!=subprocess.check_output(['git','show',pin+':pipeline/local_continuum.py'],cwd=orch):raise ValueError('Dirty pinned continuum')
 check(native,fingerprint(native));os.environ['CODEX_CONT_REF']='iag';spec=importlib.util.spec_from_file_location('H_kitt_pinned_continuum',native);lc=importlib.util.module_from_spec(spec);sys.modules[spec.name]=lc;spec.loader.exec_module(lc)
 if(lc.MAX_SHIFT_MODEL,lc.MAX_SHIFT,lc.MODEL_Q,lc.ENV_PERCENTILE)!=(.15,.03,90.,95.):raise ValueError('Unexpected ORCH settings')
 comparison_path=ROOT/plan['comparison'];comparison=json.loads(comparison_path.read_text());check(comparison_path,fingerprint(comparison_path));b=comparison['broadening'];observation=require_codex_path('data.kurucz2005_residual');check(observation,comparison['observation'])
 holding='solar_kpno_kurucz2005_corrected';allowed,reason=gate_holding(holding,'kpno_solar_atlas')
 if not allowed:raise ValueError(reason)
 assert_not_renormalising(holding,pre_normalised=True,local_window_envelope=True)
 profile_path=ROOT/plan['iag_receipt'];profile=json.loads(profile_path.read_text());iag=require_codex_path('data.solar_iag_baker2020');check(iag,profile['iag']);check(profile_path,fingerprint(profile_path))
 sys.path.insert(0,str(require_codex_path('engines.ispec')));import ispec
 import inspect
 effects=Path(inspect.getfile(ispec.apply_post_fundamental_effects));check(effects,comparison['ispec_effects'])
 raw={};reports={};source_receipts=[];wave=None
 for item in plan['sources']:
  folder=ROOT/item['path'];r=json.loads((folder/'report.json').read_text());check(folder/'report.json',fingerprint(folder/'report.json'));check(folder/'flux.csv',r['flux'])
  if r['case']!=item['case']:raise ValueError('Wrong named source')
  if item['kind']=='coupled':
   if r['independent']['physical_rate_residual']>=r['controls']['rtol']or r['independent']['particle_conservation_error']>1e-10:raise ValueError('Unadmitted coupled parent')
   check(folder/'populations.npz',r['coupled']['populations'])
  elif item['kind']=='fixed':
   if r['SE_updated']or r['SE_balance_claimed']or r['conservation']>1e-10:raise ValueError('Wrong fixed-population contract')
   for key in ['radiation','continuum_radiation']:
    if max(r[key]['final_J_change'],r[key]['final_boundary_change'])>=1e-8:raise ValueError('Unadmitted fixed radiation')
   check(folder/'populations.npz',r['populations'])
  else:raise ValueError('Unknown source kind')
  data=np.genfromtxt(folder/'flux.csv',delimiter=',',names=True)
  if not all(np.isfinite(data[k]).all()for k in data.dtype.names)or np.any(data['normalized_flux']<=0):raise ValueError('Invalid source flux')
  if wave is None:wave=data['wavelength_air_A'].copy()
  else:np.testing.assert_array_equal(wave,data['wavelength_air_A'])
  raw[item['label']]=data;reports[item['label']]=r;source_receipts.append(dict(**item,report=fingerprint(folder/'report.json'),flux=r['flux']))
 baseline=plan['baseline'];base=reports[baseline]
 for item in plan['sources']:
  r=reports[item['label']]
  for key in ['source_model','atom','controls','cap_length_m','context_helper','backend','LW_constants']:
   if r[key]!=base[key]:raise ValueError('Different shared physical input '+key)
  matched=reports['coupled_'+item['case']]
  if r['state']!=matched['state']:raise ValueError('Different named physical state')
 for c in plan['centres_air_A']:
  if wave[0]>c-5 or wave[-1]<c+5:raise ValueError('Incomplete model continuum coverage')
 w,obs=read_kurucz2005_residual(observation,7765,7782);keep=(w>=wave[0])&(w<=wave[-1]);w=w[keep];obs=obs[keep];original_obs=obs.copy()
 if not np.all(np.diff(w)>0)or not np.isfinite(obs).all()or np.any(obs<=0)or w[0]>min(plan['centres_air_A'])-5 or w[-1]<max(plan['centres_air_A'])+5:raise ValueError('Invalid/short observed AIR coverage')
 segment=np.array([(wave[0]/10,wave[-1]/10)],dtype=[('wave_base',float),('wave_top',float)]);models={};padding=[];arrays=dict(observed_wave_air_A=w,observed_shipped_flux=obs,model_wave_air_A=wave);rv_values=dict(no_gravity=0.,plus_gravity=SOLAR_GRAVITATIONAL_REDSHIFT)
 for rvlabel in plan['RVs']:
  for label,data in raw.items():
   flux=ispec.apply_post_fundamental_effects(wave/10,data['normalized_flux'].copy(),segment,macroturbulence=b['vmac_km_s'],vsini=b['vsini_km_s'],limb_darkening_coeff=b['limb_darkening'],R=b['R'],vrad=(rv_values[rvlabel],))
   zero_indices=validate_processed_flux(wave,flux,plan['centres_air_A']);padding.append(dict(rv=rvlabel,label=label,zero_indices=zero_indices,zero_wavelengths_air_A=wave[zero_indices].tolist()))
   models[rvlabel,label]=flux;arrays['model_'+rvlabel+'_'+label]=flux
 fits={};iag_cache={};support=[]
 for rvlabel in plan['RVs']:
  for ci,centre in enumerate(plan['centres_air_A']):
   for reference in plan['references']:
    if reference=='iag':
     if centre not in iag_cache:iag_cache[centre]=lc.solar_atlas_reference(centre-5,centre+5,b['R'])
     rw,rf=iag_cache[centre];reader=sys.modules['measure_band_ew']
     if reader.IAG_FITS.resolve()!=iag.resolve():raise ValueError('Wrong IAG reader holding')
    for label in raw:
     if reference=='synthesis':rw,rf=wave,models[rvlabel,label]
     if rw[0]>centre-5 or rw[-1]<centre+5:raise ValueError('Truncated continuum reference')
     ref=np.interp(w,rw,rf);rec,cont=lc.fit_model_guided(w,obs,centre,w,ref,exclude_half_width_A=plan['half_width_A'],reference=reference)
     record={k:(None if isinstance(v,float)and not np.isfinite(v)else v)for k,v in rec.__dict__.items()}
     if cont is not None and(not np.isfinite(cont).all()or np.any(cont<=0)):raise ValueError('Invalid ORCH continuum')
     fits[rvlabel,ci,reference,label]=(record,cont)
     if cont is not None:arrays[f'continuum_{rvlabel}_{ci}_{reference}_{label}']=cont
     arrays[f'reference_{rvlabel}_{ci}_{reference}_{label}']=ref
     mask,threshold=selection(w,obs,centre,ref);trace,_=trace_fit(w,obs,centre,ref,mask)
     # Legacy helper's bound is NOT used for application; only identical fit/support arithmetic.
     bound=512*np.finfo(float).eps*trace['geometry']['design_condition_2']
     if abs(trace['level']-rec.level_at_centre)>bound or abs(trace['slope_per_A']-rec.slope_per_A)>bound:raise ValueError('Support trace differs from pinned coefficients')
     support.append(dict(label=label,rv=rvlabel,centre_air_A=centre,reference=reference,selector_threshold=threshold,selected_indices=np.flatnonzero(mask).tolist(),pinned_record=record,trace_scope='Coefficient/support arithmetic only; legacy helper application bound ignored',**{k:trace[k]for k in ['level','slope_per_A','geometry','fit_indices','surviving_indices','coefficient_support_count','postfit_survivor_count','clipping_stop_reason','clipping_limit_reached']}))
 for rvlabel in plan['RVs']:
  for ci in range(len(plan['centres_air_A'])):
   first=fits[rvlabel,ci,'iag',baseline]
   for label in raw:
    rec,cont=fits[rvlabel,ci,'iag',label]
    if rec!=first[0]or(cont is None)!=(first[1]is None):raise ValueError('IAG continuum changed with model label')
    if cont is not None:np.testing.assert_array_equal(cont,first[1])
 rows=[];records=[]
 for rvlabel in plan['RVs']:
  for ci,centre in enumerate(plan['centres_air_A']):
   lo=centre-plan['half_width_A'];hi=centre+plan['half_width_A'];mask=(w>=lo)&(w<=hi);area=lambda flux:1000*window_integral(w,flux,lo,hi)
   for reference in plan['references']:
    common=fits[rvlabel,ci,reference,baseline][1]
    for label in raw:
     rec,cont=fits[rvlabel,ci,reference,label];model=models[rvlabel,label];mf=np.interp(w,wave,model);model_area=1000*window_integral(wave,model,lo,hi)
     row=dict(label=label,rv=rvlabel,reference=reference,centre_air_A=centre,n_observed_window_pixels=int(mask.sum()),model_area_mA=model_area,shipped_observed_area_mA=area(obs),continuum_applied=rec['applied'],continuum_level=rec['level_at_centre'],continuum_slope_per_A=rec['slope_per_A'],continuum_n_bins=rec['n_bins'],continuum_reason=rec['reason'])
     for method,divisor in [('shipped',np.ones_like(obs)),('common',common),('own',cont)]:
      if divisor is None:
       for key in ['observed_area_mA','area_residual_mA','mean_residual','rms_residual']:row[method+'_'+key]=None
      else:
       conditioned=obs/divisor;residual=(mf-conditioned)[mask];oa=area(conditioned);row.update({method+'_observed_area_mA':oa,method+'_area_residual_mA':model_area-oa,method+'_mean_residual':float(np.mean(residual)),method+'_rms_residual':float(np.sqrt(np.mean(residual**2)))})
     rows.append(row);records.append(dict(label=label,rv=rvlabel,reference=reference,centre_air_A=centre,continuum=rec))
 contrasts=[]
 for rvlabel in plan['RVs']:
  for ci,centre in enumerate(plan['centres_air_A']):
   for reference in plan['references']:
    by={r['label']:r for r in rows if r['rv']==rvlabel and r['reference']==reference and r['centre_air_A']==centre}
    for anchor,prefix in [('coupled','coupled_'),('model_anchor','model_anchor_'),('native_anchor','native_anchor_')]:
     A=by[prefix+'model_mass_H'];B=by[prefix+'native_network_H'];contrasts.append(dict(anchor=anchor,rv=rvlabel,reference=reference,centre_air_A=centre,model_area_contrast_mA=B['model_area_mA']-A['model_area_mA'],own_observation_area_contrast_mA=None if A['own_observed_area_mA']is None or B['own_observed_area_mA']is None else B['own_observed_area_mA']-A['own_observed_area_mA'],common_area_residual_contrast_mA=None if A['common_area_residual_mA']is None or B['common_area_residual_mA']is None else B['common_area_residual_mA']-A['common_area_residual_mA'],own_area_residual_contrast_mA=None if A['own_area_residual_mA']is None or B['own_area_residual_mA']is None else B['own_area_residual_mA']-A['own_area_residual_mA']))
 np.testing.assert_array_equal(obs,original_obs)
 for path,receipt in checks:
  if fingerprint(path)!=receipt:raise ValueError('Source changed during replay')
 a.output.mkdir(parents=True,exist_ok=False)
 for name,data in [('summary.csv',rows),('contrasts.csv',contrasts)]:
  with(a.output/name).open('w',newline='')as stream:writer=csv.DictWriter(stream,fieldnames=list(data[0]),lineterminator='\n');writer.writeheader();writer.writerows(data)
 np.savez(a.output/'profiles_and_continua.npz',**arrays)
 (a.output/'support.json').write_text(json.dumps(support,indent=2,allow_nan=False)+'\n')
 import matplotlib
 matplotlib.use('Agg');import matplotlib.pyplot as plt
 fig,axes=plt.subplots(3,3,figsize=(15,11));rv=plan['plot_rv'];ref=plan['plot_reference']
 for ci,centre in enumerate(plan['centres_air_A']):
  mask=abs(w-centre)<=.65;mk=abs(wave-centre)<=.65;cont=fits[rv,ci,ref,baseline][1]
  if cont is not None:axes[0,ci].plot(w[mask],(obs/cont)[mask],color='black',label='Kitt Peak / fixed baseline continuum')
  else:axes[0,ci].text(.1,.5,'Continuum refused',transform=axes[0,ci].transAxes)
  for label,style in [('coupled_model_mass_H','-'),('coupled_native_network_H','--')]:
   model=models[rv,label];axes[0,ci].plot(wave[mk],model[mk],style,label=label)
   if cont is not None:axes[1,ci].plot(w[mask],(np.interp(w,wave,model)-obs/cont)[mask],style,label=label)
  for anchor in ['coupled','model_anchor','native_anchor']:
   delta=models[rv,anchor+'_native_network_H']-models[rv,anchor+'_model_mass_H'];axes[2,ci].plot(wave[mk],delta[mk],label=anchor)
  axes[0,ci].set_title(f'{centre:.3f} AIR A');axes[0,ci].set_ylabel('Normalized flux');axes[1,ci].set_ylabel('Model minus conditioned observation');axes[2,ci].set_ylabel('Native-H minus model-H');
  for ax in axes[:,ci]:ax.set_xlabel('AIR A');ax.legend(fontsize=7)
 fig.suptitle('Bare-background O diagnostics: fixed broadening, plus-gravity, pinned IAG continuum\nNo abundance/continuum/RV tuning; no Fe/CN replacement');fig.tight_layout();fig.savefig(a.output/'kitt_comparison.png',dpi=140);plt.close(fig)
 result=dict(schema='codex.bride_H_kitt_orch_analysis/1',scientific_admission='HOLD',scope=plan['scope'],generator=fingerprint(__file__),plan=fingerprint(plan_path),orch_plan=fingerprint(orch_plan),orch_commit=pin,current_native=fingerprint(native),comparison=fingerprint(comparison_path),observation=fingerprint(observation),holding=holding,policy_reason=reason,frame='Delivered observed vacuum nm -> canonical AIR Angstrom once by established reader; model AIR sampled output is already AIR. No further observation wavelength/telluric correction.',iag=fingerprint(iag),iag_reader=fingerprint(reader.__file__),ispec_effects=fingerprint(effects),broadening=b,RV_km_s=rv_values,sources=source_receipts,continuum_records=records,rows=len(rows),contrasts=len(contrasts),applied=sum(r['continuum']['applied']for r in records),refused=sum(not r['continuum']['applied']for r in records),common_continuum_baseline=baseline,model_dependent_continuum_is_diagnostic=True,statistic='Unweighted sampled residuals/finite AIR areas; no noise model, chi-square, abundance fit or significance claim.',limitations=['Bare-background O column has no added Fe/CN opacity; not a replacement product or direct blanketed-baseline comparison.','Only VIS triplet output/continuum windows, not all-band validation; frozen populations are not SE solutions.','Pinned native local-continuum replay and existing broadening only, not a full ORCH pipeline execution.','Both RV controls retained, no optimization; no preferred H/EOS/abundance or observed normalization adoption.','Common and per-model continuum paths retained separately; input observation stays byte-protected.'])
 result['preserved_ispec_zero_padding_outside_comparison_support']=padding
 result['support_trace_helper']=fingerprint(ROOT/'pipeline/bride_continuum_audit.py')
 (a.output/'report.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(json.dumps({k:result[k]for k in ['rows','contrasts','applied','refused','orch_commit']}))

if __name__=='__main__':main()
