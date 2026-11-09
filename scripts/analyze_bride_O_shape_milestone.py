"""Complete O777 shape/broadening diagnostic package, without inverse fitting."""
import argparse,csv,importlib.util,inspect,json,os,subprocess,sys,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def shape_metrics(wave,flux,centre,*,core_half_A,feature_half_A,shoulder_inner_A,shoulder_outer_A):
 from pipeline.bride_feii_profile import feature_metrics
 from pipeline.bride_observables import window_integral
 out=feature_metrics(wave,flux,centre,core_half_A=core_half_A,shoulder_inner_A=shoulder_inner_A,shoulder_outer_A=shoulder_outer_A)
 out['relative_feature_area_mA']=1000*window_integral(wave,np.asarray(flux)/out['shoulder_level'],centre-feature_half_A,centre+feature_half_A)
 out['shoulder_asymmetry']=(out['right_shoulder']-out['left_shoulder'])/out['shoulder_level']
 ranges=[]
 for sign in [-1,1]:
  lo,hi=sorted([centre+sign*shoulder_inner_A,centre+sign*shoulder_outer_A]);mask=(wave>=lo)&(wave<=hi)
  if mask.sum()<2:raise ValueError('Undersampled shoulder diagnostic')
  ranges.append(float(np.ptp(np.asarray(flux)[mask])/out['shoulder_level']))
 out['maximum_shoulder_fractional_range']=max(ranges)
 return out


def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();started=time.perf_counter()
 from astropy.io import fits
 from config.constants import require_codex_path
 from pipeline.bride_reference import fingerprint
 from pipeline.wavelength_util import vac_to_air
 from pipeline.bride_feii_profile import feature_metrics
 from pipeline.telluric_policy import gate_holding
 from pipeline.prenormalised_guard import assert_not_renormalising
 from scripts.analyze_bride_H_kitt_orch import validate_processed_flux
 plan_path=ROOT/'config/experiments/rya1281_O_shape_milestone.json';plan=json.loads(plan_path.read_text());checks=[]
 def check(path,receipt):
  if fingerprint(path)!=receipt:raise ValueError('Changed milestone input '+str(path))
  checks.append((path,receipt))
 for name,receipt in plan['source_files'].items():check(ROOT/name,receipt)
 parent=json.loads((ROOT/'data/results/rya1281/H_kitt_orch_v1/report.json').read_text());labels=[x['label']for x in parent['sources']]
 with np.load(ROOT/'data/results/rya1281/H_kitt_orch_v1/profiles_and_continua.npz')as data:archived={k:data[k].copy()for k in data.files}
 atlas=require_codex_path('data.solar_iag_baker2020');kitt=require_codex_path('data.kurucz2005_residual');check(atlas,parent['iag']);check(kitt,parent['observation']);policies={}
 for holding,instrument in [('solar_kpno_kurucz2005_corrected','kpno_solar_atlas'),('solar_iag','iag_fts_solar_atlas')]:
  allowed,reason=gate_holding(holding,instrument)
  if not allowed:raise ValueError(reason)
  policies[holding]=reason
 assert_not_renormalising('solar_kpno_kurucz2005_corrected',pre_normalised=True,local_window_envelope=True)
 orch=require_codex_path('experiments.bride_current_orch_checkout');native=orch/'pipeline/local_continuum.py';pin=parent['orch_commit'];check(native,parent['current_native'])
 if subprocess.check_output(['git','rev-parse','HEAD'],cwd=orch,text=True).strip()!=pin or native.read_bytes()!=subprocess.check_output(['git','show',pin+':pipeline/local_continuum.py'],cwd=orch):raise ValueError('Changed pinned ORCH')
 os.environ['CODEX_CONT_REF']='iag';spec=importlib.util.spec_from_file_location('O_shape_pinned_continuum',native);lc=importlib.util.module_from_spec(spec);sys.modules[spec.name]=lc;spec.loader.exec_module(lc)
 sys.path.insert(0,str(require_codex_path('engines.ispec')));import ispec
 check(Path(inspect.getfile(ispec.apply_post_fundamental_effects)),parent['ispec_effects'])
 with fits.open(atlas)as h:
  iw=vac_to_air(1e8/np.asarray(h[1].data['v'],float));iff=np.asarray(h[1].data['s'],float).copy();flags=np.asarray(h[1].data['flags'],float).copy();order=np.argsort(iw);iw,iff,flags=iw[order],iff[order],flags[order]
 kw=archived['observed_wave_air_A'];kf=archived['observed_shipped_flux'];mw=archived['model_wave_air_A'];observations={};flag_records=[]
 for ci,centre in enumerate(plan['centres_air_A']):
  gw,gf=lc.solar_atlas_reference(centre-5,centre+5,parent['broadening']['R']);reader=sys.modules['measure_band_ew']
  if reader.IAG_FITS.resolve()!=atlas.resolve():raise ValueError('Wrong IAG holding')
  np.testing.assert_array_equal(np.interp(kw,gw,gf),archived[f'reference_plus_gravity_{ci}_iag_coupled_model_mass_H'])
  for label,wave,flux in [('kitt_shipped',kw,kf),('iag_native',iw,iff),('iag_pinned_gaussian',gw,gf)]:
   keep=abs(wave-centre)<=1.1;sw=np.array(wave[keep]);sf=np.array(flux[keep])
   if sw[0]>centre-1.0 or sw[-1]<centre+1.0 or not np.isfinite(sf).all()or np.any(sf<=0):raise ValueError('Invalid/truncated shape observation')
   observations[ci,label]=(sw,sf)
  keep=abs(iw-centre)<=1.1
  if not keep.any()or not np.isfinite(flags[keep]).all()or np.any(flags[keep]!=0):raise ValueError('IAG nonzero/missing feature flags')
  flag_records.append(dict(centre_air_A=centre,half_width_A=1.1,count=int(keep.sum()),all_exact_zero=True))
 raw={}
 for source in parent['sources']:
  folder=ROOT/source['path'];check(folder/'report.json',source['report']);check(folder/'flux.csv',source['flux']);data=np.genfromtxt(folder/'flux.csv',delimiter=',',names=True);np.testing.assert_array_equal(data['wavelength_air_A'],mw);raw[source['label']]=data['normalized_flux'].copy()
 b=parent['broadening'];segment=np.array([(mw[0]/10,mw[-1]/10)],dtype=[('wave_base',float),('wave_top',float)]);models={};padding=[];variant_parameters={}
 for variant in plan['broadening_variants']:
  params={k:variant.get(k,b[k])for k in ['vmac_km_s','vsini_km_s','R','limb_darkening']};name=variant['name'];variant_parameters[name]=params
  for rv in plan['RVs']:
   for label in labels:
    flux=ispec.apply_post_fundamental_effects(mw/10,raw[label].copy(),segment,macroturbulence=params['vmac_km_s'],vsini=params['vsini_km_s'],limb_darkening_coeff=params['limb_darkening'],R=params['R'],vrad=(parent['RV_km_s'][rv],));zeros=validate_processed_flux(mw,flux,plan['centres_air_A'])
    if name=='nominal':np.testing.assert_array_equal(flux,archived['model_'+rv+'_'+label])
    models[name,rv,label]=flux;padding.append(dict(variant=name,rv=rv,label=label,zero_padding_indices=zeros))
  print(json.dumps(dict(broadening_variant=name,elapsed_seconds=time.perf_counter()-started)),flush=True)
 source_metrics=[];comparisons=[];inventory=[];profiles={};metrics={};canary_contrast=[];canary_area=[];grid_deltas=[]
 intake=json.loads((ROOT/'data/results/rya1281/anchor_influence_v1/candidate_intake.json').read_text())['rows']
 def measure(wave,flux,centre,scheme):
  args=dict(core_half_A=plan['core_half_A'],feature_half_A=plan['feature_half_A'],**plan['shoulders'][scheme]);m=shape_metrics(wave,flux,centre,**args)
  for scale in plan['scale_canaries']:
   other=shape_metrics(wave,flux*scale,centre,**args);dc=abs(other['contrast']-m['contrast']);da=abs(other['relative_feature_area_mA']-m['relative_feature_area_mA']);canary_contrast.append(dc);canary_area.append(da)
   if dc>plan['scale_contrast_atol']or da>plan['scale_area_atol_mA']:raise ValueError('Multiplicative-level invariance failed')
  return m
 for rv in plan['RVs']:
  beta=parent['RV_km_s'][rv]*1000/299792458.
  for ci,rest in enumerate(plan['centres_air_A']):
   centre=float(rest*np.sqrt((1+beta)/(1-beta)))
   for scheme,bands in plan['shoulders'].items():
    for side,sign in [('blue',-1),('red',1)]:
     lo,hi=sorted([centre+sign*bands['shoulder_inner_A'],centre+sign*bands['shoulder_outer_A']]);candidates=[dict(line_id=r['canonical_row']['line_id'],species=r['canonical_row']['species'],wave_air_A=r['canonical_row']['wavelength_air_A'],treatment=r['treatment'])for r in intake if lo<=float(r['canonical_row']['wavelength_air_A'])<=hi];inventory.append(dict(rv=rv,rest_air_A=rest,centre_air_A=centre,scheme=scheme,side=side,lo_air_A=lo,hi_air_A=hi,candidates=candidates))
    for observation in plan['observation_legs']:
     ow,of=observations[ci,observation];m=measure(ow,of,centre,scheme);metrics[rv,ci,scheme,observation]=m;source_metrics.append(dict(kind='observation',variant='shipped',rv=rv,label=observation,rest_air_A=rest,centre_air_A=centre,scheme=scheme,**{k:v for k,v in m.items()if k!='scope'}))
    for variant in variant_parameters:
     for label in labels:
      flux=models[variant,rv,label];m=measure(mw,flux,centre,scheme);metrics[variant,rv,ci,scheme,label]=m;source_metrics.append(dict(kind='model',variant=variant,rv=rv,label=label,rest_air_A=rest,centre_air_A=centre,scheme=scheme,**{k:v for k,v in m.items()if k!='scope'}))
      for step in plan['grids_air_A']:
       offset=np.linspace(-plan['plot_half_A'],plan['plot_half_A'],int(round(2*plan['plot_half_A']/step))+1);grid=centre+offset;core=abs(offset)<=plan['core_half_A']+1e-12;mf=np.interp(grid,mw,flux);sampled=shape_metrics(grid,mf,centre,core_half_A=plan['core_half_A'],feature_half_A=plan['feature_half_A'],**bands);grid_deltas.append(dict(kind='model',variant=variant,rv=rv,label=label,rest_air_A=rest,scheme=scheme,step_air_A=step,contrast_delta=sampled['contrast']-m['contrast'],feature_area_delta_mA=sampled['relative_feature_area_mA']-m['relative_feature_area_mA']))
       if step==plan['grids_air_A'][0]:profiles[f'model_{variant}_{rv}_{ci}_{label}']=mf;profiles[f'wave_{rv}_{ci}']=grid
       for observation in plan['observation_legs']:
        ow,of=observations[ci,observation];df=np.interp(grid,ow,of);om=metrics[rv,ci,scheme,observation];osample=shape_metrics(grid,df,centre,core_half_A=plan['core_half_A'],feature_half_A=plan['feature_half_A'],**bands);comparisons.append(dict(variant=variant,rv=rv,label=label,rest_air_A=rest,centre_air_A=centre,scheme=scheme,observation=observation,step_air_A=step,n_core_grid_points=int(core.sum()),model_contrast=m['contrast'],observed_contrast=om['contrast'],contrast_residual=m['contrast']-om['contrast'],model_relative_core_area_mA=m['relative_core_area_mA'],observed_relative_core_area_mA=om['relative_core_area_mA'],model_relative_feature_area_mA=m['relative_feature_area_mA'],observed_relative_feature_area_mA=om['relative_feature_area_mA'],relative_feature_area_residual_mA=m['relative_feature_area_mA']-om['relative_feature_area_mA'],raw_core_rms=float(np.sqrt(np.mean((mf[core]-df[core])**2))),shoulder_ratio_core_rms=float(np.sqrt(np.mean((mf[core]/m['shoulder_level']-df[core]/om['shoulder_level'])**2))),observed_sampled_contrast_delta=osample['contrast']-om['contrast'],observed_sampled_feature_area_delta_mA=osample['relative_feature_area_mA']-om['relative_feature_area_mA']))
        if step==plan['grids_air_A'][0]:profiles[f'observed_{rv}_{ci}_{observation}']=df
 envelopes=[]
 for rv in plan['RVs']:
  for ci,rest in enumerate(plan['centres_air_A']):
   for scheme in plan['shoulders']:
    for label in labels:
     values=[metrics[variant,rv,ci,scheme,label]for variant in variant_parameters]
     for observation in plan['observation_legs']:
      obs=metrics[rv,ci,scheme,observation];lo=min(v['contrast']for v in values);hi=max(v['contrast']for v in values);al=min(v['relative_feature_area_mA']for v in values);ah=max(v['relative_feature_area_mA']for v in values);envelopes.append(dict(rv=rv,label=label,rest_air_A=rest,scheme=scheme,observation=observation,model_contrast_min=lo,model_contrast_max=hi,observed_contrast=obs['contrast'],observed_contrast_inside_tested_interval=lo<=obs['contrast']<=hi,model_feature_area_min_mA=al,model_feature_area_max_mA=ah,observed_feature_area_mA=obs['relative_feature_area_mA'],observed_area_inside_tested_interval=al<=obs['relative_feature_area_mA']<=ah))
 for path,receipt in checks:
  if fingerprint(path)!=receipt:raise ValueError('Protected source changed during milestone')
 a.output.mkdir(parents=True,exist_ok=False)
 for name,rows in [('metrics.csv',source_metrics),('comparisons.csv',comparisons),('sampling.csv',grid_deltas),('broadening_intervals.csv',envelopes)]:
  with(a.output/name).open('w',newline='')as stream:writer=csv.DictWriter(stream,fieldnames=list(rows[0]),lineterminator='\n');writer.writeheader();writer.writerows(rows)
 (a.output/'shoulder_inventory.json').write_text(json.dumps(inventory,indent=2,allow_nan=False)+'\n');np.savez(a.output/'profiles.npz',**profiles)
 import matplotlib
 matplotlib.use('Agg');import matplotlib.pyplot as plt
 fig,axes=plt.subplots(3,3,figsize=(15,11));rv=plan['plot_rv'];scheme=plan['plot_shoulder'];label='coupled_model_mass_H'
 for ci,rest in enumerate(plan['centres_air_A']):
  grid=profiles[f'wave_{rv}_{ci}'];centre=next(r['centre_air_A']for r in source_metrics if r['rv']==rv and r['rest_air_A']==rest);offset=grid-centre
  for observation in plan['observation_legs']:
   f=profiles[f'observed_{rv}_{ci}_{observation}'];m=metrics[rv,ci,scheme,observation];axes[0,ci].plot(offset,f,label=observation,lw=1);axes[1,ci].plot(offset,f/m['shoulder_level'],label=observation,lw=1)
  f=profiles[f'model_nominal_{rv}_{ci}_{label}'];m=metrics['nominal',rv,ci,scheme,label];axes[0,ci].plot(offset,f,color='black',label='Model-H nominal');axes[1,ci].plot(offset,f/m['shoulder_level'],color='black',label='Model-H nominal')
  for variant in variant_parameters:
   f=profiles[f'model_{variant}_{rv}_{ci}_{label}'];m=metrics[variant,rv,ci,scheme,label];axes[2,ci].plot(offset,f/m['shoulder_level'],label=variant,lw=1)
  data=profiles[f'observed_{rv}_{ci}_kitt_shipped'];om=metrics[rv,ci,scheme,'kitt_shipped'];axes[2,ci].plot(offset,data/om['shoulder_level'],color='black',ls='--',label='Kitt')
  axes[0,ci].set_title(f'{rest:.3f} AIR A');axes[0,ci].set_ylabel('Shipped normalized flux');axes[1,ci].set_ylabel('Fixed-shoulder ratio');axes[2,ci].set_ylabel('Broadening OAT controls / own shoulders')
  for ax in axes[:,ci]:ax.set_xlabel('AIR offset from fixed +gravity centre (A)');ax.axvspan(-plan['core_half_A'],plan['core_half_A'],alpha=.08,color='gray');ax.legend(fontsize=6)
 fig.suptitle('O777 shape milestone: no optimized widths/RV/continuum/abundance\nAlternative shoulders and both RVs retained in tables; repeated references are correlated');fig.tight_layout();fig.savefig(a.output/'shape_milestone.png',dpi=140);plt.close(fig)
 result=dict(schema='codex.bride_O_shape_milestone/1',scientific_admission='HOLD',generator=fingerprint(__file__),metric_helper=fingerprint(ROOT/'pipeline/bride_feii_profile.py'),plan=fingerprint(plan_path),source_files=plan['source_files'],scope=plan['scope'],milestone_acceptance=plan['milestone_acceptance'],orch_commit=pin,policies=policies,frame=parent['frame'],observation=parent['observation'],iag=parent['iag'],iag_reader=fingerprint(reader.__file__),broadening_variants=variant_parameters,RV_km_s=parent['RV_km_s'],iag_flags=flag_records,metrics=len(source_metrics),comparisons=len(comparisons),sampling_records=len(grid_deltas),broadening_interval_records=len(envelopes),nominal_replays=12,scale_canaries=len(canary_contrast),max_scale_contrast_error=max(canary_contrast),max_scale_feature_area_error_mA=max(canary_area),padding=padding,telemetry=dict(seconds=time.perf_counter()-started,threads=1),limitations=['No inverse fit or selection of an optimal broadening/RV/atlas/shoulder; tested intervals are not physical uncertainty bounds.','Shoulders contain neighboring absorption and possibly other triplet wings; candidate positions are not blend identifications.','Common-grid residuals are unweighted and correlated; finer interpolation is not more resolving power.','Shape ratios remove a common multiplicative level only, not wavelength-dependent continua, missing opacity or all blend/atmosphere/source errors.','Native and Gaussian IAG are correlated views of one holding; source normalization-error independence unestablished.','No full3D/all-band/abundance/model/default/product/production/publication/merge admission.'])
 (a.output/'report.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(json.dumps({k:result[k]for k in ['metrics','comparisons','sampling_records','nominal_replays','scale_canaries','max_scale_contrast_error','max_scale_feature_area_error_mA','telemetry']}))

if __name__=='__main__':main()
