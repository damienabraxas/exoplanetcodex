"""Fixed support/reference controls on archived paired-H Kitt Peak diagnostics."""
import argparse,csv,importlib.util,json,subprocess,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
 from config.constants import require_codex_path
 from pipeline.bride_reference import fingerprint
 from pipeline.bride_fixed_support import shared_masks,fixed_fit
 from pipeline.bride_observables import window_integral
 from pipeline.bride_crossed_continuum import finite_accounting
 plan_path=ROOT/'config/experiments/rya1281_H_bracket_support.json';plan=json.loads(plan_path.read_text());source=ROOT/plan['source_archive'];checks=[]
 def check(path,receipt):
  if fingerprint(path)!=receipt:raise ValueError('Changed source '+str(path))
  checks.append((path,receipt))
 for name,receipt in plan['source_files'].items():check(source/name,receipt)
 parent=json.loads((source/'report.json').read_text());supports=json.loads((source/'support.json').read_text());orch=require_codex_path('experiments.bride_current_orch_checkout');pin=parent['orch_commit'];native=orch/'pipeline/local_continuum.py'
 if subprocess.check_output(['git','rev-parse','HEAD'],cwd=orch,text=True).strip()!=pin or native.read_bytes()!=subprocess.check_output(['git','show',pin+':pipeline/local_continuum.py'],cwd=orch):raise ValueError('Changed pinned ORCH')
 check(native,parent['current_native']);check(require_codex_path('data.kurucz2005_residual'),parent['observation']);check(require_codex_path('data.solar_iag_baker2020'),parent['iag'])
 spec=importlib.util.spec_from_file_location('bracket_pinned_continuum',native);lc=importlib.util.module_from_spec(spec);sys.modules[spec.name]=lc;spec.loader.exec_module(lc)
 if(lc.MAX_SHIFT_MODEL,lc.MODEL_Q,lc.DEGREE)!=(.15,90.,1)or plan['clip']or plan['selector_percentile']!=lc.MODEL_Q:raise ValueError('Changed control contract')
 with np.load(source/'profiles_and_continua.npz')as data:arrays={k:data[k].copy()for k in data.files}
 w=arrays['observed_wave_air_A'];obs=arrays['observed_shipped_flux'];mw=arrays['model_wave_air_A'];before={k:v.copy()for k,v in arrays.items()};baseline=plan['baseline'];labels=[s['label']for s in parent['sources']];lookup={(r['rv'],r['centre_air_A'],r['reference'],r['label']):r for r in supports};records=[];rows=[];continua={};replay=[]
 for rv in plan['RVs']:
  for ci,centre in enumerate(plan['centres_air_A']):
   references={ref:arrays[f'reference_{rv}_{ci}_{ref}_{baseline}']for ref in plan['references']};shared=shared_masks(w,centre,references['synthesis'],references['iag'],inner=plan['half_width_A'],outer=plan['continuum_half_width_A'],percentile=plan['selector_percentile'])
   eligible=(abs(w-centre)>plan['half_width_A'])&(abs(w-centre)<=plan['continuum_half_width_A'])
   for label in labels:
    syn=lookup[rv,centre,'synthesis',label];iag=lookup[rv,centre,'iag',label];masks={}
    for name,indices in [('native_selected_synthesis',syn['selected_indices']),('native_coeff_synthesis',syn['fit_indices']),('native_coeff_iag',iag['fit_indices'])]:
     mask=np.zeros(len(w),bool);mask[indices]=True;masks[name]=mask
    masks.update(shared);model=arrays['model_'+rv+'_'+label];model_area=1000*window_integral(mw,model,centre-plan['half_width_A'],centre+plan['half_width_A']);mf=np.interp(w,mw,model);window=abs(w-centre)<=plan['half_width_A']
    for name in plan['supports']:
     mask=masks[name]
     if np.any(mask&~eligible):raise ValueError('Support outside eligible sidebands')
     for ref in plan['references']:
      denominator=arrays[f'reference_{rv}_{ci}_{ref}_{label}'];is_shared=name in shared;fit,cont=fixed_fit(w,obs,centre,denominator,mask,maximum_shift=lc.MAX_SHIFT_MODEL,minimum_pixels=lc.MIN_MODEL_PIX,minimum_each_side=plan['minimum_each_side'],require_bracket=is_shared)
      identity=dict(label=label,rv=rv,centre_air_A=centre,support=name,reference=ref)
      if name=='native_coeff_'+ref:
       original=lookup[rv,centre,ref,label];bound=512*np.finfo(float).eps*original['geometry']['design_condition_2'];dl=abs(fit['level']-original['pinned_record']['level_at_centre']);ds=abs(fit['slope_per_A']-original['pinned_record']['slope_per_A'])
       if dl>bound or ds>bound or fit['applied']!=original['pinned_record']['applied']:raise ValueError('Native coefficient replay changed')
       replay.append(dict(**identity,level_difference=dl,slope_difference=ds,bound=bound))
      if cont is not None:
       np.testing.assert_array_equal(fit['indices'],np.flatnonzero(mask));conditioned=obs/cont;area=1000*window_integral(w,conditioned,centre-plan['half_width_A'],centre+plan['half_width_A']);rms=float(np.sqrt(np.mean((mf[window]-conditioned[window])**2)))
      else:area=None;rms=None
      row=dict(**identity,applied=fit['applied'],reason=fit['reason'],level=fit['level'],slope_per_A=fit['slope_per_A'],count=fit['count'],blue_count=fit['blue_count'],red_count=fit['red_count'],require_bracket=is_shared,centre_outside_hull_A=None if fit['geometry']is None else fit['geometry']['centre_outside_hull_A'],centre_weight_l1=None if fit['geometry']is None else fit['geometry']['centre_weight_l1'],observed_area_mA=area,model_area_mA=model_area,area_residual_mA=None if area is None else model_area-area,rms_residual=rms)
      rows.append(row);records.append(dict(**identity,fit=fit));continua[rv,ci,name,ref,label]=cont
   for name in shared:
    values=[continua[rv,ci,name,'iag',label]for label in labels]
    for cont in values[1:]:
     if(cont is None)!=(values[0]is None):raise ValueError('Shared IAG refusal depends on model label')
     if cont is not None:np.testing.assert_array_equal(cont,values[0])
 for k,v in arrays.items():np.testing.assert_array_equal(v,before[k])
 contrasts=[]
 for rv in plan['RVs']:
  for centre in plan['centres_air_A']:
   for name in plan['supports']:
    by={(r['label'],r['reference']):r for r in rows if r['rv']==rv and r['centre_air_A']==centre and r['support']==name}
    for label in labels:
     A=by[label,'synthesis'];B=by[label,'iag'];contrasts.append(dict(kind='reference',label=label,rv=rv,centre_air_A=centre,support=name,reference='iag_minus_synthesis',area_contrast_mA=None if A['observed_area_mA']is None or B['observed_area_mA']is None else B['observed_area_mA']-A['observed_area_mA']))
    for anchor in ['coupled','model_anchor','native_anchor']:
     for ref in plan['references']:
      A=by[anchor+'_model_mass_H',ref];B=by[anchor+'_native_network_H',ref];contrasts.append(dict(kind='H_residual',label=anchor,rv=rv,centre_air_A=centre,support=name,reference=ref,area_contrast_mA=None if A['area_residual_mA']is None or B['area_residual_mA']is None else B['area_residual_mA']-A['area_residual_mA']))
 accounting=[]
 for rv in plan['RVs']:
  for centre in plan['centres_air_A']:
   for label in labels:
    by={(r['support'],r['reference']):r for r in rows if r['rv']==rv and r['centre_air_A']==centre and r['label']==label}
    corners={support+'_'+ref:dict(area=by[name,reference]['observed_area_mA'])for support,name in [('before','native_coeff_synthesis'),('after','balanced_iag_per_side')]for ref,reference in [('before','synthesis'),('after','iag')]}
    change=finite_accounting(corners,'area')
    if change['available']and abs(change['closure'])>1e-10:raise ValueError('Finite accounting closure failed')
    accounting.append(dict(label=label,rv=rv,centre_air_A=centre,**change))
 for path,receipt in checks:
  if fingerprint(path)!=receipt:raise ValueError('Source changed during fixed-support analysis')
 a.output.mkdir(parents=True,exist_ok=False)
 for name,data in [('summary.csv',rows),('contrasts.csv',contrasts),('accounting.csv',accounting)]:
  with(a.output/name).open('w',newline='')as stream:writer=csv.DictWriter(stream,fieldnames=list(data[0]),lineterminator='\n');writer.writeheader();writer.writerows(data)
 (a.output/'fits.json').write_text(json.dumps(records,indent=2,allow_nan=False)+'\n')
 import matplotlib
 matplotlib.use('Agg');import matplotlib.pyplot as plt
 fig,axes=plt.subplots(2,3,figsize=(15,8));rv=plan['plot_rv'];name=plan['plot_support']
 for ci,centre in enumerate(plan['centres_air_A']):
  mask=abs(w-centre)<=.65;mk=abs(mw-centre)<=.65;model=arrays['model_'+rv+'_'+baseline];axes[0,ci].plot(mw[mk],model[mk],color='black',label='Coupled model-H')
  for ref in plan['references']:
   cont=continua[rv,ci,name,ref,baseline]
   if cont is not None:axes[0,ci].plot(w[mask],(obs/cont)[mask],label='Kitt / fixed '+ref+' denominator')
   else:axes[0,ci].text(.1,.2,'Refused '+ref,transform=axes[0,ci].transAxes)
   r=next(r for r in records if r['label']==baseline and r['rv']==rv and r['centre_air_A']==centre and r['support']==name and r['reference']==ref);idx=r['fit']['indices'];axes[1,ci].scatter(w[idx]-centre,obs[idx]/arrays[f'reference_{rv}_{ci}_{ref}_{baseline}'][idx],s=8,label=ref+' ratios')
  axes[0,ci].set_title(f'{centre:.3f} AIR A');axes[0,ci].set_xlabel('AIR A');axes[0,ci].set_ylabel('Normalized flux');axes[1,ci].set_xlabel('Offset from line centre (AIR A)');axes[1,ci].set_ylabel('Observed/reference on fixed shared pixels')
  for ax in axes[:,ci]:ax.legend(fontsize=8)
 fig.suptitle('Diagnostic balanced IAG-per-side support, no clipping; fixed plus-gravity models\nNo continuum/abundance adoption; original native supports retained separately');fig.tight_layout();fig.savefig(a.output/'fixed_support.png',dpi=140);plt.close(fig)
 result=dict(schema='codex.bride_H_bracket_support_analysis/1',scientific_admission='HOLD',scope=plan['scope'],generator=fingerprint(__file__),helper=fingerprint(ROOT/'pipeline/bride_fixed_support.py'),plan=fingerprint(plan_path),source_files=plan['source_files'],orch_commit=pin,pinned_native=fingerprint(native),maximum_shift=lc.MAX_SHIFT_MODEL,minimum_pixels=lc.MIN_MODEL_PIX,observation=parent['observation'],iag=parent['iag'],frame=parent['frame'],broadening=parent['broadening'],RV_km_s=parent['RV_km_s'],rows=len(rows),applied=sum(r['applied']for r in rows),refused=sum(not r['applied']for r in rows),native_coefficient_replays=replay,by_support={name:dict(applied=sum(r['applied']for r in rows if r['support']==name),refused=sum(not r['applied']for r in rows if r['support']==name))for name in plan['supports']},limitations=['Fixed support is an intervention, not a new validated continuum or ratified estimator.','Shared bracketing criteria are prospective diagnostic guards, not physical error bounds.','No clipping; original native selected/exact coefficient supports retained as controls.','No new broadening, RV, observation conditioning outside saved local diagnostics, RT/EOS/provider/abundance or product admission.','Finite areas/unweighted residuals have no noise/significance/abundance interpretation.'])
 result['finite_accounting']=dict(quartets=len(accounting),source_support='native_coeff_synthesis',target_support='balanced_iag_per_side',scope='Finite path-average changes in conditioned observed area only; correlated interventions, not unique causal/error-budget terms',helper=fingerprint(ROOT/'pipeline/bride_crossed_continuum.py'))
 (a.output/'report.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(json.dumps({k:result[k]for k in ['rows','applied','refused','by_support']}))

if __name__=='__main__':main()
