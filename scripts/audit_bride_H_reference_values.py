"""Same-pixel reference-value accounting and held-source lineage; no rescaling."""
import argparse,csv,hashlib,importlib.util,json,os,subprocess,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
 from astropy.io import fits
 from scipy.ndimage import gaussian_filter1d
 from config.constants import require_codex_path
 from pipeline.bride_reference import fingerprint
 from pipeline.bride_continuum_audit import support_geometry
 from pipeline.bride_crossed_continuum import finite_accounting
 from pipeline.bride_observables import window_integral
 from pipeline.wavelength_util import vac_to_air
 from scripts.audit_bride_kitt_residuals import integral
 plan_path=ROOT/'config/experiments/rya1281_H_reference_values.json';plan=json.loads(plan_path.read_text());checks=[]
 def check(path,receipt):
  if fingerprint(path)!=receipt:raise ValueError('Changed held source '+str(path))
  checks.append((path,receipt))
 for name,receipt in plan['source_files'].items():check(ROOT/name,receipt)
 for name,receipt in plan['papers'].items():check(require_codex_path('refs.library')/name,receipt)
 source=ROOT/'data/results/rya1281/H_kitt_orch_v1';parent=json.loads((source/'report.json').read_text());fixed=json.loads((ROOT/'data/results/rya1281/H_bracket_support_v1/fits.json').read_text());lookup={(r['label'],r['rv'],r['centre_air_A'],r['reference']):r['fit']for r in fixed if r['support']==plan['support']};labels=[r['label']for r in parent['sources']]
 with np.load(source/'profiles_and_continua.npz')as data:arrays={k:data[k].copy()for k in data.files}
 original={k:v.copy()for k,v in arrays.items()};w=arrays['observed_wave_air_A'];obs=arrays['observed_shipped_flux'];atlas=require_codex_path('data.solar_iag_baker2020');kitt=require_codex_path('data.kurucz2005_residual');check(atlas,parent['iag']);check(kitt,parent['observation'])
 orch=require_codex_path('experiments.bride_current_orch_checkout');native=orch/'pipeline/local_continuum.py';pin=parent['orch_commit'];check(native,parent['current_native'])
 if subprocess.check_output(['git','rev-parse','HEAD'],cwd=orch,text=True).strip()!=pin or native.read_bytes()!=subprocess.check_output(['git','show',pin+':pipeline/local_continuum.py'],cwd=orch):raise ValueError('Changed pinned helper')
 os.environ['CODEX_CONT_REF']='iag';spec=importlib.util.spec_from_file_location('reference_values_pinned_continuum',native);lc=importlib.util.module_from_spec(spec);sys.modules[spec.name]=lc;spec.loader.exec_module(lc)
 with fits.open(atlas)as hdus:
  headers=[str(h.header)for h in hdus];data=hdus[1].data;v=np.array(data['v'],float);s=np.array(data['s'],float);flags=np.array(data['flags'],float);columns=list(data.names)
 if not np.isfinite(v).all()or np.any(v<=0)or not np.isfinite(s).all()or not np.isfinite(flags).all():raise ValueError('Invalid held IAG arrays')
 order=np.argsort(vac_to_air(1e8/v));iw=vac_to_air(1e8/v)[order];fl=s[order];q=flags[order];source_replays=[]
 for ci,centre in enumerate(plan['centres_air_A']):
  rw,rf=lc.solar_atlas_reference(centre-5,centre+5,parent['broadening']['R']);reader=sys.modules['measure_band_ew']
  if reader.IAG_FITS.resolve()!=atlas.resolve():raise ValueError('Wrong pinned IAG holding')
  keep=(iw>=centre-6)&(iw<=centre+6);pw=iw[keep];pf=fl[keep];step=float(np.median(np.diff(pw)));fwhm=.5*((centre-5)+(centre+5))/float(parent['broadening']['R']);rebuilt=gaussian_filter1d(pf,fwhm/2.3548/step)
  np.testing.assert_array_equal(pw,rw);np.testing.assert_array_equal(rebuilt,rf)
  sampled=np.interp(w,pw,rebuilt)
  for rv in plan['RVs']:
   for label in labels:np.testing.assert_array_equal(sampled,arrays[f'reference_{rv}_{ci}_iag_{label}'])
  values,counts=np.unique(q[keep],return_counts=True);source_replays.append(dict(centre_air_A=centre,native_count=len(pw),native_median_step_air_A=step,gaussian_sigma_pixels=fwhm/2.3548/step,flag_counts={str(float(x)):int(n)for x,n in zip(values,counts)},rebuilt_helper_and_archived_samples='exact array identity',scope='Same-source format/reader replay, not independent physical reference or calibrated FTS LSF'))
 rows=[];pixels=[];plot_data=[]
 for rv in plan['RVs']:
  for ci,centre in enumerate(plan['centres_air_A']):
   for label in labels:
    A=lookup[label,rv,centre,'synthesis'];B=lookup[label,rv,centre,'iag']
    if not A['applied']or not B['applied']or A['indices']!=B['indices']:raise ValueError('Unmatched accepted balanced support')
    idx=np.array(A['indices'],int);x=w-centre;syn=arrays[f'reference_{rv}_{ci}_synthesis_{label}'];iag=arrays[f'reference_{rv}_{ci}_iag_{label}'];
    if any(not np.isfinite(z[idx]).all()or np.any(z[idx]<=0)for z in [obs,syn,iag]):raise ValueError('Invalid selected reference values')
    ratio_s=obs[idx]/syn[idx];ratio_i=obs[idx]/iag[idx];point=obs[idx]*(syn[idx]-iag[idx])/(iag[idx]*syn[idx])
    np.testing.assert_allclose(ratio_i-ratio_s,point,rtol=0,atol=16*np.finfo(float).eps)
    geometry=support_geometry(w,centre,idx);weights=np.array(geometry['centre_weights']);cs=np.polyfit(x[idx],ratio_s,1);ct=np.polyfit(x[idx],ratio_i,1);diff=np.polyfit(x[idx],point,1);bound=1024*np.finfo(float).eps*geometry['design_condition_2']
    np.testing.assert_array_equal(cs,np.array([A['slope_per_A'],A['level']]));np.testing.assert_array_equal(ct,np.array([B['slope_per_A'],B['level']]))
    np.testing.assert_allclose(ct-cs,diff,rtol=0,atol=bound);weighted=float(weights@point)
    if abs(weighted-(B['level']-A['level']))>bound:raise ValueError('Level influence reconstruction failed')
    area=lambda level,slope:1000*window_integral(w,obs/(level+slope*x),centre-plan['half_width_A'],centre+plan['half_width_A'])
    corners={lev+'_'+sl:dict(area=area(level,slope))for lev,level in [('before',A['level']),('after',B['level'])]for sl,slope in [('before',A['slope_per_A']),('after',B['slope_per_A'])]};account=finite_accounting(corners,'area');ca=A['level']+A['slope_per_A']*x;cb=B['level']+B['slope_per_A']*x;direct=1000*integral(w,obs*(1/ca-1/cb),centre-plan['half_width_A'],centre+plan['half_width_A'])
    if not account['available']or abs(account['closure'])>1e-10 or abs(direct-account['total'])>1e-10:raise ValueError('Finite observed-area identity failed')
    base=dict(label=label,rv=rv,centre_air_A=centre,count=len(idx));row=dict(**base,synthesis_level=A['level'],iag_level=B['level'],level_difference=B['level']-A['level'],weighted_ratio_difference=weighted,slope_difference_per_A=B['slope_per_A']-A['slope_per_A'],observed_area_difference_mA=account['total'],direct_area_difference_mA=direct,level_path_average_mA=account['support_path_average'],slope_path_average_mA=account['denominator_path_average'],level_slope_interaction_mA=account['interaction'],closure_mA=account['closure'])
    for key,values in [('kitt',obs[idx]),('synthesis',syn[idx]),('iag',iag[idx]),('iag_minus_synthesis',(iag-syn)[idx]),('kitt_over_synthesis',ratio_s),('kitt_over_iag',ratio_i)]:
     row[key+'_median']=float(np.median(values));row[key+'_mean']=float(np.mean(values));row[key+'_min']=float(np.min(values));row[key+'_max']=float(np.max(values))
    rows.append(row)
    for pos,(j,weight)in enumerate(zip(idx,weights)):pixels.append(dict(**base,index=int(j),wave_air_A=float(w[j]),offset_air_A=float(x[j]),side='blue'if x[j]<0 else'red',kitt=float(obs[j]),synthesis=float(syn[j]),iag=float(iag[j]),kitt_over_synthesis=float(ratio_s[pos]),kitt_over_iag=float(ratio_i[pos]),ratio_difference=float(ratio_i[pos]-ratio_s[pos]),algebraic_ratio_difference=float(point[pos]),centre_influence_weight=float(weight)))
    if label==plan['baseline']and rv=='plus_gravity':plot_data.append((centre,x[idx],obs[idx],syn[idx],iag[idx],ratio_s,ratio_i))
 for key,value in arrays.items():np.testing.assert_array_equal(value,original[key])
 for path,receipt in checks:
  if fingerprint(path)!=receipt:raise ValueError('Source changed during value audit')
 a.output.mkdir(parents=True,exist_ok=False)
 for name,data in [('summary.csv',rows),('pixels.csv',pixels)]:
  with(a.output/name).open('w',newline='')as stream:writer=csv.DictWriter(stream,fieldnames=list(data[0]),lineterminator='\n');writer.writeheader();writer.writerows(data)
 with kitt.open()as stream:kitt_header=[next(stream).rstrip()for _ in range(3)]
 lineage=dict(schema='codex.bride_H_reference_lineage/1',papers=plan['papers'],kitt=dict(receipt=parent['observation'],header=kitt_header,interpretation='Shipped residual irradiance, vacuum nm including gravitational redshift; local canonical AIR conversion once; no own continuum or new correction.'),iag=dict(receipt=parent['iag'],fits_headers=headers,columns=columns,declared_column_units_in_FITS=False,reader=fingerprint(reader.__file__),reader_interpretation='v is vacuum wavenumber cm^-1, s is solar transmission/normalized flux; canonical 1e8/v then vac_to_air once; source construction includes continuum fitting; no absolute irradiance units assigned.',publisher_reader_present=(atlas.parent/'open_data.py').exists(),format_limit='FITS column names carry no TUNIT metadata; exact source-file lineage and primary paper/reader convention retained, not invented unit headers.'),source_notes=['Baker2020 p1-2: disk-integrated solar flux and approximate resolving power1e6, not disk-centre intensity.','Baker2020 p5-7 sections3.1.4/4.1/4.2: solar-template normalization and fitted continuum corrections are part of construction; final spectra combined after ephemeris/iodine/zero-point shifts. This does not prove the local offset cause.','Baker2020 p7-8: some telluric residual regions replaced with solar spline and flagged; exact zero O777 flags do not certify physical accuracy.','Reiners2016 p4 section5.2: normalized flux atlases have distinct normalization procedures; illustrative agreement elsewhere is not O777 calibration.','Historical RYA929 intake manifest describes absolute irradthu, not this residual irradrelwl; old unresolved medium/absolute units are not reassigned to the current residual file.'],resolution_limit='Retained helper GaussianR500000 is numerical reference smoothing; approximate native R is not an independently calibrated FTS kernel. No deconvolution, new kernel or RV calibration.',source_replays=source_replays)
 (a.output/'lineage.json').write_text(json.dumps(lineage,indent=2,allow_nan=False)+'\n')
 import matplotlib
 matplotlib.use('Agg');import matplotlib.pyplot as plt
 fig,axes=plt.subplots(2,3,figsize=(15,8))
 for ci,(centre,offset,k,s,i,rs,ri)in enumerate(plot_data):
  for values,label in [(k,'Kitt'),(s,'Bare synthetic reference'),(i,'IAG reference')]:axes[0,ci].scatter(offset,values,s=9,label=label)
  for values,label in [(rs,'Kitt / synthesis'),(ri,'Kitt / IAG')]:axes[1,ci].scatter(offset,values,s=9,label=label)
  axes[0,ci].set_title(f'{centre:.3f} AIR A');axes[0,ci].set_ylabel('Existing reference/observed normalized values');axes[1,ci].set_ylabel('Same-pixel ratios')
  for ax in axes[:,ci]:ax.set_xlabel('AIR offset from fixed centre (A)');ax.legend(fontsize=8)
 fig.suptitle('Unchanged balanced pixels: reference values and ratios\nNo atlas rescaling, new continuum, physical attribution or abundance adoption');fig.tight_layout();fig.savefig(a.output/'reference_values.png',dpi=140);plt.close(fig)
 result=dict(schema='codex.bride_H_reference_values/1',scientific_admission='HOLD',scope=plan['scope'],generator=fingerprint(__file__),plan=fingerprint(plan_path),source_files=plan['source_files'],papers=plan['papers'],orch_commit=pin,rows=len(rows),pixels=len(pixels),source_replays=source_replays,max_level_identity_error=float(max(abs(r['level_difference']-r['weighted_ratio_difference'])for r in rows)),max_area_identity_error_mA=float(max(abs(r['observed_area_difference_mA']-r['direct_area_difference_mA'])for r in rows)),max_accounting_closure_mA=float(max(abs(r['closure_mA'])for r in rows)),limitations=['Repeated spectra/RV labels reuse the same observation, not independent evidence.','Exact algebra/source-reader parity does not establish which reference continuum or model is physically correct.','Level/slope finite paths are correlated bookkeeping, not independent uncertainty or causal terms.','No source/normalization/RV/broadening/physics/default/product change, abundance or all-band admission.'])
 (a.output/'report.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(json.dumps({k:result[k]for k in ['rows','pixels','max_level_identity_error','max_area_identity_error_mA','max_accounting_closure_mA']}))

if __name__=='__main__':main()
