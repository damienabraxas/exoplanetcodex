"""Replay actual pinned RYA1315 ORCH continuum on frozen Bride inputs."""
import argparse,csv,json,os,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--orch-worktree',type=Path,required=True);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    from config.constants import require_codex_path
    from pipeline.bride_reference import fingerprint
    from pipeline.bride_continuum_audit import selection
    from pipeline.kp_atlas_integrity import read_kurucz2005_residual
    from pipeline.reflected_solar_rv import SOLAR_GRAVITATIONAL_REDSHIFT
    from scripts.audit_bride_kitt_residuals import integral
    os.environ['CODEX_CONT_REF']='iag'
    from pipeline import local_continuum as lc
    sys.path.insert(0,str(require_codex_path('engines.ispec')));import ispec
    plan_path=ROOT/'config/experiments/rya1281_current_orch.json';plan=json.loads(plan_path.read_text())
    if plan['schema']!='codex.bride_current_orch_plan/1':raise ValueError('Wrong crossed plan')
    import importlib.util,subprocess
    sha=subprocess.check_output(['git','rev-parse','HEAD'],cwd=args.orch_worktree,text=True).strip()
    if sha!=plan['orch_commit']:raise ValueError('Wrong pinned ORCH checkout')
    native_path=args.orch_worktree/'pipeline/local_continuum.py'
    blob=subprocess.check_output(['git','show',sha+':pipeline/local_continuum.py'],cwd=args.orch_worktree)
    if native_path.read_bytes()!=blob:raise ValueError('Changed current ORCH source')
    branch_blob=subprocess.check_output(['git','show',plan['orch_branch']+':pipeline/local_continuum.py'],cwd=ROOT)
    if branch_blob!=blob:raise ValueError('RYA1315 branch continuum differs from pinned main')
    spec=importlib.util.spec_from_file_location('rya1315_native_continuum',native_path)
    current=importlib.util.module_from_spec(spec);sys.modules[spec.name]=current;spec.loader.exec_module(current)
    if (current.MAX_SHIFT,current.MAX_SHIFT_MODEL,current.MODEL_Q,current.ENV_PERCENTILE)!=(.03,.15,90.,95.):raise ValueError('Unexpected current ORCH controls')
    x=np.linspace(-5,5,1001);unit=np.ones_like(x);synthetic=[]
    for level in (.92,0.):
        before,_=lc.fit_model_guided(x,level*unit,0.,x,unit,exclude_half_width_A=.5)
        after,_=current.fit_model_guided(x,level*unit,0.,x,unit,exclude_half_width_A=.5)
        synthetic.append(dict(level=level,historical=before.__dict__,current=after.__dict__))
    envelope,_=current.fit(x,.92*unit,0.)
    if synthetic[0]['historical']['applied'] or not synthetic[0]['current']['applied'] or synthetic[1]['current']['applied'] or envelope.applied:raise ValueError('Current ORCH bound canary failed')
    comparison_path=ROOT/'data/results/rya1281/solar_bound_v1/kitt_peak/report.json'
    comparison=json.loads(comparison_path.read_text());b=comparison['broadening']
    observation=require_codex_path('data.kurucz2005_residual')
    if fingerprint(observation)!=comparison['observation']:raise ValueError('Changed delivered observation')
    old=json.loads((ROOT/'data/results/rya1281/continuum_support_v1/report.json').read_text())
    if fingerprint(lc.__file__)!=old['native'] or lc.MAX_SHIFT!=.03:raise ValueError('Changed historical continuum rule')
    w,obs=read_kurucz2005_residual(observation,7765,7782)
    from scripts.measure_band_ew import IAG_FITS
    iag=require_codex_path('data.solar_iag_baker2020')
    profile=json.loads((ROOT/'data/results/rya1281/feii_profile_source_v1/profile/report.json').read_text())
    if IAG_FITS.resolve()!=iag.resolve() or fingerprint(iag)!=profile['iag']:raise ValueError('Changed IAG holding')
    result=dict(schema='codex.bride_current_orch/1',scientific_admission='HOLD',generator=fingerprint(__file__),
        trace_kernel=fingerprint(ROOT/'pipeline/bride_continuum_audit.py'),
        plan=fingerprint(plan_path),comparison=fingerprint(comparison_path),native=fingerprint(lc.__file__),
        observation=fingerprint(observation),iag=fingerprint(iag),iag_reader=fingerprint(ROOT/'scripts/measure_band_ew.py'),
        scope=plan['scope'],orch_commit=sha,current_native=fingerprint(native_path),synthetic_controls=synthetic,envelope_sag_control=envelope.__dict__,sources=[],records=[])
    rows=[];replays=0
    for mode in plan['stark_modes']:
        transfer=require_codex_path('experiments.bride_feii_width_transfer')/mode
        analysis_root=require_codex_path('experiments.bride_feii_width_analysis')/mode
        anchor_path=require_codex_path('experiments.bride_feii_width_anchors')/mode/'report.json'
        analysis_path=analysis_root/'report.json';analysis=json.loads(analysis_path.read_text())
        anchor=json.loads(anchor_path.read_text())
        if anchor['analysis']!=fingerprint(analysis_path) or anchor['adapter']!=result['trace_kernel']:raise ValueError('Changed source trace')
        if analysis['stark_mode']!=mode:raise ValueError('Wrong Stark alternative')
        result['sources'].append(dict(stark_mode=mode,analysis=fingerprint(analysis_path),anchors=fingerprint(anchor_path)))
        for name in plan['cases']:
            source=next(c for c in analysis['cases'] if c['name']==name)
            saved=next(c['traces'] for c in anchor['cases'] if c['name']==name)
            run_path=transfer/name/'report.json';run=json.loads(run_path.read_text())
            if fingerprint(run_path)!=source['source']:raise ValueError('Changed source transfer')
            models={}
            for identity,leg in [('before','zero'),('after','feii')]:
                path=transfer/name/leg/'flux.csv'
                if fingerprint(path)!=next(r['flux'] for r in run['legs'] if r['mode']==leg):raise ValueError('Changed source spectrum')
                raw=np.genfromtxt(path,delimiter=',',names=True);mw=raw['wavelength_air_A']
                segment=np.array([(mw[0]/10,mw[-1]/10)],dtype=[('wave_base',float),('wave_top',float)])
                mf=ispec.apply_post_fundamental_effects(mw/10,raw['normalized_flux'].copy(),segment,
                    macroturbulence=b['vmac_km_s'],vsini=b['vsini_km_s'],limb_darkening_coeff=b['limb_darkening'],R=b['R'],vrad=(SOLAR_GRAVITATIONAL_REDSHIFT,))
                models[identity]=(mw,mf)
            for reference in plan['references']:
                for centre in plan['centres_air_A']:
                    refs={};selected={};coefficients={};native={}
                    for identity,leg in [('before','zero'),('after','feii')]:
                        apath=analysis_root/name/leg/reference/'report.json'
                        if fingerprint(apath)!=source['audit_receipts'][leg][reference]:raise ValueError('Changed native audit')
                        rw,rf=models[identity] if reference=='synthesis' else lc.solar_atlas_reference(centre-5,centre+5,b['R'])
                        keep=(rw>=centre-5.01)&(rw<=centre+5.01)
                        refs[identity]=np.interp(w,rw[keep],rf[keep])
                        selected[identity],_=selection(w,obs,centre,refs[identity])
                        native[identity]=saved[leg][str(centre)+'_'+reference]
                        np.testing.assert_array_equal(np.flatnonzero(selected[identity]),native[identity]['selection_indices'])
                        mask=np.zeros(len(w),bool);mask[native[identity]['fit_indices']]=True;coefficients[identity]=mask
                    for identity,leg in [('before','zero'),('after','feii')]:
                        old_record,old_cont=lc.fit_model_guided(w,obs,centre,w,refs[identity],exclude_half_width_A=.5,reference=reference)
                        new_record,new_cont=current.fit_model_guided(w,obs,centre,w,refs[identity],exclude_half_width_A=.5,reference=reference)
                        expected=native[identity]
                        old_area=None if old_cont is None else 1000*integral(w,1-obs/old_cont,centre-.5,centre+.5)
                        if (old_record.level_at_centre,old_record.slope_per_A,old_record.n_bins,old_record.applied,old_area)!=(expected['level'],expected['slope_per_A'],expected['reported_n_bins'],expected['applied'],expected.get('observed_area_mA')):raise ValueError('Historical native replay differs')
                        replays+=1
                        new_area=None if new_cont is None else 1000*integral(w,1-obs/new_cont,centre-.5,centre+.5)
                        exact=old_record.__dict__==new_record.__dict__ and ((old_cont is None and new_cont is None) or (old_cont is not None and new_cont is not None and np.array_equal(old_cont,new_cont)))
                        item=dict(stark_mode=mode,case=name,reference=reference,centre_air_A=centre,leg=leg,
                            historical=old_record.__dict__,current=new_record.__dict__,historical_area_mA=old_area,current_area_mA=new_area,
                            exact_record_and_continuum=exact,area_delta_mA=None if old_area is None or new_area is None else new_area-old_area)
                        result['records'].append(item)
                        rows.append({k:v for k,v in item.items()if k not in ('historical','current')})
    result['exact_native_endpoint_replays']=replays
    args.output.mkdir(parents=True,exist_ok=False)
    (args.output/'report.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    fields=list(dict.fromkeys(key for row in rows for key in row))
    with (args.output/'summary.csv').open('x')as stream:
        writer=csv.DictWriter(stream,fieldnames=fields,lineterminator='\n');writer.writeheader();writer.writerows(rows)
    print(json.dumps(dict(records=len(rows),exact_native_endpoint_replays=replays,exact_current_matches=sum(r['exact_record_and_continuum']for r in result['records']),output=str(args.output))))


if __name__=='__main__':main()
