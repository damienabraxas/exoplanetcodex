"""Cross immutable FeII width references and actual continuum-fit supports."""
import argparse,csv,json,os,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    from config.constants import require_codex_path
    from pipeline.bride_reference import fingerprint
    from pipeline.bride_crossed_continuum import crossed_fits,finite_accounting
    from pipeline.bride_continuum_audit import selection
    from pipeline.kp_atlas_integrity import read_kurucz2005_residual
    from pipeline.reflected_solar_rv import SOLAR_GRAVITATIONAL_REDSHIFT
    from scripts.audit_bride_kitt_residuals import integral
    os.environ['CODEX_CONT_REF']='iag'
    from pipeline import local_continuum as lc
    sys.path.insert(0,str(require_codex_path('engines.ispec')));import ispec
    plan_path=ROOT/'config/experiments/rya1281_crossed_continuum.json';plan=json.loads(plan_path.read_text())
    if plan['schema']!='codex.bride_crossed_continuum_plan/1':raise ValueError('Wrong crossed plan')
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
    result=dict(schema='codex.bride_crossed_continuum/1',scientific_admission='HOLD',generator=fingerprint(__file__),
        kernel=fingerprint(ROOT/'pipeline/bride_crossed_continuum.py'),trace_kernel=fingerprint(ROOT/'pipeline/bride_continuum_audit.py'),
        plan=fingerprint(plan_path),comparison=fingerprint(comparison_path),native=fingerprint(lc.__file__),
        observation=fingerprint(observation),iag=fingerprint(iag),iag_reader=fingerprint(ROOT/'scripts/measure_band_ew.py'),
        scope=plan['scope'],sources=[],records=[])
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
                    for design in plan['designs']:
                        masks=coefficients if design=='coefficient_support_no_clipping' else selected
                        fits=crossed_fits(w,obs,centre,refs,masks,clip=design=='selected_with_clipping')
                        corners={}
                        for corner,(trace,continuum) in fits.items():
                            area=None if continuum is None else 1000*integral(w,1-obs/continuum,centre-.5,centre+.5)
                            item=dict(applied=trace['applied'],reason=trace['reason'],level=trace.get('level'),slope_per_A=trace.get('slope_per_A'),
                                observed_area_mA=area,fit_indices=trace.get('fit_indices'),count=trace['coefficient_support_count'],
                                stop_reason=trace['clipping_stop_reason'],passes=trace['clipping_pass_count'])
                            corners[corner]=item
                            if corner in ('before_before','after_after') and design!='selected_no_clipping':
                                identity=corner.split('_')[0];expected=native[identity]
                                if (item['level'],item['slope_per_A'],area,item['fit_indices'])!=(expected['level'],expected['slope_per_A'],expected.get('observed_area_mA'),expected['fit_indices']):
                                    raise ValueError('Native coefficient/area endpoint replay differs')
                                replays+=1
                        accounting={key:finite_accounting(corners,key) for key in ('level','observed_area_mA')}
                        result['records'].append(dict(stark_mode=mode,case=name,reference=reference,centre_air_A=centre,design=design,corners=corners,accounting=accounting))
                        rows.append(dict(stark_mode=mode,case=name,reference=reference,centre_air_A=centre,design=design,
                            **{key:value for key,value in accounting['observed_area_mA'].items() if key not in ('available','reason')}))
    result['exact_native_endpoint_replays']=replays
    args.output.mkdir(parents=True,exist_ok=False)
    (args.output/'report.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    fields=list(dict.fromkeys(key for row in rows for key in row))
    with (args.output/'summary.csv').open('x')as stream:
        writer=csv.DictWriter(stream,fieldnames=fields,lineterminator='\n');writer.writeheader();writer.writerows(rows)
    print(json.dumps(dict(records=len(rows),exact_native_endpoint_replays=replays,output=str(args.output))))


if __name__=='__main__':main()
