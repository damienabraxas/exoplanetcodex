"""Replay actual coefficient supports for the fixed-O FeII differential."""
import argparse,json,os,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--transfer',type=Path,required=True);p.add_argument('--analysis',type=Path,required=True);p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    from config.constants import require_codex_path
    from pipeline.bride_reference import fingerprint
    from pipeline.bride_continuum_audit import selection,trace_fit
    from pipeline.kp_atlas_integrity import read_kurucz2005_residual
    from pipeline.reflected_solar_rv import SOLAR_GRAVITATIONAL_REDSHIFT
    os.environ['CODEX_CONT_REF']='iag'
    from pipeline import local_continuum as lc
    from scripts.audit_bride_kitt_residuals import integral
    sys.path.insert(0,str(require_codex_path('engines.ispec')));import ispec
    comparison=json.loads((ROOT/'data/results/rya1281/solar_bound_v1/kitt_peak/report.json').read_text());b=comparison['broadening']
    observation=require_codex_path('data.kurucz2005_residual')
    if fingerprint(observation)!=comparison['observation']:raise ValueError('Changed observation')
    expected=json.loads((ROOT/'data/results/rya1281/continuum_support_v1/report.json').read_text())
    if fingerprint(lc.__file__)!=expected['native'] or lc.MAX_SHIFT!=.03:raise ValueError('Changed historical estimator')
    w,obs=read_kurucz2005_residual(observation,7765,7782)
    analysis=json.loads((args.analysis/'report.json').read_text())
    result=dict(schema='codex.bride_feii_anchor_trace/1',scientific_admission='HOLD',generator=fingerprint(__file__),
        analysis=fingerprint(args.analysis/'report.json'),native=fingerprint(lc.__file__),adapter=fingerprint(ROOT/'pipeline/bride_continuum_audit.py'),
        observation=fingerprint(observation),scope='Exact native-record/area replay, frozen historical3percent rule. No new continuum/adoption.',cases=[])
    for case in analysis['cases']:
        name=case['name'];pairs={}
        for mode in ('zero','feii'):
            path=args.transfer/name/mode/'flux.csv';receipt=json.loads((path.parent/'report.json').read_text())
            if fingerprint(path)!=receipt['flux']:raise ValueError('Changed spectrum')
            raw=np.genfromtxt(path,delimiter=',',names=True);mw=raw['wavelength_air_A']
            segment=np.array([(mw[0]/10,mw[-1]/10)],dtype=[('wave_base',float),('wave_top',float)])
            mf=ispec.apply_post_fundamental_effects(mw/10,raw['normalized_flux'].copy(),segment,
                macroturbulence=b['vmac_km_s'],vsini=b['vsini_km_s'],limb_darkening_coeff=b['limb_darkening'],R=b['R'],vrad=(SOLAR_GRAVITATIONAL_REDSHIFT,))
            pairs[mode]={}
            for reference in ('synthesis','iag'):
                apath=args.analysis/name/mode/reference/'report.json'
                if fingerprint(apath)!=case['audit_receipts'][mode][reference]:raise ValueError('Changed native audit')
                audit=json.loads(apath.read_text())
                for li,line in enumerate(audit['lines']):
                    centre=line['centre_air_A'];rw,rf=(mw,mf) if reference=='synthesis' else lc.solar_atlas_reference(centre-5,centre+5,b['R'])
                    keep=(rw>=centre-5.01)&(rw<=centre+5.01);ref=np.interp(w,rw[keep],rf[keep]);mask,_=selection(w,obs,centre,ref)
                    trace,cont=trace_fit(w,obs,centre,ref,mask);native=line['model_guided']
                    if (trace['level'],trace['slope_per_A'],trace['reported_n_bins'],trace['applied'])!=(native['level_at_centre'],native['slope_per_A'],native['n_bins'],native['applied']):raise ValueError('Native trace replay differs')
                    trace['selection_indices']=np.flatnonzero(mask).tolist()
                    if cont is not None:
                        area=1000*integral(w,1-obs/cont,centre-.5,centre+.5)
                        saved=next(r['observed_area_mA'] for r in line['variants'] if r['name']=='model_guided_'+reference+'_diagnostic')
                        if area!=saved:raise ValueError('Native area replay differs')
                        trace['observed_area_mA']=area
                    pairs[mode][str(centre)+'_'+reference]=trace
        changes=[]
        for key,before in pairs['zero'].items():
            after=pairs['feii'][key]
            changes.append(dict(key=key,selection_symmetric_difference=len(set(before['selection_indices'])^set(after['selection_indices'])),
                coefficient_support_symmetric_difference=len(set(before['fit_indices'])^set(after['fit_indices'])),
                before_actual_count=before['coefficient_support_count'],after_actual_count=after['coefficient_support_count'],
                before_blue=before['geometry']['blue_count'],before_red=before['geometry']['red_count'],
                after_blue=after['geometry']['blue_count'],after_red=after['geometry']['red_count']))
        result['cases'].append(dict(name=name,changes=changes,traces=pairs))
    args.output.mkdir(parents=True,exist_ok=False);(args.output/'report.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps([dict(name=c['name'],changes=c['changes']) for c in result['cases']]))


if __name__=='__main__':main()
