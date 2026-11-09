"""Native-matching support replay for molecular LTE continuum sensitivities."""
import argparse,importlib.util,json,os,subprocess,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--analysis',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    from config.constants import require_codex_path
    from pipeline.bride_reference import fingerprint
    from pipeline.bride_continuum_audit import selection,trace_fit
    from pipeline.kp_atlas_integrity import read_kurucz2005_residual
    from pipeline.reflected_solar_rv import SOLAR_GRAVITATIONAL_REDSHIFT
    from scripts.audit_bride_kitt_residuals import integral
    os.environ['CODEX_CONT_REF']='iag'
    analysis=json.loads((a.analysis/'report.json').read_text());run=json.loads((a.source/'report.json').read_text())
    if fingerprint(a.source/'report.json')!=analysis['source']:raise ValueError('Changed native run')
    orchard=require_codex_path('experiments.bride_current_orch_checkout');path=orchard/'pipeline/local_continuum.py'
    if fingerprint(path)!=analysis['current_orch']:raise ValueError('Changed current native continuum')
    spec=importlib.util.spec_from_file_location('native_molecules_trace_current',path);lc=importlib.util.module_from_spec(spec);sys.modules[spec.name]=lc;spec.loader.exec_module(lc)
    sys.path.insert(0,str(require_codex_path('engines.ispec')));import ispec
    comparison_path=ROOT/'data/results/rya1281/solar_bound_v1/kitt_peak/report.json'
    if fingerprint(comparison_path)!=analysis['comparison']:raise ValueError('Changed comparison')
    b=json.loads(comparison_path.read_text())['broadening'];observation=require_codex_path('data.kurucz2005_residual')
    if fingerprint(observation)!=analysis['observation']:raise ValueError('Changed observation')
    w,obs=read_kurucz2005_residual(observation,7765,7782)
    result=dict(schema='codex.bride_native_molecule_trace/1',scientific_admission='HOLD',generator=fingerprint(__file__),
        analysis=fingerprint(a.analysis/'report.json'),native=fingerprint(path),adapter=fingerprint(ROOT/'pipeline/bride_continuum_audit.py'),
        scope='Exact record/continuum-array/area replay of actual current native estimator; historical trace bound inactive on all admitted records. Reconstructed coefficient supports, not postfit survivor counts.',records=[],changes=[])
    for case in run['cases']:
        flux_path=a.source/case['name']/'flux.csv'
        if fingerprint(flux_path)!=case['flux']:raise ValueError('Changed native spectrum')
        raw=np.genfromtxt(flux_path,delimiter=',',names=True);mw=raw['wavelength_air_A'];segment=np.array([(mw[0]/10,mw[-1]/10)],dtype=[('wave_base',float),('wave_top',float)])
        for rvlabel,rv in [('no_gravity',0.),('plus_gravity',SOLAR_GRAVITATIONAL_REDSHIFT)]:
            mf=ispec.apply_post_fundamental_effects(mw/10,raw['normalized_flux'].copy(),segment,
                macroturbulence=b['vmac_km_s'],vsini=b['vsini_km_s'],limb_darkening_coeff=b['limb_darkening'],R=b['R'],vrad=(rv,))
            for saved in [r for r in analysis['continuum']if r['case']==case['name']and r['rv']==rvlabel]:
                centre=saved['centre_air_A'];reference=saved['reference']
                rw,rf=(mw,mf)if reference=='synthesis'else lc.solar_atlas_reference(centre-5,centre+5,b['R'])
                keep=(rw>=centre-5.01)&(rw<=centre+5.01);ref=np.interp(w,rw[keep],rf[keep]);mask,_=selection(w,obs,centre,ref)
                native,continuum=lc.fit_model_guided(w,obs,centre,w,ref,exclude_half_width_A=.5,reference=reference)
                if native.__dict__!=saved['record']:raise ValueError('Current native record changed')
                if not native.applied or abs(native.level_at_centre-1)>.03:raise ValueError('Historical trace bound cannot replay this current-native case')
                trace,cont=trace_fit(w,obs,centre,ref,mask)
                if (trace['level'],trace['slope_per_A'],trace['reported_n_bins'],trace['applied'])!=(native.level_at_centre,native.slope_per_A,native.n_bins,native.applied):raise ValueError('Trace/current native mismatch')
                np.testing.assert_array_equal(cont,continuum)
                area=1000*integral(w,1-obs/cont,centre-.5,centre+.5)
                if area!=saved['observed_area_mA']:raise ValueError('Native area mismatch')
                trace['selection_indices']=np.flatnonzero(mask).tolist();trace['observed_area_mA']=area
                result['records'].append(dict(case=case['name'],rv=rvlabel,centre_air_A=centre,reference=reference,trace=trace))
    for item in result['records']:
        base_name='off_refined'if item['case'].endswith('_refined')else'off'
        before=next(r['trace']for r in result['records']if(r['case'],r['rv'],r['centre_air_A'],r['reference'])==(base_name,item['rv'],item['centre_air_A'],item['reference']))
        after=item['trace'];g0=before['geometry'];g1=after['geometry']
        result['changes'].append(dict(case=item['case'],rv=item['rv'],centre_air_A=item['centre_air_A'],reference=item['reference'],
            selection_difference=len(set(before['selection_indices'])^set(after['selection_indices'])),
            coefficient_support_difference=len(set(before['fit_indices'])^set(after['fit_indices'])),
            before_count=before['coefficient_support_count'],after_count=after['coefficient_support_count'],
            before_blue=g0['blue_count'],before_red=g0['red_count'],after_blue=g1['blue_count'],after_red=g1['red_count'],
            before_outside_hull_A=g0['centre_outside_hull_A'],after_outside_hull_A=g1['centre_outside_hull_A'],
            before_weight_l1=g0['centre_weight_l1'],after_weight_l1=g1['centre_weight_l1']))
    a.output.mkdir(parents=True,exist_ok=False);(a.output/'report.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(json.dumps(dict(exact_replays=len(result['records']))))


if __name__=='__main__':main()
