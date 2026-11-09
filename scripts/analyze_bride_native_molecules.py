"""Native LTE molecular-file sensitivities and current-ORCH reference diagnostics."""
import argparse,csv,importlib.util,json,os,subprocess,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    from config.constants import require_codex_path
    from pipeline.bride_reference import fingerprint
    from pipeline.kp_atlas_integrity import read_kurucz2005_residual
    from pipeline.reflected_solar_rv import SOLAR_GRAVITATIONAL_REDSHIFT
    from scripts.audit_bride_kitt_residuals import integral
    os.environ['CODEX_CONT_REF']='iag'
    run=json.loads((a.source/'report.json').read_text());plan=json.loads((ROOT/'config/experiments/rya1281_native_molecules.json').read_text())
    if run['plan']!=fingerprint(ROOT/'config/experiments/rya1281_native_molecules.json')or run['current_orch_commit']!=plan['orch_commit']:raise ValueError('Changed native plan')
    orchard=require_codex_path('experiments.bride_current_orch_checkout');native=orchard/'pipeline/local_continuum.py'
    if native.read_bytes()!=subprocess.check_output(['git','show',plan['orch_commit']+':pipeline/local_continuum.py'],cwd=orchard):raise ValueError('Changed native ORCH continuum')
    spec=importlib.util.spec_from_file_location('molecular_current_continuum',native);lc=importlib.util.module_from_spec(spec);sys.modules[spec.name]=lc;spec.loader.exec_module(lc)
    if (lc.MAX_SHIFT_MODEL,lc.MAX_SHIFT,lc.MODEL_Q)!=(.15,.03,90.):raise ValueError('Changed continuum settings')
    sys.path.insert(0,str(require_codex_path('engines.ispec')));import ispec
    comparison_path=ROOT/'data/results/rya1281/solar_bound_v1/kitt_peak/report.json';comparison=json.loads(comparison_path.read_text());b=comparison['broadening']
    observation=require_codex_path('data.kurucz2005_residual')
    if fingerprint(observation)!=comparison['observation']:raise ValueError('Changed Kitt product')
    ow,obs=read_kurucz2005_residual(observation,7765,7782)
    profile=json.loads((ROOT/'data/results/rya1281/feii_profile_source_v1/profile/report.json').read_text())
    if fingerprint(require_codex_path('data.solar_iag_baker2020'))!=profile['iag']:raise ValueError('Changed IAG source')
    result=dict(schema='codex.bride_native_molecule_analysis/1',scientific_admission='HOLD',generator=fingerprint(__file__),
        source=fingerprint(a.source/'report.json'),current_orch=fingerprint(native),comparison=fingerprint(comparison_path),
        observation=fingerprint(observation),iag=profile['iag'],intrinsic=[],refinement=[],continuum=[],
        limitations='Legacy-source LTE native5750/4.5, not canonical solar5772/4.438 or Bride NLTE correction. EPS1e-4 and native5decimal output retained; correlated finite window/reference diagnostics, no abundance or adopted continuum.')
    centres=(7771.944,7774.166,7775.388);spectra={};areas={}
    for case in run['cases']:
        path=a.source/case['name']/'flux.csv'
        if fingerprint(path)!=case['flux']:raise ValueError('Changed native flux')
        raw=np.genfromtxt(path,delimiter=',',names=True);spectra[case['name']]=raw
        areas[case['name']]={str(c):1000*integral(raw['wavelength_air_A'],1-raw['normalized_flux'],c-.5,c+.5)for c in centres}
    for case,raw in spectra.items():
        baseline='off_refined'if case.endswith('_refined')else'off';base=spectra[baseline]
        np.testing.assert_array_equal(raw['wavelength_air_A'],base['wavelength_air_A'])
        for centre in centres:
            # Two independently rounded profiles: max paired 1A area error0.01mA.
            result['intrinsic'].append(dict(case=case,centre_air_A=centre,area_mA=areas[case][str(centre)],
                off_area_mA=areas[baseline][str(centre)],delta_mA=areas[case][str(centre)]-areas[baseline][str(centre)],
                paired_native_rounding_area_bound_mA=.01,
                max_full_region_flux_delta=float(np.max(abs(raw['normalized_flux']-base['normalized_flux'])))))
    for centre in centres:
        coarse=areas['all'][str(centre)]-areas['off'][str(centre)];fine=areas['all_refined'][str(centre)]-areas['off_refined'][str(centre)]
        summed=sum(areas[k][str(centre)]-areas['off'][str(centre)]for k in ('CH','CN','C2','NH','OH'))
        result['refinement'].append(dict(centre_air_A=centre,coarse_delta_mA=coarse,fine_delta_mA=fine,difference_mA=fine-coarse,
            sum_family_deltas_mA=summed,all_minus_family_sum_mA=coarse-summed,
            caveat='Overlap and native cutoff/rounding can contribute; family effects are not independent uncertainty terms.'))
    for name,raw in spectra.items():
        wave=raw['wavelength_air_A'];segment=np.array([(wave[0]/10,wave[-1]/10)],dtype=[('wave_base',float),('wave_top',float)])
        for rvlabel,rv in [('no_gravity',0.),('plus_gravity',SOLAR_GRAVITATIONAL_REDSHIFT)]:
            flux=ispec.apply_post_fundamental_effects(wave/10,raw['normalized_flux'].copy(),segment,
                macroturbulence=b['vmac_km_s'],vsini=b['vsini_km_s'],limb_darkening_coeff=b['limb_darkening'],R=b['R'],vrad=(rv,))
            for centre in centres:
                model_area=1000*integral(wave,1-flux,centre-.5,centre+.5)
                for reference in ('synthesis','iag'):
                    rw,rf=(wave,flux)if reference=='synthesis'else lc.solar_atlas_reference(centre-5,centre+5,b['R'])
                    if reference=='iag':
                        reader=sys.modules['measure_band_ew']
                        if reader.IAG_FITS.resolve()!=require_codex_path('data.solar_iag_baker2020').resolve():raise ValueError('Wrong actual IAG helper input')
                        result['iag_reader']=fingerprint(reader.__file__)
                    keep=(rw>=centre-5.01)&(rw<=centre+5.01);ref=np.interp(ow,rw[keep],rf[keep])
                    rec,cont=lc.fit_model_guided(ow,obs,centre,ow,ref,exclude_half_width_A=.5,reference=reference)
                    area=None if cont is None else 1000*integral(ow,1-obs/cont,centre-.5,centre+.5)
                    result['continuum'].append(dict(case=name,rv=rvlabel,centre_air_A=centre,reference=reference,record=rec.__dict__,
                        observed_area_mA=area,model_area_mA=model_area))
    a.output.mkdir(parents=True,exist_ok=False)
    with (a.output/'intrinsic.csv').open('x')as stream:
        w=csv.DictWriter(stream,fieldnames=list(result['intrinsic'][0]),lineterminator='\n');w.writeheader();w.writerows(result['intrinsic'])
    rows=[]
    for item in result['continuum']:
        base_name='off_refined'if item['case'].endswith('_refined')else'off'
        base=next(r for r in result['continuum']if r['case']==base_name and (r['rv'],r['centre_air_A'],r['reference'])==(item['rv'],item['centre_air_A'],item['reference']))
        rows.append(dict(case=item['case'],rv=item['rv'],centre_air_A=item['centre_air_A'],reference=item['reference'],
            applied=item['record']['applied'],level=item['record']['level_at_centre'],observed_area_mA=item['observed_area_mA'],
            observed_delta_mA=None if item['observed_area_mA']is None or base['observed_area_mA']is None else item['observed_area_mA']-base['observed_area_mA'],
            model_area_mA=item['model_area_mA'],model_delta_mA=item['model_area_mA']-base['model_area_mA']))
    with (a.output/'continuum.csv').open('x')as stream:
        w=csv.DictWriter(stream,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
    (a.output/'report.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(intrinsic=len(result['intrinsic']),continuum=len(rows),refinement=result['refinement'])))


if __name__=='__main__':main()
