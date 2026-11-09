"""Matched frozen/coupled spectra and actual current-ORCH continuum, both RVs."""
import argparse,csv,json,os,subprocess,sys,importlib.util
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    from config.constants import require_codex_path
    from pipeline.bride_reference import fingerprint
    from pipeline.kp_atlas_integrity import read_kurucz2005_residual
    from pipeline.reflected_solar_rv import SOLAR_GRAVITATIONAL_REDSHIFT
    from scripts.audit_bride_kitt_residuals import integral
    os.environ['CODEX_CONT_REF']='iag'
    orchard=require_codex_path('experiments.bride_current_orch_checkout')
    pin=json.loads((ROOT/'config/experiments/rya1281_current_orch.json').read_text())['orch_commit']
    if subprocess.check_output(['git','rev-parse','HEAD'],cwd=orchard,text=True).strip()!=pin:raise ValueError('Changed current ORCH checkout')
    native=orchard/'pipeline/local_continuum.py'
    if native.read_bytes()!=subprocess.check_output(['git','show',pin+':pipeline/local_continuum.py'],cwd=orchard):raise ValueError('Dirty ORCH source')
    spec=importlib.util.spec_from_file_location('coupled_current_continuum',native);lc=importlib.util.module_from_spec(spec);sys.modules[spec.name]=lc;spec.loader.exec_module(lc)
    if (lc.MAX_SHIFT_MODEL,lc.MAX_SHIFT,lc.MODEL_Q)!=(.15,.03,90.):raise ValueError('Unexpected native continuum settings')
    sys.path.insert(0,str(require_codex_path('engines.ispec')));import ispec
    comparison_path=ROOT/'data/results/rya1281/solar_bound_v1/kitt_peak/report.json'
    comparison=json.loads(comparison_path.read_text());b=comparison['broadening']
    observation=require_codex_path('data.kurucz2005_residual')
    if fingerprint(observation)!=comparison['observation']:raise ValueError('Changed Kitt holding')
    w,obs=read_kurucz2005_residual(observation,7765,7782)
    profile=json.loads((ROOT/'data/results/rya1281/feii_profile_source_v1/profile/report.json').read_text())
    if fingerprint(require_codex_path('data.solar_iag_baker2020'))!=profile['iag']:raise ValueError('Changed IAG holding')
    result=dict(schema='codex.bride_feii_coupled_analysis/1',scientific_admission='HOLD',generator=fingerprint(__file__),
        orch_commit=pin,current_native=fingerprint(native),observation=fingerprint(observation),iag=profile['iag'],
        comparison=fingerprint(comparison_path),scope='Fresh O SE/radiation at fixed atmosphere/EOS; both fixed RVs/references. Diagnostic finite areas, not isolated EWs or adopted continuum/abundance.',cases=[])
    rows=[];pops={}
    for mode in ('natural','cowley','provided'):
        origin=args.source/mode;run=json.loads((origin/'report.json').read_text())
        if run['mode']!=mode or run['name']!='canonical_refined':raise ValueError('Wrong coupled case')
        if fingerprint(origin/'flux.csv')!=run['flux'] or fingerprint(origin/'populations.npz')!=run['coupled']['populations']:raise ValueError('Changed coupled flux/populations')
        if run['independent']['physical_rate_residual']>=run['controls']['rtol']:raise ValueError('Unaccepted physical balance')
        fixed=(require_codex_path('experiments.bride_feii_transfer')/'canonical_refined/feii') if mode=='natural' else (require_codex_path('experiments.bride_feii_width_transfer')/mode/'canonical_refined/feii')
        prior=json.loads((fixed/'report.json').read_text())
        if fingerprint(fixed/'report.json')!=run['frozen_comparison'] or fingerprint(fixed/'flux.csv')!=prior['flux']:raise ValueError('Changed matched frozen case')
        with np.load(origin/'populations.npz')as data:pops[mode]=data['populations'].copy();initial=data['initial_populations'].copy()
        if mode=='natural':initial_reference=initial
        else:np.testing.assert_array_equal(initial,initial_reference)
        spectra={key:np.genfromtxt(path,delimiter=',',names=True)for key,path in [('fixed',fixed/'flux.csv'),('coupled',origin/'flux.csv')]}
        np.testing.assert_array_equal(spectra['fixed']['wavelength_air_A'],spectra['coupled']['wavelength_air_A'])
        mw=spectra['fixed']['wavelength_air_A'];segment=np.array([(mw[0]/10,mw[-1]/10)],dtype=[('wave_base',float),('wave_top',float)])
        cases=[]
        for rvlabel,rv in [('no_gravity',0.),('plus_gravity',SOLAR_GRAVITATIONAL_REDSHIFT)]:
            models={key:ispec.apply_post_fundamental_effects(mw/10,raw['normalized_flux'].copy(),segment,
                macroturbulence=b['vmac_km_s'],vsini=b['vsini_km_s'],limb_darkening_coeff=b['limb_darkening'],R=b['R'],vrad=(rv,))for key,raw in spectra.items()}
            for centre in (7771.944,7774.166,7775.388):
                model_areas={key:1000*integral(mw,1-flux,centre-.5,centre+.5)for key,flux in models.items()}
                for reference in ('synthesis','iag'):
                    pair={}
                    for key,flux in models.items():
                        rw,rf=(mw,flux)if reference=='synthesis'else lc.solar_atlas_reference(centre-5,centre+5,b['R'])
                        if reference=='iag':
                            reader=sys.modules['measure_band_ew']
                            if reader.IAG_FITS.resolve()!=require_codex_path('data.solar_iag_baker2020').resolve():raise ValueError('Current helper read wrong IAG holding')
                            result['iag_reader']=fingerprint(reader.__file__)
                        keep=(rw>=centre-5.01)&(rw<=centre+5.01);ref=np.interp(w,rw[keep],rf[keep])
                        rec,cont=lc.fit_model_guided(w,obs,centre,w,ref,exclude_half_width_A=.5,reference=reference)
                        area=None if cont is None else 1000*integral(w,1-obs/cont,centre-.5,centre+.5)
                        pair[key]=dict(record=rec.__dict__,observed_area_mA=area)
                    delta=None if any(v['observed_area_mA']is None for v in pair.values())else pair['coupled']['observed_area_mA']-pair['fixed']['observed_area_mA']
                    item=dict(rv=rvlabel,centre_air_A=centre,reference=reference,pair=pair,model_areas_mA=model_areas,
                        model_area_delta_mA=model_areas['coupled']-model_areas['fixed'],observed_area_delta_mA=delta)
                    cases.append(item);rows.append(dict(mode=mode,rv=rvlabel,centre_air_A=centre,reference=reference,
                        model_area_delta_mA=item['model_area_delta_mA'],observed_area_delta_mA=delta,
                        fixed_observed_area_mA=pair['fixed']['observed_area_mA'],coupled_observed_area_mA=pair['coupled']['observed_area_mA']))
        result['cases'].append(dict(mode=mode,source=fingerprint(origin/'report.json'),cases=cases,
            max_relative_population_change_from_initial=float(np.max(abs(pops[mode]/initial-1))),independent=run['independent']))
    result['population_comparison_to_natural']={mode:float(np.max(abs(pop/pops['natural']-1)))for mode,pop in pops.items()}
    result['caution']='Pairwise changes below retained1e-6joint tolerance are not a physical no-feedback bound; no tighter convergence or abundance accuracy inferred.'
    args.output.mkdir(parents=True,exist_ok=False)
    (args.output/'report.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    with (args.output/'summary.csv').open('x')as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]),lineterminator='\n');writer.writeheader();writer.writerows(rows)
    print(json.dumps(dict(rows=len(rows),populations=result['population_comparison_to_natural'],max_model_area_delta_mA=max(abs(r['model_area_delta_mA'])for r in rows),max_observed_area_delta_mA=max(abs(r['observed_area_delta_mA'])for r in rows if r['observed_area_delta_mA']is not None))))


if __name__=='__main__':main()
