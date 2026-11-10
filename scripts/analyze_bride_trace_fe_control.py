"""Paired fixed-population Fe I trace-line control; unchanged Kitt AIR comparison."""
import argparse,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--none',type=Path,required=True);p.add_argument('--natural',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    from config.constants import require_codex_path
    from pipeline.bride_reference import fingerprint
    from pipeline.bride_observables import window_integral
    from pipeline.kp_atlas_integrity import read_kurucz2005_residual
    from pipeline.telluric_policy import gate_holding
    from pipeline.reflected_solar_rv import SOLAR_GRAVITATIONAL_REDSHIFT
    sys.path.insert(0,str(require_codex_path('engines.ispec')));import ispec
    import inspect,matplotlib
    matplotlib.use('Agg');import matplotlib.pyplot as plt
    receipts=[json.loads((d/'report.json').read_text())for d in (args.none,args.natural)]
    if [r['mode']for r in receipts]!=['none','natural']:raise ValueError('Wrong paired modes')
    for key in ('generator','provider','line_coefficients','lines','partition','source_manifest','source_populations','binding','wide_baseline'):
        if receipts[0][key]!=receipts[1][key]:raise ValueError('Mixed control inputs: '+key)
    arrays=[]
    for d,r in zip((args.none,args.natural),receipts):
        if fingerprint(d/'flux.csv')!=r['flux']:raise ValueError('Changed control flux')
        arrays.append(np.genfromtxt(d/'flux.csv',delimiter=',',names=True))
    wave=arrays[0]['wavelength_air_A'];np.testing.assert_array_equal(wave,arrays[1]['wavelength_air_A'])
    np.testing.assert_array_equal(arrays[0]['continuum_Fnu_W_m2_Hz'],arrays[1]['continuum_Fnu_W_m2_Hz'])
    previous=ROOT/'data/results/rya1281/solar_bound_v1/kitt_peak/report.json';prior=json.loads(previous.read_text())
    allowed,reason=gate_holding('solar_kpno_kurucz2005_corrected','kpno_solar_atlas')
    if not allowed:raise ValueError(reason)
    path=require_codex_path('data.kurucz2005_residual')
    if fingerprint(path)!=prior['observation']:raise ValueError('Changed Kitt observation')
    ow,obs=read_kurucz2005_residual(path,7763,7783)
    b=prior['broadening'];segments=np.array([(wave[0]/10,wave[-1]/10)],dtype=[('wave_base',float),('wave_top',float)])
    result=dict(schema='codex.bride_trace_fe_control_review/1',scientific_admission='HOLD',generator=fingerprint(__file__),
        controls=[fingerprint(d/'report.json')for d in (args.none,args.natural)],observation=fingerprint(path),policy=reason,
        broadening=b,ispec_effects=fingerprint(inspect.getsourcefile(ispec.apply_post_fundamental_effects)),cases=[])
    profiles=[wave,np.interp(wave,ow,obs)];labels=['wavelength_air_A','shipped_observation']
    args.output.mkdir(parents=True,exist_ok=False);fig,axes=plt.subplots(3,1,figsize=(12,9))
    for ax in axes:ax.plot(ow,obs,color='black',lw=.9,label='Kitt shipped continuum')
    for name,data in zip(('none','natural'),arrays):
        case=dict(mode=name,intrinsic=[],variants=[])
        for c in (7771.944,7774.166,7775.388):case['intrinsic'].append(dict(centre_air_A=c,area_mA=1000*window_integral(wave,data['normalized_flux'],c-.5,c+.5)))
        for label,rv in [('canonical_broadened',0.),('plus_gravitational_shift',SOLAR_GRAVITATIONAL_REDSHIFT)]:
            f=ispec.apply_post_fundamental_effects(wave/10,data['normalized_flux'].copy(),segments,
                macroturbulence=b['vmac_km_s'],vsini=b['vsini_km_s'],limb_darkening_coeff=b['limb_darkening'],R=b['R'],vrad=(rv,))
            profiles.append(f);labels.append(name+'_'+label)
            metrics=[]
            for c in (7771.944,7774.166,7775.388):metrics.append(dict(centre_air_A=c,model_area_mA=1000*window_integral(wave,f,c-.5,c+.5),observed_area_mA=1000*window_integral(ow,obs,c-.5,c+.5)))
            case['variants'].append(dict(name=label,rv_km_s=rv,metrics=metrics))
            if rv!=0:
                for ax in axes:ax.plot(wave,f,lw=.9,label=('O only' if name=='none' else 'O + 41 Fe I natural-Voigt lines'))
        result['cases'].append(case)
    result['added_intrinsic_areas_mA']=[n['area_mA']-z['area_mA']for z,n in zip(result['cases'][0]['intrinsic'],result['cases'][1]['intrinsic'])]
    result['max_intrinsic_profile_difference']=float(np.max(abs(arrays[1]['normalized_flux']-arrays[0]['normalized_flux'])))
    ranges=[(7770.5,7776.5), (7777.8,7781.2), (7766.9,7770.7)]
    for ax,limits in zip(axes,ranges):
        ax.set_xlim(*limits)
        inside=(wave>=limits[0])&(wave<=limits[1])
        floor=min(float(np.min(a[inside]))for a in profiles[1:])
        ax.set_ylim(max(0,floor-.02),1.005)
        ax.legend(fontsize=8);ax.grid(alpha=.2);ax.set_ylabel('Normalized flux')
    axes[-1].set_xlabel('AIR wavelength (Angstrom)')
    fig.suptitle('Fixed O populations + 41 source-matched natural-Voigt Fe I lines\nIncomplete trace control; no pressure/molecular opacity or fit; canonical gravitational shift')
    fig.tight_layout();figure=args.output/'trace_fe_control.png';fig.savefig(figure,dpi=150);plt.close(fig)
    path=args.output/'profiles.csv';np.savetxt(path,np.column_stack(profiles),delimiter=',',header=','.join(labels),comments='')
    result.update(profiles=fingerprint(path),figure=fingerprint(figure),caveats='FeI subset, natural Voigt only; distinct Kurucz LTE neutral reservoir reported. No complete blend/continuum/abundance admission. Both RV variants retained; no optimization.')
    (args.output/'report.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print('Added intrinsic window area mA',result['added_intrinsic_areas_mA']);print('Max profile delta',result['max_intrinsic_profile_difference'])


if __name__=='__main__':main()
