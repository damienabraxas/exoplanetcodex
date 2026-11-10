"""Paired H2 collider-only FeI pressure control; no fit or full-EOS claim."""
import argparse,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--runs',type=Path,nargs='+',required=True)
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
    receipts=[json.loads((d/'report.json').read_text())for d in args.runs]
    for d,r in zip(args.runs,receipts):
        if r['schema']!='codex.bride_h2_collider_flux/1' or r['input_manifest']!=fingerprint(d/'input_manifest.json'):raise ValueError('Wrong/changed H2 control report input')
    identities=[(r['mode'],r['population_mode'])for r in receipts]
    expected=[('off','fixed'),('molecular','fixed')]
    if identities!=expected:raise ValueError('Require all declared controls in order; no missing/best-case filter')
    names=['H2_off','H2_molecular']
    for key in ('generator','stark_source','stark_kernel','pressure_source','pressure_provider','provider','line_coefficients','lines','partition','source_manifest','source_populations','binding','wide_baseline','chemistry_source','chemistry_adapter','collider_provider','stark_mode'):
        if receipts[0][key]!=receipts[1][key]:raise ValueError('Mixed control inputs: '+key)
    arrays=[]
    for d,r in zip(args.runs,receipts):
        if fingerprint(d/'flux.csv')!=r['flux']:raise ValueError('Changed control flux')
        arrays.append(np.genfromtxt(d/'flux.csv',delimiter=',',names=True))
    wave=arrays[0]['wavelength_air_A']
    for array in arrays[1:]:np.testing.assert_array_equal(wave,array['wavelength_air_A'])
    for array in arrays[1:]:np.testing.assert_array_equal(arrays[0]['continuum_Fnu_W_m2_Hz'],array['continuum_Fnu_W_m2_Hz'])
    previous=ROOT/'data/results/rya1281/solar_bound_v1/kitt_peak/report.json';prior=json.loads(previous.read_text())
    allowed,reason=gate_holding('solar_kpno_kurucz2005_corrected','kpno_solar_atlas')
    if not allowed:raise ValueError(reason)
    path=require_codex_path('data.kurucz2005_residual')
    if fingerprint(path)!=prior['observation']:raise ValueError('Changed Kitt observation')
    ow,obs=read_kurucz2005_residual(path,7763,7783)
    b=prior['broadening'];segments=np.array([(wave[0]/10,wave[-1]/10)],dtype=[('wave_base',float),('wave_top',float)])
    result=dict(schema='codex.bride_h2_collider_review/1',scientific_admission='HOLD',generator=fingerprint(__file__),
        controls=[fingerprint(d/'report.json')for d in args.runs],observation=fingerprint(path),policy=reason,
        broadening=b,ispec_effects=fingerprint(inspect.getsourcefile(ispec.apply_post_fundamental_effects)),cases=[])
    profiles=[wave,np.interp(wave,ow,obs)];labels=['wavelength_air_A','shipped_observation']
    args.output.mkdir(parents=True,exist_ok=False);fig,axes=plt.subplots(3,1,figsize=(12,9))
    for ax in axes:ax.plot(ow,obs,color='black',lw=.9,label='Kitt shipped continuum')
    for name,data in zip(names,arrays):
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
                for ax in axes:ax.plot(wave,f,lw=.9,label=name.replace('_',' '))
        result['cases'].append(case)
    widths=[]
    for directory,receipt in zip(args.runs,receipts):
        if fingerprint(directory/'widths.csv')!=receipt['widths']:raise ValueError('Changed width decomposition')
        widths.append(np.genfromtxt(directory/'widths.csv',delimiter=',',names=True,dtype=None,encoding=None))
    for key in ('line_id','depth_index','gamma_rad_s_inv','gamma_native_HHe_s_inv','gamma_Stark_s_inv'):
        np.testing.assert_array_equal(widths[0][key],widths[1][key])
    w=widths[1];den=w['gamma_native_HHe_s_inv']
    net=(w['gamma_H2_added_s_inv']+w['gamma_HI_depletion_s_inv'])/den
    result['width_decomposition']=dict(
        H2_added_relative_native_HHe_range=[float(min(w['gamma_H2_added_s_inv']/den)),float(max(w['gamma_H2_added_s_inv']/den))],
        HI_depletion_relative_native_HHe_range=[float(min(w['gamma_HI_depletion_s_inv']/den)),float(max(w['gamma_HI_depletion_s_inv']/den))],
        net_relative_native_HHe_range=[float(min(net)),float(max(net))])
    result['oxygen_aperture_max_flux_delta']=[float(max(abs(arrays[1]['normalized_flux'][abs(wave-c)<=.5]-arrays[0]['normalized_flux'][abs(wave-c)<=.5])))for c in (7771.944,7774.166,7775.388)]
    result['area_changes_vs_H2_off_mA']={c['mode']:[v['area_mA']-base['area_mA']for v,base in zip(c['intrinsic'],result['cases'][0]['intrinsic'])]for c in result['cases']}
    result['max_profile_change_vs_H2_off']={name:float(np.max(abs(array['normalized_flux']-arrays[0]['normalized_flux'])))for name,array in zip(names,arrays)}
    ranges=[(7770.5,7776.5), (7777.8,7781.2), (7766.9,7770.7)]
    for ax,limits in zip(axes,ranges):
        ax.set_xlim(*limits)
        inside=(wave>=limits[0])&(wave<=limits[1])
        floor=min(float(np.min(a[inside]))for a in profiles[1:])
        ax.set_ylim(max(0,floor-.02),1.005)
        ax.legend(fontsize=8);ax.grid(alpha=.2);ax.set_ylabel('Normalized flux')
    axes[-1].set_xlabel('AIR wavelength (Angstrom)')
    fig.suptitle('41 source-matched Fe I lines: paired H2 pressure colliders\nFixed O/native atomic continuum and collisions; depleted HI plus H2, no fit')
    fig.tight_layout();figure=args.output/'h2_collider_control.png';fig.savefig(figure,dpi=150);plt.close(fig)
    fig,axes=plt.subplots(3,1,figsize=(10,8))
    delta=arrays[1]['normalized_flux']-arrays[0]['normalized_flux']
    for ax,c in zip(axes,(7771.944,7774.166,7775.388)):
        ax.plot(wave,delta*1e6);ax.axhline(0,color='black',lw=.5);ax.set_xlim(c-.5,c+.5)
        inside=(wave>=c-.5)&(wave<=c+.5);limit=float(np.max(abs(delta[inside])))*1.1e6
        ax.set_ylim(-limit,limit);ax.ticklabel_format(axis='x',useOffset=False)
        ax.set_ylabel('Flux change (ppm)');ax.grid(alpha=.2)
    axes[-1].set_xlabel('AIR wavelength (Angstrom)')
    fig.suptitle('H2 collider-only intervention minus atomic control; fixed O populations')
    fig.tight_layout();difference_figure=args.output/'h2_collider_difference.png';fig.savefig(difference_figure,dpi=150);plt.close(fig)
    result['difference_figure']=fingerprint(difference_figure)
    path=args.output/'profiles.csv';np.savetxt(path,np.column_stack(profiles),delimiter=',',header=','.join(labels),comments='')
    result.update(profiles=fingerprint(path),figure=fingerprint(figure),caveats='FeI pressure collider intervention only, depleted HI plus H2, fixed Cowley Stark/O populations/native continuum and O collisions. No full EOS, molecular opacity or abundance admission. Both RV variants retained; no optimization.')
    (args.output/'report.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print('Area changes vs H2 off',result['area_changes_vs_H2_off_mA']);print('Max flux change',result['max_profile_change_vs_H2_off'])


if __name__=='__main__':main()
