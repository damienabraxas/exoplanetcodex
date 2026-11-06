"""Compare frozen FeII controls to both atlases without gf/width/RV fitting."""
import argparse,csv,json,os,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    from config.constants import require_codex_path
    from pipeline.bride_reference import fingerprint
    from pipeline.bride_feii_profile import feature_metrics
    from pipeline.kp_atlas_integrity import read_kurucz2005_residual
    from pipeline.reflected_solar_rv import SOLAR_GRAVITATIONAL_REDSHIFT
    from pipeline.telluric_policy import gate_holding
    from pipeline.prenormalised_guard import assert_not_renormalising
    os.environ['CODEX_CONT_REF']='iag'
    from pipeline import local_continuum as lc
    from scripts.measure_band_ew import iag_atlas
    from astropy.io import fits
    from pipeline.uv_conditioning import vac_to_air
    sys.path.insert(0,str(require_codex_path('engines.ispec')));import ispec
    import matplotlib
    matplotlib.use('Agg');import matplotlib.pyplot as plt
    plan_path=ROOT/'config/experiments/rya1281_feii_profile.json';plan=json.loads(plan_path.read_text())
    if plan['schema']!='codex.bride_feii_profile_plan/1':raise ValueError('Unknown profile plan')
    transfer=require_codex_path('experiments.bride_feii_transfer')
    canonical=json.loads((transfer/'canonical_refined/report.json').read_text());centre=float(canonical['evidence']['canonical']['wavelength_air_A'])
    kw={k:plan[k] for k in ('core_half_A','shoulder_inner_A','shoulder_outer_A')}
    comparison=json.loads((ROOT/'data/results/rya1281/solar_bound_v1/kitt_peak/report.json').read_text());b=comparison['broadening']
    kp=require_codex_path('data.kurucz2005_residual');iag=require_codex_path('data.solar_iag_baker2020')
    if fingerprint(kp)!=comparison['observation']:raise ValueError('Changed Kitt holding')
    policies={}
    for holding,instrument in [('solar_kpno_kurucz2005_corrected','kpno_solar_atlas'),('solar_iag','iag_fts_solar_atlas')]:
        ok,why=gate_holding(holding,instrument)
        if not ok:raise ValueError(why)
        policies[holding]=why
    assert_not_renormalising('solar_kpno_kurucz2005_corrected',pre_normalised=True,local_window_envelope=True)
    half=plan['plot_half_A'];step=plan['grid_step_A'];grid=centre+np.arange(-half,half+step/2,step)
    kwave,kflux=read_kurucz2005_residual(kp,centre-half-.1,centre+half+.1)
    iw,iflux=iag_atlas();gw,gflux=lc.solar_atlas_reference(centre-half,centre+half,b['R'])
    for wave in (kwave,iw,gw):
        if wave[0]>grid[0] or wave[-1]<grid[-1]:raise ValueError('Incomplete feature coverage')
    data={'kitt_shipped':np.interp(grid,kwave,kflux),'iag_native':np.interp(grid,iw,iflux),'iag_shipped_gaussian':np.interp(grid,gw,gflux)}
    with fits.open(iag) as h:
        w=vac_to_air(1e8/np.asarray(h[1].data['v'],float));flags=np.asarray(h[1].data['flags'],float);keep=abs(w-centre)<=half
        if not keep.any() or np.any(flags[keep]!=0):raise ValueError('Nonzero/missing IAG feature quality flags')
        flag_count=int(keep.sum())
    result=dict(schema='codex.bride_feii_profile_comparison/1',scientific_admission='HOLD',generator=fingerprint(__file__),
        adapter=fingerprint(ROOT/'pipeline/bride_feii_profile.py'),plan=fingerprint(plan_path),kitt=fingerprint(kp),iag=fingerprint(iag),
        native_iag_reader=fingerprint(ROOT/'scripts/measure_band_ew.py'),reference_helper=fingerprint(lc.__file__),policies=policies,
        iag_feature_flag0_count=flag_count,interpolation_step_A=step,rest_air_A=centre,broadening=b,
        caveats='Raw shipped and fixed-shoulder shape views retained. Contrast is not an isolated EW or adopted continuum. No gf/width/RV/LSF fit, pixel covariance/noise significance, extra telluric correction or preferred atlas. IAG native and Gaussian legs are correlated views of one holding. Interpolation adds no resolution.',cases=[])
    rows=[];fig,axes=plt.subplots(2,2,figsize=(13,9))
    velocities={'canonical_broadened':0.,'plus_gravitational_shift':SOLAR_GRAVITATIONAL_REDSHIFT}
    for ci,name in enumerate(('canonical_refined','lower_gravity_refined')):
        directory=transfer/name;run=json.loads((directory/'report.json').read_text())
        if run['line']!=canonical['line']:raise ValueError('Changed line identity')
        models={};receipts={}
        for mode in ('zero','feii'):
            path=directory/mode/'flux.csv';receipt=json.loads((path.parent/'report.json').read_text())
            if fingerprint(path)!=receipt['flux']:raise ValueError('Changed frozen spectrum')
            raw=np.genfromtxt(path,delimiter=',',names=True);mw=raw['wavelength_air_A'];segment=np.array([(mw[0]/10,mw[-1]/10)],dtype=[('wave_base',float),('wave_top',float)])
            receipts[mode]=fingerprint(path);models[mode]={}
            for rv,velocity in velocities.items():
                flux=ispec.apply_post_fundamental_effects(mw/10,raw['normalized_flux'].copy(),segment,macroturbulence=b['vmac_km_s'],vsini=b['vsini_km_s'],limb_darkening_coeff=b['limb_darkening'],R=b['R'],vrad=(velocity,))
                models[mode][rv]=np.interp(grid,mw,flux)
        out=dict(name=name,source=fingerprint(directory/'report.json'),spectra=receipts,metrics=[])
        for rv,velocity in velocities.items():
            beta=velocity*1000/299792458.;c=centre*np.sqrt((1+beta)/(1-beta));core=abs(grid-c)<=plan['core_half_A']
            observed={label:feature_metrics(grid,f,c,**kw) for label,f in data.items()}
            for mode in models:
                f=models[mode][rv];m=feature_metrics(grid,f,c,**kw)
                for label,g in data.items():
                    d=observed[label]
                    row=dict(case=name,rv=rv,model=mode,observation=label,centre_air_A=c,model_contrast=m['contrast'],observed_contrast=d['contrast'],
                        contrast_difference=m['contrast']-d['contrast'],model_relative_core_area_mA=m['relative_core_area_mA'],observed_relative_core_area_mA=d['relative_core_area_mA'],
                        raw_core_rms=float(np.sqrt(np.mean((f[core]-g[core])**2))),
                        shoulder_scaled_core_rms=float(np.sqrt(np.mean((f[core]/m['shoulder_level']-g[core]/d['shoulder_level'])**2))))
                    rows.append(row);out['metrics'].append(row)
        c=centre*np.sqrt((1+SOLAR_GRAVITATIONAL_REDSHIFT*1000/299792458.)/(1-SOLAR_GRAVITATIONAL_REDSHIFT*1000/299792458.))
        for label,f in data.items():
            axes[0,ci].plot(grid-c,f,label=label,lw=1)
            m=feature_metrics(grid,f,c,**kw);axes[1,ci].plot(grid-c,f/m['shoulder_level'],label=label,lw=1)
        for mode in models:
            f=models[mode]['plus_gravitational_shift'];m=feature_metrics(grid,f,c,**kw)
            axes[0,ci].plot(grid-c,f,label='model_'+mode,lw=1.5)
            axes[1,ci].plot(grid-c,f/m['shoulder_level'],label='model_'+mode,lw=1.5)
        axes[0,ci].set_title(name+': shipped levels, fixed gravity convention')
        axes[1,ci].set_title('Fixed symmetric-shoulder shape diagnostic')
        for ax in axes[:,ci]:
            ax.axvspan(-plan['core_half_A'],plan['core_half_A'],color='grey',alpha=.1)
            for sign in (-1,1):ax.axvspan(sign*plan['shoulder_inner_A'],sign*plan['shoulder_outer_A'],color='green',alpha=.07)
            ax.set(xlabel='AIR A from fixed observer FeII centre',ylabel='Normalized flux / diagnostic ratio');ax.grid(alpha=.2);ax.legend(fontsize=7)
        result['cases'].append(out)
    fig.suptitle('FeII7777.099: no fitted gf/width/RV or adopted continuum; scientific HOLD');fig.tight_layout()
    args.output.mkdir(parents=True,exist_ok=False);fig.savefig(args.output/'profile_comparison.png',dpi=140);plt.close(fig)
    with (args.output/'summary.csv').open('x') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]),lineterminator='\n');writer.writeheader();writer.writerows(rows)
    (args.output/'report.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps([r for r in rows if r['rv']=='plus_gravitational_shift' and r['model']=='feii']))


if __name__=='__main__':main()
