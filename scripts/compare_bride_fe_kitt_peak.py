"""Fixed Fe-background controls against unchanged canonical Kitt Peak AIR data.

Dedicated physics audit, not the fixed-physics parameter-grid admission path.
No abundance, continuum, RV, LSF or parameter optimization.
"""
import argparse
import inspect
import json
from pathlib import Path
import sys

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--receipt',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    from dataclasses import asdict
    import csv
    from config.constants import require_codex_path,get_star_params
    from pipeline.bride_reference import fingerprint
    from pipeline.bride_observables import window_integral
    from pipeline.kp_atlas_integrity import read_kurucz2005_residual
    from pipeline.telluric_policy import gate_holding,holding_basis,applied_state
    from pipeline.prenormalised_guard import assert_data_matches_declaration
    from pipeline.reflected_solar_rv import SOLAR_GRAVITATIONAL_REDSHIFT
    holding='solar_kpno_kurucz2005_corrected';instrument='kpno_solar_atlas'
    allowed,reason=gate_holding(holding,instrument)
    if not allowed:raise ValueError(reason)
    observation=require_codex_path('data.kurucz2005_residual')
    prior=json.loads((ROOT/'data/results/rya1281/kitt_peak_fixed_v1/manifest.json').read_text())
    if fingerprint(observation)!=prior['observed']:raise ValueError('Changed corrected observation')
    wave_obs,flux_obs=read_kurucz2005_residual(observation,7768,7780)
    evidence=assert_data_matches_declaration(holding,wave_obs,flux_obs,declared=True,where='Fe background physics diagnostic')
    sys.path.insert(0,str(require_codex_path('engines.ispec')));import ispec
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    solar=get_star_params('solar')
    with (ROOT/'data/catalog/instrument_catalog.csv').open()as stream:
        catalog=next(r for r in csv.DictReader(stream)if r['instrument_id']==instrument)
    R=float(catalog['resolving_power_max'])
    receipt=json.loads(args.receipt.read_text());store=require_codex_path('experiments.bride_mac_grid')
    result=dict(schema='codex.bride_fe_kitt_peak_comparison/1',scientific_admission='HOLD',
        generator=fingerprint(__file__),receipt=fingerprint(args.receipt),observation=fingerprint(observation),
        holding=holding,policy=dict(reason=reason,basis=holding_basis(holding,instrument),applied_state=applied_state(holding),
            normalization=asdict(evidence),action='Serve corrected pre-normalized observation unchanged'),
        frame='Observed delivered vacuum nm -> canonical AIR Angstrom; no extra observed shift/correction/continuum',
        broadening=dict(vmac_km_s=solar['vmac'],vsini_km_s=solar['vsini'],R=R,limb_darkening=.6,
            variants=['canonical_solar_broadening','plus_canonical_solar_gravitational_shift']),
        ispec_effects=fingerprint(inspect.getsourcefile(ispec.apply_post_fundamental_effects)),
        caveats='Native5750/4.5/A8.70 differs from canonical Sun; incomplete blends/background, approximate Gaussian catalog R instead of measured FTS LSF, no local error vector or chi-square/adoption. No optimizer or best-case ranking.',cases=[])
    args.output.mkdir(parents=True,exist_ok=False);fig,axes=plt.subplots(2,1,figsize=(12,7),sharex=True)
    view_obs=(wave_obs>=7770.5)&(wave_obs<=7776.8)
    axes[0].plot(wave_obs[view_obs],flux_obs[view_obs],color='black',lw=1,label='Kitt Peak corrected2005 unchanged')
    for entry in receipt['cases']:
        if entry['name']=='raw_pair_finite_column_cap':continue
        path=store/entry['source_run']/entry['name']
        if fingerprint(path/'manifest.json')!=entry['source_manifest'] or fingerprint(path/'flux.csv')!=entry['flux']:
            raise ValueError('Changed control input')
        manifest=json.loads((path/'manifest.json').read_text())
        if manifest['status']!='CONVERGED_EXPERIMENT':raise ValueError('Unconverged control')
        data=np.genfromtxt(path/'flux.csv',delimiter=',',names=True);wave=data['wavelength_air_A'];native=data['normalized_flux']
        segments=np.array([(wave[0]/10,wave[-1]/10)],dtype=[('wave_base',float),('wave_top',float)])
        case=dict(name=entry['name'],source_manifest=entry['source_manifest'],flux=entry['flux'],
            experimental_fe_background=manifest.get('experimental_fe_background'),variants=[])
        columns=[wave,native,np.interp(wave,wave_obs,flux_obs)]
        for label,rv in [('canonical_solar_broadening',0.),('plus_canonical_solar_gravitational_shift',SOLAR_GRAVITATIONAL_REDSHIFT)]:
            flux=ispec.apply_post_fundamental_effects(wave/10,native.copy(),segments,
                macroturbulence=solar['vmac'],vsini=solar['vsini'],limb_darkening_coeff=.6,R=R,vrad=(rv,))
            sample=[]
            for centre in (7771.944,7774.166,7775.388):
                mask=(wave>=centre-.5)&(wave<=centre+.5)
                sample.append(dict(centre_air_A=centre,model_area_mA=1000*window_integral(wave,flux,centre-.5,centre+.5),
                    observed_area_mA=1000*window_integral(wave_obs,flux_obs,centre-.5,centre+.5),
                    unweighted_RMS=float(np.sqrt(np.mean((flux[mask]-np.interp(wave[mask],wave_obs,flux_obs))**2)))))
            case['variants'].append(dict(name=label,rv_km_s=rv,metrics=sample));columns.append(flux)
            if rv==0.:
                view=(wave>=7770.5)&(wave<=7776.8)
                axes[0].plot(wave[view],flux[view],lw=.9,label=entry['name'])
                axes[1].plot(wave[view],flux[view]-np.interp(wave[view],wave_obs,flux_obs),lw=.9,label=entry['name'])
        output=args.output/(entry['name']+'_profiles.csv')
        np.savetxt(output,np.column_stack(columns),delimiter=',',header='wavelength_air_A,intrinsic,observation,canonical_broadened,plus_gravitational_shift',comments='')
        case['profiles']=fingerprint(output);result['cases'].append(case)
    axes[0].set_ylabel('Normalized disk flux');axes[1].set_ylabel('Model - observed')
    axes[1].set_xlabel('AIR wavelength (Angstrom)');axes[1].axhline(0,color='black',lw=.5)
    for ax in axes:ax.set_xlim(7770.5,7776.8);ax.legend(fontsize=7);ax.grid(alpha=.2)
    fig.suptitle('Fixed O8.70 / native5750,4.5 — Fe background diagnostics, no fit\nPlotted model RV=0; canonical gravitational-shift variant reported separately')
    fig.tight_layout();figure=args.output/'kitt_peak_fe_controls.png';fig.savefig(figure,dpi=160);plt.close(fig)
    result['figure']=fingerprint(figure)
    (args.output/'report.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print('Wrote unchanged-observation comparison:',args.output)


if __name__=='__main__':main()
