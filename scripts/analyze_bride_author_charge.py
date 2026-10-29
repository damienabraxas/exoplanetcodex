"""Compare sourced charge controls with frozen reference and prior paired case."""
import argparse
import json
from pathlib import Path
import sys

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference',type=Path,required=True)
    parser.add_argument('--runs',nargs='+',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--record-failed-controls',action='store_true',
                        help='Retain failed manifests/histories as failure receipts, without population curves')
    args=parser.parse_args()
    from config.constants import require_codex_path
    from pipeline.bride_reference import fingerprint,read_marcs_structure
    from pipeline.bride_observables import window_integral
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    if fingerprint(args.reference)['sha256']!='c3113efdb21e3f40d79c1836b0ddb0ebb41dbceecaf04d622af6b7a8d8c080ec':
        raise ValueError('Changed reference')
    reference=json.loads(args.reference.read_text())
    atmosphere=require_codex_path('grids.marcs_standard')/reference['atmosphere_file']['name']
    if fingerprint(atmosphere)!=reference['atmosphere_file']:raise ValueError('Changed atmosphere')
    structure=read_marcs_structure(atmosphere);tau=structure['logtau_500']
    target=np.array(reference['reference']['departures']);subset=(tau>=-3)&(tau<=0)
    args.output.mkdir(parents=True,exist_ok=False)
    result=dict(schema='codex.bride_author_charge_comparison/1',scientific_admission='HOLD',
        generator=fingerprint(__file__),reference=fingerprint(args.reference),cases=[],failed_controls=[],
        scope='logtau[-3,0] is a declared diagnostic subset, not a contribution function; no abundance fit')
    fig,axes=plt.subplots(1,3,figsize=(14,4));fluxfig,fluxax=plt.subplots(figsize=(10,4))
    for ax,level,label in zip(axes,(5,9,39),('Lower triplet','Upper triplet','O II ground')):
        ax.plot(tau,target[level],color='black',label='Independent reference')
        ax.set(xlabel='log tau500',ylabel='Departure b',title=label)
        if level==39:ax.set_yscale('log')
    oldreceipt=ROOT/'data/results/rya1281/native_charge_v1/analysis/summary.json'
    old=json.loads(oldreceipt.read_text())
    prior=next(x for x in old['cases']if x['name']=='native_tautsp_lte_balanced_CH0')
    store=require_codex_path('experiments.bride_mac_grid')
    priorpath=store/prior['source_run']/prior['name']
    if fingerprint(priorpath/'populations.npz')!=prior['populations']:
        raise ValueError('Changed prior paired populations')
    result['prior_receipt']=fingerprint(oldreceipt)
    previous=np.load(priorpath/'populations.npz');b=previous['populations']/previous['nstar']
    for ax,level in zip(axes,(5,9,39)):ax.plot(tau,b[level],ls='--',label='Prior raw CH0 LTE pair')
    for directory in args.runs:
        summary=json.loads((directory/'summary.json').read_text())
        if summary['status']!='EXPERIMENTS_FINISHED' or summary['atom']!=reference['atom']:
            raise ValueError('Changed/unfinished source run')
        for entry in summary['cases']:
            casepath=directory/entry['name'];manifest_path=casepath/'manifest.json'
            if fingerprint(manifest_path)!=entry['manifest']:raise ValueError('Changed manifest')
            manifest=json.loads(manifest_path.read_text());case=manifest['case']
            if manifest['status']=='FAILED' and args.record_failed_controls:
                result['failed_controls'].append(dict(name=entry['name'],source_run=directory.name,
                    source_summary=fingerprint(directory/'summary.json'),source_manifest=fingerprint(manifest_path),
                    history=fingerprint(casepath/'history.jsonl'),error=manifest['error'],
                    elapsed_seconds=manifest['elapsed_seconds']))
                continue
            if (manifest['status']!='CONVERGED_EXPERIMENT' or manifest['atmosphere']!=reference['atmosphere_file']
                    or (case['xi'],case['abundance'],case['drop_bottom'])!=(1.,8.7,0)
                    or manifest['experimental_collision_convention'] not in ('barklem2018_MN_LTE_pair','native_tautsp_lte_balanced_CH0')
                    or manifest['diagnostic_collision_scales']):
                raise ValueError('Unexpected native-node physics')
            data=np.load(casepath/'populations.npz')
            np.testing.assert_array_equal(data['nstar'],previous['nstar'])
            np.testing.assert_array_equal(data['temperature_K'],structure['temperature_K'])
            b=data['populations']/data['nstar'];delta=b-target
            fluxpath=casepath/'flux.csv'
            if fingerprint(fluxpath)!=manifest['flux']:raise ValueError('Changed flux')
            flux=np.genfromtxt(fluxpath,delimiter=',',names=True)
            curve=np.column_stack((tau,b[5],b[9],b[39],target[5],target[9],target[39]))
            curvepath=args.output/(entry['name']+'_departures.csv')
            np.savetxt(curvepath,curve,delimiter=',',header='logtau500,lower_b,upper_b,OII_b,reference_lower_b,reference_upper_b,reference_OII_b',comments='')
            sampled=int(np.argmin(abs(tau+1.2621)))
            cap_diagnostic=None
            if manifest.get('finite_column_top_cap'):
                cap_path=casepath/'top_cap.npz';cap=np.load(cap_path)
                length=float(cap['effective_length_m'])
                if length!=manifest['finite_column_top_cap']['effective_length_m']:
                    raise ValueError('Changed cap column metadata')
                selected=[]
                for wavelength in (121.57,130.20,135.56,277.2477,500.,777.4):
                    k=int(np.argmin(abs(cap['wavelength_vac_nm']-wavelength)))
                    selected.append(dict(wavelength_vac_nm=float(cap['wavelength_vac_nm'][k]),
                        vertical_optical_depth=float(cap['chi_m_inv'][k]*length)))
                cap_diagnostic=dict(snapshot=fingerprint(cap_path),selected_optical_depths=selected,
                                    gridpoint_flux=fingerprint(casepath/'gridpoint_flux.csv'))
            result['cases'].append(dict(name=entry['name'],source_run=directory.name,
                source_summary=fingerprint(directory/'summary.json'),source_manifest=fingerprint(manifest_path),
                populations=fingerprint(casepath/'populations.npz'),flux=fingerprint(fluxpath),
                departure_curve=fingerprint(curvepath),iterations=manifest['iterations'],
                elapsed_seconds=manifest['elapsed_seconds'],final_diagnostics=manifest['final_diagnostics'],
                upper_max_abs_b_difference_photospheric_subset=float(np.max(abs(delta[9,subset]))),
                upper_max_abs_b_difference_all_depths=float(np.max(abs(delta[9]))),
                upper_max_abs_b_difference_excluding_top_point=float(np.max(abs(delta[9,1:]))),
                lower_max_abs_b_difference_excluding_top_point=float(np.max(abs(delta[5,1:]))),
                all_levels_depths_max_abs_b_difference=float(np.max(abs(delta))),
                top_lower_upper_ion_b=b[[5,9,39],0].tolist(),
                sample_logtau500=float(tau[sampled]),sample_lower_upper_ion_b=b[[5,9,39],sampled].tolist(),
                triplet_fixed_window_mA=[float(window_integral(flux['wavelength_air_A'],flux['normalized_flux'],w-.5,w+.5)*1000)
                                          for w in (7771.94,7774.166,7775.388)],
                experimental_collision_convention=manifest['experimental_collision_convention'],
                author_charge_control=manifest.get('author_charge_control'),
                finite_column_top_cap=manifest.get('finite_column_top_cap'),
                top_cap_diagnostic=cap_diagnostic,
                lte_uv_background_lines=manifest.get('lte_uv_background_lines')))
            for ax,level in zip(axes,(5,9,39)):ax.plot(tau,b[level],label=entry['name'])
            fluxax.plot(flux['wavelength_air_A'],flux['normalized_flux'],label=entry['name'])
    for ax in axes:ax.legend(fontsize=6)
    fig.tight_layout();fig.savefig(args.output/'departures.png',dpi=150);plt.close(fig)
    fluxax.set(xlabel='Air wavelength (Angstrom)',ylabel='Intrinsic normalized flux')
    fluxax.legend(fontsize=7);fluxfig.tight_layout();fluxfig.savefig(args.output/'flux.png',dpi=150);plt.close(fluxfig)
    result['status']='DIAGNOSTIC_COMPLETE'
    (args.output/'summary.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print([(c['name'],c['upper_max_abs_b_difference_all_depths'])for c in result['cases']])


if __name__=='__main__':main()
