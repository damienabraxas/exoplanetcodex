"""Canonical Kitt Peak policy plus fixed iSpec broadening/input-response diagnostics.

Run in the project Python on the Mac (pandas/iSpec); RT/SE remains in the isolated
LW environment. No corrected observation is re-normalized or telluric-corrected.
"""
from __future__ import annotations
import argparse
from dataclasses import asdict
import inspect
from importlib.metadata import version
import json
from pathlib import Path
import platform
import subprocess
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--runs',type=Path,required=True)
    p.add_argument('--additional-runs',type=Path,nargs='*',default=[])
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    if platform.system()!='Darwin':raise RuntimeError('Mac-only analysis request')
    from config.constants import require_codex_path,get_star_params
    from pipeline.bride_reference import fingerprint
    from pipeline.bride_observables import window_integral
    from pipeline.kp_atlas_integrity import read_kurucz2005_residual
    from pipeline.telluric_policy import gate_holding,holding_basis,applied_state
    from pipeline.band_policy import check_intake
    from pipeline.prenormalised_guard import assert_data_matches_declaration
    from pipeline.reflected_solar_rv import SOLAR_GRAVITATIONAL_REDSHIFT
    combined=[];summaries=[];seen=set();common_physics=None
    for directory in [args.runs]+args.additional_runs:
        run=json.loads((directory/'summary.json').read_text())
        if run['status']!='EXPERIMENTS_FINISHED':raise ValueError('Experiments not finished')
        identity=(run['atom'],run['context_builder'])
        if common_physics is None:common_physics=identity
        elif identity!=common_physics:raise ValueError('Cannot combine different atoms or context physics into one grid')
        summaries.append(dict(directory_id=directory.name,summary=fingerprint(directory/'summary.json')))
        for entry in run['cases']:
            if entry['name']in seen:raise ValueError('Duplicate case identity across runs')
            seen.add(entry['name'])
            if fingerprint(directory/entry['name']/'manifest.json')!=entry['manifest']:
                raise ValueError('Changed case manifest')
            combined.append((directory,entry))
    holding='solar_kpno_kurucz2005_corrected';instrument='kpno_solar_atlas'
    permitted,reason=gate_holding(holding,instrument)
    if not permitted:raise ValueError(reason)
    band=check_intake(7774.166,'synthesis',instrument=instrument)
    observed_path=require_codex_path('data.kurucz2005_residual')
    prior=json.loads((ROOT/'data/results/rya1281/kitt_peak_fixed_v1/manifest.json').read_text())
    if fingerprint(observed_path)!=prior['observed']:raise ValueError('Corrected Kitt Peak product changed since baseline')
    with observed_path.open()as f:header=''.join(next(f)for _ in range(5))
    if 'vacuum wavelength including gravitational red shift'not in header:
        raise ValueError('Wrong source wavelength medium: expected delivered vacuum nm')
    wobs,fobs=read_kurucz2005_residual(observed_path,7768,7780)
    evidence=assert_data_matches_declaration(holding,wobs,fobs,declared=True,where='Bride Mac experiment')
    ispec_dir=require_codex_path('engines.ispec');sys.path.insert(0,str(ispec_dir))
    import ispec
    effects=Path(inspect.getsourcefile(ispec.apply_post_fundamental_effects))
    stars=get_star_params('solar')
    import csv
    with (ROOT/'data/catalog/instrument_catalog.csv').open()as f:
        inst=next(r for r in csv.DictReader(f)if r['instrument_id']==instrument)
    rlow=float(inst['resolving_power_min']);rhigh=float(inst['resolving_power_max'])
    macro=float(stars['vmac']);rotation=float(stars['vsini'])
    broad={
        'intrinsic':(0.,0.,0.,0.),
        'instrument_only':(0.,0.,rhigh,0.),
        'solar_broadening':(macro,rotation,rhigh,0.),
        'macro_minus1':(max(0.,macro-1),rotation,rhigh,0.),
        'macro_plus1':(macro+1,rotation,rhigh,0.),
        'rotation_minus1':(macro,max(0.,rotation-1),rhigh,0.),
        'rotation_plus1':(macro,rotation+1,rhigh,0.),
        'resolution_low':(macro,rotation,rlow,0.),
        'solar_gravitational_shift':(macro,rotation,rhigh,SOLAR_GRAVITATIONAL_REDSHIFT),
        'rv_minus_half':(macro,rotation,rhigh,-.5),
        'rv_plus_half':(macro,rotation,rhigh,.5)}
    args.output.mkdir(parents=True,exist_ok=False)
    result=dict(schema='codex.bride_mac_sensitivity_analysis/1',
        scientific_admission='HOLD: no best-fit selection, no calibrated solar abundance or full uncertainty budget',
        generator=fingerprint(__file__),code_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        input_summaries=summaries,host=run['host'],
        observation=fingerprint(observed_path),holding=holding,
        frame_contract=dict(observation_native='vacuum nm including solar gravitational redshift',
            comparison_and_plots='air Angstrom',engine_internal='vacuum nm',
            boundary='One canonical vac_to_air conversion on observed input; one air_to_vac conversion for the RT call',
            no_observed_velocity_shift_or_second_telluric_correction=True),
        orchestrator_policy=dict(gate_permitted=permitted,gate_reason=reason,holding_basis=holding_basis(holding,instrument),
            applied_state=applied_state(holding),band=asdict(band),normalization_evidence=asdict(evidence),
            action='Serve corrected pre-normalized product unchanged; no second correction or continuum placement'),
        ispec_effects=fingerprint(effects),ispec_signature=str(inspect.signature(ispec.apply_post_fundamental_effects)),
        solar_parameters=stars,instrument_resolution_interval=[rlow,rhigh],
        broadening_controls={k:dict(vmac_km_s=v[0],vsini_km_s=v[1],R=v[2],rv_km_s=v[3],limb_darkening=.6)for k,v in broad.items()},
        interpretation=['No observational continuum or parameter optimizer; all responses reported, not ranked as a best star model',
          'True Teff changes use independent MARCS structures; xi changes line/SE microturbulence at fixed t01 atmosphere',
          '5500/6000 K controls use common one-layer bottom cut, with separate 5750 K depth controls',
          'R interval is catalog capability, not a measured local FTS sinc LSF; iSpec Gaussian resolution is an approximation',
          'Solar broadening/rotation use canonical inputs; +/-1 km/s and RV +/-0.5 are diagnostic probes, not inferred uncertainties',
          'Published Kitt Peak continuum and corrected telluric product remain frozen; no local pixel telluric covariance is manufactured'],
        cases=[],failed_cases=[])
    result['capability_contract']=dict(geometry='1D plane-parallel',horizontal_radiative_transfer=False,
        full_3D=False,full_spectrum_validated=False,output_air_A=[7768.,7780.],
        observational_comparison='O I777nm triplet, normalized disk-integrated flux',
        microturbulence_scope='Line/SE microturbulence at fixed native MARCS t01 atmospheric structure',
        abundance_scope='Trace oxygen offsets; source background composition retained apart from explicit O override',
        background_scope='Native LW passive atoms/default continuum, no complete weak-blend/molecular spectrum',
        interpolation_scope='Native nodes and explicitly labeled experimental full-field interpolants; no fresh atmosphere/EOS equilibrium solve',
        raw_model_flux='Stellar surface Fnu; not absolute Earth irradiance. Distance dilution cancels only in normalized comparison.',
        physics_admission='HOLD: experimental native/interpolated atmospheres and full41-level atom; original MULTI discrepancy open')
    science_python=require_codex_path('engines.bride_eval')/'.venv/bin/python'
    query="import json,platform;from importlib.metadata import version;print(json.dumps(dict(python=platform.python_version(),packages={n:version(n) for n in ['lightweaver','numpy','scipy','numba','matplotlib']})))"
    result['science_environment']=json.loads(subprocess.check_output([str(science_python),'-c',query],text=True))
    result['analysis_environment']={n:version(n)for n in ('numpy','scipy','pandas','matplotlib')}
    centres=(7771.944,7774.166,7775.388)
    plotdata={};native_data={};basewave=None
    def metrics(wave,flux):
        values=[]
        for c in centres:
            lo,hi=c-.5,c+.5;mask=(wave>=lo)&(wave<=hi)
            ongrid=np.interp(wave[mask],wobs,fobs)
            values.append(dict(centre_air_A=c,
                model_area_mA=1000*window_integral(wave,flux,lo,hi),
                observed_area_mA=1000*window_integral(wobs,fobs,lo,hi),
                min_flux=float(flux[mask].min()),minimum_at_air_A=float(wave[mask][flux[mask].argmin()]),
                unweighted_RMS=float(np.sqrt(np.mean((flux[mask]-ongrid)**2)))))
        return values
    for directory,case in combined:
        path=directory/case['name'];m=json.loads((path/'manifest.json').read_text())
        if m.get('diagnostic_collision_scales') or m.get('experimental_collision_convention') or m.get('lte_h_background_lines') or m.get('lte_uv_background_lines'):
            raise ValueError('Diagnostic rate or background changes belong to their dedicated audit, not the fixed-atom parameter grid')
        if m['status']!='CONVERGED_EXPERIMENT':
            result['failed_cases'].append(dict(name=case['name'],reason=m.get('error','unknown')));continue
        if fingerprint(path/'flux.csv')!=m['flux']:raise ValueError('Changed case spectrum')
        d=np.genfromtxt(path/'flux.csv',delimiter=',',names=True)
        wave=d['wavelength_air_A'];native=d['normalized_flux'];segments=np.array([(wave[0]/10,wave[-1]/10)],dtype=[('wave_base',float),('wave_top',float)])
        native_data[case['name']]=(wave.copy(),native.copy())
        variants=broad if case['name']=='baseline' else {'solar_broadening':broad['solar_broadening'],'intrinsic':broad['intrinsic']}
        arrays=[];names=[];records={}
        for name,(vmac,vsini,R,rv)in variants.items():
            flux=native.copy()if name=='intrinsic'else ispec.apply_post_fundamental_effects(wave/10,native.copy(),segments,
                macroturbulence=vmac,vsini=vsini,limb_darkening_coeff=.6,R=R,vrad=(rv,))
            mask=(wave>=7770.5)&(wave<=7776.8)
            if not np.isfinite(flux[mask]).all()or np.any(flux[mask]<=0):raise ValueError('Invalid broadened analysis window')
            records[name]=metrics(wave,flux);arrays.append(flux);names.append(name)
            if name=='solar_broadening':plotdata[case['name']]=(wave.copy(),flux.copy())
        np.savetxt(args.output/(case['name']+'_profiles.csv'),np.column_stack([wave,native,np.interp(wave,wobs,fobs)]+arrays),
            delimiter=',',header=','.join(['wavelength_air_A','native_normalized_flux','observed_flux']+names),comments='')
        result['cases'].append(dict(name=case['name'],inputs=m['case'],source_manifest=fingerprint(path/'manifest.json'),
            source_directory_id=directory.name,atmosphere=m['atmosphere'],
            original_population_file=fingerprint(path/'populations.npz'),original_history=fingerprint(path/'history.jsonl'),
            final_diagnostics=m['final_diagnostics'],elapsed_seconds=m['elapsed_seconds'],variants=records))
    if 'baseline'not in plotdata:raise ValueError('Missing converged baseline')
    basewave,base=plotdata['baseline']
    case_names=list(plotdata)
    for name in case_names:
        if not np.array_equal(plotdata[name][0],basewave):raise ValueError('Case grid identity mismatch')
    np.savez_compressed(args.output/'profile_grid.npz',wavelength_air_A=basewave,
        case_names=np.asarray(case_names),native_normalized_flux=np.stack([native_data[n][1]for n in case_names]),
        solar_broadened_normalized_flux=np.stack([plotdata[n][1]for n in case_names]),
        observed_residual_flux=np.interp(basewave,wobs,fobs))
    # Broadening cannot create missing total absorption: also report a common wide
    # analysis-span integral, separate from the narrower per-component apertures.
    result['broad_solar_analysis_span_mA']=1000*window_integral(basewave,base,7770.5,7776.8)
    result['observed_analysis_span_mA']=1000*window_integral(wobs,fobs,7770.5,7776.8)
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    groups=[('Depth controls',['baseline','depth_cut1','depth_cut2']),
        ('Microturbulence (fixed atmosphere)',['baseline','xi_minus','xi_plus','xi_low','xi_high']),
        ('Teff (common bottom cut)',['depth_cut1','teff_low_cut1','teff_high_cut1']),
        ('Gravity',['baseline','logg_low','logg_high']),
        ('Metal mixture (common bottom cut; native alpha pattern)',['depth_cut1','metal_low_cut1','metal_high_cut1']),
        ('Abundance probes (not fits)',['baseline','abundance_minus','abundance_plus'])]
    fig,axes=plt.subplots(len(groups),1,figsize=(11,14),sharex=True)
    for ax,(title,names)in zip(axes,groups):
        ax.plot(wobs,fobs,color='black',lw=1,label='Kitt Peak corrected2005')
        minima=[float(fobs[(wobs>=7770.5)&(wobs<=7776.8)].min())]
        for name in names:
            if name in plotdata:
                w,f=plotdata[name];ax.plot(w,f,lw=1,label=name)
                minima.append(float(f[(w>=7770.5)&(w<=7776.8)].min()))
        ax.set_xlim(7770.5,7776.8);ax.set_ylim(max(0,min(minima)-.03),1.015);ax.set_title(title,fontsize=10);ax.legend(fontsize=7,ncol=3)
    axes[-1].set_xlabel('Air wavelength (Angstrom)');fig.tight_layout();fig.savefig(args.output/'parameter_responses.png',dpi=160);plt.close(fig)
    d=np.genfromtxt(args.output/'baseline_profiles.csv',delimiter=',',names=True)
    fig,ax=plt.subplots(figsize=(11,4));ax.plot(wobs,fobs,color='black',label='Kitt Peak corrected2005')
    for name in ('intrinsic','solar_broadening','solar_gravitational_shift'):
        ax.plot(d['wavelength_air_A'],d[name],label=name)
    ax.set_xlim(7770.5,7776.8);ax.set_ylim(.65,1.015);ax.legend(fontsize=8);ax.set_xlabel('Air wavelength (Angstrom)')
    ax.set_title('Canonical broadening and known gravitational-shift control; no fit')
    fig.tight_layout();fig.savefig(args.output/'broadening_alignment.png',dpi=160);plt.close(fig)
    by_name={c['name']:c for c in result['cases']}
    families=[('fine_T','Linear thermodynamic fields'),('logfine_T','Log pressure/density'),
              ('meshfine_T','Log pressure/density and positive layer thickness')]
    fig,axes=plt.subplots(3,1,figsize=(9,10),sharex=True)
    result['fine_temperature_responses']={}
    for ax,(prefix,title)in zip(axes,families):
        fine=sorted([c for c in result['cases']if c['name'].startswith(prefix)],key=lambda c:c['inputs']['teff'])
        if not fine:continue
        central=by_name[prefix+'5772'];table=[]
        for c in fine:
            table.append(dict(teff=c['inputs']['teff'],window_areas_mA=[r['model_area_mA']for r in c['variants']['solar_broadening']]))
        result['fine_temperature_responses'][prefix]=table
        for i,centre in enumerate(centres):
            baseline_area=central['variants']['solar_broadening'][i]['model_area_mA']
            ys=[c['variants']['solar_broadening'][i]['model_area_mA']-baseline_area for c in fine]
            ax.plot([c['inputs']['teff']for c in fine],ys,marker='o',label=f'{centre:.3f} A')
        ax.axhline(0,color='grey',lw=.5);ax.axvline(stars['teff'],color='grey',lw=.5)
        ax.set_ylabel('Area change (mA)');ax.set_title(title,fontsize=10);ax.legend(fontsize=8)
    axes[-1].set_xlabel('Interpolated Teff (K);10K sampling, interpolation error assessed separately')
    fig.tight_layout();fig.savefig(args.output/'fine_temperature_response.png',dpi=160);plt.close(fig)
    native=by_name['depth_cut1']['variants']['solar_broadening']
    definitions=[('temperature','linear','holdout_interp_5750'),('temperature','log_positive','holdout_log_5750'),
        ('temperature','log_positive_layers','holdout_mesh_5750'),('gravity','linear','holdout_logg_linear'),
        ('gravity','log_positive','holdout_logg_log'),('gravity','log_positive_layers','holdout_logg_mesh')]
    result['interpolation_holdouts']=[]
    from pipeline.bride_reference import read_marcs_structure
    from pipeline.bride_parameter_checks import interpolate_marcs_structure,marcs_filename
    store=require_codex_path('grids.marcs_standard');native_structure=read_marcs_structure(store/marcs_filename(5750,4.5))
    fieldmask=(native_structure['logtau_500']>=-3)&(native_structure['logtau_500']<=0)
    for axis,method,name in definitions:
        if name not in by_name:continue
        c=by_name[name];held=c['variants']['solar_broadening'];inputs=c['inputs']
        nodes=[]
        for t in inputs['bracket_teff']:
            for g in inputs['bracket_logg']:nodes.append(((t,g),read_marcs_structure(store/marcs_filename(t,g))))
        reconstructed,_=interpolate_marcs_structure(nodes,(5750,4.5),thermodynamic_scale=method)
        item=dict(axis=axis,method=method,name=name,
            area_difference_mA=[h['model_area_mA']-n['model_area_mA']for h,n in zip(held,native)],
            area_relative_error_percent=[100*(h['model_area_mA']/n['model_area_mA']-1)for h,n in zip(held,native)],
            structure_diagnostics_over_native_logtau500_minus3_to0=dict(
                max_temperature_error_K=float(np.max(abs(reconstructed['temperature_K'][fieldmask]-native_structure['temperature_K'][fieldmask]))),
                max_pressure_relative_error=float(np.max(abs(reconstructed['electron_pressure_Pa'][fieldmask]/native_structure['electron_pressure_Pa'][fieldmask]-1))),
                max_density_relative_error=float(np.max(abs(reconstructed['rho_kg_m3'][fieldmask]/native_structure['rho_kg_m3'][fieldmask]-1))),
                max_layer_thickness_relative_error=float(np.max(abs(np.diff(reconstructed['height_m'])[fieldmask[:-1]]/np.diff(native_structure['height_m'])[fieldmask[:-1]]-1)))))
        result['interpolation_holdouts'].append(item)
    fig,axes=plt.subplots(1,2,figsize=(10,4))
    for ax,axis in zip(axes,('temperature','gravity')):
        records=[r for r in result['interpolation_holdouts']if r['axis']==axis]
        for i,centre in enumerate(centres):
            ax.plot(range(len(records)),[abs(r['area_relative_error_percent'][i])for r in records],marker='o',label=f'{centre:.3f} A')
        ax.set_xticks(range(len(records)),[r['method'].replace('log_positive_layers','log P/rho/layers').replace('log_positive','log P/rho')for r in records],rotation=15)
        ax.set_yscale('log');ax.set_title(axis+' held-out native node');ax.set_ylabel('Absolute area error (%)');ax.legend(fontsize=7)
    fig.tight_layout();fig.savefig(args.output/'interpolation_holdouts.png',dpi=160);plt.close(fig)
    selected=[c for c in result['cases']if c['name'].startswith('meshfine_T')]
    if selected:
        fig,ax=plt.subplots(figsize=(11,4));ax.plot(wobs,fobs,color='black',label='Kitt Peak corrected2005')
        for c in selected:
            w,f=plotdata[c['name']];ax.plot(w,f,label=f"{c['inputs']['teff']}K")
        ax.set_xlim(7770.5,7776.8);ax.set_ylim(.65,1.015);ax.legend(ncol=3,fontsize=8)
        ax.set_title('Near-solar10K sequence, logg4.438, fixed O8.70; no fitting');ax.set_xlabel('Air wavelength (Angstrom)')
        fig.tight_layout();fig.savefig(args.output/'fine_solar_profiles.png',dpi=160);plt.close(fig)
    with (args.output/'analysis.json').open('x')as f:json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
    grid=dict(schema='codex.bride_experimental_grid_index/1',capabilities=result['capability_contract'],
        atom=run['atom'],data_store_key='experiments.bride_mac_grid',host=run['host'],nodes=[])
    for c in result['cases']:
        grid['nodes'].append(dict(name=c['name'],coordinates=c['inputs'],atmosphere=c['atmosphere'],
            population_relative_path=f"{c['source_directory_id']}/{c['name']}/populations.npz",
            populations=c['original_population_file'],diagnostics=c['final_diagnostics']))
    with (args.output/'grid_index.json').open('x')as f:json.dump(grid,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps(dict(successful_cases=len(result['cases']),failed_cases=result['failed_cases'],policy=result['orchestrator_policy']),indent=2))


if __name__=='__main__':main()
