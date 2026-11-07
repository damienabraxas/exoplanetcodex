"""Wide fixed-population O777 flux replay with explicit provenance and overlap control."""
import argparse
import json
from pathlib import Path
import platform
import sys
import time
import subprocess
from importlib.metadata import version
from dataclasses import replace

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case-directory',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if platform.system()!='Darwin' or version('lightweaver')!='0.16.1':raise RuntimeError('Mac/LW0.16.1 bounded replay only')
    from config.constants import require_codex_path
    from pipeline.bride_reference import fingerprint,read_marcs_structure,make_reference_context
    from pipeline.bride_parameter_checks import interpolate_marcs_structure,marcs_filename
    from pipeline.bride_solar_control import density_column_structure,bound_fe_factory
    from pipeline.bride_opacity_binding import load_bound_coefficients
    from pipeline.bride_top_cap import FiniteColumnCap,converge_cap_radiation
    from pipeline.multi_atom import read_multi_atom
    from lightweaver.atomic_table import PeriodicTable
    from pipeline.bride_observables import flux_quadrature, integrate_flux
    from pipeline.bride_top_cap import emergent_rays
    from pipeline.wavelength_util import air_to_vac
    import resource
    import os
    origin=args.case_directory;summary=json.loads((origin.parent/'summary.json').read_text())
    entry=next(x for x in summary['cases']if x['name']==origin.name)
    if fingerprint(origin/'manifest.json')!=entry['manifest']:raise ValueError('Changed run manifest')
    m=json.loads((origin/'manifest.json').read_text());case=m['case']
    if fingerprint(origin/'flux.csv')!=m['flux']:raise ValueError('Changed original flux')
    if m['status']!='CONVERGED_EXPERIMENT' or case.get('experimental_column_contract')!='density_integral_log_native_top' or m['experimental_collision_convention']!='barklem2018_MN_LTE_pair':
        raise ValueError('Unsupported or incomplete bound experiment')
    if fingerprint(origin/'experimental_input.json')!=m['experimental_bound_input']['receipt']:raise ValueError('Changed experiment binding')
    if fingerprint(ROOT/'pipeline/bride_reference.py')!=summary['context_builder'] or fingerprint(ROOT/'pipeline/bride_solar_control.py')!=m['experimental_bound_input']['adapter']:
        raise ValueError('Changed experimental context/provider code')
    for name in ('pipeline/bride_fe_background.py','pipeline/bride_top_cap.py'):
        recorded=subprocess.check_output(['git','show',summary['code_commit']+':'+name],cwd=ROOT)
        if recorded!=(ROOT/name).read_bytes():raise ValueError('Changed inherited provider/cap equation')
    atom_path=require_codex_path('experiments.bride_oxygen_atom')
    if fingerprint(atom_path)!=summary['atom']:raise ValueError('Changed atom')
    source=replace(read_multi_atom(atom_path),abundance=case['abundance'])
    store=require_codex_path('grids.marcs_standard');nodes=[]
    for (t,g),receipt in zip([(t,g)for t in case['bracket_teff']for g in case['bracket_logg']],m['atmosphere']['sources']):
        path=store/marcs_filename(t,g)
        if fingerprint(path)!=receipt:raise ValueError('Changed atmosphere source')
        nodes.append(((t,g),read_marcs_structure(path)))
    base,interpolation=interpolate_marcs_structure(nodes,(case['teff'],case['logg']),thermodynamic_scale=case['thermodynamic_scale'])
    s,contract=density_column_structure(base,nodes,interpolation)
    certificate=require_codex_path('experiments.bride_bound_point_coefficients')/case['bound_point_coefficients']
    provider,bound=bound_fe_factory(certificate,base,s,source.abundance,column_contract=contract)
    if bound!=json.loads((origin/'experimental_input.json').read_text()):raise ValueError('Reconstructed binding differs')
    data,_=load_bound_coefficients(certificate,base,source.abundance)
    prior_rate=json.loads((ROOT/'data/results/rya1281/solar_bound_v1/rates/5772/report.json').read_text())
    if fingerprint(origin/'manifest.json')!=prior_rate['source_manifest'] or fingerprint(origin/'populations.npz')!=prior_rate['source_populations']:
        raise ValueError('Saved nominal populations differ from archived fresh-rate audit')
    with np.load(origin/'populations.npz')as saved:pop=saved['populations'].copy();nstar=saved['nstar'].copy()
    args.output.mkdir(parents=True,exist_ok=False)
    start=time.perf_counter();load_start=os.getloadavg()
    length=s['column_mass_kg_m2'][0]/s['rho_kg_m3'][0]
    ctx=make_reference_context(s,source,rays=m['controls']['rays'],qscale=m['controls']['qscale'],
        fixed_populations=pop,fe_background=provider,extra_wavelengths_vac_nm=data['wavelength_vac_nm'],
        experimental_column_cap=FiniteColumnCap(length))
    np.testing.assert_allclose(ctx.kwargs['eqPops'].atomicPops[PeriodicTable['O']].nStar,nstar,rtol=64*np.finfo(float).eps,atol=0)
    radiation=converge_cap_radiation(ctx,rtol=1e-8,maxiter=1000,deadline=start+180)
    np.testing.assert_array_equal(ctx.kwargs['eqPops'].atomicPops[PeriodicTable['O']].n,pop)
    continuum=make_reference_context(s,source,rays=m['controls']['rays'],qscale=m['controls']['qscale'],
        fixed_populations=pop,continuum_only=True,fe_background=provider,
        extra_wavelengths_vac_nm=data['wavelength_vac_nm'],experimental_column_cap=FiniteColumnCap(length))
    continuum_radiation=converge_cap_radiation(continuum,rtol=1e-8,maxiter=1000,deadline=start+180)
    # Explicit wide grid includes >=2.6 A padding beyond every +/-5 A neighborhood.
    wave=7763.+np.arange(4001)*.005
    vac=air_to_vac(wave)/10
    mu,weight=flux_quadrature(m['controls']['flux_rays'])
    _,out=emergent_rays(ctx,vac,mu,weight)
    _,cout=emergent_rays(continuum,vac,mu,weight)
    line=integrate_flux(out,mu,weight);cont=integrate_flux(cout,mu,weight)
    normalized=line/cont
    old=np.genfromtxt(origin/'flux.csv',delimiter=',',names=True)
    # Recompute on the original coordinates, avoiding interpolation as a replay test.
    _,repeat=emergent_rays(ctx,air_to_vac(old['wavelength_air_A'])/10,mu,weight)
    _,crepeat=emergent_rays(continuum,air_to_vac(old['wavelength_air_A'])/10,mu,weight)
    replay=integrate_flux(repeat,mu,weight)/integrate_flux(crepeat,mu,weight)
    difference=float(np.max(abs(replay-old['normalized_flux'])))
    if difference>1e-6:raise ValueError('Original-grid replay exceeds original outer tolerance')
    np.testing.assert_array_equal(ctx.kwargs['eqPops'].atomicPops[PeriodicTable['O']].n,pop)
    np.testing.assert_array_equal(continuum.kwargs['eqPops'].atomicPops[PeriodicTable['O']].n,pop)
    np.savetxt(args.output/'flux.csv',np.column_stack((wave,line,cont,normalized)),delimiter=',',
        header='wavelength_air_A,Fnu_W_m2_Hz,continuum_Fnu_W_m2_Hz,normalized_flux',comments='')
    result=dict(schema='codex.bride_wide_flux_replay/1',scientific_admission='HOLD',
        generator=fingerprint(__file__),context_builder=fingerprint(ROOT/'pipeline/bride_reference.py'),
        source_manifest=fingerprint(origin/'manifest.json'),source_populations=fingerprint(origin/'populations.npz'),
        source_flux=fingerprint(origin/'flux.csv'),atom=summary['atom'],binding=fingerprint(origin/'experimental_input.json'),
        wave_air_A=[float(wave[0]),float(wave[-1])],step_A=.005,flux_rays=len(mu),
        original_grid_max_normalized_flux_difference=difference,overlap_threshold=1e-6,
        radiation=radiation,continuum_radiation=continuum_radiation,flux=fingerprint(args.output/'flux.csv'),
        telemetry=dict(elapsed_seconds=time.perf_counter()-start,process_peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            peak_rss_units='Darwin bytes; process lifetime, includes imports/setup',
            host_load_start=load_start,host_load_end=os.getloadavg(),context_threads=1),
        caveats='Fixed saved populations, no SE update. Same O bound-free continuum and Fe intervention; weak neighboring bound-bound/molecular opacity absent. Not full spectrum or continuum adoption.')
    (args.output/'report.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(result['telemetry']));print('Replay overlap',difference)


if __name__=='__main__':main()
