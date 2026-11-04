"""Fresh source-matched solar Teff/gravity O solve with optional sampling refinement."""
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
    parser.add_argument('--refine-numerics',action='store_true')
    parser.add_argument('--warm-start-run',type=Path)
    parser.add_argument('--ng-control',choices=('record','off'))
    parser.add_argument('--case-directory',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--lines',type=Path,required=True)
    parser.add_argument('--plan',type=Path,required=True)
    parser.add_argument('--case',required=True)
    parser.add_argument('--coefficients',type=Path,required=True)
    parser.add_argument('--pressure-source',type=Path,required=True)
    parser.add_argument('--stark-source',type=Path,required=True)
    args=parser.parse_args()
    if args.warm_start_run is not None and not args.refine_numerics:raise ValueError('Matched initial guess is restricted to numerical refinement')
    if args.ng_control is not None and (not args.refine_numerics or args.warm_start_run is None):raise ValueError('Ng comparison requires refined matched initialization')
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
    m=json.loads((origin/'manifest.json').read_text())
    from pipeline.bride_thermal_sweep import validate_plan,reconstruct,warm_start,effective_controls,validate_spline_receipt
    plan=json.loads(args.plan.read_text());validate_plan(plan,json.loads((ROOT/'config/experiments/rya1281_mac_solar_bound_control.json').read_text()))
    selected=[c for c in plan['cases']if c['name']==args.case]
    if len(selected)!=1:raise ValueError('Missing/ambiguous requested sweep case')
    case=selected[0]
    if plan['controls']!=m['controls']:raise ValueError('Changed reference solve controls')
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
    base,s,contract,atmosphere_sources=reconstruct(case)
    certificate=args.coefficients/case['bound_point_coefficients']
    provider,bound=bound_fe_factory(certificate,base,s,source.abundance,column_contract=contract)
    m['controls']=effective_controls(m['controls'],args.refine_numerics)
    ng_experiment=None
    if args.ng_control is not None:
        ng_path=ROOT/'config/experiments/rya1281_ng_control.json';ng_plan=json.loads(ng_path.read_text())
        if (ng_plan['schema']!='codex.bride_ng_experiment/1' or case['name']!=ng_plan['case']
                or args.ng_control not in ng_plan['modes'] or args.plan.resolve()!=ROOT/ng_plan['physical_plan']):raise ValueError('Outside declared Ng experiment')
        m['controls']['ng']=args.ng_control=='record';ng_experiment=fingerprint(ng_path)
    data,cert=load_bound_coefficients(certificate,base,source.abundance)
    producer=json.loads((certificate/cert['source_report']['name']).read_text())
    producer_case=next(c for c in producer['cases']if c['name']==case['name'])
    if producer_case['requested_xi_km_s']!=case['xi']:raise ValueError('Relabelled producer microturbulence')
    np.testing.assert_array_equal(data['requested_vturb_m_s'],np.full(s['ndep'],case['xi']*1000))
    np.testing.assert_array_equal(data['producer_vturb_m_s'],producer_case['producer_xi_m_s'])
    prior_rate=json.loads((ROOT/'data/results/rya1281/solar_bound_v1/rates/5772/report.json').read_text())
    if fingerprint(origin/'manifest.json')!=prior_rate['source_manifest'] or fingerprint(origin/'populations.npz')!=prior_rate['source_populations']:
        raise ValueError('Saved nominal populations differ from archived fresh-rate audit')
    with np.load(origin/'populations.npz')as saved:pop=saved['populations'].copy();nstar=saved['nstar'].copy()
    from pipeline.bride_trace_fe_lines import trace_fe_factory
    from lightweaver.atomic_table import KuruczPfTable
    from lightweaver.utils import get_data_path
    from scripts.prepare_bride_fe_lines import read_packet
    packet=read_packet(args.lines)
    if packet['schema']!='codex.bride_trace_fe_lines/1':raise ValueError('Wrong trace source schema')
    ges=require_codex_path('engines.ispec')/'input/linelists/transitions/GESv6_atom_hfs_iso.420_920nm/atomic_lines.tsv'
    if fingerprint(ges)!=packet['source'] or fingerprint(ROOT/'data/linelists/canonical_gf.csv')!=packet['canonical']:
        raise ValueError('Changed trace source data')
    # Dense VIS wavelengths resolve the added opacity in the radiation/cap solve.
    wave=7763.+np.arange(4001)*.005;vac=air_to_vac(wave)/10
    extra=np.unique(np.r_[data['wavelength_vac_nm'],vac])
    args.output.mkdir(parents=True,exist_ok=False);start=time.perf_counter();load_start=os.getloadavg()
    length=s['column_mass_kg_m2'][0]/s['rho_kg_m3'][0]
    chosen=packet['lines']
    pressure=json.loads(args.pressure_source.read_text());pfolder=args.pressure_source.parent
    if pressure['schema']!='codex.bride_fe_pressure_source/1' or pressure['lines']!=fingerprint(args.lines):raise ValueError('Wrong pressure source/line packet')
    if pressure['pressure_kernel']!=fingerprint(ROOT/'pipeline/bride_fe_pressure.py'):raise ValueError('Changed pressure kernel')
    for relative,receipt in pressure['sources'].items():
        if fingerprint(require_codex_path('engines.ispec')/relative)!=receipt:raise ValueError('Changed original pressure source')
    for field,name in [('binary','pressure_oracle'),('driver','driver.f'),('input','input.txt'),('native_stdout','stdout.txt')]:
        if fingerprint(pfolder/name)!=pressure[field]:raise ValueError('Changed native pressure evidence')
    native={}
    for row in (pfolder/'stdout.txt').read_text().splitlines():
        parts=row.split()
        if parts and parts[0]=='BRIDE_RECIPE':native[int(parts[1])]=(parts[2],float(parts[3]),float(parts[4]))
    if len(pressure['recipes'])!=len(chosen):raise ValueError('Incomplete native recipes')
    for i,(r,line,e) in enumerate(zip(pressure['recipes'],chosen,packet['evidence']),1):
        if (r['line_id']!=line['line_id'] or native[i]!=(r['recipe'],r['sigma_au'],r['alpha'])
                or r['log_gamma6_10000']!=float(e['source_row']['turbospectrum_fdamp'])):raise ValueError('Pressure recipe differs from native/source row')
    from pipeline.bride_fe_pressure import pressure_fe_factory,pressure_gamma
    from pipeline.bride_fe_stark import stark_fe_factory,stark_gamma
    from scripts.prepare_bride_fe_stark import source_rows
    stark=json.loads(args.stark_source.read_text())
    if stark['schema']!='codex.bride_fe_stark_source/1'or stark['rows']!=source_rows():raise ValueError('Changed Stark source rows')
    if stark['kernel']!=fingerprint(ROOT/'pipeline/bride_fe_stark.py'):raise ValueError('Changed Stark kernel')
    for relative,receipt in stark['sources'].items():
        if fingerprint(require_codex_path('engines.ispec')/relative)!=receipt:raise ValueError('Changed Stark original source')
    line_provider=stark_fe_factory(provider,chosen,packet['support_vac_nm'],pressure['recipes'],stark['rows'],'cowley')
    input_manifest=dict(schema='codex.bride_thermal_input/1',numeric_refinement=args.refine_numerics,atmosphere_sources=atmosphere_sources,case=case,plan=fingerprint(args.plan),sweep_adapter=fingerprint(ROOT/'pipeline/bride_thermal_sweep.py'),certificate=fingerprint(certificate/'manifest.json'),source_manifest=fingerprint(origin/'manifest.json'),
        source_binding=fingerprint(origin/'experimental_input.json'),local_binding=bound,initial_populations=fingerprint(origin/'populations.npz'),
        lines=fingerprint(args.lines),pressure_source=fingerprint(args.pressure_source),stark_source=fingerprint(args.stark_source),
        generator=fingerprint(__file__),mode='cowley',population_mode='coupled',
        H2='OMITTED: atomic baseline retained; H2 EOS/continuum/O-collision consistency not admitted',controls=m['controls'],
        code_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),scientific_admission='HOLD')
    if args.ng_control is not None:
        input_manifest.update(ng_control=args.ng_control,ng_experiment=ng_experiment,
            iteration_recorder=fingerprint(ROOT/'pipeline/bride_iteration_trace.py'))
    ctx=make_reference_context(s,source,rays=m['controls']['rays'],qscale=m['controls']['qscale'],ng=m['controls']['ng'],
        fixed_populations=None,fe_background=line_provider,extra_wavelengths_vac_nm=extra,
        experimental_column_cap=FiniteColumnCap(length))
    from pipeline.bride_coupled_control import solve_coupled
    target_nstar=ctx.activeAtoms[0].nStar.copy()
    initialization=dict(mode='original_solar_departures',populations=fingerprint(origin/'populations.npz'))
    if args.warm_start_run is not None:
        from pipeline.bride_thermal_initial import matched_start
        expected=dict(input_manifest,baseline_controls=plan['controls'])
        initial,initialization=matched_start(args.warm_start_run,expected,target_nstar,ctx.activeAtoms[0].nTotal)
        ctx.activeAtoms[0].n[:]=initial
    else:ctx.activeAtoms[0].n[:]=warm_start(pop,nstar,target_nstar,ctx.activeAtoms[0].nTotal)
    input_manifest['initialization']=initialization
    (args.output/'input_manifest.json').write_text(json.dumps(input_manifest,indent=2,allow_nan=False)+'\n')
    from pipeline.bride_native_spline import NativeTautsp,install_native_collision_control
    from pipeline.bride_charge_contract import install_author_mn_pair
    oracle=NativeTautsp(require_codex_path('experiments.bride_multi_spline_oracle'))
    local_spline=install_native_collision_control(ctx,source,oracle,source_directions=True,reversible_ch0=False)
    validate_spline_receipt(local_spline,m['native_spline_control'],source,s['ndep'])
    local_mn=install_author_mn_pair(ctx,source,require_codex_path('experiments.bride_barklem_oxygen'),oracle)
    if local_mn!=m['author_charge_control']:raise ValueError('Changed author collision source/prescription')
    m['native_spline_control']=local_spline
    if args.ng_control is None:pop,coupled=solve_coupled(ctx,source,m,args.output,start)
    else:
        from pipeline.bride_iteration_trace import StatEquilTrace
        with StatEquilTrace(ctx,args.output/'iteration_updates.jsonl',args.output/'history.jsonl')as observed:
            pop,coupled=solve_coupled(observed,source,m,args.output,start)
        coupled.update(iteration_updates=fingerprint(args.output/'iteration_updates.jsonl'),
            recorder=fingerprint(ROOT/'pipeline/bride_iteration_trace.py'),recorder_overhead_seconds=observed.overhead_seconds)
    radiation=coupled['fresh_radiation']
    np.testing.assert_array_equal(ctx.kwargs['eqPops'].atomicPops[PeriodicTable['O']].nStar,target_nstar)
    continuum=make_reference_context(s,source,rays=m['controls']['rays'],qscale=m['controls']['qscale'],
        fixed_populations=pop,continuum_only=True,fe_background=provider,extra_wavelengths_vac_nm=extra,
        experimental_column_cap=FiniteColumnCap(length))
    continuum_radiation=converge_cap_radiation(continuum,rtol=1e-8,maxiter=1000,deadline=start+m['controls']['seconds_per_case'])
    mu,weight=flux_quadrature(m['controls']['flux_rays'])
    _,out=emergent_rays(ctx,vac,mu,weight);_,cout=emergent_rays(continuum,vac,mu,weight)
    line=integrate_flux(out,mu,weight);cont=integrate_flux(cout,mu,weight)
    for context in (ctx,continuum):np.testing.assert_array_equal(context.kwargs['eqPops'].atomicPops[PeriodicTable['O']].n,pop)
    old=json.loads((ROOT/'data/results/rya1281/wide_continuum_v1/replay.json').read_text())
    prior=require_codex_path('experiments.bride_mac_grid')/'rya1281-wide-flux-v4'/'flux.csv'
    if fingerprint(prior)!=old['flux']:raise ValueError('Changed no-blend wide baseline')
    comparison=np.genfromtxt(prior,delimiter=',',names=True)
    np.testing.assert_array_equal(wave,comparison['wavelength_air_A'])
    difference=float(np.max(abs(line/cont-comparison['normalized_flux'])))
    path=args.output/'flux.csv'
    np.savetxt(path,np.column_stack((wave,line,cont,line/cont)),delimiter=',',header='wavelength_air_A,Fnu_W_m2_Hz,continuum_Fnu_W_m2_Hz,normalized_flux',comments='')
    eq=ctx.kwargs['eqPops'];a=ctx.atmos.pyAtmos
    partition=KuruczPfTable(atomicAbundance=eq.abundance)[PeriodicTable['Fe']]
    ion=partition.lte_ionisation(a)[0]
    native=eq.atomicPops[PeriodicTable['Fe']]
    native_neutral=native.n[[l.stage==0 for l in native.model.levels]].sum(axis=0)
    ratio=ion/native_neutral
    np.savetxt(args.output/'neutral_reservoir.csv',np.column_stack((s['logtau_500'],native_neutral,ion,ratio)),
        delimiter=',',header='logtau500,native_FeI_m3,Kurucz_trace_FeI_m3,trace_over_native',comments='')
    neutral_colliders=[]
    for element in ('H','He'):
        state=eq.atomicPops[PeriodicTable[element]]
        neutral_colliders.append(state.n[[l.stage==0 for l in state.model.levels]].sum(axis=0))
    nh,nhe=neutral_colliders
    np.savetxt(args.output/'neutral_colliders.csv',np.column_stack((s['logtau_500'],a.nHTot,nh,nhe)),delimiter=',',
        header='logtau500,total_H_nuclei_m3,neutral_H_m3,neutral_He_m3',comments='')
    import csv
    widthfile=args.output/'widths.csv'
    with widthfile.open('w') as f:
        writer=csv.writer(f,lineterminator="\n");writer.writerow(['line_id','depth_index','gamma_rad_s_inv','gamma_HHe_s_inv','gamma_Stark_s_inv','gamma_used_s_inv'])
        for line,r,sr in zip(chosen,pressure['recipes'],stark['rows']):
            widths=pressure_gamma(a.temperature,nh,nhe,np.zeros_like(nh),r)
            st=stark_gamma(a.temperature,a.ne,sr,'cowley')
            for d,g in enumerate(widths):writer.writerow([line['line_id'],d,line['gamma_rad_s_inv'],g,st[d],line['gamma_rad_s_inv']+g+st[d]])
    result=dict(schema='codex.bride_thermal_flux/1',numeric_refinement=args.refine_numerics,atmosphere_sources=atmosphere_sources,case=case,plan=fingerprint(args.plan),sweep_adapter=fingerprint(ROOT/'pipeline/bride_thermal_sweep.py'),certificate=fingerprint(certificate/'manifest.json'),scientific_admission='HOLD',mode='cowley',
        generator=fingerprint(__file__),input_manifest=fingerprint(args.output/'input_manifest.json'),stark_source=fingerprint(args.stark_source),stark_kernel=fingerprint(ROOT/'pipeline/bride_fe_stark.py'),population_mode='coupled',coupled=coupled,pressure_source=fingerprint(args.pressure_source),pressure_provider=fingerprint(ROOT/'pipeline/bride_fe_pressure.py'),provider=fingerprint(ROOT/'pipeline/bride_trace_fe_lines.py'),
        line_coefficients=fingerprint(ROOT/'pipeline/bride_lte_blend.py'),lines=fingerprint(args.lines),
        partition=fingerprint(Path(get_data_path())/'pf_Kurucz.input'),selected_line_count=len(chosen),
        excluded_component_count=len(packet['excluded']),source_manifest=fingerprint(origin/'manifest.json'),
        source_populations=fingerprint(origin/'populations.npz'),source_flux=fingerprint(origin/'flux.csv'),
        source_binding=fingerprint(origin/'experimental_input.json'),local_binding=bound,wide_baseline=fingerprint(prior),
        no_blend_baseline_max_flux_difference=difference,flux=fingerprint(path),
        neutral_reservoir=fingerprint(args.output/'neutral_reservoir.csv'),neutral_colliders=fingerprint(args.output/'neutral_colliders.csv'),widths=fingerprint(widthfile),
        trace_to_native_FeI_range=[float(ratio.min()),float(ratio.max())],radiation=radiation,continuum_radiation=continuum_radiation,
        telemetry=dict(elapsed_seconds=time.perf_counter()-start,process_peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            host_load_start=load_start,host_load_end=os.getloadavg(),context_threads=1),
        caveats='41 source-matched FeI lines, Explicit neutral-H/He widths plus NONE/supplied/Cowley Stark prescription. H2 omitted to retain atomic EOS/continuum/O collision baseline; no HFS/isotope/molecular/other-species opacity; Teff/logg use interpolated source fields with new local densities and LTE populations, not a fresh atmosphere-equilibrium solution. Kurucz LTE reservoir is separate from limited native Fe atom, ratio reported. Same transfer equation, no spectrum multiplication or opacity scale; Population mode explicitly recorded, optional bounded O SE; no Fe abundance measurement or opacity/continuum adoption.')
    result['initialization']=initialization
    if args.ng_control is not None:result.update(ng_control=args.ng_control,ng_experiment=ng_experiment)
    (args.output/'report.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(mode='cowley',delta=difference,neutral_ratio=result['trace_to_native_FeI_range'],telemetry=result['telemetry'])))


if __name__=='__main__':main()
