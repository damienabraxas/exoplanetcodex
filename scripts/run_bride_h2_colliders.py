"""Paired fixed-O FeI pressure intervention using qualified full-column H2 tables."""
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
    parser.add_argument('--lines',type=Path,required=True)
    parser.add_argument('--mode',choices=('off','molecular'),required=True)
    parser.add_argument('--pressure-source',type=Path,required=True)
    parser.add_argument('--stark-source',type=Path,required=True)
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
    from pipeline.bride_h2_tables import H2Tables
    from pipeline.bride_h2_colliders import h2_collider_factory
    tables=H2Tables(require_codex_path('experiments.bride_molecular_tables'))
    factory=h2_collider_factory(provider,chosen,packet['support_vac_nm'],pressure['recipes'],stark['rows'],'cowley',tables,args.mode)
    captured=[]
    def line_provider(*args,**kwargs):
        backend=factory(*args,**kwargs);captured.append(backend);return backend
    input_manifest=dict(schema='codex.bride_h2_collider_input/1',source_manifest=fingerprint(origin/'manifest.json'),
        binding=fingerprint(origin/'experimental_input.json'),initial_populations=fingerprint(origin/'populations.npz'),
        lines=fingerprint(args.lines),pressure_source=fingerprint(args.pressure_source),stark_source=fingerprint(args.stark_source),
        generator=fingerprint(__file__),mode=args.mode,population_mode='fixed',stark_mode='cowley',
        H2='Colliders only: depleted HI plus H2; native EOS/continuum/O collisions unchanged',
        chemistry_source=tables.receipt,chemistry_adapter=fingerprint(ROOT/'pipeline/bride_h2_tables.py'),
        collider_provider=fingerprint(ROOT/'pipeline/bride_h2_colliders.py'),controls=m['controls'],
        code_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),scientific_admission='HOLD')
    (args.output/'input_manifest.json').write_text(json.dumps(input_manifest,indent=2,allow_nan=False)+'\n')
    ctx=make_reference_context(s,source,rays=m['controls']['rays'],qscale=m['controls']['qscale'],ng=False,
        fixed_populations=pop,fe_background=line_provider,extra_wavelengths_vac_nm=extra,
        experimental_column_cap=FiniteColumnCap(length))
    native_H_before=ctx.kwargs['eqPops'].atomicPops[PeriodicTable['H']].n.copy()
    native_Hminus_before=ctx.kwargs['eqPops'].HminPops.copy()
    radiation=converge_cap_radiation(ctx,rtol=1e-8,maxiter=1000,deadline=start+180)
    np.testing.assert_allclose(ctx.kwargs['eqPops'].atomicPops[PeriodicTable['O']].nStar,nstar,rtol=64*np.finfo(float).eps,atol=0)
    continuum=make_reference_context(s,source,rays=m['controls']['rays'],qscale=m['controls']['qscale'],
        fixed_populations=pop,continuum_only=True,fe_background=provider,extra_wavelengths_vac_nm=extra,
        experimental_column_cap=FiniteColumnCap(length))
    continuum_radiation=converge_cap_radiation(continuum,rtol=1e-8,maxiter=1000,deadline=start+180)
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
    if args.mode=='off':
        prior_stark=np.genfromtxt(ROOT/'data/results/rya1281/stark_coupled_v1/cowley/flux.csv',delimiter=',',names=True)
        np.testing.assert_array_equal(line/cont,prior_stark['normalized_flux'])
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
    np.testing.assert_array_equal(eq.atomicPops[PeriodicTable['H']].n,native_H_before)
    np.testing.assert_array_equal(eq.HminPops,native_Hminus_before)
    nh,nhe=neutral_colliders
    usedH,usedHe,usedH2=captured[0].colliders
    for backend in captured:
        for actual,expected in zip(backend.colliders,(usedH,usedHe,usedH2)):np.testing.assert_array_equal(actual,expected)
    np.testing.assert_array_equal(usedHe,nhe)
    np.savetxt(args.output/'neutral_colliders.csv',np.column_stack((s['logtau_500'],a.nHTot,nh,nhe,usedH,usedH2)),delimiter=',',
        header='logtau500,total_H_nuclei_m3,native_neutral_H_m3,neutral_He_m3,used_neutral_H_m3,used_H2_m3',comments='')
    import csv
    widthfile=args.output/'widths.csv'
    with widthfile.open('w') as f:
        writer=csv.writer(f,lineterminator="\n");writer.writerow(['line_id','depth_index','gamma_rad_s_inv','gamma_native_HHe_s_inv','gamma_H2_added_s_inv','gamma_HI_depletion_s_inv','gamma_Stark_s_inv','gamma_used_s_inv'])
        for line,r,sr in zip(chosen,pressure['recipes'],stark['rows']):
            widths=pressure_gamma(a.temperature,nh,nhe,np.zeros_like(nh),r)
            addition=pressure_gamma(a.temperature,np.zeros_like(nh),np.zeros_like(nh),usedH2,r)
            depletion=pressure_gamma(a.temperature,nh-usedH,np.zeros_like(nh),np.zeros_like(nh),r)
            st=stark_gamma(a.temperature,a.ne,sr,'cowley')
            used=pressure_gamma(a.temperature,usedH,usedHe,usedH2,r)+st+line['gamma_rad_s_inv']
            np.testing.assert_allclose(used,captured[0].line_width(line,a),rtol=2e-15,atol=0)
            for d,g in enumerate(widths):writer.writerow([line['line_id'],d,line['gamma_rad_s_inv'],g,addition[d],-depletion[d],st[d],used[d]])
    h2screen=[]
    for (_,node),receipt in zip(nodes,m['atmosphere']['sources']):
        sourcepath=require_codex_path('grids.marcs_standard')/receipt['name']
        text=sourcepath.read_text().splitlines();i=next(i for i,l in enumerate(text)if l.startswith(' k  lgPgas'))
        table=np.array([[float(x)for x in l.split()]for l in text[i+1:i+1+node['ndep']]])
        ratio_h2=10**(table[:,4]-table[:,2])
        h2screen.append(dict(source=receipt,HI_H2_header=text[i],max_H2_over_HI=float(ratio_h2.max()),max_H2_fraction_of_HI_width_term=float(.85*ratio_h2.max())))
    result=dict(schema='codex.bride_h2_collider_flux/1',scientific_admission='HOLD',mode=args.mode,
        generator=fingerprint(__file__),input_manifest=fingerprint(args.output/'input_manifest.json'),stark_source=fingerprint(args.stark_source),stark_kernel=fingerprint(ROOT/'pipeline/bride_fe_stark.py'),population_mode='fixed',stark_mode='cowley',chemistry_source=tables.receipt,chemistry_adapter=fingerprint(ROOT/'pipeline/bride_h2_tables.py'),collider_provider=fingerprint(ROOT/'pipeline/bride_h2_colliders.py'),pressure_source=fingerprint(args.pressure_source),pressure_provider=fingerprint(ROOT/'pipeline/bride_fe_pressure.py'),provider=fingerprint(ROOT/'pipeline/bride_trace_fe_lines.py'),
        line_coefficients=fingerprint(ROOT/'pipeline/bride_lte_blend.py'),lines=fingerprint(args.lines),
        partition=fingerprint(Path(get_data_path())/'pf_Kurucz.input'),selected_line_count=len(chosen),
        excluded_component_count=len(packet['excluded']),source_manifest=fingerprint(origin/'manifest.json'),
        source_populations=fingerprint(origin/'populations.npz'),source_flux=fingerprint(origin/'flux.csv'),
        binding=fingerprint(origin/'experimental_input.json'),wide_baseline=fingerprint(prior),
        no_blend_baseline_max_flux_difference=difference,flux=fingerprint(path),
        neutral_reservoir=fingerprint(args.output/'neutral_reservoir.csv'),neutral_colliders=fingerprint(args.output/'neutral_colliders.csv'),widths=fingerprint(widthfile),source_corner_H2_screen=h2screen,
        trace_to_native_FeI_range=[float(ratio.min()),float(ratio.max())],radiation=radiation,continuum_radiation=continuum_radiation,
        telemetry=dict(elapsed_seconds=time.perf_counter()-start,process_peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            host_load_start=load_start,host_load_end=os.getloadavg(),context_threads=1),
        caveats='Paired FeI pressure-collider intervention only. H2 formation depletes HI; native H populations, ne, continuum and O collisions unchanged. No coupled O SE, complete EOS or molecular opacity. Fixed source Cowley Stark. No fitted opacity/damping/abundance; scientific HOLD.')
    (args.output/'report.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(mode=args.mode,delta=difference,neutral_ratio=result['trace_to_native_FeI_range'],telemetry=result['telemetry'])))


if __name__=='__main__':main()
