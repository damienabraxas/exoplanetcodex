"""Fresh coupled-O controls for natural and source-width FeII; fixed atmosphere/EOS."""
import argparse,csv,json,os,platform,resource,subprocess,sys,time
from dataclasses import replace
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def selected_line(plan):
    from config.constants import require_codex_path
    from pipeline.bride_reference import fingerprint
    from pipeline.wavelength_util import air_to_vac
    from pipeline.bride_lte_screen import source_eligibility
    intake_path=ROOT/'data/results/rya1281/anchor_influence_v1/candidate_intake.json'
    screen_path=ROOT/'data/results/rya1281/lte_opacity_screen_v1/report.json';screen=json.loads(screen_path.read_text())
    if fingerprint(intake_path)!=screen['intake']:raise ValueError('Changed strength intake')
    item=next(r for r in json.loads(intake_path.read_text())['rows'] if r['canonical_row']['line_id']==plan['line_id'])
    row=item['canonical_row']
    if row['species']!='Fe II' or source_eligibility(item)!='SOURCE_ELIGIBLE_FOR_NATURAL_ONLY_SCREEN':raise ValueError('Unqualified FeII row')
    for case in screen['cases']:
        if next(r for r in case['dispositions'] if r['line_id']==plan['line_id'])['reason']!='SCREENED_NATURAL_ONLY':raise ValueError('Strength screen refused line')
    with (ROOT/'data/linelists/canonical_gf.csv').open() as f:
        current=next(r for r in csv.DictReader(f) if r['line_id']==plan['line_id'])
    if current!=row:raise ValueError('Changed frozen canonical line')
    import io
    main=subprocess.check_output(['git','show','origin/main:data/linelists/canonical_gf.csv'],cwd=ROOT,text=True)
    if next(r for r in csv.DictReader(io.StringIO(main)) if r['line_id']==plan['line_id'])!=row:raise ValueError('Current main changed FeII line')
    ges=require_codex_path('engines.ispec')/'input/linelists/transitions/GESv6_atom_hfs_iso.420_920nm/atomic_lines.tsv'
    with ges.open() as f:matches=[r for r in csv.DictReader(f,delimiter='\t') if all(r[k]==v for k,v in item['source_rows'][0].items())]
    if len(matches)!=1 or int(matches[0]['spectrum_synthe_isotope'])!=0:raise ValueError('Missing/ambiguous/component source')
    src=matches[0]
    line=dict(line_id=plan['line_id'],lambda0_vac_nm=float(air_to_vac(float(row['wavelength_air_A'])))/10,
        excitation_eV=float(row['excitation_potential_eV']),loggf=float(row['log_gf']),gamma_rad_s_inv=float(src['turbospectrum_rad']))
    return line,dict(canonical=row,source=src,source_file=fingerprint(ges),intake=fingerprint(intake_path),screen=fingerprint(screen_path))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--mode',required=True,choices=('natural','cowley','provided'));p.add_argument('--case',required=True);p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    from importlib.metadata import version
    if platform.system()!='Darwin' or version('lightweaver')!='0.16.1':raise RuntimeError('Mac/LW0.16.1 control only')
    from config.constants import require_codex_path
    from pipeline.bride_reference import fingerprint,make_reference_context
    from pipeline.bride_thermal_sweep import reconstruct
    from pipeline.bride_solar_control import bound_fe_factory
    from pipeline.bride_opacity_binding import load_bound_coefficients
    from pipeline.bride_top_cap import FiniteColumnCap,converge_cap_radiation,emergent_rays
    from pipeline.bride_observables import flux_quadrature,integrate_flux
    from pipeline.multi_atom import read_multi_atom
    from pipeline.wavelength_util import air_to_vac
    from pipeline.bride_fe_stark import stark_fe_factory
    from pipeline.bride_feii_control import feii_factory
    from scripts.prepare_bride_fe_lines import read_packet
    from scripts.prepare_bride_fe_stark import source_rows
    from lightweaver.atomic_table import PeriodicTable
    plan_path=ROOT/'config/experiments/rya1281_feii_coupled.json';plan=json.loads(plan_path.read_text())
    if plan['schema']!='codex.bride_feii_coupled_plan/1' or args.case not in plan['cases']:raise ValueError('Outside declared plan')
    line,evidence=selected_line(plan)
    from pipeline.bride_feii_width_control import feii_width_factory
    width_path=ROOT/plan['width_source'];width=json.loads(width_path.read_text())
    if width['schema']!='codex.bride_feii_broadening/1' or width['line']!=line or width['evidence']!=evidence:raise ValueError('Changed verified FeII width identity')
    for relative,receipt in width['sources'].items():
        if fingerprint(require_codex_path('engines.ispec')/relative)!=receipt:raise ValueError('Changed FeII native width source')
    native_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=require_codex_path('engines.ispec'),text=True).strip()
    if native_commit!=width['source_commit']:raise ValueError('Changed native width source commit')
    recipe=dict(width['recipe'],line_id=line['line_id']);stark_row=dict(width['stark_row'],line_id=line['line_id'])
    natural_origin=require_codex_path('experiments.bride_feii_transfer')/args.case
    natural_run=json.loads((natural_origin/'report.json').read_text())
    natural_path=natural_origin/'feii/flux.csv'
    if fingerprint(natural_path)!=next(r for r in natural_run['legs'] if r['mode']=='feii')['flux']:raise ValueError('Changed natural baseline flux')
    natural=np.genfromtxt(natural_path,delimiter=',',names=True)
    frozen=json.loads((ROOT/'config/experiments/rya1281_continuum_anchor.json').read_text())
    case=next(c for c in frozen['cases'] if c['name']==args.case)
    origin=require_codex_path(case['resource'])/case['relative'];run=json.loads((origin/'report.json').read_text())
    previous=json.loads((ROOT/'data/results/rya1281/continuum_support_v1/report.json').read_text())
    if fingerprint(origin/'report.json')!=next(c for c in previous['cases'] if c['name']==args.case)['source_report']:raise ValueError('Changed frozen run')
    if fingerprint(origin/'flux.csv')!=run['flux'] or fingerprint(origin/'populations.npz')!=run['coupled']['populations']:raise ValueError('Changed source flux/populations')
    raw=np.genfromtxt(origin/'flux.csv',delimiter=',',names=True)
    with np.load(origin/'populations.npz') as data:pop=data['populations'].copy()
    base,s,contract,sources=reconstruct(run['case'])
    if sources!=run['atmosphere_sources']:raise ValueError('Changed atmosphere sources')
    atom_path=require_codex_path('experiments.bride_oxygen_atom');atom=replace(read_multi_atom(atom_path),abundance=run['case']['abundance'])
    folder=require_codex_path('experiments.bride_thermal_coefficients')/run['case']['bound_point_coefficients']
    bound_provider,bound=bound_fe_factory(folder,base,s,atom.abundance,column_contract=contract)
    if bound!=run['local_binding']:raise ValueError('Changed bound background')
    data,_=load_bound_coefficients(folder,base,atom.abundance)
    packet_path=ROOT/'data/results/rya1281/trace_fe_v1/lines.json';packet=read_packet(packet_path)
    pressure_path=ROOT/'data/results/rya1281/fe_pressure_v1/source/report.json'
    stark_path=ROOT/'data/results/rya1281/stark_coupled_v1/source/report.json'
    pressure=json.loads(pressure_path.read_text());stark=json.loads(stark_path.read_text())
    if fingerprint(packet_path)!=run['lines'] or fingerprint(pressure_path)!=run['pressure_source'] or fingerprint(stark_path)!=run['stark_source']:raise ValueError('Changed inherited line data')
    for key,name in [('pressure_provider','bride_fe_pressure.py'),('stark_kernel','bride_fe_stark.py'),('provider','bride_trace_fe_lines.py'),('line_coefficients','bride_lte_blend.py')]:
        if fingerprint(ROOT/'pipeline'/name)!=run[key]:raise ValueError('Changed inherited kernel')
    for source in (pressure,stark):
        for relative,receipt in source['sources'].items():
            if fingerprint(require_codex_path('engines.ispec')/relative)!=receipt:raise ValueError('Changed inherited native source')
    if stark['rows']!=source_rows():raise ValueError('Changed source Stark rows')
    provider=stark_fe_factory(bound_provider,packet['lines'],packet['support_vac_nm'],pressure['recipes'],stark['rows'],'cowley')
    wave=raw['wavelength_air_A'];vac=air_to_vac(wave)/10;extra=np.unique(np.r_[data['wavelength_vac_nm'],vac])
    controls=run['coupled']['controls'];length=s['column_mass_kg_m2'][0]/s['rho_kg_m3'][0]
    def construct(factory,*,fixed=None,continuum_only=False):
        return make_reference_context(s,atom,rays=controls['rays'],qscale=controls['qscale'],ng=controls['ng'] and fixed is None,fixed_populations=fixed,continuum_only=continuum_only,
            fe_background=factory,extra_wavelengths_vac_nm=extra,experimental_column_cap=FiniteColumnCap(length))
    if controls['seconds_per_case']!=plan['seconds_per_case']:raise ValueError('Changed coupled budget')
    if fingerprint(ROOT/'pipeline/bride_coupled_control.py')!=run['coupled']['helper']:raise ValueError('Changed inherited coupled solver')
    from pipeline.bride_coupled_control import solve_coupled
    from pipeline.bride_rate_audit import collision_classes,process_flow
    from pipeline.bride_numerical_validation import physical_rate_matrix
    result=dict(schema='codex.bride_feii_coupled/1',scientific_admission='HOLD',name=args.case,mode=args.mode,
        generator=fingerprint(__file__),plan=fingerprint(plan_path),line=line,evidence=evidence,atom=fingerprint(atom_path),
        width_source=fingerprint(width_path),source_run=fingerprint(origin/'report.json'),source_manifest=run['source_manifest'],
        source_flux=run['source_flux'],source_populations=fingerprint(origin/'populations.npz'),atmosphere_sources=sources,
        local_binding=bound,controls=controls,feii_provider=fingerprint(ROOT/'pipeline/bride_feii_control.py'),
        width_provider=fingerprint(ROOT/'pipeline/bride_feii_width_control.py'),width_kernel=fingerprint(ROOT/'pipeline/bride_feii_profile.py'),
        inherited_pressure=fingerprint(pressure_path),inherited_stark=fingerprint(stark_path),
        scope=plan['scope'],initialization='Exact same frozen converged O populations in every mode; not sequential warm starts')
    args.output.mkdir(parents=True,exist_ok=False);(args.output/'input.json').write_text(json.dumps(result,indent=2)+'\n')
    start=time.perf_counter()
    chosen=feii_factory(provider,line,packet['support_vac_nm']) if args.mode=='natural' else feii_width_factory(provider,line,packet['support_vac_nm'],recipe,stark_row,args.mode)
    ctx=construct(chosen);active=ctx.activeAtoms[0];target_nstar=active.nStar.copy()
    if pop.shape!=active.n.shape or not np.allclose(pop.sum(axis=0),active.nTotal,rtol=1e-10,atol=0):raise ValueError('Changed initial population identity')
    active.n[:]=pop
    original=dict(controls=controls,experimental_collision_convention='barklem2018_MN_LTE_pair',
        native_spline_control=run['coupled']['source_spline'],author_charge_control=run['coupled']['source_mn'])
    final,coupled=solve_coupled(ctx,atom,original,args.output,start)
    np.testing.assert_array_equal(active.nStar,target_nstar)
    groups=collision_classes(active.atomicModel,ctx.atmos.pyAtmos,ctx.kwargs['eqPops'])
    collision=sum(groups.values());active.compute_collisions(fillDiagonal=False)
    np.testing.assert_allclose(collision,active.C,rtol=256*np.finfo(float).eps,atol=0)
    radiative=np.zeros_like(collision)
    for transition in active.trans:
        radiative[transition.j,transition.i]+=transition.Rij
        radiative[transition.i,transition.j]+=transition.Rji
    net,gross=process_flow(collision+radiative,final)
    residual=float(np.max(abs(net)/np.maximum(gross,np.finfo(float).tiny)))
    conservation=float(np.max(abs(final.sum(axis=0)/active.nTotal-1)))
    if residual>=controls['rtol'] or conservation>1e-10:raise ValueError('Independent fresh physical balance fails original accuracy')
    independent_path=args.output/'independent_rates.npz'
    np.savez_compressed(independent_path,collision_rates_s_inv=collision,radiative_rates_s_inv=radiative,net_flow_m3_s=net,gross_flow_m3_s=gross)
    result['independent']=dict(physical_rate_residual=residual,particle_conservation_error=conservation,
        classes=list(groups),rates=fingerprint(independent_path),kernel=fingerprint(ROOT/'pipeline/bride_rate_audit.py'),
        scope='Fresh individual collision classes plus transition radiative rates; no accelerated Gamma used')
    continuum=construct(bound_provider,fixed=final,continuum_only=True)
    continuum_radiation=converge_cap_radiation(continuum,rtol=1e-8,maxiter=1000,deadline=start+controls['seconds_per_case'])
    mu,weights=flux_quadrature(controls['flux_rays'])
    _,out=emergent_rays(ctx,vac,mu,weights);flux=integrate_flux(out,mu,weights)
    _,cout=emergent_rays(continuum,vac,mu,weights);cont=integrate_flux(cout,mu,weights)
    for context in (ctx,continuum):np.testing.assert_array_equal(context.kwargs['eqPops'].atomicPops[PeriodicTable['O']].n,final)
    flux_path=args.output/'flux.csv'
    np.savetxt(flux_path,np.column_stack((wave,flux,cont,flux/cont)),delimiter=',',header='wavelength_air_A,Fnu_W_m2_Hz,continuum_Fnu_W_m2_Hz,normalized_flux',comments='')
    frozen_leg=require_codex_path('experiments.bride_feii_transfer')/args.case/'feii' if args.mode=='natural' else require_codex_path('experiments.bride_feii_width_transfer')/args.mode/args.case/'feii'
    frozen_receipt=json.loads((frozen_leg/'report.json').read_text())
    if fingerprint(frozen_leg/'flux.csv')!=frozen_receipt['flux']:raise ValueError('Changed frozen comparison')
    frozen=np.genfromtxt(frozen_leg/'flux.csv',delimiter=',',names=True)
    result.update(coupled=coupled,continuum_radiation=continuum_radiation,flux=fingerprint(flux_path),
        frozen_comparison=fingerprint(frozen_leg/'report.json'),max_normalized_flux_change_from_fixed=float(np.max(abs(flux/cont-frozen['normalized_flux']))),
        max_continuum_relative_change=float(np.max(abs(cont/frozen['continuum_Fnu_W_m2_Hz']-1))),
        telemetry=dict(seconds=time.perf_counter()-start,threads=1,host=platform.node(),load=os.getloadavg(),peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
    (args.output/'report.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(mode=args.mode,iterations=coupled['iterations'],independent=result['independent'],flux_change=result['max_normalized_flux_change_from_fixed'],telemetry=result['telemetry'])))


if __name__=='__main__':main()
