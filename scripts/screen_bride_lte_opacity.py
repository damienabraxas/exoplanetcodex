"""Local LTE opacity screen at frozen anchors; no formal RT/SE or adoption."""
import argparse,csv,io,json,os,platform,subprocess,sys,time,resource
from dataclasses import replace
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    from importlib.metadata import version
    if platform.system()!='Darwin' or version('lightweaver')!='0.16.1':raise RuntimeError('Declared Mac/LW0.16.1 screen only')
    import lightweaver.constants as c
    from lightweaver.atomic_table import KuruczPfTable,PeriodicTable
    from lightweaver.utils import get_data_path
    from config.constants import require_codex_path
    from pipeline.bride_reference import fingerprint,make_reference_context
    from pipeline.bride_thermal_sweep import reconstruct
    from pipeline.bride_solar_control import bound_fe_factory
    from pipeline.bride_top_cap import FiniteColumnCap
    from pipeline.multi_atom import read_multi_atom
    from pipeline.wavelength_util import air_to_vac
    from pipeline.reflected_solar_rv import SOLAR_GRAVITATIONAL_REDSHIFT
    from pipeline.bride_lte_blend import lte_line_coefficients
    from pipeline.bride_lte_screen import source_eligibility,independent_ion_populations
    plan_path=ROOT/'config/experiments/rya1281_lte_screen.json';plan=json.loads(plan_path.read_text())
    if plan['schema']!='codex.bride_lte_screen_plan/1':raise ValueError('Unknown plan')
    influence_path=ROOT/plan['influence_report'];influence=json.loads(influence_path.read_text())
    intake_path=ROOT/plan['intake'];intake=json.loads(intake_path.read_text())['rows']
    if fingerprint(ROOT/'data/linelists/canonical_gf.csv')!=influence['canonical']:raise ValueError('Changed frozen catalog')
    main_text=subprocess.check_output(['git','show',plan['main_commit']+':data/linelists/canonical_gf.csv'],text=True,cwd=ROOT)
    current={r['line_id']:r for r in csv.DictReader(io.StringIO(main_text))}
    changed=[r['canonical_row']['line_id'] for r in intake if current.get(r['canonical_row']['line_id'])!=r['canonical_row']]
    ges=require_codex_path('engines.ispec')/'input/linelists/transitions/GESv6_atom_hfs_iso.420_920nm/atomic_lines.tsv'
    if fingerprint(ges)!=influence['original']:raise ValueError('Changed original line source')
    with ges.open() as stream:original=[r for r in csv.DictReader(stream,delimiter='\t') if 7763<=float(r['wave_A'])<=7783]
    pixels_path=ROOT/'data/results/rya1281/continuum_anchor_v1/pixel_diagnostics.csv'
    if fingerprint(pixels_path)!=influence['derived_pixels']:raise ValueError('Changed anchor grid')
    with pixels_path.open() as stream:pixels=list(csv.DictReader(stream))
    frozen=json.loads((ROOT/'config/experiments/rya1281_continuum_anchor.json').read_text())
    prior=json.loads((ROOT/'data/results/rya1281/continuum_support_v1/report.json').read_text())
    result=dict(schema='codex.bride_lte_opacity_screen/1',scientific_admission='HOLD',
        generator=fingerprint(__file__),adapter=fingerprint(ROOT/'pipeline/bride_lte_screen.py'),plan=fingerprint(plan_path),
        kernel=fingerprint(ROOT/'pipeline/bride_lte_blend.py'),context=fingerprint(ROOT/'pipeline/bride_reference.py'),
        intake=fingerprint(intake_path),influence=fingerprint(influence_path),partition=fingerprint(Path(get_data_path())/'pf_Kurucz.input'),
        main_commit=plan['main_commit'],main_catalog_changed_ids=changed,lightweaver=version('lightweaver'),
        frame='Observer AIR -> inverse canonical relativistic gravitational shift in AIR -> vacuum nm once. Local thermal+xi profiles; macro/rotation/instrument convolution not applied to local coefficients.',
        limitations='Coefficient screen only, not a flux/EW prediction or abundance inference. Fixed-ne atomic LTE chemistry, no pressure/Stark/molecular opacity, charge/EOS/SE feedback or new opacity installed. Historical3percent anchors retained; main15percent model-guided bound unexercised.',cases=[])
    start=time.perf_counter();depth_rows=[];line_rows=[]
    for case,oldcase in zip(frozen['cases'],prior['cases']):
        if time.perf_counter()-start>plan['seconds']:raise TimeoutError('Declared screen budget exhausted')
        directory=require_codex_path(case['resource'])/case['relative'];run=json.loads((directory/'report.json').read_text())
        if fingerprint(directory/'report.json')!=oldcase['source_report']:raise ValueError('Changed frozen run')
        popfile=directory/'populations.npz'
        if fingerprint(popfile)!=run['coupled']['populations']:raise ValueError('Changed fixed oxygen populations')
        with np.load(popfile) as data:pop=data['populations'].copy()
        base,s,contract,receipts=reconstruct(run['case'])
        if receipts!=run['atmosphere_sources']:raise ValueError('Changed atmosphere sources')
        atompath=require_codex_path('experiments.bride_oxygen_atom');atom=replace(read_multi_atom(atompath),abundance=run['case']['abundance'])
        certificate=require_codex_path('experiments.bride_thermal_coefficients')/run['case']['bound_point_coefficients']
        provider,bound=bound_fe_factory(certificate,base,s,atom.abundance,column_contract=contract)
        if bound!=run['local_binding'] or result['partition']!=run['partition']:raise ValueError('Changed frozen coefficient/population contract')
        rr=[r for r in pixels if r['case']==case['name'] and r['coefficient_support']=='True']
        observer=np.array(sorted({float(r['centre_air_A'])+float(r['offset_air_A']) for r in rr}))
        beta=SOLAR_GRAVITATIONAL_REDSHIFT*1000/c.CLight
        wave=air_to_vac(observer*np.sqrt((1-beta)/(1+beta)))/10
        ctx=make_reference_context(s,atom,rays=1,qscale=1,fixed_populations=pop,continuum_only=True,fe_background=provider,
            extra_wavelengths_vac_nm=wave,experimental_column_cap=FiniteColumnCap(s['column_mass_kg_m2'][0]/s['rho_kg_m3'][0]))
        a=ctx.atmos.pyAtmos;eq=ctx.kwargs['eqPops'];grid=np.asarray(ctx.spect.wavelength);indices=np.searchsorted(grid,wave)
        np.testing.assert_array_equal(grid[indices],wave)
        background=np.asarray(ctx.background.chi)[indices]-np.asarray(ctx.background.sca)[indices]
        if np.any(background<=0) or not np.isfinite(background).all():raise ValueError('Invalid baseline absorption')
        table=KuruczPfTable(atomicAbundance=eq.abundance);cache={};dispositions=[];groups={};footprints={}
        reference_depth=int(np.argmin(abs(s['logtau_500']-plan['reference_logtau500'])))
        for centre in (7771.944,7774.166,7775.388):
            for mode in ('synthesis','iag'):
                rows=[r for r in rr if float(r['centre_air_A'])==centre and r['reference']==mode]
                ix=np.array([int(np.searchsorted(observer,float(r['centre_air_A'])+float(r['offset_air_A']))) for r in rows])
                weights=np.abs([float(r['centre_weight']) for r in rows]);weights=weights/weights.sum()
                footprints[f'{centre}_{mode}']=(ix,weights)
        for item in intake:
            if time.perf_counter()-start>plan['seconds']:raise TimeoutError('Declared screen budget exhausted')
            row=item['canonical_row'];line_id=row['line_id'];reason=source_eligibility(item,changed_on_main=line_id in changed)
            if reason!='SOURCE_ELIGIBLE_FOR_NATURAL_ONLY_SCREEN':dispositions.append(dict(line_id=line_id,reason=reason));continue
            src=item['source_rows'][0];matches=[r for r in original if all(r[k]==src[k] for k in src)]
            if len(matches)!=1:raise ValueError('Intake source fragment no longer unique')
            if int(matches[0]['spectrum_synthe_isotope'])!=0:
                dispositions.append(dict(line_id=line_id,reason='ISOTOPE_COMPONENT_UNQUALIFIED'));continue
            element=PeriodicTable[int(row['key_z'])];stage=int(float(row['ion']))-1
            if element not in cache:
                part=table[element]
                if a.temperature.min()<part.Tpf[0] or a.temperature.max()>part.Tpf[-1]:
                    cache[element]=None
                else:
                    ions=part.lte_ionisation(a);independent=independent_ion_populations(part,a);total=part.abundance*a.nHTot
                    error=float(np.max(abs(ions-independent)/total))
                    closure=float(np.max(abs(ions.sum(axis=0)/total-1)))
                    if max(error,closure)>1024*np.finfo(float).eps:raise ValueError('Independent LTE population check failed')
                    cache[element]=(part,ions,error,closure)
            if cache[element] is None:dispositions.append(dict(line_id=line_id,reason='PARTITION_COLUMN_DOMAIN'));continue
            part,ions,error,closure=cache[element]
            if not 0<=stage<min(len(ions),len(part.pf)):
                dispositions.append(dict(line_id=line_id,reason='ION_STAGE_UNAVAILABLE'));continue
            lam=float(air_to_vac(float(row['wavelength_air_A'])))/10
            upper=float(row['excitation_potential_eV'])*c.EV+c.HPlanck*c.CLight/(lam*1e-9)
            if upper>=part.ionPot[stage]:dispositions.append(dict(line_id=line_id,reason='UPPER_AT_OR_ABOVE_IONIZATION'));continue
            U=np.exp(np.interp(a.temperature,part.Tpf,part.pf[stage]))
            chi,_=lte_line_coefficients(wave,a.temperature,ions[stage]/U,a.vturb,lambda0_nm=lam,mass_amu=element.mass,
                excitation_eV=float(row['excitation_potential_eV']),loggf=float(row['log_gf']),gamma_s_inv=float(src['turbospectrum_rad']))
            if not np.isfinite(chi).all() or np.any(chi<0):raise ValueError('Invalid line extinction')
            group='modeled_FeI' if item['treatment']=='trace_FeI' else row['species']
            if group not in groups:groups[group]=np.zeros_like(background)
            groups[group]+=chi
            dispositions.append(dict(line_id=line_id,reason='SCREENED_NATURAL_ONLY',species=row['species'],group=group,
                ion_population_fraction_error=error,ion_conservation_error=closure))
            for name,(ix,h) in footprints.items():
                ratio=chi[ix]/background[ix];loc=np.unravel_index(np.argmax(ratio),ratio.shape)
                line_rows.append(dict(case=case['name'],footprint=name,line_id=line_id,species=row['species'],group=group,
                    max_line_to_absorption=float(ratio[loc]),peak_logtau500=float(s['logtau_500'][loc[1]]),
                    peak_observer_air_A=float(observer[ix[loc[0]]]),
                    at_reference_depth_max_ratio=float(np.max(ratio[:,reference_depth]))))
        for group,chi in groups.items():
            for name,(ix,h) in footprints.items():
                mean=h@chi[ix];bg=h@background[ix]
                for depth in range(len(a.temperature)):
                    depth_rows.append(dict(case=case['name'],footprint=name,group=group,depth_index=depth,logtau500=float(s['logtau_500'][depth]),
                        mean_chi_m_inv=float(mean[depth]),mean_absorption_m_inv=float(bg[depth]),mean_ratio=float(mean[depth]/bg[depth])))
        np.testing.assert_array_equal(eq.atomicPops[PeriodicTable['O']].n,pop)
        result['cases'].append(dict(name=case['name'],case=run['case'],run=fingerprint(directory/'report.json'),populations=fingerprint(popfile),
            atmosphere_sources=receipts,atom=fingerprint(atompath),local_binding=bound,
            dispositions=dispositions,reason_counts=dict(__import__('collections').Counter(r['reason'] for r in dispositions)),
            sampled_anchor_count=len(wave),reference_depth_index=reference_depth,reference_logtau500_actual=float(s['logtau_500'][reference_depth]),
            oxygen_populations_unchanged=True))
    result['telemetry']=dict(seconds=time.perf_counter()-start,host=platform.node(),threads=1,load=os.getloadavg(),peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    args.output.mkdir(parents=True,exist_ok=False)
    for name,rows in [('line_ranks.csv',line_rows),('species_depths.csv',depth_rows)]:
        with (args.output/name).open('x') as stream:
            writer=csv.DictWriter(stream,fieldnames=list(rows[0]),lineterminator='\n');writer.writeheader();writer.writerows(rows)
    (args.output/'report.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(telemetry=result['telemetry'],cases=[dict(name=x['name'],reasons=x['reason_counts']) for x in result['cases']])))


if __name__=='__main__':main()
