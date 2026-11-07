"""Source-qualified inventory of omitted atomic and CNO molecular candidates."""
import argparse,csv,io,json,re,shlex,subprocess,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    from config.constants import require_codex_path
    from pipeline.bride_reference import fingerprint
    from pipeline.molecular_lists import VENDORED_DIR,MANIFEST_PATH,load_manifest
    from pipeline.molecular_canonical import canonical_rows,PROVENANCE_FIELDS
    from scripts.migrate_molecular_transitions_rya1130 import metadata_for
    from pipeline.reflected_solar_rv import SOLAR_GRAVITATIONAL_REDSHIFT
    plan_path=ROOT/'config/experiments/rya1281_missing_opacity.json';plan=json.loads(plan_path.read_text())
    if plan['schema']!='codex.bride_missing_opacity_plan/1':raise ValueError('Wrong missing-opacity plan')
    manifest=load_manifest();engine=require_codex_path('engines.ispec');native=engine/'input/linelists/turbospectrum/molecules'
    result=dict(schema='codex.bride_missing_opacity/1',scientific_admission='HOLD',generator=fingerprint(__file__),plan=fingerprint(plan_path),
        molecular_parser=fingerprint(ROOT/'pipeline/molecular_canonical.py'),manifest=fingerprint(MANIFEST_PATH),
        wrapper=fingerprint(engine/'ispec/synth/turbospectrum.py'),air_evidence=fingerprint(engine/'synthesizer/turbospectrum/source/babsma.f'),
        native_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=engine,text=True).strip(),
        scope=plan['scope'],frame='Stored molecule coordinates interpreted as AIR A by TS consumer; native bsyn AIR documented in babsma/hydropac. This does not promote unknown primary-release frame metadata.',files=[])
    lo,hi=plan['range_consumer_air_A'];molecular=[]
    for molecule in plan['molecules']:
        entry=manifest['molecules'][molecule];metadata=metadata_for(entry)
        for filename in entry['files']:
            match=re.fullmatch(r'(.*)_(\d+)-(\d+)\.bsyn',filename)
            if match is None:raise ValueError('Unrecognized routed source filename')
            lower,upper=map(float,match.groups()[1:])
            if not ((lower>=lo/10 and upper<=hi/10)or(lower<=lo/10<=upper)or(lower<=hi/10<=upper)):continue
            path=VENDORED_DIR/entry['vendored_subdir']/filename;installed=native/filename
            if path.read_bytes()!=installed.read_bytes():raise ValueError('Installed molecular source differs from canonical vendor')
            with path.open()as stream:header=shlex.split(stream.readline());label=stream.readline().strip()
            expected=int(header[-1]);count=0;in_region=0
            for row in canonical_rows(path,molecule,metadata):
                count+=1;wave=row['wavelength']['value_A'];energy=row['lower_energy_eV'];loggf=row['transition_probability']['log_gf']
                if not np.isfinite([wave,energy,loggf]).all()or wave<=0:raise ValueError('Invalid molecular numeric core')
                if not lo<=wave<=hi:continue
                in_region+=1
                molecular.append(dict(transition_id=row['transition_id'],molecule=molecule,isotopologue=row['isotopologue'],
                    filename=filename,consumer_air_A=wave,lower_energy_eV=energy,stored_log_gf=loggf,
                    raw_label=row['transition_identity']['raw_label'],intake_status=row['intake_status'],
                    primary_medium=row['wavelength']['medium'],primary_frame=row['wavelength']['frame'],
                    probability_kind=row['transition_probability']['kind']))
            if count!=expected:raise ValueError('Molecular header/data count differs')
            result['files'].append(dict(molecule=molecule,filename=filename,source=fingerprint(path),installed=fingerprint(installed),
                declared_rows=expected,parsed_rows=count,region_rows=in_region,label=label,
                missing_provenance=[field for field in PROVENANCE_FIELDS if not metadata.get(field)]))
    # Geometry is a membership census, not an identity match or strength calculation.
    geometry=[];c_light_km_s=299792.458
    for rvlabel,rv in [('no_gravity',0.),('plus_gravity',SOLAR_GRAVITATIONAL_REDSHIFT)]:
        factor=np.sqrt((1+rv/c_light_km_s)/(1-rv/c_light_km_s))
        windows=[('O_'+str(c),c-plan['core_half_A'],c+plan['core_half_A'])for c in plan['centres_air_A']]
        centre=plan['feii_air_A']*factor
        windows += [('FeII_blue_shoulder',centre-plan['shoulder_outer_A'],centre-plan['shoulder_inner_A']),
                    ('FeII_red_shoulder',centre+plan['shoulder_inner_A'],centre+plan['shoulder_outer_A'])]
        for molecule in plan['molecules']:
            for label,left,right in windows:
                selected=[r for r in molecular if r['molecule']==molecule and left<=r['consumer_air_A']*factor<=right]
                geometry.append(dict(rv=rvlabel,molecule=molecule,window=label,left_observer_air_A=left,right_observer_air_A=right,
                    source_rows=len(selected),isotopologues=';'.join(sorted({r['isotopologue']for r in selected}))))
    intake_path=ROOT/'data/results/rya1281/anchor_influence_v1/candidate_intake.json'
    intake={r['canonical_row']['line_id']:r for r in json.loads(intake_path.read_text())['rows']}
    screen_path=ROOT/'data/results/rya1281/lte_opacity_screen_v1/report.json';screen=json.loads(screen_path.read_text())
    if screen['intake']!=fingerprint(intake_path):raise ValueError('Changed atomic intake')
    current_text=subprocess.check_output(['git','show','origin/main:data/linelists/canonical_gf.csv'],cwd=ROOT,text=True)
    current={r['line_id']:r for r in csv.DictReader(io.StringIO(current_text))}
    ranks_path=ROOT/'data/results/rya1281/lte_opacity_screen_v1/line_ranks.csv'
    original_screen=require_codex_path('experiments.bride_lte_screen')
    if (original_screen/'report.json').read_bytes()!=screen_path.read_bytes() or (original_screen/'line_ranks.csv').read_bytes()!=ranks_path.read_bytes():raise ValueError('Atomic rank archive differs from immutable source')
    # All original ranks retained separately; do not interpret these coefficient ratios as flux.
    with ranks_path.open()as stream:ranks=list(csv.DictReader(stream))
    aggregated={}
    for row in ranks:
        if row['group']=='modeled_FeI' or row['line_id']=='gf_079010':continue
        key=(row['case'],row['line_id']);value=float(row['at_reference_depth_max_ratio'])
        if key not in aggregated or value>float(aggregated[key]['at_reference_depth_max_ratio']):aggregated[key]=row
    atomic=[]
    for (case,line_id),row in sorted(aggregated.items(),key=lambda x:(x[0][0],-float(x[1]['at_reference_depth_max_ratio']))):
        source=intake[line_id]['canonical_row'];fresh=current.get(line_id)
        atomic.append(dict(case=case,line_id=line_id,species=row['species'],wavelength_air_A=source['wavelength_air_A'],
            excitation_eV=source['excitation_potential_eV'],loggf=source['log_gf'],reference=source['loggf_reference'],
            current_main_matches=fresh==source,strongest_frozen_footprint=row['footprint'],
            at_reference_depth_max_ratio=row['at_reference_depth_max_ratio'],max_column_ratio_same_footprint=row['max_line_to_absorption'],
            scope='Historical natural-LTE local coefficient ratio, not flux/formation/EW prediction or physical-source admission'))
    result.update(atomic_intake=fingerprint(intake_path),atomic_screen=fingerprint(screen_path),atomic_ranks=fingerprint(ranks_path),
        current_main=subprocess.check_output(['git','rev-parse','origin/main'],cwd=ROOT,text=True).strip(),
        molecule_region_rows=len(molecular),molecule_intake_statuses={s:sum(r['intake_status']==s for r in molecular)for s in sorted({r['intake_status']for r in molecular})},
        atomic_rank_rows=len(atomic),molecular_strength_admission='HOLD: primary transition identity/probability, component/isotope weights, partition convention, equilibrium depletion/EOS and widths not established by source-row counts')
    args.output.mkdir(parents=True,exist_ok=False)
    for name,rows in [('molecular_rows.csv',molecular),('geometry.csv',geometry),('atomic_ranks.csv',atomic)]:
        with (args.output/name).open('x')as stream:
            writer=csv.DictWriter(stream,fieldnames=list(rows[0]),lineterminator='\n');writer.writeheader();writer.writerows(rows)
    (args.output/'report.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(files=len(result['files']),molecular_rows=len(molecular),atomic_rows=len(atomic),statuses=result['molecule_intake_statuses'])))


if __name__=='__main__':main()
