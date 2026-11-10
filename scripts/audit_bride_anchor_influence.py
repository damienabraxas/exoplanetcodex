"""Account for fixed-support continuum differences and source availability."""
import argparse,csv,json,sys
from collections import Counter
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    from config.constants import require_codex_path
    from pipeline.bride_reference import fingerprint
    from pipeline.bride_anchor_influence import fixed_support_blocks
    import matplotlib
    matplotlib.use('Agg');import matplotlib.pyplot as plt
    plan_path=ROOT/'config/experiments/rya1281_anchor_influence.json';plan=json.loads(plan_path.read_text())
    if plan['schema']!='codex.bride_anchor_influence_plan/1':raise ValueError('Unknown plan')
    source_path=ROOT/plan['source_report'];source=json.loads(source_path.read_text())
    pixels_path=ROOT/plan['derived_pixels']
    if fingerprint(pixels_path)!=source['derived_pixels']:raise ValueError('Changed derived pixel input')
    packet_path=ROOT/plan['trace_packet'];packet=json.loads(packet_path.read_text())
    canonical=ROOT/'data/linelists/canonical_gf.csv'
    original=require_codex_path('engines.ispec')/'input/linelists/transitions/GESv6_atom_hfs_iso.420_920nm/atomic_lines.tsv'
    if fingerprint(canonical)!=packet['canonical'] or fingerprint(original)!=packet['source']:raise ValueError('Changed canonical/source intake')
    with canonical.open() as stream:candidates=[r for r in csv.DictReader(stream) if 7763<=float(r['wavelength_air_A'])<=7783]
    with original.open() as stream:original_rows=[r for r in csv.DictReader(stream,delimiter='\t') if 7763<=float(r['wave_A'])<=7783]
    included={r['line_id'] for r in packet['lines']};excluded={r['canonical_row']['line_id'] for r in packet['excluded']}
    if included!={r['canonical_row']['line_id'] for r in packet['evidence']}:raise ValueError('Trace identity mismatch')
    intake=[]
    for row in candidates:
        symbol=row['species'].split()[0];element=symbol+' '+str(int(float(row['ion'])))
        matches=[r for r in original_rows if r['element']==element
                 and abs(float(r['wave_A'])-float(row['wavelength_air_A']))<1e-6
                 and abs(float(r['lower_state_eV'])-float(row['excitation_potential_eV']))<1e-6]
        treatment=('trace_FeI' if row['line_id'] in included else 'excluded_combined_FeI' if row['line_id'] in excluded
                   else 'active_OI_candidate' if row['species']=='O I' else 'not_in_trace_subset')
        intake.append(dict(canonical_row=row,treatment=treatment,exact_identity_matches=len(matches),
            source_rows=[{k:r[k] for k in ('element','wave_A','lower_state_eV','loggf','rad','turbospectrum_rad','stark','waals')} for r in matches],
            canonical_gf_matches_any_source=any(abs(float(r['loggf'])-float(row['log_gf']))<1e-6 for r in matches)))
    with pixels_path.open() as stream:pixels=list(csv.DictReader(stream))
    result=dict(schema='codex.bride_anchor_influence/1',scientific_admission='HOLD',
        generator=fingerprint(__file__),adapter=fingerprint(ROOT/'pipeline/bride_anchor_influence.py'),
        plan=fingerprint(plan_path),source_report=fingerprint(source_path),derived_pixels=fingerprint(pixels_path),
        trace_packet=fingerprint(packet_path),canonical=fingerprint(canonical),original=fingerprint(original),
        bin_width_air_A=plan['bin_width_air_A'],candidate_accounting=dict(Counter(r['treatment'] for r in intake)),
        source_accounting=dict(Counter(str(r['exact_identity_matches']) for r in intake)),
        caveats='Fixed-support linear centre accounting, not independent errors or a physical opacity attribution. Candidate centres within bins are an inventory only, not blend identities, wing coverage, predicted strengths or source/abundance admission. Original frozen spectra/reference choices unchanged.',cases=[])
    blocks_csv=[];figures={}
    for case in source['cases']:
        out=dict(name=case['name'],physical_case=case['physical_case'],lines=[])
        fig,axes=plt.subplots(3,2,figsize=(14,11))
        for li,line in enumerate(case['lines']):
            centre=line['centre_air_A'];refs={}
            for ref in ('synthesis','iag'):
                rr=[r for r in pixels if r['case']==case['name'] and float(r['centre_air_A'])==centre and r['reference']==ref]
                refs[ref]=sorted(rr,key=lambda r:float(r['offset_air_A']))
            offsets=np.array([float(r['offset_air_A']) for r in refs['synthesis']])
            np.testing.assert_array_equal(offsets,[float(r['offset_air_A']) for r in refs['iag']])
            ratios={m:np.array([float(r['relative_ratio']) for r in rr]) for m,rr in refs.items()}
            lo=dict(centre_air_A=centre,masks={})
            for mi,maskref in enumerate(('synthesis','iag')):
                rr=refs[maskref];mask=np.array([r['coefficient_support']=='True' for r in rr])
                x=offsets[mask];h=np.array([float(r['centre_weight']) for r in rr])[mask]
                a,b=ratios['synthesis'][mask],ratios['iag'][mask]
                analysis=fixed_support_blocks(x,a,b,h,plan['bin_width_air_A'])
                fixed=line['frozen_support'][maskref]
                expected=fixed['iag']['level']-fixed['synthesis']['level']
                independent=float(np.polyfit(x,b,1)[-1]-np.polyfit(x,a,1)[-1])
                bound=128*np.finfo(float).eps*max(1.,float(np.sum(abs(h)))*max(float(np.max(abs(a))),float(np.max(abs(b)))))
                if abs(analysis['centre_delta']-expected)>bound or abs(independent-expected)>bound:raise ValueError('Fixed-support reconstruction differs')
                analysis.update(expected_level_delta=expected,independent_fit_delta=independent,roundoff_bound=bound,
                    coefficient_support_count=int(mask.sum()),level_reconstruction_error=analysis['centre_delta']-expected,
                    observed_area_delta_mA=fixed['iag']['observed_area_mA']-fixed['synthesis']['observed_area_mA'])
                for block in analysis['blocks']:
                    block['lo_air_A']=centre+block['lo_offset_A'];block['hi_air_A']=centre+block['hi_offset_A']
                    nearby=[r for r in intake if block['lo_air_A']<=float(r['canonical_row']['wavelength_air_A'])<block['hi_air_A']]
                    block['candidate_line_ids']=[r['canonical_row']['line_id'] for r in nearby]
                    block['candidate_treatments']=dict(Counter(r['treatment'] for r in nearby))
                    blocks_csv.append(dict(case=case['name'],centre_air_A=centre,mask_reference=maskref,
                        **{k:v for k,v in block.items() if k not in ('candidate_line_ids','candidate_treatments')},
                        candidate_line_ids=';'.join(block['candidate_line_ids'])))
                lo['masks'][maskref]=analysis
                top=sorted(analysis['blocks'],key=lambda r:abs(r['centre_delta']),reverse=True)[:12]
                ax=axes[li,mi];values=[r['centre_delta'] for r in top]
                ax.bar(np.arange(len(top)),values,color=['tab:blue' if v>=0 else 'tab:orange' for v in values])
                ax.set_xticks(np.arange(len(top)),[f"{r['lo_air_A']:.2f}" for r in top],rotation=70,fontsize=7)
                ax.set(title=f'{centre:.3f} AIR A; frozen {maskref} support',ylabel='Signed centre-level contribution',xlabel='0.1A bin lower AIR edge; 12 largest absolute contributions')
                ax.axhline(0,color='grey',lw=.7);ax.grid(axis='y',alpha=.2)
            out['lines'].append(lo)
        fig.suptitle(case['name']+': fixed-mask IAG-minus-synthesis ratio accounting; HOLD');fig.tight_layout();figures[case['name']]=fig;result['cases'].append(out)
    args.output.mkdir(parents=True,exist_ok=False)
    for name,fig in figures.items():fig.savefig(args.output/(name+'_influence.png'),dpi=140);plt.close(fig)
    with (args.output/'blocks.csv').open('x') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(blocks_csv[0]),lineterminator='\n');writer.writeheader();writer.writerows(blocks_csv)
    (args.output/'candidate_intake.json').write_text(json.dumps(dict(scope=result['caveats'],rows=intake),indent=2,allow_nan=False)+'\n')
    (args.output/'report.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(candidate_accounting=result['candidate_accounting'],source_accounting=result['source_accounting'])))


if __name__=='__main__':main()
