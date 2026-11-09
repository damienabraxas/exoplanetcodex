"""Join exact-manifest solar spectra, fresh rate replays and photoionization audits."""
import argparse
import csv
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--analysis',type=Path,required=True)
    parser.add_argument('--rate-reports',type=Path,nargs='+',required=True)
    parser.add_argument('--photo-reports',type=Path,nargs='+',required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    from pipeline.bride_reference import fingerprint
    analysis=json.loads(args.analysis.read_text())
    rates=[(p,json.loads(p.read_text()))for p in args.rate_reports]
    photos=[(p,json.loads(p.read_text()))for p in args.photo_reports]
    cases=sorted([c for c in analysis['cases']if c['experimental_bound_input']is not None],key=lambda c:c['parameters']['teff'])
    if len(cases)!=5 or len(rates)!=5 or len(photos)!=5:raise ValueError('Requires all five declared solar controls/replays/photo audits')
    rows=[];used_rates=set();used_photos=set()
    for case in cases:
        selected=[(p,r)for p,r in rates if r['source_manifest']==case['source_manifest']]
        if len(selected)!=1:raise ValueError('Missing/ambiguous exact source rate replay')
        rate_path,rate=selected[0]
        if rate['status']!='DIAGNOSTIC_COMPLETE':raise ValueError('Incomplete rate replay')
        selected=[(p,r)for p,r in photos if r['source_report']==fingerprint(rate_path)]
        if len(selected)!=1:raise ValueError('Missing/ambiguous exact source photo audit')
        photo_path,photo=selected[0]
        if photo['arrays']!=[m['array']for m in rate['modes']]:raise ValueError('Photo/rate array identities differ')
        top=min(photo['triplet_lower_selected'],key=lambda r:r['logtau500'])
        diag=case['final_diagnostics'];mode=rate['modes'][0]
        row=dict(teff=case['parameters']['teff'],logg=case['parameters']['logg'],abundance=case['parameters']['abundance'],
            iterations=case['iterations'],seconds=case['elapsed_seconds'],effective_cap_length_m=case['effective_cap_length_m'],
            final_population_change=diag['population_change'],final_J_change=diag['J_change'],
            in_iteration_rate_residual=diag['physical_rate_residual'],fresh_rate_residual=mode['max_physical_rate_residual'],
            fresh_radiation_iterations=mode['radiation']['iterations'],top_lower_photo_s_inv=top['independent_photo_s_inv'],
            top_lower_J_weighted_over_B=top['J_weighted_over_B'],independent_photo_over_backend=top['independent_over_backend'],
            lower_top_b=case['top_departures']['lower_triplet'],upper_top_b=case['top_departures']['upper_triplet'],OII_top_b=case['top_departures']['OII_ground'])
        row.update({f'window_{i+1}_area_mA':m['intrinsic_window_area_mA']for i,m in enumerate(case['metrics'])})
        rows.append(row);used_rates.add(str(rate_path));used_photos.add(str(photo_path))
    if len(used_rates)!=5 or len(used_photos)!=5:raise ValueError('Reused source report across temperature cases')
    args.output.mkdir(parents=True,exist_ok=False)
    with (args.output/'summary.csv').open('x')as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]),lineterminator='\n');writer.writeheader();writer.writerows(rows)
    result=dict(schema='codex.bride_solar_diagnostics_join/1',scientific_admission='HOLD',generator=fingerprint(__file__),
        analysis=fingerprint(args.analysis),rate_reports=[fingerprint(p)for p,r in rates],photo_reports=[fingerprint(p)for p,r in photos],
        summary=fingerprint(args.output/'summary.csv'),cases=rows,temperature_steps=analysis['temperature_steps'],
        interpretation='Complete predeclared five-point series; each rate/photo report joined by exact source-manifest and array identity. No fit, best-case selection or observed-data adjustment. Fixed Teff sampling is not physical interpolation accuracy.')
    (args.output/'report.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print('Joined all five exact solar spectra/rate/photo diagnostics')


if __name__=='__main__':main()
