"""Held-out MARCS column/interpolation diagnostics; no RT/SE or grid adoption."""
import argparse
import csv
import json
from pathlib import Path
import sys

import numpy as np

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    from config.constants import require_codex_path
    from pipeline.bride_reference import read_marcs_structure,fingerprint
    from pipeline.bride_parameter_checks import marcs_filename,interpolate_marcs_structure
    from pipeline.bride_column_checks import integrated_column,candidate_column
    store=require_codex_path('grids.marcs_standard');models={};sources=[];native=[]
    for T in (5500,5750,6000):
        for g in (4.,4.5,5.):
            path=store/marcs_filename(T,g);s=read_marcs_structure(path)
            models[(T,g)]=s;sources.append(dict(teff=T,logg=g,source=fingerprint(path)))
            m=s['column_mass_kg_m2'];check=integrated_column(s['height_m'],s['rho_kg_m3'],m[0])
            native.append(dict(teff=T,logg=g,max_fractional_mass_closure=float(abs(check/m-1).max())))
    nodes=[((T,g),models[(T,g)])for T in (5500,6000)for g in (4.,5.)]
    targets=[(5500,4.5),(5750,4.),(5750,4.5),(5750,5.),(6000,4.5)]
    rows=[];curves=[]
    for target in targets:
        exact=models[target]
        for thermo in ('linear','log_positive','log_positive_layers'):
            s,receipt=interpolate_marcs_structure(nodes,target,thermodynamic_scale=thermo)
            if 'column_mass_kg_m2'in s:raise ValueError('Production interpolation unexpectedly admits columns')
            indices=receipt['simplex_source_indices'];weights=receipt['weights']
            for column_scale in ('linear','log','density_integrated_log_top'):
                m=candidate_column([nodes[i][1]['column_mass_kg_m2']for i in indices],weights,
                    scale='linear'if column_scale=='linear'else 'log')
                if column_scale=='density_integrated_log_top':
                    m=integrated_column(s['height_m'],s['rho_kg_m3'],m[0])
                closure=integrated_column(s['height_m'],s['rho_kg_m3'],m[0])
                rows.append(dict(teff=target[0],logg=target[1],thermodynamic_scale=thermo,column_scale=column_scale,
                    max_temperature_error_K=float(abs(s['temperature_K']-exact['temperature_K']).max()),
                    max_fractional_Pe_error=float(abs(s['electron_pressure_Pa']/exact['electron_pressure_Pa']-1).max()),
                    max_fractional_rho_error=float(abs(s['rho_kg_m3']/exact['rho_kg_m3']-1).max()),
                    max_fractional_column_error=float(abs(m/exact['column_mass_kg_m2']-1).max()),
                    top_fractional_column_error=float(m[0]/exact['column_mass_kg_m2'][0]-1),
                    max_fractional_mass_closure=float(abs(closure/m-1).max()),
                    weights=json.dumps(weights),simplex_indices=json.dumps(indices)))
                for d in range(s['ndep']):
                    curves.append(dict(teff=target[0],logg=target[1],thermodynamic_scale=thermo,column_scale=column_scale,
                        depth_1based=d+1,logtauRoss=exact['logtau_ross'][d],native_column_kg_m2=exact['column_mass_kg_m2'][d],
                        candidate_column_kg_m2=m[d],integrated_density_column_kg_m2=closure[d]))
    args.output.mkdir(parents=True,exist_ok=False)
    result=dict(schema='codex.bride_column_interpolation_audit/1',scientific_admission='HOLD',
        generator=fingerprint(__file__),column_checks=fingerprint(ROOT/'pipeline/bride_column_checks.py'),
        interpolator=fingerprint(ROOT/'pipeline/bride_parameter_checks.py'),sources=sources,native_closure=native,
        method='Four5500/6000K,g4/5corners; five distinct native held-out targets; all3 existing thermodynamic and3 candidate column conventions retained. The integrated candidate enforces mass closure by construction; only held-out native error independently tests its accuracy. No hydrostatic/radiative equilibrium certification.',
        production_column_admission=False,cases=rows)
    for name,data in [('summary.csv',rows),('columns.csv',curves)]:
        with (args.output/name).open('x')as stream:
            writer=csv.DictWriter(stream,fieldnames=list(data[0]),lineterminator='\n');writer.writeheader();writer.writerows(data)
        result[name]=fingerprint(args.output/name)
    (args.output/'report.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps([r for r in rows if (r['teff'],r['logg'])==(5750,4.5)],indent=2))


if __name__=='__main__':main()
