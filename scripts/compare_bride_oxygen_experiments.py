"""Summarize isolated oxygen numerical experiments; never admit an abundance."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))


def compare(paths):
    from pipeline.bride_reference import fingerprint
    records=[]; arrays=[]; spectra=[]
    for path in paths:
        m=json.loads((path/'manifest.json').read_text())
        if m['status']!='CONVERGED_EXPERIMENT':
            raise ValueError(f'Unconverged input: {path}')
        a=np.load(path/'populations.npz')
        b=a['departures']; ref=a['reference_departures']
        if not np.isfinite(b).all() or np.any(b<=0):
            raise ValueError('Invalid departure coefficients')
        tau=np.log10(a['reference_tau500']); selection=(tau>=-3)&(tau<=0)
        if not selection.any():raise ValueError('No diagnostic depth window')
        s=np.genfromtxt(path/'triplet_profiles.csv',delimiter=',',names=True)
        triplet={}
        for key in m['reference_comparison']:
            index=int(key)-1
            triplet[key]=dict(
                calculated_b_range=[float(b[index,selection].min()),float(b[index,selection].max())],
                reference_b_range=[float(ref[index,selection].min()),float(ref[index,selection].max())],
                max_absolute_reference_difference=float(np.max(abs(b[index,selection]-ref[index,selection]))))
        records.append(dict(name=path.name,source_directory=str(path),
            inputs={name:fingerprint(path/name) for name in
                ['manifest.json','history.jsonl','populations.npz','triplet_profiles.csv']},
            atom=m['atom'],atmosphere=m['atmosphere'],fixed_abundance=m['fixed_abundance'],
            controls=m['controls'],iterations=m['iterations'],final_diagnostics=m['final_diagnostics'],
            total_wall_seconds=m['total_wall_seconds'],
            nlte_over_lte_intensity_range=[float(s[s.dtype.names[-1]].min()),float(s[s.dtype.names[-1]].max())],
            triplet_departures_in_diagnostic_window=triplet))
        arrays.append(b);spectra.append(s)
    for r in records[1:]:
        for key in ['atom','atmosphere','fixed_abundance']:
            if r[key]!=records[0][key]:raise ValueError(f'Mismatched experiment {key}')
    differences=[]
    for index in range(1,len(records)):
        s0,s1=spectra[0],spectra[index]
        if s0.dtype.names!=s1.dtype.names or not np.array_equal(s0[s0.dtype.names[0]],s1[s1.dtype.names[0]]):
            raise ValueError('Spectrum wavelength mismatch')
        differences.append(dict(compared=records[index]['name'],baseline=records[0]['name'],
            max_absolute_departure_difference=float(np.max(abs(arrays[index]-arrays[0]))),
            max_relative_intensity_difference={name:float(np.max(abs(s1[name]/s0[name]-1)))
                for name in s0.dtype.names[2:]}))
    return dict(schema='codex.bride_oxygen_numerical_comparison/1',
        scientific_admission='HOLD: numerical experiments, no solar abundance or full-3D result',
        diagnostic_depth_window='-3 <= reference log10(tau500) <= 0; not a calculated contribution function',
        interpretation='Differences are reported, not fitted acceptance budgets. Simultaneous control changes do not isolate individual errors.',
        experiments=records,differences_from_first=differences)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--runs',type=Path,nargs='+',required=True)
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    from config.constants import assert_on_sirius
    assert_on_sirius('Bride oxygen numerical-result comparison')
    result=compare(args.runs)
    with args.output.open('x') as f:json.dump(result,f,indent=2,allow_nan=False);f.write('\n')


if __name__=='__main__':main()
