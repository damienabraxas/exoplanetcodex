"""Hash-checked fixed-orbital sensitivities and a synthetic numerical control."""
import hashlib
import json
from pathlib import Path
import re
from rya1354_two_method_summary import capsule_files,NUM

ROOT=Path(__file__).resolve().parents[1]


def summarize(destination=None):
    base=ROOT/'data/audit/rya1354'
    capsule=base/'rci_sensitivities_v1.tar.gz'
    files=capsule_files(capsule)
    root='runs/rci_sensitivity_v1'
    control='runs/rci_warning_control_v1/DC_dummy_mass_recoil_off'
    if files[f'{root}/DC_E1/lower.upper.ct']!=files[f'{control}_E1/lower.upper.ct']:
        raise ValueError('Synthetic recoil-disabled warning control changes native DC table')
    for tag in ('lower','upper'):
        if 'IEEE_DIVIDE_BY_ZERO' not in files[f'{root}/DC_{tag}/rci.stderr'].decode():
            raise ValueError('Expected stationary-mass numerical diagnostic absent')
        if 'IEEE_DIVIDE_BY_ZERO' in files[f'{control}_{tag}/rci.stderr'].decode():
            raise ValueError('Positive-mass warning control still divides by zero')
    pattern=(rf'f2\s+\d+\s+([123])\s+\+\s+f1\s+1\s+2\s+-\s+([\d.]+)\s+C\s+({NUM})\s+({NUM})\s+({NUM})'
             rf'\s+B\s+({NUM})\s+({NUM})\s+({NUM})')
    cases={}
    for case in ('DC','Breit_static_limit','VP_self_energy','CV_2s_one_hole_DC'):
        table=files[f'{root}/{case}_E1/lower.upper.ct'].decode()
        rows=[]
        for j,en,ac,gc,sc,ab,gb,sb in re.findall(pattern,table):
            conv=lambda v:float(v.replace('D','E'))
            ab,ac,gb,gc=map(conv,(ab,ac,gb,gc))
            rows.append(dict(upper_J=int(j),native_energy_cm1=float(en),A_B_s1=ab,A_C_s1=ac,
                             gf_B=gb,gf_C=gc,gauge_fraction_A=abs(ab-ac)/max(ab,ac)))
        if len(rows)!=3 or {r['upper_J'] for r in rows}!={1,2,3}:
            raise ValueError('Missing or duplicate target components')
        cases[case]=dict(rows=rows,summed_gf_B=sum(r['gf_B'] for r in rows))
    baseline=cases['DC']
    for case,value in cases.items():
        value['fractional_gf_change_from_DC']=value['summed_gf_B']/baseline['summed_gf_B']-1
        value['energy_changes_from_DC_cm1']={str(r['upper_J']):r['native_energy_cm1']-next(
            b['native_energy_cm1'] for b in baseline['rows'] if b['upper_J']==r['upper_J']) for r in value['rows']}
    result=dict(schema='rya1354.rci-sensitivities.v1',scientific_admission='HOLD',
                prediction_origin='THIS_UNPUBLISHED_THEORETICAL_CALCULATION',candidate_doi=None,
                capsule_sha256=hashlib.sha256(capsule.read_bytes()).hexdigest(),cases=cases,
                numerical_warning=dict(observed='IEEE_DIVIDE_BY_ZERO plus underflow/denormal in stationary-nucleus RCI',
                                       source='Pinned setham_gg.f90 line110 evaluates 1/EMN unconditionally; EMN=0. Recoil uses are guarded by disabled LNMS/LSMS.',
                                       positive_control='Synthetic 1 amu input with both recoil terms disabled removes divide-by-zero; full DC E1 table is byte-identical. Underflow/denormal flags remain.',
                                       disposition='No solver patch. Synthetic mass is solely a numerical control, not an oxygen isotope or physical candidate.'),
                interpretation='Coarse n3 fixed-orbital sensitivities; RCI uses a different interpolation grid from MCDHF. Effects are not independent Gaussian errors and are not combined into sigma.',
                uncertainty=None,adoption_allowed=False)
    (destination or base/'source_curation_v2/rci_sensitivities.json').write_text(json.dumps(result,indent=2)+'\n')
    print({k:v['fractional_gf_change_from_DC'] for k,v in cases.items()})


if __name__=='__main__':
    summarize()
