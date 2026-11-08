"""Deduplicate identical-state records; keep all k<=4 one-photon ground decays."""
import csv
import hashlib
import json
import math
from pathlib import Path
import re
from rya1354_two_method_summary import capsule_files,scf_diagnostics,NUM

ROOT=Path(__file__).resolve().parents[1]
LEVELS={(1,2):('3P2',0.),(1,1):('3P1',158.265),(1,0):('3P0',226.977),
        (2,2):('1D2',15867.862),(2,0):('1S0',33792.583)}


def expected_keys():
    return {(('M' if k%2 else 'E')+str(k),uo,uj,lo,lj)
            for (uo,uj),(_,ue) in LEVELS.items() for (lo,lj),(_,le) in LEVELS.items()
            if ue>le for k in range(1,5) if abs(uj-lj)<=k<=uj+lj}


def parse(text):
    sections=re.split(r'(Electric|Magnetic) 2\*\*\(\s*(\d+)\)-pole transitions',text)
    unique={}
    duplicates={}
    for i in range(1,len(sections),3):
        kind,k,body=sections[i:i+3]
        operator=('E' if kind=='Electric' else 'M')+k
        pattern=(rf'f[12]\s+(\d+)\s+(\d+)\s+\+\s+f[12]\s+(\d+)\s+(\d+)\s+\+\s+([\d.]+)\s+([CM])\s+({NUM})\s+({NUM})\s+({NUM})')
        matches=list(re.finditer(pattern,body))
        if len(matches)!=len(re.findall(r'^\s*f[12]\s+',body,re.M)):
            raise ValueError('Unparsed native forbidden record')
        for m in matches:
            u,uj,l,lj,en,gauge,a,gf,s=m.groups()
            conv=lambda v:float(v.replace('D','E'))
            ac=conv(a)
            following=body[m.end():]
            b=re.match(rf'\s+B\s+({NUM})\s+({NUM})\s+({NUM})',following)
            if kind=='Electric' and b is None:
                raise ValueError('Electric gauge missing')
            ab=conv(b.group(1)) if b else ac
            key=(operator,int(u),int(uj),int(l),int(lj))
            if not all(math.isfinite(v) and v>=0 for v in (ab,ac)) or float(en)<=0:
                raise ValueError('Invalid native forbidden rate')
            row=dict(operator=operator,upper_ordinal=int(u),upper_J=int(uj),lower_ordinal=int(l),lower_J=int(lj),
                     upper=LEVELS[(int(u),int(uj))][0],lower=LEVELS[(int(l),int(lj))][0],
                     native_energy_cm1=float(en),A_B_s1=ab,A_C_s1=ac,
                     gauge_fraction_A=None if ab==ac==0 else abs(ab-ac)/max(ab,ac),
                     computed_zero_in_this_model=(ab==ac==0))
            if key in unique and unique[key]!=row:
                raise ValueError('Duplicate directions disagree; refuse averaging or summation')
            unique[key]=row
            duplicates[key]=duplicates.get(key,0)+1
    if set(unique)!=expected_keys() or any(n!=2 for n in duplicates.values()):
        raise ValueError('Missing/extra physical channel or unexpected duplicate inventory')
    return list(unique.values())


def summarize(destination=None,qualified=False):
    base=ROOT/'data/audit/rya1354'
    capsule=base/'forbidden_ladder_v1.tar.gz'
    files=capsule_files(capsule)
    root='runs/forbidden_ladder_v1'
    receipt=json.loads(files[root+'/receipt.json'])
    if receipt['status']!='FORBIDDEN_MODEL_LADDER_COMPLETED_REVIEW_REQUIRED':
        raise ValueError('Forbidden ladder did not complete')
    controls=base/'forbidden_numerical_controls_v5.tar.gz'
    if qualified:
        files.update(capsule_files(controls))
        resumed=json.loads(files['runs/forbidden_ladder_v5/receipt.json'])
        if resumed['status']!='FORBIDDEN_MODEL_LADDER_COMPLETED_REVIEW_REQUIRED' or resumed['fixed_radial_damping']!=-.9:
            raise ValueError('Explicit numerical continuation did not complete with the declared control')
    references=list(csv.DictReader((base/'nonblind_v1/nist_OI_3000_10000.csv').open()))
    layers=[]
    for n in (3,4,5):
        selected_root='runs/forbidden_ladder_v5' if qualified and n>=4 else root
        hit_limit='MAXIMUM ITERATIONS IN SCF EXCEEDED' in files[f'{selected_root}/n{n}_ground/rmcdhf.stderr'].decode().upper()
        if qualified and hit_limit:
            raise ValueError('Iteration-limit stage cannot enter qualified SCF sequence')
        numerical_status=('ITERATION_LIMIT_HIT_UNQUALIFIED' if hit_limit else
                          'FROZEN_SEED_FROM_UNCONVERGED_N4_UNQUALIFIED' if n==5 and not qualified else 'NO_ITERATION_LIMIT_OBSERVED')
        channels=parse(files[f'{selected_root}/n{n}_multipoles/lower.upper.t'].decode())
        lifetimes={}
        for upper in ('1D2','1S0'):
            decay=[r for r in channels if r['upper']==upper]
            lifetimes[upper]={g:1/sum(r[f'A_{g}_s1'] for r in decay) for g in ('B','C')}
        comparisons=[]
        for row in channels:
            if row['upper']!='1D2' or row['lower'] not in ('3P2','3P1') or row['operator'] not in ('M1','E2'):
                continue
            matches=[r for r in references if r['Type']==row['operator']
                     and r['Lower level'].split('|')[0].strip()=='2s2.2p4'
                     and r['Lower level'].split('|')[1].strip()=='3P'
                     and r['Lower level'].split('|')[2].strip()==str(row['lower_J'])
                     and r['Upper level'].split('|')[0].strip()=='2s2.2p4'
                     and r['Upper level'].split('|')[1].strip()=='1D'
                     and r['Upper level'].split('|')[2].strip()=='2']
            reference,=matches
            comparisons.append(dict(upper=row['upper'],lower=row['lower'],operator=row['operator'],
                                    native_A_B_s1=row['A_B_s1'],reference_A_s1=float(reference['Aki']),
                                    fractional_A_residual=row['A_B_s1']/float(reference['Aki'])-1,
                                    reference_accuracy_letter=reference['Acc.'],reference_TP=reference['TP'],
                                    reference_origin='EVALUATED_THEORY_NOT_ASSUMED_LAB',
                                    uncertainty=None,experimental_accuracy_score=None))
        if len(comparisons)!=4:
            raise ValueError('Separate reference M1/E2 inventory incomplete')
        diagnostics=scf_diagnostics(files[f'{selected_root}/n{n}_ground/rmcdhf.stdout'].decode(),
                                    files[f'{selected_root}/n{n}_ground/rmcdhf.stderr'].decode())
        layers.append(dict(n=n,selected_run=selected_root,numerical_status=numerical_status,
                           native_SCF_diagnostics=diagnostics,channels=channels,model_one_photon_lifetime_s=lifetimes,
                           reference_comparisons=comparisons))
    result=dict(schema='rya1354.forbidden-convergence.v2' if qualified else 'rya1354.forbidden-convergence.v1',scientific_admission='HOLD',
                prediction_origin='THIS_UNPUBLISHED_THEORETICAL_CALCULATION',candidate_doi=None,
                capsule_sha256=hashlib.sha256(capsule.read_bytes()).hexdigest(),layers=layers,
                level_inventory='Five 2p4 levels from retained NIST handbook inventory; same parity, J<=2, all applicable one-photon operators through k=4',
                duplicate_policy='Two equal f1/f2 directions checked and retained once per physical channel; never summed',
                ladder_qualified_for_convergence=False,
                native_stopping_sequence_completed=qualified,
                full_orbital_residual_convergence_validated=False,
                numerical_failure='Original n4 hits the 100-cycle SCF limit despite exit zero; original n5 freezes those unqualified n4 orbitals. Preserve as diagnostics; numerical continuation is a separate run, never a silent replacement.',
                warning='Fixed n<3 orbitals; VV-only Dirac-Coulomb models. Core, Hamiltonian, radial and independent relativistic-method checks and calibrated uncertainty remain owed. Higher-order multiphoton rates not certified.',
                uncertainty=None,adoption_allowed=False)
    if qualified:
        result['numerical_controls_capsule_sha256']=hashlib.sha256(controls.read_bytes()).hexdigest()
        result['numerical_failure']+=' Explicit resumed sequence uses n3 original plus n4/n5 from the completed -0.9 damping-only run; this establishes native stopping only, not full orbital-residual convergence or active-space accuracy.'
    name='forbidden_convergence_v2.json' if qualified else 'forbidden_convergence.json'
    (destination or base/'source_curation_v2'/name).write_text(json.dumps(result,indent=2)+'\n')
    print([(layer['n'],layer['model_one_photon_lifetime_s']) for layer in layers])


if __name__=='__main__':
    import argparse
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--qualified',action='store_true',help='Require the explicit completed numerical-control capsule; never fall back to the unqualified sequence')
    summarize(qualified=ap.parse_args().qualified)
