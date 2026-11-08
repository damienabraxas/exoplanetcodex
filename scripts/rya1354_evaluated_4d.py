"""Score physically identified exposed NIST 4d rows without forcing grid aliases.

These evaluated/theoretical references are not laboratory accuracy tests. Use
native computed energies and native gf, with no measured-energy substitution.
"""
import csv
import hashlib
import json
import math
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]


def evaluate(destination=None):
    base=ROOT/'data/audit/rya1354'
    nist=base/'nonblind_v1/nist_OI_3000_10000.csv'
    native=base/'nonblind_v1/native_summary.json'
    rows=list(csv.DictReader(nist.open()))
    predictions=json.loads(native.read_text())['coverage']['optical_4d_E1']
    references=[]
    for row in rows:
        lc,lt,lj=map(str.strip,row['Lower level'].split('|'))
        uc,ut,uj=map(str.strip,row['Upper level'].split('|'))
        if (lc,lt,uc,ut)!=('2s2.2p3.(4S*).3p','5P','2s2.2p3.(4S*).4d','5D*'):
            continue
        if row['Type'] or row['TP']!='T5803,LS':
            raise ValueError('Unexpected multipole or source: review reference scope')
        lower_j,upper_j=int(lj),int(uj)
        gi,gk=map(int,row['gi   gk'].split('-'))
        if (gi,gk)!=(2*lower_j+1,2*upper_j+1):
            raise ValueError('Source J/degeneracy contradiction')
        prediction,=[p for p in predictions if (p['lower_J'],p['upper_J'])==(lower_j,upper_j)]
        observed_lower,observed_upper=map(float,row['Ei           Ek'].split('  -'))
        reference_gf=gi*float(row['fik'])
        if reference_gf<=0 or prediction['gf_B']<=0:
            raise ValueError('Nonpositive oscillator strength')
        references.append(dict(physical_reference_id=f'NIST_OI_2p3-4S-3p_5P_J{lj}__2p3-4S-4d_5Do_J{uj}',
                               reference=row,native_prediction=prediction,
                               reference_gf=reference_gf,reference_origin='EVALUATED_THEORY_HIBBERT1991_LS',
                               delta_loggf_dex=math.log10(prediction['gf_B']/reference_gf),
                               native_energy_residual_cm1=prediction['native_energy_cm1']-(observed_upper-observed_lower),
                               uncertainty=None,covariance=None,experimental_accuracy_score=None,
                               isotope_hfs='UNSPECIFIED',observed_energy_substitution=False,
                               dependence_groups=['4d5D_multiplet','Hibbert1991_common_theory_scale',f'4d5D_J{uj}'],
                               grid_alias_status='NOT_ASSIGNED_TO_AN_INDIVIDUAL_ROUNDED_615X_GRID_ROW'))
    if len(references)!=9 or {(r['native_prediction']['lower_J'],r['native_prediction']['upper_J']) for r in references} != {
            (1,0),(1,1),(1,2),(2,1),(2,2),(2,3),(3,2),(3,3),(3,4)}:
        raise ValueError('Incomplete or duplicate physical multiplet inventory')
    result=dict(schema='rya1354.evaluated-4d-validation.v1',scientific_admission='HOLD',
                prediction_origin='THIS_UNPUBLISHED_THEORETICAL_CALCULATION',candidate_doi=None,
                exposure='NONBLIND_EXTERNAL_THEORY_COMPARISON',
                inputs={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (nist,native)},
                criterion='Source-qualified configurations/parent/term/J/parity and operator, not closest calculated wavelength or gf fit',
                rows=references,experimental_gf_accuracy_rows=0,
                warning='All nine independently identified NIST components scored; rounded grid aliases still unresolved. Accuracy letters are not calibrated sigmas. No stellar grade inferred.')
    (destination or base/'source_curation_v2/evaluated_4d_validation.json').write_text(json.dumps(result,indent=2)+'\n')
    print(f'{len(references)} physical evaluated/theoretical references scored; no laboratory accuracy claim')


if __name__=='__main__':
    evaluate()
