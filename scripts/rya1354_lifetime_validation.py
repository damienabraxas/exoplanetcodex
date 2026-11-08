"""Hash-checked exploratory 3s 3S1 model lifetime versus one primary experiment.

No sigma, chi-square or adopted uncertainty is inferred. All one-photon angular
channels in the six-lower-level inventory are retained, including computed zeros.
"""
from pathlib import Path
import hashlib
import json
import math
import re
import tarfile

ROOT=Path(__file__).resolve().parents[1]
NUM=r'[\d.D+\-]+'


def parse_table(text,destination):
    sections=re.split(r'(Electric|Magnetic) 2\*\*\(\s*(\d+)\)-pole transitions',text)
    rows=[]
    for i in range(1,len(sections),3):
        kind,k,body=sections[i:i+3]
        operator=('E' if kind=='Electric' else 'M')+k
        pattern=(rf'f2\s+(\d+)\s+(\d+)\s+-\s+f1\s+(\d+)\s+(\d+)\s+([+-])\s+([\d.]+)\s+([CM])\s+({NUM})\s+({NUM})\s+({NUM})')
        matches=list(re.finditer(pattern,body))
        if len(matches)!=len(re.findall(r'^\s*f[12]\s+',body,re.M)):
            raise ValueError('Unexpected or unparsed native transition')
        for match in matches:
            u,uj,l,lj,parity,energy,gauge,a,gf,s=match.groups()
            conv=lambda v:float(v.replace('D','E'))
            ac=conv(a)
            following=body[match.end():]
            b=re.match(rf'\s+B\s+({NUM})\s+({NUM})\s+({NUM})',following)
            if kind=='Electric' and b is None:
                raise ValueError('Electric row missing second native gauge')
            ab=conv(b.group(1)) if b else ac
            if u!='1' or uj!='1' or not all(math.isfinite(v) and v>=0 for v in (ab,ac)) or float(energy)<=0:
                raise ValueError('Wrong upper state, nonfinite/negative rate or downward energy')
            rows.append(dict(destination=destination,operator=operator,lower_ordinal=int(l),lower_J=int(lj),
                             lower_parity=parity,native_energy_cm1=float(energy),
                             A_B_s1=ab,A_C_s1=ac,magnetic_gauge='NATIVE_M' if kind=='Magnetic' else None,
                             computed_zero_in_this_model=(ab==ac==0)))
    return rows


def model_lifetimes(rows):
    expected={('to_ground','E1',1,0),('to_ground','E1',2,0),('to_ground','E1',1,1),
              ('to_ground','E1',1,2),('to_ground','E1',2,2),
              ('to_ground','M2',1,1),('to_ground','M2',1,2),('to_ground','M2',2,2),
              ('to_ground','E3',1,2),('to_ground','E3',2,2),
              ('to_5S','M1',1,2),('to_5S','E2',1,2),('to_5S','M3',1,2)}
    keys=[(r['destination'],r['operator'],r['lower_ordinal'],r['lower_J']) for r in rows]
    if len(keys)!=len(set(keys)) or set(keys)!=expected:
        raise ValueError('Missing, extra or duplicate decay channel; refuse lifetime')
    result={}
    for gauge in ('B','C'):
        total=sum(r[f'A_{gauge}_s1'] for r in rows)
        if not math.isfinite(total) or total<=0:
            raise ValueError('Invalid decay total')
        result[gauge]=dict(total_A_s1=total,model_one_photon_lifetime_ns=1e9/total)
        for row in rows:
            row[f'model_branching_fraction_{gauge}']=row[f'A_{gauge}_s1']/total
    return result


def evaluate(destination=None):
    base=ROOT/'data/audit/rya1354'
    capsule=base/'radial_decay_controls_v1.tar.gz'
    with tarfile.open(capsule) as archive:
        files={m.name:archive.extractfile(m).read() for m in archive.getmembers() if m.isfile()}
    root='runs/decay_inventory_v1'
    receipt=json.loads(files[root+'/receipt.json'])
    if receipt['status']!='MODEL_SPACE_OPERATORS_COMPLETED_NOT_ACCURACY_VALIDATED':
        raise ValueError('Native decay run did not complete')
    for relative,record in receipt['files'].items():
        if hashlib.sha256(files[root+'/'+relative]).hexdigest()!=record['sha256']:
            raise ValueError('Retained artifact hash mismatch')
    rows=[]
    for case in ('to_ground','to_5S'):
        rows+=parse_table(files[f'{root}/{case}/lower.upper.t'].decode(),case)
    baseline=json.loads((base/'nonblind_v1/native_summary.json').read_text())['coverage']['UV_E1']
    if sorted(r['A_B_s1'] for r in rows if r['operator']=='E1')!=sorted(r['A_B_s1'] for r in baseline):
        raise ValueError('Fixed-wavefunction E1 positive control changed')
    model=model_lifetimes(rows)
    source_receipt=base/'source_curation_v2/source_receipt.json'
    level_receipt=base/'source_curation_v2/nist_low_levels.raw.html'
    curation=json.loads(source_receipt.read_text())
    primary,=[s for s in curation['sources'] if s.get('archive_relative_path','').endswith('smith1971.pdf')]
    experiment=dict(primary_source='Smith et al. 1971 ApJ165 217, Table1 first row and Table2 present work',
                    primary_pdf_sha256=primary['sha256'],primary_url=primary['url'],
                    upper='O I 2s2 2p3(4S)3s 3S(o), J=1',quantity='lifetime',value_ns=1.79,
                    reported_uncertainty_ns=.17,method='83 keV beam-foil; velocity measured after foil; cascade-corrected two-exponential fit',
                    uncertainty_convention='Reported combination of statistical rms fitting errors and beam-velocity error (4%); combination rule/confidence not specified sufficiently for calibrated Gaussian scoring',
                    isotope='NOT_EXPLICITLY_IDENTIFIED',
                    dependence='One upper lifetime; shared velocity calibration and cascade/background fit. Brooks1977 1302A is borrowed Smith1971, not independent. Recommended f=0.048±0.003 combines other studies and is not an additional observation.',
                    exposure='EXPOSED_NONBLIND; value already visible in Li2023 comparison before primary acquisition',
                    review='PRIMARY_CONTENT_INSPECTED; DISTINCT_CURATOR_AND_HUMAN_REVIEW_OWED')
    result=dict(schema='rya1354.nonblind-lifetime-validation.v1',scientific_admission='HOLD',
                prediction_origin='THIS_UNPUBLISHED_THEORETICAL_CALCULATION',candidate_doi=None,
                capsule_sha256=hashlib.sha256(capsule.read_bytes()).hexdigest(),
                source_receipts={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
                                 for p in (source_receipt,level_receipt)},
                scope='EXPLORATORY_MODEL_ONE_PHOTON_RATE_INVENTORY_NOT_CALIBRATED_PHYSICAL_LIFETIME',
                lower_level_inventory='NIST handbook: ground 3P2/1/0, 1D2, 1S0, 3s5So2, all below 3s3So1; k<=3 by angular triangle',
                excluded_physics='Configuration/radial/Hamiltonian convergence and higher-order multiphoton contributions not validated',
                model=model,channels=rows,experiment=experiment,
                raw_residual_ns={g:model[g]['model_one_photon_lifetime_ns']-experiment['value_ns'] for g in ('B','C')},
                calibrated_sigma=None,chi_square=None,adoption_allowed=False)
    (destination or base/'source_curation_v2/lifetime_validation.json').write_text(json.dumps(result,indent=2)+'\n')
    print(model)


if __name__=='__main__':
    evaluate()
