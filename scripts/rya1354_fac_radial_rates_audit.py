"""Portable full-evidence provenance and independent radial decay-rate replay."""
import hashlib
import json
from pathlib import Path

from rya1354_native_parts import load_parts
from rya1354_fac_potential_audit import verify_directory
from rya1354_fac_radial_rates import rate_input
from rya1354_fac_radial_recouple import CASES
from rya1354_fac_radial_angular_remaining_metadata import hold
from rya1354_fac_no_correction_audit import configuration_keys,verify_config_inventory
from rya1354_fac_expanded_rate_audit import review as rate_review
from rya1354_fac_rate_audit import constants


def basis_keys(data):
    marker=b'Mixing Coefficients';end=data.find(marker)
    if end<0:raise ValueError('Native basis coefficient boundary missing')
    # Retain the recognized heading, including its leading native decoration.
    # Cutting immediately before the words leaves a bogus separator-only row.
    return configuration_keys(data[:end+len(marker)].decode())


def summarize(base):
    directory=base/'fac_radial_rates_complete_v41';verify_directory(directory)
    files,_=load_parts(directory);root='runs/fac_radial_rates_v41/'
    receipt=json.loads(files[root+'receipt.json']);hold(receipt)
    index_raw=files[root+'index.json'];index=json.loads(index_raw);hold(index)
    plan=json.loads((base.parents[2]/'config/atomic/rya1354/lifetime_plan.v41.json').read_text())
    digest=lambda data:hashlib.sha256(data).hexdigest()
    if (receipt['status']!='FAC_RADIAL_RATES_COMPLETED_REVIEW_REQUIRED' or set(index['cases'])!=set(CASES)
        or digest(index_raw)!=receipt['index_sha256'] or digest(index_raw)!=plan['index_sha256']
        or receipt['probe']!=index['probe'] or receipt['fac_pin']!=index['fac_pin']):
        raise ValueError('Radial rates/index not qualified')
    bundles={key:json.loads((base/name/'manifest.json').read_text()) for key,name in (
        ('first','fac_radial_ci_first_failed_v33'),('identity','fac_radial_identity_complete_v35'),
        ('remaining','fac_radial_ci_remaining_complete_v36'),('angular_first','fac_radial_angular_first_complete_v38'),
        ('angular_remaining','fac_radial_angular_remaining_complete_v40'))}
    for key,record in index['parent_receipts'].items():
        if bundles[key]['files'][record['path']]['sha256']!=record['sha256']:
            raise ValueError('Original parent receipt differs')
    if set(index['parent_receipts'])!=set(bundles):raise ValueError('Incomplete original parent chain')
    c=constants(files[root+'source_documents/consts.h'].decode())
    original_constants=bundles['first']['files']['runs/fac_radial_ci_first_v33/source_documents/faclib/consts.h']
    if digest(files[root+'source_documents/consts.h'])!=original_constants['sha256']:
        raise ValueError('Pinned source constants differ')
    expected={case+'_'+g for case in CASES for g in ('B','C')}
    if set(receipt['cases'])!=expected:raise ValueError('Missing or extra radial gauge case')
    result=dict(schema='rya1354.fac-radial-rates-audit.v41',scientific_admission='HOLD',exposure='EXPOSED_NONBLIND',
        adoption_allowed=False,accepted_lifetime=None,uncertainty=None,covariance=None,calibrated_sigma=None,cases={})
    for case in CASES:
        entry=index['cases'][case];mapping=entry['mapping'];scale=100 if case.endswith('100') else 10
        first=case=='ground_coherent10';source=bundles['first' if first else 'remaining']['files']
        angular=bundles['angular_first' if first else 'angular_remaining']['files']
        if (angular[entry['angular']+'/report.json']['sha256']!=entry['report_sha256']
            or angular[entry['angular']+'/receipt.json']['sha256']!=entry['receipt_sha256']):
            raise ValueError('Case-specific mapping hash differs')
        # Index contains the complete report, so also check its exact JSON bytes.
        report_bytes=(json.dumps(mapping,indent=2)+'\n').encode()
        if digest(report_bytes)!=entry['report_sha256']:raise ValueError('Embedded full LS report differs')
        for suffix,gauge in [('B',2),('C',1)]:
            name=case+'_'+suffix;p=root+name+'/'
            text,requests,lower=rate_input(scale,gauge,mapping)
            if files[p+'input.sf'].decode()!=text:raise ValueError('Declared radial rate input differs')
            for fn in ('o_export_000.c','o_export_000.cm','source.pot','ci.pot'):
                record=source[entry['source']+'/'+fn]
                if digest(files[p+fn])!=record['sha256'] or len(files[p+fn])!=record['bytes']:
                    raise ValueError('Original state/potential differs: '+name+'/'+fn)
            keys=basis_keys(files[p+'o.basis'])
            configs={n:verify_config_inventory(files[p+'configs.'+n+'.tsv'].decode(),keys) for n in ('before','after')}
            qualified=entry['qualified_CI']
            if configs!=qualified['configuration_state'] or receipt['cases'][name]['qualified_CI']!=qualified:
                raise ValueError('Original zero-correction/radial state differs')
            audit=rate_review(files[p+'raw.tsv'].decode(),files[p+'rates.tr'].decode(),requests,gauge,mapping,c)
            if (audit!=receipt['cases'][name]['audit'] or requests!=receipt['cases'][name]['requests']
                or lower!=receipt['cases'][name]['lower_states'] or gauge!=receipt['cases'][name]['gauge']):
                raise ValueError('Independent raw/table/coverage replay differs')
            result['cases'][name]=dict(audit=audit,lower_states=lower,configuration_state=configs)
    previous_path=base/'source_curation_v2/fac_full_rates_v30.json'
    previous=json.loads(previous_path.read_text())
    result['previous_report_sha256']=digest(previous_path.read_bytes())
    result['paired_n3_radial_comparison']={case:{g:dict(
        n3_lifetime_ns=previous['cases'][case+'_'+g]['audit']['exploratory_model_lifetime_ns'],
        radial_lifetime_ns=result['cases'][case+'_'+g]['audit']['exploratory_model_lifetime_ns'])
        for g in ('B','C')} for case in CASES}
    return result


if __name__=='__main__':print(json.dumps(summarize(Path(__file__).resolve().parents[1]/'data/audit/rya1354'),indent=2))
