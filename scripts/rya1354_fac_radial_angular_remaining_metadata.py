"""All four expanded LS cases: portable provenance and native-label metadata only."""
from decimal import Decimal
import hashlib
import json
import math
from pathlib import Path
import struct
from rya1354_native_parts import load_parts
from rya1354_fac_potential_audit import verify_directory
from rya1354_fac_radial_angular_metadata import summarize as first_summary
from rya1354_fac_radial_recouple import qualified_canary
from rya1354_fac_radial_recouple_batch import CASES
from rya1354_fac_radial_ls_audit import ls_blocks,label_rows,stdout_ids,EXTENDED
from rya1354_fac_recoupling_audit import records

FIELDS=('levels','full_LS_coefficients_reconstructed','label_ID_overflows','targets','target_metadata','orth_diagnostic','signed_rounding_fallbacks','weight_rounding_fallbacks','derived_weight_binary_differences','derived_weight_portability_rule')


def native_metadata(files,prefix,report):
    cm=records(files[prefix+'state.cm']);dims=struct.unpack('<6i',cm[1])
    if dims!=(8,6600,12,6600,4354454,16):raise ValueError('Actual CM dimensions differ')
    expected=[];energies={};block_names={};ls=ls_blocks(files[prefix+'state.lsj.c'].decode(),EXTENDED)
    for b in range(dims[-1]):
        _,nb,nv,jp,p=struct.unpack('<5i',cm[2+4*b]);parity='+' if p==1 else '-';ids=struct.unpack('<'+'i'*nv,cm[3+4*b]);e=struct.unpack('<'+'d'*(nv+1),cm[4+4*b])
        if len(ls[b])!=nb or len({row['name'] for row in ls[b]})!=nb:raise ValueError('Actual LS block names incomplete/ambiguous')
        names={row['name'] for row in ls[b]}
        for i,ident in enumerate(ids):expected.append((ident,jp-1,parity));energies[ident]=e[0]+e[i+1];block_names[ident]=names
    labels=label_rows(files[prefix+'state.lsj.lbl'].decode(),6600,expected,stdout_ids(files[prefix+'jj2lsj.stdout'].decode()))
    if set(report['levels'])!={str(i) for i in range(6600)} or report['full_LS_coefficients_reconstructed']!=4354454 or report['label_ID_overflows']!=5600:raise ValueError('Complete angular report scope differs')
    derived_weight_binary_differences=[]
    for ident,j,p in expected:
        item=report['levels'][str(ident)];label=labels[ident]
        if (item['two_J'],item['parity'],item['energy_au'])!=(j,p,energies[ident]) or label['energy']!=format(energies[ident],'.9f') or not math.isfinite(energies[ident]):raise ValueError('Original label/CM/report identity differs')
        if (not label['components'] or item['dominant']!=label['components'][0][2] or item['native_component_count']!=len(label['components'])
            or Decimal(format(item['dominant_printed_weight'],'.8f'))!=Decimal(format(float(label['components'][0][0])**2,'.8f')) or item['label_display_id_overflow']!=label['display_overflow']
            or any(name not in block_names[ident] for _,_,name in label['components'])):raise ValueError('Full native composition metadata differs')
        local_weight=float(label['components'][0][0])**2
        if item['dominant_printed_weight']!=local_weight:derived_weight_binary_differences.append(dict(ID=ident,Sirius=item['dominant_printed_weight'],Mac=local_weight,native_eight_decimal_match=True))
    targets={term:[int(i) for i,row in report['levels'].items() if row['dominant']=='1s(2).2s(2).2p(3)4S.3s_'+term] for term in ('3S','5S')}
    if targets!=report['targets'] or any(len(v)!=1 for v in targets.values()):raise ValueError('Target parent/term identity ambiguous')
    return dict(levels=6600,full_LS_coefficients_reconstructed=report['full_LS_coefficients_reconstructed'],label_ID_overflows=report['label_ID_overflows'],targets=targets,target_metadata={term:report['levels'][str(ids[0])] for term,ids in targets.items()},orth_diagnostic=report['max_transformation_orthonormality_residual'],signed_rounding_fallbacks=report['signed_rounding_fallbacks'],weight_rounding_fallbacks=report['weight_rounding_fallbacks'],derived_weight_binary_differences=derived_weight_binary_differences,derived_weight_portability_rule='Native eight-decimal printed precision for derived metadata only; raw Sirius values retained; signed angular replay unchanged')


def hold(record):
    if record['scientific_admission']!='HOLD' or record['exposure']!='EXPOSED_NONBLIND' or record['adoption_allowed'] is not False or any(record[k] is not None for k in ('accepted_lifetime','uncertainty','covariance','calibrated_sigma')):raise ValueError('Scientific hold fields differ')


def summarize(base):
    first=first_summary(base)
    first_manifest=json.loads((base/'fac_radial_angular_first_complete_v38/manifest.json').read_text())
    directory=base/'fac_radial_angular_remaining_complete_v40';verify_directory(directory);files,_=load_parts(directory)
    root='runs/fac_radial_angular_remaining_v40/';batch=json.loads(files[root+'receipt.json']);hold(batch)
    digest=lambda data:hashlib.sha256(data).hexdigest()
    corrected,_=load_parts(base/'fac_ls_fast_corrected_v39');cp='runs/fac_ls_fast_canary_v39/';canary=json.loads(corrected[cp+'report.json']);qualified_canary(canary)
    prior,_=load_parts(base/'fac_full_angular_complete_v29');op='runs/fac_full_recoupling_v29/';old=json.loads(prior[op+'receipt.json']);old_digest=digest(prior[op+'receipt.json']);support=prior[op+'ground_coherent10/isodata'];del prior
    if (batch['status']!='FAC_RADIAL_ANGULAR_REMAINING_COMPLETED_REVIEW_REQUIRED' or set(batch['cases'])!=set(CASES)
        or batch['canary_report_sha256']!=digest(corrected[cp+'report.json'])
        or batch['first_angular_receipt_sha256']!=first_manifest['files']['runs/fac_radial_angular_first_v38/receipt.json']['sha256']
        or files[root+'drivers/rya1354_fac_radial_ls_audit.py']!=corrected[cp+'rya1354_fac_radial_ls_audit.py']):raise ValueError('Angular batch provenance/status/scope differs')
    source_manifest=json.loads((base/'fac_radial_ci_remaining_complete_v36/manifest.json').read_text())
    other_manifests={name:json.loads((base/bundle/'manifest.json').read_text()) for name,bundle in [('first','fac_radial_ci_first_failed_v33'),('identity','fac_radial_identity_complete_v35')]}
    qualified=json.loads((base/'source_curation_v2/fac_radial_ci_v36.json').read_text())
    result=dict(schema='rya1354.fac-radial-angular-four-case-metadata.v40',scientific_admission='HOLD',exposure='EXPOSED_NONBLIND',adoption_allowed=False,accepted_lifetime=None,uncertainty=None,covariance=None,calibrated_sigma=None,scope='Independent complete archive/provenance/native-label/CM/composition metadata only; signed numerical products were verified on Sirius',cases={'ground_coherent10':{k:first[k] for k in FIELDS}})
    for case in CASES:
        child=root+case+'/';receipt=json.loads(files[child+'receipt.json']);report=json.loads(files[child+'report.json']);hold(receipt)
        if (receipt['status']!='FAC_RADIAL_ANGULAR_REPLAY_COMPLETED_REVIEW_REQUIRED' or receipt['case']!=case
            or batch['cases'][case]['receipt_sha256']!=digest(files[child+'receipt.json']) or batch['cases'][case]['report_sha256']!=digest(files[child+'report.json'])
            or receipt['report_sha256']!=digest(files[child+'report.json']) or receipt['canary_report_sha256']!=batch['canary_report_sha256']
            or receipt['qualified_CI']!=qualified['cases'][case] or receipt['old_converter_receipt_sha256']!=old_digest or receipt['grasp_pin']!=old['grasp_pin']
            or {row['sha256'] for row in old['steps']}!={receipt['binary_sha256']} or receipt['python_version']!=canary['python_version'] or receipt['numpy_version']!=canary['numpy_version']
            or files[child+'drivers/rya1354_fac_radial_ls_audit.py']!=corrected[cp+'rya1354_fac_radial_ls_audit.py']):raise ValueError('Angular child source/runtime/qualification differs: '+case)
        refs={'first':(other_manifests['first'],'runs/fac_radial_ci_first_v33/receipt.json'),'identity':(other_manifests['identity'],'runs/fac_radial_identity_replay_v35/receipt.json'),'remaining':(source_manifest,'runs/fac_radial_ci_remaining_v36/receipt.json')}
        if any(receipt[name+'_receipt_sha256']!=manifest['files'][path]['sha256'] for name,(manifest,path) in refs.items()):raise ValueError('Original CI receipt chain differs')
        prefix=child+case+'/'
        for original,dest in [('o_export_000.c','state.c'),('o_export_000.cm','state.cm')]:
            record=source_manifest['files']['runs/fac_radial_ci_remaining_v36/'+case+'/'+original]
            if digest(files[prefix+dest])!=record['sha256'] or len(files[prefix+dest])!=record['bytes']:raise ValueError('Original C/CM transfer differs: '+case)
        if files[prefix+'isodata']!=support or receipt['isodata_sha256']!=digest(support):raise ValueError('Support input differs')
        view=native_metadata(files,prefix,report)
        if any(batch['cases'][case][k]!=view[k] for k in ('targets','full_LS_coefficients_reconstructed','label_ID_overflows')):raise ValueError('Batch result summary differs')
        if any(view['target_metadata'][term]['two_J']!=j or view['target_metadata'][term]['parity']!='-' for term,j in [('3S',2),('5S',4)]):raise ValueError('Target J/parity differs')
        result['cases'][case]=view
    return result

if __name__=='__main__':print(json.dumps(summarize(Path(__file__).resolve().parents[1]/'data/audit/rya1354'),indent=2))
