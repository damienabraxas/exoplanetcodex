"""Portable archive/provenance/label checks; angular products stay on Sirius."""
from decimal import Decimal
import hashlib
import json
import math
from pathlib import Path
import struct
from rya1354_native_parts import load_parts
from rya1354_fac_potential_audit import verify_directory
from rya1354_fac_radial_recouple import qualified_canary
from rya1354_fac_radial_ls_audit import ls_blocks,label_rows,stdout_ids,EXTENDED
from rya1354_fac_recoupling_audit import records


def summarize(base):
    directory=base/'fac_ls_fast_qualified_v37';verify_directory(directory);canary_files,_=load_parts(directory)
    old_cp='runs/fac_ls_fast_canary_v37/';initial_canary=json.loads(canary_files[old_cp+'report.json'])
    directory=base/'fac_ls_fast_corrected_v39';verify_directory(directory);corrected_files,_=load_parts(directory)
    cp='runs/fac_ls_fast_canary_v39/';canary=json.loads(corrected_files[cp+'report.json']);qualified_canary(canary)
    if canary['prior_canary_report_sha256']!=hashlib.sha256(canary_files[old_cp+'report.json']).hexdigest() or initial_canary['negative_controls']['vector_sign']!='Missing LS vector identity':raise ValueError('Canary correction provenance differs')
    prior,_=load_parts(base/'fac_full_angular_complete_v29');old_report=prior['evidence/fac_full_replay_v29b/report.json'];old=json.loads(old_report)
    digest=lambda value:hashlib.sha256(value).hexdigest()
    if digest(old_report)!=canary['reference_sha256'] or canary['converter_receipt_sha256']!=digest(prior['runs/fac_full_recoupling_v29/receipt.json']):raise ValueError('Fast reference provenance differs')
    for name,item in canary['cases'].items():
        reference=old['cases'][name]
        if (item['levels']!=len(reference['levels']) or item['targets']!=reference['targets'] or item['old_orth_diagnostic']!=reference['max_transformation_orthonormality_residual']
            or not math.isfinite(item['new_orth_diagnostic']) or item['signed_rounding_fallbacks']<0 or item['weight_rounding_fallbacks']<0):raise ValueError('Fast replay reference metadata differs')
    directory=base/'fac_radial_angular_first_complete_v38';verify_directory(directory);files,_=load_parts(directory)
    root='runs/fac_radial_angular_first_v38/';receipt=json.loads(files[root+'receipt.json']);report=json.loads(files[root+'report.json'])
    if (receipt['status']!='FAC_RADIAL_ANGULAR_REPLAY_COMPLETED_REVIEW_REQUIRED' or receipt['case']!='ground_coherent10'
        or receipt['canary_report_sha256']!=digest(canary_files[old_cp+'report.json']) or receipt['report_sha256']!=digest(files[root+'report.json'])
        or files[root+'drivers/rya1354_fac_radial_ls_audit.py']!=corrected_files[cp+'rya1354_fac_radial_ls_audit.py'] or corrected_files[cp+'rya1354_fac_radial_ls_audit.py']!=canary_files[old_cp+'rya1354_fac_radial_ls_audit.py']):raise ValueError('Angular replay status/module/canary provenance differs')
    if receipt['scientific_admission']!='HOLD' or receipt['exposure']!='EXPOSED_NONBLIND' or receipt['adoption_allowed'] is not False or any(receipt[k] is not None for k in ('accepted_lifetime','uncertainty','covariance','calibrated_sigma')):raise ValueError('Scientific hold fields differ')
    manifests={name:json.loads((base/bundle/'manifest.json').read_text()) for name,bundle in [('first','fac_radial_ci_first_failed_v33'),('identity','fac_radial_identity_complete_v35'),('remaining','fac_radial_ci_remaining_complete_v36')]}
    for name,path in [('first','runs/fac_radial_ci_first_v33/receipt.json'),('identity','runs/fac_radial_identity_replay_v35/receipt.json'),('remaining','runs/fac_radial_ci_remaining_v36/receipt.json')]:
        if receipt[name+'_receipt_sha256']!=manifests[name]['files'][path]['sha256']:raise ValueError('Qualified CI receipt provenance differs')
    qualified=json.loads((base/'source_curation_v2/fac_radial_ci_v36.json').read_text())
    if receipt['qualified_CI']!=qualified['cases']['ground_coherent10']:raise ValueError('Qualified original CI result differs')
    prefix=root+'ground_coherent10/';native='runs/fac_radial_ci_first_v33/ground_coherent10/'
    for original,dest in [('o_export_000.c','state.c'),('o_export_000.cm','state.cm')]:
        record=manifests['first']['files'][native+original]
        if len(files[prefix+dest])!=record['bytes'] or digest(files[prefix+dest])!=record['sha256']:raise ValueError('Native input transfer differs')
    old_receipt=json.loads(prior['runs/fac_full_recoupling_v29/receipt.json'])
    if receipt['old_converter_receipt_sha256']!=digest(prior['runs/fac_full_recoupling_v29/receipt.json']) or receipt['grasp_pin']!=old_receipt['grasp_pin'] or {x['sha256'] for x in old_receipt['steps']}!={receipt['binary_sha256']} or digest(files[prefix+'isodata'])!=receipt['isodata_sha256'] or files[prefix+'isodata']!=prior['runs/fac_full_recoupling_v29/ground_coherent10/isodata']:raise ValueError('Original converter/support input differs')
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
    return dict(schema='rya1354.fac-radial-angular-metadata.v38',scientific_admission='HOLD',exposure='EXPOSED_NONBLIND',adoption_allowed=False,accepted_lifetime=None,uncertainty=None,covariance=None,calibrated_sigma=None,scope='Independent archive/provenance/native-label/CM metadata checks only; full signed matrix/vector numerical replay performed on Sirius',derived_weight_binary_differences=derived_weight_binary_differences,derived_weight_portability_rule='Native eight-decimal printed precision for derived metadata only; raw Sirius values retained; signed angular replay unchanged',prior_invalid_vector_control_rejection=initial_canary['negative_controls']['vector_sign'],canary=canary,case=receipt['case'],levels=6600,full_LS_coefficients_reconstructed=report['full_LS_coefficients_reconstructed'],label_ID_overflows=report['label_ID_overflows'],targets=targets,target_metadata={term:report['levels'][str(ids[0])] for term,ids in targets.items()},orth_diagnostic=report['max_transformation_orthonormality_residual'],signed_rounding_fallbacks=report['signed_rounding_fallbacks'],weight_rounding_fallbacks=report['weight_rounding_fallbacks'])

if __name__=='__main__':print(json.dumps(summarize(Path(__file__).resolve().parents[1]/'data/audit/rya1354'),indent=2))
