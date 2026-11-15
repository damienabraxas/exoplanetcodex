"""Independent archived radial readiness replay; no atomic solver on Mac."""
import hashlib
import json
from pathlib import Path
from rya1354_native_parts import load_parts
from rya1354_fac_potential_audit import verify_directory
from rya1354_fac_potential_fields import layout
from rya1354_fac_radial_preflight import review


def summarize(base):
    directory=base/'fac_radial_preflight_complete_v32';verify_directory(directory);files,_=load_parts(directory);root='runs/fac_radial_preflight_v32/'
    parent,_=load_parts(base/'fac_full_valence_complete_v28');receipt=json.loads(files[root+'receipt.json'])
    if receipt['status']!='FAC_RADIAL_PREFLIGHT_COMPLETED_REVIEW_REQUIRED' or receipt['baseline_receipt_sha256']!=hashlib.sha256(parent['runs/fac_full_valence_v28/receipt.json']).hexdigest():raise ValueError('Radial preflight provenance/incomplete status differs')
    if receipt['binary_sha256']!=json.loads(parent['runs/fac_full_valence_v28/receipt.json'])['probe']['upstream_binary_sha256']:raise ValueError('Radial preflight native binary differs')
    spec=layout(files[root+'source_documents/faclib/radial.c'].decode(),files[root+'source_documents/faclib/orbital.c'].decode(),files[root+'source_documents/faclib/consts.h'].decode())
    if spec!=receipt['potential_layout']:raise ValueError('Radial layout receipt differs')
    result=dict(schema='rya1354.fac-radial-preflight-audit.v32',scientific_admission='HOLD',exposure='EXPOSED_NONBLIND',adoption_allowed=False,accepted_lifetime=None,uncertainty=None,covariance=None,calibrated_sigma=None,cases={})
    for policy in ('ground','equal'):
        for label,scale in [('coherent10',10),('coherent100',100)]:
            name=policy+'_'+label;p=root+name+'/';data={k.removeprefix(p):v for k,v in files.items() if k.startswith(p)}
            if data['source.pot']!=parent['runs/fac_full_valence_v28/'+name+'/o.pot']:raise ValueError('Radial probe substituted paired potential')
            reviewed=review(data,scale,spec)
            if reviewed!=receipt['cases'][name]:raise ValueError('Archived radial replay differs from native receipt')
            result['cases'][name]=reviewed
    return result

if __name__=='__main__':print(json.dumps(summarize(Path(__file__).resolve().parents[1]/'data/audit/rya1354'),indent=2))
