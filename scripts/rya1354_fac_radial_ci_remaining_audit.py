"""Independent complete four-case fixed-potential CI metadata replay, v36."""
import hashlib
import json
from pathlib import Path
import tempfile
from rya1354_native_parts import load_parts
from rya1354_fac_potential_audit import verify_directory
from rya1354_fac_radial_ci_audit import summarize as first_summary
from rya1354_fac_radial_ci_v36 import review,CASES
from rya1354_fac_potential_fields import layout
from rya1354_fac_radial_identity_replay_v35 import source_proof


def summarize(base):
    first=first_summary(base)  # Fully verify old large evidence before loading the new capsule.
    first_manifest=json.loads((base/'fac_radial_ci_first_failed_v33/manifest.json').read_text())
    old_first='runs/fac_radial_ci_first_v33/receipt.json'
    directory=base/'fac_radial_ci_remaining_complete_v36';verify_directory(directory);files,_=load_parts(directory)
    root='runs/fac_radial_ci_remaining_v36/';r=json.loads(files[root+'receipt.json'])
    if r['scientific_admission']!='HOLD' or r['exposure']!='EXPOSED_NONBLIND' or r['adoption_allowed'] is not False or any(r[k] is not None for k in ('accepted_lifetime','uncertainty','covariance','calibrated_sigma')):raise ValueError('Scientific hold fields differ')
    parent,_=load_parts(base/'fac_full_valence_complete_v28');old='runs/fac_full_valence_v28/';pr=json.loads(parent[old+'receipt.json'])
    preflight,_=load_parts(base/'fac_radial_preflight_complete_v32');pfr='runs/fac_radial_preflight_v32/'
    identities,_=load_parts(base/'fac_radial_identity_complete_v35');ir='runs/fac_radial_identity_replay_v35/receipt.json'
    digest=lambda data:hashlib.sha256(data).hexdigest()
    if (r['status']!='FAC_RADIAL_CI_REMAINING_COMPLETED_REVIEW_REQUIRED' or set(r['cases'])!={n for n,_ in CASES}
        or r['baseline_receipt_sha256']!=digest(parent[old+'receipt.json']) or r['preflight_receipt_sha256']!=digest(preflight[pfr+'receipt.json'])
        or r['first_native_receipt_sha256']!=first_manifest['files'][old_first]['sha256'] or r['first_identity_receipt_sha256']!=digest(identities[ir])
        or r['probe']!=pr['probe'] or r['fac_pin']!=pr['fac_pin']):raise ValueError('Remaining CI provenance/scope/status differs')
    source_proof(files[root+'source_documents/faclib/structure.c'].decode())
    spec=layout(files[root+'source_documents/faclib/radial.c'].decode(),files[root+'source_documents/faclib/orbital.c'].decode(),files[root+'source_documents/faclib/consts.h'].decode())
    if spec!=r['potential_layout']:raise ValueError('Potential source layout differs')
    result=dict(schema='rya1354.fac-radial-ci-four-case-audit.v36',scientific_admission='HOLD',exposure='EXPOSED_NONBLIND',adoption_allowed=False,accepted_lifetime=None,uncertainty=None,covariance=None,calibrated_sigma=None,retained_v33_status=first['retained_v33_status'],retained_v34_status=first['retained_v34_status'],cases=dict(first['cases']))
    for case,scale in CASES:
        prefix=root+case+'/'
        if files[prefix+'source.pot']!=parent[old+case+'/o.pot'] or files[prefix+'source.pot']!=preflight[pfr+case+'/source.pot']:raise ValueError('Paired potential differs: '+case)
        with tempfile.TemporaryDirectory(prefix='rya1354-ci-metadata-') as temp:
            out=Path(temp)
            for path,data in files.items():
                if path.startswith(prefix):
                    name=path.removeprefix(prefix)
                    if Path(name).name!=name:raise ValueError('Unexpected nested evidence')
                    (out/name).write_bytes(data)
            reviewed=review(out,spec,scale)
        if reviewed!=r['cases'][case]:raise ValueError('Remaining CI archived replay differs: '+case)
        result['cases'][case]=reviewed
    return result

if __name__=='__main__':print(json.dumps(summarize(Path(__file__).resolve().parents[1]/'data/audit/rya1354'),indent=2))
