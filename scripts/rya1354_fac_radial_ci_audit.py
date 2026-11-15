"""Independent lossless v33/v34/v35 identity replay; no atomic solver."""
import hashlib
import json
from pathlib import Path
import tempfile
from rya1354_native_parts import load_parts
from rya1354_fac_potential_audit import verify_directory
from rya1354_fac_radial_ci_v35 import review
from rya1354_fac_potential_fields import layout
from rya1354_fac_radial_identity_replay_v35 import source_proof


def summarize(base):
    directory=base/'fac_radial_ci_first_failed_v33';verify_directory(directory);files,_=load_parts(directory)
    directory=base/'fac_radial_identity_complete_v35';verify_directory(directory);replays,_=load_parts(directory)
    root='runs/fac_radial_ci_first_v33/';native=json.loads(files[root+'receipt.json'])
    v34='runs/fac_radial_identity_replay_v34/';v35='runs/fac_radial_identity_replay_v35/'
    previous=json.loads(replays[v34+'receipt.json']);r=json.loads(replays[v35+'receipt.json'])
    parent,_=load_parts(base/'fac_full_valence_complete_v28');old='runs/fac_full_valence_v28/';pr=json.loads(parent[old+'receipt.json'])
    preflight,_=load_parts(base/'fac_radial_preflight_complete_v32');pfr='runs/fac_radial_preflight_v32/'
    digest=lambda data:hashlib.sha256(data).hexdigest()
    if native['status']!='FAILED_RETAINED' or native['error']!='ValueError: Native coupled basis duplicate/identity differs' or previous['status']!='FAILED_RETAINED' or previous['error']!='ValueError: CSF identity permutation differs':raise ValueError('Retained failure status differs')
    if r['status']!='FAC_RADIAL_IDENTITY_REPLAY_COMPLETED_REVIEW_REQUIRED' or r['source_native_receipt_sha256']!=digest(files[root+'receipt.json']) or r['previous_replay_receipt_sha256']!=digest(replays[v34+'receipt.json']) or previous['source_native_receipt_sha256']!=digest(files[root+'receipt.json']):raise ValueError('Replay chain provenance differs')
    if native['baseline_receipt_sha256']!=digest(parent[old+'receipt.json']) or native['preflight_receipt_sha256']!=digest(preflight[pfr+'receipt.json']) or native['probe']!=pr['probe']:raise ValueError('CI provenance differs')
    if r['native_recalculation'] or set(r['cases'])!={'ground_coherent10'} or r['canaries']!={k:v['transport'] for k,v in pr['cases'].items()}:raise ValueError('Replay scope/canaries differ')
    if source_proof(files[root+'source_documents/faclib/structure.c'].decode())!=r['source_proof']:raise ValueError('Primary identity source proof differs')
    spec=layout(files[root+'source_documents/faclib/radial.c'].decode(),files[root+'source_documents/faclib/orbital.c'].decode(),files[root+'source_documents/faclib/consts.h'].decode())
    if spec!=native['potential_layout']:raise ValueError('Potential source layout differs')
    prefix=root+'ground_coherent10/'
    if files[prefix+'source.pot']!=parent[old+'ground_coherent10/o.pot'] or files[prefix+'source.pot']!=preflight[pfr+'ground_coherent10/source.pot']:raise ValueError('Paired potential differs')
    with tempfile.TemporaryDirectory(prefix='rya1354-ci-metadata-') as temp:
        out=Path(temp)
        for path,data in files.items():
            if path.startswith(prefix):
                name=path.removeprefix(prefix)
                if Path(name).name!=name:raise ValueError('Unexpected nested evidence')
                (out/name).write_bytes(data)
        result=review(out,spec)
    if result!=r['cases']['ground_coherent10']:raise ValueError('CI archived replay differs')
    return dict(schema='rya1354.fac-radial-ci-audit.v35',scientific_admission='HOLD',exposure='EXPOSED_NONBLIND',adoption_allowed=False,accepted_lifetime=None,uncertainty=None,covariance=None,calibrated_sigma=None,retained_v33_status=native['status'],retained_v34_status=previous['status'],native_recalculation=False,cases={'ground_coherent10':result})

if __name__=='__main__':print(json.dumps(summarize(Path(__file__).resolve().parents[1]/'data/audit/rya1354'),indent=2))
