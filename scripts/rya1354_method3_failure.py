"""Replay the rejected n4 method3 attempt; never emit an accepted lifetime."""
import argparse
import hashlib
import json
from pathlib import Path

from rya1354_native_parts import load_parts
from rya1354_two_method_summary import scf_diagnostics

ROOT=Path(__file__).resolve().parents[1]


def summarize(bundle):
    files,manifest=load_parts(bundle)
    prefix='runs/lifetime_relax_n4_method3_v2/'
    receipt=json.loads(files[prefix+'receipt.json'])
    if receipt['status']!='FAILED_RETAINED' or receipt.get('two_s_method3') is not True or receipt['n']!=4:
        raise ValueError('Unexpected method3 failure receipt')
    stdout=files[prefix+'ground/rmcdhf.stdout'].decode()
    stderr=files[prefix+'ground/rmcdhf.stderr'].decode()
    if ('continuum state' not in stdout or 'dimensional limit (currently 439)' not in stdout
            or any(prefix+state+'/state.lsj.lbl' in files for state in ('ground','quintet','triplet'))):
        raise ValueError('Native failure or non-promotion evidence differs')
    diagnostics=scf_diagnostics(stdout,stderr)
    if diagnostics['all_reported_functions_below_ACC']:
        raise ValueError('Failed method3 attempt cannot be orbital-qualified')
    return dict(schema='rya1354.method3-failure.v7',scientific_admission='HOLD',
        exposure='EXPOSED_NONBLIND',adoption_allowed=False,accepted_lifetime=None,
        uncertainty=None,covariance=None,calibrated_sigma=None,
        generator_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        parts=manifest['parts'],native_receipt_sha256=hashlib.sha256(files[prefix+'receipt.json']).hexdigest(),
        native_status=receipt['status'],native_error=receipt['error'],diagnostics=diagnostics,
        disposition='Native IN refuses radial extent beyond439 and warns of possible continuum state. Optimization diverged; no ground/fresh restart/lifetime promoted. Method3 is not a successful remedy or a precision result.')


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--bundle',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args()
    args.out.write_text(json.dumps(summarize(args.bundle),indent=2)+'\n')
