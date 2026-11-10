"""Audit terminal n4 MR evidence and retain the failed relaxation separately."""
import argparse
import hashlib
import json
from pathlib import Path

from rya1354_native_parts import load_parts
from rya1354_mr_summary import fixed_cases
from rya1354_lifetime_ladder import STATES

ROOT=Path(__file__).resolve().parents[1]


def summarize(mr_bundle,relax_bundle):
    files,mr=load_parts(mr_bundle)
    root='runs/lifetime_mr_n4_v1'
    receipt=json.loads(files[root+'/receipt.json'])
    dc='n\nn\nn\nn\nn\nn\n';breit='y\ny\n1.d-6\nn\nn\nn\nn\n'
    specs=[dict(name=name,core=core,multireference=mr,physics=physics)
           for name,core,mr,physics in [('DC',False,False,dc),('VV_MR_DC',False,True,dc),
           ('CV_2s_MR_DC',True,True,dc),('CV_2s_MR_Breit',True,True,breit)]]
    if (receipt['status']!='BALANCED_CORE_LIFETIMES_COMPLETED_REVIEW_REQUIRED'
            or receipt['n']!=4 or receipt['multireference'] is not True
            or receipt['case_specs']!=specs):
        raise ValueError('Unexpected completed MR model inventory')
    for spec in specs:
        for state,state_spec in STATES.items():
            roots=state_spec['roots'].removesuffix('5\n') if state=='ground' else state_spec['roots']
            expected=f'y\nstate\n{spec["physics"]}{roots}'
            if files[f'{root}/{spec["name"]}_{state}/rci.stdin'].decode()!=expected:
                raise ValueError('Native Hamiltonian/root input differs from declared MR model')
    cases=fixed_cases(files,receipt['case_specs'],root)
    failed,relax=load_parts(relax_bundle)
    prefix='runs/lifetime_relax_n4_v1/'
    failure=json.loads(failed[prefix+'receipt.json'])
    stdout=failed[prefix+'ground/rmcdhf.stdout'].decode()
    if (failure['status']!='FAILED_RETAINED' or failure['n']!=4
            or 'Convergence not obtained' not in stdout
            or any(prefix+state+'/state.lsj.lbl' in failed for state in ('ground','quintet','triplet'))):
        raise ValueError('Failed relaxation disposition does not match native evidence')
    return dict(schema='rya1354.terminal-lifetime-audit.v6',exposure='EXPOSED_NONBLIND',
        scientific_admission='HOLD',adoption_allowed=False,uncertainty=None,covariance=None,
        calibrated_sigma=None,chi_square=None,candidate_doi=None,
        generator_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        fixed_orbital_MR=dict(cases=cases,parts=mr['parts'],same_space_DC_reproduced=True,
            scope='Fixed qualified n4 VV orbitals; MR CSFs verified against native generator; all13channels. No MR orbital optimization or isotope precision; stationary-nucleus RCI warnings retained.'),
        failed_relaxation=dict(status=failure['status'],error=failure['error'],parts=relax['parts'],
            accepted_lifetime=None,native_stdout_sha256=hashlib.sha256(failed[prefix+'ground/rmcdhf.stdout']).hexdigest(),
            disposition='Ground2s integration methods1/2 failed node selection. No optimized state promoted; no subsequent state or lifetime accepted.'))


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--mr-bundle',type=Path,required=True)
    ap.add_argument('--relax-bundle',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args()
    args.out.write_text(json.dumps(summarize(args.mr_bundle,args.relax_bundle),indent=2)+'\n')
