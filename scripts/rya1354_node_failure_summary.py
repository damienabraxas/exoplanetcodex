"""Replay the terminal n4 node failures without admitting a relaxed lifetime."""
import argparse
import hashlib
import json
from pathlib import Path
import re

from rya1354_native_parts import load_parts
from rya1354_scf_history import history
from rya1354_lifetime_relax import scf_controls

ROOT=Path(__file__).resolve().parents[1]


def attempts(text):
    pattern=(r'Method (\d) unable to solve for\s+(\d+)([spdfgh]-?)\s+orbital'
             r'.*?Iteration number:\s+(\d+), limit:\s+(\d+)'
             r'.*?Lower bound on energy:\s*([\d.D+\-]+), upper bound:\s*([\d.D+\-]+)'
             r'.*?Number of nodes counted:\s*(\d+), Correct number:\s*(\d+)')
    result=[]
    for m in re.finditer(pattern,text,re.S):
        method,n,shell,used,limit,lower,upper,nodes,target=m.groups()
        width=float(upper.replace('D','E'))-float(lower.replace('D','E'))
        result.append(dict(method=int(method),orbital=n+shell,iterations=int(used),limit=int(limit),
            counted_nodes=int(nodes),required_nodes=int(target),native_energy_bracket_width=width,
            native_width_floor=1/int(n)**3,iteration_limit_reached=int(used)>=int(limit),
            width_guard_triggered=width<1/int(n)**3))
    if not result:
        raise ValueError('Native failed radial-solve inventory missing')
    return result


def summarize(base):
    prior,_=load_parts(base/'lifetime_relax_failed_v6')
    old_history=history(prior['runs/lifetime_relax_n4_v1/ground/rmcdhf.stdout'].decode())
    seed,_=load_parts(base/'lifetime_aux_v4b')
    cases={}
    for name,bundle,both in [('nodebudget','nodebudget_failed_v9',False),('pcomponents','pcomponents_failed_v9',True)]:
        files,manifest=load_parts(base/bundle)
        run='lifetime_relax_n4_nodebudget_v3' if not both else 'lifetime_relax_n4_pcomponents_v4'
        prefix='runs/'+run+'/'
        receipt=json.loads(files[prefix+'receipt.json'])
        if (receipt['status']!='FAILED_RETAINED' or receipt['n']!=4 or receipt['radial_iteration_floor']!=24
                or receipt.get('both_2p_spectroscopic',False)!=both or receipt['two_s_method3']):
            raise ValueError('Unexpected node-failure model/receipt')
        if files[prefix+'ground/rmcdhf.stdin'].decode()!=scf_controls('ground',4,1500,.5,False,
                radial_iteration_floor=24,both_2p_spectroscopic=both):
            raise ValueError('Native node-control stdin differs')
        for field in ('state.c','state.w','isodata'):
            if files[prefix+'ground/'+field]!=seed['runs/lifetime_core_n4_v1/CV_2s_DC_ground/'+field]:
                raise ValueError('Qualified initializer bytes differ')
        if any(prefix+state+'/state.lsj.lbl' in files for state in ('ground','quintet','triplet')):
            raise ValueError('Failed ground unexpectedly promoted')
        stdout=files[prefix+'ground/rmcdhf.stdout'].decode()
        if 'Convergence not obtained' not in stdout:
            raise ValueError('Native node-failure termination absent')
        failed=attempts(stdout);rows=history(stdout)
        complete=[r for r in rows if r['complete_printed_cycle']]
        if (failed[-1]['orbital']!='2s' or failed[-1]['counted_nodes']!=3
                or failed[-1]['required_nodes']!=1 or failed[-1]['iteration_limit_reached']
                or not failed[-1]['width_guard_triggered']):
            raise ValueError('Final refusal does not demonstrate pre-cap node/bracket failure')
        cases[name]=dict(native_status=receipt['status'],native_error=receipt['error'],
            native_receipt_sha256=hashlib.sha256(files[prefix+'receipt.json']).hexdigest(),parts=manifest['parts'],
            last_printed_outer_iteration=rows[-1]['iteration'],last_complete_cycle=complete[-1],
            failed_attempts=failed,accepted_lifetime=None,qualified_initializer_preserved=True)
        if not both:
            if rows[:-1]!=old_history[:-1]:
                raise ValueError('Claimed unchanged pre-failure trajectory differs')
            cases[name]['first_193_cycles_exactly_reproduce_default_budget']=len(complete)==193
    return dict(schema='rya1354.node-failures.v9',exposure='EXPOSED_NONBLIND',scientific_admission='HOLD',
        adoption_allowed=False,uncertainty=None,covariance=None,calibrated_sigma=None,
        generator_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),cases=cases,
        conclusion='Both final2s refusals occur before24iterations with bracket width below native1/n^3 and3nodes instead of1. Increasing the cap or constraining both2pcomponents did not yield qualified states. Physical cause remains unresolved; no accepted relaxed lifetime.')


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--base',type=Path,default=ROOT/'data/audit/rya1354')
    ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args();args.out.write_text(json.dumps(summarize(args.base),indent=2)+'\n')
