"""Replay terminal NSIC0 preconditioning, preserving non-admission and load baseline."""
import argparse
import hashlib
import json
from pathlib import Path

from rya1354_native_parts import load_parts
from rya1354_node_failure_summary import attempts
from rya1354_scf_history import history
from rya1354_staged_failure import update_counts
from rya1354_lifetime_relax import scf_controls
from rya1354_wave_probe import decode

ROOT=Path(__file__).resolve().parents[1]


def summarize(base):
    files,manifest=load_parts(base/'staged_nsic0_failed_v11')
    root='runs/lifetime_relax_n4_staged_nsic0_v6/';stage=root+'fixed_2s_ground/'
    receipt=json.loads(files[root+'receipt.json'])
    if (receipt['status']!='FAILED_RETAINED' or receipt['optimization_nsic']!=0 or receipt['diagnostics']
            or receipt.get('last_completed_phase') or any('released_2s_' in name or name.endswith('state.lsj.lbl') for name in files)):
        raise ValueError('Failed NSIC0 campaign unexpectedly qualified/released')
    if files[stage+'rmcdhf.stdin'].decode()!=scf_controls('ground',4,1500,.5,False,
            radial_iteration_floor=24,both_2p_spectroscopic=True,freeze_2s=True,optimization_nsic=0):
        raise ValueError('Native NSIC0 input differs')
    for name in ('isodata','rcsf.inp','rwfn.inp'):
        if files[stage+name]!=files[stage+'fixed_load/'+name]:
            raise ValueError('All-fixed baseline input differs')
    if files[stage+'fixed_load/rmcdhf.stdin'].decode()!=scf_controls('ground',4,1,.5,True,
            radial_iteration_floor=24,both_2p_spectroscopic=True,all_fixed=True):
        raise ValueError('All-fixed baseline control differs')
    a=decode(files[stage+'fixed_load/rwfn.out']);b=decode(files[stage+'rwfn.out'])
    fixed_records=[dict(n=k[0],kappa=k[1],exact_native_baseline_match=a[k]==b[k],
        max_P_change=max(abs(x-y) for x,y in zip(a[k]['p'],b[k]['p']))) for k in ((1,-1),(2,-1))]
    text=files[stage+'rmcdhf.stdout'].decode();rows=history(text);complete=[r for r in rows if r['complete_printed_cycle']]
    counts=update_counts(text)
    if any(counts[r['iteration']]!=14 for r in complete):
        raise ValueError('Unexpected extra orbital updates in NSIC0 control')
    failed=attempts(text)
    if failed[-1]['orbital']!='2p' or failed[-1]['counted_nodes']!=1 or failed[-1]['required_nodes']!=0:
        raise ValueError('Unexpected final node refusal')
    return dict(schema='rya1354.nsic0-failure.v11',scientific_admission='HOLD',exposure='EXPOSED_NONBLIND',
        adoption_allowed=False,accepted_lifetime=None,uncertainty=None,covariance=None,calibrated_sigma=None,
        generator_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        parts=manifest['parts'],native_receipt_sha256=hashlib.sha256(files[root+'receipt.json']).hexdigest(),
        native_status=receipt['status'],native_error=receipt['error'],last_printed_cycle=rows[-1]['iteration'],
        last_complete_cycle=complete[-1],failed_attempts=failed,
        all_completed_cycles_have_14_regular_updates=True,fixed_records=fixed_records,
        native_loading_baseline_exact_for_fixed1s2s=all(r['exact_native_baseline_match'] for r in fixed_records),
        conclusion='NSIC0 still fails2p node selection before any stage qualifies; removing extra updates is insufficient. Fixed1s matches native loading baseline; fixed2s differs by roundoff. Upstream failure recovery calls global ORTHSC, which can renormalize fixed functions; the exact bit-change timing is not established. Exact frozen-record guard remains unchanged; no preconditioner/restart/release/lifetime accepted. Physical instability unresolved.')


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--base',type=Path,default=ROOT/'data/audit/rya1354')
    ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args();args.out.write_text(json.dumps(summarize(args.base),indent=2)+'\n')
