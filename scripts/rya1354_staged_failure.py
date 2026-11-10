"""Audit rejected fixed2s preconditioning and native fixed-load normalization."""
import argparse
import hashlib
import json
from pathlib import Path
import re

from rya1354_native_parts import load_parts
from rya1354_node_failure_summary import attempts
from rya1354_scf_history import history
from rya1354_lifetime_relax import scf_controls
from rya1354_wave_probe import decode

ROOT=Path(__file__).resolve().parents[1]


def update_counts(text):
    blocks=re.split(r'^\s*Iteration number\s+(\d+)\s*$',text,flags=re.M)
    counts={}
    for i in range(1,len(blocks),2):
        rows=[]
        for line in blocks[i+1].splitlines():
            fields=line.split()
            if len(fields)>=11 and re.fullmatch(r'\d+[spdfgh]-?',fields[0]) and fields[2] in ('1','2','3','4'):
                rows.append(fields[0])
        counts[int(blocks[i])]=len(rows)
    return counts


def summarize(base):
    files,parts=load_parts(base/'staged_failed_v10');prefix='runs/lifetime_relax_n4_staged_v5/'
    receipt=json.loads(files[prefix+'receipt.json'])
    if receipt['status']!='FAILED_RETAINED' or receipt['diagnostics'] or receipt.get('last_completed_phase'):
        raise ValueError('Failed preconditioner must not qualify any stage')
    if files[prefix+'fixed_2s_ground/rmcdhf.stdin'].decode()!=scf_controls('ground',4,1500,.5,False,
            radial_iteration_floor=24,both_2p_spectroscopic=True,freeze_2s=True):
        raise ValueError('Staged native input differs')
    if any('released_2s_' in name or name.endswith('state.lsj.lbl') for name in files):
        raise ValueError('Failed preconditioner unexpectedly promoted/released')
    stdout=files[prefix+'fixed_2s_ground/rmcdhf.stdout'].decode()
    rows=history(stdout);complete=[r for r in rows if r['complete_printed_cycle']];failed=attempts(stdout)
    if failed[-1]['orbital']!='2p-' or failed[-1]['counted_nodes']!=2 or failed[-1]['required_nodes']!=0:
        raise ValueError('Native staged node failure differs')
    control,cp=load_parts(base/'fixed_load_control_v10');root='validation/fixed_load_control_v10/'
    cr=json.loads(control[root+'receipt.json'])
    if (cr['status']!='NATIVE_LOADING_NORMALIZATION_VERIFIED_NO_STATE_ADMISSION' or
            cr['failed_receipt_sha256']!=hashlib.sha256(files[prefix+'receipt.json']).hexdigest()):
        raise ValueError('Native fixed-load control source/status differs')
    for name in ('isodata','rcsf.inp','rwfn.inp'):
        if control[root+name]!=files[prefix+'fixed_2s_ground/'+name]:
            raise ValueError('Loading-control input identity differs')
    a=decode(control[root+'fixed_load/rwfn.out']);b=decode(files[prefix+'fixed_2s_ground/rwfn.out'])
    if any(a[k]!=b[k] for k in ((1,-1),(2,-1))):
        raise ValueError('Native loading baseline does not preserve actual fixed records')
    counts=update_counts(stdout)
    return dict(schema='rya1354.staged-failure.v10',scientific_admission='HOLD',exposure='EXPOSED_NONBLIND',
        adoption_allowed=False,accepted_lifetime=None,uncertainty=None,covariance=None,calibrated_sigma=None,
        generator_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        parts=parts['parts'],native_receipt_sha256=hashlib.sha256(files[prefix+'receipt.json']).hexdigest(),
        native_status=receipt['status'],native_error=receipt['error'],
        last_printed_cycle=rows[-1]['iteration'],last_complete_cycle=complete[-1],failed_attempts=failed,
        first_three_cycle_update_counts={str(i):counts[i] for i in (1,2,3)},
        fixed_loading_control=dict(parts=cp['parts'],fixed_records=cr['fixed_records'],
            exact_after_native_loading=True,raw_preload_identity_is_not_the_correct_baseline=True),
        conclusion='Fixed2s preconditioning diverged and failed2p-; no stage/lifetime qualified. NSIC20 adds up to20updates inside everycycle; causal impact remains a separate control. Fixed1s/2s exactly match same-input native loading baseline; raw pre-load serialization differs by native normalization roundoff, not measured orbital evolution.')


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--base',type=Path,default=ROOT/'data/audit/rya1354')
    ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args();args.out.write_text(json.dumps(summarize(args.base),indent=2)+'\n')
