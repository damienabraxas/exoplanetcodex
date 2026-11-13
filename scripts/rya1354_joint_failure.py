"""Retain timed-out joint4*/5* ground optimization; never accept its wave."""
import argparse
import hashlib
import json
from pathlib import Path

from rya1354_native_parts import load_parts
from rya1354_extension_audit import summarize as n5_summary
from rya1354_joint_outer import controls
from rya1354_frozen_extend import controls as n5_controls
from rya1354_scf_history import history
from rya1354_wave_probe import decode

ROOT=Path(__file__).resolve().parents[1]
RUN='runs/lifetime_joint45_n5_v1/'


def summarize(base):
    files,manifest=load_parts(base/'joint_failed_v16')
    prior=n5_summary(base,base/'n5_complete_v15');prior['parts']={}
    source,_=load_parts(base/'n5_complete_v15')
    receipt=json.loads(files[RUN+'receipt.json'])
    if (receipt['status']!='FAILED_RETAINED' or receipt['n']!=5 or receipt['diagnostics']
            or not receipt['error'].startswith('TimeoutExpired:') or '43200 seconds' not in receipt['error']):
        raise ValueError('Unexpected joint failure disposition')
    if receipt['source_audit']!=prior or receipt['source_receipt_sha256']!=hashlib.sha256(source['runs/lifetime_frozen_prior_n5_v1/receipt.json']).hexdigest():
        raise ValueError('Joint initializer qualification differs')
    path=RUN+'joint45_ground/'
    for name in ('state.c','state.w','isodata'):
        if files[path+name]!=source['runs/lifetime_frozen_prior_n5_v1/outer_n5_ground_verify/'+name]:
            raise ValueError('Joint failed ground source differs')
    for name in ('rcsf.inp','rwfn.inp','isodata'):
        if files[path+name]!=files[path+'fixed_load/'+name]:
            raise ValueError('Joint loading inputs differ')
    if (files[path+'rmcdhf.stdin'].decode()!=controls('ground')
            or files[path+'fixed_load/rmcdhf.stdin'].decode()!=n5_controls('ground',True,True)):
        raise ValueError('Joint declared radial controls differ')
    if files[path+'state.c']!=files[path+'rcsf.inp'] or files[path+'state.w']!=files[path+'rwfn.inp']:
        raise ValueError('Joint initializer aliases differ')
    rows=history(files[path+'rmcdhf.stdout'].decode())
    complete=[r for r in rows if r['complete_printed_cycle']]
    if not complete or any(r['max_reported_residual']<=.0625**6 for r in complete):
        raise ValueError('Unexpected joint residual trajectory; inspect before disposition')
    a,b=decode(files[path+'fixed_load/rwfn.out']),decode(files[path+'rwfn.out'])
    fixed=[k for k in a if k[0]<=3]
    if set(a)!=set(b) or len(a)!=25 or len(fixed)!=9:
        raise ValueError('Retained joint wave inventory differs')
    return dict(schema='rya1354.joint-failure-audit.v16',status=receipt['status'],error=receipt['error'],
        native_receipt_sha256=hashlib.sha256(files[RUN+'receipt.json']).hexdigest(),
        completed_cycle_count=len(complete),last_complete_cycle=complete[-1],last_printed_cycle=rows[-1],
        best_complete_cycle=min(complete,key=lambda x:x['max_reported_residual']),
        original_ACC=.0625**6,exact_fixed_n_le_3_records_in_retained_wave=sum(a[k]==b[k] for k in fixed),
        retained_wave_timing='Terminal ORBOUT file; exact cycle-to-file correspondence not asserted',
        scientific_admission='HOLD',exposure='EXPOSED_NONBLIND',adoption_allowed=False,
        accepted_lifetime=None,uncertainty=None,covariance=None,calibrated_sigma=None,
        qualified_states=0,parts=manifest['parts'],
        generator_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        disposition='12-hour ground timeout with unresolved oscillatory orbital residuals. No state/rate qualification, no failed-wave initializer and no root-cause or accuracy claim.')


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--base',type=Path,default=ROOT/'data/audit/rya1354')
    ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args();args.out.write_text(json.dumps(summarize(args.base),indent=2)+'\n')
