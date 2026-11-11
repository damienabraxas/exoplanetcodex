"""Audit closed outer-layer ground stages; no final relaxed state or lifetime."""
import argparse
import hashlib
import json
from pathlib import Path

from rya1354_native_parts import load_parts
from rya1354_lifetime_summary import scf_state
from rya1354_lifetime_relax import scf_controls
from rya1354_wave_probe import decode

ROOT=Path(__file__).resolve().parents[1]


def summarize(base):
    files,manifest=load_parts(base/'from_n3_ground_partial_v11')
    checkpoint=json.loads(files['evidence/from_n3_ground_closed_v11.json'])
    cases=['outer_n4_ground','outer_n4_ground_verify']
    if checkpoint['cases']!=cases or checkpoint['accepted_lifetime'] is not None or checkpoint['adoption_allowed']:
        raise ValueError('Closed-stage checkpoint scope differs')
    for name,record in checkpoint['files'].items():
        if name not in files:
            if Path(name).name.startswith('mcp.') and name in manifest['omitted_angular_caches']:continue
            raise ValueError('Closed-stage native artifact missing')
        if hashlib.sha256(files[name]).hexdigest()!=record['sha256'] or len(files[name])!=record['bytes']:
            raise ValueError('Closed-stage snapshot hash differs')
    target,_=load_parts(base/'lifetime_aux_v4b')
    initial,_=load_parts(base/'from_n3_checkpoint_v11')
    root='runs/lifetime_relax_n4_from_relaxed_n3_v7/';states={}
    for case in cases:
        path=root+case;verify=case.endswith('_verify')
        diagnostics=scf_state(files,path,'ground')
        if any(diagnostics[k]!=v for k,v in checkpoint['diagnostics'][case].items()):
            raise ValueError('Native diagnostic differs from checkpoint')
        if files[path+'/rmcdhf.stdin'].decode()!=scf_controls('ground',4,1500,.5,verify,
                radial_iteration_floor=24,both_2p_spectroscopic=True,optimization_nsic=0,outer_layer_only=True):
            raise ValueError('Outer-layer native input differs')
        for name in ('isodata','rcsf.inp','rwfn.inp'):
            if files[path+'/'+name]!=files[path+'/fixed_load/'+name]:
                raise ValueError('Loading baseline input differs')
        if files[path+'/fixed_load/rmcdhf.stdin'].decode()!=scf_controls('ground',4,1,.5,True,
                radial_iteration_floor=24,both_2p_spectroscopic=True,all_fixed=True):
            raise ValueError('All-fixed baseline control differs')
        if files[path+'/state.c']!=target['runs/lifetime_core_n4_v1/CV_2s_DC_ground/state.c']:
            raise ValueError('Target full2s-SD configuration space differs')
        for adopted,native in (('state.c','rcsf.inp'),('state.w','rwfn.out'),('state.m','rmix.out')):
            if files[path+'/'+adopted]!=files[path+'/'+native]:
                raise ValueError('Promoted state differs from native output')
        if files[path+'/isodata']!=target['runs/lifetime_core_n4_v1/CV_2s_DC_ground/isodata']:
            raise ValueError('Nucleus identity differs')
        a=decode(files[path+'/fixed_load/rwfn.out']);b=decode(files[path+'/state.w'])
        if set(a)!=set(b) or len(a)!=16:
            raise ValueError('Native orbital inventory differs')
        fixed=[k for k in a if k[0]<=3]
        if len(fixed)!=9 or any(a[k]!=b[k] for k in fixed):
            raise ValueError('Frozen n<=3 record changed')
        source=files[root+'outer_n4_ground/state.w'] if verify else initial[root+'initial_ground/state.w']
        if files[path+'/rwfn.inp']!=source:
            raise ValueError('Initializer/fresh-restart wave identity differs')
        states[case]=dict(diagnostics=diagnostics,exact_frozen_n_le_3_records=9,
                          physical_labels_verified=True,source_wave_preserved=True)
    return dict(schema='rya1354.closed-ground-preconditioner-audit.v11',scientific_admission='HOLD',
        exposure='EXPOSED_NONBLIND',adoption_allowed=False,accepted_lifetime=None,
        uncertainty=None,covariance=None,calibrated_sigma=None,
        generator_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        captured_utc=checkpoint['captured_utc'],parts=manifest['parts'],states=states,
        scope='Only restricted outer4* ground optimization and fresh restart numerically qualified, with n<=3 fixed. Odd states/full release/campaign/lifetime remain pending. No final relaxed-state/accuracy/grade admission.')


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--base',type=Path,default=ROOT/'data/audit/rya1354')
    ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args();args.out.write_text(json.dumps(summarize(args.base),indent=2)+'\n')
