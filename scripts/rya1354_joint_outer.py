"""Matched n5 sensitivity: jointly optimize4*/5*, retain core-adapted n<=3."""
import argparse
from datetime import datetime,timezone
import hashlib
import json
from pathlib import Path
import platform
import shutil

from rya1354_frozen_extend import source_files,qualify_source,controls as n5_controls,native_text
from rya1354_extension_audit import audit_files,RUN,CORE,SOURCE
from rya1354_lifetime_relax import validate_core_layer
from rya1354_frozen_layer import STATES,STATUS
from rya1354_lifetime_staged import require_frozen
from rya1354_lifetime_summary import scf_state,decay
from rya1354_solver_smoke import run,sha
from rya1354_triplet_correlation import transitions,require_triplet_quintet_order
from rya1354_wave_probe import decode


def controls(state,verify=False):
    text=n5_controls(state,verify)
    if text.count('\n5*\n')!=1:
        raise ValueError('Unexpected native outer-shell control layout')
    return text.replace('\n5*\n','\n4*,5*\n')


def qualified_source(args):
    prior=qualify_source(args.n4_source,args.parent,args.n4_core)
    validate_core_layer(args.core,5)
    source=source_files(args.n4_source,SOURCE)
    target={CORE+str(p.relative_to(args.core)):p.read_bytes() for p in
            [args.core/'receipt.json',*[args.core/('CV_2s_DC_'+s)/n for s in STATES for n in ('state.c','state.w','isodata')]]}
    files=source_files(args.source,RUN)
    receipt=json.loads(files[RUN+'receipt.json'])
    for name,record in receipt['files'].items():
        if Path(name).name.startswith('mcp.'):
            continue
        value=files[RUN+name]
        if len(value)!=record['bytes'] or hashlib.sha256(value).hexdigest()!=record['sha256']:
            raise ValueError('N5 source artifact differs: '+name)
    return audit_files(files,target,source,dict(source_audit=prior['source_audit'],
                       decay=decay(source,SOURCE)),{})


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    for name in ('source','n4-source','parent','n4-core','core','grasp-bin','strict-bin','out'):
        ap.add_argument('--'+name,type=Path,required=True)
    args=ap.parse_args()
    if platform.node().lower().split('.')[0]!='sirius':
        raise RuntimeError('Atomic solvers belong on Sirius')
    for name,value in vars(args).items():
        setattr(args,name,value.resolve())
    audit=qualified_source(args)
    args.out.mkdir(parents=True,exist_ok=False)
    drivers=args.out/'drivers';drivers.mkdir()
    for p in Path(__file__).parent.glob('rya1354_*.py'):
        shutil.copyfile(p,drivers/p.name)
    receipt=dict(schema='rya1354.joint-outer-n5.v15',status='RUNNING',n=5,
        source_receipt_sha256=sha(args.source/'receipt.json'),source_audit=audit,
        started_utc=datetime.now(timezone.utc).isoformat(),steps=[],diagnostics={},
        scientific_admission='HOLD',exposure='EXPOSED_NONBLIND',adoption_allowed=False,
        accepted_lifetime=None,uncertainty=None,covariance=None,calibrated_sigma=None,
        model='Same n5 full2s-SD/DC CSFs/nucleus and qualified sequential-layer n5 initializer. Jointly vary4*,5*; freeze9 n<=3 records against native same-input loading control. Original grid/ACCY, NSIC0, fresh restart. Historical n3 component limitation and rejected full-release failure retained. No observed energy substitution or accuracy claim.')
    try:
        selected={}
        for state in STATES:
            seed=args.source/('outer_n5_'+state+'_verify')
            for verify in (False,True):
                out=args.out/('joint45_'+state+('_verify' if verify else ''));out.mkdir()
                for name in ('state.c','state.w','isodata'):
                    shutil.copyfile(seed/name,out/name)
                shutil.copyfile(out/'state.c',out/'rcsf.inp');shutil.copyfile(out/'state.w',out/'rwfn.inp')
                receipt['steps'].append(run(args.grasp_bin/'rangular','y\n',out,1800))
                baseline=out/'fixed_load';baseline.mkdir()
                for name in ('rcsf.inp','rwfn.inp','isodata'):
                    shutil.copyfile(out/name,baseline/name)
                receipt['steps'].append(run(args.grasp_bin/'rangular','y\n',baseline,1800))
                receipt['steps'].append(run(args.strict_bin/'rmcdhf',n5_controls(state,True,True),baseline,43200))
                if 'All subshell radial wavefunctions are fixed' not in native_text(baseline,'rmcdhf'):
                    raise ValueError('Native all-fixed loading control refused')
                loaded=decode((baseline/'rwfn.out').read_bytes())
                if (set(loaded)!=set(decode((out/'rwfn.inp').read_bytes())) or len(loaded)!=25
                        or sum(k[0]<=3 for k in loaded)!=9):
                    raise ValueError('Joint45 native loading inventory differs')
                receipt['steps'].append(run(args.strict_bin/'rmcdhf',controls(state,verify),out,43200))
                if any(s['semantic_errors'] for s in receipt['steps']):
                    raise ValueError('Native semantic refusal')
                require_frozen(baseline/'rwfn.out',out/'rwfn.out',False,fixed_max_n=3)
                for name in ('state.c','isodata'):
                    if (out/name).read_bytes()!=(args.source/('outer_n5_'+state+'_verify')/name).read_bytes():
                        raise ValueError('Joint45 CSF/nucleus differs')
                for native,extension in (('rwfn.out','w'),('rmix.out','m')):
                    shutil.copyfile(out/native,out/('state.'+extension))
                receipt['steps'].append(run(args.grasp_bin/'jj2lsj','state\nn\nn\ny\n',out,600))
                diagnostic=scf_state(source_files(out,'stage/'),'stage',state)
                receipt['diagnostics'][out.name]=diagnostic;receipt['last_qualified_stage']=out.name
                (args.out/'progress.json').write_text(json.dumps(receipt,indent=2)+'\n')
                seed=out
            selected[state]=seed
        require_triplet_quintet_order((selected['triplet']/'state.lsj.lbl').read_text(),(selected['quintet']/'state.lsj.lbl').read_text())
        for name,lower,operators,same in [('to_ground','ground','E1,M2,E3',False),('to_5S','quintet','M1,E2,M3',True)]:
            receipt['steps']+=transitions(args.grasp_bin,selected[lower],selected['triplet'],args.out/name,operators,same)
        receipt['status']=STATUS
    except Exception as exc:
        receipt.update(status='FAILED_RETAINED',error=type(exc).__name__+': '+str(exc))
        raise
    finally:
        receipt['files']={str(p.relative_to(args.out)):dict(sha256=sha(p),bytes=p.stat().st_size)
                          for p in sorted(args.out.rglob('*')) if p.is_file()}
        (args.out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')


if __name__=='__main__':
    main()
