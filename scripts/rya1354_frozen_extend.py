"""Declared layerwise n4-to-n5 frozen-prior-layer sensitivity, Sirius only."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import platform
import shutil

from rya1354_frozen_layer import validate, digest, PARENT, TARGET, STATES, STATUS
from rya1354_lifetime_relax import validate_core_layer, scf_controls
from rya1354_lifetime_summary import scf_state
from rya1354_lifetime_ladder import identity
from rya1354_lifetime_staged import require_frozen
from rya1354_radial_controls import initialization_input
from rya1354_solver_smoke import run, sha
from rya1354_triplet_correlation import transitions, require_triplet_quintet_order
from rya1354_wave_probe import decode
from rya1354_mr_inventory import inventory


def source_files(root,prefix):
    return {prefix+str(p.relative_to(root)):p.read_bytes() for p in root.rglob('*')
            if p.is_file() and not p.name.startswith('mcp.')}


def qualify_source(source,parent,core):
    """Recheck original six stages and immutable frozen-model selected states."""
    original=source_files(parent,PARENT)
    target={TARGET+str(p.relative_to(core)):p.read_bytes() for p in
            [core/'receipt.json',*[core/('CV_2s_DC_'+s)/n for s in STATES for n in ('state.c','isodata')]]}
    audit=validate(original,target)
    receipt=json.loads((source/'receipt.json').read_text())
    if receipt['status']!=STATUS or receipt['source_audit']!=audit:
        raise ValueError('Unqualified frozen-layer source')
    for name,record in receipt['files'].items():
        p=source/name
        if not p.is_file() or sha(p)!=record['sha256'] or p.stat().st_size!=record['bytes']:
            raise ValueError('Frozen-layer source artifact differs: '+name)
    for state in STATES:
        for name in ('state.c','state.w','state.m','state.lsj.lbl','isodata'):
            if (source/state/name).read_bytes()!=original[PARENT+'outer_n4_'+state+'_verify/'+name]:
                raise ValueError('Frozen-layer selected state differs')
    return dict(receipt_sha256=sha(source/'receipt.json'),source_audit=audit)


def controls(state,verify=False,all_fixed=False):
    return scf_controls(state,5,1 if all_fixed else 1500,.5,verify,
        radial_iteration_floor=24,both_2p_spectroscopic=True,optimization_nsic=0,
        outer_layer_only=not all_fixed,all_fixed=all_fixed)


def native_text(root,exe):
    return (root/(exe+'.stdout')).read_text()+(root/(exe+'.stderr')).read_text()


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    for name in ('source','parent','source-core','core','grasp-bin','strict-bin','out'):
        ap.add_argument('--'+name,type=Path,required=True)
    args=ap.parse_args()
    if platform.node().lower().split('.')[0]!='sirius':
        raise RuntimeError('Atomic solvers belong on Sirius')
    for name,value in vars(args).items():
        setattr(args,name,value.resolve())
    source=qualify_source(args.source,args.parent,args.source_core)
    target=validate_core_layer(args.core,5)
    args.out.mkdir(parents=True,exist_ok=False)
    drivers=args.out/'drivers';drivers.mkdir()
    for p in Path(__file__).parent.glob('rya1354_*.py'):
        shutil.copyfile(p,drivers/p.name)
    receipt=dict(schema='rya1354.frozen-extension.v13',status='RUNNING',n=5,
        source=source,target_core_receipt_sha256=target,steps=[],diagnostics={},
        started_utc=datetime.now(timezone.utc).isoformat(),scientific_admission='HOLD',
        exposure='EXPOSED_NONBLIND',adoption_allowed=False,accepted_lifetime=None,
        uncertainty=None,covariance=None,calibrated_sigma=None,
        model='Layerwise full2s-SD/DC n5: qualified frozen-layer n4 shapes fixed, only new n5 optimized. Original grid/ACCY and fresh restart. Historical n3 component limitation retained. No full core relaxation, observed energy substitution or accuracy claim.')
    def progress():
        (args.out/'progress.json').write_text(json.dumps(receipt,indent=2)+'\n')
    try:
        selected={}
        for state in STATES:
            seed=args.out/('initial_'+state);seed.mkdir()
            for name in ('state.c','isodata'):
                shutil.copyfile(args.core/('CV_2s_DC_'+state)/name,seed/name)
            inventory((seed/'state.c').read_text(),'even' if state=='ground' else 'odd',True)
            if (seed/'isodata').read_bytes()!=(args.source/state/'isodata').read_bytes():
                raise ValueError('Source/target nucleus differs')
            shutil.copyfile(seed/'state.c',seed/'rcsf.inp')
            shutil.copyfile(args.source/state/'state.w',seed/'previous.w')
            receipt['steps'].append(run(args.grasp_bin/'rwfnestimate',initialization_input(3,439,.03125),seed,1800))
            if ('All required subshell radial wavefunctions' not in native_text(seed,'rwfnestimate')
                    or receipt['steps'][-1]['semantic_errors']):
                raise ValueError('Native n5 initialization refused')
            wave=decode((seed/'rwfn.inp').read_bytes())
            if set(wave)!=set(decode((args.core/('CV_2s_DC_'+state)/'state.w').read_bytes())) or len(wave)!=25:
                raise ValueError('Initialized n5 inventory differs')
            shutil.copyfile(seed/'rwfn.inp',seed/'state.w')
            for verify in (False,True):
                out=args.out/('outer_n5_'+state+('_verify' if verify else ''));out.mkdir()
                for name in ('state.c','state.w','isodata'):
                    shutil.copyfile(seed/name,out/name)
                shutil.copyfile(out/'state.c',out/'rcsf.inp');shutil.copyfile(out/'state.w',out/'rwfn.inp')
                receipt['steps'].append(run(args.grasp_bin/'rangular','y\n',out,1800))
                baseline=out/'fixed_load';baseline.mkdir()
                for name in ('rcsf.inp','rwfn.inp','isodata'):
                    shutil.copyfile(out/name,baseline/name)
                receipt['steps'].append(run(args.grasp_bin/'rangular','y\n',baseline,1800))
                receipt['steps'].append(run(args.strict_bin/'rmcdhf',controls(state,True,True),baseline,43200))
                if 'All subshell radial wavefunctions are fixed' not in native_text(baseline,'rmcdhf'):
                    raise ValueError('Native loading baseline refused')
                receipt['steps'].append(run(args.strict_bin/'rmcdhf',controls(state,verify),out,43200))
                if any(s['semantic_errors'] for s in receipt['steps']):
                    raise ValueError('Native semantic refusal')
                require_frozen(baseline/'rwfn.out',out/'rwfn.out',False,fixed_max_n=4)
                a=decode((baseline/'rwfn.out').read_bytes())
                if len(a)!=25 or sum(k[0]<=4 for k in a)!=16:
                    raise ValueError('Native loading inventory differs')
                for name in ('state.c','isodata'):
                    if (out/name).read_bytes()!=(args.core/('CV_2s_DC_'+state)/name).read_bytes():
                        raise ValueError('Target CSF/nucleus differs')
                for native,extension in (('rwfn.out','w'),('rmix.out','m')):
                    shutil.copyfile(out/native,out/('state.'+extension))
                receipt['steps'].append(run(args.grasp_bin/'jj2lsj','state\nn\nn\ny\n',out,600))
                identity((out/'state.lsj.lbl').read_text(),state)
                diagnostic=scf_state(source_files(out,'stage/'),'stage',state)
                receipt['diagnostics'][out.name]=diagnostic
                receipt['last_qualified_stage']=out.name;progress();seed=out
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
