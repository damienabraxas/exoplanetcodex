"""Fixed2s preconditioning, then full2s release; no lifetime from intermediate states."""
import argparse
from datetime import datetime,timezone
import json
from pathlib import Path
import platform
import shutil
import sys

from rya1354_lifetime_relax import validate_core_layer,scf_controls
from rya1354_lifetime_ladder import STATES,identity
from rya1354_solver_smoke import run,sha
from rya1354_two_method_summary import scf_diagnostics
from rya1354_triplet_correlation import transitions,require_triplet_quintet_order
from rya1354_wave_probe import decode


def require_frozen(source,target,freeze_2s):
    """Fixed functions must retain their native records, not just their labels."""
    a,b=decode(source.read_bytes()),decode(target.read_bytes())
    if set(a)!=set(b):
        raise ValueError('Orbital inventory changed within matched staged model')
    for key in ((1,-1),(2,-1)) if freeze_2s else ((1,-1),):
        if a[key]!=b[key]:
            raise ValueError('Declared fixed radial function changed: '+str(key))


def optimize(args,phase,state,seed,receipt):
    freeze=phase=='fixed_2s';selected=None
    for verify in (False,True):
        out=args.out/(phase+'_'+state+('_verify' if verify else ''));out.mkdir()
        source=selected if verify else seed
        for name in ('isodata','state.c','state.w'):
            shutil.copyfile(source/name,out/name)
        shutil.copyfile(out/'state.c',out/'rcsf.inp');shutil.copyfile(out/'state.w',out/'rwfn.inp')
        receipt['steps'].append(run(args.grasp_bin/'rangular','y\n',out,1800))
        controls=scf_controls(state,4,args.cycles,args.damping,verify,radial_iteration_floor=24,
                              both_2p_spectroscopic=True,freeze_2s=freeze)
        receipt['steps'].append(run(args.strict_bin/'rmcdhf',controls,out,43200))
        if any(step['semantic_errors'] for step in receipt['steps']):
            raise ValueError('Native refusal; no staged-state promotion')
        diagnostics=scf_diagnostics((out/'rmcdhf.stdout').read_text(),(out/'rmcdhf.stderr').read_text())
        if diagnostics['iteration_limit_hit'] or not diagnostics['all_reported_functions_below_ACC']:
            raise ValueError('Staged orbital criterion unmet')
        require_frozen(out/'rwfn.inp',out/'rwfn.out',freeze)
        if (out/'rcsf.inp').read_bytes()!=(seed/'state.c').read_bytes() or (out/'isodata').read_bytes()!=(seed/'isodata').read_bytes():
            raise ValueError('Staged CSF/nucleus identity changed')
        for native,extension in (('rwfn.out','w'),('rmix.out','m')):
            shutil.copyfile(out/native,out/('state.'+extension))
        receipt['steps'].append(run(args.grasp_bin/'jj2lsj','state\nn\nn\ny\n',out,600))
        identity((out/'state.lsj.lbl').read_text(),state)
        receipt['diagnostics'][out.name]=diagnostics
        receipt['last_qualified_stage']=out.name
        (args.out/'progress.json').write_text(json.dumps(receipt,indent=2)+'\n')
        selected=out
    return selected


def optimize_stages(args,receipt):
    previous={state:args.core/('CV_2s_DC_'+state) for state in STATES}
    for phase in receipt['phases']:
        selected={state:optimize(args,phase,state,previous[state],receipt) for state in STATES}
        require_triplet_quintet_order((selected['triplet']/'state.lsj.lbl').read_text(),(selected['quintet']/'state.lsj.lbl').read_text())
        receipt['last_completed_phase']=phase;previous=selected
        (args.out/'progress.json').write_text(json.dumps(receipt,indent=2)+'\n')
    return selected


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--grasp-bin',type=Path,required=True)
    ap.add_argument('--strict-bin',type=Path,required=True)
    ap.add_argument('--core',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--cycles',type=int,default=1500)
    ap.add_argument('--damping',type=float,default=.5)
    args=ap.parse_args()
    if platform.node().lower().split('.')[0]!='sirius':
        raise RuntimeError('Atomic solvers belong on Sirius')
    if args.cycles<1 or not 0<args.damping<1:
        raise ValueError('Invalid staged numerical controls')
    for attr in ('grasp_bin','strict_bin','core','out'):
        setattr(args,attr,getattr(args,attr).resolve())
    initializer=validate_core_layer(args.core,4)
    args.out.mkdir(parents=True,exist_ok=False)
    drivers=args.out/'drivers';drivers.mkdir();shutil.copyfile(__file__,drivers/'driver.py')
    for name in ('rya1354_lifetime_relax','rya1354_lifetime_ladder','rya1354_solver_smoke',
                 'rya1354_two_method_summary','rya1354_triplet_correlation','rya1354_wave_probe'):
        shutil.copyfile(sys.modules[name].__file__,drivers/(name+'.py'))
    receipt=dict(schema='rya1354.staged-core-relaxation.v9',status='RUNNING',n=4,
        scientific_admission='HOLD',adoption_allowed=False,started_utc=datetime.now(timezone.utc).isoformat(),
        core_receipt_sha256=initializer,steps=[],diagnostics={},cycles=args.cycles,damping=args.damping,
        phases=['fixed_2s','released_2s'],accepted_lifetime=None,
        model='Same n4 full2s-SD/DC CSFs. All3states first vary2p*/3*/4* with1s/2s fixed, strict/fresh verification; then vary2*/3*/4* with1s fixed, strict/fresh verification. Explicit2p-/2p spectroscopic, NSOLV24, original439grid/ACCY. No observed/failed-wave inputs; no lifetime from preconditioner.')
    try:
        selected=optimize_stages(args,receipt)
        for name,lower,operators,same in [('to_ground',selected['ground'],'E1,M2,E3',False),
                                         ('to_5S',selected['quintet'],'M1,E2,M3',True)]:
            receipt['steps']+=transitions(args.grasp_bin,lower,selected['triplet'],args.out/name,operators,same)
        receipt['status']='CORE_RELAXED_MODEL_LIFETIME_COMPLETED_REVIEW_REQUIRED'
        # Terminal native completion still requires independent offline lifetime audit.
    except Exception as exc:
        receipt.update(status='FAILED_RETAINED',error=type(exc).__name__+': '+str(exc));raise
    finally:
        receipt['files']={str(p.relative_to(args.out)):dict(sha256=sha(p),bytes=p.stat().st_size)
                          for p in sorted(args.out.rglob('*')) if p.is_file()}
        (args.out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')


if __name__=='__main__':
    main()
