"""Reoptimize matched n3/n4/n5 core-expanded states; preserve fixed-orbital control."""
import argparse
from datetime import datetime,timezone
import json
import math
from pathlib import Path
import platform
import shutil

from rya1354_solver_smoke import run,sha
from rya1354_lifetime_ladder import STATES,identity
from rya1354_two_method_summary import scf_diagnostics
from rya1354_triplet_correlation import transitions,require_triplet_quintet_order


def variational_orbitals(n):
    """Same radial scope for all states; keep only 1s fixed."""
    if n not in (3, 4, 5):
        raise ValueError('Supported core layers are n3/n4/n5')
    return ','.join(f'{layer}*' for layer in range(2, n + 1))


def validate_core_layer(core, n):
    """Refuse a mislabeled/incomplete fixed-orbital initializer campaign."""
    receipt = json.loads((core / 'receipt.json').read_text())
    if (receipt.get('n') != n or receipt.get('status') !=
            'BALANCED_CORE_LIFETIMES_COMPLETED_REVIEW_REQUIRED'):
        raise ValueError('Core initializer layer/status differs from requested model')
    for state in STATES:
        seed = core / ('CV_2s_DC_' + state)
        identity((seed / 'state.lsj.lbl').read_text(), state)
        for name in ('state.c', 'state.w', 'isodata'):
            path = seed / name
            record = receipt['files'][str(path.relative_to(core))]
            if sha(path) != record['sha256'] or path.stat().st_size != record['bytes']:
                raise ValueError('Core initializer differs from native receipt: ' + str(path))
    return sha(core / 'receipt.json')


def scf_controls(state,n,cycles,damping,verify,*,two_s_method3=False,radial_iteration_floor=None,both_2p_spectroscopic=False,freeze_2s=False,optimization_nsic=20,all_fixed=False):
    """Try method3 for2s optimization; fresh restart restores node-checked defaults."""
    spectroscopic=('2s,2p-,2p' if both_2p_spectroscopic else '2s,2p')+(',3s' if state!='ground' else '')
    grid=f'n\nn\ny\nn\ny\n{math.exp(-65/16)/8:.17g}\n0.03125\n0\n439\ny\n{.0625**6:.17g}\n'
    if optimization_nsic<0:
        raise ValueError('Additional orbital updates must be nonnegative')
    nsic='n\n' if verify else f'y\n{optimization_nsic}\n'
    methods='y\n\n\n2s\n\n' if two_s_method3 and not verify else 'n\n'
    if radial_iteration_floor is not None and radial_iteration_floor<1:
        raise ValueError('Native radial iteration floor must be positive')
    # Pinned SOLVE uses MAX(NSOLV,3*n), despite GETSCD calling NSOLV a factor.
    radial=f'y\n{radial_iteration_floor}\n' if radial_iteration_floor is not None else 'n\n'
    varied=variational_orbitals(n)
    if freeze_2s:
        varied=varied.replace('2*','2p*',1)
    if all_fixed:
        varied=''
    return (f'{grid}{STATES[state]["roots"]}{varied}\n{spectroscopic}\n{cycles}\n'
            f'y\nn\n{methods}n\ny\nn\n{-damping}\nn\n{nsic}{radial}n\n1\n')


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--grasp-bin',type=Path,required=True)
    ap.add_argument('--strict-bin',type=Path,required=True)
    ap.add_argument('--core',type=Path,required=True)
    ap.add_argument('--n',type=int,choices=(3,4,5),default=3)
    ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--damping',type=float,default=.9)
    ap.add_argument('--resume-ground',type=Path,help='Explicit same-CSF failed ground stage as an unqualified initializer only')
    ap.add_argument('--resume-quintet',type=Path,help='Explicit same-CSF quintet continuation initializer')
    ap.add_argument('--resume-triplet',type=Path,help='Explicit same-CSF triplet continuation initializer')
    ap.add_argument('--cycles',type=int,default=500)
    ap.add_argument('--two-s-method3',action='store_true',help='Method3 for2s optimization only; mandatory fresh default-method node-checked restart before promotion')
    ap.add_argument('--radial-iteration-floor',type=int,help='Explicit NSOLV: pinned SOLVE permits MAX(NSOLV,3*n) inner iterations; node and bracket criteria unchanged')
    ap.add_argument('--both-2p-spectroscopic',action='store_true',help='Explicitly constrain both physical2p-/2p components; preserve old selection as a separate control')
    args=ap.parse_args()
    if not 0<args.damping<1 or args.cycles<1:
        raise ValueError('Damping must lie strictly between zero and one')
    if args.radial_iteration_floor is not None and args.radial_iteration_floor<1:
        raise ValueError('Native radial iteration floor must be positive')
    if platform.node().lower().split('.')[0]!='sirius':
        raise RuntimeError('Atomic solvers belong on Sirius')
    args.grasp_bin=args.grasp_bin.resolve();args.strict_bin=args.strict_bin.resolve()
    core_receipt_sha256=validate_core_layer(args.core,args.n)
    varied=variational_orbitals(args.n)
    args.out.mkdir(parents=True,exist_ok=False)
    shutil.copyfile(__file__,args.out/'driver.py')
    receipt=dict(schema='rya1354.core-relaxed-lifetime.v1',status='RUNNING',scientific_admission='HOLD',
                 started_utc=datetime.now(timezone.utc).isoformat(),steps=[],diagnostics={},
                 n=args.n,varied_orbitals=varied,core_receipt_sha256=core_receipt_sha256,
                 two_s_method3=args.two_s_method3,
                 radial_iteration_floor=args.radial_iteration_floor,
                 both_2p_spectroscopic=args.both_2p_spectroscopic,
                 resume_triplet=str(args.resume_triplet) if args.resume_triplet else None,
                 damping=args.damping,resume_ground=str(args.resume_ground) if args.resume_ground else None,
                 cycles=args.cycles,resume_quintet=str(args.resume_quintet) if args.resume_quintet else None,
                 model=f'Existing n{args.n} full2s SD CSFs and orbitals as initializer; reoptimize {varied} radial functions with fixed1s only, fine439grid/original ACCY; strict orbital stop/fresh restart; DC only, no observed inputs')
    selected={}
    try:
        for state,spec in STATES.items():
            seed=args.core/('CV_2s_DC_'+state)
            for verify in (False,True):
                out=args.out/(state+('_verify' if verify else ''))
                out.mkdir()
                source=args.out/state if verify else seed
                for file in ('isodata','state.w','state.c'):
                    shutil.copyfile(source/file,out/file)
                resumed=args.resume_ground if state=='ground' else args.resume_quintet if state=='quintet' else args.resume_triplet
                if not verify and resumed:
                    if ((resumed/'rcsf.inp').read_bytes()!=(out/'state.c').read_bytes() or
                            (resumed/'isodata').read_bytes()!=(out/'isodata').read_bytes()):
                        raise ValueError('Same-layer resume input identity differs')
                    shutil.copyfile(resumed/'rwfn.out',out/'state.w')
                shutil.copyfile(out/'state.c',out/'rcsf.inp');shutil.copyfile(out/'state.w',out/'rwfn.inp')
                receipt['steps'].append(run(args.grasp_bin/'rangular','y\n',out,1800))
                controls=scf_controls(state,args.n,args.cycles,args.damping,verify,
                                      two_s_method3=args.two_s_method3,
                                      radial_iteration_floor=args.radial_iteration_floor,
                                      both_2p_spectroscopic=args.both_2p_spectroscopic)
                receipt['steps'].append(run(args.strict_bin/'rmcdhf',controls,out,43200))
                if any(s['semantic_errors'] for s in receipt['steps']):
                    raise ValueError('Native numerical refusal; no state promotion')
                d=scf_diagnostics((out/'rmcdhf.stdout').read_text(),(out/'rmcdhf.stderr').read_text())
                if d['iteration_limit_hit'] or not d['all_reported_functions_below_ACC']:
                    raise ValueError('Core-relaxed orbital criterion unmet')
                receipt['diagnostics'][out.name]=d
                for src,ext in [('rwfn.out','w'),('rmix.out','m')]:
                    shutil.copyfile(out/src,out/f'state.{ext}')
                receipt['steps'].append(run(args.grasp_bin/'jj2lsj','state\nn\nn\ny\n',out,600))
                identity((out/'state.lsj.lbl').read_text(),state)
            selected[state]=out
        require_triplet_quintet_order((selected['triplet']/'state.lsj.lbl').read_text(),(selected['quintet']/'state.lsj.lbl').read_text())
        for name,lower,operators,same in [('to_ground',selected['ground'],'E1,M2,E3',False),
                                          ('to_5S',selected['quintet'],'M1,E2,M3',True)]:
            receipt['steps']+=transitions(args.grasp_bin,lower,selected['triplet'],args.out/name,operators,same)
        receipt['status']='CORE_RELAXED_MODEL_LIFETIME_COMPLETED_REVIEW_REQUIRED'
    except Exception as exc:
        receipt.update(status='FAILED_RETAINED',error=f'{type(exc).__name__}: {exc}')
        raise
    finally:
        receipt['files']={str(p.relative_to(args.out)):dict(sha256=sha(p),bytes=p.stat().st_size)
                          for p in sorted(args.out.rglob('*')) if p.is_file()}
        (args.out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')


if __name__=='__main__':
    main()
