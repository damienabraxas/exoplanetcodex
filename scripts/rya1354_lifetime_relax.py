"""Reoptimize matched n3 core-expanded states; preserve fixed-orbital control."""
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


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--grasp-bin',type=Path,required=True)
    ap.add_argument('--strict-bin',type=Path,required=True)
    ap.add_argument('--core',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--damping',type=float,default=.9)
    ap.add_argument('--resume-ground',type=Path,help='Explicit same-CSF failed ground stage as an unqualified initializer only')
    ap.add_argument('--resume-quintet',type=Path,help='Explicit same-CSF quintet continuation initializer')
    ap.add_argument('--cycles',type=int,default=500)
    args=ap.parse_args()
    if not 0<args.damping<1 or args.cycles<1:
        raise ValueError('Damping must lie strictly between zero and one')
    if platform.node().lower().split('.')[0]!='sirius':
        raise RuntimeError('Atomic solvers belong on Sirius')
    args.grasp_bin=args.grasp_bin.resolve();args.strict_bin=args.strict_bin.resolve()
    args.out.mkdir(parents=True,exist_ok=False)
    shutil.copyfile(__file__,args.out/'driver.py')
    receipt=dict(schema='rya1354.core-relaxed-lifetime.v1',status='RUNNING',scientific_admission='HOLD',
                 started_utc=datetime.now(timezone.utc).isoformat(),steps=[],diagnostics={},
                 damping=args.damping,resume_ground=str(args.resume_ground) if args.resume_ground else None,
                 cycles=args.cycles,resume_quintet=str(args.resume_quintet) if args.resume_quintet else None,
                 model='Existing n3 full2s SD CSFs and orbitals as initializer; reoptimize n2/n3 radial functions with fixed1s only, fine439grid/original ACCY; strict orbital stop/fresh restart; DC only, no observed inputs')
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
                resumed=args.resume_ground if state=='ground' else args.resume_quintet if state=='quintet' else None
                if not verify and resumed:
                    if ((resumed/'rcsf.inp').read_bytes()!=(out/'state.c').read_bytes() or
                            (resumed/'isodata').read_bytes()!=(out/'isodata').read_bytes()):
                        raise ValueError('Same-layer resume input identity differs')
                    shutil.copyfile(resumed/'rwfn.out',out/'state.w')
                shutil.copyfile(out/'state.c',out/'rcsf.inp');shutil.copyfile(out/'state.w',out/'rwfn.inp')
                receipt['steps'].append(run(args.grasp_bin/'rangular','y\n',out,1800))
                spectroscopic='2s,2p'+(',3s' if state!='ground' else '')
                grid=f'n\nn\ny\nn\ny\n{math.exp(-65/16)/8:.17g}\n0.03125\n0\n439\ny\n{.0625**6:.17g}\n'
                nsic='n\n' if verify else 'y\n20\n'
                controls=f'{grid}{spec["roots"]}2*,3*\n{spectroscopic}\n{args.cycles}\ny\nn\nn\nn\ny\nn\n{-args.damping}\nn\n{nsic}n\nn\n1\n'
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
