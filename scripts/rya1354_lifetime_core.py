"""Balanced fixed-orbital core/Hamiltonian controls for the triplet decay."""
import argparse
from datetime import datetime,timezone
import json
from pathlib import Path
import platform
import shutil

from rya1354_solver_smoke import run,sha
from rya1354_lifetime_ladder import STATES,identity
from rya1354_triplet_correlation import transitions,require_triplet_quintet_order


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--grasp-bin',type=Path,required=True)
    ap.add_argument('--ladder',type=Path,required=True)
    ap.add_argument('--n',type=int,choices=(3,4,5),required=True)
    ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args()
    if platform.node().lower().split('.')[0]!='sirius':
        raise RuntimeError('Atomic solvers belong on Sirius')
    args.grasp_bin=args.grasp_bin.resolve()
    args.out.mkdir(parents=True,exist_ok=False)
    shutil.copyfile(__file__,args.out/'driver.py')
    receipt=dict(schema='rya1354.balanced-lifetime-core.v1',status='RUNNING',scientific_admission='HOLD',
                 started_utc=datetime.now(timezone.utc).isoformat(),steps=[],n=args.n,
                 model='Fixed balanced strict orbitals; native RCI interpolation grid; identical VV DC control, full 2s-hole SD/DC, same core CSFs plus static-limit Breit. 1s frozen, point stationary Z8, recoil off. No observed-energy substitution.')
    try:
        for case,core,physics in [('DC',False,'n\nn\nn\nn\nn\nn\n'),
                                  ('CV_2s_DC',True,'n\nn\nn\nn\nn\nn\n'),
                                  ('CV_2s_Breit',True,'y\ny\n1.d-6\nn\nn\nn\nn\n')]:
            selected={}
            for state,spec in STATES.items():
                seed=args.ladder/f'n{args.n}_{state}_verify'
                out=args.out/(case+'_'+state)
                out.mkdir()
                identity((seed/'state.lsj.lbl').read_text(),state)
                for file in ('isodata','state.w'):
                    shutil.copyfile(seed/file,out/file)
                if core:
                    active=','.join(f'{args.n}{l}' for l in 'spdfg'[:args.n])
                    config=spec['config'].format(occupation='*').replace('2s(2,i)','2s(2,*)')
                    receipt['steps'].append(run(args.grasp_bin/'rcsfgenerate',f'*\n0\n{config}\n\n{active}\n{spec["js"]}\n2\nn\n',out,1800))
                    shutil.copyfile(out/'rcsf.out',out/'state.c')
                else:
                    shutil.copyfile(seed/'state.c',out/'state.c')
                shutil.copyfile(out/'state.c',out/'rcsf.inp')
                roots=spec['roots'].removesuffix('5\n') if state=='ground' else spec['roots']
                for exe,stdin in [('rangular','y\n'),('rci',f'y\nstate\n{physics}{roots}'),('jj2lsj','state\ny\nn\ny\n')]:
                    receipt['steps'].append(run(args.grasp_bin/exe,stdin,out,43200))
                    if receipt['steps'][-1]['semantic_errors']:
                        raise ValueError('Native scientific refusal')
                if not (out/'state.cm').is_file():
                    raise ValueError('Native CI mixing file absent')
                identity((out/'state.lsj.lbl').read_text(),state)
                selected[state]=out
            require_triplet_quintet_order((selected['triplet']/'state.lsj.lbl').read_text(),(selected['quintet']/'state.lsj.lbl').read_text())
            for name,lower,operators,same in [('to_ground',selected['ground'],'E1,M2,E3',False),
                                              ('to_5S',selected['quintet'],'M1,E2,M3',True)]:
                receipt['steps']+=transitions(args.grasp_bin,lower,selected['triplet'],args.out/(case+'_'+name),operators,same,ci=True)
            receipt['last_completed_case']=case
            (args.out/'progress.json').write_text(json.dumps(receipt,indent=2)+'\n')
        receipt['status']='BALANCED_CORE_LIFETIMES_COMPLETED_REVIEW_REQUIRED'
    except Exception as exc:
        receipt.update(status='FAILED_RETAINED',error=f'{type(exc).__name__}: {exc}')
        raise
    finally:
        receipt['files']={str(p.relative_to(args.out)):dict(sha256=sha(p),bytes=p.stat().st_size)
                          for p in sorted(args.out.rglob('*')) if p.is_file()}
        (args.out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')


if __name__=='__main__':
    main()
