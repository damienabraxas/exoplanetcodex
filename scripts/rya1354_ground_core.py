"""Fixed-orbital ground-term core/Hamiltonian controls, no external fitting."""
import argparse
from datetime import datetime,timezone
import json
from pathlib import Path
import platform
import shutil

from rya1354_solver_smoke import run,sha


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--grasp-bin',type=Path,required=True)
    ap.add_argument('--seed',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--combined-only',action='store_true')
    ap.add_argument('--hamiltonian-only',action='store_true',help='Identical-space DC and static Breit only; no new core CSFs')
    args=ap.parse_args()
    if args.combined_only and args.hamiltonian_only:
        raise ValueError('Choose one control subset')
    if platform.node().lower().split('.')[0]!='sirius':
        raise RuntimeError('Atomic solvers belong on Sirius')
    args.grasp_bin=args.grasp_bin.resolve()
    args.out.mkdir(parents=True,exist_ok=False)
    shutil.copyfile(__file__,args.out/'driver.py')
    receipt=dict(schema='rya1354.ground-core.v1',status='RUNNING',scientific_admission='HOLD',
                 started_utc=datetime.now(timezone.utc).isoformat(),steps=[],
                 model='Fixed strict_n3 orbitals. Native RCI interpolation grid. Point stationary Z8; recoil off. No observed energies or rates fitted.')
    cases=[('DC',None,'n\nn\nn\nn\nn\nn\n'),
           ('Breit_static_limit',None,'y\ny\n1.d-6\nn\nn\nn\nn\n'),
           ('CV_2s_one_hole','1s(2,i)2s(2,1)2p(4,*)','n\nn\nn\nn\nn\nn\n'),
           ('CV_2s_two_holes','1s(2,i)2s(2,*)2p(4,*)','n\nn\nn\nn\nn\nn\n')]
    if args.combined_only:
        cases=[('CV_2s_two_holes_Breit_static_limit','1s(2,i)2s(2,*)2p(4,*)','y\ny\n1.d-6\nn\nn\nn\nn\n')]
    if args.hamiltonian_only:
        cases=cases[:2]
        receipt['model']='Fixed declared seed orbitals and identical CSFs; DC/static Breit only, no core expansion. Native RCI interpolation grid; point stationary Z8, recoil off. No observed energies/rates fitted.'
    try:
        for name,config,physics in cases:
            out=args.out/name
            out.mkdir()
            for filename in ('isodata','state.w'):
                shutil.copyfile(args.seed/filename,out/filename)
            if config:
                receipt['steps'].append(run(args.grasp_bin/'rcsfgenerate',f'*\n0\n{config}\n\n3s,3p,3d\n0,4\n2\nn\n',out,600))
                shutil.copyfile(out/'rcsf.out',out/'state.c')
            else:
                shutil.copyfile(args.seed/'state.c',out/'state.c')
            shutil.copyfile(out/'state.c',out/'rcsf.inp')
            for exe,stdin in [('rangular','y\n'),('rci',f'y\nstate\n{physics}1,2\n1\n1,2\n'),
                              ('jj2lsj','state\ny\nn\ny\n')]:
                receipt['steps'].append(run(args.grasp_bin/exe,stdin,out,7200))
                if receipt['steps'][-1]['semantic_errors']:
                    raise ValueError('Native scientific refusal')
            # Preserve CI ordinal/LS assignments; no state selection by proximity.
            import re
            blocks=re.findall(r'^\s*\d+\s+(\d+)\s+([+-])\s+[^\n]+\n([^\n]+)',(out/'state.lsj.lbl').read_text(),re.M)
            expected=[('0','3P'),('0','1S'),('1','3P'),('2','3P'),('2','1D')]
            if len(blocks)!=5 or any(j!=ej or p!='+' or '2p(4)' not in line or not line.rstrip().endswith('_'+term)
                                     for (j,p,line),(ej,term) in zip(blocks,expected)):
                raise ValueError('Ground LS inventory changed')
            transition=args.out/(name+'_multipoles')
            transition.mkdir()
            for tag in ('lower','upper'):
                for ext in ('c','cm','w'):
                    shutil.copyfile(out/f'state.{ext}',transition/f'{tag}.{ext}')
                shutil.copyfile(out/'state.lsj.lbl',transition/f'{tag}.lsj.lbl')
            shutil.copyfile(out/'isodata',transition/'isodata')
            for exe,stdin in [('rbiotransform','y\ny\nlower\nupper\ny\n'),
                              ('rtransition','y\ny\nlower\nupper\nM1,E2,M3,E4\nn\nn\n')]:
                receipt['steps'].append(run(args.grasp_bin/exe,stdin,transition,1800))
            if not (transition/'lower.upper.ct').is_file() or any(s['semantic_errors'] for s in receipt['steps']):
                raise ValueError('CI multipole table missing or refused')
        receipt['status']='GROUND_CORE_CONTROLS_COMPLETED_REVIEW_REQUIRED'
    except Exception as exc:
        receipt.update(status='FAILED_RETAINED',error=f'{type(exc).__name__}: {exc}')
        raise
    finally:
        receipt['files']={str(p.relative_to(args.out)):dict(sha256=sha(p),bytes=p.stat().st_size)
                          for p in sorted(args.out.rglob('*')) if p.is_file()}
        (args.out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')


if __name__=='__main__':
    main()
