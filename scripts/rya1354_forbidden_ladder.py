"""Ground-configuration VV-SD convergence with native M1/E2/M3/E4 operators.

All five lowest 2p4 levels share an EOL orbital set. Fixed n<3 radial functions
start from coverage_v8 ground HF; optimize new layers only. No observed energies
or transition strengths enter the optimization. Same-file duplicate directions
must not be summed. This is not a certified forbidden-line data product.
"""
import argparse
from datetime import datetime,timezone
import json
from pathlib import Path
import platform
import re
import shutil
from rya1354_solver_smoke import run,sha


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--grasp-bin',type=Path,required=True)
    ap.add_argument('--seed',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--n-max',type=int,choices=(3,4,5),default=5)
    ap.add_argument('--n-min',type=int,choices=(3,4,5),default=3)
    ap.add_argument('--resume',type=Path,help='Same-layer retained state directory for numerical SCF continuation')
    ap.add_argument('--scf-cycles',type=int,default=500)
    ap.add_argument('--fixed-damping',action='store_true',
                    help='Numerical stabilization only: native ODAMP=-damping-factor for varied orbitals; all other settings and tolerances unchanged')
    ap.add_argument('--damping-factor',type=float,default=.5)
    args=ap.parse_args()
    if platform.node().lower().split('.')[0]!='sirius':
        raise RuntimeError('Solver calculations belong on Sirius')
    args.grasp_bin=args.grasp_bin.resolve()
    if args.n_min>args.n_max or (args.n_min>3 and args.resume is None) or args.scf_cycles<=0:
        raise ValueError('Invalid layer range, missing resume state or nonpositive iteration count')
    if not 0<args.damping_factor<1:
        raise ValueError('Native fixed damping must lie strictly between zero and one')
    args.out.mkdir(parents=True,exist_ok=False)
    shutil.copyfile(Path(__file__),args.out/'driver.py')
    receipt=dict(schema='rya1354.forbidden-vv-sd.v1',status='RUNNING',scientific_admission='HOLD',
                 started_utc=datetime.now(timezone.utc).isoformat(),steps=[],
                 fixed_radial_damping=(-args.damping_factor if args.fixed_damping else None),
                 model='Point static Z8, Dirac-Coulomb; closed 1s2/2s2 occupations and fixed n<3 radial functions; valence SD; EOL five ground-configuration levels; new-layer optimization; no observed-level substitution')
    try:
        for n in range(args.n_min,args.n_max+1):
            out=args.out/f'n{n}_ground'
            out.mkdir()
            active=','.join(f'{n}{l}' for l in 'spdfg'[:n])
            for exe,stdin in [('rnucleus','8\n0\n0\n0\n0\n0\n'),
                              ('rcsfgenerate',f'*\n0\n1s(2,i)2s(2,i)2p(4,*)\n\n{active}\n0,4\n2\nn\n')]:
                receipt['steps'].append(run(args.grasp_bin/exe,stdin,out,600))
            shutil.copyfile(out/'rcsf.out',out/'rcsf.inp')
            resumed=(n==args.n_min and args.resume is not None)
            if resumed and (out/'rcsf.inp').read_bytes()!=(args.resume/'state.c').read_bytes():
                raise ValueError('Resume configuration bytes differ; numerical continuation refused')
            previous=(args.resume/'state.w' if resumed else
                      args.seed/'state.w' if n==3 else args.out/f'n{n-1}_ground'/'state.w')
            shutil.copyfile(previous,out/'previous.w')
            init='y\n1\nprevious.w\n*\n'+('' if resumed else '3\n*\n')
            scf=f'y\n1,2\n1\n1,2\n5\n{n}*\n\n{args.scf_cycles}\n'
            if args.fixed_damping:
                # getscd.f90: unchanged grid/ACCY, roots and standard weights;
                # modify only native ODAMP, retaining methods/node/orthogonality controls.
                scf=(f'n\nn\nn\nn\n1,2\n1\n1,2\n5\n{n}*\n\n{args.scf_cycles}\n'
                     f'y\nn\nn\nn\ny\nn\n{-args.damping_factor}\nn\nn\nn\nn\n1\n')
            for exe,stdin in [('rangular','y\n'),('rwfnestimate',init),
                              ('rmcdhf',scf)]:
                receipt['steps'].append(run(args.grasp_bin/exe,stdin,out,3600))
                if receipt['steps'][-1]['semantic_errors']:
                    raise ValueError(f'{exe}: native scientific refusal or SCF limit')
            for src,ext in [('rcsf.inp','c'),('rwfn.out','w'),('rmix.out','m')]:
                shutil.copyfile(out/src,out/f'state.{ext}')
            receipt['steps'].append(run(args.grasp_bin/'jj2lsj','state\nn\nn\ny\n',out,600))
            # Leading labels, not wavelength proximity, identify the five ASFs.
            blocks=re.findall(r'^\s*\d+\s+(\d+)\s+([+-])\s+[^\n]+\n([^\n]+)',
                              (out/'state.lsj.lbl').read_text(),re.M)
            expected=[('0','3P'),('0','1S'),('1','3P'),('2','3P'),('2','1D')]
            if len(blocks)!=5 or any(j!=ej or p!='+' or '2p(4)' not in line or not line.rstrip().endswith('_'+term)
                                     for (j,p,line),(ej,term) in zip(blocks,expected)):
                raise ValueError('Configured five-level identity inventory not established')
            transition=args.out/f'n{n}_multipoles'
            transition.mkdir()
            for tag in ('lower','upper'):
                for ext in ('c','w','m'):
                    shutil.copyfile(out/f'state.{ext}',transition/f'{tag}.{ext}')
            shutil.copyfile(out/'isodata',transition/'isodata')
            for exe,stdin in [('rbiotransform','y\nn\nlower\nupper\ny\n'),
                              ('rtransition','y\nn\nlower\nupper\nM1,E2,M3,E4\nn\nn\n')]:
                receipt['steps'].append(run(args.grasp_bin/exe,stdin,transition,1800))
            if not (transition/'lower.upper.t').is_file() or any(s['semantic_errors'] for s in receipt['steps']):
                raise ValueError('Native multipole result missing or input refused')
        receipt['status']='FORBIDDEN_MODEL_LADDER_COMPLETED_REVIEW_REQUIRED'
    except Exception as exc:
        receipt.update(status='FAILED_RETAINED',error=f'{type(exc).__name__}: {exc}')
        raise
    finally:
        receipt['files']={str(p.relative_to(args.out)):{'sha256':sha(p),'bytes':p.stat().st_size}
                          for p in sorted(args.out.rglob('*')) if p.is_file()}
        (args.out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')


if __name__=='__main__':
    main()
