"""Sirius-only zero-substitution UV/IR and native M1/E2 operator exploration.

These are reference-configuration diagnostics, not converged rates. Preserve all
native state labels for review; do not identify transitions from wavelengths.
"""
from pathlib import Path
import argparse
import json
import platform
import shutil
from datetime import datetime, timezone
from rya1354_solver_smoke import run, sha


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--grasp-bin',type=Path,required=True)
    ap.add_argument('--seed',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args()
    args.grasp_bin=args.grasp_bin.resolve()
    args.seed=args.seed.resolve()
    if platform.node().lower().split('.')[0]!='sirius':
        raise RuntimeError('Solver calculations belong on Sirius')
    args.out.mkdir(parents=True,exist_ok=False)
    receipt=dict(schema='rya1354.coverage-smoke.v1',status='RUNNING',scientific_admission='HOLD',
                 started_utc=datetime.now(timezone.utc).isoformat(),script_sha256=sha(Path(__file__)),
                 model='Zero substitutions; point static Z=8; Dirac-Coulomb; all orbitals varied except triplet 3p freezes n<3 at quintet seed; native energies',
                 steps=[],warning='Native operator demonstration only; LS identities and numerical convergence require review; no lifetimes or uncertainty inferred')
    cases=[('ground','1s(2,i)2s(2,i)2p(4,i)','2p','0,4','1,2\n1\n1,2\n5\n'),
           ('triplet_3s','1s(2,i)2s(2,i)2p(3,i)3s(1,i)','3s,3p','2,2','1\n'),
           ('triplet_3p','1s(2,i)2s(2,i)2p(3,i)3p(1,i)','3s,3p','0,4','1\n2\n2\n5\n'),
           ('quintet_3d','1s(2,i)2s(2,i)2p(3,i)3d(1,i)','3s,3p,3d','0,8','1\n1\n1\n1\n1\n5\n'),
           ('quintet_4d','1s(2,i)2s(2,i)2p(3,i)4d(1,i)','4s,4p,4d','0,8','1\n1\n1\n1\n1\n5\n')]
    try:
        for name,config,active,js,states in cases:
            out=args.out/name
            out.mkdir()
            for exe,stdin in [('rnucleus','8\n0\n0\n0\n0\n0\n'),
                              ('rcsfgenerate',f'*\n0\n{config}\n\n{active}\n{js}\n0\nn\n')]:
                receipt['steps'].append(run(args.grasp_bin/exe,stdin,out,3600))
            shutil.copyfile(out/'rcsf.out',out/'rcsf.inp')
            init='y\n2\n*\n'
            if name in ('triplet_3s','triplet_3p'):
                seed_case='grasp_3s_J2' if name=='triplet_3s' else 'grasp_3p_J123'
                shutil.copyfile(args.seed/seed_case/'state.w',out/'previous.w')
                init='y\n1\nprevious.w\n*\n'
            varied='3*' if name=='triplet_3p' else '*'
            for exe,stdin in [('rangular','y\n'),('rwfnestimate',init),
                              ('rmcdhf',f'y\n{states}{varied}\n{varied}\n100\n')]:
                receipt['steps'].append(run(args.grasp_bin/exe,stdin,out,3600))
            for src,ext in [('rcsf.inp','c'),('rmix.out','m'),('rwfn.out','w')]:
                shutil.copyfile(out/src,out/f'state.{ext}')
            receipt['steps'].append(run(args.grasp_bin/'jj2lsj','state\nn\nn\ny\n',out,3600))
        for name,lower,upper,operators in [
                ('ground_M1_E2',args.out/'ground',args.out/'ground','M1,E2'),
                ('UV_E1',args.out/'ground',args.out/'triplet_3s','E1'),
                ('triplet_8446_E1',args.out/'triplet_3s',args.out/'triplet_3p','E1'),
                ('IR_3d_E1',args.seed/'grasp_3p_J123',args.out/'quintet_3d','E1'),
                ('optical_4d_E1',args.seed/'grasp_3p_J123',args.out/'quintet_4d','E1')]:
            out=args.out/name
            out.mkdir()
            for tag,source in [('lower',lower),('upper',upper)]:
                for ext in ('c','m','w'):
                    shutil.copyfile(source/f'state.{ext}',out/f'{tag}.{ext}')
            shutil.copyfile(lower/'isodata',out/'isodata')
            for exe,stdin in [('rbiotransform','y\nn\nlower\nupper\ny\n'),
                              ('rtransition',f'y\nn\nlower\nupper\n{operators}\n'+('n\nn\n' if 'M1' in operators else ''))]:
                receipt['steps'].append(run(args.grasp_bin/exe,stdin,out,3600))
        receipt['status']='NATIVE_OPERATORS_COMPLETED_REVIEW_REQUIRED'
    except Exception as exc:
        receipt.update(status='FAILED_RETAINED',error=f'{type(exc).__name__}: {exc}')
        raise
    finally:
        receipt['finished_utc']=datetime.now(timezone.utc).isoformat()
        receipt['files']={str(p.relative_to(args.out)):{'sha256':sha(p),'bytes':p.stat().st_size}
                          for p in sorted(args.out.rglob('*')) if p.is_file()}
        (args.out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')


if __name__=='__main__':
    main()
