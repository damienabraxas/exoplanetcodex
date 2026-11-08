"""Complete angularly permitted one-photon operators for configured O I 3s 3S1.

Downward bound-level inventory: ground 2p4 (five levels), then 3s5S2.
Jupper=1 and Jlower<=2 limit k<=3. Parity permits E1/M2/E3 to ground
and M1/E2/M3 to metastable 5S2. This is a model-space inventory, not an
accuracy claim or certification of the physical completeness of the atom.
"""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import platform
import shutil
from rya1354_solver_smoke import run,sha


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--grasp-bin',type=Path,required=True)
    ap.add_argument('--coverage',type=Path,required=True)
    ap.add_argument('--seed',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args()
    if platform.node().lower().split('.')[0]!='sirius':
        raise RuntimeError('Solver calculations belong on Sirius')
    args.grasp_bin=args.grasp_bin.resolve()
    args.out.mkdir(parents=True,exist_ok=False)
    shutil.copyfile(Path(__file__),args.out/'driver.py')
    receipt=dict(schema='rya1354.3s3S1-decay-inventory.v1',status='RUNNING',scientific_admission='HOLD',
                 started_utc=datetime.now(timezone.utc).isoformat(),steps=[],
                 note='Fixed retained coverage_v8 wavefunctions; no orbitals fit to experiment; model-space one-photon completeness through k=3')
    upper=args.coverage/'triplet_3s'
    if '3s_3S' not in (upper/'state.lsj.lbl').read_text():
        raise ValueError('Upper LS identity absent')
    try:
        for name,lower,operators,extra in [
                ('to_ground',args.coverage/'ground','E1,M2,E3',''),
                ('to_5S',args.seed/'grasp_3s_J2','M1,E2,M3','n\nn\n')]:
            out=args.out/name
            out.mkdir()
            for tag,source in [('lower',lower),('upper',upper)]:
                for ext in ('c','m','w'):
                    shutil.copyfile(source/f'state.{ext}',out/f'{tag}.{ext}')
                shutil.copyfile(source/'state.lsj.lbl',out/f'{tag}.lsj.lbl')
            shutil.copyfile(lower/'isodata',out/'isodata')
            for exe,stdin in [('rbiotransform','y\nn\nlower\nupper\ny\n'),
                              ('rtransition',f'y\nn\nlower\nupper\n{operators}\n{extra}')]:
                receipt['steps'].append(run(args.grasp_bin/exe,stdin,out,600))
        if any(s['semantic_errors'] for s in receipt['steps']):
            raise ValueError('Native semantic refusal')
        receipt['status']='MODEL_SPACE_OPERATORS_COMPLETED_NOT_ACCURACY_VALIDATED'
    except Exception as exc:
        receipt.update(status='FAILED_RETAINED',error=f'{type(exc).__name__}: {exc}')
        raise
    finally:
        receipt['files']={str(p.relative_to(args.out)):{'sha256':sha(p),'bytes':p.stat().st_size}
                          for p in sorted(args.out.rglob('*')) if p.is_file()}
        (args.out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')


if __name__=='__main__':
    main()
