"""Sirius-only initialization controls; no SCF or transition accuracy claim.

Replay the valid n5 CSF/previous-n4 orbitals in new directories, varying only
initializer and radial grid. Upstream getinfo.f90/genrwf.f90 own prompt ordering.
"""
import argparse
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import platform
import shutil
from rya1354_solver_smoke import run, sha


def initialization_input(method, points, spacing):
    rnt=math.exp(-65/16)/8
    # Nondefault settings, no debug, default summary filename; change grid only,
    # retain c; read old orbitals then initialize all missing orbitals explicitly.
    return (f'n\nn\n\ny\nn\ny\n{rnt:.17g}\n{spacing:.17g}\n0\n{points}\n'
            f'1\nprevious.w\n*\n{method}\n*\nn\n')


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--grasp-bin',type=Path,required=True)
    ap.add_argument('--prepared',type=Path,required=True,help='Retained valid n5 lower CSF/isodata/previous.w directory')
    ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args()
    if platform.node().lower().split('.')[0]!='sirius':
        raise RuntimeError('Atomic calculations belong on Sirius')
    args.grasp_bin=args.grasp_bin.resolve()
    args.out.mkdir(parents=True,exist_ok=False)
    shutil.copyfile(Path(__file__),args.out/'driver.py')
    receipt=dict(schema='rya1354.radial-initialization-controls.v1',scientific_admission='HOLD',
                 started_utc=datetime.now(timezone.utc).isoformat(),cases=[],
                 note='Initialization only, no SCF; no parameter chosen from gf agreement')
    for name,method,points,spacing in [('TF_220',2,220,.0625),('H_220',3,220,.0625),
                                      ('H_260',3,260,.0625),('H_440',3,440,.03125)]:
        out=args.out/name
        out.mkdir()
        for file in ('rcsf.inp','isodata','previous.w'):
            shutil.copyfile(args.prepared/file,out/file)
        case=dict(name=name,method=method,points=points,spacing=spacing,status='RUNNING')
        try:
            case['step']=run(args.grasp_bin/'rwfnestimate',initialization_input(method,points,spacing),out,180)
            output=(out/'rwfnestimate.stderr').read_text()+(out/'rwfnestimate.stdout').read_text()
            if case['step']['semantic_errors'] or 'All required subshell radial wavefunctions' not in output:
                raise ValueError('Native completion marker or scientific input refused')
            if not (out/'rwfn.inp').stat().st_size:
                raise ValueError('Empty wavefunction artifact')
            case['status']='INITIALIZATION_COMPLETED_NOT_SCF_VALIDATED'
        except (OSError,RuntimeError,ValueError) as exc:
            case.update(status='FAILED_RETAINED',error=f'{type(exc).__name__}: {exc}')
        receipt['cases'].append(case)
    receipt['files']={str(p.relative_to(args.out)):{'sha256':sha(p),'bytes':p.stat().st_size}
                      for p in sorted(args.out.rglob('*')) if p.is_file()}
    (args.out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print([(c['name'],c['status']) for c in receipt['cases']])


if __name__=='__main__':
    main()
