"""Require actual native orbital criterion before radiation postprocessing."""
import argparse
from datetime import datetime,timezone
import json
from pathlib import Path
import platform
import shutil

from rya1354_solver_smoke import sha
from rya1354_triplet_correlation import transitions
from rya1354_two_method_summary import scf_diagnostics


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--grasp-bin',type=Path,required=True)
    ap.add_argument('--runs',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--source',type=Path,help='Explicit source for a single selected n')
    ap.add_argument('--n',type=int,nargs='+',choices=(3,4,5),default=[3,4,5])
    args=ap.parse_args()
    if platform.node().lower().split('.')[0]!='sirius':
        raise RuntimeError('Atomic solvers belong on Sirius')
    args.grasp_bin=args.grasp_bin.resolve()
    if args.source and len(args.n)!=1:
        raise ValueError('Explicit source requires exactly one layer')
    args.out.mkdir(parents=True,exist_ok=False)
    shutil.copyfile(__file__,args.out/'driver.py')
    receipt=dict(schema='rya1354.strict-multipoles.v1',status='RUNNING',scientific_admission='HOLD',
                 started_utc=datetime.now(timezone.utc).isoformat(),steps=[],diagnostics={})
    try:
        for n in args.n:
            source=args.source or args.runs/f'strict_n{n}_v1'
            source_receipt=json.loads((source/'receipt.json').read_text())
            if source_receipt['status']!='NATIVE_STOP_COMPLETED_REVIEW_REQUIRED':
                raise ValueError('SCF run failed')
            d=scf_diagnostics((source/'rmcdhf.stdout').read_text(),(source/'rmcdhf.stderr').read_text(),
                              function_accy=source_receipt.get('native_ACC',.0625**6))
            if d['iteration_limit_hit'] or not d['all_reported_functions_below_ACC']:
                raise ValueError('Native orbital criterion not satisfied')
            receipt['diagnostics'][str(n)]=d
            receipt['steps']+=transitions(args.grasp_bin,source,source,args.out/f'n{n}_multipoles','M1,E2,M3,E4',True)
        receipt['status']='ORBITAL_QUALIFIED_MULTIPOLES_COMPLETED_REVIEW_REQUIRED'
    except Exception as exc:
        receipt.update(status='FAILED_RETAINED',error=f'{type(exc).__name__}: {exc}')
        raise
    finally:
        receipt['files']={str(p.relative_to(args.out)):dict(sha256=sha(p),bytes=p.stat().st_size)
                          for p in sorted(args.out.rglob('*')) if p.is_file()}
        (args.out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')


if __name__=='__main__':
    main()
