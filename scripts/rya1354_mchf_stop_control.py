"""Compare original/strict MCHF updates on the native upstream Si fixture."""
import argparse
import json
from pathlib import Path
import platform
import shutil

from rya1354_solver_smoke import run,sha


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--atsp',type=Path,required=True)
    ap.add_argument('--strict-bin',type=Path,required=True)
    ap.add_argument('--fixture',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args()
    if platform.node().lower().split('.')[0]!='sirius':
        raise RuntimeError('Native solvers belong on Sirius')
    args.atsp=args.atsp.resolve()
    args.strict_bin=args.strict_bin.resolve()
    args.out.mkdir(parents=True,exist_ok=False)
    shutil.copyfile(__file__,args.out/'driver.py')
    receipt=dict(schema='rya1354.mchf-stop-control.v1',scope='ONE_CYCLE_SOFTWARE_CONTROL_NOT_CONVERGED_ATOM',steps=[],status='RUNNING')
    controls=(args.fixture/'mchf.stdin').read_text()
    if controls.count('500,0')!=1:
        raise ValueError('Unexpected original fixture input')
    controls=controls.replace('500,0','1,0')
    try:
        for case,binpath in [('original',args.atsp/'bin'),('strict',args.strict_bin)]:
            out=args.out/case
            out.mkdir()
            for file in ('cfg.inp','cfg.h','c.lst','ico.lst','ih.01.lst','yint.lst','wfn.inp'):
                shutil.copyfile(args.fixture/file,out/file)
            receipt['steps'].append(run(binpath/'mchf',controls,out,600))
        matches={file:(args.out/'original'/file).read_bytes()==(args.out/'strict'/file).read_bytes()
                 for file in ('wfn.out','1De.l')}
        if not all(matches.values()):
            raise ValueError('Strict variant changes updates before stopping divergence')
        receipt.update(status='ONE_CYCLE_ORBITAL_AND_MIXING_BYTES_IDENTICAL',matches=matches,
                       warning='Deliberate one-cycle cap, not convergence validation. Native old maximum-cycle wording is misleading; no atom accepted from these controls.')
    except Exception as exc:
        receipt.update(status='FAILED_RETAINED',error=f'{type(exc).__name__}: {exc}')
        raise
    finally:
        receipt['files']={str(p.relative_to(args.out)):dict(sha256=sha(p),bytes=p.stat().st_size)
                          for p in sorted(args.out.rglob('*')) if p.is_file()}
        (args.out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')


if __name__=='__main__':
    main()
