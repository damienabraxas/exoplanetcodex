"""Replay pinned serial MCHF Si 1D n4 fixture, preserving all native evidence."""
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
    ap.add_argument('--atsp',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args()
    if platform.node().lower().split('.')[0]!='sirius':
        raise RuntimeError('Solver calculations belong on Sirius')
    args.atsp=args.atsp.resolve()
    args.out.mkdir(parents=True,exist_ok=False)
    shutil.copyfile(Path(__file__),args.out/'driver.py')
    receipt=dict(schema='rya1354.modern-mchf-native-control.v1',scientific_admission='HOLD',status='RUNNING',
                 started_utc=datetime.now(timezone.utc).isoformat(),steps=[])
    try:
        for source,target in [('examples/files_c/1De1.4.c','cfg.inp'),('examples/e1/e1.14_3.w','wfn.inp'),
                              ('examples/e1/e1.14_4_1D.l','expected.l')]:
            shutil.copyfile(args.atsp/source,args.out/target)
        receipt['steps'].append(run(args.atsp/'bin/nonh','',args.out,600))
        stdin='1D,14.\n1,1\nall\n1s 2s 2p 3s 3p\ny\nn\ny\nn\n.f.,1.e-9,1.e-9\nn\n500,0\ny\n'
        receipt['steps'].append(run(args.atsp/'bin/mchf',stdin,args.out,1800))
        output=args.out/'1De.l'
        if not output.is_file():
            raise ValueError('Native eigenvector output absent')
        energy=lambda p:float(re.search(r'^\s*1\s+(-\d+\.\d+)',p.read_text(),re.M).group(1))
        receipt.update(status='NATIVE_CONTROL_COMPLETED_NUMERICAL_COMPARISON_FOR_REVIEW',
                       reference_energy_hartree=energy(args.out/'expected.l'),
                       computed_energy_hartree=energy(output),
                       energy_difference_hartree=energy(output)-energy(args.out/'expected.l'),
                       byte_identical_to_reference=output.read_bytes()==(args.out/'expected.l').read_bytes())
    except Exception as exc:
        receipt.update(status='FAILED_RETAINED',error=f'{type(exc).__name__}: {exc}')
        raise
    finally:
        receipt['files']={str(p.relative_to(args.out)):{'sha256':sha(p),'bytes':p.stat().st_size}
                          for p in sorted(args.out.rglob('*')) if p.is_file()}
        (args.out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')


if __name__=='__main__':
    main()
