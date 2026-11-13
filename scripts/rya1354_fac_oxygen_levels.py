"""FAC neutral-O level-identity feasibility only; no transition/lifetime output."""
import argparse
from collections import Counter
from datetime import datetime,timezone
import hashlib
import json
from pathlib import Path
import platform
import re
import shutil
import subprocess

from rya1354_fac_bootstrap import PIN,level_inventory,sha


def input_text():
    return """SetAtom('O')
SetBreit(0)
SetSE(0)
SetVP(0)
SetMS(0,0)
Closed('1s')
Config('ground', '2s2 2p4')
Config('odd3s', '2s2 2p3 3s1')
ConfigEnergy(0)
SetOptimizePrint(1)
OptimizeRadial(['ground'])
ConfigEnergy(1)
Structure('o.lev.b', ['ground', 'odd3s'])
MemENTable('o.lev.b')
PrintTable('o.lev.b', 'o.lev', 1)
"""


def inspect_inventory(text):
    rows=level_inventory(text)
    expected=Counter({(0,0):2,(0,2):1,(0,4):2,(1,0):1,(1,2):4,(1,4):4,(1,6):1})
    if (re.search(r'^O Z\s*=\s*8\.0\s*$',text,re.M) is None
            or set(re.findall(r'^NELE\s*=\s*(\d+)',text,re.M))!={'8'}
            or Counter((r['parity'],r['two_J']) for r in rows)!=expected):
        raise ValueError('Neutral-O electron/parity/J inventory differs')
    for row in rows:
        # Pinned ConstructLevelName/IsClosedShellNR may omit a full2s shell.
        label=row['labels'][1].removeprefix('2s2.')
        if (row['parity']==0 and label!='2p4') or (row['parity']==1 and label!='2p3.3s1'):
            raise ValueError('Neutral-O dominant configuration label differs')
    return rows


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    for name in ('source','bootstrap','out'):
        ap.add_argument('--'+name,type=Path,required=True)
    args=ap.parse_args()
    if platform.node().lower().split('.')[0]!='sirius':
        raise RuntimeError('Atomic calculations belong on Sirius')
    args.source=args.source.resolve();args.bootstrap=args.bootstrap.resolve();args.out=args.out.resolve()
    head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=args.source,text=True).strip()
    build=json.loads((args.bootstrap/'receipt.json').read_text());binary=args.source/'sfac/sfac'
    if (head!=PIN or build['pin']!=PIN or build['status']!='FAC_BUILD_AND_UPSTREAM_INVENTORY_COMPLETED_NOT_O_VALIDATED'
            or sha(binary)!=build['binary']['sha256'] or build['tracked_source_diff']):
        raise ValueError('FAC build/pin/executable qualification differs')
    args.out.mkdir(parents=True,exist_ok=False)
    shutil.copyfile(__file__,args.out/'driver.py')
    shutil.copyfile(Path(__file__).with_name('rya1354_fac_bootstrap.py'),args.out/'rya1354_fac_bootstrap.py')
    (args.out/'input.sf').write_text(input_text())
    receipt=dict(schema='rya1354.fac-oxygen-levels.v17',status='RUNNING',pin=PIN,
        started_utc=datetime.now(timezone.utc).isoformat(),binary_sha256=sha(binary),
        bootstrap_receipt_sha256=sha(args.bootstrap/'receipt.json'),scientific_admission='HOLD',
        exposure='EXPOSED_NONBLIND',adoption_allowed=False,accepted_lifetime=None,
        uncertainty=None,covariance=None,calibrated_sigma=None,
        scope='Minimal 1s2 2s2 2p4 plus 1s2 2s2 2p3 3s neutral-O level inventory only. No transition/gf/lifetime, LS-term admission, GRASP matching or calibrated accuracy.',
        physics='FAC local-central-potential relativistic CI; ground-only OptimizeRadial with native defaults/diagnostic printing; theory-only ConfigEnergy workflow. Breit/SE/VP and normal/specific mass shifts explicitly disabled. Native standard nuclear parameters retained; no point-nucleus matching claim.')
    try:
        with (args.out/'sfac.stdout').open('xb') as out,(args.out/'sfac.stderr').open('xb') as err:
            process=subprocess.run([str(binary),'input.sf'],cwd=args.out,stdout=out,stderr=err,timeout=300,check=False)
        receipt['returncode']=process.returncode
        if process.returncode:
            raise RuntimeError('FAC neutral-O exit '+str(process.returncode))
        rows=inspect_inventory((args.out/'o.lev').read_text())
        receipt['levels']=rows
        receipt['status']='FAC_O_LEVEL_INVENTORY_COMPLETED_REVIEW_REQUIRED'
    except Exception as exc:
        receipt.update(status='FAILED_RETAINED',error=type(exc).__name__+': '+str(exc))
        raise
    finally:
        receipt['files']={str(p.relative_to(args.out)):dict(sha256=sha(p),bytes=p.stat().st_size)
                          for p in sorted(args.out.rglob('*')) if p.is_file()}
        (args.out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')


if __name__=='__main__':
    main()
