"""Pinned FAC stopping-rule replay and cold-start level/basis controls; no rates."""
import argparse
from datetime import datetime,timezone
from decimal import Decimal
import json
import math
from pathlib import Path
import platform
import shutil
import subprocess

from rya1354_fac_bootstrap import PIN,sha
from rya1354_fac_oxygen_levels import input_text as original_input,inspect_inventory

CASES=(('cap_canary',3.0,1),('explicit_default',3.0,256),('tight10',.3,256),('tight100',.03,256))


def input_text(multiplier,limit):
    return original_input().replace('SetOptimizePrint(1)',
        f'SetOptimizeTolerance({multiplier})\nSetOptimizeMaxIter({limit})\nSetOptimizePrint(1)')+"BasisTable('o.basis',0)\nBasisTable('o_export',10)\nSavePotential('o.pot')\n"


def diagnostics(text,multiplier,limit):
    """Replay printed energy/potential conditions; positive native marker required."""
    rows=[]
    names=('tol','atol','potential_change','tol0','tol1','atol0','ahx','dph0','dph1','sta0','sta1')
    for line in text.splitlines():
        if not line.startswith('optimize loop:'):
            continue
        f=line.split()
        if len(f)!=16:
            raise ValueError('FAC optimization row layout differs')
        row=dict(iteration=int(f[2]),native_miter=int(f[3]),muconv=int(f[4]))
        values=list(map(float,f[5:]))
        if not all(math.isfinite(x) for x in values):
            raise ValueError('Nonfinite FAC optimization quantity')
        row.update(zip(names,values))
        for key,scale in [('tol0','1e-4'),('tol1','1e-8'),('atol0','1e-5')]:
            value=f[5+names.index(key)]
            if Decimal(value)!=Decimal(str(multiplier))*Decimal(scale):
                raise ValueError('FAC threshold multiplier differs')
        rows.append(row)
    if len(rows)<2 or [r['iteration'] for r in rows]!=list(range(len(rows))):
        raise ValueError('FAC optimization cycle inventory differs')
    def numerical(row):
        sta=min(row['sta0'],row['sta1'])
        return ((row['tol']<=row['tol0']*sta and row['atol']<=row['atol0']*sta)
                or row['tol']<=row['tol1']*sta)
    previous,last=rows[-2:]
    qualified=(last['iteration']<limit and last['native_miter']==last['iteration'] and last['native_miter']>0
        and last['muconv']==1 and previous['muconv']==1 and previous['iteration']>7
        and previous['ahx']<=1e-5 and numerical(previous) and numerical(last))
    return dict(iterations=len(rows),last=last,previous=previous,printed_rule_verified=qualified,
        native_limit=limit,multiplier=multiplier,
        scope='Pinned FAC printed local-potential/energy stopping conditions plus native positive miter marker; not GRASP orbital ACCY or physical accuracy.')


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    for name in ('source','bootstrap','out'):
        ap.add_argument('--'+name,type=Path,required=True)
    args=ap.parse_args()
    if platform.node().lower().split('.')[0]!='sirius':
        raise RuntimeError('Atomic calculations belong on Sirius')
    args.source=args.source.resolve();args.bootstrap=args.bootstrap.resolve();args.out=args.out.resolve()
    build=json.loads((args.bootstrap/'receipt.json').read_text());binary=args.source/'sfac/sfac'
    if (subprocess.check_output(['git','rev-parse','HEAD'],cwd=args.source,text=True).strip()!=PIN
            or subprocess.check_output(['git','diff','--name-only'],cwd=args.source,text=True).strip()
            or build['status']!='FAC_BUILD_AND_UPSTREAM_INVENTORY_COMPLETED_NOT_O_VALIDATED'
            or sha(binary)!=build['binary']['sha256']):
        raise ValueError('FAC source/executable qualification differs')
    args.out.mkdir(parents=True,exist_ok=False)
    for name in ('rya1354_fac_controls.py','rya1354_fac_bootstrap.py','rya1354_fac_oxygen_levels.py'):
        shutil.copyfile(Path(__file__).with_name(name),args.out/name)
    for name in ('faclib/consts.h','faclib/radial.c','faclib/structure.c','sfac/sfac.c','doc/manual.tex'):
        dest=args.out/'source_documents'/name;dest.parent.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(args.source/name,dest)
    receipt=dict(schema='rya1354.fac-numerical-controls.v18',status='RUNNING',pin=PIN,
        started_utc=datetime.now(timezone.utc).isoformat(),binary_sha256=sha(binary),cases={},
        scientific_admission='HOLD',exposure='EXPOSED_NONBLIND',adoption_allowed=False,
        accepted_lifetime=None,uncertainty=None,covariance=None,calibrated_sigma=None,
        scope='Explicit cold-start FAC potential stopping controls and basis/mixing export only. No rates, LS-term admission, independent converged benchmark or accuracy calibration.')
    try:
        for name,multiplier,limit in CASES:
            out=args.out/name;out.mkdir();(out/'input.sf').write_text(input_text(multiplier,limit))
            with (out/'sfac.stdout').open('xb') as stdout,(out/'sfac.stderr').open('xb') as stderr:
                result=subprocess.run([str(binary),'input.sf'],cwd=out,stdout=stdout,stderr=stderr,timeout=300,check=False)
            text=(out/'sfac.stdout').read_text()+(out/'sfac.stderr').read_text()
            diagnostic=diagnostics(text,multiplier,limit)
            record=dict(returncode=result.returncode,diagnostics=diagnostic)
            if name=='cap_canary':
                if diagnostic['printed_rule_verified']:
                    raise ValueError('Capped FAC canary wrongly qualified')
                record.update(expected_rejection=True,levels_written=(out/'o.lev').is_file())
            else:
                if result.returncode or not diagnostic['printed_rule_verified']:
                    raise ValueError('FAC declared stopping control failed: '+name)
                record['levels']=inspect_inventory((out/'o.lev').read_text())
                record['exports']=sorted(p.name for p in out.glob('o_export*'))
                if not record['exports'] or not (out/'o.basis').is_file() or not (out/'o.pot').is_file():
                    raise ValueError('FAC basis/potential export missing')
            receipt['cases'][name]=record
            (args.out/'progress.json').write_text(json.dumps(receipt,indent=2)+'\n')
        receipt['status']='FAC_NUMERICAL_CONTROLS_COMPLETED_REVIEW_REQUIRED'
    except Exception as exc:
        receipt.update(status='FAILED_RETAINED',error=type(exc).__name__+': '+str(exc))
        raise
    finally:
        receipt['files']={str(p.relative_to(args.out)):dict(sha256=sha(p),bytes=p.stat().st_size)
                          for p in sorted(args.out.rglob('*')) if p.is_file()}
        (args.out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')


if __name__=='__main__':
    main()
