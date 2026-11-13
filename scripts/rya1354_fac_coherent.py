"""Separately declared coherent inner/outer FAC numerical controls; no rates."""
import argparse
from datetime import datetime,timezone
import json
from pathlib import Path
import platform
import shutil
import subprocess

from rya1354_fac_bootstrap import PIN,sha
from rya1354_fac_controls import input_text as outer_input,diagnostics
from rya1354_fac_oxygen_levels import inspect_inventory

CASES=(('coherent10',.3,10),('coherent100',.03,100))


def input_text(multiplier,scale):
    inner=''.join(f"SetOption('orbital:{key}',{value/scale:.17g})\n" for key,value in
                  [('enerelerr',1e-4),('enerelerr1',1e-8),('eneabserr',1e-5)])
    text=outer_input(multiplier,256)
    return text.replace('SetOptimizeTolerance',inner+'SetOptimizeTolerance',1)


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
        raise ValueError('FAC pinned source/executable differs')
    args.out.mkdir(parents=True,exist_ok=False)
    for name in ('rya1354_fac_coherent.py','rya1354_fac_controls.py','rya1354_fac_bootstrap.py','rya1354_fac_oxygen_levels.py'):
        shutil.copyfile(Path(__file__).with_name(name),args.out/name)
    for name in ('faclib/orbital.c','faclib/radial.c','faclib/consts.h'):
        dest=args.out/'source_documents'/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(args.source/name,dest)
    receipt=dict(schema='rya1354.fac-coherent-controls.v18',status='RUNNING',pin=PIN,
        binary_sha256=sha(binary),started_utc=datetime.now(timezone.utc).isoformat(),cases={},
        scientific_admission='HOLD',exposure='EXPOSED_NONBLIND',adoption_allowed=False,
        accepted_lifetime=None,uncertainty=None,covariance=None,calibrated_sigma=None,
        scope='Separate cold-start inner/outer tolerance sensitivity, no failed-wave input, source patch, rates, LS-term mapping or calibrated accuracy.')
    try:
        for name,multiplier,scale in CASES:
            out=args.out/name;out.mkdir();(out/'input.sf').write_text(input_text(multiplier,scale))
            with (out/'sfac.stdout').open('xb') as stdout,(out/'sfac.stderr').open('xb') as stderr:
                result=subprocess.run([str(binary),'input.sf'],cwd=out,stdout=stdout,stderr=stderr,timeout=300,check=False)
            diagnostic=diagnostics((out/'sfac.stdout').read_text(),multiplier,256)
            if result.returncode or not diagnostic['printed_rule_verified']:
                raise ValueError('Coherent FAC numerical case failed: '+name)
            levels=inspect_inventory((out/'o.lev').read_text())
            exports=sorted(p.name for p in out.glob('o_export*'))
            if not exports or not (out/'o.basis').is_file() or not (out/'o.pot').is_file():
                raise ValueError('Coherent FAC basis/potential export absent')
            receipt['cases'][name]=dict(returncode=result.returncode,diagnostics=diagnostic,levels=levels,
                inner_scale=scale,exports=exports)
            (args.out/'progress.json').write_text(json.dumps(receipt,indent=2)+'\n')
        receipt['status']='FAC_COHERENT_CONTROLS_COMPLETED_REVIEW_REQUIRED'
    except Exception as exc:
        receipt.update(status='FAILED_RETAINED',error=type(exc).__name__+': '+str(exc))
        raise
    finally:
        receipt['files']={str(p.relative_to(args.out)):dict(sha256=sha(p),bytes=p.stat().st_size)
                          for p in sorted(args.out.rglob('*')) if p.is_file()}
        (args.out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')


if __name__=='__main__':
    main()
