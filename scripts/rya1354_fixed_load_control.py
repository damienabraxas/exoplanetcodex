"""Prove native loading normalization separately from optimization; Sirius only."""
import argparse
from datetime import datetime,timezone
import json
from pathlib import Path
import platform
import shutil
import sys
from types import SimpleNamespace

from rya1354_lifetime_staged import load_baseline,require_frozen
from rya1354_solver_smoke import sha
from rya1354_wave_probe import decode


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--grasp-bin',type=Path,required=True)
    ap.add_argument('--strict-bin',type=Path,required=True)
    ap.add_argument('--failed-run',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args()
    if platform.node().lower().split('.')[0]!='sirius':
        raise RuntimeError('Native loading controls belong on Sirius')
    for attr in ('grasp_bin','strict_bin','failed_run','out'):
        setattr(args,attr,getattr(args,attr).resolve())
    failed=json.loads((args.failed_run/'receipt.json').read_text())
    if failed['status']!='FAILED_RETAINED' or failed['n']!=4:
        raise ValueError('Unexpected failed staged source')
    source=args.failed_run/'fixed_2s_ground'
    for name in ('isodata','rcsf.inp','rwfn.inp','rwfn.out'):
        path=source/name;record=failed['files'][str(path.relative_to(args.failed_run))]
        if sha(path)!=record['sha256'] or path.stat().st_size!=record['bytes']:
            raise ValueError('Native failed-source identity differs')
    args.out.mkdir(parents=True,exist_ok=False)
    for name in ('isodata','rcsf.inp','rwfn.inp'):
        shutil.copyfile(source/name,args.out/name)
    shutil.copyfile(source/'rwfn.out',args.out/'failed_rwfn.out')
    shutil.copyfile(__file__,args.out/'driver.py')
    for name in ('rya1354_lifetime_staged','rya1354_lifetime_relax'):
        shutil.copyfile(sys.modules[name].__file__,args.out/(name+'.py'))
    receipt=dict(schema='rya1354.native-fixed-load-control.v10',status='RUNNING',steps=[],
        captured_utc=datetime.now(timezone.utc).isoformat(),scientific_admission='HOLD',
        adoption_allowed=False,accepted_lifetime=None,failed_receipt_sha256=sha(args.failed_run/'receipt.json'))
    try:
        baseline=load_baseline(SimpleNamespace(grasp_bin=args.grasp_bin,strict_bin=args.strict_bin,damping=.5),
                               args.out,'ground',receipt)
        require_frozen(baseline,args.out/'failed_rwfn.out',True)
        raw=decode((args.out/'rwfn.inp').read_bytes());loaded=decode(baseline.read_bytes())
        receipt['fixed_records']=[dict(n=k[0],kappa=k[1],raw_record_identical=raw[k]==loaded[k],
            max_P_change=max(abs(a-b) for a,b in zip(raw[k]['p'],loaded[k]['p'])),
            normalized_record_identical_to_failed_output=True) for k in ((1,-1),(2,-1))]
        receipt.update(status='NATIVE_LOADING_NORMALIZATION_VERIFIED_NO_STATE_ADMISSION',
            interpretation='Same-grid native LODRWF/INTRPQ load renormalizes incoming orbitals by roundoff. All-fixed native baseline exactly reproduces fixed1s/2s records in rejected optimization. Compare against that baseline, not raw pre-load serialization; no tolerance or SCF-criterion change.')
    except Exception as exc:
        receipt.update(status='FAILED_RETAINED',error=type(exc).__name__+': '+str(exc));raise
    finally:
        receipt['files']={str(p.relative_to(args.out)):dict(sha256=sha(p),bytes=p.stat().st_size)
                          for p in sorted(args.out.rglob('*')) if p.is_file()}
        (args.out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')


if __name__=='__main__':
    main()
