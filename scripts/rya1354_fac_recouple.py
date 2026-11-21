"""Unmodified GRASP angular recoupling of exact FAC exports; no new Hamiltonian."""
import argparse
from datetime import datetime,timezone
import json
from pathlib import Path
import platform
import shutil
import subprocess

from rya1354_solver_smoke import run,sha
from rya1354_fac_controls import diagnostics
from rya1354_fac_mapping import verify_transport

CONTROLS='state\ny\nn\nn\ny\n0\n0\ny\n1\n'
GRASP_PIN='9006157730a82ac839f2b4ff4e938bcba63a539e'


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    for name in ('source','grasp','isodata','out'):
        ap.add_argument('--'+name,type=Path,required=True)
    args=ap.parse_args()
    if platform.node().lower().split('.')[0]!='sirius':
        raise RuntimeError('Atomic representation transforms belong on Sirius')
    for name,value in vars(args).items():setattr(args,name,value.resolve())
    parent=json.loads((args.source/'receipt.json').read_text())
    if parent['status']!='FAC_COHERENT_CONTROLS_COMPLETED_REVIEW_REQUIRED':
        raise ValueError('Unqualified FAC numerical source')
    for name,record in parent['files'].items():
        p=args.source/name
        if sha(p)!=record['sha256'] or p.stat().st_size!=record['bytes']:
            raise ValueError('FAC numerical source artifact changed')
    if (subprocess.check_output(['git','rev-parse','HEAD'],cwd=args.grasp,text=True).strip()!=GRASP_PIN
            or subprocess.check_output(['git','diff','--name-only'],cwd=args.grasp,text=True).strip()):
        raise ValueError('GRASP angular converter source differs')
    audited={}
    for case,multiplier in [('coherent10',.3),('coherent100',.03)]:
        source=args.source/case
        d=diagnostics((source/'sfac.stdout').read_text(),multiplier,256)
        if not d['printed_rule_verified']:raise ValueError('FAC source stopping guard failed')
        audited[case]=verify_transport((source/'o.basis').read_text(),(source/'o_export_000.c').read_text(),
                                       (source/'o_export_000.cm').read_bytes())
    args.out.mkdir(parents=True,exist_ok=False)
    drivers=args.out/'drivers';drivers.mkdir()
    for p in Path(__file__).parent.glob('rya1354_*.py'):shutil.copyfile(p,drivers/p.name)
    reader=args.grasp/'src/appl/jj2lsj90'
    for name in ('getmixblock.f90','jj2lsj_code.f90','jj2lsj2K.f90'):
        dest=args.out/'converter_source'/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(reader/name,dest)
    receipt=dict(schema='rya1354.fac-recoupling.v19',status='RUNNING',started_utc=datetime.now(timezone.utc).isoformat(),
        source_receipt_sha256=sha(args.source/'receipt.json'),grasp_pin=GRASP_PIN,transport=audited,steps=[],
        scientific_admission='HOLD',exposure='EXPOSED_NONBLIND',adoption_allowed=False,
        accepted_lifetime=None,uncertainty=None,covariance=None,calibrated_sigma=None,
        scope='Exact FAC C/CM angular-representation conversion only, no orbital/Hamiltonian calculation or rate. Unmodified converter, all levels/full composition/transformation output, unique labeling disabled. Z8 isodata is converter support only, not a replacement/match of FAC physical nucleus.')
    try:
        for case in audited:
            source=args.source/case;out=args.out/case;out.mkdir()
            for original,dest in [('o_export_000.c','state.c'),('o_export_000.cm','state.cm')]:
                shutil.copyfile(source/original,out/dest)
            shutil.copyfile(args.isodata,out/'isodata')
            receipt['steps'].append(run(args.grasp/'build/bin/jj2lsj',CONTROLS,out,300))
            if receipt['steps'][-1]['semantic_errors']:
                raise ValueError('Native recoupling semantic refusal')
            for name in ('state.lsj.lbl','state.lsj.c','state.lsj.j','state.lsj.T'):
                if not (out/name).is_file() or not (out/name).stat().st_size:
                    raise ValueError('Full native recoupling artifact missing: '+name)
            if (out/'state.c').read_bytes()!=(source/'o_export_000.c').read_bytes() or (out/'state.cm').read_bytes()!=(source/'o_export_000.cm').read_bytes():
                raise ValueError('Recoupling changed original FAC C/CM input')
        receipt['status']='FAC_ANGULAR_TRANSFORM_COMPLETED_REVIEW_REQUIRED'
    except Exception as exc:
        receipt.update(status='FAILED_RETAINED',error=type(exc).__name__+': '+str(exc));raise
    finally:
        receipt['files']={str(p.relative_to(args.out)):dict(sha256=sha(p),bytes=p.stat().st_size)
                          for p in sorted(args.out.rglob('*')) if p.is_file()}
        (args.out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')


if __name__=='__main__':main()
