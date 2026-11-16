"""Sequential remaining LS cases; unchanged qualified child, isolated process groups."""
import argparse
from datetime import datetime,timezone
import json
import os
from pathlib import Path
import platform
import shutil
import signal
import subprocess
import sys
from rya1354_fac_bootstrap import sha
from rya1354_fac_rates import verify_receipt
from rya1354_fac_radial_recouple import qualified_canary

CASES=('ground_coherent100','equal_coherent10','equal_coherent100')


def command(a,case):
    if case not in CASES:raise ValueError('Undeclared remaining case')
    cmd=[sys.executable,str(Path(__file__).with_name('rya1354_fac_radial_recouple.py'))]
    for key in ('grasp','first','identity','remaining','old_converted','canary'):cmd.extend(['--'+key,str(getattr(a,key))])
    return cmd+['--case',case,'--out',str(a.out/case)]


def stop_owned_group(pid):
    try:os.killpg(pid,signal.SIGTERM)
    except ProcessLookupError:pass


def main():
    ap=argparse.ArgumentParser()
    for name in ('grasp','first','identity','remaining','old_converted','canary','first_angular','out'):ap.add_argument('--'+name,type=Path,required=True)
    a=ap.parse_args()
    if platform.node().lower().split('.')[0]!='sirius':raise RuntimeError('Native angular work belongs on Sirius')
    for key,value in vars(a).items():setattr(a,key,value.resolve())
    first=verify_receipt(a.first_angular);canary=json.loads((a.canary/'report.json').read_text());qualified_canary(canary)
    if first['status']!='FAC_RADIAL_ANGULAR_REPLAY_COMPLETED_REVIEW_REQUIRED' or first['case']!='ground_coherent10' or first['report_sha256']!=sha(a.first_angular/'report.json') or sha(a.first_angular/'drivers/rya1354_fac_radial_ls_audit.py')!=sha(a.canary/'rya1354_fac_radial_ls_audit.py'):raise ValueError('Qualified first angular case differs')
    a.out.mkdir(parents=True,exist_ok=False);drivers=a.out/'drivers';drivers.mkdir()
    for p in Path(__file__).parent.glob('rya1354_*.py'):shutil.copyfile(p,drivers/p.name)
    receipt=dict(schema='rya1354.fac-radial-angular-remaining.v40',status='RUNNING',started_utc=datetime.now(timezone.utc).isoformat(),first_angular_receipt_sha256=sha(a.first_angular/'receipt.json'),canary_report_sha256=sha(a.canary/'report.json'),cases={},scientific_admission='HOLD',exposure='EXPOSED_NONBLIND',adoption_allowed=False,accepted_lifetime=None,uncertainty=None,covariance=None,calibrated_sigma=None,scope='Three remaining original-CM angular conversions and exact signed replay only; unchanged current child/math/controls; no Hamiltonian/orbitals/rates/lifetime.')
    try:
        for case in CASES:
            cmd=command(a,case);(a.out/(case+'.command.json')).write_text(json.dumps(cmd,indent=2)+'\n')
            with (a.out/(case+'.stdout')).open('xb') as stdout,(a.out/(case+'.stderr')).open('xb') as stderr:
                p=subprocess.Popen(cmd,stdout=stdout,stderr=stderr,start_new_session=True)
                try:p.wait(timeout=2400)
                except subprocess.TimeoutExpired:
                    stop_owned_group(p.pid)
                    try:p.wait(timeout=10)
                    except subprocess.TimeoutExpired:
                        try:os.killpg(p.pid,signal.SIGKILL)
                        except ProcessLookupError:pass
                        p.wait()
                    raise ValueError('Owned angular case timeout: '+case)
            if p.returncode:
                stop_owned_group(p.pid)
                raise ValueError('Angular case exit '+str(p.returncode)+': '+case)
            child=verify_receipt(a.out/case)
            if child['status']!='FAC_RADIAL_ANGULAR_REPLAY_COMPLETED_REVIEW_REQUIRED' or child['case']!=case or child['canary_report_sha256']!=receipt['canary_report_sha256']:raise ValueError('Angular child incomplete/unqualified: '+case)
            report=json.loads((a.out/case/'report.json').read_text())
            receipt['cases'][case]=dict(receipt_sha256=sha(a.out/case/'receipt.json'),report_sha256=sha(a.out/case/'report.json'),targets=report['targets'],full_LS_coefficients_reconstructed=report['full_LS_coefficients_reconstructed'],label_ID_overflows=report['label_ID_overflows'])
            (a.out/'progress.json').write_text(json.dumps(receipt,indent=2)+'\n')
        receipt['status']='FAC_RADIAL_ANGULAR_REMAINING_COMPLETED_REVIEW_REQUIRED'
    except Exception as exc:
        receipt.update(status='FAILED_RETAINED',error=type(exc).__name__+': '+str(exc));raise
    finally:
        receipt['files']={str(p.relative_to(a.out)):dict(sha256=sha(p),bytes=p.stat().st_size) for p in sorted(a.out.rglob('*')) if p.is_file()}
        (a.out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')

if __name__=='__main__':main()
