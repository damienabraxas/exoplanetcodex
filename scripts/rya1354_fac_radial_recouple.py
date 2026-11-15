"""Single selected 6600-level native LS conversion and fully qualified replay."""
import argparse
from datetime import datetime,timezone
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
from rya1354_fac_bootstrap import sha
from rya1354_fac_rates import verify_receipt
from rya1354_fac_full_recouple import CONTROLS,GRASP_PIN
from rya1354_fac_radial_ci_v36 import review as ci_review
from rya1354_fac_radial_ls_audit import review,EXTENDED

CASES=('ground_coherent10','ground_coherent100','equal_coherent10','equal_coherent100')


def qualified_canary(report):
    location=report.get('vector_control_location',{})
    if location.get('field')!='actual_first_vector_coefficient' or float(location.get('original',0))!=-float(location.get('mutated',0)) or abs(float(location.get('original',0)))<=.01:raise ValueError('Fast replay canary unqualified: actual coefficient location missing')
    if (report.get('schema')!='rya1354.fac-ls-fast-canary.v39' or report['status']!='EXACT_LEGACY_EQUIVALENCE_AND_NEGATIVE_CONTROLS_PASS' or set(report['cases'])!=set(CASES)
        or any(not x['all_nondiagnostic_report_fields_exact'] or x['full_LS_coefficients']!=329198 for x in report['cases'].values())
        or report['negative_controls']!={'d4_seniority':'Unsupported NK d-shell seniority','vector_sign':'Expanded LS coefficient/sign replay differs','matrix_sign':'Expanded LS coefficient/sign replay differs','explicit_stdout_ID':'Explicit stdout/binary ID order differs'} or not report['forced_proposal_recovery']['pass_exact']):raise ValueError('Fast replay canary unqualified')


def main():
    ap=argparse.ArgumentParser()
    for n in ('grasp','first','identity','remaining','old_converted','canary','out'):ap.add_argument('--'+n,type=Path,required=True)
    ap.add_argument('--case',choices=CASES,required=True);a=ap.parse_args()
    if platform.node().lower().split('.')[0]!='sirius':raise RuntimeError('Angular transforms and replay belong on Sirius')
    import numpy as np
    for k,v in vars(a).items():
        if isinstance(v,Path):setattr(a,k,v.resolve())
    first=verify_receipt(a.first);identity=verify_receipt(a.identity);remaining=verify_receipt(a.remaining);old=verify_receipt(a.old_converted)
    canary=json.loads((a.canary/'report.json').read_text());qualified_canary(canary)
    binary=a.grasp/'build/bin/jj2lsj';isodata=a.old_converted/'ground_coherent10/isodata'
    if (identity['status']!='FAC_RADIAL_IDENTITY_REPLAY_COMPLETED_REVIEW_REQUIRED' or first['status']!='FAILED_RETAINED'
        or identity['source_native_receipt_sha256']!=sha(a.first/'receipt.json') or remaining['status']!='FAC_RADIAL_CI_REMAINING_COMPLETED_REVIEW_REQUIRED'
        or remaining['first_native_receipt_sha256']!=sha(a.first/'receipt.json') or remaining['first_identity_receipt_sha256']!=sha(a.identity/'receipt.json')
        or old['status']!='FAC_ANGULAR_TRANSFORM_COMPLETED_REVIEW_REQUIRED' or old['grasp_pin']!=GRASP_PIN
        or {x['sha256'] for x in old['steps']}!={sha(binary)}
        or subprocess.check_output(['git','rev-parse','HEAD'],cwd=a.grasp,text=True).strip()!=GRASP_PIN
        or subprocess.check_output(['git','diff','--name-only'],cwd=a.grasp,text=True).strip()
        or canary['converter_receipt_sha256']!=sha(a.old_converted/'receipt.json')
        or canary['python_version']!=platform.python_version() or canary['numpy_version']!=np.__version__
        or sha(a.canary/'rya1354_fac_radial_ls_audit.py')!=sha(Path(__file__).with_name('rya1354_fac_radial_ls_audit.py'))
        or any(os.environ.get(k)!='1' for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS'))):raise ValueError('Angular source/runtime/canary provenance differs')
    first_case=a.case=='ground_coherent10';source=(a.first if first_case else a.remaining)/a.case
    source_receipt=first if first_case else remaining;scale=100 if a.case.endswith('coherent100') else 10
    expected=identity['cases'][a.case] if first_case else remaining['cases'][a.case]
    if ci_review(source,source_receipt['potential_layout'],scale)!=expected:raise ValueError('Qualified original CI replay differs')
    a.out.mkdir(parents=True,exist_ok=False);drivers=a.out/'drivers';drivers.mkdir()
    for p in Path(__file__).parent.glob('rya1354_*.py'):shutil.copyfile(p,drivers/p.name)
    for n in ('getmixblock.f90','jj2lsj_code.f90','jj2lsj2K.f90','packLS.f90'):
        p=a.out/'converter_source'/n;p.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(a.grasp/'src/appl/jj2lsj90'/n,p)
    out=a.out/a.case;out.mkdir()
    for src,dest in [('o_export_000.c','state.c'),('o_export_000.cm','state.cm')]:shutil.copyfile(source/src,out/dest)
    shutil.copyfile(isodata,out/'isodata');(out/'jj2lsj.stdin').write_text(CONTROLS)
    receipt=dict(schema='rya1354.fac-radial-recoupling.v38',status='RUNNING',started_utc=datetime.now(timezone.utc).isoformat(),case=a.case,grasp_pin=GRASP_PIN,binary_sha256=sha(binary),isodata_sha256=sha(isodata),first_receipt_sha256=sha(a.first/'receipt.json'),identity_receipt_sha256=sha(a.identity/'receipt.json'),remaining_receipt_sha256=sha(a.remaining/'receipt.json'),old_converter_receipt_sha256=sha(a.old_converted/'receipt.json'),canary_report_sha256=sha(a.canary/'report.json'),python_version=platform.python_version(),numpy_version=np.__version__,qualified_CI=expected,scientific_admission='HOLD',exposure='EXPOSED_NONBLIND',adoption_allowed=False,accepted_lifetime=None,uncertainty=None,covariance=None,calibrated_sigma=None,scope='Original FAC C/CM exact angular representation conversion only; full unmodified native output and complete signed replay. Z8 isodata converter support only, not a change/match of FAC nucleus. No new Hamiltonian/orbitals/rates/lifetime.')
    try:
        with (out/'jj2lsj.stdout').open('xb') as stdout,(out/'jj2lsj.stderr').open('xb') as stderr:
            result=subprocess.run(['/usr/bin/time','-v','-o','native.time',str(binary)],input=CONTROLS.encode(),cwd=out,stdout=stdout,stderr=stderr,timeout=1800,check=False)
        if result.returncode:raise ValueError('Native LS conversion exit '+str(result.returncode))
        output=((out/'jj2lsj.stdout').read_text(errors='replace')+(out/'jj2lsj.stderr').read_text(errors='replace')).upper()
        if any(x in output for x in ('TOO MANY ELECTRONS FOR THIS CODE','CONFIGURATION DATA WRONG','PROGRAM HALTED')) or ('THE ORDER OF' in output and 'IS INCONSISTENT' in output):raise ValueError('Native LS semantic refusal')
        for src,dest in [('o_export_000.c','state.c'),('o_export_000.cm','state.cm')]:
            if sha(source/src)!=sha(out/dest):raise ValueError('Converter changed original FAC C/CM')
        s={'o_export_000.cm':(source/'o_export_000.cm').read_bytes()}
        c={n:(out/n).read_bytes() for n in ('state.lsj.c','state.lsj.j','state.lsj.lbl','state.lsj.T','jj2lsj.stdout')}
        report=review(s,c,EXTENDED);(a.out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
        receipt['report_sha256']=sha(a.out/'report.json');receipt['status']='FAC_RADIAL_ANGULAR_REPLAY_COMPLETED_REVIEW_REQUIRED'
    except Exception as exc:
        receipt.update(status='FAILED_RETAINED',error=type(exc).__name__+': '+str(exc));raise
    finally:
        receipt['files']={str(p.relative_to(a.out)):dict(sha256=sha(p),bytes=p.stat().st_size) for p in sorted(a.out.rglob('*')) if p.is_file()}
        (a.out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')

if __name__=='__main__':main()
