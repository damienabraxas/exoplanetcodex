"""Hash-qualified expanded radial states and complete both-gauge decay rates."""
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

from rya1354_fac_bootstrap import PIN,sha
from rya1354_fac_rates import verify_receipt
from rya1354_fac_radial_ci_v36 import input_text as ci_input,review as ci_review
from rya1354_fac_radial_recouple import CASES,qualified_canary
from rya1354_fac_radial_angular_remaining_metadata import native_metadata,hold
from rya1354_fac_expanded_rates import rate_input as operator_input
from rya1354_fac_space import input_text as old_input
from rya1354_fac_expanded_rate_audit import review as rate_review
from rya1354_fac_rate_audit import constants


def rate_input(scale,gauge,mapping):
    """Keep the exact qualified radial CI protocol; append existing operator calls."""
    if scale not in (10,100) or gauge not in (1,2):
        raise ValueError('Undeclared radial rate controls')
    text,requests,lower=operator_input(2,3/scale,scale,gauge,mapping)
    prefix=old_input(2,3/scale,scale)
    if not text.startswith(prefix):raise ValueError('Operator input prefix differs')
    return ci_input(scale)+text[len(prefix):],requests,lower


def source_path(root,case):
    parent='fac_radial_ci_first_v33' if case=='ground_coherent10' else 'fac_radial_ci_remaining_v36'
    return root/'runs'/parent/case


def prepare(root):
    """Verify sealed original receipts and full native metadata, without new math."""
    runs=root/'runs'
    parents={key:verify_receipt(runs/name) for key,name in (
        ('first','fac_radial_ci_first_v33'),('identity','fac_radial_identity_replay_v35'),
        ('remaining','fac_radial_ci_remaining_v36'),('angular_first','fac_radial_angular_first_v38'),
        ('angular_remaining','fac_radial_angular_remaining_v40'))}
    first,identity,remaining=(parents[k] for k in ('first','identity','remaining'))
    if (first['status']!='FAILED_RETAINED' or identity['status']!='FAC_RADIAL_IDENTITY_REPLAY_COMPLETED_REVIEW_REQUIRED'
        or remaining['status']!='FAC_RADIAL_CI_REMAINING_COMPLETED_REVIEW_REQUIRED'
        or identity['source_native_receipt_sha256']!=sha(runs/'fac_radial_ci_first_v33/receipt.json')
        or remaining['first_native_receipt_sha256']!=sha(runs/'fac_radial_ci_first_v33/receipt.json')
        or remaining['first_identity_receipt_sha256']!=sha(runs/'fac_radial_identity_replay_v35/receipt.json')
        or first['probe']!=remaining['probe'] or first['potential_layout']!=remaining['potential_layout']):
        raise ValueError('Qualified radial CI chain differs')
    canary_path=runs/'fac_ls_fast_canary_v39/report.json'
    qualified_canary(json.loads(canary_path.read_text()))
    batch=parents['angular_remaining'];hold(batch)
    if (batch['status']!='FAC_RADIAL_ANGULAR_REMAINING_COMPLETED_REVIEW_REQUIRED'
        or set(batch['cases'])!=set(CASES)-{'ground_coherent10'}
        or batch['first_angular_receipt_sha256']!=sha(runs/'fac_radial_angular_first_v38/receipt.json')
        or batch['canary_report_sha256']!=sha(canary_path)):
        raise ValueError('Qualified angular batch differs')
    result=dict(schema='rya1354.fac-radial-rate-index.v41',cases={},probe=first['probe'],fac_pin=first['fac_pin'],
        potential_layout=first['potential_layout'],parent_receipts={},canary_report_sha256=sha(canary_path),
        scientific_admission='HOLD',exposure='EXPOSED_NONBLIND',adoption_allowed=False,
        accepted_lifetime=None,uncertainty=None,covariance=None,calibrated_sigma=None)
    for key,name in [('first','fac_radial_ci_first_v33'),('identity','fac_radial_identity_replay_v35'),
                     ('remaining','fac_radial_ci_remaining_v36'),('angular_first','fac_radial_angular_first_v38'),
                     ('angular_remaining','fac_radial_angular_remaining_v40')]:
        result['parent_receipts'][key]=dict(path='runs/'+name+'/receipt.json',sha256=sha(runs/name/'receipt.json'))
    for case in CASES:
        first_case=case=='ground_coherent10';scale=100 if case.endswith('100') else 10
        source=source_path(root,case)
        angular=runs/'fac_radial_angular_first_v38' if first_case else runs/'fac_radial_angular_remaining_v40'/case
        receipt=parents['angular_first'] if first_case else verify_receipt(angular)
        hold(receipt)
        report_path=angular/'report.json';report=json.loads(report_path.read_text())
        qualified=identity['cases'][case] if first_case else remaining['cases'][case]
        if (receipt['status']!='FAC_RADIAL_ANGULAR_REPLAY_COMPLETED_REVIEW_REQUIRED' or receipt['case']!=case
            or receipt['report_sha256']!=sha(report_path) or receipt['qualified_CI']!=qualified
            or receipt['first_receipt_sha256']!=result['parent_receipts']['first']['sha256']
            or receipt['identity_receipt_sha256']!=result['parent_receipts']['identity']['sha256']
            or receipt['remaining_receipt_sha256']!=result['parent_receipts']['remaining']['sha256']):
            raise ValueError('Case-specific angular qualification differs: '+case)
        if sha(angular/'drivers/rya1354_fac_radial_ls_audit.py')!=sha(runs/'fac_ls_fast_canary_v39/rya1354_fac_radial_ls_audit.py'):
            raise ValueError('Qualified numerical LS module differs')
        if not first_case and receipt['canary_report_sha256']!=sha(canary_path):
            raise ValueError('Corrected angular canary hash differs')
        if not first_case and (batch['cases'][case]['receipt_sha256']!=sha(angular/'receipt.json')
                              or batch['cases'][case]['report_sha256']!=sha(report_path)):
            raise ValueError('Batch child hash differs')
        if ci_review(source,result['potential_layout'],scale)!=qualified:
            raise ValueError('Original CI replay differs: '+case)
        converted=angular/case
        if any(sha(source/src)!=sha(converted/dest) for src,dest in [('o_export_000.c','state.c'),('o_export_000.cm','state.cm')]):
            raise ValueError('Mapped original state differs')
        files={n:(converted/n).read_bytes() for n in ('state.cm','state.lsj.c','state.lsj.lbl','jj2lsj.stdout')}
        metadata=native_metadata(files,'',report)
        for term,j in [('3S',2),('5S',4)]:
            if metadata['target_metadata'][term]['two_J']!=j or metadata['target_metadata'][term]['parity']!='-':
                raise ValueError('Target J/parity differs')
        result['cases'][case]=dict(mapping=report,qualified_CI=qualified,source=str(source.relative_to(root)),
            angular=str(angular.relative_to(root)),report_sha256=sha(report_path),receipt_sha256=sha(angular/'receipt.json'))
    return result


def run_native(binary,out,timeout=1800):
    """Bound the owned native process group, retaining all files on failure."""
    with (out/'sfac.stdout').open('xb') as stdout,(out/'sfac.stderr').open('xb') as stderr:
        process=subprocess.Popen(['/usr/bin/time','-v','-o','native.time',str(binary),'input.sf'],cwd=out,
            stdout=stdout,stderr=stderr,start_new_session=True)
        try:
            process.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            os.killpg(process.pid,signal.SIGTERM)
            try:process.wait(timeout=10)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid,signal.SIGKILL);process.wait()
            raise ValueError('Owned radial rate process timed out')
    if process.returncode:raise ValueError('Native radial rate exit '+str(process.returncode))


def main():
    ap=argparse.ArgumentParser()
    for key in ('root','source','index','out'):ap.add_argument('--'+key,type=Path,required=True)
    ap.add_argument('--index-sha256',required=True)
    a=ap.parse_args()
    if platform.node().lower().split('.')[0]!='sirius':raise RuntimeError('Atomic rate work belongs on Sirius')
    a.root=a.root.resolve();a.source=a.source.resolve();a.out=a.out.resolve()
    if any(os.environ.get(k)!='1' for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS')):raise ValueError('Single-thread controls differ')
    if sha(a.index)!=a.index_sha256:raise ValueError('Declared index hash differs')
    index=json.loads(a.index.read_text());hold(index)
    if index!=prepare(a.root):raise ValueError('Frozen radial rate index differs')
    probe=index['probe'];binary=Path(probe['path'])
    if (index['fac_pin']!=PIN or subprocess.check_output(['git','rev-parse','HEAD'],cwd=a.source,text=True).strip()!=PIN
        or subprocess.check_output(['git','diff','--name-only'],cwd=a.source,text=True).strip()
        or sha(binary)!=probe['sha256'] or sha(a.source/'libfac.a')!=probe['libfac_sha256']
        or sha(a.source/'sfac/sfac')!=probe['upstream_binary_sha256']):raise ValueError('Rate solver/source/frontend provenance differs')
    a.out.mkdir(parents=True,exist_ok=False);drivers=a.out/'drivers';drivers.mkdir()
    for path in Path(__file__).parent.glob('rya1354_*.py'):shutil.copyfile(path,drivers/path.name)
    shutil.copyfile(a.index,a.out/'index.json')
    source_docs=a.out/'source_documents';source_docs.mkdir()
    for name in ('consts.h','transition.c'):
        shutil.copyfile(a.source/'faclib'/name,source_docs/name)
    c=constants((source_docs/'consts.h').read_text())
    receipt=dict(schema='rya1354.fac-radial-rates.v41',status='RUNNING',started_utc=datetime.now(timezone.utc).isoformat(),
        index_sha256=a.index_sha256,probe=probe,fac_pin=PIN,cases={},scientific_admission='HOLD',exposure='EXPOSED_NONBLIND',
        adoption_allowed=False,accepted_lifetime=None,uncertainty=None,covariance=None,calibrated_sigma=None)
    try:
        for case in CASES:
            entry=index['cases'][case];mapping=entry['mapping'];scale=100 if case.endswith('100') else 10
            source=a.root/entry['source']
            for suffix,gauge in [('B',2),('C',1)]:
                name=case+'_'+suffix;out=a.out/name;out.mkdir()
                shutil.copyfile(source/'source.pot',out/'source.pot')
                text,requests,lower=rate_input(scale,gauge,mapping);(out/'input.sf').write_text(text)
                run_native(binary,out)
                state=ci_review(out,index['potential_layout'],scale,declared_input=text)
                if state!=entry['qualified_CI']:raise ValueError('Rate CI/radial replay differs: '+name)
                for fn in ('o_export_000.c','o_export_000.cm','source.pot','ci.pot'):
                    if sha(out/fn)!=sha(source/fn):raise ValueError('Regenerated state/potential bytes differ: '+name+'/'+fn)
                audit=rate_review((out/'raw.tsv').read_text(),(out/'rates.tr').read_text(),requests,gauge,mapping,c)
                receipt['cases'][name]=dict(gauge=gauge,requests=requests,lower_states=lower,qualified_CI=state,audit=audit)
                (a.out/'progress.json').write_text(json.dumps(receipt,indent=2)+'\n')
        receipt['status']='FAC_RADIAL_RATES_COMPLETED_REVIEW_REQUIRED'
    except Exception as exc:
        receipt.update(status='FAILED_RETAINED',error=type(exc).__name__+': '+str(exc));raise
    finally:
        receipt['files']={str(p.relative_to(a.out)):dict(sha256=sha(p),bytes=p.stat().st_size)
                          for p in sorted(a.out.rglob('*')) if p.is_file()}
        (a.out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')


if __name__=='__main__':main()
