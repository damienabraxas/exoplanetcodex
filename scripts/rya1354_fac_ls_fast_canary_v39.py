"""Sirius-only exact old-report equivalence and full signed negative controls."""
import argparse
from datetime import datetime,timezone
import json
import os
from pathlib import Path
import platform
import re
import shutil
import struct
from rya1354_fac_bootstrap import sha
from rya1354_fac_rates import verify_receipt
from rya1354_fac_radial_ls_audit import review
from rya1354_fac_recoupling_audit import records

EXTRA={'max_transformation_orthonormality_residual','signed_rounding_fallbacks','weight_rounding_fallbacks','orthonormality_arithmetic'}


def inputs(source,converted,case):
    p=source/case;c=converted/case
    return ({'o_export_000.cm':(p/'o_export_000.cm').read_bytes()},
            {name:(c/name).read_bytes() for name in ('state.lsj.c','state.lsj.j','state.lsj.lbl','state.lsj.T','jj2lsj.stdout')})


def flip_actual_vector_coefficient(text):
    lines=text.splitlines(keepends=True)
    start=next(i for i,line in enumerate(lines) if 'Ssms=' in line)+2
    for i in range(start,len(lines)):
        if 'Ssms=' in lines[i] or '2*J =' in lines[i]:break
        if not re.fullmatch(r'\s*(?:[-+]?\d+\.\d{8}\s*)+',lines[i]):continue
        for match in re.finditer(r'[-+]?\d+\.\d{8}',lines[i]):
            if abs(float(match[0]))>.01:
                before=match[0];after=before[1:] if before.startswith('-') else '-'+before
                lines[i]=lines[i][:match.start()]+after+lines[i][match.end():]
                return ''.join(lines),dict(line=i,column=match.start(),original=before,mutated=after,field='actual_first_vector_coefficient')
    raise ValueError('No nonzero coefficient inside actual first vector')


def main():
    ap=argparse.ArgumentParser()
    for name in ('source','converted','reference','previous','out'):ap.add_argument('--'+name,type=Path,required=True)
    a=ap.parse_args()
    if platform.node().lower().split('.')[0]!='sirius':raise RuntimeError('Angular canaries belong on Sirius')
    import numpy as np
    sr=verify_receipt(a.source);cr=verify_receipt(a.converted)
    if sr['status']!='FAC_FULL_VALENCE_LEVELS_COMPLETED_REVIEW_REQUIRED' or cr['status']!='FAC_ANGULAR_TRANSFORM_COMPLETED_REVIEW_REQUIRED' or cr['source_receipt_sha256']!=sha(a.source/'receipt.json'):raise ValueError('Old complete source/converter provenance differs')
    if sha(a.reference)!='741f0199af1dc6793b866d38c1fef38e058c4f2d14c01b6ee141c4a345ed4be6':raise ValueError('Qualified full signed LS reference hash differs')
    previous=json.loads((a.previous/'report.json').read_text())
    if previous['schema']!='rya1354.fac-ls-fast-canary.v37' or previous['negative_controls']['vector_sign']!='Missing LS vector identity':raise ValueError('Prior control disposition differs')
    baseline=json.loads(a.reference.read_text());a.out.mkdir(parents=True,exist_ok=False)
    for p in Path(__file__).parent.glob('rya1354_*.py'):shutil.copyfile(p,a.out/p.name)
    result=dict(prior_canary_report_sha256=sha(a.previous/'report.json'),python_version=platform.python_version(),numpy_version=np.__version__,thread_environment={k:os.environ.get(k) for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS')},schema='rya1354.fac-ls-fast-canary.v39',prior_control_disposition='v37 vector-sign control modified an energy identity; prior report retained, coefficient-sign coverage qualified only by v39',status='RUNNING',started_utc=datetime.now(timezone.utc).isoformat(),reference_sha256=sha(a.reference),source_receipt_sha256=sha(a.source/'receipt.json'),converter_receipt_sha256=sha(a.converted/'receipt.json'),cases={},negative_controls={},scientific_admission='HOLD',exposure='EXPOSED_NONBLIND',adoption_allowed=False,accepted_lifetime=None,uncertainty=None,covariance=None,calibrated_sigma=None)
    try:
        for case in baseline['cases']:
            s,c=inputs(a.source,a.converted,case);r=review(s,c);old=baseline['cases'][case]
            if {k:v for k,v in r.items() if k not in EXTRA}!={k:v for k,v in old.items() if k not in EXTRA}:raise ValueError('Exact legacy scientific/index/signed replay differs: '+case)
            result['cases'][case]=dict(full_LS_coefficients=r['full_LS_coefficients_reconstructed'],levels=len(r['levels']),targets=r['targets'],old_orth_diagnostic=old['max_transformation_orthonormality_residual'],new_orth_diagnostic=r['max_transformation_orthonormality_residual'],signed_rounding_fallbacks=r['signed_rounding_fallbacks'],weight_rounding_fallbacks=r['weight_rounding_fallbacks'],all_nondiagnostic_report_fields_exact=True)
        s,c=inputs(a.source,a.converted,'ground_coherent100')
        import rya1354_fac_radial_ls_audit as audit
        original=audit.product
        def perturbed(values,matrix):
            proposed=original(values,matrix);proposed[0,0]+=.01;return proposed
        audit.product=perturbed
        try:recovered=review(s,c)
        finally:audit.product=original
        old=baseline['cases']['ground_coherent100']
        if recovered['signed_rounding_fallbacks']<1 or {k:v for k,v in recovered.items() if k not in EXTRA}!={k:v for k,v in old.items() if k not in EXTRA}:raise ValueError('Strict legacy fallback recovery control failed')
        result['forced_proposal_recovery']=dict(pass_exact=True,fallbacks=recovered['signed_rounding_fallbacks'])
        for mode in ('d4_seniority','vector_sign','matrix_sign','explicit_stdout_ID'):
            changed=dict(c)
            if mode=='d4_seniority':
                lines=c['state.lsj.c'].decode().splitlines();done=False
                for i,line in enumerate(lines[:-1]):
                    if re.search(r'3d\(\s*4\)',line) and '3P2' in lines[i+1]:lines[i+1]=lines[i+1].replace('3P2','3P9',1);done=True;break
                if not done:raise ValueError('No d4 seniority control')
                changed['state.lsj.c']='\n'.join(lines).encode()
            elif mode=='vector_sign':
                mutated,location=flip_actual_vector_coefficient(c['state.lsj.j'].decode())
                changed['state.lsj.j']=mutated.encode();result['vector_control_location']=location
            elif mode=='matrix_sign':
                rr=records(c['state.lsj.T']);_,_,_,nb=struct.unpack('<4i',rr[1]);done=False
                for i in range(3+nb,len(rr)):
                    if len(rr[i])==16:
                        jj,l,v=struct.unpack('<iid',rr[i])
                        if abs(v)>.01:rr[i]=struct.pack('<iid',jj,l,-v);done=True;break
                if not done:raise ValueError('No matrix-sign control')
                changed['state.lsj.T']=b''.join(struct.pack('<i',len(r))+r+struct.pack('<i',len(r)) for r in rr)
            else:
                lines=c['jj2lsj.stdout'].decode().splitlines();waiting=False;done=False
                for i,line in enumerate(lines):
                    if 'Weights of major contributors to ASF in LS-coupling:' in line:waiting=True;continue
                    if waiting and re.match(r'^\s*\d+\s+\d+\s+[+-]\s+.*\bof\b',line):lines[i]=re.sub(r'^(\s*)\d+',r'\g<1>9999',line,count=1);done=True;break
                if not done:raise ValueError('No stdout identity control')
                changed['jj2lsj.stdout']='\n'.join(lines).encode()
            try:review(s,changed)
            except ValueError as exc:
                expected={'d4_seniority':'Unsupported NK d-shell seniority','vector_sign':'Expanded LS coefficient/sign replay differs','matrix_sign':'Expanded LS coefficient/sign replay differs','explicit_stdout_ID':'Explicit stdout/binary ID order differs'}[mode]
                if str(exc)!=expected:raise ValueError('Wrong rejection mechanism for '+mode+': '+str(exc))
                result['negative_controls'][mode]=str(exc)
            else:raise ValueError('Corruption admitted: '+mode)
        result['status']='EXACT_LEGACY_EQUIVALENCE_AND_NEGATIVE_CONTROLS_PASS'
    except Exception as exc:
        result.update(status='FAILED_RETAINED',error=type(exc).__name__+': '+str(exc));raise
    finally:
        (a.out/'report.json').write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':main()
