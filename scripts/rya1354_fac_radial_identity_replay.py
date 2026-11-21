"""v34 source-mapped identity replay of retained v33; no native recalculation."""
import argparse
from datetime import datetime,timezone
import json
from pathlib import Path
import platform
import shutil
import subprocess
from rya1354_fac_bootstrap import PIN,sha
from rya1354_fac_rates import verify_receipt
from rya1354_fac_potential_fields import layout
from rya1354_fac_radial_ci_v34 import review
from rya1354_fac_stream_transport_v34 import transport


def source_proof(text):
    func=text.split('int ConstructLevelName(',1)[1].split('int CompareStatePointer(',1)[0]
    for token in ('if (nc)', 'if (n == n0)', 'nq0 += nq;', '!IsClosedComplex(n0, nq0)', 'sprintf(ashell, ".%1d*%1d", n0, nq0)', '!IsClosedShellNR(n0, kl0, nq0, full_name)'):
        if token not in func:raise ValueError('Pinned principal-complex source semantics differ')
    return dict(native_principal_complex='ConstructLevelName nc sums occupations by n, omits closed complexes; sname/name may omit filled shells',identity='symmetry + principal-complex nc + NR sname + coupled jj name; exact exported-occupation reconstruction and bijection, no energy/order guess')


def main():
    ap=argparse.ArgumentParser()
    for n in ('source','native','baseline','out'):ap.add_argument('--'+n,type=Path,required=True)
    a=ap.parse_args()
    if platform.node().lower().split('.')[0]!='sirius':raise RuntimeError('Large native evidence replay belongs on Sirius')
    for k,v in vars(a).items():setattr(a,k,v.resolve())
    native=verify_receipt(a.native);parent=verify_receipt(a.baseline)
    if native['status']!='FAILED_RETAINED' or native['error']!='ValueError: Native coupled basis duplicate/identity differs' or native['baseline_receipt_sha256']!=sha(a.baseline/'receipt.json'):raise ValueError('Retained native failure/provenance differs')
    if subprocess.check_output(['git','rev-parse','HEAD'],cwd=a.source,text=True).strip()!=PIN or subprocess.check_output(['git','diff','--name-only'],cwd=a.source,text=True).strip():raise ValueError('Primary source differs')
    proof=source_proof((a.source/'faclib/structure.c').read_text())
    a.out.mkdir(parents=True,exist_ok=False);drivers=a.out/'drivers';drivers.mkdir()
    for p in Path(__file__).parent.glob('rya1354_*.py'):shutil.copyfile(p,drivers/p.name)
    receipt=dict(schema='rya1354.fac-radial-identity-replay.v34',status='RUNNING',started_utc=datetime.now(timezone.utc).isoformat(),source_native_receipt_sha256=sha(a.native/'receipt.json'),native_status_retained='FAILED_RETAINED',source_proof=proof,cases={},canaries={},scientific_admission='HOLD',exposure='EXPOSED_NONBLIND',adoption_allowed=False,accepted_lifetime=None,uncertainty=None,covariance=None,calibrated_sigma=None,native_recalculation=False)
    try:
        for name,item in parent['cases'].items():
            p=a.baseline/name;r=transport(p/'o.basis',p/'o_export_000.c',p/'o_export_000.cm')
            if r!=item['transport']:raise ValueError('Principal-complex canary differs: '+name)
            receipt['canaries'][name]=r
        spec=layout((a.source/'faclib/radial.c').read_text(),(a.source/'faclib/orbital.c').read_text(),(a.source/'faclib/consts.h').read_text())
        if spec!=native['potential_layout']:raise ValueError('Potential source layout differs')
        receipt['cases']['ground_coherent10']=review(a.native/'ground_coherent10',spec)
        if (a.native/'ground_coherent10/source.pot').read_bytes()!=(a.baseline/'ground_coherent10/o.pot').read_bytes():raise ValueError('Paired potential differs')
        receipt['status']='FAC_RADIAL_IDENTITY_REPLAY_COMPLETED_REVIEW_REQUIRED'
    except Exception as exc:
        receipt.update(status='FAILED_RETAINED',error=type(exc).__name__+': '+str(exc));raise
    finally:
        # Verify every original output remains unchanged after the replay.
        verify_receipt(a.native)
        receipt['files']={str(p.relative_to(a.out)):dict(sha256=sha(p),bytes=p.stat().st_size) for p in sorted(a.out.rglob('*')) if p.is_file()}
        (a.out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')

if __name__=='__main__':main()
