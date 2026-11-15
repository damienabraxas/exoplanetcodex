"""Full fixed-orbital valence CI rates from hash-qualified full LS mapping."""
import argparse
from datetime import datetime,timezone
import json
from pathlib import Path
import platform
import shutil
import subprocess
from rya1354_fac_full_space import input_text,inspect
from rya1354_fac_space import input_text as sd_ground_input
from rya1354_fac_expanded_rates import rate_input as sd_rate_input
from rya1354_fac_no_correction_audit import configuration_keys,verify_config_inventory
from rya1354_fac_controls import diagnostics
from rya1354_fac_rates import verify_receipt
from rya1354_fac_bootstrap import PIN,sha
from rya1354_fac_recouple import GRASP_PIN
from rya1354_fac_expanded_rate_audit import review as rate_review
from rya1354_fac_rate_audit import constants

AUDIT_SHA='741f0199af1dc6793b866d38c1fef38e058c4f2d14c01b6ee141c4a345ed4be6'


def rate_input(policy,mult,scale,gauge,mapping):
    text,requests,lower=sd_rate_input(2,mult,scale,gauge,mapping);old=sd_ground_input(2,mult,scale)
    if not text.startswith(old):raise ValueError('Full rate operator prefix differs')
    return input_text(policy,mult,scale)+text[len(old):],requests,lower


def main():
    ap=argparse.ArgumentParser()
    for n in ('source','states','mapped','baseline','audit','out'):ap.add_argument('--'+n,type=Path,required=True)
    a=ap.parse_args()
    if platform.node().lower().split('.')[0]!='sirius':raise RuntimeError('Full valence rate calculations belong on Sirius')
    for n,v in vars(a).items():setattr(a,n,v.resolve())
    sr=verify_receipt(a.states);cr=verify_receipt(a.mapped);br=verify_receipt(a.baseline)
    if sha(a.audit)!=AUDIT_SHA:raise ValueError('Declared full LS audit hash differs')
    mapping=json.loads(a.audit.read_text());probe=sr['probe'];binary=Path(probe['path'])
    if (sr['status']!='FAC_FULL_VALENCE_LEVELS_COMPLETED_REVIEW_REQUIRED' or cr['status']!='FAC_ANGULAR_TRANSFORM_COMPLETED_REVIEW_REQUIRED'
        or cr['source_receipt_sha256']!=sha(a.states/'receipt.json') or mapping['source_receipt_sha256']!=sha(a.states/'receipt.json') or mapping['converter_receipt_sha256']!=sha(a.mapped/'receipt.json')
        or sr['baseline_receipt_sha256']!=sha(a.baseline/'receipt.json') or sr['probe']!=br['probe'] or cr['grasp_pin']!=GRASP_PIN
        or sr['fac_pin']!=PIN or subprocess.check_output(['git','rev-parse','HEAD'],cwd=a.source,text=True).strip()!=PIN
        or subprocess.check_output(['git','diff','--name-only'],cwd=a.source,text=True).strip()
        or sha(binary)!=probe['sha256'] or sha(a.source/'libfac.a')!=probe['libfac_sha256'] or sha(a.source/'sfac/sfac')!=probe['upstream_binary_sha256']):raise ValueError('Full rate source/mapping/frontend qualification differs')
    a.out.mkdir(parents=True,exist_ok=False);drivers=a.out/'drivers';drivers.mkdir()
    for p in Path(__file__).parent.glob('rya1354_*.py'):shutil.copyfile(p,drivers/p.name)
    receipt=dict(schema='rya1354.fac-full-valence-rates.v30',status='RUNNING',started_utc=datetime.now(timezone.utc).isoformat(),
        source_receipt_sha256=sha(a.states/'receipt.json'),mapping_receipt_sha256=sha(a.mapped/'receipt.json'),audit_sha256=AUDIT_SHA,probe=probe,
        cases={},scientific_admission='HOLD',exposure='EXPOSED_NONBLIND',adoption_allowed=False,accepted_lifetime=None,uncertainty=None,covariance=None,calibrated_sigma=None)
    try:
        for policy in ('ground','equal'):
            for label,mult,scale in [('coherent10',.3,10),('coherent100',.03,100)]:
                name=policy+'_'+label;s=a.states/name;m=mapping['cases'][name]
                if (s/'input.sf').read_text()!=input_text(policy,mult,scale) or not diagnostics((s/'sfac.stdout').read_text(),mult,256)['printed_rule_verified']:raise ValueError('Full rate source input/stopping differs')
                inspect((s/'o.lev').read_text());keys=configuration_keys((s/'o.basis').read_text())
                for n in ('before','after'):verify_config_inventory((s/('configs.'+n+'.tsv')).read_text(),keys)
                for suffix,gauge in [('B',2),('C',1)]:
                    case=name+'_'+suffix;out=a.out/case;out.mkdir();text,requests,lower=rate_input(policy,mult,scale,gauge,m);(out/'input.sf').write_text(text)
                    with (out/'sfac.stdout').open('xb') as stdout,(out/'sfac.stderr').open('xb') as stderr:
                        p=subprocess.run([str(binary),'input.sf'],cwd=out,stdout=stdout,stderr=stderr,timeout=900,check=False)
                    if p.returncode:raise ValueError(f'Full rate native exit{p.returncode}: '+case)
                    d=diagnostics((out/'sfac.stdout').read_text(),mult,256)
                    if not d['printed_rule_verified']:raise ValueError('Full rate stopping failure')
                    for fn in ('o_export_000.c','o_export_000.cm','o.pot'):
                        if (out/fn).read_bytes()!=(s/fn).read_bytes():raise ValueError('Full rate regenerated mapped state/potential differs')
                    configs={n:verify_config_inventory((out/('configs.'+n+'.tsv')).read_text(),keys) for n in ('before','after')}
                    audit=rate_review((out/'raw.tsv').read_text(),(out/'rates.tr').read_text(),requests,gauge,m,constants((a.source/'faclib/consts.h').read_text()))
                    receipt['cases'][case]=dict(diagnostics=d,requests=requests,lower_states=lower,gauge=gauge,configuration_state=configs,audit=audit)
                    (a.out/'progress.json').write_text(json.dumps(receipt,indent=2)+'\n')
        receipt['status']='FAC_FULL_VALENCE_RATES_COMPLETED_REVIEW_REQUIRED'
    except Exception as exc:
        receipt.update(status='FAILED_RETAINED',error=type(exc).__name__+': '+str(exc));raise
    finally:
        receipt['files']={str(p.relative_to(a.out)):dict(sha256=sha(p),bytes=p.stat().st_size) for p in sorted(a.out.rglob('*')) if p.is_file()}
        (a.out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')

if __name__=='__main__':main()
