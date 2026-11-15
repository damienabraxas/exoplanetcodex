"""Paired 4s/4p CI protocol; only orbital-control scale varies, v35 identity unchanged."""
import argparse
from collections import Counter
from datetime import datetime,timezone
import json
from pathlib import Path
import platform
import shutil
import subprocess
from rya1354_fac_bootstrap import PIN,sha,level_inventory
from rya1354_fac_radial_preflight import input_text as radial_input,OLD,NEW,name,wave,review as radial_review
from rya1354_fac_stream_transport_v35 import EXTENDED,configurations,inventory,transport
from rya1354_fac_potential_fields import layout,compare
from rya1354_fac_no_correction_audit import verify_config_inventory
from rya1354_fac_rates import verify_receipt


def input_text(scale=10):
    text=radial_input(scale)+"Closed('1s')\n";groups=[]
    for i,q in enumerate(configurations(EXTENDED)):
        group='ext'+str(i);groups.append(group)
        cfg='2s2 '+' '.join(shell+str(n) for (shell,_,_),n in zip(EXTENDED,q) if n)
        text+=f"Config('{group}', '{cfg}')\n"
    text+="ConfigEnergy(0)\nConfigDiagnostics('configs.before.tsv')\n"
    text+=f"Structure('o.lev.b', {groups!r})\nMemENTable('o.lev.b')\nPrintTable('o.lev.b','o.lev',1)\nBasisTable('o.basis',0)\nBasisTable('o_export',10)\nConfigDiagnostics('configs.after.tsv')\n"
    for n,k in OLD+NEW:text+=f"WaveFuncTable('ci_{name(n,k)}.wav',{n},{k})\n"
    return text+"SavePotential('ci.pot')\n"


def review(out,spec,scale=10):
    text=(out/'input.sf').read_text()
    if text!=input_text(scale):raise ValueError('Declared radial CI input differs')
    files={p.name:p.read_bytes() for p in out.iterdir() if p.is_file() and (p.suffix in ('.pot','.wav','.overlap') or p.name=='input.sf')}
    files['input.sf']=radial_input(scale).encode()
    preflight=radial_review(files,scale,spec)
    for n,k in OLD+NEW:
        prefix='before_' if (n,k) in OLD else 'new_'
        if wave(files[prefix+name(n,k)+'.wav'].decode(),n,k)!=wave(files['ci_'+name(n,k)+'.wav'].decode(),n,k):raise ValueError('CI changed radial orbital')
    potential=compare(files['source.pot'],files['ci.pot'],spec)
    rows=level_inventory((out/'o.lev').read_text());expected=inventory(EXTENDED,configurations(EXTENDED))
    if Counter((r['parity'],r['two_J']) for r in rows)!=expected or set(r['index'] for r in rows)!=set(range(6600)):raise ValueError('Extended CI level identity/J/parity inventory differs')
    # Header checks are independent of the coefficient export's electron count.
    import re
    lt=(out/'o.lev').read_text()
    if not re.search(r'^O Z\s*=\s*8\.0\s*$',lt,re.M) or set(re.findall(r'^NELE\s*=\s*(\d+)',lt,re.M))!={'8'}:raise ValueError('Native element/electron inventory differs')
    keys=set()
    with (out/'o.basis').open() as src:
        for line in src:
            if 'Mixing Coefficients' in line:break
            if line.strip() and not line.startswith('#') and 'Basis Table' not in line:
                f=line.split()
                if len(f)!=10:raise ValueError('Basis layout differs')
                keys.add((int(f[4]),int(f[5])))
    configs={v:verify_config_inventory((out/('configs.'+v+'.tsv')).read_text(),keys) for v in ('before','after')}
    if any(v['groups']!=114 for v in configs.values()):raise ValueError('Extended configuration group count differs')
    tr=transport(out/'o.basis',out/'o_export_000.c',out/'o_export_000.cm',EXTENDED)
    return dict(radial_preflight=preflight,post_CI_orbitals_unchanged=True,potential=potential,configuration_state=configs,transport=tr)


CASES=(('ground_coherent100',100),('equal_coherent10',10),('equal_coherent100',100))


def main():
    ap=argparse.ArgumentParser()
    for n in ('source','baseline','preflight','first','identity','out'):ap.add_argument('--'+n,type=Path,required=True)
    a=ap.parse_args()
    if platform.node().lower().split('.')[0]!='sirius':raise RuntimeError('CI runs belong on Sirius')
    for k,v in vars(a).items():setattr(a,k,v.resolve())
    baseline=verify_receipt(a.baseline);preflight=verify_receipt(a.preflight)
    first=verify_receipt(a.first);identity=verify_receipt(a.identity);probe=baseline['probe'];binary=Path(probe['path'])
    if (baseline['status']!='FAC_FULL_VALENCE_LEVELS_COMPLETED_REVIEW_REQUIRED' or preflight['status']!='FAC_RADIAL_PREFLIGHT_COMPLETED_REVIEW_REQUIRED'
        or identity['status']!='FAC_RADIAL_IDENTITY_REPLAY_COMPLETED_REVIEW_REQUIRED' or identity['source_native_receipt_sha256']!=sha(a.first/'receipt.json')
        or first['status']!='FAILED_RETAINED' or first['probe']!=probe or first['baseline_receipt_sha256']!=sha(a.baseline/'receipt.json')
        or preflight['baseline_receipt_sha256']!=sha(a.baseline/'receipt.json') or baseline['fac_pin']!=PIN
        or subprocess.check_output(['git','rev-parse','HEAD'],cwd=a.source,text=True).strip()!=PIN
        or subprocess.check_output(['git','diff','--name-only'],cwd=a.source,text=True).strip()
        or sha(binary)!=probe['sha256'] or sha(a.source/'libfac.a')!=probe['libfac_sha256'] or sha(a.source/'sfac/sfac')!=probe['upstream_binary_sha256']):raise ValueError('Paired CI provenance differs')
    if (a.first/'ground_coherent10/input.sf').read_text()!=input_text(10):raise ValueError('First-case protocol changed')
    for case,scale in CASES:
        if (a.preflight/case/'source.pot').read_bytes()!=(a.baseline/case/'o.pot').read_bytes():raise ValueError('Paired preflight potential differs')
    a.out.mkdir(parents=True,exist_ok=False);drivers=a.out/'drivers';drivers.mkdir()
    for p in Path(__file__).parent.glob('rya1354_*.py'):shutil.copyfile(p,drivers/p.name)
    for n in ('faclib/structure.c','faclib/radial.c','faclib/orbital.c','faclib/config.c','faclib/consts.h','sfac/sfac.c','doc/manual.tex'):
        p=a.out/'source_documents'/n;p.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(a.source/n,p)
    spec=layout((a.source/'faclib/radial.c').read_text(),(a.source/'faclib/orbital.c').read_text(),(a.source/'faclib/consts.h').read_text())
    receipt=dict(schema='rya1354.fac-radial-ci-remaining.v36',status='RUNNING',started_utc=datetime.now(timezone.utc).isoformat(),fac_pin=PIN,probe=probe,baseline_receipt_sha256=sha(a.baseline/'receipt.json'),preflight_receipt_sha256=sha(a.preflight/'receipt.json'),first_native_receipt_sha256=sha(a.first/'receipt.json'),first_identity_receipt_sha256=sha(a.identity/'receipt.json'),potential_layout=spec,cases={},scientific_admission='HOLD',exposure='EXPOSED_NONBLIND',adoption_allowed=False,accepted_lifetime=None,uncertainty=None,covariance=None,calibrated_sigma=None,scope='Remaining three paired 4s/4p fixed-potential CI cases;114NRgroups/6600levels each;no freshSCF,LS,rates,lifetime or full-n4 convergence claim.')
    try:
        for case,scale in CASES:
            out=a.out/case;out.mkdir();shutil.copyfile(a.baseline/case/'o.pot',out/'source.pot');(out/'input.sf').write_text(input_text(scale))
            with (out/'sfac.stdout').open('xb') as stdout,(out/'sfac.stderr').open('xb') as stderr:
                p=subprocess.run(['/usr/bin/time','-v','-o','native.time',str(binary),'input.sf'],cwd=out,stdout=stdout,stderr=stderr,timeout=1800,check=False)
            if p.returncode:raise ValueError('Extended native exit '+str(p.returncode)+': '+case)
            receipt['cases'][case]=review(out,spec,scale)
            (a.out/'progress.json').write_text(json.dumps(receipt,indent=2)+'\n')
        receipt['status']='FAC_RADIAL_CI_REMAINING_COMPLETED_REVIEW_REQUIRED'
    except Exception as exc:
        receipt.update(status='FAILED_RETAINED',error=type(exc).__name__+': '+str(exc));raise
    finally:
        receipt['files']={str(p.relative_to(a.out)):dict(sha256=sha(p),bytes=p.stat().st_size) for p in sorted(a.out.rglob('*')) if p.is_file()}
        (a.out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')

if __name__=='__main__':main()
