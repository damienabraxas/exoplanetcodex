"""Complete four-valence-electron CI closure in fixed n3 orbital set; no rates."""
import argparse
from collections import Counter
from datetime import datetime,timezone
from itertools import product
import ast
import json
from pathlib import Path
import platform
import re
import shutil
import subprocess
from rya1354_fac_space import SHELLS,configurations as sd_configurations,config_text,expected_inventory
from rya1354_fac_no_correction import input_text as sd_input,config_state
from rya1354_fac_rates import verify_receipt,read_files
from rya1354_fac_bootstrap import PIN,sha,level_inventory
from rya1354_fac_controls import diagnostics


def configurations():
    return [q for q in product(*(range(min(cap,4)+1) for _,_,cap in SHELLS)) if sum(q)==4]


def input_text(policy,mult,scale):
    text=sd_input(policy,mult,scale);old=set(sd_configurations(2));extra='';groups=[]
    for i,q in enumerate(configurations()):
        if q in old:continue
        name='closure'+str(i);groups.append(name);extra+=f"Config('{name}', '{config_text(q)}')\n"
    marker='SetConfigEnergyMode(-1)'
    if text.count(marker)!=1:raise ValueError('Closure configuration insertion point differs')
    text=text.replace(marker,extra+marker,1)
    pattern=r"^Structure\('o\.lev\.b', (\[.*\])\)$"
    matches=list(re.finditer(pattern,text,re.M))
    if len(matches)!=1:raise ValueError('Closure structure group list differs')
    m=matches[0];old_groups=ast.literal_eval(m[1]);new=f"Structure('o.lev.b', {old_groups+groups!r})"
    return text[:m.start()]+new+text[m.end():]


def inspect(text):
    rows=level_inventory(text);expected=expected_inventory(configurations())
    if (not re.search(r'^O Z\s*=\s*8\.0\s*$',text,re.M) or set(re.findall(r'^NELE\s*=\s*(\d+)',text,re.M))!={'8'}
        or Counter((r['parity'],r['two_J']) for r in rows)!=expected or set(r['index'] for r in rows)!=set(range(len(rows)))):
        raise ValueError('Full-space electron/parity/J/original-ID inventory differs')
    return rows


def main():
    ap=argparse.ArgumentParser()
    for n in ('source','baseline','out'):ap.add_argument('--'+n,type=Path,required=True)
    a=ap.parse_args()
    if platform.node().lower().split('.')[0]!='sirius':raise RuntimeError('Atomic CI closure runs belong on Sirius')
    for n,v in vars(a).items():setattr(a,n,v.resolve())
    parent=verify_receipt(a.baseline);probe=parent['probe'];binary=Path(probe['path'])
    if (parent['status']!='FAC_NO_CORRECTION_COMPLETED_REVIEW_REQUIRED' or parent['fac_pin']!=PIN
        or subprocess.check_output(['git','rev-parse','HEAD'],cwd=a.source,text=True).strip()!=PIN
        or subprocess.check_output(['git','diff','--name-only'],cwd=a.source,text=True).strip()
        or sha(binary)!=probe['sha256'] or sha(a.source/'libfac.a')!=probe['libfac_sha256'] or sha(a.source/'sfac/sfac')!=probe['upstream_binary_sha256']):raise ValueError('Closure pinned source/frontend/library differs')
    for policy in ('ground','equal'):
        for label,mult,scale in [('coherent10',.3,10),('coherent100',.03,100)]:
            p=a.baseline/(policy+'_'+label)/'levels'
            if (p/'input.sf').read_text()!=sd_input(policy,mult,scale) or not diagnostics((p/'sfac.stdout').read_text(),mult,256)['printed_rule_verified']:raise ValueError('Closure paired baseline differs')
    a.out.mkdir(parents=True,exist_ok=False);drivers=a.out/'drivers';drivers.mkdir()
    for p in Path(__file__).parent.glob('rya1354_*'):
        if p.is_file():shutil.copyfile(p,drivers/p.name)
    for name in ('doc/manual.tex','faclib/config.c','faclib/structure.c','faclib/radial.c','faclib/consts.h'):
        dest=a.out/'source_documents'/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(a.source/name,dest)
    receipt=dict(schema='rya1354.fac-full-valence-space.v28',status='RUNNING',started_utc=datetime.now(timezone.utc).isoformat(),fac_pin=PIN,
        baseline_receipt_sha256=sha(a.baseline/'receipt.json'),probe=probe,cases={},scientific_admission='HOLD',exposure='EXPOSED_NONBLIND',adoption_allowed=False,
        accepted_lifetime=None,uncertainty=None,covariance=None,calibrated_sigma=None,scope='All4valence electrons in fixed2p/3s/3p/3d set;1s2/2s2 fixed;zero config corrections;both potentials;paired saved potential byteidentical toSD;levels/basis only, no LS/rates.')
    try:
        from rya1354_fac_full_space_audit import transport
        from rya1354_fac_no_correction_audit import configuration_keys,verify_config_inventory
        for policy in ('ground','equal'):
            for label,mult,scale in [('coherent10',.3,10),('coherent100',.03,100)]:
                name=policy+'_'+label;out=a.out/name;out.mkdir();(out/'input.sf').write_text(input_text(policy,mult,scale))
                with (out/'sfac.stdout').open('xb') as stdout,(out/'sfac.stderr').open('xb') as stderr:
                    p=subprocess.run([str(binary),'input.sf'],cwd=out,stdout=stdout,stderr=stderr,timeout=900,check=False)
                if p.returncode:raise ValueError(f'Closure native exit{p.returncode}: '+name)
                d=diagnostics((out/'sfac.stdout').read_text(),mult,256)
                if not d['printed_rule_verified']:raise ValueError('Closure stopping guard failed: '+name)
                rows=inspect((out/'o.lev').read_text());s=read_files(out);keys=configuration_keys(s['o.basis'].decode())
                configs={n:verify_config_inventory(s['configs.'+n+'.tsv'].decode(),keys) for n in ('before','after')}
                if any(v['groups']!=31 for v in configs.values()):raise ValueError('Closure groups differ')
                if s['o.pot']!=(a.baseline/name/'levels/o.pot').read_bytes():raise ValueError('Closure changed paired saved potential: '+name)
                tr=transport(s['o.basis'].decode(),s['o_export_000.c'].decode(),s['o_export_000.cm'])
                receipt['cases'][name]=dict(diagnostics=d,configuration_state=configs,transport=tr,paired_saved_potential_byte_identical=True,paired_potential_sha256=sha(a.baseline/name/'levels/o.pot'))
                (a.out/'progress.json').write_text(json.dumps(receipt,indent=2)+'\n')
        receipt['status']='FAC_FULL_VALENCE_LEVELS_COMPLETED_REVIEW_REQUIRED'
    except Exception as exc:
        receipt.update(status='FAILED_RETAINED',error=type(exc).__name__+': '+str(exc));raise
    finally:
        receipt['files']={str(p.relative_to(a.out)):dict(sha256=sha(p),bytes=p.stat().st_size) for p in sorted(a.out.rglob('*')) if p.is_file()}
        (a.out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')

if __name__=='__main__':main()
