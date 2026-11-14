"""Both-gauge expanded-space FAC raw rates after full native LS replay; Sirius."""
import argparse
from datetime import datetime,timezone
import json
from pathlib import Path
import platform
import shutil
import subprocess
from rya1354_fac_bootstrap import PIN,sha,level_inventory
from rya1354_fac_rates import verify_receipt,read_files
from rya1354_fac_expanded_ls_audit import review
from rya1354_fac_space import input_text,inspect
from rya1354_fac_controls import diagnostics
from rya1354_fac_recouple import CONTROLS,GRASP_PIN


def rate_input(rank,mult,scale,gauge,mapping):
    upper,=mapping['targets']['3S'];states=mapping['levels'];u=states[str(upper)]
    lower=sorted(int(i) for i,r in states.items() if r['energy_au']<u['energy_au'])
    if not lower:raise ValueError('No expanded lower states')
    requests=[];groups={}
    for low in lower:
        r=states[str(low)];ju=u['two_J'];jl=r['two_J']
        for k2 in range(max(2,abs(ju-jl)),ju+jl+1,2):
            k=k2//2
            electric=((r['parity']!=u['parity'])==(k%2==1))
            m=-k if electric else k
            groups.setdefault(m,[]).append(low)
            requests.append(dict(upper=upper,lower=low,multipole=m,operator=('E' if electric else 'M')+str(k)))
    text=input_text(rank,mult,scale)+f'SetTransitionOptions({gauge},0)\nSetTransitionCut(0,0)\n'
    for m,lows in sorted(groups.items(),key=lambda x:(abs(x[0]),x[0])):
        text+=f"TransitionTable('rates.tr.b',{lows},[{upper}],{m})\n"
        for low in lows:text+=f"RawTransition('raw.tsv',{low},{upper},{m})\n"
    text+="PrintTable('rates.tr.b','rates.tr',1)\n"
    return text,requests,lower


def main():
    ap=argparse.ArgumentParser()
    for n in ('source','states','mapped','rate_parent','out'):ap.add_argument('--'+n,type=Path,required=True)
    a=ap.parse_args()
    if platform.node().lower().split('.')[0]!='sirius':raise RuntimeError('Atomic rates belong on Sirius')
    for n,v in vars(a).items():setattr(a,n,v.resolve())
    sreceipt=verify_receipt(a.states);mreceipt=verify_receipt(a.mapped);rreceipt=verify_receipt(a.rate_parent)
    if (sreceipt['status']!='FAC_EXPANDED_LEVELS_COMPLETED_REVIEW_REQUIRED' or mreceipt['status']!='FAC_ANGULAR_TRANSFORM_COMPLETED_REVIEW_REQUIRED'
        or mreceipt['source_receipt_sha256']!=sha(a.states/'receipt.json') or mreceipt['grasp_pin']!=GRASP_PIN
        or rreceipt['status']!='FAC_MODEL_RATES_COMPLETED_REVIEW_REQUIRED' or sreceipt['pin']!=PIN
        or subprocess.check_output(['git','rev-parse','HEAD'],cwd=a.source,text=True).strip()!=PIN
        or subprocess.check_output(['git','diff','--name-only'],cwd=a.source,text=True).strip()):raise ValueError('Expanded provenance differs')
    probe=rreceipt['probe'];binary=Path(probe['path'])
    if sha(binary)!=probe['sha256'] or sha(a.source/'libfac.a')!=probe['libfac_sha256'] or sha(a.source/'sfac/sfac')!=probe['upstream_binary_sha256']:raise ValueError('Unchanged diagnostic frontend/library differs')
    prepared={}
    for rank in (1,2):
        for label,mult,scale in [('coherent10',.3,10),('coherent100',.03,100)]:
            name=f'n3_r{rank}_{label}';s=read_files(a.states/name);c=read_files(a.mapped/name)
            if (s['input.sf'].decode()!=input_text(rank,mult,scale) or not diagnostics(s['sfac.stdout'].decode(),mult,256)['printed_rule_verified']
                or c['state.c']!=s['o_export_000.c'] or c['state.cm']!=s['o_export_000.cm'] or c['jj2lsj.stdin'].decode()!=CONTROLS):raise ValueError('Expanded source input differs')
            inspect(s['o.lev'].decode(),rank);prepared[name]=(rank,mult,scale,s,review(s,c,rank))
    a.out.mkdir(parents=True,exist_ok=False);drivers=a.out/'drivers';drivers.mkdir()
    for p in Path(__file__).parent.glob('rya1354_*.py'):shutil.copyfile(p,drivers/p.name)
    receipt=dict(schema='rya1354.fac-expanded-rates.v25',status='RUNNING',started_utc=datetime.now(timezone.utc).isoformat(),
        source_receipt_sha256=sha(a.states/'receipt.json'),mapping_receipt_sha256=sha(a.mapped/'receipt.json'),rate_parent_receipt_sha256=sha(a.rate_parent/'receipt.json'),probe=probe,
        cases={},scientific_admission='HOLD',exposure='EXPOSED_NONBLIND',adoption_allowed=False,accepted_lifetime=None,uncertainty=None,covariance=None,calibrated_sigma=None,
        scope='Fresh all-lower-state/operator enumeration from full LS mapping; bothgauges; native status never zero-filled; no complete lifetime until independent raw/table audit.')
    try:
        for name,(rank,mult,scale,s,mapping) in prepared.items():
            for suffix,gauge in [('B',2),('C',1)]:
                case=name+'_'+suffix;out=a.out/case;out.mkdir();text,requests,lower=rate_input(rank,mult,scale,gauge,mapping);(out/'input.sf').write_text(text)
                with (out/'sfac.stdout').open('xb') as stdout,(out/'sfac.stderr').open('xb') as stderr:
                    process=subprocess.run([str(binary),'input.sf'],cwd=out,stdout=stdout,stderr=stderr,timeout=900,check=False)
                if process.returncode:raise ValueError(f'Expanded rate native exit{process.returncode}: '+case)
                d=diagnostics((out/'sfac.stdout').read_text(),mult,256)
                if not d['printed_rule_verified']:raise ValueError('Expanded rates stopping guard failed')
                for f in ('o_export_000.c','o_export_000.cm'):
                    if (out/f).read_bytes()!=s[f]:raise ValueError('Expanded rate state differs from mapping')
                raw=[l for l in (out/'raw.tsv').read_text().splitlines() if l.strip() and not l.startswith('#')]
                if len(raw)!=len(requests):raise ValueError('Incomplete expanded raw queries')
                receipt['cases'][case]=dict(gauge=gauge,mode=0,lower_states=lower,requests=requests,diagnostics=d,raw_records=len(raw))
                (a.out/'progress.json').write_text(json.dumps(receipt,indent=2)+'\n')
        receipt['status']='FAC_EXPANDED_RATES_COMPLETED_REVIEW_REQUIRED'
    except Exception as exc:
        receipt.update(status='FAILED_RETAINED',error=type(exc).__name__+': '+str(exc));raise
    finally:
        receipt['files']={str(p.relative_to(a.out)):dict(sha256=sha(p),bytes=p.stat().st_size) for p in sorted(a.out.rglob('*')) if p.is_file()}
        (a.out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')

if __name__=='__main__':main()
