"""Declared balanced n3 valence spaces; determinant-count guard; Sirius runner."""
import argparse
from collections import Counter
from datetime import datetime, timezone
from itertools import product,combinations
import json
from pathlib import Path
import platform
import re
import shutil
import subprocess
from rya1354_fac_bootstrap import PIN,sha,level_inventory
from rya1354_fac_coherent import input_text as minimal_input
from rya1354_fac_controls import diagnostics

SHELLS=(('2p',1,6),('3s',0,2),('3p',1,6),('3d',2,10))
REFERENCES=((4,0,0,0),(3,1,0,0))


def configurations(rank):
    if rank not in (1,2):raise ValueError('Only declared S/SD spaces')
    return [q for q in product(*(range(min(cap,4)+1) for _,_,cap in SHELLS))
            if sum(q)==4 and min(sum(abs(a-b) for a,b in zip(q,r))//2 for r in REFERENCES)<=rank]


def config_text(q):
    return '2s2 '+' '.join(name+str(n) for (name,_,_),n in zip(SHELLS,q) if n)


def expected_inventory(configs):
    """Independent magnetic determinant counting: N(J)=N(M=J)-N(M=J+1)."""
    result=Counter()
    for q in configs:
        states=Counter({0:1});parity=0
        for (_,ell,_),n in zip(SHELLS,q):
            parity+=ell*n
            ms=[2*m+s for m in range(-ell,ell+1) for s in (-1,1)]
            shell=Counter(sum(c) for c in combinations(ms,n))
            states=Counter({total:sum(a*b for x,a in states.items() for y,b in shell.items() if x+y==total)
                            for total in {x+y for x in states for y in shell}})
        for j in range(0,max(states)+1,2):
            count=states[j]-states[j+2]
            if count:result[(parity%2,j)]+=count
    return result


def input_text(rank,multiplier,scale):
    configs=configurations(rank);groups=['ground','odd3s'];extra=''
    for i,q in enumerate(configs):
        if q in REFERENCES:continue
        group='corr'+str(i);groups.append(group)
        extra+=f"Config('{group}', '{config_text(q)}')\n"
    text=minimal_input(multiplier,scale).replace('ConfigEnergy(0)',extra+'ConfigEnergy(0)',1)
    text=text.replace("Structure('o.lev.b', ['ground', 'odd3s'])",f"Structure('o.lev.b', {groups!r})")
    return text


def inspect(text,rank):
    rows=level_inventory(text)
    if (not re.search(r'^O Z\s*=\s*8\.0\s*$',text,re.M)
            or set(re.findall(r'^NELE\s*=\s*(\d+)',text,re.M))!={'8'}
            or Counter((r['parity'],r['two_J']) for r in rows)!=expected_inventory(configurations(rank))):
        raise ValueError('Expanded O electron/parity/J inventory differs from determinant count')
    if set(r['index'] for r in rows)!=set(range(len(rows))):raise ValueError('Expanded level IDs incomplete')
    return rows


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    for name in ('source','bootstrap','out'):ap.add_argument('--'+name,type=Path,required=True)
    args=ap.parse_args()
    if platform.node().lower().split('.')[0]!='sirius':raise RuntimeError('Atomic calculations belong on Sirius')
    args.source=args.source.resolve();args.bootstrap=args.bootstrap.resolve();args.out=args.out.resolve()
    build=json.loads((args.bootstrap/'receipt.json').read_text());binary=args.source/'sfac/sfac'
    if (subprocess.check_output(['git','rev-parse','HEAD'],cwd=args.source,text=True).strip()!=PIN
            or subprocess.check_output(['git','diff','--name-only'],cwd=args.source,text=True).strip()
            or sha(binary)!=build['binary']['sha256']):raise ValueError('FAC source/executable differs')
    args.out.mkdir(parents=True,exist_ok=False);drivers=args.out/'drivers';drivers.mkdir()
    for p in Path(__file__).parent.glob('rya1354_fac*.py'):shutil.copyfile(p,drivers/p.name)
    for name in ('doc/manual.tex','faclib/consts.h','faclib/structure.c'):
        dest=args.out/'source_documents'/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(args.source/name,dest)
    receipt=dict(schema='rya1354.fac-space.v23',status='RUNNING',pin=PIN,binary_sha256=sha(binary),
        started_utc=datetime.now(timezone.utc).isoformat(),cases={},scientific_admission='HOLD',exposure='EXPOSED_NONBLIND',
        adoption_allowed=False,accepted_lifetime=None,uncertainty=None,covariance=None,calibrated_sigma=None,
        scope='Balanced n3 S/SD from both references;1s/2s fixed, ground-only potential optimization; no rates/LS identity/accuracy admission.')
    try:
        for rank in (1,2):
            for label,mult,scale in [('coherent10',.3,10),('coherent100',.03,100)]:
                name=f'n3_r{rank}_{label}';out=args.out/name;out.mkdir();(out/'input.sf').write_text(input_text(rank,mult,scale))
                with (out/'sfac.stdout').open('xb') as stdout,(out/'sfac.stderr').open('xb') as stderr:
                    p=subprocess.run([str(binary),'input.sf'],cwd=out,stdout=stdout,stderr=stderr,timeout=900,check=False)
                if p.returncode:raise ValueError(f'Native FAC exit{p.returncode}: '+name)
                d=diagnostics((out/'sfac.stdout').read_text(),mult,256)
                if not d['printed_rule_verified']:raise ValueError('Expanded stopping guard failed: '+name)
                rows=inspect((out/'o.lev').read_text(),rank)
                if any(not (out/f).is_file() for f in ('o.basis','o_export_000.c','o_export_000.cm','o.pot')):
                    raise ValueError('Expanded basis/potential absent')
                receipt['cases'][name]=dict(rank=rank,configurations=[list(q) for q in configurations(rank)],
                    diagnostics=d,levels=rows,returncode=p.returncode)
                (args.out/'progress.json').write_text(json.dumps(receipt,indent=2)+'\n')
        receipt['status']='FAC_EXPANDED_LEVELS_COMPLETED_REVIEW_REQUIRED'
    except Exception as exc:
        receipt.update(status='FAILED_RETAINED',error=type(exc).__name__+': '+str(exc));raise
    finally:
        receipt['files']={str(p.relative_to(args.out)):dict(sha256=sha(p),bytes=p.stat().st_size)
                          for p in sorted(args.out.rglob('*')) if p.is_file()}
        (args.out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')

if __name__=='__main__':main()
