"""v47: size the Li et al. 2023 (arXiv 2304.02310, Table 1) O I computational scheme on
Sirius, CSF GENERATION ONLY (no orbitals, energies, rates or grades).

Two stages per parity and active-set layer n_max:
  mcdhf: Table 1 "MR in MCDHF", 1s and 2s closed, SD from 2p and the outer valence electron;
  rci:   Table 1 "MR in RCI" (adds np to n=8, ns to n=7 and the 2s2p4 ns/nd families),
         1s closed, SD from every n=2 subshell and the outer valence electron.
MR orbitals above n_max are omitted (recorded per case). Active sets are n_max{s,p,d,f,g,h}
truncated at l<=n_max-1; Li's final set is {12s,12p,11d,11f,9g,7h}.
"""
import argparse
from collections import Counter
from datetime import datetime, timezone
import json
from pathlib import Path
import platform
import re
import shutil
import time

from rya1354_solver_smoke import run, sha

L='spdfgh'
TWOJ={'even':'0,4','odd':'2,4'}   # 2p4 3P0,1,2/1D2/1S0; 2p3 3s 3S1 and 5S2
LAYERS=(6,7,8,9)


def references(stage,parity,n_max):
    """Native rcsfgenerate reference lines for one stage/parity, and the omitted MR configurations."""
    if stage not in ('mcdhf','rci') or parity not in ('even','odd'):raise ValueError('Unknown stage/parity')
    s2='1s(2,i)2s(2,'+('i' if stage=='mcdhf' else '*')+')'
    if parity=='even':
        fam=[(s2+'2p(4,*)',2)]
        fam+=[(s2+f'2p(3,*){n}p(1,*)',n) for n in range(3,(6 if stage=='mcdhf' else 8)+1)]
        fam+=[(s2+f'2p(3,*){n}f(1,*)',n) for n in (4,5)]
        if stage=='rci':
            fam+=[(f'1s(2,i)2s(1,*)2p(4,*){n}s(1,*)',n) for n in range(3,8)]
            fam+=[(f'1s(2,i)2s(1,*)2p(4,*){n}d(1,*)',n) for n in range(3,7)]
    else:
        fam=[(s2+f'2p(3,*){n}s(1,*)',n) for n in range(3,(6 if stage=='mcdhf' else 7)+1)]
        fam+=[(s2+f'2p(3,*){n}d(1,*)',n) for n in range(3,6)]
    kept=[c for c,n in fam if n<=n_max];omitted=[c for c,n in fam if n>n_max]
    return kept,omitted


def active(n_max):
    return ','.join(f'{n_max}{l}' for i,l in enumerate(L) if i<n_max)


def stdin_for(stage,parity,n_max):
    kept,_=references(stage,parity,n_max)
    return '*\n0\n'+'\n'.join(kept)+'\n\n'+active(n_max)+'\n'+TWOJ[parity]+'\n2\nn\n'


ORB=re.compile(r'(\d+)([spdfgh])(-?)\s*\(\s*(\d+)\)')


def inventory_stream(path,parity,stage):
    """Streaming version of rya1354_mr_inventory.inventory: native 3-record CSFs, 8 electrons,
    closed 1s (and 2s for mcdhf), block parity; counts per J block. Duplicates are checked by
    the SHA-1 of each record to bound memory."""
    import hashlib
    blocks=Counter();seen=set();buf=[]
    with open(path) as f:
        head=[next(f).rstrip('\n') for _ in range(5)]
        if head[0]!='Core subshells:' or head[1].strip() or head[4]!='CSF(s):':
            raise ValueError('Unexpected native CSF header/core representation')
        for line in f:
            line=line.rstrip('\n')
            if line.startswith('*') or line.startswith(' *'):
                if buf:raise ValueError('Incomplete native CSF record at block end')
                continue
            buf.append(line)
            if len(buf)<3:continue
            key=hashlib.sha1('\n'.join(buf).encode()).digest()
            if key in seen:raise ValueError('Duplicate native CSF')
            seen.add(key)
            occ={}
            for n,l,minus,k in ORB.findall(buf[0]):occ[n+l+minus]=int(k)
            tot=sum(occ.values())
            s1=occ.get('1s',0);s2=occ.get('2s',0)
            if tot!=8 or s1!=2 or (stage=='mcdhf' and s2!=2):
                raise ValueError('Electron/core occupation differs from declared model')
            odd=sum(L.index(o[1] if o[1] in L else re.search('[spdfgh]',o)[0])*k for o,k in occ.items())%2
            label=buf[2].split()[-1]
            if odd!=(parity=='odd') or not label.endswith('+' if parity=='even' else '-'):
                raise ValueError('Wrong native occupation/block parity')
            blocks[label]+=1;buf=[]
    if buf:raise ValueError('Incomplete native CSF records')
    return dict(CSFs=sum(blocks.values()),blocks=dict(sorted(blocks.items())),duplicate_CSFs=0)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--grasp-bin',type=Path,required=True);ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--timeout',type=int,default=3600);a=ap.parse_args()
    if platform.node().lower().split('.')[0]!='sirius':raise RuntimeError('Atomic native generators belong on Sirius')
    a.grasp_bin=a.grasp_bin.resolve();a.out.mkdir(parents=True,exist_ok=False)
    shutil.copyfile(__file__,a.out/'driver.py')
    receipt=dict(schema='rya1354.li-scheme-sizing.v47',status='RUNNING',scientific_admission='HOLD',
        started_utc=datetime.now(timezone.utc).isoformat(),steps=[],cases={},
        source='Li et al. 2023 arXiv2304.02310 Table 1 MR/AS; reduced layers n_max 6-9; no fitted or observed inputs',
        scope='CSF_GENERATION_ONLY_NO_WAVEFUNCTION_ENERGY_RATE_OR_GRADE')
    try:
        for n_max in LAYERS:
            for stage in ('mcdhf','rci'):
                for parity in ('even','odd'):
                    name=f'{stage}_{parity}_n{n_max}';out=a.out/name;out.mkdir()
                    t0=time.monotonic()
                    step=run(a.grasp_bin/'rcsfgenerate',stdin_for(stage,parity,n_max),out,a.timeout)
                    step['wall_s']=round(time.monotonic()-t0,3);receipt['steps'].append(step)
                    if step['semantic_errors']:raise ValueError('Native generator refusal: '+name)
                    kept,omitted=references(stage,parity,n_max)
                    receipt['cases'][name]=dict(references=kept,omitted_references=omitted,active=active(n_max),
                        twoJ=TWOJ[parity],rcsf_bytes=(out/'rcsf.out').stat().st_size,rcsf_sha256=sha(out/'rcsf.out'),
                        wall_s=step['wall_s'],**inventory_stream(out/'rcsf.out',parity,stage))
                    (a.out/'progress.json').write_text(json.dumps(receipt,indent=2)+'\n')
        receipt['status']='SIZING_COMPLETED_NO_ATOMIC_PREDICTION'
    except Exception as exc:
        receipt.update(status='FAILED_RETAINED',error=f'{type(exc).__name__}: {exc}');raise
    finally:
        receipt['finished_utc']=datetime.now(timezone.utc).isoformat()
        (a.out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')


if __name__=='__main__':main()
