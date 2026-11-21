"""v53: GRASP MCDHF-EOL for O I in the Li et al. 2023 (arXiv 2304.02310) scheme, Sirius only.

Separate even and odd orbital sets (transitions through rbiotransform), as in Li 2023 Sect. 2.3:
  stage s0: Dirac-Fock EOL on the Table 1 "MR in MCDHF" configurations alone (0 excitations);
            all MR orbitals varied and spectroscopic.
  layer n (n = 3..6): MR-SD, 1s and 2s closed, substitutions from 2p and the outer electron into
            the MR orbitals plus every orbital with principal number <= n (v47 'mcdhf' families);
            only the NEW non-MR orbitals of shell n are varied (correlation orbitals); everything
            earlier is frozen.
EOL targets (Li 2023 Sect. 2.3 / Sect. 3.1): the 45 lowest even states (2p4, 2p3(4S)np n=3-6,
2p3(4S)nf n=4,5) and the 36 lowest odd states (2p3(4S)ns n=3-6, 2p3(4S)nd n=3-5, 2p3(2D)3s),
standard (2J+1) weights. Leading jj2lsj labels of every targeted level must reproduce those
families per J block, at every stage, or the run stops.

After each stage: rbiotransform + rtransition E1 (Babushkin and Coulomb) between all targeted
even and odd states. No observed energy, gf or lifetime is an input. No result is adopted.

A step exceeding its declared budget, the address-space cap or the disk guard is a MEASURED
OUTCOME that ends the run; any refusal, nonzero exit, non-converged orbital or failed check is
FAILED_RETAINED. The receipt carries a files inventory.
"""
import argparse
from collections import Counter
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import platform
import re
import shutil

from rya1354_li_scheme_sizing import inventory_stream, references
from rya1354_rci_feasibility import BudgetExceeded, MEM_CAP_BYTES, RWFN_STDIN, disk_guard, dir_bytes, timed
from rya1354_solver_smoke import sha
from rya1354_triplet_correlation import transitions
from rya1354_two_method_summary import scf_diagnostics

L='spdfgh'
LAYERS=(3,4,5,6)
TWOJ={'even':'0,10','odd':'0,8'}
BUDGET_S={'rcsfgenerate':3600,'rnucleus':600,'rangular':14400,'rwfnestimate':1800,'rmcdhf':86400,'jj2lsj':3600}
RMCDHF_ITERATIONS=200
LAYER_ESTIMATE='y\n1\nprevious.w\n*\n3\n*\n'   # previous stage's orbitals, screened hydrogenic for the new shell
JJ2LSJ_STDIN='state\nn\nn\ny\n'

# LS terms (2S+1, L, J values) of each target family; Li's 45 even / 36 odd lowest states.
TERMS={'ground':[('3P',(0,1,2)),('1D',(2,)),('1S',(0,))],
       '4S.p':[('5P',(1,2,3)),('3P',(0,1,2))],'4S.f':[('5F',(1,2,3,4,5)),('3F',(2,3,4))],
       '4S.s':[('5S',(2,)),('3S',(1,))],'4S.d':[('5D',(0,1,2,3,4)),('3D',(1,2,3))],
       '2D.s':[('3D',(1,2,3)),('1D',(2,))]}
TARGET_FAMILIES={'even':[('ground',2)]+[('4S.p',n) for n in (3,4,5,6)]+[('4S.f',n) for n in (4,5)],
                 'odd':[('4S.s',n) for n in (3,4,5,6)]+[('4S.d',n) for n in (3,4,5)]+[('2D.s',3)]}


def expected_targets(parity):
    """{J: Counter((family, n))} of the targeted states."""
    out={}
    for fam,n in TARGET_FAMILIES[parity]:
        for _,js in TERMS[fam]:
            for j in js:out.setdefault(j,Counter())[(fam,n)]+=1
    return out


def roots(parity):
    return {j:sum(c.values()) for j,c in sorted(expected_targets(parity).items())}


def mr_max(parity):
    """Highest principal number per l among the MR orbitals (from the v47 'mcdhf' families)."""
    kept,_=references('mcdhf',parity,6)
    top={}
    for line in kept:
        for n,l in re.findall(r'(\d+)([spdfgh])\(',line):top[l]=max(top.get(l,0),int(n))
    return top


def active_set(parity,n):
    top=dict(mr_max(parity))
    for i,l in enumerate(L):
        if i<n:top[l]=max(top.get(l,0),n)
    return ','.join(f'{top[l]}{l}' for l in L if l in top)


def new_orbitals(parity,n):
    top=mr_max(parity)
    # GRASP GETRSL: '3d' alone selects only 3d (j=5/2); '3d*' selects 3d- and 3d.
    return [f'{n}{l}*' for i,l in enumerate(L) if i<n and top.get(l,0)<n]


def generator_stdin(parity,n):
    kept,omitted=references('mcdhf',parity,6)
    if omitted:raise ValueError('Li MCDHF MR must be complete')
    if n==0:return '*\n0\n'+'\n'.join(kept)+'\n\n'+active_set(parity,0)+'\n'+TWOJ[parity]+'\n0\nn\n'
    return '*\n0\n'+'\n'.join(kept)+'\n\n'+active_set(parity,n)+'\n'+TWOJ[parity]+'\n2\nn\n'


def block_order(path):
    """J of each native CSF block, in file order."""
    order=[];buf=[];seen=False
    with open(path) as f:
        for _ in range(5):next(f)
        for line in f:
            if line.startswith('*') or line.startswith(' *'):seen=False;continue
            buf.append(line.rstrip('\n'))
            if len(buf)==3:
                if not seen:
                    lab=buf[2].split()[-1]
                    order.append(int(lab[:-1]) if lab[:-1].isdigit() else lab);seen=True
                buf=[]
    return order


def levels_stdin(order,parity):
    want=roots(parity)
    if sorted(order)!=sorted(want) or len(set(order))!=len(order):
        raise ValueError(f'CSF blocks {order} differ from targeted J {sorted(want)}')
    return ''.join(','.join(str(i) for i in range(1,want[j]+1))+'\n' for j in order)


LEVEL=re.compile(r'^\s*(\d+)\s+(\d+)\s+([+-])\s+(-?\d+\.\d+)\s+[\d.]+%\s*$')


def family(label):
    if '2p(4)' in label:return ('ground',2)
    m=re.search(r'2p\(3\)(4S|2D|2P)\d?\.(\d+)([spdf])_',label)
    return (f'{m[1]}.{m[3]}',int(m[2])) if m else ('other',0)


def identity(text,parity):
    """Leading-label families per J must equal the targeted Li states."""
    got={};lines=text.splitlines();i=0
    while i<len(lines):
        m=LEVEL.match(lines[i])
        if m:
            if m[3]!=('+' if parity=='even' else '-'):raise ValueError('Wrong level parity in lsj.lbl')
            got.setdefault(int(m[2]),Counter())[family(lines[i+1].split()[-1])]+=1
        i+=1
    if got!=expected_targets(parity):
        raise ValueError(f'Leading physical state identity differs ({parity}): {dict(got)}')
    return {j:{f'{k[0]}{k[1]}':v for k,v in c.items()} for j,c in sorted(got.items())}


def rmcdhf_stdin(levels,varied,spectroscopic):
    return f'y\n{levels}5\n{varied}\n{spectroscopic}\n{RMCDHF_ITERATIONS}\n'


def stage(a,receipt,parity,n,previous):
    name=f's0_{parity}' if n==0 else f'n{n}_{parity}'
    out=a.out/name;out.mkdir();case=dict(parity=parity,layer=n,steps={})
    receipt['stages'][name]=case
    case['disk_free_before']=disk_guard(a.out)
    def step(exe,stdin,bin_=None):
        case['steps'][exe]=timed((bin_ or a.grasp_bin)/exe,stdin,out,BUDGET_S[exe],MEM_CAP_BYTES)
    step('rcsfgenerate',generator_stdin(parity,n))
    case.update(active_set=active_set(parity,n),**inventory_stream(out/'rcsf.out',parity,'mcdhf'))
    order=block_order(out/'rcsf.out')
    levels=levels_stdin(order,parity)
    for j,k in roots(parity).items():
        if case['blocks'][f'{j}{"+" if parity=="even" else "-"}']<k:raise ValueError(f'Block J={j} has fewer CSFs than targets')
    if n==6 and a.v47:
        for j,c in a.v47[f'mcdhf_{parity}_n6']['blocks'].items():
            if case['blocks'].get(j)!=c:raise ValueError(f'n6 {parity} block {j} differs from v47: {case["blocks"].get(j)} != {c}')
        case['v47_n6_blocks_match']=True
    os.link(out/'rcsf.out',out/'rcsf.inp')
    step('rnucleus','8\n0\n0\n0\n0\n0\n')
    step('rangular','y\n')
    case['angular_bytes']=dir_bytes(out,'mcp.')
    case['disk_free_after_rangular']=disk_guard(a.out)
    if n==0:
        step('rwfnestimate',RWFN_STDIN)
        varied,spec='*','*'
    else:
        os.link(previous/'state.w',out/'previous.w')
        step('rwfnestimate',LAYER_ESTIMATE)
        varied,spec=','.join(new_orbitals(parity,n)),''
    case['varied']=varied
    step('rmcdhf',rmcdhf_stdin(levels,varied,spec),a.strict_bin)
    diag=scf_diagnostics((out/'rmcdhf.stdout').read_text(),(out/'rmcdhf.stderr').read_text())
    case['scf']=diag
    if diag['iteration_limit_hit'] or not diag['all_reported_functions_below_ACC']:
        raise ValueError(f'Native orbital criterion unmet in {name}')
    for src,ext in (('rcsf.inp','c'),('rwfn.out','w'),('rmix.out','m')):os.link(out/src,out/f'state.{ext}')
    step('jj2lsj',JJ2LSJ_STDIN)
    case['leading_families']=identity((out/'state.lsj.lbl').read_text(),parity)
    (a.out/'progress.json').write_text(json.dumps({k:v for k,v in receipt.items() if k!='files'},indent=2)+'\n')
    return out


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--grasp-bin',type=Path,required=True)
    ap.add_argument('--strict-bin',type=Path,required=True)
    ap.add_argument('--v47-receipt',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    a=ap.parse_args()
    if platform.node().lower().split('.')[0]!='sirius':raise RuntimeError('Atomic solvers belong on Sirius')
    for k,v in vars(a).items():
        if isinstance(v,Path):setattr(a,k,v.resolve())
    a.out.mkdir(parents=True,exist_ok=False)
    shutil.copyfile(__file__,a.out/'driver.py')
    a.v47=json.loads(a.v47_receipt.read_text())['cases']
    receipt=dict(schema='rya1354.mcdhf-eol.v53',status='RUNNING',scientific_admission='HOLD',exposure='EXPOSED_NONBLIND',
        adoption_allowed=False,started_utc=datetime.now(timezone.utc).isoformat(),stages={},transitions={},
        v47_receipt_sha256=sha(a.v47_receipt),budgets_s=BUDGET_S,mem_cap_bytes=MEM_CAP_BYTES,
        targets={p:{str(j):k for j,k in roots(p).items()} for p in ('even','odd')},
        source='Li et al. 2023 arXiv2304.02310 Table 1 MR in MCDHF, Sect. 2.3 EOL targets; no fitted or observed inputs',
        model='Point static Z=8 nucleus, Dirac-Coulomb MCDHF-EOL, standard weights; separate even/odd orbitals; strict orbital-convergence rmcdhf',
        stopped=None)
    try:
        prev={}
        for n in (0,*LAYERS):
            for parity in ('even','odd'):prev[parity]=stage(a,receipt,parity,n,prev.get(parity))
            key=f'{"s0" if n==0 else f"n{n}"}_E1';t=a.out/key
            receipt['transitions'][key]=transitions(a.grasp_bin,prev['even'],prev['odd'],t,'E1',False)
            receipt['last_completed_stage']=key
        receipt['status']='MCDHF_EOL_COMPLETED_REVIEW_REQUIRED'
    except BudgetExceeded as exc:
        receipt.update(status='MEASURED_STOP_REVIEW_REQUIRED',stopped=str(exc))
    except Exception as exc:
        receipt.update(status='FAILED_RETAINED',error=f'{type(exc).__name__}: {exc}');raise
    finally:
        receipt['finished_utc']=datetime.now(timezone.utc).isoformat()
        seen={};files={}
        for p in sorted(a.out.rglob('*')):
            if p.is_file() and p.name!='receipt.json':
                key=(p.stat().st_dev,p.stat().st_ino)
                if key not in seen:seen[key]=dict(sha256=sha(p),bytes=p.stat().st_size)
                files[str(p.relative_to(a.out))]=seen[key]
        receipt['files']=files
        (a.out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')


if __name__=='__main__':
    main()
