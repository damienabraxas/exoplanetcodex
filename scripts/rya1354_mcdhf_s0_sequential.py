"""v58: sequential spectroscopic Dirac-Fock for the Li 2023 O I MCDHF MR (Sirius only).

v53-v57: varying every MR orbital together leaves rotations within {1s,2s} and within each complete
Rydberg manifold free (the s0 CSF space is closed under them), and the strict orbital criterion
never converges; leading labels are not invariant under those rotations either (v57 even_A).
v58 removes the freedom the standard way: per parity, a core step, then ONE new spectroscopic
orbital per step with everything earlier frozen, each step on the MR configurations included so
far (all parents), EOL on the targeted Li states those configurations carry (lowest K per J,
standard weights). The last step of each parity holds the full s0 MR with all 45 even / 36 odd
targets.

  even: core 2p4 (vary 1s,2s,2p) -> 3p -> 4p -> 5p -> 6p -> 4f -> 5f
  odd:  core 2p3 3s (vary 1s,2s,2p,3s) -> 4s -> 5s -> 6s -> 3d -> 4d -> 5d

Every step: strict rmcdhf, method 3 for the varied orbitals (v56 non-default path), 500 cycles; the
strict criterion must hold for the varied orbitals; leading jj2lsj families per J must equal the
included targets. Then E1 Babushkin/Coulomb between the final even and odd sets. Refusals,
non-convergence, native stops, prompt mismatches, identity failures: FAILED_RETAINED (stop).
No observed input; nothing adopted.
"""
import argparse
from collections import Counter
import math
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import platform
import re
import shutil

from rya1354_li_scheme_sizing import inventory_stream, references
from rya1354_mcdhf_eol import BUDGET_S, TERMS, TARGET_FAMILIES, block_order, family, LEVEL
from rya1354_mcdhf_s0_controls import prompt_check
from rya1354_rci_feasibility import BudgetExceeded, MEM_CAP_BYTES, disk_guard, timed
from rya1354_solver_smoke import sha
from rya1354_triplet_correlation import transitions
from rya1354_two_method_summary import scf_diagnostics

NSCF=500
# v59: the lifetime ladder's declared half grid (v58 FAILED_RETAINED: 1s floor 2.16e-7 > ACCY on the default
# H=1/16, N=220 grid in every v53-v58 run). Same RNT and ACCY, H=0.03125, N=439, in rwfnestimate and rmcdhf.
RNT=math.exp(-65/16)/8
H=0.03125
N=439
ACCY=.0625**6
GRID=f'{RNT:.17g}\n{H:.17g}\n0\n{N}\n'
FIRST_ESTIMATE=f'n\nn\n\ny\nn\ny\n{GRID}2\n*\n3\n*\nn\n'               # TF, then screened hydrogenic
STEP_ESTIMATE=f'n\nn\n\ny\nn\ny\n{GRID}1\nprevious.w\n*\n2\n*\n3\n*\nn\n'  # frozen earlier orbitals first


def step_estimate(n_previous):
    """v60: read frozen orbitals from every earlier step of the parity (newest first; frozen orbitals are carried
    unchanged), then TF, then screened hydrogenic for the new orbital(s)."""
    return f'n\nn\n\ny\nn\ny\n{GRID}'+''.join(f'1\nprev{i}.w\n*\n' for i in range(n_previous))+'2\n*\n3\n*\nn\n'
# v60: each Rydberg orbital is optimised on its OWN l-series only (2p3 n'l, n'<=n, all parents; targets the
# 4S-parent states of that series), with the core fixed; lower orbitals of the same l stay in the list so the
# new orbital is kept orthogonal to them. v59 FAILED_RETAINED: with the 2p4 ground states in the 4f step's
# EOL, 4f (j=7/2) collapsed to a compact correlation-like function (e=0.62 Ha, 0 nodes).
# (label, target families of this step's list, varied orbitals); the final step is a fixed-orbital CI.
def _series(fam,ns,first=()):
    return [(f'{n}{fam[-1]}',list(first)+[(fam,m) for m in ns[:i+1]],[f'{n}{fam[-1]}']) for i,n in enumerate(ns)]
STEPS={'even':[('core',[('ground',2)],['1s','2s','2p'])]+_series('4S.p',(3,4,5,6))+_series('4S.f',(4,5))
              +[('ci',list(TARGET_FAMILIES['even']),[])],
       'odd':[('core',[('4S.s',3),('2D.s',3)],['1s','2s','2p','3s'])]
              +[(f'{n}s',[('4S.s',3),('2D.s',3)]+[('4S.s',m) for m in range(4,n+1)],[f'{n}s']) for n in (4,5,6)]
              +_series('4S.d',(3,4,5))+[('ci',list(TARGET_FAMILIES['odd']),[])]}


def included(parity,k):
    """Target families of step k (0-based)."""
    return STEPS[parity][k][1]


def expected(families):
    out={}
    for fam,n in families:
        for _,js in TERMS[fam]:
            for j in js:out.setdefault(j,Counter())[(fam,n)]+=1
    return out


def config_lines(parity,families):
    """MR configuration lines (all parents) carrying the included families, in MR order."""
    kept,_=references('mcdhf',parity,6)
    want=set()
    for fam,n in families:
        want.add('2p(4,*)' if fam=='ground' else f'2p(3,*){n}{fam[-1]}(1,*)')
    lines=[c for c in kept if any(c.endswith(w) for w in want)]
    if len(lines)!=len(want):raise ValueError(f'MR lines for {families} not found')
    return lines


def twoj_range(families):
    js=sorted(expected(families))
    if js!=list(range(js[0],js[-1]+1)):raise ValueError(f'Targeted J values are not contiguous: {js}')
    return f'{2*js[0]},{2*js[-1]}'


def generator_stdin(parity,k):
    fams=included(parity,k);lines=config_lines(parity,fams)
    return '*\n0\n'+'\n'.join(lines)+'\n\n'+active(lines)+'\n'+twoj_range(fams)+'\n0\nn\n'


def active(lines):
    top={}
    for line in lines:
        for n,l in re.findall(r'(\d+)([spdfgh])\(',line):top[l]=max(top.get(l,0),int(n))
    return ','.join(f'{top[l]}{l}' for l in 'spdfgh' if l in top)


def levels_stdin(order,families):
    want={j:sum(c.values()) for j,c in expected(families).items()}
    if sorted(order)!=sorted(want) or len(set(order))!=len(order):
        raise ValueError(f'CSF blocks {order} differ from targeted J {sorted(want)}')
    return ''.join(','.join(str(i) for i in range(1,want[j]+1))+'\n' for j in order)


def rmcdhf_stdin(levels,varied):
    """v56/v57 non-default path (method 3 for the varied orbitals); varied orbitals are spectroscopic."""
    v=','.join(f'{o}*' for o in varied)   # [] -> empty lists: every orbital fixed (CI)
    return (f'n\nn\ny\nn\ny\n{GRID}y\n{ACCY:.17g}\n{levels}5\n{v}\n{v}\n{NSCF}\n'
            f'y\nn\ny\n\n\n{v}\n\n'
            'n\nn\nn\nn\nn\nn\n1\n')


def identity(text,parity,families):
    got={}
    lines=text.splitlines()
    for i,line in enumerate(lines):
        m=LEVEL.match(line)
        if m:
            if m[3]!=('+' if parity=='even' else '-'):raise ValueError('Wrong level parity in lsj.lbl')
            got.setdefault(int(m[2]),Counter())[family(lines[i+1].split()[-1])]+=1
    if got!=expected(families):raise ValueError(f'Leading physical state identity differs ({parity}): {dict(got)}')
    return {j:{f'{k[0]}{k[1]}':v for k,v in c.items()} for j,c in sorted(got.items())}


def step(a,receipt,parity,k,previous):
    """previous: earlier step directories of this parity, newest first."""
    label,_,varied=STEPS[parity][k];fams=included(parity,k)
    name=f'{parity}_{k}_{label}';out=a.out/name;out.mkdir()
    case=dict(parity=parity,step=k,label=label,varied=varied,families=[f'{f}{n}' for f,n in fams],steps={})
    receipt['steps'][name]=case;case['disk_free_before']=disk_guard(a.out)
    def run(exe,stdin,bin_=None):
        case['steps'][exe]=timed((bin_ or a.grasp_bin)/exe,stdin,out,BUDGET_S[exe],MEM_CAP_BYTES)
    run('rcsfgenerate',generator_stdin(parity,k))
    case.update(**inventory_stream(out/'rcsf.out',parity,'mcdhf'))
    os.link(out/'rcsf.out',out/'rcsf.inp')
    run('rnucleus','8\n0\n0\n0\n0\n0\n');run('rangular','y\n')
    if not previous:run('rwfnestimate',FIRST_ESTIMATE)
    else:
        for i,d in enumerate(previous):os.link(d/'rwfn.out',out/f'prev{i}.w')
        case['estimate_sources']=[d.name for d in previous]
        run('rwfnestimate',step_estimate(len(previous)))
    levels=levels_stdin(block_order(out/'rcsf.out'),fams)
    if varied:
        run('rmcdhf',rmcdhf_stdin(levels,varied),a.strict_bin)
        case['prompts']=prompt_check(out,'method3_all')
        diag=scf_diagnostics((out/'rmcdhf.stdout').read_text(),(out/'rmcdhf.stderr').read_text())
        case['scf']=diag
        if diag['iteration_limit_hit'] or not diag['all_reported_functions_below_ACC']:
            raise ValueError(f'Strict orbital criterion unmet in {name}')
    else:   # fixed-orbital CI over the full s0 list: pinned standard rmcdhf, nothing varied
        run('rmcdhf',rmcdhf_stdin(levels,[]))
        case['prompts']=prompt_check(out,'method3_all')
        if 'performing CI calculation' not in (out/'rmcdhf.stderr').read_text():
            raise ValueError(f'{name}: rmcdhf did not report a fixed-orbital CI')
    for src,ext in (('rcsf.inp','c'),('rwfn.out','w'),('rmix.out','m')):os.link(out/src,out/f'state.{ext}')
    run('jj2lsj','state\nn\nn\ny\n')
    case['leading_families']=identity((out/'state.lsj.lbl').read_text(),parity,fams)
    (a.out/'progress.json').write_text(json.dumps({k_:v for k_,v in receipt.items() if k_!='files'},indent=2)+'\n')
    return out


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--grasp-bin',type=Path,required=True)
    ap.add_argument('--strict-bin',type=Path,required=True)
    ap.add_argument('--v56-receipt',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    a=ap.parse_args()
    if platform.node().lower().split('.')[0]!='sirius':raise RuntimeError('Atomic solvers belong on Sirius')
    for k,v in vars(a).items():
        if isinstance(v,Path):setattr(a,k,v.resolve())
    a.out.mkdir(parents=True,exist_ok=False);shutil.copyfile(__file__,a.out/'driver.py')
    v56=json.loads(a.v56_receipt.read_text())['inputs']
    receipt=dict(schema='rya1354.mcdhf-s0-sequential.v60',grid=dict(RNT=RNT,H=H,N=N,ACCY=ACCY),status='RUNNING',scientific_admission='HOLD',exposure='EXPOSED_NONBLIND',
        adoption_allowed=False,started_utc=datetime.now(timezone.utc).isoformat(),steps={},transitions=None,
        v56_receipt_sha256=sha(a.v56_receipt),nscf=NSCF,mem_cap_bytes=MEM_CAP_BYTES,
        plan={p:[dict(label=l,families=[f'{f}{n}' for f,n in fams],varied=v) for l,fams,v in s] for p,s in STEPS.items()},stopped=None)
    try:
        last={}
        for parity in ('even','odd'):
            done=[]
            for k in range(len(STEPS[parity])):done.insert(0,step(a,receipt,parity,k,list(done)))
            prev=done[0]
            final=receipt['steps'][prev.name]
            if final['blocks']!=v56[parity]['blocks']:raise ValueError(f'final {parity} list differs from the v56 s0 list')
            final['equals_full_s0_list']=True;last[parity]=prev
        receipt['transitions']=transitions(a.grasp_bin,last['even'],last['odd'],a.out/'s0_E1','E1',False)
        receipt['status']='S0_SEQUENTIAL_COMPLETED_REVIEW_REQUIRED'
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
