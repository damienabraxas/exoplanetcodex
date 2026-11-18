"""v57: is the v53-v56 s0 non-convergence a redundant-rotation invariance? (Sirius only)

Hypothesis (from v56): in the Li MCDHF MR without excitations, the CSF space of each J block is
closed under rotations among {1s,2s} (both closed) and within each complete Rydberg manifold
(even: 3p..6p, 4f..5f; odd: 3s..6s, 3d..5d). Those rotations leave every EOL eigenvalue and every
many-electron state unchanged, so an orbital-change criterion cannot converge for them while the
physics is fixed. The strict rmcdhf (no energy stop) then never terminates.

Test: per parity, two runs of the identical s0 problem with the PINNED STANDARD rmcdhf (native
energy stop), method 3 for every orbital (v56 non-default path), 500 cycles, started from
  A: Thomas-Fermi then screened hydrogenic (v53 start)
  B: screened hydrogenic for every orbital
Measured: target eigenvalues (jj2lsj), leading LS families, native orbital residuals, and E1
line strengths (Babushkin and Coulomb) for every even-odd target pair from A and from B.

Declared verdict:
  INVARIANCE_SUPPORTED if, for both parities, both runs stop on the native criterion with the
    targeted families, every orbital whose final residual exceeds ACCY lies in a predicted set,
    every target eigenvalue agrees A vs B within 5e-7 Hartree, and every E1 pair with S >= 1e-6 au
    agrees A vs B within 1e-3 relative in both gauges;
  INVARIANCE_REFUTED if both runs converge with the families but any of the agreement tests fails;
  INCONCLUSIVE otherwise (a run did not stop natively or lost the target families).
Refusals, prompt mismatches or missing files are FAILED_RETAINED. Nothing is adopted.
"""
import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import platform
import re
import shutil

from rya1354_li_scheme_sizing import inventory_stream
from rya1354_mcdhf_eol import BUDGET_S, RWFN_STDIN, generator_stdin, block_order, levels_stdin, identity
from rya1354_mcdhf_s0_controls import rmcdhf_stdin, prompt_check, classify_stop
from rya1354_rci_feasibility import BudgetExceeded, MEM_CAP_BYTES, disk_guard, parse_time, timed
from rya1354_solver_smoke import sha
from rya1354_triplet_correlation import transitions
from rya1354_two_method_summary import scf_diagnostics

STARTS={'A':RWFN_STDIN,'B':'y\n3\n*\n'}
INVARIANT_SETS={'even':[{'1s','2s'},{'3p','4p','5p','6p'},{'4f','5f'}],
                'odd':[{'1s','2s'},{'3s','4s','5s','6s'},{'3d','4d','5d'}]}
E_TOL=5e-7
S_FLOOR=1e-6
S_RTOL=1e-3
LEVEL=re.compile(r'^\s*(\d+)\s+(\d+)\s+([+-])\s+(-?\d+\.\d+)\s+[\d.]+%\s*$')
E1=re.compile(r'^\s*(f[12])\s+(\d+)\s+(\d+)\s+([+-])\s+(f[12])\s+(\d+)\s+(\d+)\s+([+-])\s+([\d.]+)\s+C\s+(\S+)\s+(\S+)\s+(\S+)\s*\n'
              r'\s+B\s+(\S+)\s+(\S+)\s+(\S+)\s*$',re.M)


def in_invariant_set(orbital,parity):
    nl=orbital.rstrip('-')
    return any(nl in s for s in INVARIANT_SETS[parity])


def energies(lbl):
    return {(int(m[2]),i):float(m[4]) for i,m in enumerate(filter(None,(LEVEL.match(l) for l in lbl.splitlines())))}


def line_strengths(table):
    f=lambda x:float(x.replace('D','E'))
    return {(m[1],int(m[2]),int(m[3]),m[5],int(m[6]),int(m[7])):dict(E=float(m[9]),S_C=f(m[12]),S_B=f(m[15]))
            for m in E1.finditer(table)}


def compare(ea,eb,sa,sb):
    if set(ea)!=set(eb):raise ValueError('Target level sets differ between runs')
    if set(sa)!=set(sb) or not sa:raise ValueError('E1 pair sets differ between runs or are empty')
    de=max(abs(ea[k]-eb[k]) for k in ea)
    worst={}
    for g in ('S_B','S_C'):
        rel=[(abs(sa[k][g]-sb[k][g])/max(abs(sa[k][g]),abs(sb[k][g])),k) for k in sa if max(abs(sa[k][g]),abs(sb[k][g]))>=S_FLOOR]
        worst[g]=max(rel)[0] if rel else 0.0
        worst[g+'_pairs']=len(rel)
    ok=de<=E_TOL and worst['S_B']<=S_RTOL and worst['S_C']<=S_RTOL
    return dict(max_abs_dE_hartree=de,max_rel_dS=worst,n_pairs=len(sa),agree=ok)


def run(a,receipt,parity,start):
    name=f'{parity}_{start}';out=a.out/name;out.mkdir();case=dict(parity=parity,start=start,steps={})
    receipt['runs'][name]=case;case['disk_free_before']=disk_guard(a.out)
    def step(exe,stdin,bin_=None):
        case['steps'][exe]=timed((bin_ or a.grasp_bin)/exe,stdin,out,BUDGET_S[exe],MEM_CAP_BYTES)
    step('rcsfgenerate',generator_stdin(parity,0))
    case.update(**inventory_stream(out/'rcsf.out',parity,'mcdhf'))
    if a.v56 and a.v56[parity]['blocks']!=case['blocks']:raise ValueError(f'{parity} s0 list differs from v56')
    os.link(out/'rcsf.out',out/'rcsf.inp')
    step('rnucleus','8\n0\n0\n0\n0\n0\n');step('rangular','y\n');step('rwfnestimate',STARTS[start])
    levels=levels_stdin(block_order(out/'rcsf.out'),parity)
    try:
        step('rmcdhf',rmcdhf_stdin('method3_all',parity,levels))
    except ValueError as exc:
        kind=classify_stop(str(exc),(out/'rmcdhf.stdout').read_text(errors='replace'))
        if kind is None:raise
        case['steps']['rmcdhf']=dict(stop=kind,**parse_time(out/'rmcdhf.time'))
        prompt_check(out,'method3_all');case['outcome']=f'NOT_CONVERGED_{kind.upper()}';return None
    case['prompts']=prompt_check(out,'method3_all')
    diag=scf_diagnostics((out/'rmcdhf.stdout').read_text(),(out/'rmcdhf.stderr').read_text())
    case['scf']=diag
    above=sorted(k for k,v in diag['reported_orbital_self_consistency'].items() if v>diag['native_function_ACC'])
    case['orbitals_above_ACC']=above
    case['above_ACC_all_in_predicted_sets']=all(in_invariant_set(o,parity) for o in above)
    for src,ext in (('rcsf.inp','c'),('rwfn.out','w'),('rmix.out','m')):os.link(out/src,out/f'state.{ext}')
    step('jj2lsj','state\nn\nn\ny\n')
    lbl=(out/'state.lsj.lbl').read_text()
    try:
        case['leading_families']=identity(lbl,parity)
    except ValueError as exc:
        case.update(outcome='CONVERGED_IDENTITY_DIFFERS',identity_error=str(exc));return None
    case['energies']={f'{j}:{i}':e for (j,i),e in energies(lbl).items()}
    case['outcome']='NATIVE_STOP_IDENTITY_OK';return out


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--grasp-bin',type=Path,required=True)
    ap.add_argument('--v56-receipt',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    a=ap.parse_args()
    if platform.node().lower().split('.')[0]!='sirius':raise RuntimeError('Atomic solvers belong on Sirius')
    for k,v in vars(a).items():
        if isinstance(v,Path):setattr(a,k,v.resolve())
    a.out.mkdir(parents=True,exist_ok=False);shutil.copyfile(__file__,a.out/'driver.py')
    a.v56=json.loads(a.v56_receipt.read_text())['inputs']
    receipt=dict(schema='rya1354.mcdhf-s0-invariance.v57',status='RUNNING',scientific_admission='HOLD',exposure='EXPOSED_NONBLIND',
        adoption_allowed=False,started_utc=datetime.now(timezone.utc).isoformat(),runs={},transitions={},comparison=None,verdict=None,
        v56_receipt_sha256=sha(a.v56_receipt),invariant_sets={p:[sorted(s) for s in v] for p,v in INVARIANT_SETS.items()},
        tolerances=dict(E_hartree=E_TOL,S_floor_au=S_FLOOR,S_rel=S_RTOL),mem_cap_bytes=MEM_CAP_BYTES,stopped=None)
    try:
        got={}
        for parity in ('even','odd'):
            for start in STARTS:got[(parity,start)]=run(a,receipt,parity,start)
        if all(got.values()):
            tables={}
            for start in STARTS:
                t=a.out/f'E1_{start}'
                receipt['transitions'][start]=transitions(a.grasp_bin,got[('even',start)],got[('odd',start)],t,'E1',False)
                tables[start]=line_strengths((t/'lower.upper.t').read_text())
            ea={f'{p}:{k}':v for p in ('even','odd') for k,v in receipt['runs'][f'{p}_A']['energies'].items()}
            eb={f'{p}:{k}':v for p in ('even','odd') for k,v in receipt['runs'][f'{p}_B']['energies'].items()}
            cmp_=compare(ea,eb,tables['A'],tables['B'])
            cmp_['above_ACC_all_in_predicted_sets']=all(receipt['runs'][n]['above_ACC_all_in_predicted_sets'] for n in receipt['runs'])
            receipt['comparison']=cmp_
            receipt['verdict']='INVARIANCE_SUPPORTED' if cmp_['agree'] and cmp_['above_ACC_all_in_predicted_sets'] else 'INVARIANCE_REFUTED'
        else:
            receipt['verdict']='INCONCLUSIVE'
        receipt['status']='S0_INVARIANCE_TEST_COMPLETED_REVIEW_REQUIRED'
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
