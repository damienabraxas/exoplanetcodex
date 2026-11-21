"""v54: declared integration-method controls for the v53 s0 Dirac-Fock EOL stage (Sirius only).

v53 FAILED_RETAINED: with GRASP's default method 1 for spectroscopic orbitals, 4f/5f could not
be solved; the method-2 fallback left 5f (j=7/2) frozen at self-consistency 0.166 and 2s at
1.7e-5 while the EOL energy was stationary. v54 runs, for each parity, the identical s0 problem
(Li 2023 MCDHF MR, 45 even / 36 odd EOL targets, standard weights, all orbitals varied and
spectroscopic, strict rmcdhf, 500 cycles) under three integration-method variants:

  default     GRASP default mode (method 1 for spectroscopic orbitals; v53 settings)
  method3_f   non-default mode, method 3 for the f orbitals only (none in the odd MR)
  method3_all non-default mode, method 3 for every varied orbital (the repo's reference_method3)

Each variant is a MEASURED outcome (converged under the strict criterion or not). Selection rule
(declared): per parity, the first variant in the order above that converges with the targeted
leading LS families is selected. If both parities select, E1 Babushkin/Coulomb rates between the
selected s0 even and odd states are computed. Refusals, prompt mismatches, nonzero exits or
missing files are FAILED_RETAINED. No observed input; nothing adopted.
"""
import argparse
from datetime import datetime, timezone
import json
import os
import re
from pathlib import Path
import platform
import shutil

from rya1354_li_scheme_sizing import inventory_stream
from rya1354_mcdhf_eol import BUDGET_S, RWFN_STDIN, generator_stdin, block_order, levels_stdin, identity, mr_max
from rya1354_rci_feasibility import BudgetExceeded, MEM_CAP_BYTES, disk_guard, parse_time, timed
from rya1354_solver_smoke import sha
from rya1354_triplet_correlation import transitions
from rya1354_two_method_summary import scf_diagnostics

VARIANTS=('default','method3_f','method3_all')
SCHEMA='rya1354.mcdhf-s0-controls.v56'
NSCF=500
PROMPT_ERRORS=('GETRSL: UNABLE','INPUT NOT CORRECT','REDO','MUST EXCEED','VALUE OUT OF RANGE','REPEATED SUBSHELL','SUBSHELL NOT OCCUPIED')


def method3_list(variant,parity):
    if variant=='method3_all':return '*'
    return ','.join(f'{n}f*' for n in range(4,mr_max(parity).get('f',0)+1))


def rmcdhf_stdin(variant,parity,levels):
    if variant=='default':
        return f'y\n{levels}5\n*\n*\n{NSCF}\n'
    m3=method3_list(variant,parity)
    if not m3:return None   # identical to default for this parity
    # RSCFVU default? n; SETDBG debug? n (v54 FAILED_RETAINED: this prompt was missing); GETSCD light/grid? n;
    # ACCY? n; [GETOLD: levels, weights, varied, spectroscopic]; NSCF; modify other defaults? y; node threshold? n;
    # select methods? y; methods 1..4; first-oscillation? n; subshell acceleration? n; eigenvector acceleration? n;
    # NSIC? n; NSOLV? n; Schmidt order? n; SCF orthonormalization order 1 (update order = the default-mode LSORT)
    return (f'n\nn\nn\nn\n{levels}5\n*\n*\n{NSCF}\n'
            f'y\nn\ny\n\n\n{m3}\n\n'
            'n\nn\nn\nn\nn\nn\n1\n')


def classify_stop(message,stdout):
    """'iteration_limit' | 'improv' (measured non-convergence) | None (refusal)."""
    if message.endswith("rc=0 ['MAXIMUM ITERATIONS IN SCF EXCEEDED']"):return 'iteration_limit'
    if (message.endswith("rc=1 ['ERROR STOP']") and 'Error in SUBROUTINE IMPROV' in stdout
            and 'Convergence not obtained' in stdout):return 'improv'
    return None


def prompt_check(out,variant):
    err=(out/'rmcdhf.stderr').read_text(errors='replace')
    bad=[m for m in PROMPT_ERRORS if m in err.upper()]
    if bad:raise ValueError(f'rmcdhf prompt mismatch in {out.name}: {bad}')
    if variant!='default':
        for need in (r'Modify other defaults\?',r'Select a different integration method',r'Method\s+3\s*:',r'Schmidt orthogonalised'):
            if not re.search(need,err):raise ValueError(f'rmcdhf non-default prompt not reached in {out.name}: {need!r}')
    return 'checked'


def run_variant(a,receipt,parity,variant,base):
    name=f'{parity}_{variant}';out=a.out/name
    case=dict(parity=parity,variant=variant,steps={})
    order=block_order(base/'rcsf.out');levels=levels_stdin(order,parity)
    stdin=rmcdhf_stdin(variant,parity,levels)
    if stdin is None:
        case.update(outcome='NOT_APPLICABLE_IDENTICAL_TO_DEFAULT');receipt['runs'][name]=case;return None
    out.mkdir();receipt['runs'][name]=case
    for f in ('rcsf.inp','isodata','rwfn.inp',*sorted(p.name for p in base.glob('mcp.*'))):os.link(base/f,out/f)
    try:
        case['steps']['rmcdhf']=timed(a.strict_bin/'rmcdhf',stdin,out,BUDGET_S['rmcdhf'],MEM_CAP_BYTES)
        limited=False
    except ValueError as exc:
        # Measured non-convergence: a clean (rc=0) iteration-limit stop, or (v56) the native IMPROV stop
        # "Convergence not obtained" after every integration method failed for an orbital (rc=1 ERROR STOP).
        # Anything else is a refusal.
        stdout=(out/'rmcdhf.stdout').read_text(errors='replace')
        kind=classify_stop(str(exc),stdout)
        if kind is None:raise
        improv=kind=='improv'
        case['steps']['rmcdhf']=dict(iteration_limit=not improv,improv_stop=improv,**parse_time(out/'rmcdhf.time'))
        if improv:
            case['prompts']=prompt_check(out,variant)
            failed=re.findall(r'Failure; equation for orbital\s+(\S+)\s+could not be solved using method (\d)',stdout)
            case.update(outcome='NOT_CONVERGED_IMPROV_STOP',unsolved=failed[-2:],
                        iterations=len(re.findall(r'Iteration number\s+\d+',stdout)))
            return None
        limited=True
    case['prompts']=prompt_check(out,variant)
    diag=scf_diagnostics((out/'rmcdhf.stdout').read_text(),(out/'rmcdhf.stderr').read_text())
    case['scf']=diag
    worst=sorted(diag['reported_orbital_self_consistency'].items(),key=lambda kv:-kv[1])[:4]
    case['worst_orbitals']=worst
    if limited or diag['iteration_limit_hit'] or not diag['all_reported_functions_below_ACC']:
        case['outcome']='NOT_CONVERGED';return None
    for src,ext in (('rcsf.inp','c'),('rwfn.out','w'),('rmix.out','m')):os.link(out/src,out/f'state.{ext}')
    case['steps']['jj2lsj']=timed(a.grasp_bin/'jj2lsj','state\nn\nn\ny\n',out,BUDGET_S['jj2lsj'],MEM_CAP_BYTES)
    try:
        case['leading_families']=identity((out/'state.lsj.lbl').read_text(),parity)
    except ValueError as exc:
        case.update(outcome='CONVERGED_IDENTITY_DIFFERS',identity_error=str(exc));return None
    case['outcome']='CONVERGED_IDENTITY_OK';return out


def base_inputs(a,receipt,parity):
    base=a.out/f'{parity}_inputs';base.mkdir();rec=dict(steps={})
    receipt['inputs'][parity]=rec
    rec['disk_free_before']=disk_guard(a.out)
    for exe,stdin in (('rcsfgenerate',generator_stdin(parity,0)),):
        rec['steps'][exe]=timed(a.grasp_bin/exe,stdin,base,BUDGET_S[exe],MEM_CAP_BYTES)
    rec.update(**inventory_stream(base/'rcsf.out',parity,'mcdhf'))
    os.link(base/'rcsf.out',base/'rcsf.inp')
    for exe,stdin in (('rnucleus','8\n0\n0\n0\n0\n0\n'),('rangular','y\n'),('rwfnestimate',RWFN_STDIN)):
        rec['steps'][exe]=timed(a.grasp_bin/exe,stdin,base,BUDGET_S[exe],MEM_CAP_BYTES)
    if a.v53:
        v53=a.v53.get(f's0_{parity}')
        if v53 and v53['blocks']!=rec['blocks']:raise ValueError(f'{parity} s0 list differs from v53')
        rec['v53_blocks_match']=bool(v53)
    return base


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--grasp-bin',type=Path,required=True)
    ap.add_argument('--strict-bin',type=Path,required=True)
    ap.add_argument('--v53-receipt',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    a=ap.parse_args()
    if platform.node().lower().split('.')[0]!='sirius':raise RuntimeError('Atomic solvers belong on Sirius')
    for k,v in vars(a).items():
        if isinstance(v,Path):setattr(a,k,v.resolve())
    a.out.mkdir(parents=True,exist_ok=False)
    shutil.copyfile(__file__,a.out/'driver.py')
    a.v53=json.loads(a.v53_receipt.read_text())['stages']
    receipt=dict(schema=SCHEMA,status='RUNNING',scientific_admission='HOLD',exposure='EXPOSED_NONBLIND',
        adoption_allowed=False,started_utc=datetime.now(timezone.utc).isoformat(),inputs={},runs={},selected={},transitions=None,
        v53_receipt_sha256=sha(a.v53_receipt),variants=VARIANTS,nscf=NSCF,mem_cap_bytes=MEM_CAP_BYTES,
        selection_rule='per parity, first variant in declared order that converges (strict) with the targeted leading LS families',
        stopped=None)
    try:
        for parity in ('even','odd'):
            base=base_inputs(a,receipt,parity)
            for variant in VARIANTS:
                ok=run_variant(a,receipt,parity,variant,base)
                if ok and parity not in receipt['selected']:receipt['selected'][parity]=dict(variant=variant,path=ok.name)
            (a.out/'progress.json').write_text(json.dumps(receipt,indent=2)+'\n')
        if set(receipt['selected'])=={'even','odd'}:
            ev,od=(a.out/receipt['selected'][p]['path'] for p in ('even','odd'))
            receipt['transitions']=transitions(a.grasp_bin,ev,od,a.out/'s0_E1',f'E1',False)
            receipt['status']='S0_CONTROLS_COMPLETED_SELECTED_REVIEW_REQUIRED'
        else:
            receipt['status']='S0_CONTROLS_COMPLETED_NO_SELECTION'
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
