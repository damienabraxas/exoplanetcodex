"""v48-v52: RCI feasibility levers for the Li et al. 2023 (arXiv 2304.02310, Table 1) O I scheme on
Sirius. CSF lists and COST only: no physical orbitals, no adopted energy, rate or grade.

Section A (reduction): for n_max 6 and 7, both parities, regenerate the v47 "MR in RCI" SD list
(sha256 must equal the v47 receipt), generate the MR itself (0 excitations, same active set and
2J range) and run native rcsfinteract against it with the Dirac-Coulomb (1) and the
Dirac-Coulomb-Breit (2) reduction. Every reduced list is inventoried, must be a subset of the
full list, and must reproduce rcsfinteract's own per-block table.

Section B (cost): odd J=1 block (the 2p3 3s 3S1 upper level) of the same scheme, DCB-reduced, at
n_max 4, 5, 6: rnucleus (point static Z=8), rangular, rwfnestimate (Thomas-Fermi, NOT physical
orbitals), rci (Breit static limit, lowest state), each under GNU time -v. Energies printed by rci
from Thomas-Fermi orbitals are NOT results; only wall time, peak RSS and disk are.

A step exceeding its declared budget (or the disk guard) is a MEASURED OUTCOME: it is recorded
and ends its section. Any refusal, nonzero exit or failed cross-check stops the run
(FAILED_RETAINED). The receipt carries a files inventory (sha256, bytes) of every output.
"""
import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import shutil
import signal
import subprocess
import time

from rya1354_li_scheme_sizing import TWOJ, active, inventory_stream, references
from rya1354_solver_smoke import sha

BUDGET_S={'rcsfgenerate':3600,'rcsfinteract':7200,'rnucleus':600,'rangular':14400,'rwfnestimate':1800,'rci':21600}
DISK_GUARD_BYTES=40*1024**3
MARKERS=('TOO MANY ELECTRONS FOR THIS CODE','CONFIGURATION DATA WRONG','PROGRAM HALTED',
         'MAXIMUM ITERATIONS IN SCF EXCEEDED','ERROR STOP','DIFFERENT CLOSE SHELLS',
         'NOT A CONFIGURATION SYMMETRY LIST','ERROR WHEN OPENING','DIFFERENT ORDER OF PEEL ORBITALS')
REDUCTION_LAYERS=(6,7)
COST_LAYERS=(4,5,6)
COST_TWOJ='2,2'
ESTIMATE_OUTPUT='rwfn.inp'   # v50 FAILED_RETAINED: rwfnestimate writes rwfn.inp; rwfn.out is rmcdhf's
RCI_BREIT_STATIC='y\ny\n1.d-6\nn\nn\nn\nn\n'   # transverse on, frequencies scaled to the static limit; no VP/NMS/SMS/SE


class BudgetExceeded(Exception):
    pass


class MemoryCapReached(BudgetExceeded):
    """v52: a section-B step stopped at the declared address-space cap (a measured outcome)."""


MEM_CAP_BYTES=12*1024**3   # prlimit --as on section-B native steps (15 GB host, shared)
MEMORY_MARKERS=('UNABLE TO ALLOCATE','ALLOCATION WOULD EXCEED MEMORY LIMIT','CANNOT ALLOCATE MEMORY','OUT OF MEMORY')
RWFN_STDIN='y\n2\n*\n3\n*\n'   # v51 FAILED_RETAINED: Thomas-Fermi cannot build 5g; remaining subshells -> screened hydrogenic


def generator_stdin(stage,parity,n_max,excitations,twoj=None):
    kept,_=references(stage,parity,n_max)
    return '*\n0\n'+'\n'.join(kept)+'\n\n'+active(n_max)+'\n'+(twoj or TWOJ[parity])+f'\n{excitations}\nn\n'


def parse_time(path):
    text=path.read_text()
    rss=re.search(r'Maximum resident set size \(kbytes\): (\d+)',text)
    return dict(max_rss_kib=int(rss[1]) if rss else None)


def timed(executable,stdin,cwd,budget,mem_cap=None):
    """One native step under GNU time -v in its own process group; killed as a group at budget."""
    name=executable.name;cwd=Path(cwd).resolve()   # v48: a relative -o path resolved against the child's cwd
    if not executable.is_file():raise FileNotFoundError(executable)
    (cwd/f'{name}.stdin').write_text(stdin)
    t0=time.monotonic()
    with (cwd/f'{name}.stdout').open('xb') as out,(cwd/f'{name}.stderr').open('xb') as err:
        cap=['/usr/bin/prlimit',f'--as={mem_cap}','--'] if mem_cap else []
        p=subprocess.Popen(['/usr/bin/time','-v','-o',str(cwd/f'{name}.time')]+cap+[str(executable)],
                           stdin=subprocess.PIPE,stdout=out,stderr=err,cwd=cwd,start_new_session=True)
        try:
            p.communicate(stdin.encode(),timeout=budget)
        except subprocess.TimeoutExpired:
            os.killpg(p.pid,signal.SIGKILL);p.wait()
            raise BudgetExceeded(f'{name} exceeded {budget} s in {cwd.name}')
    record=dict(executable=str(executable),sha256=sha(executable),returncode=p.returncode,
                wall_s=round(time.monotonic()-t0,3),budget_s=budget,mem_cap_bytes=mem_cap)
    record.update(parse_time(cwd/f'{name}.time'))
    output=((cwd/f'{name}.stdout').read_text(errors='replace')+(cwd/f'{name}.stderr').read_text(errors='replace')).upper()
    record['semantic_errors']=[m for m in MARKERS if m in output]
    memory=[m for m in MEMORY_MARKERS if m in output]
    if mem_cap and memory:
        raise MemoryCapReached(f'{name} stopped at the {mem_cap} B address-space cap in {cwd.name}: {memory} rc={p.returncode} rss={record["max_rss_kib"]} KiB')
    if p.returncode or record['semantic_errors']:
        raise ValueError(f'Native refusal in {cwd.name}/{name}: rc={p.returncode} {record["semantic_errors"]}')
    return record


def records(path):
    """Native 3-line CSF records, trailing blanks stripped (rcsfinteract TRIMs); block separators skipped."""
    with open(path) as f:
        for _ in range(5):next(f)
        buf=[]
        for line in f:
            if line.startswith('*') or line.startswith(' *'):continue
            buf.append(line.rstrip())
            if len(buf)==3:yield '\n'.join(buf);buf=[]
    if buf:raise ValueError('Incomplete native CSF records')


def digest_set(path):
    return {hashlib.sha1(r.encode()).digest() for r in records(path)}


def peel_header(path):
    with open(path) as f:return [next(f) for _ in range(4)][3].split()


def mr_input(full,mr,dest):
    """v49 FAILED_RETAINED: rcsfgenerate with 0 excitations lists only the orbitals the MR occupies,
    so its peel header is not a prefix of the full list's (rcsfinteract set_CSF_list requires that).
    rcsfinteract parses MR records against the FULL list's orbital table (PRSRSL on rcsf.inp); the MR
    header is used only for the prefix check. Derive rcsfmr.inp = full header + the native MR body,
    byte-identical, after checking the MR orbitals are an order-preserving subset of the full header."""
    with open(full) as f:fh=[next(f) for _ in range(4)]
    text=Path(mr).read_text();lines=text.split('\n')
    if lines[:3]!=[x.rstrip('\n') for x in fh[:3]]:raise ValueError('MR and full list header records differ')
    mo,fo=lines[3].split(),fh[3].split();it=iter(fo)
    if not all(o in it for o in mo):raise ValueError('MR orbitals are not an ordered subset of the full peel list')
    body='\n'.join(lines[4:])
    Path(dest).write_text(''.join(fh)+body)
    if Path(dest).read_text().split('\n')[4:]!=lines[4:]:raise ValueError('Derived MR body differs from native MR')
    return dict(native_mr_orbitals=mo,header_from='full rcsf.out',body_sha256=sha_text(body))


def sha_text(text):
    return hashlib.sha256(text.encode()).hexdigest()


def interact_table(stdout):
    rows=re.findall(r'^\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s*$',stdout,re.M)
    if not rows:raise ValueError('rcsfinteract block table absent')
    return [dict(block=int(b),MR=int(m),before=int(x),after=int(y)) for b,m,x,y in rows]


def disk_guard(path):
    free=shutil.disk_usage(path).free
    if free<DISK_GUARD_BYTES:raise BudgetExceeded(f'disk guard: {free} B free < {DISK_GUARD_BYTES}')
    return free


def dir_bytes(path,prefix=''):
    return sum(p.stat().st_size for p in path.iterdir() if p.is_file() and p.name.startswith(prefix))


def generate(bindir,out,stage,parity,n_max,excitations,twoj=None):
    out.mkdir()
    step=timed(bindir/'rcsfgenerate',generator_stdin(stage,parity,n_max,excitations,twoj),out,BUDGET_S['rcsfgenerate'])
    return step


def reduce(bindir,base,parity,n_max,mode,full,mr):
    """rcsfinteract in its own directory; inputs hard-linked (no copy, identical inode)."""
    out=base/mode;out.mkdir()
    os.link(full,out/'rcsf.inp');derivation=mr_input(full,mr,out/'rcsfmr.inp')
    step=timed(bindir/'rcsfinteract',{'DC':'1\n','DCB':'2\n'}[mode],out,BUDGET_S['rcsfinteract'])
    table=interact_table((out/'rcsfinteract.stdout').read_text())
    inv=inventory_stream(out/'rcsf.out',parity,'rci')
    step['rcsfmr_derivation']=derivation
    return step,table,inv


def check_reduction(table,inv_full,inv_red,inv_mr,full,red,mr):
    labels=sorted(inv_full['blocks'],key=lambda s:(s[-1],int(s[:-1])))
    if [t['before'] for t in table]!=[inv_full['blocks'][k] for k in labels]:
        raise ValueError('rcsfinteract Before counts differ from the full-list inventory')
    if [t['after'] for t in table]!=[inv_red['blocks'].get(k,0) for k in labels]:
        raise ValueError('rcsfinteract After counts differ from the reduced-list inventory')
    if [t['MR'] for t in table]!=[inv_mr['blocks'].get(k,0) for k in labels]:
        raise ValueError('rcsfinteract MR counts differ from the MR inventory')
    red_set=digest_set(red)
    if not red_set<=digest_set(full):raise ValueError('Reduced list is not a subset of the full list')
    if not digest_set(mr)<=red_set:raise ValueError('Reduced list omits an MR CSF')
    return labels


def section_reduction(a,receipt,v47):
    for n_max in REDUCTION_LAYERS:
        for parity in ('even','odd'):
            name=f'reduce_{parity}_n{n_max}';base=a.out/name;base.mkdir()
            case=receipt['cases'][name]=dict(section='A',references=references('rci',parity,n_max)[0],
                active=active(n_max),twoJ=TWOJ[parity],steps={})
            case['disk_free_before']=disk_guard(a.out)
            case['steps']['full']=generate(a.grasp_bin,base/'full','rci',parity,n_max,2)
            full=base/'full'/'rcsf.out'
            case['full_sha256']=sha(full);case['full_bytes']=full.stat().st_size
            ref=v47['cases'][f'rci_{parity}_n{n_max}']
            if (case['full_sha256'],case['full_bytes'])!=(ref['rcsf_sha256'],ref['rcsf_bytes']):
                raise ValueError('Regenerated full list differs from the v47 receipt: '+name)
            case['steps']['mr']=generate(a.grasp_bin,base/'mr','rci',parity,n_max,0)
            mr=base/'mr'/'rcsf.out'
            inv_full=dict(CSFs=sum(ref['blocks'].values()),blocks=ref['blocks'])
            inv_mr=inventory_stream(mr,parity,'rci');case['MR']=inv_mr
            case['full']=inv_full
            for mode in ('DC','DCB'):
                step,table,inv=reduce(a.grasp_bin,base,parity,n_max,mode,full,mr)
                labels=check_reduction(table,inv_full,inv,inv_mr,full,base/mode/'rcsf.out',mr)
                case['steps'][mode]=step
                case[mode]=dict(table=table,labels=labels,**inv,rcsf_sha256=sha(base/mode/'rcsf.out'),
                    rcsf_bytes=(base/mode/'rcsf.out').stat().st_size,
                    kept_fraction={k:inv['blocks'].get(k,0)/inv_full['blocks'][k] for k in labels})
            (a.out/'progress.json').write_text(json.dumps(receipt,indent=2)+'\n')


def section_cost(a,receipt,section_a=None):
    for n_max in COST_LAYERS:
        name=f'cost_odd_J1_n{n_max}';base=a.out/name;base.mkdir()
        case=receipt['cases'][name]=dict(section='B',references=references('rci','odd',n_max)[0],
            active=active(n_max),twoJ=COST_TWOJ,steps={},
            scope='COST_ONLY: Thomas-Fermi/screened-hydrogenic orbitals; rci energies are not results')
        case['disk_free_before']=disk_guard(a.out)
        case['steps']['full']=generate(a.grasp_bin,base/'full','rci','odd',n_max,2,COST_TWOJ)
        case['steps']['mr']=generate(a.grasp_bin,base/'mr','rci','odd',n_max,0,COST_TWOJ)
        full,mr=base/'full'/'rcsf.out',base/'mr'/'rcsf.out'
        inv_full=inventory_stream(full,'odd','rci');inv_mr=inventory_stream(mr,'odd','rci')
        step,table,inv=reduce(a.grasp_bin,base,'odd',n_max,'DCB',full,mr)
        case['steps']['DCB']=step
        check_reduction(table,inv_full,inv,inv_mr,full,base/'DCB'/'rcsf.out',mr)
        case.update(full=inv_full,MR=inv_mr,DCB=dict(table=table,**inv))
        if n_max in REDUCTION_LAYERS:
            other=receipt['cases'].get(f'reduce_odd_n{n_max}') or (section_a or {}).get('cases',{}).get(f'reduce_odd_n{n_max}')
            if other and other['DCB']['blocks'].get('1-')!=inv['blocks'].get('1-'):
                raise ValueError('J=1 block count differs from the section A reduction')
        run=base/'run';run.mkdir()
        os.link(base/'DCB'/'rcsf.out',run/'rcsf.inp')
        try:
            cap=MEM_CAP_BYTES
            case['steps']['rnucleus']=timed(a.grasp_bin/'rnucleus','8\n0\n0\n0\n0\n0\n',run,BUDGET_S['rnucleus'],cap)
            case['steps']['rwfnestimate']=timed(a.grasp_bin/'rwfnestimate',RWFN_STDIN,run,BUDGET_S['rwfnestimate'],cap)
            case['steps']['rangular']=timed(a.grasp_bin/'rangular','y\n',run,BUDGET_S['rangular'],cap)
            case['angular_bytes']=dir_bytes(run,'mcp.');case['disk_free_after_rangular']=disk_guard(a.out)
            os.link(run/'rcsf.inp',run/'state.c');os.link(run/ESTIMATE_OUTPUT,run/'state.w')
            case['steps']['rci']=timed(a.grasp_bin/'rci','y\nstate\n'+RCI_BREIT_STATIC+'1\n',run,BUDGET_S['rci'],cap)
            case['run_bytes']=dir_bytes(run)
        except BudgetExceeded as exc:
            case['budget_stop']=str(exc);receipt['budget_stops'].append(str(exc))
            (a.out/'progress.json').write_text(json.dumps(receipt,indent=2)+'\n')
            return
        (a.out/'progress.json').write_text(json.dumps(receipt,indent=2)+'\n')


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--grasp-bin',type=Path,required=True);ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--v47-receipt',type=Path,required=True)
    ap.add_argument('--sections',choices=('AB','B'),default='AB')
    ap.add_argument('--section-a-receipt',type=Path,help='completed section A receipt (v51) for the B cross-check');a=ap.parse_args()
    section_a=None
    if a.sections=='B':
        if not a.section_a_receipt:raise ValueError('Section B alone needs the section A receipt')
        section_a=json.loads(a.section_a_receipt.read_text())
        if sorted(k for k in section_a['cases'] if k.startswith('reduce_'))!=sorted(f'reduce_{p}_n{n}' for p in ('even','odd') for n in REDUCTION_LAYERS):
            raise ValueError('Section A receipt is incomplete')
    if platform.node().lower().split('.')[0]!='sirius':raise RuntimeError('Atomic native runs belong on Sirius')
    a.grasp_bin=a.grasp_bin.resolve();a.out=a.out.resolve();a.out.mkdir(parents=True,exist_ok=False)
    v47=json.loads(a.v47_receipt.read_text())
    if v47.get('status')!='SIZING_COMPLETED_NO_ATOMIC_PREDICTION':raise ValueError('v47 receipt is not the completed sizing')
    shutil.copyfile(__file__,a.out/'driver.py')
    receipt=dict(schema='rya1354.rci-feasibility.v52',status='RUNNING',scientific_admission='HOLD',
        started_utc=datetime.now(timezone.utc).isoformat(),cases={},budget_stops=[],
        v47_receipt_sha256=sha(a.v47_receipt),budgets_s=BUDGET_S,disk_guard_bytes=DISK_GUARD_BYTES,
        source='Li et al. 2023 arXiv2304.02310 Table 1 MR in RCI; rcsfinteract (Gaigalas & Froese Fischer); no fitted or observed inputs',
        scope='CSF_REDUCTION_AND_COST_ONLY_NO_PHYSICAL_ORBITALS_ENERGY_RATE_OR_GRADE')
    try:
        if a.sections=='AB':
            try:section_reduction(a,receipt,v47)
            except BudgetExceeded as exc:receipt['budget_stops'].append('A: '+str(exc))
        else:
            receipt.update(sections='B',section_a_receipt_sha256=sha(a.section_a_receipt))
        section_cost(a,receipt,section_a)
        receipt['status']='FEASIBILITY_COMPLETED_NO_ATOMIC_PREDICTION'
    except Exception as exc:
        receipt.update(status='FAILED_RETAINED',error=f'{type(exc).__name__}: {exc}');raise
    finally:
        receipt['finished_utc']=datetime.now(timezone.utc).isoformat()
        seen={};files={}
        for p in sorted(a.out.rglob('*')):
            if not p.is_file() or p.name=='receipt.json':continue
            st=p.stat();key=(st.st_dev,st.st_ino)
            if key not in seen:seen[key]=dict(sha256=sha(p),bytes=st.st_size)
            files[str(p.relative_to(a.out))]=dict(seen[key])
        receipt['files']=files
        (a.out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')


if __name__=='__main__':main()
