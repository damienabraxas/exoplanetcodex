"""Task-local pinned FAC build and upstream Fe example; no O prediction."""
import argparse
from datetime import datetime,timezone
import hashlib
import json
import math
from pathlib import Path
import platform
import shutil
import subprocess

PIN='9ab6defb6af62e1bc8251a35ff995809abef717c'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def level_inventory(text):
    """Read verbose FAC level records without importing its compiled interface."""
    rows=[]
    for line in text.splitlines():
        f=line.split()
        if len(f)<9 or not f[0].isdigit():
            continue
        try:
            index,ibase=int(f[0]),int(f[1]);energy=float(f[2]);parity,vnl,j=map(int,f[3:6])
        except ValueError:
            continue
        if parity not in (0,1) or j<0 or not math.isfinite(energy):
            raise ValueError('Invalid FAC level quantum numbers')
        rows.append(dict(index=index,ibase=ibase,energy_eV=energy,parity=parity,vnl=vnl,two_J=j,
                         labels=f[6:]))
    if not rows or len({r['index'] for r in rows})!=len(rows):
        raise ValueError('Missing/duplicate FAC level inventory')
    return rows


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args()
    if platform.node().lower().split('.')[0]!='sirius':
        raise RuntimeError('Atomic compilation/calculations belong on Sirius')
    args.source=args.source.resolve();args.out=args.out.resolve()
    head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=args.source,text=True).strip()
    if head!=PIN or subprocess.check_output(['git','diff','--name-only'],cwd=args.source,text=True).strip():
        raise ValueError('FAC pinned tracked source differs')
    args.out.mkdir(parents=True,exist_ok=False)
    shutil.copyfile(__file__,args.out/'driver.py')
    for name in ('README.md','LICENSE','configure.ac','doc/manual.tex','demo/README','demo/structure/fe17_structure.sf','demo/structure/ref/ne.lev','demo/structure/ref/ne.tr'):
        dest=args.out/'source_documents'/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(args.source/name,dest)
    receipt=dict(schema='rya1354.fac-bootstrap.v17',status='RUNNING',pin=PIN,
        source_url='https://github.com/flexible-atomic-code/fac',started_utc=datetime.now(timezone.utc).isoformat(),
        steps=[],scientific_admission='HOLD',adoption_allowed=False,accepted_lifetime=None,
        uncertainty=None,covariance=None,calibrated_sigma=None,
        scope='Task-local build and upstream Fe example execution/inventory only. No neutral-O calculation, independent oxygen benchmark, physical accuracy or grade admission.')
    def execute(name,command,cwd,timeout):
        with (args.out/(name+'.stdout')).open('xb') as stdout,(args.out/(name+'.stderr')).open('xb') as stderr:
            result=subprocess.run(command,cwd=cwd,stdout=stdout,stderr=stderr,timeout=timeout,check=False)
        receipt['steps'].append(dict(name=name,command=command,cwd=str(cwd),returncode=result.returncode))
        if result.returncode:
            raise RuntimeError(name+': exit '+str(result.returncode))
    try:
        execute('gcc_version',['gcc','--version'],args.source,30)
        execute('gfortran_version',['gfortran','--version'],args.source,30)
        execute('configure',['./configure','--prefix='+str(args.out/'unused_install_prefix'),
                '--enable-copt=-O2 -std=gnu17','--enable-fopt=-O2 -fallow-argument-mismatch',
                'CC=gcc','FC=gfortran','F77=gfortran'],args.source,180)
        execute('make',['make','-j1'],args.source,1800)
        binary=args.source/'sfac/sfac'
        if not binary.is_file():
            raise ValueError('FAC executable absent after build')
        receipt['binary']=dict(path=str(binary),sha256=sha(binary),bytes=binary.stat().st_size)
        smoke=args.out/'upstream_structure';smoke.mkdir()
        shutil.copyfile(args.source/'demo/structure/fe17_structure.sf',smoke/'input.sf')
        execute('upstream_example',[str(binary),'input.sf'],smoke,600)
        actual=level_inventory((smoke/'ne.lev').read_text())
        expected=level_inventory((args.source/'demo/structure/ref/ne.lev').read_text())
        fields=('index','ibase','parity','vnl','two_J','labels')
        if len(actual)!=len(expected) or any(any(a[k]!=b[k] for k in fields) for a,b in zip(actual,expected)):
            raise ValueError('Upstream example level identity inventory differs')
        receipt['upstream_level_inventory']=dict(levels=len(actual),exact_quantum_inventory=True,
            max_absolute_energy_difference_eV=max(abs(a['energy_eV']-b['energy_eV']) for a,b in zip(actual,expected)),
            numerical_accuracy_acceptance=False,
            note='Difference reported without invented tolerance. This is not the full PFAC-dependent upstream checker or an oxygen accuracy test.')
        if not (smoke/'ne.tr').is_file() or not (smoke/'ne.tr.b').stat().st_size:
            raise ValueError('Upstream radiative output missing')
        receipt['status']='FAC_BUILD_AND_UPSTREAM_INVENTORY_COMPLETED_NOT_O_VALIDATED'
    except Exception as exc:
        receipt.update(status='FAILED_RETAINED',error=type(exc).__name__+': '+str(exc))
        raise
    finally:
        receipt['tracked_source_diff']=subprocess.check_output(['git','diff','--name-only'],cwd=args.source,text=True).strip()
        receipt['files']={str(p.relative_to(args.out)):dict(sha256=sha(p),bytes=p.stat().st_size)
                          for p in sorted(args.out.rglob('*')) if p.is_file()}
        (args.out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')


if __name__=='__main__':
    main()
