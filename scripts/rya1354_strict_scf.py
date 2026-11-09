"""Isolated convergence-only GRASP diagnostic; preserve the pinned build untouched.

Remove only the alternative energy-stationarity stop. The original orbital-change
criterion, Hamiltonian, integrals and iteration updates are unchanged. This is a
declared local solver variant, not an upstream release or an adopted data product.
"""
import argparse
from datetime import datetime, timezone
import difflib
import json
import math
from pathlib import Path
import platform
import shlex
import shutil
import subprocess

from rya1354_solver_smoke import run, sha

ENERGY_STOP='         IF (ABS((WTAEV - WTAEV0)/WTAEV) < 0.001*ACCY) CONVG = .TRUE.'


def strict_source(text):
    if text.count(ENERGY_STOP)!=1:
        raise ValueError('Pinned energy-stop statement is absent or ambiguous')
    return text.replace(ENERGY_STOP,'! RYA-1354 diagnostic: require native orbital-change convergence; no energy-only stop.')


def build(grasp, out):
    out.mkdir(parents=True,exist_ok=False)
    source=grasp/'src/appl/rmcdhf90/scf.f90'
    original=source.read_text()
    strict=strict_source(original)
    (out/'scf.original.f90').write_text(original)
    (out/'scf.f90').write_text(strict)
    (out/'strict_scf.patch').write_text(''.join(difflib.unified_diff(original.splitlines(True),strict.splitlines(True),fromfile='scf.original.f90',tofile='scf.f90')))
    objects=grasp/'build/src/appl/rmcdhf90'
    includes=[objects,*(grasp/'build/modules'/m for m in ('dvd90','9290','mod'))]
    command=['/usr/bin/gfortran','-std=legacy','-fallow-argument-mismatch','-fno-automatic','-O3',
             *('-I'+str(p) for p in includes),'-c',str(out/'scf.f90'),'-o',str(out/'scf.f90.o')]
    with (out/'build.log').open('x') as log:
        subprocess.run(command,cwd=out,stdout=log,stderr=log,check=True)
        link=shlex.split((objects/'CMakeFiles/rmcdhf.dir/link.txt').read_text())
        link=[word for word in link if not word.startswith('-Wl,--dependency-file=')]
        # GNU Fortran's generated link includes object paths relative to objects.
        old='CMakeFiles/rmcdhf.dir/scf.f90.o'
        if link.count(old)!=1 or link.count('-o')!=1:
            raise ValueError('Unexpected pinned link command')
        link[link.index(old)]=str(out/'scf.f90.o')
        link[link.index('-o')+1]=str(out/'rmcdhf')
        subprocess.run(link,cwd=objects,stdout=log,stderr=log,check=True)
    receipt=dict(schema='rya1354.strict-scf-build.v1',upstream_source_sha256=sha(source),
                 original_executable_sha256=sha(grasp/'build/bin/rmcdhf'),
                 compile=command,link=link,change='Remove alternative energy-only convergence; retain native SCNSTY<=ACCY',
                 upstream_revision=subprocess.check_output(['git','rev-parse','HEAD'],cwd=grasp,text=True).strip(),
                 upstream_dirty=subprocess.check_output(['git','status','--porcelain'],cwd=grasp,text=True),
                 files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in out.iterdir() if p.is_file()})
    (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')


def continue_scf(args):
    out=args.out.resolve()
    out.mkdir(parents=True,exist_ok=False)
    shutil.copyfile(__file__,out/'driver.py')
    receipt=dict(schema='rya1354.strict-scf-run.v1',status='RUNNING',scientific_admission='HOLD',
                 started_utc=datetime.now(timezone.utc).isoformat(),seed=str(args.seed),
                 n=args.n,damping=args.damping,cycles=args.cycles,steps=[],
                 numerical_variant='STRICT_ORBITAL_STOP' if args.strict_bin else 'UNMODIFIED_UPSTREAM',
                 state_kind=args.state_kind,
                 grid_spacing=args.grid_spacing,grid_points=args.grid_points,
                 native_ACC=args.grid_spacing**6 if args.grid_spacing is not None else .0625**6,
                 extra_improvements=args.extra_improvements)
    try:
        for name in ('isodata','rcsf.inp',*sorted(p.name for p in args.seed.glob('mcp.*'))):
            shutil.copyfile(args.seed/name,out/name)
        shutil.copyfile(args.seed/'state.w',out/'rwfn.inp')
        # Same EOL five roots and standard weights; only the declared numerical
        # settings differ. Blank integration/node controls retain native defaults.
        roots='1,2\n1\n1,2\n5\n' if args.state_kind=='ground' else '1\n'
        spectroscopic='' if args.state_kind=='ground' or args.n!=3 else '3s'
        grid=('n\nn\nn\nn\n' if args.grid_spacing is None else
              f'n\nn\ny\nn\ny\n{math.exp(-65/16)/8:.17g}\n{args.grid_spacing}\n0\n{args.grid_points}\ny\n{args.grid_spacing**6:.17g}\n')
        controls=(f'{grid}{roots}{args.n}*\n{spectroscopic}\n{args.cycles}\n'
                  f'y\nn\nn\nn\ny\nn\n{-args.damping}\nn\n'
                  +('n\n' if not args.extra_improvements else f'y\n{args.extra_improvements}\n')
                  +'n\nn\n1\n')
        binary=(args.strict_bin or args.grasp_bin)/'rmcdhf'
        step=run(binary,controls,out,21600)
        receipt['steps'].append(step)
        if step['semantic_errors']:
            raise ValueError('Native scientific refusal or iteration-limit exhaustion')
        for src,ext in [('rcsf.inp','c'),('rwfn.out','w'),('rmix.out','m')]:
            shutil.copyfile(out/src,out/f'state.{ext}')
        receipt['steps'].append(run(args.grasp_bin/'jj2lsj','state\nn\nn\ny\n',out,600))
        receipt['status']='NATIVE_STOP_COMPLETED_REVIEW_REQUIRED'
    except Exception as exc:
        receipt.update(status='FAILED_RETAINED',error=f'{type(exc).__name__}: {exc}')
        raise
    finally:
        receipt['files']={str(p.relative_to(out)):dict(sha256=sha(p),bytes=p.stat().st_size)
                          for p in sorted(out.rglob('*')) if p.is_file()}
        (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--build-from',type=Path)
    ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--grasp-bin',type=Path)
    ap.add_argument('--strict-bin',type=Path)
    ap.add_argument('--seed',type=Path)
    ap.add_argument('--n',type=int,choices=(3,4,5))
    ap.add_argument('--cycles',type=int,default=500)
    ap.add_argument('--damping',type=float,default=.9)
    ap.add_argument('--extra-improvements',type=int,default=0)
    ap.add_argument('--state-kind',choices=('ground','quintet3s'),default='ground')
    ap.add_argument('--grid-spacing',type=float)
    ap.add_argument('--grid-points',type=int)
    args=ap.parse_args()
    if platform.node().lower().split('.')[0]!='sirius':
        raise RuntimeError('Solver build/calculations belong on Sirius')
    for attr in ('out','build_from','grasp_bin','strict_bin','seed'):
        value=getattr(args,attr)
        if value is not None:
            setattr(args,attr,value.resolve())
    if args.build_from:
        build(args.build_from,args.out)
    else:
        if not args.seed or not args.grasp_bin or args.n is None or args.cycles<1 or not 0<args.damping<1 or args.extra_improvements<0:
            raise ValueError('Invalid continuation configuration')
        if ((args.grid_spacing is None)!=(args.grid_points is None) or
                (args.grid_spacing is not None and (args.grid_spacing<=0 or not 2<=args.grid_points<=590))):
            raise ValueError('Both valid grid spacing and points are required')
        continue_scf(args)


if __name__=='__main__':
    main()
