"""Nonblind n3 valence-correlation test for the primary triplet lifetime failure.

Use the existing ground VV-SD n3 model and independently relaxed triplet seed.
Freeze n<3 functions; optimize n3 only, no observed-energy substitution. This
bounded test does not certify all atom physics or lifetime accuracy.
"""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import platform
import re
import shutil

from rya1354_solver_smoke import run, sha


def require_triplet_quintet_order(triplet,quintet):
    energies=[]
    for text,j,term in ((triplet,1,'3s_3S'),(quintet,2,'3s_5S')):
        match=re.search(rf'^\s*1\s+{j}\s+-\s+(-[\d.]+)[^\n]*\n([^\n]+)',text,re.M)
        if match is None or '2p(3)4S' not in match[2] or not match[2].rstrip().endswith(term):
            raise ValueError('Leading triplet/quintet physical identity absent')
        energies.append(float(match[1]))
    if energies[0]<=energies[1]:
        raise ValueError('Model reverses physical triplet/quintet energy ordering; refuse complete lifetime')


def transitions(bindir, lower, upper, out, operators, same_parity,*,ci=False):
    out.mkdir()
    for tag,source in [('lower',lower),('upper',upper)]:
        for ext in ('c','cm' if ci else 'm','w'):
            shutil.copyfile(source/f'state.{ext}',out/f'{tag}.{ext}')
        shutil.copyfile(source/'state.lsj.lbl',out/f'{tag}.lsj.lbl')
    shutil.copyfile(lower/'isodata',out/'isodata')
    steps=[]
    mixing='y' if ci else 'n'
    for exe,stdin in [('rbiotransform',f'y\n{mixing}\nlower\nupper\ny\n'),
                      ('rtransition',f'y\n{mixing}\nlower\nupper\n{operators}\n'+('n\nn\n' if same_parity else ''))]:
        steps.append(run(bindir/exe,stdin,out,1800))
    if any(s['semantic_errors'] for s in steps) or not (out/('lower.upper.ct' if ci else 'lower.upper.t')).is_file():
        raise ValueError('Native transition refusal or missing table')
    return steps


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--grasp-bin',type=Path,required=True)
    ap.add_argument('--strict-bin',type=Path,required=True)
    ap.add_argument('--ground',type=Path,required=True)
    ap.add_argument('--triplet-seed',type=Path,required=True)
    ap.add_argument('--quintet-seed',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--n',type=int,choices=(3,4),default=3)
    args=ap.parse_args()
    if platform.node().lower().split('.')[0]!='sirius':
        raise RuntimeError('Atomic solvers belong on Sirius')
    for attr in vars(args):
        if isinstance(getattr(args,attr),Path):
            setattr(args,attr,getattr(args,attr).resolve())
    args.out.mkdir(parents=True,exist_ok=False)
    shutil.copyfile(__file__,args.out/'driver.py')
    receipt=dict(schema='rya1354.triplet-correlation.v1',status='RUNNING',scientific_admission='HOLD',
                 started_utc=datetime.now(timezone.utc).isoformat(),steps=[],
                 model=f'n{args.n} VV-SD; fixed n<{args.n} functions; strict orbital stop; no observed-energy substitution')
    upper=args.out/'triplet_3s'
    upper.mkdir()
    try:
        active=','.join(f'{args.n}{l}' for l in 'spdf'[:args.n])
        for exe,stdin in [('rnucleus','8\n0\n0\n0\n0\n0\n'),
                          ('rcsfgenerate',f'*\n0\n1s(2,i)2s(2,i)2p(3,*)3s(1,*)\n\n{active}\n2,2\n2\nn\n')]:
            receipt['steps'].append(run(args.grasp_bin/exe,stdin,upper,600))
        shutil.copyfile(upper/'rcsf.out',upper/'rcsf.inp')
        shutil.copyfile(args.triplet_seed/'state.w',upper/'previous.w')
        for exe,stdin in [('rangular','y\n'),('rwfnestimate','y\n1\nprevious.w\n*\n3\n*\n')]:
            receipt['steps'].append(run(args.grasp_bin/exe,stdin,upper,600))
        spectroscopic='3s' if args.n==3 else ''
        controls=f'n\nn\nn\nn\n1\n{args.n}*\n{spectroscopic}\n500\ny\nn\nn\nn\ny\nn\n-0.9\nn\nn\nn\nn\n1\n'
        receipt['steps'].append(run(args.strict_bin/'rmcdhf',controls,upper,21600))
        if any(s['semantic_errors'] for s in receipt['steps']):
            raise ValueError('Native numerical or scientific refusal')
        for src,ext in [('rcsf.inp','c'),('rwfn.out','w'),('rmix.out','m')]:
            shutil.copyfile(upper/src,upper/f'state.{ext}')
        receipt['steps'].append(run(args.grasp_bin/'jj2lsj','state\nn\nn\ny\n',upper,600))
        labels=(upper/'state.lsj.lbl').read_text()
        if labels.count('3s_3S')!=1 or '2p(3)4S' not in labels:
            raise ValueError('Target 2p3(4S)3s3S identity absent')
        qlabels=(args.quintet_seed/'state.lsj.lbl').read_text()
        if qlabels.count('3s_5S')!=1:
            raise ValueError('Quintet destination identity absent')
        require_triplet_quintet_order(labels,qlabels)
        receipt['steps']+=transitions(args.grasp_bin,args.ground,upper,args.out/'to_ground','E1,M2,E3',False)
        receipt['steps']+=transitions(args.grasp_bin,args.quintet_seed,upper,args.out/'to_5S','M1,E2,M3',True)
        receipt['status']='CORRELATED_MODEL_DECAY_COMPLETED_REVIEW_REQUIRED'
    except Exception as exc:
        receipt.update(status='FAILED_RETAINED',error=f'{type(exc).__name__}: {exc}')
        raise
    finally:
        receipt['files']={str(p.relative_to(args.out)):dict(sha256=sha(p),bytes=p.stat().st_size)
                          for p in sorted(args.out.rglob('*')) if p.is_file()}
        (args.out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')


if __name__=='__main__':
    main()
