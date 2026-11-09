"""Rebuild balanced O I lifetime states from strictly qualified reference orbitals.

State-specific reference configurations, then valence SD n3/n4/n5 with frozen
earlier layers. All promotions require the native orbital criterion and correct
leading LS/parent/J/parity. A fresh NSIC=0 restart verifies each optimized state.
No observed energy, gf or lifetime is an optimization input.
"""
import argparse
from datetime import datetime,timezone
import json
import math
from pathlib import Path
import platform
import re
import shutil

from rya1354_solver_smoke import run,sha
from rya1354_two_method_summary import scf_diagnostics
from rya1354_triplet_correlation import transitions,require_triplet_quintet_order
from rya1354_radial_controls import initialization_input

STATES={
    'ground':dict(config='1s(2,i)2s(2,i)2p(4,{occupation})',js='0,4',roots='1,2\n1\n1,2\n5\n'),
    'quintet':dict(config='1s(2,i)2s(2,i)2p(3,{occupation})3s(1,{occupation})',js='4,4',roots='1\n'),
    'triplet':dict(config='1s(2,i)2s(2,i)2p(3,{occupation})3s(1,{occupation})',js='2,2',roots='1\n')}


def identity(text,state):
    blocks=re.findall(r'^\s*\d+\s+(\d+)\s+([+-])\s+[^\n]+\n([^\n]+)',text,re.M)
    expected=([('0','+','_3P'),('0','+','_1S'),('1','+','_3P'),('2','+','_3P'),('2','+','_1D')]
              if state=='ground' else [('2' if state=='quintet' else '1','-','3s_5S' if state=='quintet' else '3s_3S')])
    if len(blocks)!=len(expected) or any(j!=ej or p!=ep or not line.rstrip().endswith(term)
                                       or ('2p(4)' if state=='ground' else '2p(3)4S') not in line
                                       for (j,p,line),(ej,ep,term) in zip(blocks,expected)):
        raise ValueError('Leading physical state identity changed: '+state)


def optimize(args,out,state,n,seed,receipt,*,verify=False):
    out.mkdir()
    if n==0 or verify:
        for name in ('isodata','rcsf.inp',*sorted(p.name for p in seed.glob('mcp.*'))):
            shutil.copyfile(seed/name,out/name)
        shutil.copyfile(seed/'state.w',out/'rwfn.inp')
    else:
        spec=STATES[state]
        active=','.join(f'{n}{l}' for l in 'spdfg'[:n])
        config=spec['config'].format(occupation='*')
        for exe,stdin in [('rnucleus','8\n0\n0\n0\n0\n0\n'),
                          ('rcsfgenerate',f'*\n0\n{config}\n\n{active}\n{spec["js"]}\n2\nn\n')]:
            receipt['steps'].append(run(args.grasp_bin/exe,stdin,out,600))
        shutil.copyfile(out/'rcsf.out',out/'rcsf.inp')
        shutil.copyfile(seed/'state.w',out/'previous.w')
        estimate=initialization_input(3,439,.03125) if args.half_grid else 'y\n1\nprevious.w\n*\n3\n*\n'
        for exe,stdin in [('rangular','y\n'),('rwfnestimate',estimate)]:
            receipt['steps'].append(run(args.grasp_bin/exe,stdin,out,1800))
    varied='*' if n==0 else f'{n}*'
    spectroscopic='*' if n==0 else '3s' if n==3 and state!='ground' else ''
    damping=.5 if n==5 else .9
    nsic=0 if verify or n==0 else 20
    grid=(f'n\nn\ny\nn\ny\n{math.exp(-65/16)/8:.17g}\n0.03125\n0\n439\ny\n{.0625**6:.17g}\n'
          if args.half_grid else 'n\nn\nn\nn\n')
    controls=(f'{grid}{STATES[state]["roots"]}{varied}\n{spectroscopic}\n500\n'
              +'y\nn\n'+('y\n\n\n*\n\n' if n==0 and args.reference_method3 else 'n\n')
              +f'n\ny\nn\n{-damping}\nn\n'+('n\n' if nsic==0 else f'y\n{nsic}\n')+'n\nn\n1\n')
    step=run(args.strict_bin/'rmcdhf',controls,out,43200)
    receipt['steps'].append(step)
    if step['semantic_errors']:
        raise ValueError('Native refusal or SCF iteration limit; no layer promotion')
    diagnostic=scf_diagnostics((out/'rmcdhf.stdout').read_text(),(out/'rmcdhf.stderr').read_text())
    if diagnostic['iteration_limit_hit'] or not diagnostic['all_reported_functions_below_ACC']:
        raise ValueError('Native orbital criterion unmet; no layer promotion')
    for source,ext in [('rcsf.inp','c'),('rwfn.out','w'),('rmix.out','m')]:
        shutil.copyfile(out/source,out/f'state.{ext}')
    receipt['steps'].append(run(args.grasp_bin/'jj2lsj','state\nn\nn\ny\n',out,600))
    identity((out/'state.lsj.lbl').read_text(),state)
    receipt['qualified_states'].append(dict(path=str(out.relative_to(args.out)),state=state,n=n,
                                            fresh_restart=verify,diagnostics=diagnostic,
                                            seed_state_w_sha256=sha(seed/'state.w'),
                                            frozen_layers='NONE' if n==0 else f'n<{n}',damping=damping,NSIC=nsic))


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--grasp-bin',type=Path,required=True)
    ap.add_argument('--strict-bin',type=Path,required=True)
    ap.add_argument('--ground-seed',type=Path,required=True)
    ap.add_argument('--triplet-seed',type=Path,required=True)
    ap.add_argument('--quintet-seed',type=Path,required=True)
    ap.add_argument('--n-max',type=int,choices=(3,4,5),default=5)
    ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--reference-method3',action='store_true',help='Explicit native integration-method control for reference functions; retain node constraints')
    ap.add_argument('--half-grid',action='store_true',help='H0.03125/N439 for all stages, same explicit original ACCY; preserve endpoint')
    args=ap.parse_args()
    if platform.node().lower().split('.')[0]!='sirius':
        raise RuntimeError('Atomic solvers belong on Sirius')
    for attr,value in vars(args).items():
        if isinstance(value,Path):
            setattr(args,attr,value.resolve())
    args.out.mkdir(parents=True,exist_ok=False)
    shutil.copyfile(__file__,args.out/'driver.py')
    receipt=dict(schema='rya1354.balanced-lifetime-ladder.v1',status='RUNNING',scientific_admission='HOLD',
                 started_utc=datetime.now(timezone.utc).isoformat(),steps=[],qualified_states=[],
                 reference_method3=args.reference_method3,
                 half_grid=args.half_grid,native_ACC=.0625**6,
                 model='Point stationary Z8/DC. Strict-qualified all-orbital reference states; VV-SD n3..nmax, closed 1s2/2s2 occupations, frozen earlier radial layers. Independent state-specific cores; no observed-energy substitution or fitted strength.')
    seeds={s:getattr(args,s+'_seed') for s in STATES}
    try:
        for n in (0,*range(3,args.n_max+1)):
            for state in STATES:
                out=args.out/f'n{n}_{state}'
                optimize(args,out,state,n,seeds[state],receipt)
                verified=args.out/f'n{n}_{state}_verify'
                optimize(args,verified,state,n,out,receipt,verify=True)
                seeds[state]=verified
            require_triplet_quintet_order((seeds['triplet']/'state.lsj.lbl').read_text(),(seeds['quintet']/'state.lsj.lbl').read_text())
            for name,lower,operators,same in [('to_ground',seeds['ground'],'E1,M2,E3',False),
                                              ('to_5S',seeds['quintet'],'M1,E2,M3',True)]:
                receipt['steps']+=transitions(args.grasp_bin,lower,seeds['triplet'],args.out/f'n{n}_{name}',operators,same)
            receipt['last_completed_layer']=n
            # Checkpoint metadata only. Final receipt always rehashes all files.
            (args.out/'progress.json').write_text(json.dumps({k:v for k,v in receipt.items() if k!='files'},indent=2)+'\n')
        receipt['status']='BALANCED_LIFETIME_LADDER_COMPLETED_REVIEW_REQUIRED'
    except Exception as exc:
        receipt.update(status='FAILED_RETAINED',error=f'{type(exc).__name__}: {exc}')
        raise
    finally:
        receipt['files']={str(p.relative_to(args.out)):dict(sha256=sha(p),bytes=p.stat().st_size)
                          for p in sorted(args.out.rglob('*')) if p.is_file()}
        (args.out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')


if __name__=='__main__':
    main()
