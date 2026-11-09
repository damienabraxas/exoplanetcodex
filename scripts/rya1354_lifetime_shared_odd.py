"""Matched n3 core model with native statistical EOL orbitals shared by 3S/5S.

This is a separate physical orbital model, not a fallback relabeling of failed
state-specific SCF. RCI regenerates individual native J blocks, without editing
mixing coefficients, before the complete triplet decay inventory is evaluated.
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
from rya1354_lifetime_ladder import STATES,identity
from rya1354_two_method_summary import scf_diagnostics
from rya1354_triplet_correlation import transitions,require_triplet_quintet_order


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--grasp-bin',type=Path,required=True)
    ap.add_argument('--strict-bin',type=Path,required=True)
    ap.add_argument('--core',type=Path,required=True)
    ap.add_argument('--ground',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--cycles',type=int,default=500)
    ap.add_argument('--resume-odd',type=Path,help='Explicit same-CSF continuation initializer, never a qualified state by itself')
    args=ap.parse_args()
    if args.cycles<1:
        raise ValueError('Positive iteration allowance required')
    if platform.node().lower().split('.')[0]!='sirius':
        raise RuntimeError('Atomic solvers belong on Sirius')
    args.grasp_bin=args.grasp_bin.resolve();args.strict_bin=args.strict_bin.resolve()
    args.out.mkdir(parents=True,exist_ok=False)
    shutil.copyfile(__file__,args.out/'driver.py')
    receipt=dict(schema='rya1354.shared-odd-lifetime.v1',status='RUNNING',scientific_admission='HOLD',
                 started_utc=datetime.now(timezone.utc).isoformat(),steps=[],diagnostics={},
                 cycles=args.cycles,resume_odd=str(args.resume_odd) if args.resume_odd else None,
                 model='n3 full2s SD, fixed1s, ground5-state EOL independently qualified; odd3S1/5S2 share orbitals with native statistical3:5 EOL weights. Reoptimize n2/n3; fine439grid/original ACCY; DC only, no observed inputs')
    try:
        identity((args.ground/'state.lsj.lbl').read_text(),'ground')
        gd=scf_diagnostics((args.ground/'rmcdhf.stdout').read_text(),(args.ground/'rmcdhf.stderr').read_text())
        if gd['iteration_limit_hit'] or not gd['all_reported_functions_below_ACC']:
            raise ValueError('Ground initializer is not qualified')
        seed=args.core/'CV_2s_DC_quintet'
        for verify in (False,True):
            out=args.out/('odd_verify' if verify else 'odd')
            out.mkdir()
            if not verify:
                for exe,stdin in [('rnucleus','8\n0\n0\n0\n0\n0\n'),
                                  ('rcsfgenerate','*\n0\n1s(2,i)2s(2,*)2p(3,*)3s(1,*)\n\n3s,3p,3d\n2,4\n2\nn\n')]:
                    receipt['steps'].append(run(args.grasp_bin/exe,stdin,out,600))
                shutil.copyfile(out/'rcsf.out',out/'rcsf.inp');shutil.copyfile(seed/'state.w',out/'rwfn.inp')
                if args.resume_odd:
                    if ((args.resume_odd/'rcsf.inp').read_bytes()!=(out/'rcsf.inp').read_bytes() or
                            (args.resume_odd/'isodata').read_bytes()!=(out/'isodata').read_bytes()):
                        raise ValueError('Shared-EOL resume model identity differs')
                    shutil.copyfile(args.resume_odd/'rwfn.out',out/'rwfn.inp')
            else:
                for src,dst in [('state.c','rcsf.inp'),('state.w','rwfn.inp'),('isodata','isodata')]:
                    shutil.copyfile(args.out/'odd'/src,out/dst)
            receipt['steps'].append(run(args.grasp_bin/'rangular','y\n',out,1800))
            grid=f'n\nn\ny\nn\ny\n{math.exp(-65/16)/8:.17g}\n0.03125\n0\n439\ny\n{.0625**6:.17g}\n'
            nsic='n\n' if verify else 'y\n20\n'
            controls=f'{grid}1\n1\n5\n2*,3*\n2s,2p,3s\n{args.cycles}\ny\nn\nn\nn\ny\nn\n-0.5\nn\n{nsic}n\nn\n1\n'
            receipt['steps'].append(run(args.strict_bin/'rmcdhf',controls,out,43200))
            if any(step['semantic_errors'] for step in receipt['steps']):
                raise ValueError('Native shared-EOL refusal; no lifetime')
            d=scf_diagnostics((out/'rmcdhf.stdout').read_text(),(out/'rmcdhf.stderr').read_text())
            if d['iteration_limit_hit'] or not d['all_reported_functions_below_ACC']:
                raise ValueError('Shared odd EOL orbital criterion unmet')
            receipt['diagnostics'][out.name]=d
            for src,ext in [('rcsf.inp','c'),('rwfn.out','w'),('rmix.out','m')]:
                shutil.copyfile(out/src,out/f'state.{ext}')
            receipt['steps'].append(run(args.grasp_bin/'jj2lsj','state\nn\nn\ny\n',out,600))
            blocks=re.findall(r'^\s*\d+\s+(\d+)\s+([+-])\s+[^\n]+\n([^\n]+)',(out/'state.lsj.lbl').read_text(),re.M)
            if len(blocks)!=2 or any(j!=ej or p!='-' or '2p(3)4S' not in line or not line.rstrip().endswith(term)
                                      for (j,p,line),(ej,term) in zip(blocks,[('1','3s_3S'),('2','3s_5S')])):
                raise ValueError('Shared odd physical identities changed')
        selected={}
        for state,spec in STATES.items():
            out=args.out/('CI_'+state);out.mkdir()
            source=args.ground if state=='ground' else args.out/'odd_verify'
            for filename in ('isodata','state.w'):
                shutil.copyfile(source/filename,out/filename)
            if state=='ground':
                shutil.copyfile(source/'state.c',out/'state.c')
            else:
                config=spec['config'].format(occupation='*').replace('2s(2,i)','2s(2,*)')
                receipt['steps'].append(run(args.grasp_bin/'rcsfgenerate',f'*\n0\n{config}\n\n3s,3p,3d\n{spec["js"]}\n2\nn\n',out,600))
                shutil.copyfile(out/'rcsf.out',out/'state.c')
            shutil.copyfile(out/'state.c',out/'rcsf.inp')
            roots=spec['roots'].removesuffix('5\n') if state=='ground' else spec['roots']
            for exe,stdin in [('rangular','y\n'),('rci',f'y\nstate\nn\nn\nn\nn\nn\nn\n{roots}'),('jj2lsj','state\ny\nn\ny\n')]:
                receipt['steps'].append(run(args.grasp_bin/exe,stdin,out,1800))
            identity((out/'state.lsj.lbl').read_text(),state);selected[state]=out
        require_triplet_quintet_order((selected['triplet']/'state.lsj.lbl').read_text(),(selected['quintet']/'state.lsj.lbl').read_text())
        for name,lower,operators,same in [('to_ground',selected['ground'],'E1,M2,E3',False),
                                          ('to_5S',selected['quintet'],'M1,E2,M3',True)]:
            receipt['steps']+=transitions(args.grasp_bin,lower,selected['triplet'],args.out/name,operators,same,ci=True)
        if any(step['semantic_errors'] for step in receipt['steps']):
            raise ValueError('Native scientific refusal')
        receipt['status']='SHARED_ODD_LIFETIME_COMPLETED_REVIEW_REQUIRED'
    except Exception as exc:
        receipt.update(status='FAILED_RETAINED',error=f'{type(exc).__name__}: {exc}')
        raise
    finally:
        receipt['files']={str(p.relative_to(args.out)):dict(sha256=sha(p),bytes=p.stat().st_size)
                          for p in sorted(args.out.rglob('*')) if p.is_file()}
        (args.out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')


if __name__=='__main__':
    main()
