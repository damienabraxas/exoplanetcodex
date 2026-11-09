"""Balanced fixed-orbital core/Hamiltonian controls for the triplet decay."""
import argparse
from datetime import datetime,timezone
import json
from pathlib import Path
import platform
import shutil

from rya1354_solver_smoke import run,sha
from rya1354_mr_inventory import references
from rya1354_lifetime_ladder import STATES,identity
from rya1354_triplet_correlation import transitions,require_triplet_quintet_order


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--grasp-bin',type=Path,required=True)
    ap.add_argument('--ladder',type=Path,required=True)
    ap.add_argument('--n',type=int,choices=(3,4,5),required=True)
    ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--multireference',action='store_true',help='Explicit n4 reference-family control; includes same-space SR DC and VV/CV MR DC plus CV MR Breit')
    args=ap.parse_args()
    if args.multireference and args.n!=4:
        raise ValueError('Declared multireference inventory is restricted to n4')
    if platform.node().lower().split('.')[0]!='sirius':
        raise RuntimeError('Atomic solvers belong on Sirius')
    args.grasp_bin=args.grasp_bin.resolve()
    if args.multireference:
        parent=json.loads((args.ladder/'receipt.json').read_text())
        if parent.get('status')!='BALANCED_LIFETIME_LADDER_COMPLETED_REVIEW_REQUIRED':
            raise ValueError('Balanced initializer campaign is not qualified')
        for state in STATES:
            for name in ('state.c','state.w','isodata','state.lsj.lbl'):
                relative=f'n{args.n}_{state}_verify/{name}'
                file=args.ladder/relative
                record=parent['files'][relative]
                if sha(file)!=record['sha256'] or file.stat().st_size!=record['bytes']:
                    raise ValueError('MR initializer differs from native receipt: '+relative)
    args.out.mkdir(parents=True,exist_ok=False)
    shutil.copyfile(__file__,args.out/'driver.py')
    receipt=dict(schema='rya1354.balanced-lifetime-core.v1',status='RUNNING',scientific_admission='HOLD',
                 started_utc=datetime.now(timezone.utc).isoformat(),steps=[],n=args.n,
                 multireference=args.multireference,
                 model='Fixed balanced strict orbitals; native RCI interpolation grid; identical VV DC control, full 2s-hole SD/DC, same core CSFs plus static-limit Breit. 1s frozen, point stationary Z8, recoil off. No observed-energy substitution.')
    if args.multireference:
        receipt['model']='Fixed qualified n4 state-specific VV orbitals. Same-space SR/DC control; n<=4 Li2023 Table1 MCDHF reference-family subset expanded by VV-SD/DC, full2s-SD/DC, then same core MR CSFs plus static Breit. Existing n4 radial functions fixed; no MR orbital optimization or published-protocol reproduction. 1s inactive, no observed inputs.'
        receipt['initializer_receipt_sha256']=sha(args.ladder/'receipt.json')
    try:
        dc='n\nn\nn\nn\nn\nn\n'
        breit='y\ny\n1.d-6\nn\nn\nn\nn\n'
        cases=([('DC',False,dc,False),('VV_MR_DC',False,dc,True),
                ('CV_2s_MR_DC',True,dc,True),('CV_2s_MR_Breit',True,breit,True)]
               if args.multireference else
               [('DC',False,dc,False),('CV_2s_DC',True,dc,False),('CV_2s_Breit',True,breit,False)])
        receipt['case_specs']=[dict(name=name,core=core,multireference=mr,physics=physics)
                               for name,core,physics,mr in cases]
        for case,core,physics,mr in cases:
            selected={}
            for state,spec in STATES.items():
                seed=args.ladder/f'n{args.n}_{state}_verify'
                out=args.out/(case+'_'+state)
                out.mkdir()
                identity((seed/'state.lsj.lbl').read_text(),state)
                for file in ('isodata','state.w'):
                    shutil.copyfile(seed/file,out/file)
                if core or mr:
                    active=','.join(f'{args.n}{l}' for l in 'spdfg'[:args.n])
                    config=('\n'.join(references('even' if state=='ground' else 'odd',core,True)) if mr
                            else spec['config'].format(occupation='*').replace('2s(2,i)','2s(2,*)'))
                    receipt['steps'].append(run(args.grasp_bin/'rcsfgenerate',f'*\n0\n{config}\n\n{active}\n{spec["js"]}\n2\nn\n',out,1800))
                    shutil.copyfile(out/'rcsf.out',out/'state.c')
                else:
                    shutil.copyfile(seed/'state.c',out/'state.c')
                shutil.copyfile(out/'state.c',out/'rcsf.inp')
                roots=spec['roots'].removesuffix('5\n') if state=='ground' else spec['roots']
                for exe,stdin in [('rangular','y\n'),('rci',f'y\nstate\n{physics}{roots}'),('jj2lsj','state\ny\nn\ny\n')]:
                    receipt['steps'].append(run(args.grasp_bin/exe,stdin,out,43200))
                    if receipt['steps'][-1]['semantic_errors']:
                        raise ValueError('Native scientific refusal')
                if not (out/'state.cm').is_file():
                    raise ValueError('Native CI mixing file absent')
                identity((out/'state.lsj.lbl').read_text(),state)
                selected[state]=out
            require_triplet_quintet_order((selected['triplet']/'state.lsj.lbl').read_text(),(selected['quintet']/'state.lsj.lbl').read_text())
            for name,lower,operators,same in [('to_ground',selected['ground'],'E1,M2,E3',False),
                                              ('to_5S',selected['quintet'],'M1,E2,M3',True)]:
                receipt['steps']+=transitions(args.grasp_bin,lower,selected['triplet'],args.out/(case+'_'+name),operators,same,ci=True)
            receipt['last_completed_case']=case
            (args.out/'progress.json').write_text(json.dumps(receipt,indent=2)+'\n')
        receipt['status']='BALANCED_CORE_LIFETIMES_COMPLETED_REVIEW_REQUIRED'
    except Exception as exc:
        receipt.update(status='FAILED_RETAINED',error=f'{type(exc).__name__}: {exc}')
        raise
    finally:
        receipt['files']={str(p.relative_to(args.out)):dict(sha256=sha(p),bytes=p.stat().st_size)
                          for p in sorted(args.out.rglob('*')) if p.is_file()}
        (args.out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')


if __name__=='__main__':
    main()
