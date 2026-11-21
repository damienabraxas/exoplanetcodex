"""Native serial MCHF valence-SD comparison, no empirical gf adjustment.

Use explicit SD replacements from both open reference shells, not GENCL's
different 'type 2' complex occupation restriction. Freeze n<3 then older layers.
LS multiplet totals must be compared only with summed relativistic components.
"""
import argparse
from datetime import datetime,timezone
import json
from pathlib import Path
import platform
import shutil
from rya1354_solver_smoke import run,sha


def configuration_input(orbital,term,n):
    virtual=['2p']+[f'{k}{l}' for k in range(3,n+1) for l in 'spdfg'[:k]]
    # An explicitly empty 3s places the upper generator's subshell order in the
    # same n/l sequence as the lower state. GENCL regenerates couplings and MCHF
    # solves afresh; no existing coefficient or coupling label is reordered.
    reference='2p(3)3s(1)' if orbital=='3s' else '2p(3)3s(0)3p(1)'
    last_shell=2 if orbital=='3s' else 3
    return f'\nO_{term}\n 1s 2s\n{reference}\n\n\nsd\n{",".join(virtual)}\n1\n{last_shell}\n{term}\n\n'


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--atsp',type=Path,required=True)
    ap.add_argument('--seed',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--n-max',type=int,choices=(3,4,5),default=4)
    args=ap.parse_args()
    if platform.node().lower().split('.')[0]!='sirius':
        raise RuntimeError('Atomic calculations belong on Sirius')
    args.atsp=args.atsp.resolve()
    args.out.mkdir(parents=True,exist_ok=False)
    shutil.copyfile(Path(__file__),args.out/'driver.py')
    receipt=dict(schema='rya1354.native-mchf-ladder.v1',status='RUNNING',scientific_admission='HOLD',
                 started_utc=datetime.now(timezone.utc).isoformat(),steps=[],
                 model='Nonrelativistic MCHF; closed 1s2 2s2 occupations; explicit SD from 2p and 3s/3p; freeze n<3 radial functions then older layers; fixed SCF controls from upstream fixture, no Solar fit')
    try:
        for n in range(3,args.n_max+1):
            for tag,orbital,term,parity in [('lower','3s','5S','o'),('upper','3p','5P','e')]:
                out=args.out/f'n{n}_{tag}'
                out.mkdir()
                receipt['steps'].append(run(args.atsp/'bin/gencl',configuration_input(orbital,term,n),out,600))
                receipt['steps'].append(run(args.atsp/'bin/nonh','',out,600))
                previous=(args.seed/f'atsp_{orbital}_{term}'/'wfn.out' if n==3 else
                          args.out/f'n{n-1}_{tag}'/'wfn.out')
                shutil.copyfile(previous,out/'wfn.inp')
                varied=','.join(f'{n}{l}' for l in 'spdfg'[:n])
                spectroscopic=f'1s 2s 2p {orbital}' if n==3 else ''
                stdin=f'{term},8.\n1,1\n{varied}\n{spectroscopic}\ny\nn\ny\nn\n.f.,1.e-9,1.e-9\nn\n500,0\ny\n'
                receipt['steps'].append(run(args.atsp/'bin/mchf',stdin,out,3600))
                for src,dst in [('cfg.inp','state.c'),('wfn.out','state.w'),(f'{term}{parity}.l','state.l')]:
                    shutil.copyfile(out/src,out/dst)
            transition=args.out/f'n{n}_E1'
            transition.mkdir()
            for tag in ('lower','upper'):
                for ext in ('c','w','l'):
                    shutil.copyfile(args.out/f'n{n}_{tag}'/f'state.{ext}',transition/f'{tag}.{ext}')
            receipt['steps'].append(run(args.atsp/'bin/biotr','lower\nupper\nn\nn\nc\nE1\n*\n',transition,600))
            if receipt['steps'][-1]['semantic_errors'] or not (transition/'lower.upper.ls').is_file():
                raise ValueError('Biorthogonal transition refused or required native output absent')
        if any(s['semantic_errors'] for s in receipt['steps']):
            raise ValueError('Native scientific refusal despite zero process exit')
        receipt['status']='NATIVE_MCHF_LADDER_COMPLETED_NUMERICAL_REVIEW_REQUIRED'
    except Exception as exc:
        receipt.update(status='FAILED_RETAINED',error=f'{type(exc).__name__}: {exc}')
        raise
    finally:
        receipt['files']={str(p.relative_to(args.out)):{'sha256':sha(p),'bytes':p.stat().st_size}
                          for p in sorted(args.out.rglob('*')) if p.is_file()}
        (args.out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')


if __name__=='__main__':
    main()
