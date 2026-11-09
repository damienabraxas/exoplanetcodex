"""Independent LS MCHF check of the dominant 3s3S -> 2p4 3P E1 decay.

This is a partial E1 rate, not a certified complete relativistic lifetime.
Use native HF seeds, modern GENCL/nonh/MCHF/biotr, no empirical adjustments.
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


def validate_hf_seed(text,term):
    if 'UNABLE TO SOLVE EQUATION' in text or re.search(r'\bNaN\b',text,re.I):
        raise ValueError('Native HF radial solve failed despite possible zero exit')
    orbitals=re.search(r'There are\s+\d+ orbitals as follows:\s*\n([^\n]+)',text)
    expected=['1s','2s','2p']+(['3s'] if term=='3S' else [])
    if orbitals is None or orbitals[1].split()!=expected:
        raise ValueError('Native HF orbital identity differs from physical reference')


def configuration(term,n):
    reference='2p(4)' if term=='3P' else '2p(3)3s(1)'
    last=1 if term=='3P' else 2
    virtual=['2p']+[f'{k}{l}' for k in range(3,n+1) for l in 'spdfg'[:k]]
    return f'\nO_{term}\n 1s 2s\n{reference}\n\n\nsd\n{",".join(virtual)}\n1\n{last}\n{term}\n\n'


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--atsp',type=Path,required=True)
    ap.add_argument('--book-bin',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--n-max',type=int,choices=(3,4,5),default=5)
    ap.add_argument('--strict-bin',type=Path)
    ap.add_argument('--function-tolerance',type=float,default=1e-9)
    ap.add_argument('--qualify-reference',action='store_true',help='Native zero-replacement MCHF optimization of all reference orbitals before freezing')
    args=ap.parse_args()
    if not math.isfinite(args.function_tolerance) or args.function_tolerance<=0:
        raise ValueError('Positive finite native function tolerance required')
    if platform.node().lower().split('.')[0]!='sirius':
        raise RuntimeError('Atomic solvers belong on Sirius')
    args.atsp=args.atsp.resolve()
    args.book_bin=args.book_bin.resolve()
    if args.strict_bin:
        args.strict_bin=args.strict_bin.resolve()
    if args.qualify_reference and not args.strict_bin:
        raise ValueError('Reference qualification requires explicit strict function-stop variant')
    args.out.mkdir(parents=True,exist_ok=False)
    shutil.copyfile(__file__,args.out/'driver.py')
    receipt=dict(schema='rya1354.lifetime-mchf.v1',status='RUNNING',scientific_admission='HOLD',
                 started_utc=datetime.now(timezone.utc).isoformat(),steps=[],
                 function_tolerance=args.function_tolerance,
                 qualify_reference=args.qualify_reference,
                 scope='Nonrelativistic LS summed E1 to ground3P only, not full13-channel relativistic lifetime; valence SD/frozen1s2s and older layers; no observed energies/gf fitted')
    try:
        previous={}
        for term,reference in [('3P','2p(4)'),('3S','2p(3)3s(1)')]:
            out=args.out/('reference_'+term)
            out.mkdir()
            receipt['steps'].append(run(args.book_bin/'HF',f'O,{term},8.\n  1s  2s\n{reference}\nall\ny\ny\nn\nn\n',out,600))
            hf=(out/'HF.stdout').read_text()+(out/'HF.stderr').read_text()
            if not (out/'wfn.out').is_file():
                raise ValueError('Native HF seed failed despite possible zero exit')
            validate_hf_seed(hf,term)
            previous[term]=out/'wfn.out'
            if args.qualify_reference:
                qualified=args.out/('qualified_reference_'+term)
                qualified.mkdir()
                # Native GENCL: blank active set and replacements generate only
                # the specified reference CSFs; no hand-edited coupling records.
                generator=f'\nO_{term}\n 1s 2s\n{reference}\n\n\n\n{term}\n\n'
                receipt['steps'].append(run(args.atsp/'bin/gencl',generator,qualified,600))
                receipt['steps'].append(run(args.atsp/'bin/nonh','',qualified,600))
                shutil.copyfile(previous[term],qualified/'wfn.inp')
                orbitals='1s 2s 2p'+(' 3s' if term=='3S' else '')
                controls=f'{term},8.\n1,1\nall\n{orbitals}\ny\nn\ny\nn\n.f.,1.e-9,{args.function_tolerance:.17g}\nn\n500,0\ny\n'
                receipt['steps'].append(run(args.strict_bin/'mchf',controls,qualified,7200))
                if 'RYA1354_ORBITAL_STOP' not in (qualified/'mchf.stdout').read_text():
                    raise ValueError('Reference MCHF function criterion not established')
                previous[term]=qualified/'wfn.out'
        for n in range(3,args.n_max+1):
            states={}
            for term,parity in [('3P','e'),('3S','o')]:
                out=args.out/f'n{n}_{term}'
                out.mkdir()
                receipt['steps'].append(run(args.atsp/'bin/gencl',configuration(term,n),out,600))
                receipt['steps'].append(run(args.atsp/'bin/nonh','',out,600))
                shutil.copyfile(previous[term],out/'wfn.inp')
                varied=','.join(f'{n}{l}' for l in 'spdfg'[:n])
                spec=('1s 2s 2p'+(' 3s' if term=='3S' else '')) if n==3 else ''
                controls=f'{term},8.\n1,1\n{varied}\n{spec}\ny\nn\ny\nn\n.f.,1.e-9,{args.function_tolerance:.17g}\nn\n500,0\ny\n'
                receipt['steps'].append(run((args.strict_bin or args.atsp/'bin')/'mchf',controls,out,7200))
                if args.strict_bin and 'RYA1354_ORBITAL_STOP' not in (out/'mchf.stdout').read_text():
                    raise ValueError('MCHF orbital criterion not established')
                for src,dst in [('cfg.inp','state.c'),('wfn.out','state.w'),(term+parity+'.l','state.l')]:
                    shutil.copyfile(out/src,out/dst)
                previous[term]=out/'wfn.out'
                states[term]=out
            transition=args.out/f'n{n}_E1'
            transition.mkdir()
            for tag,term in [('lower','3P'),('upper','3S')]:
                for ext in ('c','w','l'):
                    shutil.copyfile(states[term]/('state.'+ext),transition/(tag+'.'+ext))
            receipt['steps'].append(run(args.atsp/'bin/biotr','lower\nupper\nn\nn\nc\nE1\n*\n',transition,1800))
            if not (transition/'lower.upper.ls').is_file() or any(s['semantic_errors'] for s in receipt['steps']):
                raise ValueError('Native LS transition refused or absent')
            receipt['last_completed_layer']=n
        receipt['status']='LS_DOMINANT_E1_SEQUENCE_COMPLETED_REVIEW_REQUIRED'
    except Exception as exc:
        receipt.update(status='FAILED_RETAINED',error=f'{type(exc).__name__}: {exc}')
        raise
    finally:
        receipt['files']={str(p.relative_to(args.out)):dict(sha256=sha(p),bytes=p.stat().st_size)
                          for p in sorted(args.out.rglob('*')) if p.is_file()}
        (args.out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')


if __name__=='__main__':
    main()
