"""Exploratory valence SD ladder on Sirius; native energies, no fitted gf.

The lowest even ASFs are ground-configuration states. Explicit trial ordinals
are therefore 2,3,1 for J=1,2,3; LS output MUST establish state identity before
any transition is compared. This driver never admits a result automatically.
"""
from pathlib import Path
import argparse
from datetime import datetime, timezone
import json
import platform
import shutil
import re
from rya1354_solver_smoke import run, sha


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--grasp-bin', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--n-max', type=int, choices=(3, 4, 5), default=3)
    ap.add_argument('--n-min', type=int, choices=(3, 4, 5), default=3)
    ap.add_argument('--previous-ladder', type=Path,
                    help='Explicit retained previous layer when starting above n=3')
    ap.add_argument('--new-orbitals', choices=('thomas-fermi','hydrogenic'), default='thomas-fermi')
    ap.add_argument('--seed', type=Path, required=True,
                    help='Retained single-configuration run directory; frozen starting spectroscopic orbitals')
    args = ap.parse_args()
    args.grasp_bin = args.grasp_bin.resolve()
    args.seed = args.seed.resolve()
    if args.n_min > args.n_max or (args.n_min > 3 and args.previous_ladder is None):
        raise ValueError('A resumed ladder requires an explicit previous layer and valid n range')
    if platform.node().lower().split('.')[0] != 'sirius':
        raise RuntimeError('Solver calculations belong on Sirius')
    args.out.mkdir(parents=True, exist_ok=False)
    receipt = dict(schema='rya1354.valence-sd.v1', status='RUNNING',
                   scientific_admission='HOLD', exposure='NONBLIND_EXPLORATORY',
                   started_utc=datetime.now(timezone.utc).isoformat(),
                   script_sha256=sha(Path(__file__)), steps=[],
                   new_orbital_initialization=args.new_orbitals,
                   model='Dirac-Coulomb; point static Z=8; fixed 1s2 2s2 occupations; valence SD; freeze n<3 radial orbitals, then freeze older layers; no Breit/QED',
                   target_ordinals={'odd_J2': [1], 'even_J123': [2, 3, 1]},
                   warning='Ordinal identity provisional until native LS composition review; gauge spread is not uncertainty')
    try:
        shutil.copyfile(Path(__file__),args.out/'driver.py')
        for n in range(args.n_min, args.n_max + 1):
            active = ','.join(f'{n}{l}' for l in 'spdfg'[:n])
            for tag, orbital, js, states in [('lower','3s','4,4','1\n'),
                                             ('upper','3p','2,6','2\n3\n1\n5\n')]:
                out = args.out / f'n{n}_{tag}'
                out.mkdir()
                previous_root=(args.previous_ladder if n==args.n_min and n>3 else args.out)
                previous = (previous_root / f'n{n-1}_{tag}' / 'state.w' if n > 3 else
                            args.seed / ('grasp_3s_J2' if tag == 'lower' else 'grasp_3p_J123') / 'state.w')
                shutil.copyfile(previous, out/'previous.w')
                method=2 if args.new_orbitals=='thomas-fermi' else 3
                init = f'y\n1\nprevious.w\n*\n{method}\n*\n'
                steps = [('rnucleus','8\n0\n0\n0\n0\n0\n'),
                         ('rcsfgenerate',f'*\n0\n1s(2,i)2s(2,i)2p(3,*){orbital}(1,*)\n\n{active}\n{js}\n2\nn\n')]
                for name, stdin in steps:
                    receipt['steps'].append(run(args.grasp_bin/name, stdin, out, 3600))
                shutil.copyfile(out/'rcsf.out',out/'rcsf.inp')
                for name, stdin in [('rangular','y\n'),('rwfnestimate',init),
                                    ('rmcdhf',f'y\n{states}{n}*\n{orbital if n == 3 else ""}\n100\n')]:
                    receipt['steps'].append(run(args.grasp_bin/name,stdin,out,3600))
                    if receipt['steps'][-1]['semantic_errors']:
                        raise ValueError(f'{name}: native scientific refusal or SCF limit')
                for src, ext in [('rcsf.inp','c'),('rmix.out','m'),('rwfn.out','w')]:
                    shutil.copyfile(out/src,out/f'state.{ext}')
                receipt['steps'].append(run(args.grasp_bin/'jj2lsj','state\nn\nn\ny\n',out,3600))
            # Require the leading composition of every selected ASF to identify
            # the configured target. Record compositions; this is not an accuracy gate.
            for tag, expected_j, term in [('lower',[2],'3s_5S'), ('upper',[1,2,3],'3p_5P')]:
                text = (args.out/f'n{n}_{tag}'/'state.lsj.lbl').read_text()
                blocks = re.findall(r'^\s*\d+\s+(\d+)\s+([+-])\s+[^\n]+\n([^\n]+)',text,re.M)
                parity='-' if tag=='lower' else '+'
                if ([int(b[0]) for b in blocks] != expected_j or any(term not in b[2] for b in blocks)
                        or any(b[1]!=parity for b in blocks)):
                    raise ValueError(f'n{n} {tag}: dominant LS state identity refused')
            transition = args.out/f'n{n}_E1'
            transition.mkdir()
            for tag in ('lower','upper'):
                for ext in ('c','m','w'):
                    shutil.copyfile(args.out/f'n{n}_{tag}'/f'state.{ext}',transition/f'{tag}.{ext}')
            shutil.copyfile(args.out/f'n{n}_lower'/'isodata',transition/'isodata')
            for name, stdin in [('rbiotransform','y\nn\nlower\nupper\ny\n'),
                                ('rtransition','y\nn\nlower\nupper\nE1\n')]:
                receipt['steps'].append(run(args.grasp_bin/name,stdin,transition,3600))
        if any(s['semantic_errors'] for s in receipt['steps']):
            raise ValueError('Native scientific input refusal despite zero process exit')
        receipt['status']='ORBITAL_RUNS_COMPLETED_IDENTITY_AND_NUMERICAL_REVIEW_REQUIRED'
    except Exception as exc:
        receipt.update(status='FAILED_RETAINED', error=f'{type(exc).__name__}: {exc}')
        raise
    finally:
        receipt['finished_utc']=datetime.now(timezone.utc).isoformat()
        receipt['files']={str(p.relative_to(args.out)): {'sha256':sha(p),'bytes':p.stat().st_size}
                          for p in sorted(args.out.rglob('*')) if p.is_file()}
        (args.out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')


if __name__ == '__main__':
    main()
