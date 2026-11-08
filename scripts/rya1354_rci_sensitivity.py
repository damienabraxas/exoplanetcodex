"""Fixed-orbital n3 RCI Hamiltonian/core sensitivity, never gf fitting.

RCI's pinned getcid.f90 uses its own interpolating grid (RNT=2e-6/Z,
H=.05, N=590). First run Dirac-Coulomb on the identical CSFs as a numerical
control. Isolate Breit, VP/self-energy and single-2s-hole CV within that grid.
"""
import argparse
from datetime import datetime,timezone
import json
from pathlib import Path
import platform
import shutil
from rya1354_solver_smoke import run,sha


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--grasp-bin',type=Path,required=True)
    ap.add_argument('--ladder',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--dummy-positive-mass',action='store_true',
                    help='Numerical DC warning control only: 1 amu synthetic nucleus, both recoil terms off; not an oxygen isotope candidate')
    args=ap.parse_args()
    if platform.node().lower().split('.')[0]!='sirius':
        raise RuntimeError('Atomic calculations belong on Sirius')
    args.grasp_bin=args.grasp_bin.resolve()
    args.out.mkdir(parents=True,exist_ok=False)
    shutil.copyfile(Path(__file__),args.out/'driver.py')
    receipt=dict(schema='rya1354.rci-sensitivity.v1',status='RUNNING',scientific_admission='HOLD',
                 started_utc=datetime.now(timezone.utc).isoformat(),steps=[],
                 scope='n3 only; fixed retained orbitals; point static Z8; no observed energy or gf fitting; RCI grid differs from MCDHF')
    cases=[('DC','n\nn\nn\nn\nn\nn\n',False),
           ('Breit_static_limit','y\ny\n1.d-6\nn\nn\nn\nn\n',False),
           ('VP_self_energy','n\nn\ny\nn\nn\ny\n3\n',False),
           ('CV_2s_one_hole_DC','n\nn\nn\nn\nn\nn\n',True)]
    if args.dummy_positive_mass:
        cases=[('DC_dummy_mass_recoil_off','n\nn\nn\nn\nn\nn\n',False)]
        receipt['scope']='SYNTHETIC_NUMERICAL_WARNING_CONTROL_ONLY; 1 amu input, recoil disabled; not a physical isotope or adopted candidate'
    try:
        for case,physics,core in cases:
            for tag,orbital,js,states in [('lower','3s','4,4','1\n'),('upper','3p','2,6','2\n3\n1\n')]:
                out=args.out/f'{case}_{tag}'
                out.mkdir()
                source=args.ladder/f'n3_{tag}'
                shutil.copyfile(source/'isodata',out/'isodata')
                if args.dummy_positive_mass:
                    lines=(out/'isodata').read_text().splitlines()
                    if lines.count('Mass of nucleus (in amu):')!=1:
                        raise ValueError('Nuclear mass input identity ambiguous')
                    lines[lines.index('Mass of nucleus (in amu):')+1]='1.0000000000000000'
                    (out/'isodata').write_text('\n'.join(lines)+'\n')
                shutil.copyfile(source/'state.w',out/'state.w')
                if core:
                    generator=f'*\n0\n1s(2,i)2s(2,1)2p(3,*){orbital}(1,*)\n\n3s,3p,3d\n{js}\n2\nn\n'
                    receipt['steps'].append(run(args.grasp_bin/'rcsfgenerate',generator,out,600))
                    shutil.copyfile(out/'rcsf.out',out/'state.c')
                else:
                    shutil.copyfile(source/'state.c',out/'state.c')
                shutil.copyfile(out/'state.c',out/'rcsf.inp')
                receipt['steps'].append(run(args.grasp_bin/'rangular','y\n',out,600))
                receipt['steps'].append(run(args.grasp_bin/'rci',f'y\nstate\n{physics}{states}',out,1800))
                if not (out/'state.cm').is_file():
                    raise ValueError('Native CI mixing file absent')
                receipt['steps'].append(run(args.grasp_bin/'jj2lsj','state\ny\nn\ny\n',out,600))
                text=(out/'state.lsj.lbl').read_text()
                if text.count('3s_5S' if tag=='lower' else '3p_5P')!=(1 if tag=='lower' else 3):
                    raise ValueError('Native target LS identities changed')
            out=args.out/f'{case}_E1'
            out.mkdir()
            for tag in ('lower','upper'):
                for ext in ('c','cm','w'):
                    shutil.copyfile(args.out/f'{case}_{tag}'/f'state.{ext}',out/f'{tag}.{ext}')
            shutil.copyfile(args.out/f'{case}_lower'/'isodata',out/'isodata')
            for exe,stdin in [('rbiotransform','y\ny\nlower\nupper\ny\n'),
                              ('rtransition','y\ny\nlower\nupper\nE1\n')]:
                receipt['steps'].append(run(args.grasp_bin/exe,stdin,out,1800))
            if not (out/'lower.upper.ct').is_file() and not (out/'lower.upper.t').is_file():
                raise ValueError('Native CI transition table absent')
        if any(s['semantic_errors'] for s in receipt['steps']):
            raise ValueError('Native scientific refusal despite zero process status')
        receipt['status']='FIXED_ORBITAL_SENSITIVITIES_COMPLETED_REVIEW_REQUIRED'
    except Exception as exc:
        receipt.update(status='FAILED_RETAINED',error=f'{type(exc).__name__}: {exc}')
        raise
    finally:
        receipt['files']={str(p.relative_to(args.out)):{'sha256':sha(p),'bytes':p.stat().st_size}
                          for p in sorted(args.out.rglob('*')) if p.is_file()}
        (args.out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')


if __name__=='__main__':
    main()
