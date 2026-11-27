"""RYA-1395 v71: finish the O I spectroscopic (s0) orbital stage, one parity per process (Sirius only).

Declared on RYA-1395 (comment b077e4a2) before any build or run. Reuses the RYA-1354 v60 sequential driver
(scripts/rya1354_mcdhf_s0_sequential.py) UNCHANGED except:
  1. nf steps (even 4f, 5f) vary nothing: the new orbital is screened hydrogenic (rwfnestimate option 3), not
     Thomas-Fermi, and the step is a fixed-orbital CI with the usual identity check. Evidence: v60's converged
     4f j=5/2 has eps 0.0312520 Ha vs hydrogenic 1/32 (rel 6e-5); no graded O I line involves 2p3 nf.
  2. one parity per process (--parity), so an even-parity stop cannot block odd parity.
  3. E1 between the final even and odd sets is a separate invocation (--transitions EVEN_DIR ODD_DIR).
Strict criterion, half grid, method 3, 500 cycles, memory cap, identity and final-list guards: as v60.
Report only; nothing adopted.
"""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import platform
import shutil
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
import rya1354_mcdhf_s0_sequential as S  # noqa: E402
from rya1354_rci_feasibility import BudgetExceeded  # noqa: E402
from rya1354_solver_smoke import sha  # noqa: E402
from rya1354_triplet_correlation import transitions  # noqa: E402

SCHEMA = 'rya1395.mcdhf-s0-v71'
NF = {'4f', '5f'}
S.STEPS['even'] = [(lab, fams, [] if lab in NF else varied) for lab, fams, varied in S.STEPS['even']]
_orig_estimate, _orig_step = S.step_estimate, S.step
_current = {'label': None}


def nf_estimate(n_previous):
    """Frozen orbitals from every earlier step (newest first), then screened hydrogenic for the new nf."""
    return f'n\nn\n\ny\nn\ny\n{S.GRID}' + ''.join(f'1\nprev{i}.w\n*\n' for i in range(n_previous)) + '3\n*\nn\n'


def _estimate(n_previous):
    return nf_estimate(n_previous) if _current['label'] in NF else _orig_estimate(n_previous)


def step(a, receipt, parity, k, previous):
    _current['label'] = S.STEPS[parity][k][0]
    try:
        return _orig_step(a, receipt, parity, k, previous)
    finally:
        _current['label'] = None


S.step_estimate = _estimate


def write_files(out, receipt):
    seen, files = {}, {}
    for p in sorted(out.rglob('*')):
        if p.is_file() and p.name != 'receipt.json':
            key = (p.stat().st_dev, p.stat().st_ino)
            if key not in seen:
                seen[key] = dict(sha256=sha(p), bytes=p.stat().st_size)
            files[str(p.relative_to(out))] = seen[key]
    receipt['files'] = files
    (out / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')


def run_parity(a):
    v56 = json.loads(a.v56_receipt.read_text())['inputs']
    receipt = dict(schema=SCHEMA, parity=a.parity, grid=dict(RNT=S.RNT, H=S.H, N=S.N, ACCY=S.ACCY), status='RUNNING',
                   scientific_admission='HOLD', exposure='EXPOSED_NONBLIND', adoption_allowed=False,
                   started_utc=datetime.now(timezone.utc).isoformat(), steps={}, nscf=S.NSCF, mem_cap_bytes=S.MEM_CAP_BYTES,
                   v56_receipt_sha256=sha(a.v56_receipt), nf_policy='fixed screened hydrogenic (option b)',
                   plan=[dict(label=l, families=[f'{f}{n}' for f, n in fams], varied=v) for l, fams, v in S.STEPS[a.parity]],
                   stopped=None)
    try:
        done = []
        for k in range(len(S.STEPS[a.parity])):
            done.insert(0, step(a, receipt, a.parity, k, list(done)))
        final = receipt['steps'][done[0].name]
        if final['blocks'] != v56[a.parity]['blocks']:
            raise ValueError(f'final {a.parity} list differs from the v56 s0 list')
        final['equals_full_s0_list'] = True
        receipt['final_dir'] = done[0].name
        receipt['status'] = 'S0_PARITY_COMPLETED_REVIEW_REQUIRED'
    except BudgetExceeded as exc:
        receipt.update(status='MEASURED_STOP_REVIEW_REQUIRED', stopped=str(exc))
    except Exception as exc:
        receipt.update(status='FAILED_RETAINED', error=f'{type(exc).__name__}: {exc}')
        raise
    finally:
        receipt['finished_utc'] = datetime.now(timezone.utc).isoformat()
        write_files(a.out, receipt)


def run_transitions(a):
    even, odd = (Path(p).resolve() for p in a.transitions)
    receipt = dict(schema=SCHEMA + '.e1', even=str(even), odd=str(odd), status='RUNNING',
                   started_utc=datetime.now(timezone.utc).isoformat())
    try:
        receipt['transitions'] = transitions(a.grasp_bin, even, odd, a.out / 's0_E1', 'E1', False)
        receipt['status'] = 'S0_E1_COMPLETED_REVIEW_REQUIRED'
    except Exception as exc:
        receipt.update(status='FAILED_RETAINED', error=f'{type(exc).__name__}: {exc}')
        raise
    finally:
        receipt['finished_utc'] = datetime.now(timezone.utc).isoformat()
        write_files(a.out, receipt)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--grasp-bin', type=Path, required=True)
    ap.add_argument('--strict-bin', type=Path)
    ap.add_argument('--v56-receipt', type=Path)
    ap.add_argument('--parity', choices=['even', 'odd'])
    ap.add_argument('--transitions', nargs=2, metavar=('EVEN_FINAL_DIR', 'ODD_FINAL_DIR'))
    ap.add_argument('--out', type=Path, required=True)
    a = ap.parse_args()
    if platform.node().lower().split('.')[0] != 'sirius':
        raise RuntimeError('Atomic solvers belong on Sirius')
    if bool(a.parity) == bool(a.transitions):
        raise SystemExit('give exactly one of --parity or --transitions')
    for k, v in vars(a).items():
        if isinstance(v, Path):
            setattr(a, k, v.resolve())
    a.out.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, a.out / 'driver.py')
    shutil.copyfile(S.__file__, a.out / 'rya1354_mcdhf_s0_sequential.py')
    if a.parity:
        if not (a.strict_bin and a.v56_receipt):
            raise SystemExit('--parity needs --strict-bin and --v56-receipt')
        run_parity(a)
    else:
        run_transitions(a)


if __name__ == '__main__':
    main()
