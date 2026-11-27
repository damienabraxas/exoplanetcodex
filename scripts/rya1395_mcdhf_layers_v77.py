"""RYA-1395 v77: first correlated O I MCDHF-EOL layers (n = 3..6) on top of OUR s0 finals (Sirius).

Declared on RYA-1395 before build. Reuses the v53 layer driver (rya1354_mcdhf_eol.stage: Li 2023 MCDHF MR, SD from
2p and the outer electron, 1s/2s closed; layer n adds every l<n orbital of principal number n not in the MR and
varies only those; identity gate on the 45 even / 36 odd targets) with these changes:
  * previous stage = our s0 finals (v72 even_7_ci, v75 odd_7_ci) instead of a v53 s0 stage;
  * half grid + strict ACCY (6e-8), method 3 for the varied (non-spectroscopic) correlation orbitals;
  * new orbitals: frozen previous orbitals, then Thomas-Fermi, then screened hydrogenic;
  * v75 measured-floor admission for correlation orbitals (final self-consistency <= 1e-6 AND last 50 rows
    identical); any other non-convergence or refusal is FAILED_RETAINED;
  * after each layer: E1 (B and C) + the Li 2023 Table A2 report (rya1395_s0_e1_report).
One niced process; parities in sequence.
"""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import platform
import shutil
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
import rya1354_mcdhf_eol as E  # noqa: E402
import rya1354_mcdhf_s0_sequential as S  # noqa: E402
import rya1395_s0_e1_report as REP  # noqa: E402

SCHEMA = 'rya1395.mcdhf-layers-v77'
LAYERS = (3, 4, 5, 6)
FLOOR_CEIL, STAGNANT = 1e-6, 50
LIMIT = "rc=0 ['MAXIMUM ITERATIONS IN SCF EXCEEDED']"
_state = dict(varied=set(), admitted=None)


def rmcdhf_stdin(levels, varied, spectroscopic):
    if spectroscopic:
        raise ValueError('v77 layers vary correlation orbitals only')
    _state['varied'] = {v.rstrip('*') for v in varied.split(',')}
    return (f'n\nn\ny\nn\ny\n{S.GRID}y\n{S.ACCY:.17g}\n{levels}5\n{varied}\n\n{S.NSCF}\n'
            f'y\nn\ny\n\n\n{varied}\n\n'
            'n\nn\nn\nn\nn\nn\n1\n')


def floor_check(stdout_text, varied):
    rows = {}
    for line in stdout_text.splitlines():
        f = line.split()
        if len(f) >= 11 and f[0].rstrip('-') in varied and f[2] in ('1', '2', '3', '4'):
            rows.setdefault(f[0], []).append(line.strip())
    out = {}
    for lab, r in rows.items():
        c = float(r[-1].split()[4].replace('D', 'E'))
        stag = len(r) >= STAGNANT and len(set(r[-STAGNANT:])) == 1
        out[lab] = dict(consistency=c, rows=len(r), stagnant=stag, pass_=bool(c <= FLOOR_CEIL and stag))
    return out


_timed = E.timed


def timed(exe, stdin, cwd, *args, **kw):
    try:
        return _timed(exe, stdin, cwd, *args, **kw)
    except ValueError as exc:
        if not (Path(str(exe)).name == 'rmcdhf' and str(exc).endswith(LIMIT)):
            raise
        chk = floor_check((Path(cwd) / 'rmcdhf.stdout').read_text(), _state['varied'])
        if not chk or not all(c['pass_'] for c in chk.values()):
            raise ValueError(f'{exc}; measured-floor admission refused: {chk}') from exc
        _state['admitted'] = chk
        return dict(executable=str(exe), returncode=0, semantic_errors=['MAXIMUM ITERATIONS IN SCF EXCEEDED'],
                    floor_admitted=chk, floor_ceiling=FLOOR_CEIL, stagnant_rows=STAGNANT)


_scf = E.scf_diagnostics


def scf_diagnostics(stdout, stderr):
    d = _scf(stdout, stderr)
    if _state['admitted'] is not None:
        d['strict_original'] = dict(iteration_limit_hit=d['iteration_limit_hit'],
                                    all_reported_functions_below_ACC=d['all_reported_functions_below_ACC'])
        d['floor_admitted'] = _state['admitted']
        d['iteration_limit_hit'] = False
        d['all_reported_functions_below_ACC'] = True
        _state['admitted'] = None
    return d


E.rmcdhf_stdin = rmcdhf_stdin
E.LAYER_ESTIMATE = S.STEP_ESTIMATE
E.timed = timed
E.scf_diagnostics = scf_diagnostics


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--grasp-bin', type=Path, required=True)
    ap.add_argument('--strict-bin', type=Path, required=True)
    ap.add_argument('--even-s0', type=Path, required=True)
    ap.add_argument('--odd-s0', type=Path, required=True)
    ap.add_argument('--nist-table', type=Path, required=True)
    ap.add_argument('--li-table', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--layers', type=int, nargs='+', default=list(LAYERS))
    a = ap.parse_args()
    if platform.node().lower().split('.')[0] != 'sirius':
        raise RuntimeError('Atomic solvers belong on Sirius')
    for k, v in vars(a).items():
        if isinstance(v, Path):
            setattr(a, k, v.resolve())
    a.v47 = None
    a.out.mkdir(parents=True, exist_ok=False)
    for f in (__file__, E.__file__, S.__file__, REP.__file__):
        shutil.copyfile(f, a.out / Path(f).name)
    receipt = dict(schema=SCHEMA, status='RUNNING', scientific_admission='HOLD', adoption_allowed=False,
                   started_utc=datetime.now(timezone.utc).isoformat(), stages={}, transitions={}, reports={},
                   s0=dict(even=str(a.even_s0), odd=str(a.odd_s0)), grid=dict(RNT=S.RNT, H=S.H, N=S.N, ACCY=S.ACCY),
                   targets={p: {str(j): k for j, k in E.roots(p).items()} for p in ('even', 'odd')}, stopped=None)
    prev = dict(even=a.even_s0, odd=a.odd_s0)
    try:
        for n in a.layers:
            for parity in ('even', 'odd'):
                prev[parity] = E.stage(a, receipt, parity, n, prev[parity])
                (a.out / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
            key = f'n{n}_E1'
            receipt['transitions'][key] = E.transitions(a.grasp_bin, prev['even'], prev['odd'], a.out / key, 'E1', False)
            REP.main(a.out / key, a.nist_table, a.li_table, a.out / key / 'loggf_vs_li2023.csv')
            receipt['last_completed_stage'] = key
            (a.out / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
        receipt['status'] = 'LAYERS_COMPLETED_REVIEW_REQUIRED'
    except E.BudgetExceeded as exc:
        receipt.update(status='MEASURED_STOP_REVIEW_REQUIRED', stopped=str(exc))
    except Exception as exc:
        receipt.update(status='FAILED_RETAINED', error=f'{type(exc).__name__}: {exc}')
        raise
    finally:
        receipt['finished_utc'] = datetime.now(timezone.utc).isoformat()
        (a.out / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')


if __name__ == '__main__':
    main()
