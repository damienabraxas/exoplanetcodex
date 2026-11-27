"""RYA-1395 v73: odd-parity O I s0 stage with method 1 (node counting) for 5s/6s and a node-count gate (Sirius).

Declared on RYA-1395 before any build or run. Extends v72 (own hydrogenic starts, orbital-energy gate):
  1. 5s and 6s steps: the varied orbital goes in rmcdhf's method-1 slot (outward/inward integration with node
     counting) instead of method 3. v72 FAILED_RETAINED: method 3 converged 5s strictly to a 6-node root.
  2. Node gate on every step: each varied orbital's NNP in the final rmcdhf iteration must equal n - l - 1.
     Checked before the energy gate.
Everything else is v72 (and through it v71/v60).
"""
import argparse
from pathlib import Path
import re
import shutil
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
import rya1395_mcdhf_s0_v72 as V72  # noqa: E402

V71, S = V72.V71, V72.S
SCHEMA = 'rya1395.mcdhf-s0-v73'
METHOD1 = {'5s', '6s'}
L_OF = {'s': 0, 'p': 1, 'd': 2, 'f': 3}
_orig_rmcdhf_stdin = S.rmcdhf_stdin
ROW = re.compile(r'^\s+(\d+)([spdf])(-?)\s+[-\d.D+]+\s+(\d)\s+.*\s(\d+)\s*$')


def rmcdhf_stdin(levels, varied):
    if V71._current['label'] in METHOD1 and varied:
        v = ','.join(f'{o}*' for o in varied)
        return (f'n\nn\ny\nn\ny\n{S.GRID}y\n{S.ACCY:.17g}\n{levels}5\n{v}\n{v}\n{S.NSCF}\n'
                f'y\nn\ny\n{v}\n\n\n\n'
                'n\nn\nn\nn\nn\nn\n1\n')
    return _orig_rmcdhf_stdin(levels, varied)


S.rmcdhf_stdin = rmcdhf_stdin


def final_nodes(stdout_text):
    """{orbital label: (method, NNP)} from the LAST iteration table of rmcdhf.stdout."""
    last = {}
    for line in stdout_text.splitlines():
        m = ROW.match(line)
        if m:
            last[f'{m[1]}{m[2]}{m[3]}'] = (int(m[4]), int(m[5]))
    return last


def node_gate(out, receipt, parity, k):
    varied = S.STEPS[parity][k][2]
    if not varied:
        return
    rows = final_nodes((out / 'rmcdhf.stdout').read_text())
    checks = []
    for lab, (method, nnp) in rows.items():
        nl = lab.rstrip('-')
        if nl not in varied:
            continue
        want = int(nl[:-1]) - L_OF[nl[-1]] - 1
        checks.append(dict(orbital=lab, method=method, nnp=nnp, expected=want, pass_=nnp == want))
    case = receipt['steps'][out.name]
    case['node_gate'] = checks
    if not checks:
        raise ValueError(f'{out.name}: no iteration rows found for varied orbitals {varied}')
    bad = [c for c in checks if not c['pass_']]
    if bad:
        raise ValueError(f'{out.name}: node gate failed: ' + ', '.join(f"{c['orbital']} NNP {c['nnp']} != {c['expected']}" for c in bad))


def step(a, receipt, parity, k, previous):
    lab = S.STEPS[parity][k][0]
    out = V72._v71_step(a, receipt, parity, k, previous)
    case = receipt['steps'][out.name]
    if lab in V72.HYD_STEPS:
        case['start_zeff'] = V72.start_zeff(lab)
    case['method'] = 1 if lab in METHOD1 else 3
    node_gate(out, receipt, parity, k)
    V72.gate(out, receipt, parity, k)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--grasp-bin', type=Path, required=True)
    ap.add_argument('--strict-bin', type=Path, required=True)
    ap.add_argument('--v56-receipt', type=Path, required=True)
    ap.add_argument('--parity', choices=['odd'], required=True)
    ap.add_argument('--nist-table', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    a = ap.parse_args()
    if V71.platform.node().lower().split('.')[0] != 'sirius':
        raise RuntimeError('Atomic solvers belong on Sirius')
    for key, v in vars(a).items():
        if isinstance(v, Path):
            setattr(a, key, v.resolve())
    V72.BINDING = V72.nist_binding(a.nist_table)
    a.out.mkdir(parents=True, exist_ok=False)
    for f in (__file__, V72.__file__, V71.__file__, V72.HYD.__file__, S.__file__, a.nist_table):
        shutil.copyfile(f, a.out / Path(f).name)
    V71.SCHEMA = SCHEMA
    V71.step = step
    V71.run_parity(a)


if __name__ == '__main__':
    main()
