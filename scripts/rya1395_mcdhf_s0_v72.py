"""RYA-1395 v72: O I s0 stage with our own Z_eff=1 hydrogenic nf/nd orbitals and an orbital-energy gate (Sirius).

Declared on RYA-1395 (comment 01273db5) before any build or run. Extends the v71 driver (which reuses the v60
sequential driver unchanged):
  A. nf (even, FIXED) and nd (odd, VARIED from this start) orbitals come from rya1395_hydrogenic_rwfn (Z_eff = 1,
     Pauli small component) written to hyd.w in the step directory and read with rwfnestimate option 1, instead of
     GRASP's screened-hydrogenic option (bare Z=8 for our core-free lists: v71 FAILED_RETAINED).
  D. Gate after every step: each NEW spectroscopic orbital (n >= 3) must have eps within 30% of the NIST binding
     of the lowest targeted 4S-parent level of its configuration (Li 2023 Table A1 NIST energies; O I ionization
     limit 109837.02 cm-1, NIST ASD). Every orbital's eps and <r> are reported per step.
Everything else (half grid, strict criterion, method 3, identity and final-list guards) is v71/v60.
"""
import argparse
import csv
from pathlib import Path
import shutil
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
import rya1395_mcdhf_s0_v71 as V71  # noqa: E402  (patches the v60 module: nf steps vary nothing)
import rya1395_hydrogenic_rwfn as HYD  # noqa: E402

S = V71.S
SCHEMA = 'rya1395.mcdhf-s0-v72'
HYD_STEPS = {'4f', '5f', '3d', '4d', '5d', '5s', '6s'}
FIXED_ZEFF1 = {'4f', '5f'}   # nf: FIXED at Z_eff = 1. Others (amendment C', comment 5a7f4728): VARIED from Z_eff = n/n*_NIST
GATE_REL = 0.30
IP_CM = 109837.02          # NIST ASD O I ionization energy (cm-1)
HARTREE_CM = 219474.6313632
CORE = {'1s', '2s', '2p'}
TABLE_A1 = Path(__file__).resolve().parent.parent / 'data/reference/li2023_oi/raw/tablea1.tsv'


def nist_binding(table=TABLE_A1):
    """{'3s': Ha, ...}: binding of the lowest NIST level of each 2s2 2p3 nl configuration in Li Table A1."""
    low = {}
    with open(table) as f:
        rows = [r for r in csv.reader((l for l in f if not l.startswith('#')), delimiter='\t')]
    for r in rows[3:]:
        if len(r) < 8 or '2p(3).' not in r[2]:
            continue
        nl = r[2].strip().split('.')[-1]
        e = float(r[7])
        low[nl] = min(low.get(nl, 1e99), e)
    return {nl: (IP_CM - e) / HARTREE_CM for nl, e in low.items()}


BINDING = None
_orig_timed = S.timed


def _timed(exe, stdin, cwd, *args, **kw):
    lab = V71._current['label']
    if Path(str(exe)).name == 'rwfnestimate' and lab in HYD_STEPS:
        HYD.write_rwfn(Path(cwd) / 'hyd.w', HYD.build([lab], S.RNT, S.H, S.N, zeff=start_zeff(lab)))
    return _orig_timed(exe, stdin, cwd, *args, **kw)


S.timed = _timed


def start_zeff(lab):
    if lab in FIXED_ZEFF1:
        return 1.0
    n_star = (0.5 / BINDING[lab]) ** 0.5
    return int(lab[:-1]) / n_star


def hyd_estimate(n_previous):
    """Frozen orbitals from every earlier step (newest first), then our hydrogenic file for the new orbital."""
    return f'n\nn\n\ny\nn\ny\n{S.GRID}' + ''.join(f'1\nprev{i}.w\n*\n' for i in range(n_previous)) + '1\nhyd.w\n*\nn\n'


def _estimate(n_previous):
    lab = V71._current['label']
    return hyd_estimate(n_previous) if lab in HYD_STEPS else V71._orig_estimate(n_previous)


S.step_estimate = _estimate


def gate(out, receipt, parity, k):
    lab, _, varied = S.STEPS[parity][k]
    orbs = HYD.read_rwfn(out / 'rwfn.out')
    case = receipt['steps'][out.name]
    case['orbitals'] = [dict(orbital=HYD.label(o), eps=o['eps'], mean_r=HYD.mean_r(o)) for o in orbs]
    new = set(varied) | ({lab} if lab in HYD_STEPS else set())
    checks = []
    for o in orbs:
        nl = HYD.label(o).rstrip('-')
        if nl in CORE or nl not in new:
            continue
        b = BINDING[nl]
        rel = (o['eps'] - b) / b
        checks.append(dict(orbital=HYD.label(o), eps=o['eps'], nist_binding=b, rel=rel, pass_=abs(rel) <= GATE_REL))
    case['energy_gate'] = checks
    bad = [c for c in checks if not c['pass_']]
    if bad:
        raise ValueError(f'{out.name}: orbital-energy gate failed: ' +
                         ', '.join(f"{c['orbital']} eps {c['eps']:.5f} vs NIST {c['nist_binding']:.5f}" for c in bad))


_v71_step = V71.step


def step(a, receipt, parity, k, previous):
    lab = S.STEPS[parity][k][0]
    out = _v71_step(a, receipt, parity, k, previous)
    if lab in HYD_STEPS:
        receipt['steps'][out.name]['start_zeff'] = start_zeff(lab)
    gate(out, receipt, parity, k)
    return out


def main():
    global BINDING
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--grasp-bin', type=Path, required=True)
    ap.add_argument('--strict-bin', type=Path, required=True)
    ap.add_argument('--v56-receipt', type=Path, required=True)
    ap.add_argument('--parity', choices=['even', 'odd'], required=True)
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--nist-table', type=Path, default=TABLE_A1, help='Li 2023 Table A1 (NIST energies)')
    a = ap.parse_args()
    if V71.platform.node().lower().split('.')[0] != 'sirius':
        raise RuntimeError('Atomic solvers belong on Sirius')
    for key, v in vars(a).items():
        if isinstance(v, Path):
            setattr(a, key, v.resolve())
    BINDING = nist_binding(a.nist_table)
    a.out.mkdir(parents=True, exist_ok=False)
    for f in (__file__, V71.__file__, HYD.__file__, S.__file__, a.nist_table):
        shutil.copyfile(f, a.out / Path(f).name)
    V71.SCHEMA = SCHEMA
    V71.step = step
    V71.run_parity(a)


if __name__ == '__main__':
    main()
