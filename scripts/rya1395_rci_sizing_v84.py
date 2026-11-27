"""RYA-1395 v84: size the weight-selected RCI on the converged n5 model (CSF lists and cost only).

Declared on RYA-1395 before launch. No solver runs; no energy, rate or grade comes from this stage.

1. MR selection (Li 2023 Sect. 2; Li N I uses weight > 0.05). For every target level of v77 n5 (even 45,
   odd 36), the configuration weight is the sum of the LSJ component weights in state.lsj.lbl that share a
   configuration (configuration weights are basis-invariant, jj = LSJ). A configuration with weight > 0.05
   in any level of its parity enters the MR. jj2lsj prints components above a cutoff only, so for each level
   the unprinted mass m = 1 - Comp is an upper bound on what is missing: a configuration with
   w <= 0.05 < w + m is AMBIGUOUS and is included (conservative), and listed.
2. RCI list (Li 2023 RCI scheme, as sized in v47): 1s closed; 2s opened; SD from every n=2 subshell and the
   outer electron into the v77 n5 active set of that parity (even 5s,6p,5d,5f,5g; odd 6s,5p,5d,5f,5g); the
   even MR is extended with 1s2 2s 2p4 ns/nd (orbitals present in the even n5 set); 2J range as v77.
3. rcsfinteract, Dirac-Coulomb mode (the v77 MCDHF and the planned RCI are Dirac-Coulomb), against the MR.
   Every reduced list is inventoried and checked as a subset of the full list containing every MR CSF.
4. Cost per J block, extrapolated from the v52 measured block (odd J=1, n5, 33 575 CSFs: rangular 282 s,
   rci 495 s for one eigenvalue, angular 1.12 GB) with the v52 measured exponents (time N^2.2, disk N^1.56).
   These are ESTIMATES, labelled as such; the RCI stage is declared separately with them.

Any refusal or failed check stops the run: FAILED_RETAINED.
"""
import argparse
from collections import defaultdict
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import re
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
import rya1354_rci_feasibility as F  # noqa: E402

SCHEMA = 'rya1395.rci-sizing-v84'
THRESHOLD = 0.05
ACTIVE = {'even': '5s,6p,5d,5f,5g', 'odd': '6s,5p,5d,5f,5g'}
TWOJ = {'even': '0,10', 'odd': '0,8'}
EVEN_2S_FAMILIES = {'s': range(3, 6), 'd': range(3, 6)}   # capped at the even n5 set (s <= 5s, d <= 5d)
V52 = dict(csfs=33575, rangular_s=282.0, rci_s=495.0, angular_gb=1.12, t_exp=2.2, disk_exp=1.56)
L = 'spdfgh'
TOKEN = re.compile(r'^(\d+)([spdfgh])(?:\((\d+)\))?')


def configuration(label):
    occ = {}
    for tok in label.split('.'):   # intermediate terms ride on tokens as '3p_2D'; the orbital leads each token
        m = TOKEN.match(tok)
        if not m:
            raise ValueError(f'unparsed LSJ token {tok!r} in {label!r}')
        occ[m[1] + m[2]] = occ.get(m[1] + m[2], 0) + int(m[3] or 1)
    if sum(occ.values()) != 8 or occ.get('1s') != 2:
        raise ValueError(f'not an 8-electron 1s2 configuration: {label!r}')
    return tuple(sorted(occ.items(), key=lambda kv: (int(kv[0][:-1]), L.index(kv[0][-1]))))


def levels(lbl):
    out, cur = [], None
    for line in Path(lbl).read_text().splitlines()[1:]:
        f = line.split()
        if len(f) == 5 and f[4].endswith('%'):
            cur = dict(pos=int(f[0]), J=f[1], energy=float(f[3]), comp=float(f[4][:-1]) / 100, w=defaultdict(float))
            out.append(cur)
        elif len(f) == 3:
            cur['w'][configuration(f[2])] += float(f[1])
        elif f:
            raise ValueError(f'unexpected lbl line {line!r}')
    return out


def select(lbl):
    chosen, ambiguous = {}, {}
    lv = levels(lbl)
    for lev in lv:
        missing = max(0.0, 1 - lev['comp'])
        for cfg, w in lev['w'].items():
            if w > THRESHOLD:
                chosen[cfg] = max(chosen.get(cfg, 0), w)
            elif w + missing > THRESHOLD:
                ambiguous[cfg] = max(ambiguous.get(cfg, 0), w)
    ambiguous = {c: w for c, w in ambiguous.items() if c not in chosen}
    return lv, chosen, ambiguous


def ref_line(cfg, open2s=True):
    parts = []
    for orb, k in cfg:
        parts.append(f'{orb}({k},{"i" if orb == "1s" else "*"})')
    return ''.join(parts)


def show(cfg):
    return ' '.join(f'{o}{k}' if k > 1 else o for o, k in cfg if o != '1s')


def stdin(refs, parity, excitations):
    return '*\n0\n' + '\n'.join(refs) + '\n\n' + ACTIVE[parity] + '\n' + TWOJ[parity] + f'\n{excitations}\nn\n'


def cost(n):
    r = n / V52['csfs']
    return dict(csfs=n, rangular_h=round(V52['rangular_s'] * r ** V52['t_exp'] / 3600, 2),
                rci_h_one_eigenvalue=round(V52['rci_s'] * r ** V52['t_exp'] / 3600, 2),
                angular_gb=round(V52['angular_gb'] * r ** V52['disk_exp'], 2))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--grasp-bin', type=Path, required=True)
    ap.add_argument('--v77', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    a = ap.parse_args()
    a.out.mkdir()
    receipt = dict(schema=SCHEMA, started=datetime.now(timezone.utc).isoformat(), status='RUNNING',
                   threshold=THRESHOLD, active=ACTIVE, twoj=TWOJ, v52_scaling=V52, estimates_are_not_measurements=True)
    rpath = a.out / 'receipt.json'
    rpath.write_text(json.dumps(receipt, indent=1))
    try:
        for parity in ('even', 'odd'):
            lbl = a.v77 / f'n5_{parity}' / 'state.lsj.lbl'
            lv, chosen, ambiguous = select(lbl)
            mr = sorted(set(chosen) | set(ambiguous), key=lambda c: [(int(o[:-1]), L.index(o[-1]), k) for o, k in c])
            refs = [ref_line(c) for c in mr]
            fams = []
            if parity == 'even':
                for l, ns in EVEN_2S_FAMILIES.items():
                    fams += [f'1s(2,i)2s(1,*)2p(4,*){n}{l}(1,*)' for n in ns]
            rec = dict(levels=len(lv), min_comp=min(x['comp'] for x in lv),
                       mr_selected={show(c): round(w, 4) for c, w in chosen.items()},
                       mr_ambiguous_included={show(c): round(w, 4) for c, w in ambiguous.items()},
                       references=refs + fams)
            base = a.out / parity
            base.mkdir()
            full_dir, mr_dir = base / 'full', base / 'mr'
            full_dir.mkdir()
            mr_dir.mkdir()
            rec['rcsfgenerate_full'] = F.timed(a.grasp_bin / 'rcsfgenerate', stdin(refs + fams, parity, 2), full_dir,
                                               F.BUDGET_S['rcsfgenerate'])
            rec['rcsfgenerate_mr'] = F.timed(a.grasp_bin / 'rcsfgenerate', stdin(refs + fams, parity, 0), mr_dir,
                                             F.BUDGET_S['rcsfgenerate'])
            full, mrf = full_dir / 'rcsf.out', mr_dir / 'rcsf.out'
            inv_full = F.inventory_stream(full, parity, 'rci')
            inv_mr = F.inventory_stream(mrf, parity, 'rci')
            step, table, inv_red = F.reduce(a.grasp_bin, base, parity, 5, 'DC', full, mrf)
            F.check_reduction(table, inv_full, inv_red, inv_mr, full, base / 'DC' / 'rcsf.out', mrf)
            rec['rcsfinteract_DC'] = step
            rec['blocks'] = {k: dict(full=inv_full['blocks'][k], mr=inv_mr['blocks'].get(k, 0),
                                     dc_reduced=inv_red['blocks'].get(k, 0),
                                     targets=sum(1 for x in lv if x['J'] == k[:-1]),
                                     estimate=cost(inv_red['blocks'].get(k, 0)))
                             for k in inv_full['blocks']}
            rec['total_dc_reduced'] = sum(b['dc_reduced'] for b in rec['blocks'].values())
            receipt[parity] = rec
            rpath.write_text(json.dumps(receipt, indent=1))
        est = [b['estimate'] for p in ('even', 'odd') for b in receipt[p]['blocks'].values()]
        receipt['estimate_totals'] = dict(
            rangular_h=round(sum(e['rangular_h'] for e in est), 1),
            rci_h_one_eigenvalue=round(sum(e['rci_h_one_eigenvalue'] for e in est), 1),
            angular_gb_largest_block=max(e['angular_gb'] for e in est),
            angular_gb_all=round(sum(e['angular_gb'] for e in est), 1),
            largest_block=max(e['csfs'] for e in est))
        receipt['status'] = 'COMPLETE'
    except Exception as exc:  # noqa: BLE001 - every failure is retained
        receipt['status'] = 'FAILED_RETAINED'
        receipt['error'] = f'{type(exc).__name__}: {exc}'
        raise
    finally:
        receipt['finished'] = datetime.now(timezone.utc).isoformat()
        rpath.write_text(json.dumps(receipt, indent=1))


if __name__ == '__main__':
    main()
