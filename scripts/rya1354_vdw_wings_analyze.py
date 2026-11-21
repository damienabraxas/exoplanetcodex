#!/usr/bin/env python3
"""RYA-1354 v70c analysis: error budget, harness gate, control-differenced tests, ABO vs Unsoeld.

Reads data/audit/rya1354/mg1_vdw_wings_v70c/fit_*_{MARCS,ATLAS9}.json written by
scripts/rya1354_vdw_wings.py and writes summary.csv / summary.json there.

Rules (declared on RYA-1354, comment bdf5d990):
  sigma(d) = curvature (+) |dd/dA|*0.03 (+) |dd/dA|*sigma_gf (+) |d_MARCS - d_ATLAS9|; MARCS is primary.
  Harness gate: lines in one multiplet (one ABO sigma) agree within 2 sigma (sigma of the difference
    WITHOUT the shared A term, which cancels inside a multiplet).
  Delta_ctrl = inverse-variance mean of the control lines' d (Na D, Ca triplet).
  T1/T4 pass: |d| <= 2 sigma raw AND |d - Delta_ctrl| <= 2 sqrt(sigma^2 + sigma_ctrl^2).
  T2: chi2 at the Unsoeld point (model <r2>) vs at ABO (d=0), from the chi2 curve, interpolated.
"""
import csv
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / 'data/audit/rya1354/mg1_vdw_wings_v70c'

# log gamma6 offsets (dex) of alternatives relative to ABO, from data/audit/rya1354/mg1_model_v70
# vdw_vs_vald.csv / vdw_lines.csv (Unsoeld with our model <r2>; Unsoeld hydrogenic; n*_r ABO; VALD).
ALT = {
    'Mg b':     dict(unsold_model=-7.6397 + 7.2643, unsold_hyd=-7.6282 + 7.2643, abo_r=-7.2538 + 7.2643, vald=-7.267 + 7.2643),
    'Mg 5528':  dict(unsold_model=-7.3836 + 6.9776, unsold_hyd=-7.3725 + 6.9776, abo_r=-6.9923 + 6.9776, vald=-6.979 + 6.9776),
    'Mg 8806':  dict(unsold_model=-7.7411 + 7.4103, unsold_hyd=-7.7243 + 7.4103, abo_r=-7.4180 + 7.4103, vald=-7.413 + 7.4103),
    'Mg 4s-4p': dict(unsold_model=-7.6146 + 7.1476, unsold_hyd=-7.6259 + 7.1476, abo_r=-7.1291 + 7.1476,
                     list_unsold_x2p5=-7.6259 + np.log10(2.5) + 7.1476),
    'Mg 4p-4d': dict(unsold_model=-7.3750 + 6.9227, unsold_hyd=-7.3699 + 6.9227, abo_r=-6.9255 + 6.9227,
                     list_unsold_x2p5=-7.3699 + np.log10(2.5) + 6.9227),
}
SIG_A = 0.03


def load():
    rows = {}
    for f in sorted(OUT.glob('fit_*.json')):
        r = json.loads(f.read_text())
        rows.setdefault(r['name'], {})[r['atm'].split('.')[0]] = r
    return rows


def chi2_at(r, d):
    return float(np.interp(d, r['d_grid'], r['chi2']))


def main():
    rows = load()
    out = []
    for name, by in rows.items():
        m = by.get('MARCS')
        if m is None or 'd_fit' not in m:
            out.append(dict(name=name, verdict=(m or {}).get('verdict', 'MISSING')))
            continue
        a9 = by.get('ATLAS9', {})
        s_atm = abs(m['d_fit'] - a9['d_fit']) if 'd_fit' in a9 else float('nan')
        s_noA = float(np.sqrt(m['sig_curv'] ** 2 + m['sig_gf'] ** 2 + np.nan_to_num(s_atm) ** 2))
        s_tot = float(np.hypot(s_noA, m['sig_A']))
        rec = dict(name=name, role=m['role'], group=m['group'], atlas=m['atlas'], n_fit=m['n_fit'],
                   d_marcs=m['d_fit'], d_atlas9=a9.get('d_fit'), edge=m['at_grid_edge'] or a9.get('at_grid_edge', False),
                   sig_curv=m['sig_curv'], sig_A=m['sig_A'], sig_gf=m['sig_gf'], sig_atm=s_atm,
                   sig_noA=s_noA, sigma=s_tot, dd_dA=m['dd_dA'], rms=m['rms'])
        # T2 / alternatives: chi2 relative to the ABO point (scaled to reduced chi2 = 1)
        s2 = min(m['chi2']) / max(m['n_fit'] - 2, 1)
        for k, dv in ALT.get(m['group'], {}).items():
            rec[f'{k}_offset'] = dv
            rec[f'dchi2_{k}'] = (chi2_at(m, dv) - chi2_at(m, 0.0)) / s2
            rec[f'{k}_z'] = (m['d_fit'] - dv) / s_tot
        out.append(rec)

    ok = [r for r in out if 'd_marcs' in r]
    ctrl = [r for r in ok if r['role'] == 'control']
    w = np.array([1 / r['sigma'] ** 2 for r in ctrl])
    dc = float(np.sum(w * [r['d_marcs'] for r in ctrl]) / w.sum()) if ctrl else float('nan')
    sc = float(1 / np.sqrt(w.sum())) if ctrl else float('nan')
    sc_scatter = float(np.std([r['d_marcs'] for r in ctrl], ddof=1)) if len(ctrl) > 1 else float('nan')
    sc_use = max(sc, sc_scatter) if np.isfinite(sc_scatter) else sc

    # harness gate: within-multiplet agreement
    gate = []
    for g in sorted({r['group'] for r in ok}):
        mem = [r for r in ok if r['group'] == g]
        for i in range(len(mem)):
            for j in range(i + 1, len(mem)):
                a, b = mem[i], mem[j]
                z = (a['d_marcs'] - b['d_marcs']) / np.hypot(a['sig_noA'], b['sig_noA'])
                gate.append(dict(group=g, a=a['name'], b=b['name'], diff=a['d_marcs'] - b['d_marcs'],
                                 z=float(z), pass_=bool(abs(z) <= 2)))

    for r in ok:
        r['z_raw'] = r['d_marcs'] / r['sigma']
        r['d_minus_ctrl'] = r['d_marcs'] - dc
        r['z_ctrl'] = r['d_minus_ctrl'] / np.hypot(r['sigma'], sc_use)
        if r['role'] in ('T1', 'T4'):
            r['verdict'] = 'PASS' if (abs(r['z_raw']) <= 2 and abs(r['z_ctrl']) <= 2) else 'FAIL'
        elif r['role'] == 'control':
            r['verdict'] = 'CONTROL'
        else:
            r['verdict'] = 'REPORT'
        if r['edge']:
            r['verdict'] += '(GRID_EDGE)'

    keys = sorted({k for r in out for k in r}, key=lambda k: (k not in ('name', 'role', 'group', 'verdict'), k))
    with open(OUT / 'summary.csv', 'w', newline='') as f:
        wr = csv.DictWriter(f, fieldnames=keys)
        wr.writeheader()
        for r in out:
            wr.writerow({k: (f'{v:.4f}' if isinstance(v, float) else v) for k, v in r.items()})
    summ = dict(delta_ctrl=dc, sigma_ctrl_ivw=sc, sigma_ctrl_scatter=sc_scatter, sigma_ctrl_used=sc_use,
                gate=gate, gate_pass=all(g['pass_'] for g in gate), rows=out)
    (OUT / 'summary.json').write_text(json.dumps(summ, indent=1, default=float))

    print(f"Delta_ctrl = {dc:+.3f} (ivw sigma {sc:.3f}, scatter {sc_scatter:.3f})")
    print('GATE (within-multiplet):')
    for g in gate:
        print(f"  {g['group']:10s} {g['a']:8s}-{g['b']:8s} {g['diff']:+.3f}  z {g['z']:+.2f}  {'PASS' if g['pass_'] else 'FAIL'}")
    print(f"{'line':9s} {'role':7s} {'d_MARCS':>8s} {'d_AT9':>7s} {'sigma':>6s} {'curv':>5s} {'A':>5s} {'gf':>5s} {'atm':>5s} {'z':>6s} {'z_ctrl':>6s}  verdict  | dchi2 Unsold(model)")
    for r in ok:
        print(f"{r['name']:9s} {r['role']:7s} {r['d_marcs']:+8.3f} {r['d_atlas9'] if r['d_atlas9'] is None else format(r['d_atlas9'], '+7.3f')} "
              f"{r['sigma']:6.3f} {r['sig_curv']:5.3f} {r['sig_A']:5.3f} {r['sig_gf']:5.3f} {r['sig_atm']:5.3f} "
              f"{r['z_raw']:+6.2f} {r['z_ctrl']:+6.2f}  {r['verdict']:8s} | {r.get('dchi2_unsold_model', float('nan')):8.0f}")
    for r in out:
        if 'd_marcs' not in r:
            print(f"{r['name']:9s} {r['verdict']}")


if __name__ == '__main__':
    main()
