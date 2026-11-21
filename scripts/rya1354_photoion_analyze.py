"""RYA-1354 v70b analysis: pool the box runs of rya1354_photoion.py into sigma(e) per bound Mg I term; declared
controls C1/C2 and tests T1-T5 (T4 read from alpha_T4.json); OP (Butler 1993) tables from the OpenMULTI3D atom.

Conventions fixed BEFORE the first result was read (2026-11-21):
  * energy coordinate = e above the Mg+ 3s threshold of the MODEL; lab photon energies are converted with NIST
    (e_lab = E_NIST(initial) + h nu - IP_NIST): thresholds aligned, as in NLTE practice (Lind & Amarsi 2024);
  * T1 E_r = IP_NIST + e_r(model) (primary); the model-absolute E_r (e_r - e_ground, model) is reported too;
  * sigma_model = hypot(box spread (std over R_box), |length - velocity| / 2);
  * only final symmetries with an open Mg+ 3s channel (parity == L % 2) carry a continuum below Mg+ 3p.
"""
import json
import math
import re
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

import rya1354_photoion as P
import rya1354_onevalence_model as M
import rya1354_twovalence_model as T

CM = None   # hartree in cm-1 (set in main after M.use)


def load_boxes():
    boxes = []
    for js in sorted(P.OUT.glob('box_*.json')):
        rb = int(js.stem.split('_')[1])
        if rb == 0:
            continue
        rec = json.load(open(js))
        arr = np.load(P.OUT / f'box_{rb}.npz')
        boxes.append((rb, rec, arr))
    return boxes


def open_channel(sym):
    L, S, par = sym
    return par == L % 2


def box_points(arr, key):
    """Box-normalised (e, df/de length, df/de velocity) above threshold for one initial|final-sym key."""
    e, f, fv = arr[key + '|e'], arr[key + '|f'], arr[key + '|fv']
    m = e > 0
    e, f, fv = e[m], f[m], fv[m]
    if len(e) < 3:
        return None
    o = np.argsort(e)
    e, f, fv = e[o], f[o], fv[o]
    gap = np.empty_like(e)
    gap[1:-1] = (e[2:] - e[:-2]) / 2
    gap[0], gap[-1] = e[1] - e[0], e[-1] - e[-2]
    return e[1:-1], (f / gap)[1:-1], (fv / gap)[1:-1]     # drop the edges (one-sided spacing)


def bound_points(arr, key):
    """df/de just below threshold from the unperturbed Rydberg members n* 7..BOX_N (v70 validity limit): f_n n*^3.
    Harness correction after the dry run (2026-11-21): states with n* > BOX_N are box-compressed (one of them collects
    the oscillator strength of every missing high-n member), so they cannot sample df/de."""
    e, f = arr[key + '|e'], arr[key + '|f']
    lo, hi = -1 / (2 * 7.0 ** 2), -1 / (2 * T.BOX_N ** 2)
    m = (e > lo) & (e < hi)
    return e[m], f[m] / (-2 * e[m]) ** 1.5


def sigma_curves(boxes, label, grid):
    """Per box: total sigma(grid) (length, velocity) in Mb, summed over open final syms (per-box linear interp)."""
    out = []
    for rb, rec, arr in boxes:
        ent = rec['inits'].get(label)
        if not ent:
            continue
        sl, sv = np.zeros_like(grid), np.zeros_like(grid)
        for sym, key in ent['finals'].items():
            if not open_channel(eval(sym)):
                continue
            p = box_points(arr, key)
            if p is None:
                continue
            e, dl, dv = p
            sl += np.interp(grid, e, dl, left=np.nan, right=np.nan) * P.MB_PER_DFDE
            sv += np.interp(grid, e, dv, left=np.nan, right=np.nan) * P.MB_PER_DFDE
        out.append((rb, sl, sv))
    return out


def pooled(boxes, label, sym):
    pts = []
    for rb, rec, arr in boxes:
        key = rec['inits'].get(label, {}).get('finals', {}).get(str(sym))
        if key:
            p = box_points(arr, key)
            if p:
                pts.append(np.column_stack(p))
    a = np.vstack(pts)
    return a[np.argsort(a[:, 0])]


def thr_extrap(boxes, label):
    """sigma at threshold (Mb, length): per open final sym, linear fit of pooled box points 0.002..0.012 hartree, at 0."""
    tot = 0.0
    for sym in boxes[0][1]['inits'][label]['finals']:
        if not open_channel(eval(sym)):
            continue
        pb = [box_points(a, r['inits'][label]['finals'][sym]) for _, r, a in boxes]
        pb = [p for p in pb if p]
        if not pb:
            continue
        ec = np.concatenate([p[0] for p in pb]); dc = np.concatenate([p[1] for p in pb])
        m = (ec > 0.002) & (ec < 0.012)
        if m.sum() >= 4:
            tot += np.polyval(np.polyfit(ec[m], dc[m], 1), 0.0) * P.MB_PER_DFDE
    return float(tot)


def at(curves, e0):
    L = np.array([np.interp(e0, GRID, c[1]) for c in curves])
    V = np.array([np.interp(e0, GRID, c[2]) for c in curves])
    m, mv = np.nanmean(L), np.nanmean(V)
    spread = np.nanstd(L)
    return dict(sigma_Mb=round(float(m), 2), sigma_vel_Mb=round(float(mv), 2), box_std_Mb=round(float(spread), 2),
                box_rel_spread=round(float(spread / m), 4) if m else None, per_box=[round(float(x), 2) for x in L],
                sigma_model_Mb=round(float(math.hypot(spread, abs(m - mv) / 2)), 2))


def fano_fit(pts, e0, half=2000):
    """sigma = a (q + x)^2 / (1 + x^2) + b, x = 2 (e - er) / G (e, er, G in cm-1, sigma in Mb), pooled box points."""
    e = pts[:, 0] * CM
    s = pts[:, 1] * P.MB_PER_DFDE
    m = np.abs(e - e0) < half
    e, s = e[m], s[m]
    i = np.argmax(s)

    def model(p, x):
        a, q, er, g, b = p
        y = 2 * (x - er) / g
        return a * (q + y) ** 2 / (1 + y ** 2) + b
    p0 = [s[i] / 400, 20.0, e[i], 300.0, np.median(s) * 0.5]
    r = least_squares(lambda p: model(p, e) - s, p0, bounds=([0, -200, e0 - half, 5, -1e4], [1e4, 200, e0 + half, 5000, 1e4]))
    a, q, er, g, b = r.x
    J = r.jac
    try:
        cov = np.linalg.inv(J.T @ J) * np.sum(r.fun ** 2) / max(1, len(e) - 5)
        err = np.sqrt(np.diag(cov))
    except np.linalg.LinAlgError:
        err = [np.nan] * 5
    return dict(er_cm=float(er), er_err=float(err[2]), gamma_cm=float(abs(g)), gamma_err=float(err[3]), q=float(q),
                q_err=float(err[1]), peak_Mb=float(a * (1 + q ** 2) + b), n_pts=int(len(e)),
                rms_rel=float(np.sqrt(np.mean((r.fun / np.maximum(s, 1e-9)) ** 2))))


def om_tables():
    """{our label: (lam_thr_A, [(lam_A, sigma_Mb)])} from the OpenMULTI3D Mg atom (OP / Butler 1993 data)."""
    lines = open(P.OM_ATOM).read().splitlines()
    lev, out = {}, {}
    for ln in lines:
        m = re.search(r"Level\s+(\d+)\s+=\s+(\d+[a-z])(\d)([A-Z])(\d?)(\*?)\s*'", ln)
        if m:
            nl, mult, L, star = m[2], m[3], m[4], m[6]
            term = f'{mult}{L}{star}'
            lab = '3s2 1S' if (nl == '3s' and term == '1S') else f'3s.{nl} {term}'
            lev[int(m[1])] = lab
    i = next(k for k, ln in enumerate(lines) if 'PHOTO-ION' in ln) + 1
    while i < len(lines):
        t = lines[i].split()
        if len(t) < 4 or lines[i].lstrip().startswith('*'):
            break
        lo, nq = int(t[1]), int(t[3])
        pts = [tuple(map(float, lines[i + 1 + k].split()[:2])) for k in range(nq)]
        lab = lev.get(lo)
        if lab and lab not in out:
            out[lab] = (pts[0][0], [(l, s * 1e18) for l, s in pts])
        i += 1 + nq
    return out


GRID = None


def main():
    global CM, GRID
    M.use('Mg II')
    CM = M.RY_CM
    GRID = np.linspace(1e-4, P.EMAX * 0.98, 4000)          # hartree above threshold
    boxes = load_boxes()
    if not boxes:
        raise SystemExit('no box runs')
    ip = T.SPECIES['Mg I']['ip']
    summ = dict(schema='rya1354_photoion_v70b', status=P.STATUS, boxes=[b[0] for b in boxes], npts=P.NPTS, emax_h=P.EMAX,
                channels=P.NCH, ncont=boxes[0][1]['ncont'])
    # C0 (read)
    c0 = {p.stem: json.load(open(p)) for p in P.OUT.glob('c0_*.json')}
    summ['C0'] = {k: dict(max_abs_d=v['max_abs_d'], n=v['n_lab_lines'], energy_max_rel=v['energy_max_rel'],
                          pass_=v['max_abs_d'] <= 0.01) for k, v in c0.items()}
    # every initial term: sigma curve, OP comparison
    om = om_tables()
    labels = sorted(boxes[0][1]['inits'], key=lambda l: boxes[0][1]['inits'][l]['e_nist'])
    rows, curves_all = [], {}
    for lab in labels:
        cv = sigma_curves(boxes, lab, GRID)
        if not cv:
            continue
        curves_all[lab] = cv
        L = np.nanmean([c[1] for c in cv], axis=0)
        V = np.nanmean([c[2] for c in cv], axis=0)
        sthr = thr_extrap(boxes, lab)
        a3 = at(cv, 0.003)
        r = dict(label=lab, sigma_thr_Mb=round(sthr, 3), sigma_660cm_Mb=a3['sigma_Mb'], sigma_660cm_vel_Mb=a3['sigma_vel_Mb'],
                 box_rel_spread_660cm=a3['box_rel_spread'])
        if lab in om:
            lthr, pts = om[lab]
            r['op_sigma_thr_Mb'] = round(pts[0][1], 3)
            r['op_ratio_thr'] = round(sthr / pts[0][1], 3) if pts[0][1] else None
            cmp = []
            for lam, s in pts:
                e = (1e8 / lam - 1e8 / lthr) / CM
                if 0 < e < GRID[-1]:
                    cmp.append((e * CM, s, float(np.interp(e, GRID, L))))
            r['op_points_in_range'] = len(cmp)
            if cmp:
                r['op_median_log_ratio'] = round(float(np.median([math.log10(c[2] / c[1]) for c in cmp if c[1] > 0 and c[2] > 0])), 3)
        rows.append(r)
    summ['levels'] = rows
    # sigma curves are regenerable from box_*.npz (39 MB); not written (v70b)
    # 3s3p 1P tests
    e1p = next(v['e_nist'] for v in [boxes[0][1]['inits'][P.INIT_1P]])
    E1p = ip - (-(M.IP_CM) / CM - e1p) * CM          # NIST term energy (cm-1 above Mg ground)
    tests = {}
    cv = curves_all[P.INIT_1P]
    for t in ('T2', 'T3'):
        lab = P.LAB[t]
        nu = 1e8 / T.air_to_vac(lab['lam_nm'] * 10)
        e_lab = (E1p + nu - ip) / CM
        a = at(cv, e_lab)
        z = (a['sigma_Mb'] - lab['sigma_Mb']) / math.hypot(lab['err'], a['sigma_model_Mb'])
        tests[t] = dict(lab, e_above_thr_cm=round(e_lab * CM, 1), **a, z=round(z, 2), pass_=abs(z) <= 2,
                        C1_pass=(a['box_rel_spread'] or 1) <= 0.10)
    # T1 resonance: pooled 1S final points of 3s3p 1P
    lab = P.LAB['T1']
    pts = pooled(boxes, P.INIT_1P, (0, 0, 0))
    fit = fano_fit(pts, lab['Er'] - ip)
    e_ground_model = boxes[0][1]['inits']['3s2 1S']['e_model']
    thr_model = boxes[0][1]['e_thr']
    er_abs = fit['er_cm'] + ip
    er_model_abs = (thr_model + fit['er_cm'] / CM - e_ground_model) * CM
    tests['T1'] = dict(lab, fit=fit, Er_threshold_aligned=round(er_abs, 1), Er_model_absolute=round(er_model_abs, 1),
                       pass_Er=abs(er_abs - lab['Er']) <= 300, pass_Gamma=1 / 1.5 <= fit['gamma_cm'] / lab['Gamma'] <= 1.5,
                       pass_peak=abs(fit['peak_Mb'] - lab['peak_Mb']) <= 2 * lab['peak_err'])
    # T5 continuity at threshold, per open final sym, every initial term with enough Rydberg members
    t5 = []
    for labl in labels:
        if boxes[0][1]['inits'][labl]['e_model'] - boxes[0][1]['e_thr'] > -1 / (2 * 6.0 ** 2):
            continue      # initial inside/above the Rydberg window (n* > 6): no series sample (harness correction)
        for rb, rec, arr in boxes[:1]:
            ent = rec['inits'].get(labl)
            for sym, key in ent['finals'].items():
                if not open_channel(eval(sym)):
                    continue
                eb, db = bound_points(arr, key)
                pb = [box_points(a2, r2['inits'][labl]['finals'][sym]) for _, r2, a2 in boxes]
                ec = np.concatenate([p[0] for p in pb if p]); dc = np.concatenate([p[1] for p in pb if p])
                mc = (ec > 0.002) & (ec < 0.012)          # skip the threshold box states (harness correction)
                if len(eb) < 3 or mc.sum() < 4:
                    continue
                below = np.polyval(np.polyfit(eb, db, 1), 0.0)
                above = np.polyval(np.polyfit(ec[mc], dc[mc], 1), 0.0)
                if max(abs(below), abs(above)) * P.MB_PER_DFDE < 0.01:
                    continue
                t5.append(dict(init=labl, final=sym, below_Mb=round(below * P.MB_PER_DFDE, 4), above_Mb=round(above * P.MB_PER_DFDE, 4),
                               ratio=round(above / below, 3) if below else None, n_bound=len(eb), n_cont=int(mc.sum()),
                               pass_=bool(below != 0 and abs(above / below - 1) <= 0.15)))
    tests['T5'] = dict(rows=t5, n=len(t5), n_pass=sum(r['pass_'] for r in t5))
    al = P.OUT / 'alpha_T4.json'
    if al.exists():
        tests['T4'] = json.load(open(al))
    summ['tests'] = tests
    (P.OUT / 'summary.json').write_text(json.dumps(summ, indent=1, default=str) + '\n')
    with open(P.OUT / 'levels_vs_op.csv', 'w') as fh:
        keys = sorted({k for r in rows for k in r})
        fh.write(','.join(keys) + '\n')
        for r in rows:
            fh.write(','.join(str(r.get(k, '')) for k in keys) + '\n')
    print(json.dumps({k: summ[k] for k in ('C0',)}, indent=1))
    print(json.dumps({k: v for k, v in tests.items() if k != 'T5'}, indent=1, default=str)[:4000])
    print('T5', tests['T5']['n_pass'], '/', tests['T5']['n'])
