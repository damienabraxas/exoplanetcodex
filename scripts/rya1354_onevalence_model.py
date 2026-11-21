"""RYA-1354: our own atomic model for ONE-VALENCE-ELECTRON spectra (Na I v62, Mg II v63; K I / Ca II next),
model potential with core polarization, scored against experiment and independent theory, then used to fill every
line of the species in the solar line list with a log gf and a graded sigma. Each species was declared on RYA-1354
before its run.

One valence electron in an l-dependent parametric potential (form of Marinescu, Sadeghpour & Dalgarno 1994,
PRA 49, 982) with core polarization,
    V_l(r) = -Z_l(r)/r - alpha_d/(2 r^4) (1 - exp(-(r/rc_l)^6)),
    Z_l   = zc + (Z - zc) exp(-a1 r) - r (a3 + a4 r) exp(-a2 r)      (zc = ion charge + 1),
the four a_l re-fitted HERE to the NIST term energies (n <= 9); alpha_d = static dipole polarizability of the core.
Dipole operator with the core-polarization correction  d(r) = r [1 - alpha_d/r^3 (1 - exp(-(r/rc)^3))].
gf per fine-structure line from the LS 6j fraction on NIST experimental energies.

Per-line indicators: radial cancellation CF = |int u1 d u2| / int |u1 d u2|; half the core-polarization-operator
effect (the model's analogue of a gauge difference) -> candidate sigma and class (see line_sigma).
Referees: experimental resonance-line A-values; Froese Fischer [9] B-spline CI / [16] MCHF-CV (NIST MCHF collection).
Status: CANDIDATE / HOLD (RYA-1352). Nothing here writes into a product.

Usage: rya1354_onevalence_model.py [SPECIES ...]   (default: all)
"""
import csv
import re
import json
import math
import statistics as st
from fractions import Fraction as Fr
from functools import lru_cache
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

import rya1354_na_model as V61

ROOT = Path(__file__).resolve().parents[1]
REF = ROOT / 'data/reference'
# starting values, Marinescu+1994 Table I (Na): a1, a2, a3, a4, rc (f takes rc of d)
START_NA = {0: (4.82223117, 2.45449865, -1.12255048, -1.42631393, 0.45489422),
            1: (5.08382502, 2.18226881, -1.19534623, -1.03142861, 0.45798739),
            2: (3.53324124, 2.48697936, -0.75688448, -1.27852357, 0.71875312),
            3: (1.11056646, 1.05458759, 1.73203428, -0.09265696, 0.71875312)}
R_INF, ME_AMU = 109737.31568, 5.48579909e-4
SPECIES = {
    'Na I': dict(Z=11, zc=1, alpha=0.9448, mass=22.98977, ip=41449.451, sigma=-0.7, start=START_NA,
                 levels=REF / 'na_mchf_ff2002/raw/nist_asd_na1_levels.csv',
                 ff9=REF / 'na_mchf_ff2002/raw/Na_11.9.U.lsj.db.gz', ff16=REF / 'na_mchf_ff2002/raw/Na_11.16.mchf_cv.lsj.db.gz',
                 exp={('3s', 1, '3p', 3): (6.16e7, 0.005), ('3s', 1, '3p', 1): (6.14e7, 0.005)},   # KP08 Table 4, AA
                 linelist=('Na', 'I'), out='na_model_v62', asplund=(4752, 5149, 6154, 6161, 10746), core='2p6'),
    # Mg II: alpha_d(Mg2+) 0.4814 a.u. (Safronova-type MBPT value); anchor = KP08 Table 44 3s-3p A (class A,
    # FF/Siegel/Ansbacher-lifetime/Johnson); start = the Na shape, re-fitted.
    'Mg II': dict(Z=12, zc=2, alpha=0.4814, mass=24.305, ip=121267.61, sigma=-2.2, start=START_NA,
                  levels=REF / 'onevalence_nist/MgII_levels.csv',
                  ff9=REF / 'onevalence_nist/mchf_mg2/Na_12.9.U.lsj.db.gz', ff16=REF / 'onevalence_nist/mchf_mg2/Na_12.16.mchf_cv.lsj.db.gz',
                  exp={('3s', 1, '3p', 3): (2.60e8, 0.03), ('3s', 1, '3p', 1): (2.57e8, 0.03)},
                  linelist=('Mg', 'II'), out='mg2_model_v63', asplund=(4481, 7877, 7896, 9218, 9244, 10915, 10952), core='2p6'),
    # Ca II (v64): alpha_d(Ca2+) 3.26(3) a.u. (Safronova & Safronova 2011 Sect. VI); anchors = Hettrich+2015 single-ion
    # rates (H line 2pi x 21.57(8) MHz; 8662 A 3d3/2-4p1/2 2pi x 1.482(8) MHz); referee = Safronova 2011 Table V
    # (all-order relativistic, 55 lines with stated uncertainty) in place of an FF collection (none exists for Ca II).
    'Ca II': dict(Z=20, zc=2, alpha=3.26, mass=40.078, ip=95751.87, sigma=-1.8, start=START_NA,
                  levels=REF / 'onevalence_nist/CaII_levels.csv',
                  ff9=REF / 'onevalence_nist/safronova2011_ca2/table5_ca2_A_f.csv', ff16=None,
                  exp={('4s', 1, '4p', 1): (2 * math.pi * 21.57e6, 0.0037), ('3d', 3, '4p', 1): (2 * math.pi * 1.482e6, 0.0054)},
                  linelist=('Ca', 'II'), out='ca2_model_v64', asplund=(5001, 6457, 7324, 8249, 8255), core='3p6'),
    # v69 sister-ion cores (Mg II isoelectronic). alpha_d: Al3+ 0.2649, Si4+ 0.1624 a.u. (Mitroy, Safronova & Clark 2010
    # J. Phys. B 43 202001, Ne-like ions). Core fit only (rya1354_core_fit.py); no line referee at this level.
    'Al III': dict(Z=13, zc=3, alpha=0.2649, mass=26.98154, ip=229445.71, sigma=-4.2, start=START_NA,
                   levels=REF / 'onevalence_nist/AlIII_levels.csv', ff9=None, ff16=None, exp={}, linelist=('Al', 'III'),
                   out='al3_core_v69', asplund=(), core='2p6'),
    'Si IV': dict(Z=14, zc=4, alpha=0.1624, mass=28.0855, ip=364093.1, sigma=-6.6, start=START_NA,
                  levels=REF / 'onevalence_nist/SiIV_levels.csv', ff9=None, ff16=None, exp={}, linelist=('Si', 'IV'),
                  out='si4_core_v69', asplund=(), core='2p6'),
}
SP = None   # the active species config (set by use())
NMAX = 15
X0, X1, NPTS = math.log(1e-5), math.log(1200.0), 8000
STATUS = 'CANDIDATE_HOLD_NOT_ADOPTED'


def use(name):
    """Activate one species: module-level constants the solver reads."""
    global SP, Z, ZC, ALPHA_D, IP_CM, RY_CM, START, OUT
    SP = dict(SPECIES[name], name=name)
    Z, ZC, ALPHA_D, IP_CM, START = SP['Z'], SP['zc'], SP['alpha'], SP['ip'], SP['start']
    RY_CM = 2 * R_INF / (1 + ME_AMU / SP['mass'])      # hartree in cm-1, reduced mass
    OUT = ROOT / 'data/audit/rya1354' / SP['out']
    return SP


def nist_levels():
    lev = {}
    for r in csv.DictReader(open(SP['levels'])):
        clean = {k: (v or '').strip('="[]') for k, v in r.items() if k}
        m = re.fullmatch(re.escape(SP['core']) + r'\.(\d+)([spdf])', clean['Configuration'])
        if not m or not re.fullmatch(r'\d+/2', clean['J']) or not clean['Level (cm-1)']:
            continue
        lev[(f'{m[1]}{m[2]}', int(clean['J'].split('/')[0]))] = float(clean['Level (cm-1)'])
    return lev


def parse_ff(path):
    """Referee table -> {(lo_nl, 2J, up_nl, 2J): dict(gf_L, sigma_th)} so gf_L * dE_exp / sigma_th is the referee gf on
    experimental energies. NIST MCHF .lsj.db (FF) or a transcribed A/f CSV (Safronova: gf = g_lo f, already on NIST energies)."""
    if path is None:
        return {}
    if str(path).endswith('.csv'):
        out = {}
        for r in csv.DictReader(open(path)):
            lo, up = (re.fullmatch(r'(\d+[spdf])(\d)/2', x).groups() for x in (r['lower'], r['upper']))
            out[(lo[0], int(lo[1]), up[0], int(up[1]))] = dict(gf_L=(int(lo[1]) + 1) * float(r['f']), sigma_th=1e8 / float(r['wl_vac_A']),
                                                               unc_pct=float(r['unc_pct']))
        return out
    return V61.parse_ff(path)


def nmin_by_l(levels):
    """Lowest valence n per l from the NIST levels (Na: 3,3,3,4; Ca II: 4,4,3,4)."""
    out = {}
    for nl, _ in levels:
        l = 'spdf'.index(nl[-1])
        out[l] = min(out.get(l, 99), int(nl[:-1]))
    return out


def grid():
    x = np.linspace(X0, X1, NPTS)
    return x, np.exp(x), x[1] - x[0]


def potential(r, l, p):
    a1, a2, a3, a4, rc = p
    zl = ZC + (Z - ZC) * np.exp(-a1 * r) - r * (a3 + a4 * r) * np.exp(-a2 * r)
    return -zl / r - ALPHA_D / (2 * r ** 4) * (1 - np.exp(-(r / rc) ** 6))


def solve(l, p, nstates):
    """Bound states of -y'' + (2 r^2 V + (l+1/2)^2) y = 2E r^2 y on x = ln r (u = sqrt(r) y), generalized
    sparse problem solved by shift-invert (accuracy relative to the shift; the plain tridiagonal form is
    ill-conditioned at small r). Returns energies (hartree) and z = r y with sum z^2 dx = int u^2 dr = 1,
    ordered by node count, so index k has k radial nodes (n = l + 1 + k)."""
    import scipy.sparse as sp
    from scipy.sparse.linalg import eigsh
    x, r, dx = grid()
    diag = 2 / dx ** 2 + 2 * r ** 2 * potential(r, l, p) + (l + 0.5) ** 2
    off = -np.ones(len(x) - 1) / dx ** 2
    A = sp.diags([off, diag, off], [-1, 0, 1], format='csc')
    B = sp.diags(r ** 2, format='csc')
    w, v = eigsh(A, k=nstates + 6, M=B, sigma=SP['sigma'], which='LM')
    keep = {}
    for lam, y in zip(w, v.T):
        if lam >= 0:
            continue
        z = r * y
        z = z / math.sqrt((z ** 2).sum() * dx)
        big = np.abs(z) > 1e-6 * np.abs(z).max()
        zz = z[big]
        nodes = int((np.sign(zz[1:]) != np.sign(zz[:-1])).sum())
        if z[np.argmax(np.abs(z))] < 0:
            z = -z
        keep[nodes] = (lam / 2, z)
    ks = sorted(keep)
    if ks[:nstates] != list(range(ks[0], ks[0] + min(nstates, len(ks)))):
        raise RuntimeError(f'l={l}: node ladder has gaps {ks}')
    return np.array([keep[k][0] for k in ks]), np.array([keep[k][1] for k in ks]).T, ks[0]


def nist_terms(levels):
    """J-averaged term energies (hartree, below the limit) keyed (n, l)."""
    acc = {}
    for (nl, twoj), e in levels.items():
        k = (int(nl[:-1]), 'spdf'.index(nl[-1]))
        acc.setdefault(k, []).append((twoj + 1, e))
    return {k: -(IP_CM - sum(g * e for g, e in v) / sum(g for g, _ in v)) / RY_CM for k, v in acc.items()}


def fit(levels, nfit=9):
    terms = nist_terms(levels)
    params, report = {}, {}
    for l in range(4):
        ns = [n for (n, ll) in terms if ll == l and n <= nfit]
        nmin = NMIN[l]
        target = np.array([terms[(n, l)] for n in sorted(ns)])
        idx = [n - nmin for n in sorted(ns)]

        def resid(q, rc=START[l][4]):
            w, _, k0 = solve(l, (*q, rc), max(idx) + 1)
            off = (nmin - l - 1) - k0              # valence state nmin has nmin-l-1 nodes
            return (w[[i + off for i in idx]] - target) / np.abs(target)
        r0 = resid(np.array(START[l][:4]))
        sol = least_squares(resid, np.array(START[l][:4]), bounds=([0.05, 0.05, -20, -20], [50, 50, 20, 20]), max_nfev=300)
        params[l] = (*sol.x, START[l][4])
        report['spdf'[l]] = dict(n_terms=len(ns), max_rel_err_start=float(np.max(np.abs(r0))), max_rel_err_fit=float(np.max(np.abs(sol.fun))),
                                 params=[round(float(v), 8) for v in params[l]])
    return params, report


def states(params):
    """(n, l) -> z vector. Lowest valence n: 3 for s,p,d; 4 for f. States below it are core-like and skipped."""
    out = {}
    for l in range(4):
        nmin = NMIN[l]
        w, v, k0 = solve(l, params[l], NMAX - nmin + 1)
        off = (nmin - l - 1) - k0
        for i in range(off, v.shape[1]):
            out[(nmin + i - off, l)] = (w[i], v[:, i])
    return out


def dipole(r, rc):
    return r * (1 - ALPHA_D / r ** 3 * (1 - np.exp(-(r / rc) ** 3)))


def radial(z1, z2, rc):
    _, r, dx = grid()
    t = z1 * z2 * dipole(r, rc) * dx
    s = t.sum()
    return s, abs(s) / np.abs(t).sum()


@lru_cache(None)
def _f(n):
    return math.factorial(n)


def sixj(j1, j2, j3, j4, j5, j6):
    """Racah formula, arguments as Fractions (half-integers allowed)."""
    def tri(a, b, c):
        if (a + b - c) < 0 or (a - b + c) < 0 or (-a + b + c) < 0 or (a + b + c).denominator != 1:
            return None
        return Fr(_f(int(a + b - c)) * _f(int(a - b + c)) * _f(int(-a + b + c)), _f(int(a + b + c + 1)))
    ts = [tri(j1, j2, j3), tri(j1, j5, j6), tri(j4, j2, j6), tri(j4, j5, j3)]
    if any(t is None for t in ts):
        return 0.0
    pre = math.sqrt(float(ts[0] * ts[1] * ts[2] * ts[3]))
    a = [j1 + j2 + j3, j1 + j5 + j6, j4 + j2 + j6, j4 + j5 + j3]
    b = [j1 + j2 + j4 + j5, j2 + j3 + j5 + j6, j3 + j1 + j6 + j4]
    s = 0
    for t in range(int(max(a)), int(min(b)) + 1):
        den = _f(t - int(a[0])) * _f(t - int(a[1])) * _f(t - int(a[2])) * _f(t - int(a[3])) * _f(int(b[0]) - t) * _f(int(b[1]) - t) * _f(int(b[2]) - t)
        s += (-1) ** t * _f(t + 1) / den
    return pre * s


def rc_op(l1, l2, params):
    return 0.5 * (params[min(l1, 2)][4] + params[min(l2, 2)][4])


def lines(params, levels):
    st_ = states(params)
    out = {}
    for (n1, l1), (_, z1) in st_.items():
        for (n2, l2), (_, z2) in st_.items():
            if l2 != l1 + 1:
                continue
            rc = rc_op(l1, l2, params)
            R, cf = radial(z1, z2, rc)
            r_lo, _ = radial(z1, z2, rc * 0.8)
            r_hi, _ = radial(z1, z2, rc * 1.2)
            _, rr, dx = grid()
            r_noop = (z1 * z2 * rr * dx).sum()          # bare r operator (no core-polarization correction)
            for a, b in (((n1, l1), (n2, l2)), ((n2, l2), (n1, l1))):
                for j1 in (2 * a[1] - 1, 2 * a[1] + 1):
                    for j2 in (2 * b[1] - 1, 2 * b[1] + 1):
                        if j1 < 0 or j2 < 0 or abs(j1 - j2) > 2:
                            continue
                        ka, kb = (f'{a[0]}{"spdf"[a[1]]}', j1), (f'{b[0]}{"spdf"[b[1]]}', j2)
                        if ka not in levels or kb not in levels or levels[kb] <= levels[ka]:
                            continue
                        frac = (j1 + 1) * (j2 + 1) * sixj(Fr(a[1]), Fr(j1, 2), Fr(1, 2), Fr(j2, 2), Fr(b[1]), Fr(1)) ** 2
                        de = (levels[kb] - levels[ka]) / RY_CM
                        lmax = max(a[1], b[1])

                        def gf(rr):
                            return 2 / 3 * de * lmax * rr * rr * frac
                        g = gf(R)
                        if g <= 0:
                            continue
                        out[(ka[0], j1, kb[0], j2)] = dict(gf=g, cf=cf, sens=abs(math.log10(gf(r_hi) / gf(r_lo))) / 2 if gf(r_lo) > 0 and gf(r_hi) > 0 else 9.0,
                                                          gf_noop=gf(r_noop))
    return out


def main():
    levels = nist_levels()
    global NMIN
    NMIN = nmin_by_l(levels)
    ff9, ff16 = parse_ff(SP['ff9']), parse_ff(SP['ff16'])
    params, fitrep = fit(levels)
    L = lines(params, levels)
    rows = []
    for k, g in L.items():
        lam = 1e8 / (levels[(k[2], k[3])] - levels[(k[0], k[1])])
        row = dict(lower=f'{k[0]}_{k[1]}/2', upper=f'{k[2]}_{k[3]}/2', wl_vac_A=round(lam, 3), loggf=round(math.log10(g['gf']), 4),
                   loggf_no_polarization_operator=round(math.log10(g['gf_noop']), 4) if g['gf_noop'] > 0 else None,
                   cf=round(g['cf'], 4), sens_dex=round(g['sens'], 4), cls=None)
        for tag, ff in (('ff9', ff9), ('ff16', ff16)):
            if k in ff:
                e = levels[(k[2], k[3])] - levels[(k[0], k[1])]
                row[f'{tag}_loggf'] = round(math.log10(ff[k]['gf_L'] * e / ff[k]['sigma_th']), 4)
                row[f'd_{tag}'] = round(row['loggf'] - row[f'{tag}_loggf'], 4)
        rows.append(row)
    rows.sort(key=lambda r: r['wl_vac_A'])
    global FLOOR
    FLOOR = floor(L, levels)
    sig = calibrate(rows)
    for r in rows:
        r['sigma_dex'], r['sigma_basis'] = line_sigma(r, sig)
        r['cls'] = sigma_class(r['sigma_dex'], r['cf'])
    # experiment: resonance lines
    exp = {}
    for k, (a, s) in SP['exp'].items():
        g = L[k]['gf']
        exp[f'{k[0]}_{k[1]}/2-{k[2]}_{k[3]}/2'] = dict(A_exp=a, frac_sigma=s, A_model=V61.rate(g, k, levels), dlogA=round(math.log10(V61.rate(g, k, levels) / a), 4),
                                       dlogA_without_polarization_operator=round(math.log10(V61.rate(L[k]['gf_noop'], k, levels) / a), 4))

    def stats(sel):
        d = [r['d_ff9'] for r in sel if 'd_ff9' in r]
        return dict(n=len(d), median=round(st.median(d), 4), robust_sigma=V61.robust(d),
                    frac_within_0p05=round(sum(abs(x) <= 0.05 for x in d) / len(d), 2)) if d else dict(n=0)
    by_cls = {c: stats([r for r in rows if r['cls'] == c]) for c in 'ABCDE'}
    by_cf = {b: stats([r for r in rows if lo <= r['cf'] < hi]) for b, lo, hi in (('cf>=0.5', 0.5, 9), ('0.1-0.5', 0.1, 0.5), ('0.02-0.1', 0.02, 0.1), ('<0.02', 0, 0.02))}
    OUT.mkdir(parents=True, exist_ok=True)
    with open(OUT / 'lines.csv', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=sorted({k for r in rows for k in r}), lineterminator='\n')
        w.writeheader()
        w.writerows(rows)
    # fill the solar line list, per observable FEATURE (fine-structure components closer than 0.15 A are one feature)
    comp = {}
    for r in rows:
        comp.setdefault((r['lower'].split('_')[0], r['upper'].split('_')[0]), []).append(r)
    fill = []
    for f in solar_features(levels):
        cs = [r for r in comp.get((f['lower_nl'], f['upper_nl']), []) if abs(1e8 / r['wl_vac_A'] - 1e8 / f['vac_A']) <= 0.6]
        g = sum(10 ** r['loggf'] for r in cs)
        g9 = [10 ** r['ff9_loggf'] for r in cs if 'ff9_loggf' in r]
        cf = min((r['cf'] for r in cs), default=None)
        nf = max((nmax(r) for r in cs), default=None)
        fill.append(dict(feature=f"{f['lower_nl']}-{f['upper_nl']}", components='|'.join(f"{r['lower']}-{r['upper']}" for r in cs),
                         air_A=f['air_A'], n_linelist_rows=f['n_rows'], max_depth=f['depth'], linelist_loggf=round(math.log10(f['gf']), 3),
                         linelist_sources=f['sources'], model_loggf=round(math.log10(g), 4) if g else None, model_cf=cf,
                         model_sigma_dex=round(math.sqrt(sum((10 ** r['loggf'] * r['sigma_dex']) ** 2 for r in cs)) / g, 4) if cs else None,
                         model_class=min((r['cls'] for r in cs), default=None, key=lambda c: -'ABCDE'.index(c)),
                         sigma_basis=max((r['sigma_basis'] for r in cs), default=None, key=['referee', 'upper_bound', 'no_referee'].index),
                         ff9_loggf=round(math.log10(sum(g9)), 4) if g9 and len(g9) == len(cs) else None,
                         in_nist_ff9=bool(g9) and len(g9) == len(cs),
                         d_model_minus_ff9=round(math.log10(g / sum(g9)), 4) if g9 and len(g9) == len(cs) else None,
                         d_model_minus_linelist=round(math.log10(g / f['gf']), 3) if g else None))
    with open(OUT / 'solar_fill.csv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(fill[0]), lineterminator='\n')
        w.writeheader()
        w.writerows(fill)
    meas = [x for x in fill if x['max_depth'] >= 0.01]
    summary = dict(schema='rya1354.onevalence-model.v1', species=SP['name'], status=STATUS, alpha_d=ALPHA_D, floor_dex=FLOOR, potential_fit=fitrep, experiment=exp,
                   vs_ff9_all=stats(rows), vs_ff9_by_class=by_cls, vs_ff9_by_cf=by_cf, n_lines=len(rows),
                   referee=str(SP['ff9'].relative_to(ROOT)), nmin=NMIN,
                   class_counts={c: sum(r['cls'] == c for r in rows) for c in 'ABCDE'},
                   sigma68_by_cf=sig, sigma_coverage_vs_ff9=coverage(rows),
                   fill=dict(n_lines=len(fill), n_filled=sum(x['model_loggf'] is not None for x in fill), n_measurable=len(meas),
                             n_measurable_filled=sum(x['model_loggf'] is not None for x in meas),
                             n_measurable_not_in_ff9=sum(x['model_loggf'] is not None and not x['in_nist_ff9'] for x in meas),
                             asplund_lines_vs_ff9={x['air_A']: x['d_model_minus_ff9'] for x in fill
                                                   if round(x['air_A'], 0) in SP['asplund']}))
    (OUT / 'summary.json').write_text(json.dumps(summary, indent=1) + '\n')
    return summary


CF_CANCEL = 0.1   # below this the radial integral cancels strongly (e.g. the Cooper-minimum 3s-4p)
CF_BINS = ((0.5, 9.0), (0.2, 0.5), (0.1, 0.2), (0.05, 0.1), (0.0, 0.05))


def nmax(r):
    return max(int(re.match(r'\d+', r['lower'])[0]), int(re.match(r'\d+', r['upper'])[0]))


def calibrate(rows):
    """Candidate per-line sigma: 68th percentile of |model - FF B-spline CI| in bins of our radial CF.
    Regime 'n<=7': FF converged (FF [9] and [16] agree to 0.007 dex) -> a real referee.
    Regime 'n8-9': FF begins to degrade -> sigma68 there is an UPPER bound on ours.
    n >= 10: FF [9] states are box-confined (uniform -0.47 dex offset of every n=10 transition) -> no referee;
    such lines take the n8-9 bound and are flagged. Bins with < 8 lines take the next lower-CF bin's value."""
    out = {}
    for reg, (n_lo, n_hi) in (('n<=7', (3, 7)), ('n8-9', (8, 9))):
        b = []
        for lo, hi in CF_BINS:
            d = sorted(abs(r['d_ff9']) for r in rows if 'd_ff9' in r and lo <= r['cf'] < hi and n_lo <= nmax(r) <= n_hi)
            b.append(dict(cf_lo=lo, cf_hi=hi, n=len(d), sigma68_vs_ff9=round(d[int(0.68 * (len(d) - 1))], 4) if len(d) >= 8 else None))
        for k in range(len(b) - 1, -1, -1):
            if b[k]['sigma68_vs_ff9'] is None:
                b[k]['sigma68_vs_ff9'] = b[min(k + 1, len(b) - 1)]['sigma68_vs_ff9']
        out[reg] = b
    # sparse referees (Ca II: 55 lines) leave whole regimes empty: fall back to the other regime, then to all lines
    alld = sorted(abs(r['d_ff9']) for r in rows if 'd_ff9' in r)
    overall = round(alld[int(0.68 * (len(alld) - 1))], 4) if alld else 0.1
    for reg, other in (('n<=7', 'n8-9'), ('n8-9', 'n<=7')):
        for k, b in enumerate(out[reg]):
            if b['sigma68_vs_ff9'] is None:
                b['sigma68_vs_ff9'] = out[other][k]['sigma68_vs_ff9'] or overall
                b['fallback'] = True
    # transition family (l_lower-l_upper): 68th percentile of |model - FF| on the FF-converged lines (n <= 7);
    # includes any family bias (we do NOT shift our model toward FF, which is theory, not lab)
    fam = {}
    for r in rows:
        if 'd_ff9' in r and nmax(r) <= 7 and r['cf'] >= CF_CANCEL:
            fam.setdefault(family(r), []).append(abs(r['d_ff9']))
    out['family'] = {k: dict(n=len(v), abs68=round(sorted(v)[int(0.68 * (len(v) - 1))], 4)) for k, v in sorted(fam.items()) if len(v) >= 5}
    return out


def family(r):
    return re.sub(r'\d', '', r['lower'].split('_')[0]) + '-' + re.sub(r'\d', '', r['upper'].split('_')[0])


FLOOR = None


def floor(L, levels):
    """Experimental floor: the larger of the model's mean |offset| on the resonance-line anchors and the anchors' own
    stated uncertainty, rounded UP to 0.001 dex (Na I: 0.0054 -> 0.006)."""
    off = st.mean(abs(math.log10(V61.rate(L[k]['gf'], k, levels) / a)) for k, (a, _) in SP['exp'].items())
    unc = max(math.log10(1 + s) for _, s in SP['exp'].values())
    return math.ceil(max(off, unc) * 1000 - 1e-9) / 1000


def line_sigma(r, sig):
    """Candidate per-line sigma (dex):
      n <= 7 : hypot(FLOOR, half the core-polarization-operator effect, the line's transition-family sigma vs FF)
               (family term added for Mg II, whose d-series carry 0.03-0.06 dex the operator indicator misses)
      n 8-9  : also hypot'd with the n8-9 sigma68 of model-FF in the line's CF bin (upper bound, FF degrading)
      n >= 10: same as n 8-9, flagged no_referee."""
    half = abs(r['loggf'] - r['loggf_no_polarization_operator']) / 2 if r.get('loggf_no_polarization_operator') is not None else 0.5
    fam = sig['family'].get(family(r), {}).get('abs68', 0.0)
    s = math.hypot(FLOOR, half, fam)
    if r['cf'] < CF_CANCEL:   # strong cancellation: the family term is calibrated without these, add their CF-bin term
        b = next(x for x in sig['n<=7'] if x['cf_lo'] <= r['cf'] < x['cf_hi'])
        s = math.hypot(s, b['sigma68_vs_ff9'])
    n = nmax(r)
    if n >= 8:
        b = next(x for x in sig['n8-9'] if x['cf_lo'] <= r['cf'] < x['cf_hi'])
        s = math.hypot(s, b['sigma68_vs_ff9'])
    return round(s, 4), ('referee' if n <= 7 else 'upper_bound' if n <= 9 else 'no_referee')


def sigma_class(sigma_dex, cf):
    import rya1354_theory_grade as G
    return G.accuracy_class(10 ** sigma_dex - 1, cf)


def coverage(rows):
    z = sorted(abs(r['d_ff9']) / r['sigma_dex'] for r in rows if 'd_ff9' in r and r['sigma_basis'] == 'referee')
    return dict(n=len(z), frac_within_1sigma=round(sum(x <= 1 for x in z) / len(z), 3), frac_within_2sigma=round(sum(x <= 2 for x in z) / len(z), 3),
                median_z=round(z[len(z) // 2], 3))


def solar_features(levels):
    """Na I rows of the solar line list (HFS components) -> (lower nl, upper nl) assignment by NIST energies ->
    features: rows of the same nl pair within 0.15 A are one observable feature."""
    rows = []
    for r in csv.DictReader(open(ROOT / 'data/linelists/linelist_solar.csv')):
        if (r['element'], r['ion']) != SP['linelist']:
            continue
        air = float(r['wavelength_air_A'])
        s2 = (1e4 / air) ** 2
        vac = air * (1 + 8.336624212083e-5 + 2.408926869968e-2 / (1.301065924522e2 - s2) + 1.599740894897e-4 / (3.892568793293e1 - s2))
        elo = float(r['excitation_potential_eV']) * V61.HC_EV_CM
        best = None
        for lk, el in levels.items():
            if abs(el - elo) > 2.0:
                continue
            for uk, eu in levels.items():
                dl = 'spdf'.index(uk[0][-1]) - 'spdf'.index(lk[0][-1])
                c = abs(eu - (el + 1e8 / vac)) + abs(el - elo)
                if dl in (-1, 1) and abs(uk[1] - lk[1]) <= 2 and c < 1.5 and (best is None or c < best[0]):
                    best = (c, lk[0], uk[0])
        if best:
            rows.append(dict(lower_nl=best[1], upper_nl=best[2], air=air, vac=vac, gf=10 ** float(r['log_gf']),
                             depth=float(r['central_depth'] or 0), src=r['loggf_source']))
    rows.sort(key=lambda x: (x['lower_nl'], x['upper_nl'], x['vac']))
    feats = []
    for x in rows:
        f = feats[-1] if feats else None
        if f and (f['lower_nl'], f['upper_nl']) == (x['lower_nl'], x['upper_nl']) and x['vac'] - f['_last'] <= 0.15:
            f['_m'].append(x)
            f['_last'] = x['vac']
        else:
            feats.append(dict(lower_nl=x['lower_nl'], upper_nl=x['upper_nl'], _m=[x], _last=x['vac']))
    out = []
    for f in feats:
        m = f['_m']
        w = sum(x['gf'] for x in m)
        out.append(dict(lower_nl=f['lower_nl'], upper_nl=f['upper_nl'], n_rows=len(m), gf=w,
                        air_A=round(sum(x['air'] * x['gf'] for x in m) / w, 3), vac_A=sum(x['vac'] * x['gf'] for x in m) / w,
                        depth=max(x['depth'] for x in m), sources='|'.join(sorted({x['src'] for x in m}))))
    return sorted(out, key=lambda f: f['air_A'])


if __name__ == '__main__':
    import sys
    for name in (sys.argv[1:] or list(SPECIES)):
        use(name)
        print(name, json.dumps({k: v for k, v in main().items() if k in ('experiment', 'vs_ff9_by_class', 'sigma_coverage_vs_ff9', 'fill', 'floor_dex')}, indent=1))
