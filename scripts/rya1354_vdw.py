"""RYA-1354 v70 B: van der Waals (H-collision) broadening from the model's own wavefunctions (declared 2026-11-20).

1. ABO standard: n* from NIST energies vs the ion limit; Anstee & O'Mara 1995 / Barklem & O'Mara 1997 / Barklem,
   O'Mara & Ross 1998 tables (Turbospectrum DATA/Anstee-1802-2014.dat, bilinear as in anstee.f); sigma, alpha ->
   log gamma6 (FWHM, rad s-1 cm3 per H atom) at T = 10^4 K.
2. n*_r: the active electron's <r^2> = <Psi|r1^2+r2^2|Psi> - <3s|r^2|3s>(ion), inverted through the hydrogenic
   <r^2> = n*^2 (5 n*^2 + 1 - 3 l(l+1)) / 2; same tables. Flag |n*_r - n*| > 0.1.
3. Unsold C6 with hydrogenic vs model <r^2> (report only).
Report only; nothing enters a synthesis in v70."""
import json, math, re, sys
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import rya1354_twovalence_model as T
import rya1354_onevalence_model as M

ROOT = T.ROOT
TABLES = Path.home() / 'codex/engines/Turbospectrum_NLTE/DATA/Anstee-1802-2014.dat'
A0 = 5.29177210903e-11
KB = 1.380649e-23
AMU = 1.66053906660e-27
R_MG = 109734.83  # Rydberg for 24Mg (cm-1)


def load_tables(path=TABLES):
    shapes = [(21, 18), (21, 18), (18, 18), (18, 18), (18, 18), (18, 18)]
    rows = []
    for line in open(path):
        if line.startswith('*') or not line.strip():
            continue
        try:
            vals = [float(x) for x in line.split()]
        except ValueError:
            break
        if len(vals) != 18:
            break
        rows.append(vals)
        if len(rows) == sum(s[0] for s in shapes):
            break
    out, k = [], 0
    for nr, nc in shapes:
        out.append(np.array(rows[k:k + nr]))
        k += nr
    # (row axis start, col axis start) for each pair; step 0.1
    return {('s', 'p'): (out[0], out[1], 1.0, 1.3), ('p', 'd'): (out[2], out[3], 1.3, 2.3), ('d', 'f'): (out[4], out[5], 2.3, 3.3)}


def bilin(tab, x0, y0, x, y):
    nr, nc = tab.shape
    i, j = (x - x0) / 0.1, (y - y0) / 0.1
    if not (0 <= i <= nr - 1 and 0 <= j <= nc - 1):
        return None
    i0, j0 = min(int(i), nr - 2), min(int(j), nc - 2)
    fi, fj = i - i0, j - j0
    return ((1 - fi) * (1 - fj) * tab[i0, j0] + fi * (1 - fj) * tab[i0 + 1, j0] + (1 - fi) * fj * tab[i0, j0 + 1]
            + fi * fj * tab[i0 + 1, j0 + 1])


def abo(T_, l1, n1, l2, n2):
    """sigma (a0^2), alpha for a transition between active-electron l1 (n*1) and l2 (n*2)."""
    a, b = sorted([(l1, n1), (l2, n2)], key=lambda t: 'spdfgh'.index(t[0]))
    key = (a[0], b[0])
    if key not in T_:
        return None, None
    s, al, x0, y0 = T_[key]
    sg = bilin(s, x0, y0, a[1], b[1])
    if sg is None:
        return None, None
    return sg, bilin(al, x0, y0, a[1], b[1])


def gamma6_abo(sigma, alpha, mass_amu, T=1e4):
    mu = 1 / (1 / 1.008 + 1 / mass_amu) * AMU
    vbar = math.sqrt(8 * KB * T / (math.pi * mu))
    v0 = 1e4
    w = (4 / math.pi) ** (alpha / 2) * math.gamma((4 - alpha) / 2) * v0 * sigma * A0 ** 2 * (vbar / v0) ** (1 - alpha)
    return math.log10(2 * w * 1e6)  # FWHM per perturber, cm3 s-1


def gamma6_unsold(dr2_a0, mass_amu, T=1e4):
    """Unsold: gamma6 = 8.08 vbar^(3/5) C6^(2/5) (cgs), C6 = e^2 alpha_H |d<r^2>| / hbar."""
    if dr2_a0 <= 0:
        return None
    e2 = (4.80320e-10) ** 2
    hbar = 1.054571817e-27
    aH = 6.6669e-25
    c6 = e2 * aH * dr2_a0 * (A0 * 100) ** 2 / hbar
    mu = 1 / (1 / 1.008 + 1 / mass_amu) * AMU * 1e3
    vbar = math.sqrt(8 * KB * 1e7 * T / (math.pi * mu))
    return math.log10(8.08 * vbar ** 0.6 * c6 ** 0.4)


def r2_hyd(n, l):
    return n * n * (5 * n * n + 1 - 3 * l * (l + 1)) / 2


def nstar_from_r2(r2, l):
    from scipy.optimize import brentq
    f = lambda n: r2_hyd(n, l) - r2
    lo = l + 0.5 if l > 0 else 0.5
    try:
        return brentq(f, max(lo, 0.3), 60)
    except ValueError:
        return float('nan')


def r2_expect(B, cf, L, S, C):
    """<r1^2 + r2^2> for CI vectors C (columns) in the antisymmetrized config basis cf."""
    f = {l: B.orb[l][1].T @ (B.orb[l][1] * (B.r ** 2)[:, None]) * B.dx for l in B.orb}
    n = [1 / math.sqrt(2 * (1 + ((-1) ** (S + la + lb + L)) * (la == lb and ia == ib))) for la, ia, lb, ib in cf]
    idx = {c: i for i, c in enumerate(cf)}

    def prim(la, ia, lb, ib, lc, ic, ld, id_):
        if la != lc or lb != ld:
            return 0.0
        v = 0.0
        if ib == id_:
            v += f[la][ia, ic]
        if ia == ic:
            v += f[lb][ib, id_]
        return v
    F = np.zeros((len(cf), len(cf)))
    by = {}
    for j, (lc, ic, ld, id_) in enumerate(cf):
        by.setdefault((lc, ic), []).append(j)
        by.setdefault((ld, id_), []).append(j)
    for i, (la, ia, lb, ib) in enumerate(cf):
        for j in set(by.get((la, ia), []) + by.get((lb, ib), [])):
            lc, ic, ld, id_ = cf[j]
            eta = (-1) ** (S + lc + ld + L)
            F[i, j] = 2 * n[i] * n[j] * (prim(la, ia, lb, ib, lc, ic, ld, id_) + eta * prim(la, ia, lb, ib, ld, id_, lc, ic))
    return np.einsum('ik,ij,jk->k', C, F, C)


def active(label, corepfx='3s'):
    conf = label.split(' ')[0]
    toks = conf.split('.')
    if conf == f'{corepfx}2':
        return 's', False
    if len(toks) == 2 and toks[0] == corepfx:
        return re.fullmatch(r'\d+([spdfgh])', toks[1]).group(1), False
    m = re.fullmatch(r'\d+([spdfgh])2', conf)
    return (m.group(1) if m else re.fullmatch(r'\d+([spdfgh])', toks[-1]).group(1)), True


def main(name='Mg I', out=None, mass=24.305):
    sp, B = T.build(name)
    _, S_ = T.run_ci(sp, B)
    S_ = dict(S_)
    S_.pop('_beyond_box')
    z3s = B.ion_ground
    r2_core = float((z3s ** 2 * B.r ** 2).sum() * B.dx)
    ip = sp['ip']
    lev = []
    for (L, S, par), d in S_.items():
        js = sorted(d['match'])
        if not js:
            continue
        C = d['v'][:, [d['match'][j] for j in js]]
        r2 = r2_expect(B, d['cf'], L, S, C)
        for k, j in enumerate(js):
            E, lab, _ = d['terms'][j]
            la, doubly = active(lab, sp['corepfx'] if 'corepfx' in sp and sp['corepfx'] != '2p6' else '3s')
            l = 'spdfgh'.index(la)
            nst = math.sqrt(R_MG / (ip - E)) if E < ip else float("nan")   # above the limit: autoionizing, outside ABO
            ra = float(r2[k]) - r2_core
            lev.append(dict(term=lab, E=E, l=la, doubly_excited=doubly, nstar=nst, r2_total=float(r2[k]), r2_active=ra,
                            r2_hyd=r2_hyd(nst, l), nstar_r=nstar_from_r2(ra, l)))
    lv = pd.DataFrame(lev)
    lv['dn'] = lv.nstar_r - lv.nstar
    lv['flag'] = (lv.dn.abs() > 0.1) | lv.doubly_excited
    out = Path(out or ROOT / 'data/audit/rya1354' / sp['out'])
    lv.to_csv(out / 'vdw_levels.csv', index=False)
    TT = load_tables()
    Ls = pd.read_csv(out / 'lines.csv')
    tm = {r.term.replace('*', ''): r for r in lv.itertuples()}
    rows = []
    for r in Ls.itertuples():
        a, b = tm.get(re.sub(r'\d+$', '', r.lower).replace('*', '')), tm.get(re.sub(r'\d+$', '', r.upper).replace('*', ''))
        if a is None or b is None:
            continue
        s0, a0 = abo(TT, a.l, a.nstar, b.l, b.nstar)
        s1, a1 = abo(TT, a.l, a.nstar_r, b.l, b.nstar_r)
        g0 = gamma6_abo(s0, a0, mass) if s0 else None
        g1 = gamma6_abo(s1, a1, mass) if s1 else None
        uh = gamma6_unsold(abs(b.r2_hyd - a.r2_hyd), mass)
        um = gamma6_unsold(abs(b.r2_active - a.r2_active), mass)
        rows.append(dict(lower=r.lower, upper=r.upper, wl_vac_A=r.wl_vac_A, loggf=r.loggf, type=f'{a.l}-{b.l}',
                         nstar_lo=round(a.nstar, 3), nstar_up=round(b.nstar, 3), nstar_r_lo=round(a.nstar_r, 3), nstar_r_up=round(b.nstar_r, 3),
                         sigma_abo=s0, alpha_abo=a0, log_g6_abo=g0, sigma_abo_r=s1, alpha_abo_r=a1, log_g6_abo_r=g1,
                         log_g6_unsold_hyd=uh, log_g6_unsold_model=um, flag=bool(a.flag or b.flag)))
    V = pd.DataFrame(rows)
    V.to_csv(out / 'vdw_lines.csv', index=False)
    print(lv.round(3).to_string())
    return lv, V


def v2a(w):
    s2 = (1e4 / w) ** 2
    return w / (1 + 0.0000834254 + 0.02406147 / (130 - s2) + 0.00015998 / (38.9 - s2))


def referee(out, element=('Mg', 'I')):
    """VALD log gamma6 (linelist_solar.csv) vs ours; fill candidates where VALD has none. Report only."""
    out = Path(out)
    V = pd.read_csv(out / 'vdw_lines.csv')
    d = pd.read_csv(ROOT / 'data/linelists/linelist_solar.csv', low_memory=False)
    v = d[(d.element == element[0]) & (d.ion == element[1])]
    V['air'] = np.where(V.wl_vac_A > 2000, v2a(V.wl_vac_A), V.wl_vac_A)
    m = pd.merge_asof(V.sort_values('air'), v.sort_values('wavelength_air_A')[['wavelength_air_A', 'log_gf', 'damping_vdW']],
                      left_on='air', right_on='wavelength_air_A', direction='nearest', tolerance=0.08)
    m = m[m.wavelength_air_A.notna() & ((m.loggf - m.log_gf).abs() < 1)]
    m.to_csv(out / 'vdw_vs_vald.csv', index=False)
    x = m[m.log_g6_abo.notna() & (m.damping_vdW < 0)]
    dd = x.log_g6_abo - x.damping_vdW
    r = V[V.log_g6_abo.notna() & V.log_g6_abo_r.notna()]
    dr = (r.log_g6_abo_r - r.log_g6_abo).abs()
    u = V[V.log_g6_unsold_hyd.notna() & V.log_g6_unsold_model.notna()]
    du = (u.log_g6_unsold_model - u.log_g6_unsold_hyd)
    fill = m[(m.damping_vdW == 0) & m.log_g6_abo.notna()]
    fill[['lower', 'upper', 'air', 'loggf', 'type', 'log_g6_abo', 'log_g6_abo_r', 'sigma_abo', 'alpha_abo']].to_csv(out / 'vdw_fill_candidates.csv', index=False)
    doubly = u[u.flag & (du.abs() > 0.3)]
    s = dict(n_lines=len(V), n_in_abo_tables=int(V.log_g6_abo.notna().sum()),
             control_vs_vald=dict(n=len(x), median=round(float(dd.median()), 4), abs68=round(float((dd - dd.median()).abs().quantile(.68)), 4)),
             nstar_r_effect=dict(n=len(r), p50=round(float(dr.median()), 4), p90=round(float(dr.quantile(.9)), 4), max=round(float(dr.max()), 4)),
             unsold_model_vs_hyd=dict(n=len(u), p50=round(float(du.abs().median()), 4), p90=round(float(du.abs().quantile(.9)), 4)),
             n_matched_vald=len(m), n_vald_zero=int((m.damping_vdW == 0).sum()), n_fill_candidates_in_table=len(fill),
             n_doubly_excited_unsold_shift_gt_0p3=len(doubly),
             note='report only (v70); nothing enters a synthesis')
    (out / 'vdw_summary.json').write_text(json.dumps(s, indent=1) + '\n')
    print(json.dumps(s, indent=1))
    return s


if __name__ == '__main__':
    if '--referee' in sys.argv:
        referee(sys.argv[-1])
    else:
        lv, V = main(sys.argv[1] if len(sys.argv) > 1 else 'Mg I')
        referee(ROOT / 'data/audit/rya1354' / T.SPECIES[sys.argv[1] if len(sys.argv) > 1 else 'Mg I']['out'])
