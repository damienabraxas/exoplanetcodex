"""RYA-1354 v70 A: level lifetimes and radiative damping from the model's own E1 A-values (declared 2026-11-20).

tau(level) = 1 / sum_k A_k over every model E1 decay in lines.csv (spin-changing included).
A_k = 0.6670e16 * gf / (g_u * lambda_vac[A]^2).  sigma_tau: branch-weighted quadrature of each line's sigma_dex
(lines without sigma, CF < 0.02, take the 0.05 floor); a fully correlated (linear) version is reported too.
Lab referee: Aldenius et al. 2007 (A&A 461, 767) LIF lifetimes; 3D fine structure unresolved -> the term's
g-weighted mean decay rate. PASS: |z| <= 2 for >= 5 of 6 terms. Gamma_rad(line) = Gamma_u + Gamma_l (s-1)."""
import csv, json, math, re, sys
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
LN10 = math.log(10)
FLOOR = 0.05
# Aldenius et al. 2007 Table 1 (ns, uncertainty); levels as labelled in lines.csv
ALDENIUS = {'3s.4s 3S': (11.5, 1.0), '3s.5s 3S': (29.0, 3.0), '3s.3d 3D': (5.9, 0.4), '3s.4d 3D': (17.6, 1.2),
            '3s.5d 3D': (33.0, 3.0), '3p2 3P': (2.0, 0.2)}


def gJ(label):
    return 2 * int(re.search(r'(\d+)$', label).group(1)) + 1


def term_of(label):
    return re.sub(r'\d+$', '', label).replace('*', '')


def main(model_dir):
    d = Path(model_dir)
    L = pd.read_csv(d / 'lines.csv')
    L['g_u'] = L.upper.map(gJ)
    L['A'] = 0.6670e16 * 10 ** L.loggf / (L.g_u * L.wl_vac_A ** 2)
    L['sig'] = L.sigma_dex.fillna(FLOOR).clip(lower=FLOOR)
    lev = {}
    for u, g in L.groupby('upper'):
        A = g.A.sum()
        w = g.A / A
        lev[u] = dict(A_tot=A, tau_ns=1e9 / A, rel_sig_q=math.sqrt(((w * LN10 * g.sig) ** 2).sum()),
                      rel_sig_lin=float((w * LN10 * g.sig).sum()), n_branches=len(g),
                      top_branch=g.loc[g.A.idxmax(), 'lower'], top_bf=float(w.max()))
    lv = pd.DataFrame(lev).T
    lv.index.name = 'level'
    lv['term'] = [term_of(x) for x in lv.index]
    lv['g'] = [gJ(x) for x in lv.index]
    lv.to_csv(d / 'levels_lifetimes.csv')
    # Gamma_rad per line (lower level: 0 if it is the ground state or has no model decays)
    Ag = lv.A_tot.to_dict()
    L['log_gamma_rad'] = np.log10(L.upper.map(Ag) + L.lower.map(Ag).fillna(0.0))
    L[['lower', 'upper', 'wl_vac_A', 'loggf', 'log_gamma_rad']].to_csv(d / 'lines_gamma_rad.csv', index=False)
    rows = []
    for t, (tl, el) in ALDENIUS.items():
        s = lv[lv.term == t.replace('*', '')]
        if s.empty:
            rows.append(dict(term=t, lab_ns=tl, lab_err=el, status='NOT IN MODEL'))
            continue
        A = float((s.g * s.A_tot).sum() / s.g.sum())
        rq = float((s.g * s.rel_sig_q).sum() / s.g.sum())
        rl = float((s.g * s.rel_sig_lin).sum() / s.g.sum())
        tm = 1e9 / A
        z = (tm - tl) / math.hypot(el, tm * rq)
        zl = (tm - tl) / math.hypot(el, tm * rl)
        rows.append(dict(term=t, levels=';'.join(s.index), lab_ns=tl, lab_err=el, model_ns=round(tm, 3),
                         model_err_q=round(tm * rq, 3), model_err_lin=round(tm * rl, 3), d_pct=round(100 * (tm / tl - 1), 1),
                         z=round(z, 2), z_corr=round(zl, 2), top_branch=s.top_branch.iloc[0], top_bf=round(float(s.top_bf.iloc[0]), 3)))
    R = pd.DataFrame(rows)
    R.to_csv(d / 'lifetimes_vs_aldenius2007.csv', index=False)
    npass = int((R.z.abs() <= 2).sum())
    summ = dict(rule='PASS if |z|<=2 for >=5 of 6 Aldenius 2007 terms (declared RYA-1354 v70)', n_pass=npass,
                verdict='PASS' if npass >= 5 else 'FAIL', terms=rows,
                median_abs_d_pct=float(R.d_pct.abs().median()), n_levels=len(lv))
    (d / 'lifetimes_summary.json').write_text(json.dumps(summ, indent=1, default=str) + '\n')
    print(R.to_string())
    print(summ['verdict'], npass)
    return lv, L



def vald_referee(model_dir, element=('Mg', 'I')):
    """Gamma_rad vs VALD/Kurucz damping_rad in linelist_solar.csv (report only)."""
    from rya1354_vdw import v2a
    d = Path(model_dir)
    g = pd.read_csv(d / 'lines_gamma_rad.csv')
    v = pd.read_csv(ROOT / 'data/linelists/linelist_solar.csv', low_memory=False)
    v = v[(v.element == element[0]) & (v.ion == element[1])]
    g['air'] = np.where(g.wl_vac_A > 2000, v2a(g.wl_vac_A), g.wl_vac_A)
    m = pd.merge_asof(g.sort_values('air'), v.sort_values('wavelength_air_A')[['wavelength_air_A', 'log_gf', 'damping_rad']],
                      left_on='air', right_on='wavelength_air_A', direction='nearest', tolerance=0.08)
    m = m[m.wavelength_air_A.notna() & ((m.loggf - m.log_gf).abs() < 1) & (m.damping_rad != 0)]
    m['d'] = m.log_gamma_rad - m.damping_rad
    m.to_csv(d / 'gamma_rad_vs_vald.csv', index=False)
    s = dict(n=len(m), median=round(float(m.d.median()), 4), abs68=round(float((m.d - m.d.median()).abs().quantile(.68)), 4),
             p90_abs=round(float(m.d.abs().quantile(.9)), 4))
    j = json.loads((d / 'lifetimes_summary.json').read_text())
    j['gamma_rad_vs_vald'] = s
    (d / 'lifetimes_summary.json').write_text(json.dumps(j, indent=1, default=str) + '\n')
    print(s)


if __name__ == '__main__':
    main(sys.argv[1])
    vald_referee(sys.argv[1])
