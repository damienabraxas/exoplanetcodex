"""RYA-1354 v70b: photoionization cross-sections from OUR Mg I two-valence CI model (declared on RYA-1354 before
any code or run, 2026-11-21).

Method
  * the v70 CI basis (unchanged orbitals) is enlarged by box-continuum orbitals: positive-energy eigenstates of the
    screened outer-electron potential inside a wall at R_box, projected against the core-like states and the v70
    span, re-diagonalized in the ion Hamiltonian h_l inside that complement. The one-body matrix is then no longer
    diagonal (orig-continuum coupling) -> twovalence.hamiltonian takes B.h1 (regression: identical to v70 without it);
  * configurations: every v70 pair + {Mg+ channel orbital (3s, 4s, 3p, 4p, 3d) x continuum orbital};
  * box normalisation: df/de_n = f_n * rho_n, rho_n = 2 / (e_{n+1} - e_{n-1}) over the eigenstates above the Mg+ 3s
    threshold of one final symmetry (one open channel below Mg+ 3p for every symmetry with par == L % 2);
    sigma = 2 pi^2 alpha a0^2 df/de (4.0336e-18 cm^2 per unit df/de in hartree^-1);
  * several R_box sample the discrete spectrum at shifted energies -> resonances are resolved by pooling;
  * length gauge with the core-polarization dipole operator (as v70) and velocity gauge (bare nabla).
Tests (declared): C0 bound model intact, C1 box independence, C2 gauge; T1 3p2 1S0 resonance (Er, Gamma, peak),
T2/T3 3s3p 1P absolute sigma (Rafiq 2007 354.5 nm, Madsen & Thomsen 2002 285 nm), T4 static polarizability sum rule
(Mitroy 2010), T5 bound-free continuity at threshold. Report: OP (Butler 1993) tables carried by the OpenMULTI3D atom.
Status: CANDIDATE / HOLD. Writes only data/audit/rya1354/mg1_photoion_v70b.

Usage:
  rya1354_photoion.py box R       one box radius (a0): spectra + f for every bound initial term -> box_R.npz/.json
  rya1354_photoion.py c0 R        C0 control: v70 LS gf in the enlarged basis at R (and 'c0 0' = no continuum)
  rya1354_photoion.py alpha       T4: static polarizabilities, v70 basis and enlarged basis
  rya1354_photoion.py analyze     pool boxes -> sigma(e), C1, C2, T1-T5, OP comparison, summary.json
"""
import csv
import json
import math
import re
import sys
import time
from pathlib import Path

import numpy as np
from scipy.linalg import eigh
from scipy.optimize import least_squares

sys.path.insert(0, str(Path(__file__).resolve().parent))
import rya1354_twovalence_model as T  # noqa: E402
import rya1354_onevalence_model as M  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'data/audit/rya1354/mg1_photoion_v70b'
STATUS = 'CANDIDATE_HOLD_NOT_ADOPTED'
NPTS = 6400               # finer log grid: dr <= lambda/12 at r = 430 a0 for k <= 0.45 (e <= 0.10 hartree)
EMAX = 0.10               # hartree above Mg+ 3s (21,950 cm-1; below Mg+ 3p at 35,700)
NCH = {0: 2, 1: 2, 2: 1}  # channel (ion) orbitals per l: Mg+ 3s 4s / 3p 4p / 3d
BOXES = (360, 370, 380, 390, 400, 410, 420, 430)
LMAX_SYM = 4
MB_PER_DFDE = 2 * math.pi ** 2 / 137.035999084 * 0.529177210903e-8 ** 2 * 1e18   # Mb per unit df/de (hartree^-1)
OM_ATOM = Path.home() / 'codex/openmulti3d_atoms/atom.mg86g_nq_qmax'
LAB = dict(   # Baig 2022 (Atoms 10, 39) compilation; primaries NOT HELD (IOP/APS)
    T1=dict(Er=68270.0, Er_err=5.0, Gamma=280.0, Gamma_err=5.0, q=24.5, peak_Mb=785.0, peak_err=141.0,
            src='Rafiq et al. 2007 (Bonanno 1986, Bradley 1976 agree) via Baig 2022 Table 5'),
    T2=dict(lam_nm=354.5, sigma_Mb=90.0, err=16.0, src='Rafiq et al. 2007 via Baig 2022'),
    T3=dict(lam_nm=285.0, sigma_Mb=81.0, err=23.0, src='Madsen & Thomsen 2002 (MOT) via Baig 2022'),
    T4={'3s2 1S': (71.33, 'Mitroy+2010 Table XI, RCI+MBPT [57]; CICP 71.35; sum-rule expt 74.9(2.7)'),
        '3s.3p 3P*': (101.2, 'Mitroy+2010 Table XI, 3s3p 3P0 hybrid RCI+MBPT 101.2(3); CICP 101.5')},
)
INIT_1P = '3s.3p 1P*'


def setup():
    T.GRID['NPTS'] = NPTS
    sp, B = T.build('Mg I')
    return sp, B


def add_continuum(B, rbox):
    """Append box-continuum orbitals (0 < e <= 1.15 EMAX) for every l; sets B.n0, B.h1, B.keep."""
    js = B.screen()
    dx = B.dx
    B.n0, B.h1, B.ncont = {}, {}, {}
    for l in range(T.LMAX + 1):
        eps0, Z0 = B.orb[l]
        n0 = len(eps0)
        k = int(rbox * math.sqrt(2 * 1.15 * EMAX) / math.pi) + 30
        for _ in range(4):
            e, z = B.eig(l, k, 1.15 * EMAX, extra=js, wall=rbox)      # sigma in A units (= 2E): centre 0.575 EMAX
            if e.min() < 0 and e.max() > 1.15 * EMAX:
                break
            k = int(k * 1.6)
        else:
            raise RuntimeError(f'l={l}: continuum window not covered ({e.min():.4f}..{e.max():.4f})')
        cand = z[:, (e > 0) & (e <= 1.15 * EMAX)]
        core = B.core[l]
        for _ in range(2):
            cand = cand - core @ (core.T @ cand * dx)
            cand = cand - Z0 @ (Z0.T @ cand * dx)
        S = cand.T @ cand * dx
        s, U = np.linalg.eigh(S)
        keep = s > 1e-6 * s.max()
        X = cand @ (U[:, keep] / np.sqrt(s[keep]))
        A0 = B.A(l)
        Yc, Y0 = X / B.r[:, None], Z0 / B.r[:, None]
        Hc = 0.5 * (Yc.T @ (A0 @ Yc)) * dx
        ec, C = np.linalg.eigh(0.5 * (Hc + Hc.T))
        Zc = X @ C
        Yc = Zc / B.r[:, None]
        h0c = 0.5 * (Y0.T @ (A0 @ Yc)) * dx
        h00 = 0.5 * (Y0.T @ (A0 @ Y0)) * dx
        dev = np.abs(h00 - np.diag(eps0)).max()
        if dev > 1e-6:
            raise RuntimeError(f'l={l}: v70 orbitals not diagonal in h_l on this grid ({dev:.2e})')
        h = np.zeros((n0 + len(ec),) * 2)
        h[:n0, :n0] = np.diag(eps0)
        h[n0:, n0:] = np.diag(ec)
        h[:n0, n0:] = h0c
        h[n0:, :n0] = h0c.T
        B.orb[l] = (np.concatenate([eps0, ec]), np.hstack([Z0, Zc]))
        B.n0[l], B.h1[l], B.ncont[l] = n0, h, len(ec)
    n0, nch = B.n0, NCH

    def keep(c):
        la, ia, lb, ib = c
        ca, cb = ia >= n0[la], ib >= n0[lb]
        if not ca and not cb:
            return True
        if ca and cb:
            return False
        return (ib < nch.get(lb, 0)) if ca else (ia < nch.get(la, 0))
    B.keep = keep
    return B


def syms():
    return [(L, S, p) for L in range(LMAX_SYM + 1) for S in (0, 1) for p in (0, 1)]


def solve_all(sp, B, emax_rel=EMAX * 1.05, full=False, below_ip_only=False):
    """{sym: dict(cf, e, v, terms, nist_e, match)} with every root below threshold + emax_rel (all roots if full)."""
    I = T.Integrals(B)
    terms = T.nist_terms(sp)
    e_thr = B.orb[0][0][0]
    out = {}
    for sym in syms():
        L, S, par = sym
        cf = T.configs(B, L, S, par)
        if not cf:
            continue
        t0 = time.time()
        H = T.hamiltonian(B, I, cf, L, S)
        e, v = (eigh(H) if full else eigh(H, subset_by_value=(-np.inf, e_thr + emax_rel)))
        lst = [t for t in terms.get(sym, []) if T.nmax2(t[1].split(' ')[0]) <= T.BOX_N]
        if below_ip_only:     # v70b C0 attempt 2: autoionizing NIST terms (3p2 1S, +6604 cm-1) are T1's, not bound-model
            lst = [t for t in lst if t[0] < sp['ip']]
        nist_e = np.array([-(M.IP_CM) / M.RY_CM - (sp['ip'] - t[0]) / M.RY_CM for t in lst])
        bound = int((e < e_thr).sum())
        match = T.align(e[:bound], nist_e) if len(lst) else {}
        out[sym] = dict(cf=cf, e=e, v=v, terms=lst, nist_e=nist_e, match=match)
        print(f'  sym {sym}: {len(cf)} cf, {len(e)} roots ({bound} bound), {time.time() - t0:.0f}s', flush=True)
        I.R.clear()
    return out, e_thr


def f_and_fv(B, d1, sym1, d2, sym2, i1):
    """Absorption f (length+polop, velocity) from root(s) i1 of sym1 to every root of sym2 (term-averaged LS).
    i1 may be an int or an index array (then f has shape (len(i1), n_roots2)); D is built once per sym pair."""
    (L1, S, _), (L2, _, _) = sym1, sym2
    D = T.dipole_matrix(B, d1['cf'], L1, d2['cf'], L2, S, True)
    DV = T.dipole_matrix(B, d1['cf'], L1, d2['cf'], L2, S, gauge='velocity')
    ci = d1['v'][:, i1]
    rl, rv = ci.T @ D @ d2['v'], ci.T @ DV @ d2['v']
    de = d2['e'][None, :] - np.atleast_1d(d1['e'][i1])[:, None]
    if np.ndim(i1) == 0:
        de = de[0]
    with np.errstate(divide='ignore', invalid='ignore'):
        f = 2 / 3 * de * rl ** 2 / (2 * L1 + 1)
        fv = 2 / 3 * rv ** 2 / de / (2 * L1 + 1)
    return f, fv, de


def run_box(rbox):
    OUT.mkdir(parents=True, exist_ok=True)
    sp, B = setup()
    if rbox:
        add_continuum(B, rbox)
    print(f'R_box {rbox}: continuum per l {getattr(B, "ncont", {})}', flush=True)
    st, e_thr = solve_all(sp, B)
    rec = dict(rbox=rbox, e_thr=e_thr, npts=NPTS, emax=EMAX, ncont=getattr(B, 'ncont', {}), inits={})
    arrays = {}
    for sym1, d1 in st.items():
        js = sorted(d1['match'])
        if not js:
            continue
        idx = np.array([d1['match'][j] for j in js])
        labels = [d1['terms'][j][1] for j in js]
        for j, lab, i1 in zip(js, labels, idx):
            rec['inits'][lab] = dict(sym=sym1, e_model=float(d1['e'][i1]), e_nist=float(d1['nist_e'][j]), finals={})
        for sym2, d2 in st.items():
            (L1, S1, p1), (L2, S2, p2) = sym1, sym2
            if S1 != S2 or p1 == p2 or abs(L1 - L2) > 1 or L1 + L2 == 0:
                continue
            t0 = time.time()
            f, fv, _ = f_and_fv(B, d1, sym1, d2, sym2, idx)
            for k, lab in enumerate(labels):
                key = f'{lab}|{sym2}'
                arrays[key + '|e'] = d2['e'] - e_thr
                arrays[key + '|f'] = f[k]
                arrays[key + '|fv'] = fv[k]
                rec['inits'][lab]['finals'][str(sym2)] = key
            print(f'  dipoles {sym1}->{sym2}: {len(labels)} initial terms, {time.time() - t0:.0f}s', flush=True)
    np.savez_compressed(OUT / f'box_{rbox}.npz', **arrays)
    (OUT / f'box_{rbox}.json').write_text(json.dumps(rec, indent=1, default=str) + '\n')


def run_c0(rbox):
    """C0: the v70 LS gf of the Pehlivan lab lines in the enlarged basis (R_box) or on the fine grid alone (0)."""
    OUT.mkdir(parents=True, exist_ok=True)
    sp, B = setup()
    if rbox:
        add_continuum(B, rbox)
    st, e_thr = solve_all(sp, B, emax_rel=-1e-6, below_ip_only=True)   # bound roots, bound NIST terms (attempt 2)
    mult = T.multiplets(B, st)
    rows = {(r['lower'], r['upper']): r for r in T.fs_lines(st, mult)}
    v70 = {(r['lower'], r['upper']): r for r in csv.DictReader(open(ROOT / 'data/audit/rya1354/mg1_model_v70/lines.csv'))}
    out = []
    for k, r in v70.items():
        if r['lab_loggf'] in ('', None) or k not in rows or r['loggf_pure_LS'] in ('', None):
            continue
        out.append(dict(lower=k[0], upper=k[1], wl_vac_A=r['wl_vac_A'], loggf_v70_LS=float(r['loggf_pure_LS']),
                        loggf_here_LS=rows[k]['loggf'], d=round(rows[k]['loggf'] - float(r['loggf_pure_LS']), 4)))
    erel = max(abs(d['e'][i] / d['nist_e'][j] - 1) for d in st.values() for j, i in d['match'].items())
    v70e = json.load(open(ROOT / 'data/audit/rya1354/mg1_model_v70/summary.json'))['energy_gate']
    res = dict(rbox=rbox, n_lab_lines=len(out), max_abs_d=max(abs(x['d']) for x in out), energy_max_rel=round(erel, 5),
               v70_energy_gate=v70e, lines=out)
    (OUT / f'c0_{rbox}.json').write_text(json.dumps(res, indent=1, default=str) + '\n')
    print(f'C0 R={rbox}: n={len(out)} max|d|={res["max_abs_d"]:.4f} energy max rel {erel:.5f}', flush=True)


def run_alpha():
    """T4: alpha0 = sum_f f/dE^2 (complete pseudo-spectrum of the basis, model energies) + alpha_core(Mg2+)."""
    OUT.mkdir(parents=True, exist_ok=True)
    res = {}
    for tag, rbox in (('v70_basis', 0), ('enlarged_R400', 400)):
        sp, B = setup()
        if rbox:
            add_continuum(B, rbox)
        st, e_thr = solve_all(sp, B, full=True)
        r = {}
        for label, (ref, src) in LAB['T4'].items():
            sym1 = next(s for s, d in st.items() for j in d['match'] if d['terms'][j][1] == label)
            d1 = st[sym1]
            i1 = next(i for j, i in d1['match'].items() if d1['terms'][j][1] == label)
            tot, parts, totv, totn = 0.0, {}, 0.0, 0.0
            e1n = d1['nist_e'][next(j for j, i in d1['match'].items() if i == i1)]
            for sym2, d2 in st.items():
                (L1, S1, p1), (L2, S2, p2) = sym1, sym2
                if S1 != S2 or p1 == p2 or abs(L1 - L2) > 1 or L1 + L2 == 0:
                    continue
                f, fv, de = f_and_fv(B, d1, sym1, d2, sym2, i1)
                ok = np.abs(de) > 1e-9
                a = float(np.sum(f[ok] / de[ok] ** 2))
                parts[str(sym2)] = round(a, 4)
                tot += a
                totv += float(np.sum(fv[ok] / de[ok] ** 2))
                # post-hoc variant (2026-11-21, after the declared model-energy T4 read +4.4 %): NIST energies for the
                # initial and every NIST-matched final term, model energies elsewhere (Mitroy's 'hybrid' practice);
                # f rescaled with the energy (f = 2/3 dE |D|^2 / (2L+1)) so alpha ~ |D|^2 / dE
                den = de.copy()
                for j2, i2 in d2['match'].items():
                    den[i2] = d2['nist_e'][j2] - e1n
                okn = ok & (np.abs(den) > 1e-9)
                totn += float(np.sum((f[okn] / de[okn]) / den[okn]))
            r[label] = dict(alpha_valence=round(tot, 3), alpha_core=M.ALPHA_D, alpha=round(tot + M.ALPHA_D, 3),
                            alpha_velocity=round(totv + M.ALPHA_D, 3), by_final_sym=parts, ref=ref, ref_src=src,
                            alpha_nist_energies_POSTHOC=round(totn + M.ALPHA_D, 3),
                            rel_nist_energies_POSTHOC=round((totn + M.ALPHA_D) / ref - 1, 4),
                            rel=round((tot + M.ALPHA_D) / ref - 1, 4), pass_2pct=abs((tot + M.ALPHA_D) / ref - 1) <= 0.02)
            print(f'T4 {tag} {label}: alpha {tot + M.ALPHA_D:.2f} vs {ref} ({r[label]["rel"]:+.2%})', flush=True)
        res[tag] = r
    (OUT / 'alpha_T4.json').write_text(json.dumps(res, indent=1) + '\n')


if __name__ == '__main__':
    mode = sys.argv[1]
    if mode == 'box':
        run_box(int(sys.argv[2]))
    elif mode == 'c0':
        run_c0(int(sys.argv[2]))
    elif mode == 'alpha':
        run_alpha()
    elif mode == 'analyze':
        import rya1354_photoion_analyze as A
        A.main()
    else:
        raise SystemExit(__doc__)
