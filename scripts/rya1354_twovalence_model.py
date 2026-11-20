"""RYA-1354 v65/v66: our own atomic model for TWO-VALENCE-ELECTRON spectra (Mg I, Ca I), declared on RYA-1354
before the run.

CI with core polarization (Norcross & Seaton 1976; Mitroy & Bromley 2003 style), built on the one-valence model:
  * core: the l-dependent Mg2+/Ca2+ model potential V_l already fitted to Mg II/Ca II (v63/v64), which carries the
    one-electron core polarization -alpha_d/(2 r^4) W(r);
  * two-electron Hamiltonian H = h(1) + h(2) + 1/r12 + V_d, with the dielectronic polarization
    V_d = -alpha_d (r1 . r2)/(r1^3 r2^3) sqrt(W(r1) W(r2)),  W(r) = 1 - exp(-(r/rc)^6);
  * one-electron basis per l: ion orbitals (inner electron), bound states of a screened potential V_l + J_ns
    (neutral-atom Rydberg electron) and compact small-box functions (correlation), orthonormalized and
    re-diagonalized in h_l;
  * full CI in LS coupling for each (L, S, parity); no spin-orbit (spin-changing lines flagged, not graded);
  * dipole operator with the core-polarization correction d(r) = r [1 - alpha_d/r^3 (1 - exp(-(r/rc)^3))];
  * gf per fine-structure line = (2/3) dE_NIST (2J+1)(2J'+1) {L J S; J' L' 1}^2 |<L||D||L'>|^2.
Status: CANDIDATE / HOLD (RYA-1352). Nothing here writes into a product.

Usage: rya1354_twovalence_model.py [SPECIES ...]   (default: all)
"""
import csv
import json
import math
import re
import statistics as st
import sys
from fractions import Fraction as Fr
from functools import lru_cache
from pathlib import Path

import warnings
import numpy as np
warnings.filterwarnings("ignore", message=".*encountered in matmul")
import scipy.sparse as sps
from scipy.linalg import eigh
from scipy.sparse.linalg import eigsh

import rya1354_onevalence_model as M
import rya1354_na_model as V61

ROOT = Path(__file__).resolve().parents[1]
REF = ROOT / 'data/reference'
STATUS = 'CANDIDATE_HOLD_NOT_ADOPTED'
LMAX = 4
SPECIES = {
    'Mg I': dict(core='Mg II', ip=61671.05, levels=REF / 'onevalence_nist/MgI_levels.csv', corepfx='2p6',
                 lab=REF / 'mg_gf_lab/mg1_lab_loggf_pehlivan2017.csv', linelist=('Mg', 'I'), out='mg1_model_v69', labprefix='3s',
                 intercomb_ref=(('3s2 1S0', '3s.3p 3P*1', -5.686, 'LAB: Kwong et al. 1982 lifetime (cited Bergemann 2017)'),
                                ('3s2 1S0', '3s.3p 3P*1', -5.623, 'theory: KP08 Table 40 / Tachiev & Froese Fischer BP-MCHF, class D')),
                 forbidden={2154.353: 'E2', 2089.149: 'M2', 4562.602: 'M2', 2084.537: 'E2'}),
    'Ca I': dict(core='Ca II', ip=49305.92, levels=REF / 'onevalence_nist/CaI_levels.csv', corepfx='3p6',
                 lab=REF / 'lab_gf_tables/transcribed/denhartog2021_ca1.csv', linelist=('Ca', 'I'), out='ca1_model_v66', labprefix='4s'),
    # v69 Part B: Mg I isoelectronic sister ions, BLIND (same code and rules, nothing retuned). Referee = NIST ASD graded
    # spin-allowed E1 (mostly MCHF compilations, not lab-independent); blind lab targets = the intercombination lines.
    'Al II': dict(core='Al III', ip=151862.5, levels=REF / 'onevalence_nist/AlII_levels.csv', corepfx='2p6',
                  lab=REF / 'sister_ions_nist/al2_nist_graded.csv', linelist=('Al', 'II'), out='al2_model_v69', labprefix='3s',
                  percentages=REF / 'onevalence_nist/AlII_levels.csv', floor_fallback=0.050,
                  intercomb_ref=(('3s2 1S0', '3s.3p 3P*1', -4.979, 'LAB: Traebert et al. 1999 storage-ring lifetime (NIST ASD A+, T6947)'),)),
    'Si III': dict(core='Si IV', ip=270139.3, levels=REF / 'onevalence_nist/SiIII_levels.csv', corepfx='2p6',
                   lab=REF / 'sister_ions_nist/si3_nist_graded.csv', linelist=('Si', 'III'), out='si3_model_v69', labprefix='3s',
                   percentages=REF / 'onevalence_nist/SiIII_levels.csv', floor_fallback=0.050,
                   intercomb_ref=(('3s2 1S0', '3s.3p 3P*1', -4.570, 'LAB: Kwong et al. 1983 RF ion trap, A = 1.67(10)e4 s-1 '
                                   '(as cited by A&A 2015 aa24394-14; primary paper not yet read)'),)),
}
# CI grid: same log mesh as the one-valence solver but coarser and to 900 a0 (neutral Rydberg n <= ~12)
GRID = dict(X0=math.log(1e-5), X1=math.log(900.0), NPTS=3200)
NION = 6         # ion orbitals per l (inner electron)
NRYD = 14        # screened-potential bound states per l (outer electron)
NCOMP = 14       # compact functions per l (correlation), box radius RCOMP
RCOMP = 25.0
ALPHA2_SCALE = 1.0    # v68 sensitivity: alpha_d in the TWO-electron terms only (dielectronic + dipole correction);
RC_DIEL_SCALE = 1.0   # the fitted one-electron core potential is left untouched
BOX_N = 12       # neutral Rydberg states with n above this are box-perturbed on this grid -> no_referee


# ---------------------------------------------------------------- angular algebra
@lru_cache(None)
def threej0(a, b, c):
    """(a b c; 0 0 0) for integers."""
    J = a + b + c
    if J % 2 or c < abs(a - b) or c > a + b:
        return 0.0
    g = J // 2
    f = math.factorial
    return (-1) ** g * math.sqrt(f(J - 2 * a) * f(J - 2 * b) * f(J - 2 * c) / f(J + 1)) * f(g) / (f(g - a) * f(g - b) * f(g - c))


@lru_cache(None)
def ck(l1, k, l2):
    """<l1||C^k||l2>."""
    return (-1) ** l1 * math.sqrt((2 * l1 + 1) * (2 * l2 + 1)) * threej0(l1, k, l2)


@lru_cache(None)
def sixj(a, b, c, d, e, f):
    return M.sixj(Fr(a), Fr(b), Fr(c), Fr(d), Fr(e), Fr(f))


@lru_cache(None)
def ang_coulomb(la, lb, lc, ld, L):
    """{k: A_k} with <(la lb)L| C^k(1).C^k(2) |(lc ld)L> = A_k (electron 1: a->c, electron 2: b->d)."""
    out = {}
    for k in range(max(abs(la - lc), abs(lb - ld)), min(la + lc, lb + ld) + 1):
        a = ck(la, k, lc) * ck(lb, k, ld)
        if a:
            s = sixj(la, lb, L, ld, lc, k)
            if s:
                out[k] = (-1) ** (lb + lc + L) * s * a
    return out


# ---------------------------------------------------------------- one-electron basis
class Basis:
    def __init__(self, params, nmin, alpha, rc_of_l):
        M.X0, M.X1, M.NPTS = GRID['X0'], GRID['X1'], GRID['NPTS']
        self.x, self.r, self.dx = M.grid()
        self.params, self.nmin, self.alpha, self.rc = params, nmin, alpha * ALPHA2_SCALE, rc_of_l
        sig0 = M.SP['sigma']
        M.SP['sigma'] = -4.0 * M.Z ** 2
        try:
            _, z0, _ = M.solve(0, params[0], nmin[0])
        finally:
            M.SP['sigma'] = sig0
        self.ion_ground = z0[:, nmin[0] - 1]             # ion ground orbital (Mg+ 3s, Ca+ 4s): screens the outer electron
        self.orb = {}     # l -> (eps[n], z[:, n])
        self.label = {}
        for l in range(LMAX + 1):
            self.orb[l] = self.build(l)

    def A(self, l, extra=None, wall=None):
        r, dx = self.r, self.dx
        V = M.potential(r, l, self.params[min(l, 3)])
        if extra is not None:
            V = V + extra
        if wall is not None:
            V = np.where(r > wall, 1e4, V)
        diag = 2 / dx ** 2 + 2 * r ** 2 * V + (l + 0.5) ** 2
        off = -np.ones(len(r) - 1) / dx ** 2
        return sps.diags([off, diag, off], [-1, 0, 1], format='csc')

    def eig(self, l, k, sigma, extra=None, wall=None):
        A = self.A(l, extra, wall)
        B = sps.diags(self.r ** 2, format='csc')
        w, v = eigsh(A, k=k, M=B, sigma=sigma, which='LM')
        o = np.argsort(w)
        z = (self.r[:, None] * v[:, o])
        z /= np.sqrt((z ** 2).sum(0) * self.dx)
        return w[o] / 2, z

    def build(self, l):
        """Ion orbitals + screened Rydberg states + compact functions; core-like states projected out; the span is
        re-diagonalized in the ion Hamiltonian h_l so the CI basis is orthonormal with diagonal h."""
        r, dx = self.r, self.dx
        nm = self.nmin.get(l, l + 1)
        ncore = nm - l - 1
        sig0 = M.SP['sigma']
        M.SP['sigma'] = -4.0 * M.Z ** 2                  # reach the deepest core-like state; states ordered by nodes
        try:
            _, z, k0 = M.solve(l, self.params[min(l, 3)], ncore + NION)
        finally:
            M.SP['sigma'] = sig0
        if k0 != 0:
            raise RuntimeError(f'l={l}: lowest state has {k0} nodes; core-like states missing')
        core = z[:, :ncore]
        ion = z[:, ncore:ncore + NION]
        # screening by the ion ground orbital (direct Hartree, asymptotically +1/r)
        zs = self.ion_ground
        rho = zs ** 2 * dx
        Js = np.cumsum(rho) / r + (np.cumsum(rho[::-1] / r[::-1])[::-1] - rho / r)
        wr, zr = self.eig(l, NRYD + ncore + 6, -0.5, extra=Js)
        zr = zr[:, wr < 0][:, :NRYD + ncore + 2]
        wc, zc = self.eig(l, NCOMP + ncore + 2, 2.0, wall=RCOMP)
        cand = np.hstack([ion, zr, zc])
        cand -= core @ (core.T @ cand * dx)               # Pauli: project out the core-like states
        # orthonormalize (canonical, drop near-linear-dependence)
        S = cand.T @ cand * dx
        s, U = np.linalg.eigh(S)
        keep = s > 1e-7 * s.max()
        X = cand @ (U[:, keep] / np.sqrt(s[keep]))
        Y = X / r[:, None]
        Ay = self.A(l) @ Y
        H = 0.5 * (Y.T @ Ay) * dx
        H = 0.5 * (H + H.T)
        e, C = np.linalg.eigh(H)
        Z = X @ C
        for i in range(Z.shape[1]):
            if Z[np.argmax(np.abs(Z[:, i])), i] < 0:
                Z[:, i] = -Z[:, i]
        return e, Z

    def gfun(self, la, lc):
        rc = 0.5 * (self.rc[min(la, 2)] + self.rc[min(lc, 2)]) * RC_DIEL_SCALE
        return np.sqrt(1 - np.exp(-(self.r / rc) ** 6)) / self.r ** 2

    def dip(self, la, lc, polop=True):
        r = self.r
        if not polop:
            return r
        rc = 0.5 * (self.rc[min(la, 2)] + self.rc[min(lc, 2)])
        return r * (1 - self.alpha / r ** 3 * (1 - np.exp(-(r / rc) ** 3)))


class Integrals:
    """Slater integrals R^k(ac, bd) = int int P_a P_c (r<^k / r>^{k+1}) P_b P_d by blocks of l, plus the dielectronic
    term folded into k = 1."""
    def __init__(self, B):
        self.B = B
        self.Y = {}
        self.R = {}

    def W(self, la, lc):
        za, zc = self.B.orb[la][1], self.B.orb[lc][1]
        return (za[:, :, None] * zc[:, None, :]).reshape(len(self.B.r), -1).T * self.B.dx   # (Na*Nc, npts)

    def Yk(self, la, lc, k):
        key = (la, lc, k)
        if key not in self.Y:
            r = self.B.r
            w = self.W(la, lc)
            inner = np.cumsum(w * r ** k, axis=1) / r ** (k + 1)
            outer = (np.cumsum((w / r ** (k + 1))[:, ::-1], axis=1)[:, ::-1] - w / r ** (k + 1)) * r ** k
            self.Y[key] = inner + outer
        return self.Y[key]

    def block(self, la, lc, lb, ld, k):
        """R^k[a, c, b, d] array."""
        key = (la, lc, lb, ld, k)
        if key not in self.R:
            R = self.Yk(la, lc, k) @ self.W(lb, ld).T
            if k == 1 and self.B.alpha:
                ga = self.W(la, lc) @ self.B.gfun(la, lc)
                gb = self.W(lb, ld) @ self.B.gfun(lb, ld)
                R = R - self.B.alpha * np.outer(ga, gb)
            na, nc, nb, nd = (len(self.B.orb[x][0]) for x in (la, lc, lb, ld))
            self.R[key] = R.reshape(na, nc, nb, nd)
        return self.R[key]


# ---------------------------------------------------------------- CI
def configs(B, L, S, par):
    """List of (la, ia, lb, ib) with (la, ia) <= (lb, ib), parity par, triangle, antisymmetry allowed."""
    out = []
    for la in range(LMAX + 1):
        for lb in range(la, LMAX + 1):
            if (la + lb) % 2 != par or not abs(la - lb) <= L <= la + lb:
                continue
            na, nb = len(B.orb[la][0]), len(B.orb[lb][0])
            eta = (-1) ** (S + la + lb + L)
            for ia in range(na):
                for ib in range(ia if la == lb else 0, nb):
                    if la == lb and ia == ib and eta == -1:
                        continue
                    out.append((la, ia, lb, ib))
    return out


def hamiltonian(B, I, cf, L, S):
    n = len(cf)
    H = np.zeros((n, n))
    arr = np.array(cf)
    norm = np.array([1 / math.sqrt(2 * (1 + ((-1) ** (S + la + lb + L)) * (la == lb and ia == ib))) for la, ia, lb, ib in cf])
    eta = np.array([(-1) ** (S + la + lb + L) for la, ia, lb, ib in cf])
    classes = sorted({(c[0], c[2]) for c in cf})
    idx = {cl: np.where((arr[:, 0] == cl[0]) & (arr[:, 2] == cl[1]))[0] for cl in classes}
    for (la, lb) in classes:
        i = idx[(la, lb)]
        ia, ib = arr[i, 1], arr[i, 3]
        for (lc, ld) in classes:
            j = idx[(lc, ld)]
            ic, id_ = arr[j, 1], arr[j, 3]
            blk = np.zeros((len(i), len(j)))
            if (la, lb) == (lc, ld):
                e = B.orb[la][0][ia] + B.orb[lb][0][ib]
                blk += np.diag(e) if True else 0
            # exchange one-body: a = d, b = c
            if la == ld and lb == lc:
                m = (ia[:, None] == id_[None, :]) & (ib[:, None] == ic[None, :])
                blk += m * (B.orb[la][0][ia][:, None] + B.orb[lb][0][ib][:, None]) * eta[j][None, :]
            for k, a in ang_coulomb(la, lb, lc, ld, L).items():
                R = I.block(la, lc, lb, ld, k)
                blk += a * R[ia[:, None], ic[None, :], ib[:, None], id_[None, :]]
            for k, a in ang_coulomb(la, lb, ld, lc, L).items():
                R = I.block(la, ld, lb, lc, k)
                blk += a * eta[j][None, :] * R[ia[:, None], id_[None, :], ib[:, None], ic[None, :]]
            H[np.ix_(i, j)] = 2 * norm[i][:, None] * norm[j][None, :] * blk
    return 0.5 * (H + H.T)


def velocity_radial(B, la, lc):
    """int u_a O u_c dr, O = d/dr - (lc+1)/r (la = lc+1) or d/dr + lc/r (la = lc-1); z = u sqrt(r). Identity
    RV = (E_c - E_a) RL verified to 1e-6 with a shared local potential (v68 control)."""
    za, zc = B.orb[la][1], B.orb[lc][1]
    dz = np.gradient(zc, B.x, axis=0)
    du = (dz - zc / 2) / B.r[:, None]
    k = -(lc + 1) if la == lc + 1 else lc
    return za.T @ (du + k * zc / B.r[:, None]) * B.dx


def dipole_matrix(B, cf1, L1, cf2, L2, S, polop=True, gauge='length'):
    """<cf1 L1 S || D || cf2 L2 S> (orbital reduced matrix elements, spin passive). gauge='velocity' gives the
    nabla matrix element (bare operator); divide by the model transition energy for the length-equivalent."""
    rad = {}

    def d(la, lc):
        if (la, lc) not in rad:
            if gauge == 'velocity':
                rad[(la, lc)] = velocity_radial(B, la, lc)
            else:
                za, zc = B.orb[la][1], B.orb[lc][1]
                rad[(la, lc)] = za.T @ (zc * B.dip(la, lc, polop)[:, None]) * B.dx
        return rad[(la, lc)]
    fac = math.sqrt((2 * L1 + 1) * (2 * L2 + 1))

    def prim(la, ia, lb, ib, lc, ic, ld, id_):
        v = 0.0
        if lb == ld and ib == id_ and abs(la - lc) == 1:        # electron 1 jumps
            v += (-1) ** (la + lb + L2 + 1) * fac * sixj(la, L1, lb, L2, lc, 1) * ck(la, 1, lc) * d(la, lc)[ia, ic]
        if la == lc and ia == ic and abs(lb - ld) == 1:        # electron 2 jumps
            v += (-1) ** (la + ld + L1 + 1) * fac * sixj(lb, L1, la, L2, ld, 1) * ck(lb, 1, ld) * d(lb, ld)[ib, id_]
        return v
    D = np.zeros((len(cf1), len(cf2)))
    n1 = [1 / math.sqrt(2 * (1 + ((-1) ** (S + la + lb + L1)) * (la == lb and ia == ib))) for la, ia, lb, ib in cf1]
    n2 = [1 / math.sqrt(2 * (1 + ((-1) ** (S + la + lb + L2)) * (la == lb and ia == ib))) for la, ia, lb, ib in cf2]
    # sparse: index cf2 by each orbital
    by_orb = {}
    for j, (lc, ic, ld, id_) in enumerate(cf2):
        by_orb.setdefault((lc, ic), []).append(j)
        if (ld, id_) != (lc, ic):
            by_orb.setdefault((ld, id_), []).append(j)
    for i, (la, ia, lb, ib) in enumerate(cf1):
        for j in set(by_orb.get((la, ia), []) + by_orb.get((lb, ib), [])):
            lc, ic, ld, id_ = cf2[j]
            eta = (-1) ** (S + lc + ld + L2)
            D[i, j] = 2 * n1[i] * n2[j] * (prim(la, ia, lb, ib, lc, ic, ld, id_) + eta * prim(la, ia, lb, ib, ld, id_, lc, ic))
    return D


def solve_sym(B, I, L, S, par, nroots):
    cf = configs(B, L, S, par)
    H = hamiltonian(B, I, cf, L, S)
    e, v = eigh(H, subset_by_index=[0, min(nroots, len(cf)) - 1])
    return cf, e, v


# ---------------------------------------------------------------- NIST levels
def nist_terms(sp):
    """{(L, S, par): [ (E_term cm-1 J-avg, label, {J: E}) sorted ]} for two-valence configurations above the core."""
    terms = {}
    for r in csv.DictReader(open(sp['levels'])):
        c = {k: (v or '').strip('="[]') for k, v in r.items() if k}
        conf, term, J, E = c['Configuration'], c['Term'], c['J'], c['Level (cm-1)']
        if not E or not re.fullmatch(r'\d+', J or ''):
            continue
        conf = re.sub(r'^' + re.escape(sp['corepfx']) + r'\.', '', conf)
        toks = conf.split('.')
        ne, ok = 0, True
        ls = []
        for t in toks:
            m = re.fullmatch(r'(\d+)([spdfgh])(\d?)', t)
            if not m:
                ok = False
                break
            ne += int(m[3] or 1)
            ls += ['spdfgh'.index(m[2])] * int(m[3] or 1)
        m = re.fullmatch(r'([13])([SPDFGHI])(\*?)', term)
        if not ok or ne != 2 or not m or max(ls) > LMAX:
            continue
        S = (int(m[1]) - 1) // 2
        L = 'SPDFGHI'.index(m[2])
        par = sum(ls) % 2
        if bool(m[3]) != bool(par):
            continue
        terms.setdefault((L, S, par), {}).setdefault((conf, term), {})[int(J)] = float(E.replace('?', '').replace('+x', ''))
    out = {}
    for sym, d in terms.items():
        lst = []
        for (conf, term), js in d.items():
            g = {j: 2 * j + 1 for j in js}
            lst.append((sum(g[j] * e for j, e in js.items()) / sum(g.values()), f'{conf} {term}', js))
        out[sym] = sorted(lst)
    return out


def align(ci_e, nist_e, skip_ci=0.01, skip_nist=0.05):
    """Order-preserving alignment of CI eigenvalues to NIST term energies (both hartree, negative)."""
    n, m = len(ci_e), len(nist_e)
    D = np.full((n + 1, m + 1), np.inf)
    D[0, :] = np.arange(m + 1) * skip_nist
    D[:, 0] = np.arange(n + 1) * skip_ci
    P = {}
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            c = [(D[i - 1, j - 1] + abs(ci_e[i - 1] / nist_e[j - 1] - 1), 'm'), (D[i - 1, j] + skip_ci, 'c'), (D[i, j - 1] + skip_nist, 'n')]
            D[i, j], P[(i, j)] = min(c)
    i, j, pairs = n, m, {}
    while i > 0 and j > 0:
        t = P[(i, j)]
        if t == 'm':
            pairs[j - 1] = i - 1
            i, j = i - 1, j - 1
        elif t == 'c':
            i -= 1
        else:
            j -= 1
    while i > 0 or j > 0:
        break
    return pairs


# ---------------------------------------------------------------- model
def load_core(sp):
    core = M.use(sp['core'])
    lev = M.nist_levels()
    M.NMIN = M.nmin_by_l(lev)
    fit = json.load(open(ROOT / 'data/audit/rya1354' / core['out'] / 'summary.json'))['potential_fit']
    params = {'spdf'.index(k): tuple(v['params']) for k, v in fit.items()}
    nmin = dict(M.NMIN)
    nmin[4] = 5
    return params, nmin


def build(name):
    sp = dict(SPECIES[name], name=name)
    params, nmin = load_core(sp)
    rc = {l: params[l][4] for l in range(3)}
    B = Basis(params, nmin, M.ALPHA_D, rc)
    return sp, B


def run_ci(sp, B, nroots=40):
    I = Integrals(B)
    terms = nist_terms(sp)
    e_ion = -(M.IP_CM) / M.RY_CM                      # ion ground state relative to the bare core (hartree)
    states = {}
    beyond = []
    for sym, lst in sorted(terms.items()):
        L, S, par = sym
        beyond += [t[1] for t in lst if nmax2(t[1].split(' ')[0]) > BOX_N]
        lst = [t for t in lst if nmax2(t[1].split(' ')[0]) <= BOX_N]   # outer n > BOX_N cannot live in the box
        cf, e, v = solve_sym(B, I, L, S, par, max(nroots, len(lst) + 10))
        nist_e = np.array([e_ion - (sp['ip'] - t[0]) / M.RY_CM for t in lst])
        pairs = align(e, nist_e)
        states[sym] = dict(cf=cf, e=e, v=v, terms=lst, nist_e=nist_e, match=pairs)
    states['_beyond_box'] = beyond
    return I, states


# ---------------------------------------------------------------- transitions
def multiplets(B, S_):
    """Reduced orbital dipole <T1||D||T2> for every matched NIST term pair of opposite parity, same S, |dL| <= 1,
    with and without the core-polarization dipole correction, plus the CI cancellation factor
    CF = |sum_ij c1_i D_ij c2_j| / sum_ij |c1_i D_ij c2_j|."""
    out = {}
    syms = sorted(S_)
    for a in syms:
        for b in syms:
            (L1, s1, p1), (L2, s2, p2) = a, b
            if not (a < b and s1 == s2 and p1 != p2 and abs(L1 - L2) <= 1 and L1 + L2 >= 1):
                continue
            d1, d2 = S_[a], S_[b]
            if not d1['match'] or not d2['match']:
                continue
            i1 = sorted(d1['match']); c1 = d1['v'][:, [d1['match'][j] for j in i1]]
            i2 = sorted(d2['match']); c2 = d2['v'][:, [d2['match'][j] for j in i2]]
            D = dipole_matrix(B, d1['cf'], L1, d2['cf'], L2, s1, True)
            D0 = dipole_matrix(B, d1['cf'], L1, d2['cf'], L2, s1, False)
            DV = dipole_matrix(B, d1['cf'], L1, d2['cf'], L2, s1, gauge='velocity')
            R, R0, RV = c1.T @ D @ c2, c1.T @ D0 @ c2, c1.T @ DV @ c2
            Rabs = np.abs(c1).T @ np.abs(D) @ np.abs(c2)
            e1 = [d1['e'][d1['match'][j]] for j in i1]
            e2 = [d2['e'][d2['match'][j]] for j in i2]
            for x, j1 in enumerate(i1):
                for y, j2 in enumerate(i2):
                    de = e2[y] - e1[x]          # model CI transition energy (hartree), ket - bra
                    out[(a, j1, b, j2)] = dict(rme=R[x, y], rme_noop=R0[x, y], rme_vel_eq=RV[x, y] / de if de else 0.0,
                                               cf=abs(R[x, y]) / Rabs[x, y] if Rabs[x, y] else 0.0)
    return out


# ---------------------------------------------------------------- spin-orbit / intermediate coupling (v67)
ALPHA_FS = 1 / 137.035999


def xi_radial(B, l):
    """<i|xi_l|k> on the CI orbitals of l, xi_l(r) = (alpha^2/2)(1/r) dV_l/dr of OUR fitted core potential (hartree)."""
    if l == 0:
        return None
    V = M.potential(B.r, l, B.params[min(l, 3)])
    xi = ALPHA_FS ** 2 / 2 * np.gradient(V, B.x) / B.r ** 2
    z = B.orb[l][1]
    return z.T @ (z * xi[:, None]) * B.dx


@lru_cache(None)
def sixj_h(a, b, c, d, e, f):
    """6j with half-integers given as Fractions."""
    return M.sixj(a, b, c, d, e, f)


def spin_rme(S1, S2, electron):
    """<(1/2 1/2)S1||s(i)||(1/2 1/2)S2>."""
    h = Fr(1, 2)
    six = sixj_h(h, Fr(S1), h, Fr(S2), h, Fr(1))
    ph = (-1) ** S2 if electron == 1 else (-1) ** S1
    return ph * math.sqrt((2 * S1 + 1) * (2 * S2 + 1)) * six * math.sqrt(1.5)


def so_matrix(B, xr, cf1, L1, S1, cf2, L2, S2):
    """Configuration-level <cf1 L1 S1 || sum_i xi(r_i) u^(11)(i) || cf2 L2 S2> (orbital x spin double tensor)."""
    fac = math.sqrt((2 * L1 + 1) * (2 * L2 + 1))
    sp1, sp2 = spin_rme(S1, S2, 1), spin_rme(S1, S2, 2)

    def lred(l):
        return math.sqrt(l * (l + 1) * (2 * l + 1))

    def prim(la, ia, lb, ib, lc, ic, ld, id_):
        v = 0.0
        if lb == ld and ib == id_ and la == lc and la > 0:
            v += (-1) ** (la + lb + L2 + 1) * fac * sixj(la, L1, lb, L2, lc, 1) * lred(la) * xr[la][ia, ic] * sp1
        if la == lc and ia == ic and lb == ld and lb > 0:
            v += (-1) ** (la + ld + L1 + 1) * fac * sixj(lb, L1, la, L2, ld, 1) * lred(lb) * xr[lb][ib, id_] * sp2
        return v
    n1 = [1 / math.sqrt(2 * (1 + ((-1) ** (S1 + la + lb + L1)) * (la == lb and ia == ib))) for la, ia, lb, ib in cf1]
    n2 = [1 / math.sqrt(2 * (1 + ((-1) ** (S2 + la + lb + L2)) * (la == lb and ia == ib))) for la, ia, lb, ib in cf2]
    by = {}
    for j, (lc, ic, ld, id_) in enumerate(cf2):
        by.setdefault((lc, ld, ic), []).append(j)
        by.setdefault((lc, ld, id_), []).append(j)
    W = np.zeros((len(cf1), len(cf2)))
    for i, (la, ia, lb, ib) in enumerate(cf1):
        cand = set()
        for key in ((la, lb, ia), (la, lb, ib), (lb, la, ia), (lb, la, ib)):
            cand.update(by.get(key, []))
        for j in cand:
            lc, ic, ld, id_ = cf2[j]
            eta = (-1) ** (S2 + lc + ld + L2)
            W[i, j] = 2 * n1[i] * n2[j] * (prim(la, ia, lb, ib, lc, ic, ld, id_) + eta * prim(la, ia, lb, ib, ld, id_, lc, ic))
    return W


def ion_xi_scales(B, min_fs_cm=0.1):
    """Per-l spin-orbit scale calibrated on the ION ONLY (v68): s_l = median over n of NIST/model fine structure of the
    one-valence core spectrum (Mg II for Mg I), levels with NIST splitting > min_fs_cm. No neutral-atom information."""
    lev = M.nist_levels()
    x, r, dx = M.grid()
    out = {}
    for l in range(1, 4):
        p = B.params[min(l, 3)]
        V = M.potential(r, l, p)
        xi = ALPHA_FS ** 2 / 2 * np.gradient(V, x) / r ** 2
        nm = M.NMIN.get(l, l + 1)
        _, z, _ = M.solve(l, p, nm - l + 5)
        rat = []
        for n in range(nm, nm + 5):
            lab = f'{n}{"spdf"[l]}'
            if (lab, 2 * l + 1) not in lev or (lab, 2 * l - 1) not in lev:
                continue
            nist = abs(lev[(lab, 2 * l + 1)] - lev[(lab, 2 * l - 1)])
            if nist <= min_fs_cm:
                continue
            model = (z[:, n - l - 1] ** 2 * xi * dx).sum() * (l + 0.5) * M.RY_CM
            rat.append(dict(level=lab, model_cm=round(model, 4), nist_cm=round(nist, 4), ratio=round(model / nist, 4)))
        if rat:
            out[l] = dict(scale=round(1 / st.median(x_['ratio'] for x_ in rat), 4), levels=rat)
    for l in range(4, LMAX + 1):
        out[l] = dict(scale=out.get(3, {}).get('scale', 1.0), levels=[], note='l>=4 takes the f scale (no NIST fine structure)')
    return out


def combine_so(so_l, scales):
    """Spin-orbit is linear in each l's xi: so = sum_l s_l so_l."""
    out = {}
    for l, d in so_l.items():
        for k, v in d.items():
            out[k] = out.get(k, 0.0) + scales.get(l, 1.0) * v
    return out


def ls_so_by_l(B, S_):
    """ls_so split by the l of the spin-orbit electron."""
    full = {l: xi_radial(B, l) for l in range(1, LMAX + 1)}
    out = {}
    for l in range(1, LMAX + 1):
        xr = {k: (v if k == l else np.zeros_like(v)) for k, v in full.items()}
        out[l] = ls_so(B, S_, xr)
    return out


def ls_so(B, S_, xr=None):
    """<T1 L S||V(11)||T2 L' S'> between all matched NIST terms of the same parity (|dL|,|dS| <= 1)."""
    if xr is None:
        xr = {l: xi_radial(B, l) for l in range(1, LMAX + 1)}
    out = {}
    syms = sorted(S_)
    for a in syms:
        for b in syms:
            (L1, s1, p1), (L2, s2, p2) = a, b
            if a > b or p1 != p2 or abs(L1 - L2) > 1 or abs(s1 - s2) > 1 or (s1 == 0 and s2 == 0) or (L1 == 0 and L2 == 0):
                continue
            d1, d2 = S_[a], S_[b]
            if not d1['match'] or not d2['match']:
                continue
            i1 = sorted(d1['match']); c1 = d1['v'][:, [d1['match'][j] for j in i1]]
            i2 = sorted(d2['match']); c2 = d2['v'][:, [d2['match'][j] for j in i2]]
            R = c1.T @ so_matrix(B, xr, d1['cf'], L1, s1, d2['cf'], L2, s2) @ c2
            for x, j1 in enumerate(i1):
                for y, j2 in enumerate(i2):
                    out[(a, j1, b, j2)] = R[x, y]
    return out


def ic_blocks(S_, so, xiscale=1.0):
    """Intermediate coupling per (J, parity). Diagonal = NIST LS term centroid; spin-orbit from our xi (x xiscale).
    Returns {(J, p): dict(basis=[(sym, j)], C (columns assigned to basis labels), E_model, E_nist, labels)}."""
    from scipy.optimize import linear_sum_assignment
    basis = {}
    for sym, d in S_.items():
        L, S, p = sym
        for j in sorted(d['match']):
            for J in range(abs(L - S), L + S + 1):
                basis.setdefault((J, p), []).append((sym, j))
    out = {}
    for (J, p), bs in basis.items():
        n = len(bs)
        H = np.zeros((n, n))
        for i, (a, j1) in enumerate(bs):
            H[i, i] = S_[a]['terms'][j1][0] / M.RY_CM
            for k in range(i, n):
                b, j2 = bs[k]
                key = (a, j1, b, j2) if (a, j1) <= (b, j2) else (b, j2, a, j1)
                if key not in so:
                    continue
                (La, Sa, _), (Lb, Sb, _) = (a, b) if key[0] == a and key[1] == j1 else (b, a)
                v = (-1) ** (Lb + Sa + J) * sixj(La, Sa, J, Sb, Lb, 1) * so[key] * xiscale
                H[i, k] += v
                if k != i:
                    H[k, i] += v
        e, C = np.linalg.eigh(H)
        _, col = linear_sum_assignment(-(C ** 2))           # basis label i <- eigenvector col[i]
        labels = [f"{S_[a]['terms'][j][1]}{J}" for a, j in bs]
        enist = [S_[a]['terms'][j][2].get(J) for a, j in bs]
        out[(J, p)] = dict(basis=bs, C=C[:, col], E_model=e[col] * M.RY_CM, E_nist=enist, labels=labels,
                           purity=(C[:, col] ** 2).max(0))
    return out


def term_tuned_so(S_, so, n_iter=3, min_fs_cm=0.05):
    """Energy fine-tuning per LS term (cf. Brage & Hibbert 1989, as Pehlivan Rhodin 2017 did for 4f): scale every
    spin-orbit matrix element by sqrt(f_i f_j), f_T chosen so each triplet's computed fine-structure width matches NIST.
    Uses Mg I ENERGIES only (no gf). Singlets take the factor of the triplet of the same configuration. Terms whose NIST
    splitting is unresolved, inverted or < min_fs_cm keep f = 1. Returns (so, factors)."""
    f = {}
    cur = dict(so)
    for _ in range(n_iter):
        fs = fs_report(ic_blocks(S_, cur, 1.0))
        for term, d in fs.items():
            if d['nist_cm'] >= min_fs_cm and d['model_cm'] > 0:
                f[term] = f.get(term, 1.0) * d['nist_cm'] / d['model_cm']

        def fac(sym, j):
            lab = S_[sym]['terms'][j][1]
            if lab in f:
                return f[lab]
            conf = lab.split(' ')[0]
            trip = [v for k, v in f.items() if k.split(' ')[0] == conf]
            return trip[0] if trip else 1.0
        cur = {k: v * math.sqrt(fac(k[0], k[1]) * fac(k[2], k[3])) for k, v in so.items()}
    return cur, {k: round(v, 4) for k, v in f.items()}


def ic_lines(S_, mult, blocks):
    """gf between intermediate-coupling states on NIST J-level energies; CF = radial CF of the dominant LS
    multiplet x mixing cancellation |sum| / sum|.|."""
    rows = []
    keys = sorted(blocks)
    for k1 in keys:
        for k2 in keys:
            (J1, p1), (J2, p2) = k1, k2
            if p1 >= p2 or abs(J1 - J2) > 1 or J1 + J2 == 0:
                continue
            b1, b2 = blocks[k1], blocks[k2]
            n1, n2 = len(b1['basis']), len(b2['basis'])
            D, D0, DV, Dcf = (np.zeros((n1, n2)) for _ in range(4))
            for i, (a, ja) in enumerate(b1['basis']):
                La, Sa, _ = a
                for k, (b, jb) in enumerate(b2['basis']):
                    Lb, Sb, _ = b
                    if Sa != Sb or abs(La - Lb) > 1:
                        continue
                    if (a, ja, b, jb) in mult:
                        m = mult[(a, ja, b, jb)]
                        f = (-1) ** (La + Sa + J2 + 1) * math.sqrt((2 * J1 + 1) * (2 * J2 + 1)) * sixj(La, J1, Sa, J2, Lb, 1)
                    elif (b, jb, a, ja) in mult:
                        m = mult[(b, jb, a, ja)]
                        f = (-1) ** (J1 - J2) * (-1) ** (Lb + Sa + J1 + 1) * math.sqrt((2 * J1 + 1) * (2 * J2 + 1)) * sixj(Lb, J2, Sa, J1, La, 1)
                    else:
                        continue
                    D[i, k], D0[i, k], DV[i, k], Dcf[i, k] = f * m['rme'], f * m['rme_noop'], f * m['rme_vel_eq'], m['cf']
            C1, C2 = b1['C'], b2['C']
            R, R0, RV = C1.T @ D @ C2, C1.T @ D0 @ C2, C1.T @ DV @ C2
            Rabs = np.abs(C1).T @ np.abs(D) @ np.abs(C2)
            for x in range(n1):
                for y in range(n2):
                    E1, E2 = b1['E_nist'][x], b2['E_nist'][y]
                    if E1 is None or E2 is None or E1 == E2 or not Rabs[x, y]:
                        continue
                    (lo, El, xl, bl), (up, Eu, xu, bu) = sorted(((b1['labels'][x], E1, x, b1), (b2['labels'][y], E2, y, b2)), key=lambda t: t[1])
                    de = (Eu - El) / M.RY_CM
                    gf, gf0, gfv = 2 / 3 * de * R[x, y] ** 2, 2 / 3 * de * R0[x, y] ** 2, 2 / 3 * de * RV[x, y] ** 2
                    if gf <= 1e-30:
                        continue
                    # dominant LS pair: largest |c1 c2 D| contribution
                    contrib = np.abs(C1[:, x][:, None] * D * C2[:, y][None, :])
                    i, k = np.unravel_index(contrib.argmax(), contrib.shape)
                    cf = Dcf[i, k] * abs(R[x, y]) / Rabs[x, y]
                    rows.append(dict(lower=lo, upper=up, E_lower=El, E_upper=Eu, wl_vac_A=round(1e8 / (Eu - El), 3),
                                     loggf=round(math.log10(gf), 4), loggf_no_polarization_operator=round(math.log10(gf0), 4) if gf0 > 0 else None,
                                     loggf_velocity_bare=round(math.log10(gfv), 4) if gfv > 0 else None,
                                     dT_gauge=round(abs(gf0 - gfv) / max(gf0, gfv), 4) if max(gf0, gfv) > 0 else None,
                                     cf=round(cf, 4), purity_lower=round(float(bl['purity'][xl]), 4), purity_upper=round(float(bu['purity'][xu]), 4),
                                     spin_change=lo.split(' ')[-1][0] != up.split(' ')[-1][0],
                                     nmax=max(nmax2(lo.split(' ')[0]), nmax2(up.split(' ')[0]))))
    return sorted(rows, key=lambda r: r['wl_vac_A'])


def fs_report(blocks):
    """Model vs NIST fine-structure splitting per term (max - min over J)."""
    t = {}
    for (J, p), b in blocks.items():
        for lab, em, en in zip(b['labels'], b['E_model'], b['E_nist']):
            t.setdefault(lab.rstrip('0123456789'), {})[J] = (em, en)
    out = {}
    for term, d in t.items():
        if len(d) < 2 or any(v[1] is None for v in d.values()):
            continue
        js = sorted(d)
        dm, dn = d[js[-1]][0] - d[js[0]][0], d[js[-1]][1] - d[js[0]][1]
        if abs(dn) > 0.05:
            out[term] = dict(model_cm=round(dm, 3), nist_cm=round(dn, 3), ratio=round(dm / dn, 3))
    return out


def ci_conf(label):
    """'3s.4p 3P*' -> ('3s.4p', '3P*')."""
    c, t = label.split(' ')
    return c, t


def nmax2(conf):
    return max(int(m) for m in re.findall(r'(\d+)[spdfgh]', conf))


def fs_lines(S_, mult):
    """Fine-structure lines on NIST J-level energies: gf = (2/3) dE (2J+1)(2J'+1) {L J S; J' L' 1}^2 |<L||D||L'>|^2."""
    rows = []
    for (a, j1, b, j2), m in mult.items():
        (L1, s, _), (L2, _, _) = a, b
        t1, t2 = S_[a]['terms'][j1], S_[b]['terms'][j2]
        for J1, E1 in t1[2].items():
            for J2, E2 in t2[2].items():
                if abs(J1 - J2) > 1 or J1 + J2 == 0 or E1 == E2:
                    continue
                fr = (2 * J1 + 1) * (2 * J2 + 1) * sixj(L1, J1, s, J2, L2, 1) ** 2
                if fr <= 0:
                    continue
                lo, up = ((t1, J1, E1), (t2, J2, E2)) if E1 < E2 else ((t2, J2, E2), (t1, J1, E1))
                de = (up[2] - lo[2]) / M.RY_CM
                gf, gf0 = (2 / 3 * de * fr * x ** 2 for x in (m['rme'], m['rme_noop']))
                if gf <= 0:
                    continue
                rows.append(dict(lower=f"{lo[0][1]}{lo[1]}", upper=f"{up[0][1]}{up[1]}", E_lower=lo[2], E_upper=up[2],
                                 wl_vac_A=round(1e8 / (up[2] - lo[2]), 3), loggf=round(math.log10(gf), 4),
                                 loggf_no_polarization_operator=round(math.log10(gf0), 4) if gf0 > 0 else None,
                                 cf=round(m['cf'], 4), nmax=max(nmax2(lo[0][1].split(' ')[0]), nmax2(up[0][1].split(' ')[0]))))
    return sorted(rows, key=lambda r: r['wl_vac_A'])


def air_to_vac(air):
    s2 = (1e4 / air) ** 2
    return air * (1 + 8.336624212083e-5 + 2.408926869968e-2 / (1.301065924522e2 - s2) + 1.599740894897e-4 / (3.892568793293e1 - s2))


def lab_rows(sp, conf_prefix):
    """Pehlivan Rhodin 2017: labels '3p 3P0' = 3s.3p 3P0; independent = lab lifetime (tau_theory False)."""
    out = []
    for r in csv.DictReader(l for l in open(sp['lab']) if not l.startswith('#')):
        def key(s):
            c, t = s.split(' ')
            full = '.' in c or re.fullmatch(r'\d+[spdfgh]2', c) or c == conf_prefix     # '3p2', '3p.3d', '3s2' are full configs
            conf = c if full else f'{conf_prefix}.{c}'
            odd = sum('spdfgh'.index(m[1]) * int(m[2] or 1) for m in re.finditer(r'\d+([spdfgh])(\d?)', conf)) % 2
            return conf, t[:-1] + ('*' if odd else ''), int(t[-1])
        out.append(dict(r, lo=key(r['lower']), up=key(r['upper']), loggf_exp=float(r['loggf_exp']), e=float(r['e_loggf_dex']),
                        calc=float(r['loggf_calc']), independent=r['tau_theory'] == 'False'))
    return out


def family2(r):
    """Jumping-electron family from the leading configurations, e.g. '3s.3p'->'3s.4s' = p-s."""
    a, b = r['lower'].split(' ')[0].split('.'), r['upper'].split(' ')[0].split('.')
    da = [x for x in a if x not in b] or a
    db = [x for x in b if x not in a] or b
    return f"{da[-1][-1]}-{db[-1][-1]}"


def nist_j_levels(sp):
    out = []
    for r in csv.DictReader(open(sp['levels'])):
        c = {k: (v or '').strip('="[]') for k, v in r.items() if k}
        if c['Level (cm-1)'] and re.fullmatch(r'\d+', c['J'] or ''):
            conf = re.sub(r'^' + re.escape(sp['corepfx']) + r'\.', '', c['Configuration'])
            out.append((conf, c['Term'], int(c['J']), float(c['Level (cm-1)'].replace('?', '').replace('+x', ''))))
    return out


def ls_broken(sp):
    """LS-broken configurations (pair coupling): singlet-triplet separation <= 3 x triplet fine-structure width, or
    < 1 cm-1 (unresolved). There S is not a good quantum number; the gf sum over S and J of the configuration is."""
    tl = {}
    for conf, term, j, e in nist_j_levels(sp):
        tl.setdefault(conf, {}).setdefault(term[0], {})[j] = e
    broken = set()
    for conf, d in tl.items():
        if '1' in d and '3' in d:
            e1 = st.mean(d['1'].values())
            e3 = sum((2 * j + 1) * e for j, e in d['3'].items()) / sum(2 * j + 1 for j in d['3'])
            fs = max(d['3'].values()) - min(d['3'].values())
            if abs(e1 - e3) <= max(3 * fs, 1.0):
                broken.add(conf)
    return broken


def so_variants(B, S_, sp):
    """All spin-orbit treatments (v68). Primary = 'ion' (per-l xi calibrated on the core ion only, declared).
    'abinitio' unscaled; 'tuned' = energy fine-tuning of every resolved triplet (Mg I energies, no gf);
    'hybrid' = ion everywhere + energy fine-tuning only on LS-broken configurations (post-hoc, REPORTED ONLY)."""
    so_l = ls_so_by_l(B, S_)
    sc = ion_xi_scales(B)
    ab = combine_so(so_l, {})
    ion = combine_so(so_l, {l: v['scale'] for l, v in sc.items()})
    tuned, tf = term_tuned_so(S_, ab)
    broken = ls_broken(sp)
    ion_t, _ = term_tuned_so(S_, ion)

    def conf(k, i):
        return S_[k[i]]['terms'][k[i + 1]][1].split(' ')[0]
    hyb = {k: (ion_t[k] if (conf(k, 0) in broken or conf(k, 2) in broken) else v) for k, v in ion.items()}
    return dict(ion=ion, abinitio=ab, tuned=tuned, hybrid=hyb), sc, tf


def sum_sigma(ms, gm):
    """sigma of a summed feature; components without a sigma (cancellation) carrying > 10 % of the gf -> no sigma."""
    if sum(10 ** x['loggf'] for x in ms if x['sigma_dex'] is None) > 0.1 * gm:
        return None
    return round(math.sqrt(sum((10 ** x['loggf'] * x['sigma_dex']) ** 2 for x in ms if x['sigma_dex'] is not None)) / gm, 4)


def solar_fill(sp, rows):
    """Every line-list row of the species: lower J-level by excitation energy, upper by lower + wavenumber;
    rows on the same (lower, upper) level pair are one feature."""
    levels = []
    for r in csv.DictReader(open(sp['levels'])):
        c = {k: (v or '').strip('="[]') for k, v in r.items() if k}
        if c['Level (cm-1)'] and re.fullmatch(r'\d+', c['J'] or ''):
            conf = re.sub(r'^' + re.escape(sp['corepfx']) + r'\.', '', c['Configuration'])
            par = sum('spdfgh'.index(m[1]) * int(m[2] or 1) for m in re.finditer(r'\d+([spdfgh])(\d?)', conf)) % 2
            levels.append((float(c['Level (cm-1)'].replace('?', '').replace('+x', '')), f"{conf} {c['Term']}{c['J']}", par, int(c['J'])))
    by = {(r['lower'], r['upper']): r for r in rows}
    broken = ls_broken(sp)

    def ckey(lab):
        c = lab.split(' ')[0]
        return c + ' {LSmix}' if c in broken else lab
    feats = {}
    for r in csv.DictReader(open(ROOT / 'data/linelists/linelist_solar.csv')):
        if (r['element'], r['ion']) != sp['linelist']:
            continue
        air = float(r['wavelength_air_A'])
        vac = air_to_vac(air) if air > 2000 else air
        elo = float(r['excitation_potential_eV']) * V61.HC_EV_CM
        best, cands = None, []
        for el, ll, pl, jl in levels:
            if abs(el - elo) > 20:
                continue
            for eu, lu, pu, ju in levels:
                if pu == pl or abs(ju - jl) > 1 or ju + jl == 0:       # E1 selection rules
                    continue
                c = abs(eu - el - 1e8 / vac) + abs(el - elo) / 10
                if c < 3:
                    cands.append((c, ll, lu))
                if c < 3 and (best is None or c < best[0]):
                    best = (c, ll, lu)
        if best:
            # NIST J-levels that coincide (e.g. Mg I 7f 3F J=2,3,4) cannot be told apart: compare at the TERM level
            # (sum of gf over the indistinguishable upper J), never by guessing a J
            tied = {lu for c, ll, lu in cands if ll == best[1] and c <= best[0] + 0.3 and lu.rstrip('0123456789') == best[2].rstrip('0123456789')}
            if len(tied) > 1:
                best = (best[0], best[1], best[2].rstrip('0123456789') + '{' + ','.join(sorted(x[len(best[2].rstrip('0123456789')):] for x in tied)) + '}')
        # unassigned and above the highest level the box holds (hydrogenic n = BOX_N below the limit) -> beyond the box
        lim = sp['ip'] - (M.ZC - 1) ** 2 * M.R_INF / BOX_N ** 2       # outer-electron charge: 1 (neutral), 2 (Al II), 3 (Si III)
        forb = next((t for w, t in sp.get('forbidden', {}).items() if abs(w - air) < 0.01), None)
        key = (best[1], best[2]) if best else (f'?{forb}', f'{air}') if forb else ('?>box', f'{air}') if elo + 1e8 / vac > lim else ('?', f'{air}')
        if best and (best[1].split(' ')[0] in broken or best[2].split(' ')[0] in broken):
            key = (ckey(best[1]), ckey(best[2].split('{')[0]))
        f = feats.setdefault(key, dict(lower=key[0], upper=key[1], rows=[]))
        f['rows'].append(dict(air=air, gf=10 ** float(r['log_gf']), depth=float(r['central_depth'] or 0)))
    out = []
    for (lo, up), f in feats.items():
        g = sum(x['gf'] for x in f['rows'])
        if up.endswith('{LSmix}') or lo.endswith('{LSmix}'):     # configuration-level feature: sum over S and J
            ms = [r for r in rows if ckey(r['lower']) == lo and ckey(r['upper']) == up]
            if ms:
                gm = sum(10 ** x['loggf'] for x in ms)
                m = dict(loggf=round(math.log10(gm), 4), sigma_dex=sum_sigma(ms, gm),
                         cls=max((x['cls'] for x in ms), key='ABCDE'.index), sigma_basis=ms[0]['sigma_basis'] + '+LSmix_sum')
            else:
                m = None
        elif up.endswith('}'):          # term-level feature: model gf summed over the tied upper J
            base, js = up[:-1].split('{')
            ms = [by[(lo, base + j)] for j in js.split(',') if (lo, base + j) in by]
            if ms:
                gm = sum(10 ** x['loggf'] for x in ms)
                m = dict(loggf=round(math.log10(gm), 4), sigma_dex=sum_sigma(ms, gm),
                         cls=max((x['cls'] for x in ms), key='ABCDE'.index), sigma_basis=ms[0]['sigma_basis'] + '+J_tied')
            else:
                m = None
        else:
            m = by.get((lo, up))
        spin = not lo.startswith('?') and '{LSmix}' not in lo + up and lo.split(' ')[1][0] != up.split(' ')[1][0]
        out.append(dict(lower=lo, upper=up, air_A=round(sum(x['air'] * x['gf'] for x in f['rows']) / g, 3), n_linelist_rows=len(f['rows']),
                        max_depth=max(x['depth'] for x in f['rows']), linelist_loggf=round(math.log10(g), 3),
                        model_loggf=m['loggf'] if m else None, model_sigma_dex=m['sigma_dex'] if m else None,
                        model_class=m['cls'] if m else None, sigma_basis=m['sigma_basis'] if m else None,
                        status=('graded' if m.get('sigma_dex') is not None else 'ungraded_cancellation') if m else 'unassigned' if lo == '?'
                        else 'no_model_beyond_box' if lo == '?>box' else 'non_E1_forbidden' if lo in ('?E2', '?M2')
                        else 'no_model_beyond_box' if max(nmax2(lo.split(' ')[0]), nmax2(up.split(' ')[0].split('{')[0])) > BOX_N
                        else 'no_model_spin_change' if spin
                        else 'no_model_core_excited' if lo.split(' ')[0].startswith('2p5') or up.split(' ')[0].startswith('2p5') else 'no_model_not_in_ci',
                        d_model_minus_linelist=round(m['loggf'] - math.log10(g), 3) if m else None))
    global LS_BROKEN
    LS_BROKEN = sorted(broken, key=lambda c: (c[-1], nmax2(c)))
    return sorted(out, key=lambda x: x['air_A'])


LS_BROKEN = None


def main(name, cache=None):
    import pickle
    if cache and Path(cache).exists():
        sp, B, S_ = pickle.load(open(cache, 'rb'))
        M.use(sp['core'])
    else:
        sp, B = build(name)
        _, S_ = run_ci(sp, B)
        if cache:
            pickle.dump((sp, B, S_), open(cache, 'wb'))
    S_ = dict(S_)
    beyond = S_.pop('_beyond_box')
    out = ROOT / 'data/audit/rya1354' / sp['out']
    out.mkdir(parents=True, exist_ok=True)
    # 1. energies (pre-registered: every matched term within 1% of NIST, total two-electron energy)
    erep = {}
    for sym, d in S_.items():
        for j, i in sorted(d['match'].items()):
            erep[d['terms'][j][1]] = dict(rel=round(d['e'][i] / d['nist_e'][j] - 1, 5), d_cm=round((d['e'][i] - d['nist_e'][j]) * M.RY_CM, 1))
    unmatched = [t[1] for d in S_.values() for j, t in enumerate(d['terms']) if j not in d['match']]
    # 2. lines: intermediate coupling (v67) on the LS multiplets
    mult = multiplets(B, S_)
    sov, ion_sc, tune_f = so_variants(B, S_, sp)
    so = sov['ion']                                     # declared primary (v68)
    ls_rows = {(r['lower'], r['upper']): r for r in fs_lines(S_, mult)}
    # control: xi = 0 must reproduce the pure-LS gf exactly
    ctl = {(r['lower'], r['upper']): r['loggf'] for r in ic_lines(S_, mult, ic_blocks(S_, so, 0.0))}
    ctl_dev = max(abs(ctl[k] - v['loggf']) for k, v in ls_rows.items() if k in ctl)
    ctl_missing = sum(k not in ctl for k in ls_rows)
    blocks = ic_blocks(S_, so, 1.0)
    fsr = fs_report(blocks)
    fs_all = {v: fs_report(ic_blocks(S_, sov[v], 1.0)) for v in sov}
    alt = {v: {(r['lower'], r['upper']): r['loggf'] for r in ic_lines(S_, mult, ic_blocks(S_, sov[v], 1.0))} for v in sov if v != 'ion'}
    rows = ic_lines(S_, mult, blocks)
    var, var_excl = load_variants(sp, max(abs(v['rel']) for v in erep.values()))
    for n, x in var_excl.items():
        print(f'variant {n} EXCLUDED from sigma_conv: energy {x["energy_max_rel"]} > primary + {CONV_ENERGY_SLACK}')
    for r in rows:
        k = (r['lower'], r['upper'])
        for v, d in alt.items():
            r[f'loggf_so_{v}'] = d.get(k)
        r['loggf_pure_LS'] = ls_rows[k]['loggf'] if k in ls_rows else None
        # error budget (dex), one column per term
        r['err_so_dex'] = round(abs(alt['abinitio'][k] - r['loggf']), 4) if k in alt['abinitio'] else 0.5
        for grp in ('conv', 'phys'):
            ds = [abs(d[k] - r['loggf']) for n, d in var.items() if VARIANTS[n]['group'] == grp and k in d]
            r[f'err_{grp}_dex'] = round(max(ds), 4) if ds else None
            r[f'n_{grp}_variants'] = len(ds)
        r['xi_sens_dex'] = max(r['err_so_dex'], r['err_phys_dex'] or 0.0)     # largest physics sensitivity (declared)
    by = {(r['lower'], r['upper']): r for r in rows}
    lab = lab_rows(sp, sp['labprefix'])
    for x in lab:
        k = (f'{x["lo"][0]} {x["lo"][1]}{x["lo"][2]}', f'{x["up"][0]} {x["up"][1]}{x["up"][2]}')
        x['key'] = k
        x['spin_change'] = x['lo'][1][0] != x['up'][1][0]
        r = by.get(k)
        x['model'] = r['loggf'] if r else None
        if r:
            r['lab_loggf'], r['d_lab'], r['lab_independent'] = x['loggf_exp'], round(r['loggf'] - x['loggf_exp'], 4), x['independent']
            r['ref_calc_loggf'], r['d_ref_calc'] = x['calc'], round(r['loggf'] - x['calc'], 4)
    ind = [x for x in lab if x['independent'] and x['model'] is not None]
    allowed = [x for x in lab if not x['spin_change'] and x['model'] is not None]
    spinlab = [(dict(line=x['wavelength_A'], lab=x['loggf_exp'], model=x['model'], d=None if x['model'] is None else round(x['model'] - x['loggf_exp'], 4)),
                x['key']) for x in lab if x['spin_change']]
    if ind:
        off = st.mean(abs(x['model'] - x['loggf_exp']) for x in ind)
        unc = max(x['e'] for x in ind)
    else:     # v69 sister ions: no lab-independent referee -> the Mg I measured floor carries along the sequence (not refit on theory)
        off, unc = sp['floor_fallback'], 0.0
    fl = math.ceil(max(off, unc) * 1000 - 1e-9) / 1000
    fam = {}
    for x in allowed:
        fam.setdefault(family2(by[x['key']]), []).append(abs(x['model'] - x['loggf_exp']))
    allabs = sorted(abs(x['model'] - x['loggf_exp']) for x in allowed)
    abs68 = allabs[int(0.68 * (len(allabs) - 1))]
    famsig = {k: round(sorted(v)[int(0.68 * (len(v) - 1))], 4) for k, v in fam.items() if len(v) >= 5}
    def grade(cfb):
        for r in rows:
            half = abs(r['loggf'] - r['loggf_no_polarization_operator']) / 2 if r['loggf_no_polarization_operator'] is not None else 0.5
            s = math.hypot(fl, half, famsig.get(family2(r), abs68), r['xi_sens_dex'], r['err_conv_dex'] or 0.0)
            r['sigma_basis'] = 'lab_family' if family2(r) in famsig else 'lab_pooled'
            if r['cf'] < CF_NO_GRADE:
                r['sigma_dex'], r['cls'], r['sigma_basis'] = None, 'E', 'ungraded_cancellation'
                continue
            if r['cf'] < 0.1:          # cancellation: the CF-bin term measured vs the line list (v67)
                b = next(x for x in cfb if x['cf_lo'] <= r['cf'] < x['cf_hi'])
                s = math.hypot(s, b['abs68'])
                r['sigma_basis'] += '+cf_bin'
            r['sigma_dex'] = round(s, 4)
            r['cls'] = M.sigma_class(r['sigma_dex'], r['cf'])
    grade([dict(cf_lo=lo, cf_hi=hi, abs68=abs68) for lo, hi in CF_BINS2])     # provisional
    cfb = cf_bins(solar_fill(sp, rows), rows, carry='floor_fallback' in sp)
    grade(cfb)
    z = [abs(r['d_lab']) / r['sigma_dex'] for r in rows if 'd_lab' in r and r['sigma_dex'] and not r['spin_change']]
    with open(out / 'lines.csv', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=sorted({k for r in rows for k in r}), lineterminator='\n')
        w.writeheader()
        w.writerows(rows)
    fill = solar_fill(sp, rows)
    with open(out / 'solar_fill.csv', 'w', newline='') as f:      # may be empty (Si III: no rows in the solar line list)
        w = csv.DictWriter(f, fieldnames=list(fill[0]) if fill else ['lower', 'upper', 'air_A', 'status'], lineterminator='\n')
        w.writeheader()
        w.writerows(fill)
    meas = [x for x in fill if x['max_depth'] >= 0.01]
    summary = dict(schema='rya1354.twovalence-model.v2', species=name, status=STATUS, core=sp['core'], alpha_d=M.ALPHA_D,
                   grid=GRID, basis=dict(NION=NION, NRYD=NRYD, NCOMP=NCOMP, RCOMP=RCOMP, LMAX=LMAX, n_per_l={l: len(B.orb[l][0]) for l in B.orb}),
                   energy_gate=dict(max_abs_rel=max(abs(v['rel']) for v in erep.values()), passed=all(abs(v['rel']) <= 0.01 for v in erep.values()),
                                    terms=erep, unmatched_nist_terms=unmatched, beyond_box_n_gt=BOX_N, n_beyond_box=len(beyond),
                                    note='v65 first run aligned without the box cut: n=58/59p and 24/25d NIST terms were paired '
                                         'with box-continuum CI states (max 3.3%); n>BOX_N terms are now excluded, not force-matched'),
                   lab=[dict(line=x['wavelength_A'], lower=x['lower'], upper=x['upper'], lab=x['loggf_exp'], e_lab=x['e'], ref_calc=x['calc'],
                             model=x['model'], d_model_lab=None if x['model'] is None else round(x['model'] - x['loggf_exp'], 4),
                             independent=x['independent'], spin_change=x['spin_change']) for x in lab],
                   lab_independent=dict(n=len(ind), mean_abs=round(off, 4), median=round(st.median(x['model'] - x['loggf_exp'] for x in ind), 4) if ind else None,
                                        floor_source='lab-independent lines' if ind else f"carried from Mg I ({sp['floor_fallback']})"),
                   lab_spin_allowed=dict(n=len(allowed), median=round(st.median(x['model'] - x['loggf_exp'] for x in allowed), 4), abs68=round(abs68, 4)),
                   ref_calc_vs_lab_spin_allowed=dict(median=round(st.median(x['calc'] - x['loggf_exp'] for x in allowed), 4),
                                                    abs68=round(sorted(abs(x['calc'] - x['loggf_exp']) for x in allowed)[int(0.68 * (len(allowed) - 1))], 4)),
                   floor_dex=fl, family_abs68=famsig, sigma_coverage_vs_lab_in_sample=dict(n=len(z), within_1=round(sum(v <= 1 for v in z) / len(z), 3),
                                                                                          within_2=round(sum(v <= 2 for v in z) / len(z), 3)),
                   spin_orbit=dict(primary='ion: per-l xi calibrated on the core ion only (v68 declared)', ion_scales=ion_sc,
                                   control_xi0_max_dev_vs_LS=round(ctl_dev, 6), control_xi0_missing=ctl_missing,
                                   fine_structure={v: {k: d[k] for k in sorted(d, key=lambda t: -abs(d[t]['nist_cm']))[:12]} for v, d in fs_all.items()},
                                   term_tuning_factors=tune_f,
                                   lab_spin_change_lines=[dict(x, **{f'model_{v}': (alt[v].get(k_) if v != 'ion' else x['model']) for v in sov})
                                                          for x, k_ in spinlab],
                                   ref_4571=[dict(x, **{f'model_{v}': alt[v].get((x['lower'], x['upper'])) for v in alt}) for x in intercomb(sp, by)],
                                   cf_bins=cfb, cf_no_grade=CF_NO_GRADE),
                   variants={n: dict(VARIANTS[n], n_lines=len(d)) for n, d in var.items()},
                   variants_excluded={n: dict(VARIANTS[n], energy_max_rel=x['energy_max_rel'], rule=f'conv energy > primary + {CONV_ENERGY_SLACK}',
                                              p90_abs_dloggf=p90([abs(x['g'][(r['lower'], r['upper'])] - r['loggf']) for r in rows if (r['lower'], r['upper']) in x['g']]))
                                      for n, x in var_excl.items()},
                   gauge=gauge_summary(rows),
                   composition_vs_nist=composition_test(sp, S_, B.nmin),
                   n_lines=len(rows), class_counts={c: sum(r['cls'] == c for r in rows) for c in 'ABCDE'},
                   fill=dict(n_features=len(fill), n_linelist_rows=sum(x['n_linelist_rows'] for x in fill),
                             status_counts={k: sum(x['status'] == k for x in fill) for k in sorted({x['status'] for x in fill})},
                             n_measurable=len(meas), n_measurable_graded=sum(x['status'] == 'graded' for x in meas),
                             ls_broken_configurations=LS_BROKEN))
    (out / 'summary.json').write_text(json.dumps(summary, indent=1) + '\n')
    return summary


# v68 variants: each rebuilds the full CI with ONE knob changed. group conv = numerical convergence, phys = physics input.
VARIANTS = {
    'nryd18': dict(group='conv', NRYD=18), 'ncomp18': dict(group='conv', NCOMP=18), 'rcomp35': dict(group='conv', RCOMP=35.0),
    'rcomp20': dict(group='conv', RCOMP=20.0),     # v69
    'lmax5': dict(group='conv', LMAX=5), 'npts4800': dict(group='conv', GRID_NPTS=4800), 'box1200': dict(group='conv', GRID_X1=1200.0),
    'alpha2_p5': dict(group='phys', ALPHA2_SCALE=1.05), 'alpha2_m5': dict(group='phys', ALPHA2_SCALE=0.95),
    'rcdiel_p20': dict(group='phys', RC_DIEL_SCALE=1.2), 'rcdiel_m20': dict(group='phys', RC_DIEL_SCALE=0.8),
}


def variant_dir(sp):
    return ROOT / 'data/audit/rya1354' / sp['out'] / 'variants'


def run_variant(name, vname):
    """Rebuild the CI with one knob changed; write {lower|upper: loggf} of the primary (ion-calibrated) IC model."""
    global NRYD, NCOMP, RCOMP, LMAX, ALPHA2_SCALE, RC_DIEL_SCALE
    v = VARIANTS[vname]
    for k, x in v.items():
        if k == 'GRID_NPTS':
            GRID['NPTS'] = x
        elif k == 'GRID_X1':
            GRID['X1'] = math.log(x)
        elif k != 'group':
            globals()[k] = x
    sp, B = build(name)
    _, S_ = run_ci(sp, B)
    S_ = dict(S_)
    S_.pop('_beyond_box')
    mult = multiplets(B, S_)
    sov, sc, _ = so_variants(B, S_, sp)
    rows = ic_lines(S_, mult, ic_blocks(S_, sov['ion'], 1.0))
    erel = max(abs(d['e'][i] / d['nist_e'][j] - 1) for d in S_.values() for j, i in d['match'].items())
    d = variant_dir(sp)
    d.mkdir(parents=True, exist_ok=True)
    (d / f'{vname}.json').write_text(json.dumps(dict(variant=vname, knobs=v, energy_max_rel=round(erel, 5),
                                                     ion_scales={l: x['scale'] for l, x in sc.items()},
                                                     loggf={f"{r['lower']}|{r['upper']}": r['loggf'] for r in rows}), indent=0) + '\n')


CONV_ENERGY_SLACK = 0.0010   # v69 (A3): a conv variant enters sigma_conv only if its energy max error <= primary + 0.10 %


def load_variants(sp, primary_erel=None):
    """Variant loggf by name. With primary_erel, conv variants whose energy is worse than primary + CONV_ENERGY_SLACK are
    returned separately (reported, not used in sigma: they degrade the basis, so they are not convergence steps)."""
    out, excluded = {}, {}
    d = variant_dir(sp)
    for n in VARIANTS:
        f = d / f'{n}.json'
        if f.exists():
            j = json.loads(f.read_text())
            g = {tuple(k.split('|')): v for k, v in j['loggf'].items()}
            if primary_erel is not None and VARIANTS[n]['group'] == 'conv' and j['energy_max_rel'] > primary_erel + CONV_ENERGY_SLACK:
                excluded[n] = dict(g=g, energy_max_rel=j['energy_max_rel'])
            else:
                out[n] = g
    return (out, excluded) if primary_erel is not None else out


def p90(v):
    v = sorted(v)
    return round(v[int(0.9 * (len(v) - 1))], 4) if v else None


LETTER = 'spdfghik'


def parse_conf(conf, nmin, corepfx):
    """NIST configuration ('3p2', '3p.3d', '2p6.3s2', '3s.4f') -> CI key (la, ia, lb, ib), ia = n - nmin_l (ion-orbital index)."""
    parts = [p for p in conf.split('.') if p != corepfx]
    orbs = []
    for p in parts:
        m = re.fullmatch(r'(\d+)([spdfghik])(\d?)', p)
        if not m:                      # e.g. NIST 'nd' / '(2S)' notations: not a single CI configuration
            return None
        n, l, k = int(m[1]), LETTER.index(m[2]), m[3]
        orbs += [(l, n - nmin.get(l, l + 1))] * (int(k) if k else 1)
    if len(orbs) != 2:
        return None
    (la, ia), (lb, ib) = sorted(orbs)
    return (la, ia, lb, ib)


def composition_test(sp, S_, nmin):
    """v69 (A2): CI weight of the NIST second component in every LS state NIST lists one for (MgI_levels_percentages.csv).
    Report only (no fit). Percentages are basis-dependent (ours: ion orbitals), so agreement is approximate by construction."""
    f = sp.get('percentages') or ROOT / 'data/reference/onevalence_nist' / (sp['name'].replace(' ', '') + '_levels_percentages.csv')
    if not f.exists():
        return None
    clean = lambda x: x.replace('=', '').replace('"', '').strip()
    rows, seen = [], set()
    for r in csv.DictReader(open(f)):
        conf, term, perc = clean(r['Configuration']), clean(r['Term']), clean(r['Leading percentages'])
        if ':' not in perc or (conf, term) in seen:
            continue
        seen.add((conf, term))
        lead, rest = perc.split(':', 1)
        bits = rest.split()
        if len(bits) < 3:
            continue
        p2, conf2, term2 = float(bits[0]), bits[1], bits[2]
        if term2 != term:
            continue
        L = 'SPDFGHIK'.index(term.rstrip('*')[1])
        S = (int(term[0]) - 1) // 2
        par = 1 if term.endswith('*') else 0
        lab = conf.replace(sp['corepfx'] + '.', '') + ' ' + term
        d = S_.get((L, S, par))
        if not d:
            continue
        labels = [t[1] for t in d['terms']]
        if lab not in labels or labels.index(lab) not in d['match']:
            continue
        key = parse_conf(conf2, nmin, sp['corepfx'])
        if key is None:
            continue
        v = d['v'][:, d['match'][labels.index(lab)]]
        w = float(sum(v[k] ** 2 for k, c in enumerate(d['cf']) if c == key))
        pauli = key[0] == key[2] and key[1] == key[3] and (L + S) % 2 == 1      # l^2 allows only L + S even
        rows.append(dict(level=lab, second=conf2 + ' ' + term2, nist_pct=p2, model_pct=round(100 * w, 2),
                         lead_nist_pct=float(lead.split()[0]) if lead.split() else None, nist_component_pauli_forbidden=pauli))
    d = [r['model_pct'] - r['nist_pct'] for r in rows if not r['nist_component_pauli_forbidden']]
    return dict(n=len(rows), n_pauli_forbidden_in_nist=sum(r['nist_component_pauli_forbidden'] for r in rows), median_model_minus_nist_pct=round(sorted(d)[len(d) // 2], 1) if d else None,
                rows=rows, note='report only; percentages are basis-dependent')


def gauge_summary(rows):
    """v68: velocity (bare) vs length (bare). dT is a DIAGNOSTIC: on the lab lines all dT < 1 % and Spearman(dT, |miss|)
    = -0.08; vs the line list dT tracks |d| only through CF (rho 0.30 overall, -0.11 within CF >= 0.1)."""
    d = sorted(r['dT_gauge'] for r in rows if r.get('dT_gauge') is not None)
    return dict(n=len(d), median=d[len(d) // 2], p90=d[int(0.9 * len(d))], used_in_sigma=False,
                lab_test='Spearman(dT,|model-lab|)=-0.08 (26 lab lines, all dT<0.01): no power, no prediction')


CF_NO_GRADE = 0.02
CF_BINS2 = ((0.05, 0.1), (0.02, 0.05))


MG1_CF_BINS = (0.155, 0.34)   # Mg I v68 measured abs68 for CF 0.05-0.1 / 0.02-0.05 (carried to the sister ions when a bin is thin)


def cf_bins(fill, rows, carry=False):
    """abs68 of model - line list for single-line features (n <= 8, depth >= 0.01) per CF bin (diagnostic referee:
    the line list is theory, used as FF was for Na/Mg II). Spin-changing lines are excluded: VALD's Mg I intercombination
    gf were refuted in the observed Sun in v67 (20632/11812/9822/20867/7746 A predicted 0.13-0.42 deep, observed <= 0.03)."""
    L = {(r['lower'], r['upper']): r for r in rows}
    out = []
    for lo, hi in CF_BINS2:
        d = sorted(abs(x['d_model_minus_linelist']) for x in fill if x['status'] == 'graded' and x['max_depth'] >= 0.01
                   and (x['lower'], x['upper']) in L and L[(x['lower'], x['upper'])]['nmax'] <= 8 and lo <= L[(x['lower'], x['upper'])]['cf'] < hi
                   and not L[(x['lower'], x['upper'])]['spin_change'])
        out.append(dict(cf_lo=lo, cf_hi=hi, n=len(d), abs68=round(d[int(0.68 * (len(d) - 1))], 4) if len(d) >= 5 else None))
    for k in range(len(out)):          # thin bin -> the next-higher-CF bin's value, never smaller
        if out[k]['abs68'] is None and carry:
            out[k]['abs68'], out[k]['fallback'] = MG1_CF_BINS[k], 'carried from Mg I v68'
        elif out[k]['abs68'] is None:
            out[k]['abs68'], out[k]['fallback'] = out[k - 1]['abs68'] if k else 0.15, True
    return out


def intercomb(sp, by):
    """Reference intercombination lines (theory): KP08 Table 40 = Tachiev & Froese Fischer Breit-Pauli MCHF."""
    out = []
    for lo, up, ref, src in sp.get('intercomb_ref', ()):
        r = by.get((lo, up))
        out.append(dict(lower=lo, upper=up, ref_loggf=ref, ref=src, model=r['loggf'] if r else None, d=round(r['loggf'] - ref, 4) if r else None))
    return out


if __name__ == '__main__':
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    cache = next((a.split('=', 1)[1] for a in sys.argv[1:] if a.startswith('--cache=')), None)
    vname = next((a.split('=', 1)[1] for a in sys.argv[1:] if a.startswith('--variant=')), None)
    if vname:
        run_variant(args[0], vname)
        sys.exit(0)
    for name in (args or list(SPECIES)):
        s = main(name, cache)
        print(name, json.dumps({k: v for k, v in s.items() if k not in ('lab', 'energy_gate', 'grid', 'basis')}, indent=1, default=str))
        print('energy gate', s['energy_gate']['passed'], s['energy_gate']['max_abs_rel'], 'unmatched', len(s['energy_gate']['unmatched_nist_terms']))
        for x in s['lab']:
            print(x)
