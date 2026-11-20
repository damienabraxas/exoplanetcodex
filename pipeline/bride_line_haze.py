"""Background metal-line opacity ("line haze") for the Bride, the way MULTI 2.3 adds it.

RYA-1282, 2026-11-20. The Gerber+2023 departure grid was computed with MULTI 2.3 "updated
as described in Bergemann et al. (2019)", whose Sect. 3.3 says MULTI was given background
opacity tables from MARCS/Turbospectrum with detailed LINE opacity at > 1e5 wavelengths.
In MULTI 2.3 (mul23_opacu.f OPAC, OSMET/OSLINE; Collet 2002, Carlsson 2003):

  * the line opacity is added ONLY at the frequencies of the active atom's bound-free
    transitions (KR > NLINE); the atom's own line transitions get none;
  * it is split into absorption and scattering with the van Regemorter thermalisation
    fraction  eps = ne / (C + ne),  C = sqrt(T) / (20.6 lambda_cm^3) *
    (EOSMET / P_ion(beta) + (1 - EOSMET) / P_neutral(beta)),  beta = hc / (lambda k T),
    EOSMET = 0.5 (the code default); P from PSEMI (Sobel'man tables, natural spline).

The Bride cannot read MULTI's ABSMET table (none is published), so the line opacity is
computed here directly, in LTE, from VALD "extract stellar" lists: Voigt profiles, Irwin
partition functions and ionisation energies (TS's own file, as for the bound-free terms),
Saha ionisation from the model's own electron density, model abundances (or any other
mix, e.g. solar, through `abundances_dex`). Oxygen (the active atom) and hydrogen (its own
treatment) are excluded, as MULTI's table omits the element being solved.
Molecular lines are NOT included yet (the lists carry CO, CN, CH, OH, C2, MgH...).
Units: cgs internally, the provider gets 1/m.
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np

_PI_E2_MC = 0.026540081      # pi e^2 / (m_e c), cm^2 s^-1
_K_CGS = 1.380649e-16
_H_CGS = 6.62607015e-27
_C_CGS = 2.99792458e10
_AMU = 1.66053907e-24
_A0 = 5.29177e-9
_EV = 1.602176634e-12

_ROW = re.compile(r"^'([A-Za-z][A-Za-z0-9]*) (\d+)',\s*([-+0-9.Ee]+),\s*([-+0-9.Ee]+),\s*([-+0-9.Ee]+),"
                  r"\s*([-+0-9.Ee]+),\s*([-+0-9.Ee]+),\s*([-+0-9.Ee]+),\s*([-+0-9.Ee]+),\s*([-+0-9.Ee]+),"
                  r"\s*([-+0-9.Ee]+),\s*([-+0-9.Ee]+),\s*([-+0-9.Ee]+),\s*([-+0-9.Ee]+),")


def read_vald_stellar(paths):
    """Atomic transitions from VALD 'extract stellar' long-format files.

    Returns dict of arrays: Z, ion (0 = neutral), wave_A (as VALD gives it: air above
    2000 A, vacuum below), loggf, elow_eV, eup_eV, rad, stark, waals, and the count of
    molecular rows skipped."""
    from lightweaver.atomic_table import PeriodicTable
    symbols = {}
    for z in range(1, 93):
        try:
            symbols[PeriodicTable[z].name] = z
        except Exception:                                            # noqa: BLE001
            pass
    cols = {k: [] for k in ('Z', 'ion', 'wave_A', 'loggf', 'elow_eV', 'eup_eV', 'rad', 'stark', 'waals')}
    molecular = 0
    for p in paths:
        with open(p, encoding='latin-1') as f:
            for line in f:
                if not line.startswith("'") or line.startswith("' ") or line.startswith("'_"):
                    continue
                m = _ROW.match(line)
                if m is None:
                    continue
                sp = m.group(1)
                if sp not in symbols:
                    molecular += 1
                    continue
                v = [float(m.group(i)) for i in range(3, 15)]
                cols['Z'].append(symbols[sp]); cols['ion'].append(int(m.group(2)) - 1)
                cols['wave_A'].append(v[0]); cols['loggf'].append(v[1]); cols['elow_eV'].append(v[2])
                cols['eup_eV'].append(v[4]); cols['rad'].append(v[9]); cols['stark'].append(v[10])
                cols['waals'].append(v[11])
    out = {k: np.asarray(x, dtype=int if k in ('Z', 'ion') else float) for k, x in cols.items()}
    out['molecular_rows_skipped'] = molecular
    return out


def psemi(z_ion, beta):
    """MULTI PSEMI: thermal P function for semi-empirical (van Regemorter) rates."""
    from scipy.interpolate import CubicSpline
    betref = np.array([-2.0, -1.699, -1.398, -1.0, -0.699, -0.398, 0.0, 0.301, 0.602, 1.0])
    ref = (np.array([1.160, 0.956, 0.758, 0.493, 0.331, 0.209, 0.100, 0.063, 0.040, 0.023]) if z_ion == 1 else
           np.array([1.160, 0.977, 0.788, 0.554, 0.403, 0.290, 0.214, 0.201, 0.200, 0.200]))
    b = np.asarray(beta, dtype=float)
    mid = CubicSpline(betref, ref, bc_type='natural')(np.log10(np.clip(b, 0.01, 10.)))
    lo = 0.27566 * (0.577 + np.log(np.maximum(b, 1e-300)))
    hi = 0.066 / np.sqrt(np.maximum(b, 1e-300)) if z_ion == 1 else np.full_like(b, 0.200)
    return np.where(b < 0.01, lo, np.where(b > 10., hi, mid))


def multi_line_epsilon(wave_A, T, ne_cm3, eosmet=0.5):
    """MULTI OPAC: thermalisation fraction of the background line opacity, (nwave, ndep)."""
    w = np.asarray(wave_A, dtype=float)[:, None]
    T = np.asarray(T, dtype=float)[None, :]
    beta = _H_CGS * _C_CGS * 1e8 / w / _K_CGS / T
    c1 = np.sqrt(T) / (20.6 * (w * 1e-8) ** 3)
    C = c1 * (eosmet / psemi(2, beta) + (1. - eosmet) / psemi(1, beta))
    ne = np.asarray(ne_cm3, dtype=float)[None, :]
    return ne / (C + ne)


def _stage_fractions(irwin_path, Z, T, ne_cm3, nstage=3):
    """Saha ionisation fractions (stages 0..nstage-1) and partition functions."""
    from pipeline.ts_jonabs import irwin_ionisation_eV, irwin_partition
    U = []
    for s in range(nstage):
        try:
            U.append(irwin_partition(irwin_path, Z, s, T))
        except Exception:                                            # noqa: BLE001
            break
    phi = 2. * (2 * np.pi * 9.1093837e-28 * _K_CGS * T / _H_CGS ** 2) ** 1.5
    rel = [np.ones_like(T)]
    for s in range(len(U) - 1):
        chi = irwin_ionisation_eV(irwin_path, Z, s) * _EV
        rel.append(rel[-1] * phi * U[s + 1] / U[s] * np.exp(-chi / (_K_CGS * T)) / ne_cm3)
    rel = np.array(rel)
    return rel / rel.sum(axis=0), np.array(U)


def line_haze(wave_A, T, ne_cm3, nH_cm3, nHI_cm3, vturb_cm_s, abundances_dex, lines, irwin_path,
              exclude_Z=(1, 8), window_A=3., strong_window_A=60.):
    """LTE metal-line opacity (cm^-1), (nwave, ndep), at the given wavelengths (VALD's convention).

    `abundances_dex` is A(X) = log N(X)/N(H) + 12 for Z = 1..len."""
    from scipy.special import wofz, gamma as Gamma
    from pipeline.ts_jonabs import irwin_ionisation_eV
    w = np.asarray(wave_A, dtype=float)
    T = np.asarray(T, dtype=float); ne = np.asarray(ne_cm3, dtype=float)
    nH = np.asarray(nH_cm3, dtype=float); nHI = np.asarray(nHI_cm3, dtype=float)
    vt = np.broadcast_to(np.asarray(vturb_cm_s, dtype=float), T.shape)
    keep = ~np.isin(lines['Z'], exclude_Z) & (lines['ion'] <= 1)   # I and II; 60 III rows in the solar lists
    L = {k: v[keep] for k, v in lines.items() if isinstance(v, np.ndarray)}
    out = np.zeros((len(w), len(T)))
    order = np.argsort(L['wave_A']); L = {k: v[order] for k, v in L.items()}
    lw0 = L['wave_A']
    # per-line, per-depth line strength S (cm^-1 Hz) and damping, by element/stage
    nl, nd = len(lw0), len(T)
    S = np.zeros((nl, nd)); dnuD = np.zeros((nl, nd)); gam = np.zeros((nl, nd))
    nu0 = _C_CGS / (lw0 * 1e-8)
    from lightweaver.atomic_table import PeriodicTable
    for Z in np.unique(L['Z']):
        if Z - 1 >= len(abundances_dex):
            continue
        sel_Z = L['Z'] == Z
        frac, U = _stage_fractions(irwin_path, int(Z), T, ne)
        Nel = 10. ** (abundances_dex[Z - 1] - 12.) * nH
        mass = PeriodicTable[int(Z)].mass * _AMU
        for s in range(len(U)):
            sel = sel_Z & (L['ion'] == s)
            if not sel.any():
                continue
            idx = np.where(sel)[0]
            n_over_U = (Nel * frac[s] / U[s])[None, :]
            boltz = np.exp(-L['elow_eV'][idx, None] * _EV / (_K_CGS * T[None, :]))
            stim = -np.expm1(-_H_CGS * nu0[idx, None] / (_K_CGS * T[None, :]))
            S[idx] = _PI_E2_MC * 10. ** L['loggf'][idx, None] * n_over_U * boltz * stim
            dnuD[idx] = nu0[idx, None] / _C_CGS * np.sqrt(2 * _K_CGS * T[None, :] / mass + vt[None, :] ** 2)
            # damping (FWHM, rad/s): radiative + quadratic Stark + van der Waals (VALD conventions)
            rad = L['rad'][idx]
            g = np.where(rad != 0, 10. ** rad, 2.223e15 / lw0[idx] ** 2)[:, None] * np.ones((1, nd))
            st = L['stark'][idx]
            g += np.where(st != 0, 10. ** st, 0.)[:, None] * ne[None, :] * (T[None, :] / 1e4) ** (1. / 6.)
            wa = L['waals'][idx]
            vbar = np.sqrt(8 * _K_CGS * T / np.pi * (1. / (1.008 * _AMU) + 1. / mass))
            g6 = np.zeros((len(idx), nd))
            neg = wa < 0
            g6[neg] = 10. ** wa[neg, None] * (T[None, :] / 1e4) ** 0.3
            abo = wa > 0
            if abo.any():
                sig = np.floor(wa[abo]) * _A0 ** 2; alp = wa[abo] - np.floor(wa[abo])
                v0 = 1e6
                g6[abo] = ((4. / np.pi) ** (alp[:, None] / 2.) * Gamma((4. - alp[:, None]) / 2.) * v0
                           * sig[:, None] * (vbar[None, :] / v0) ** (1. - alp[:, None])) * 2.
            uns = wa == 0
            if uns.any():
                chi = irwin_ionisation_eV(irwin_path, int(Z), s)
                eu = np.minimum(L['eup_eV'][idx][uns], chi - 0.1); el = np.minimum(L['elow_eV'][idx][uns], chi - 0.2)
                c6 = 0.3e-30 * np.abs(1. / (chi - eu) ** 2 - 1. / (chi - el) ** 2)
                g6[uns] = 17. * vbar[None, :] ** 0.6 * c6[:, None] ** 0.4
            g += g6 * nHI[None, :]
            gam[idx] = g
    # Lines whose damping wings still matter beyond `window_A` (Mg II h&k, Ca II H&K, Mg I b,
    # ...) are summed out to `strong_window_A`. Criterion: the Lorentz wing at the window edge,
    # kappa = S gam / (4 pi^2 dnu^2), exceeds 1e-9 cm^-1 at some depth (photospheric
    # continua are 1e-10..1e-6 cm^-1).
    dnu_w = _C_CGS * window_A * 1e-8 / (lw0 * 1e-8) ** 2
    wing = (S * gam / (4 * np.pi ** 2 * dnu_w[:, None] ** 2)).max(axis=1)
    strong = wing > 1e-9                                   # cm^-1 at the window edge
    lo_i = np.searchsorted(lw0, w - strong_window_A); hi_i = np.searchsorted(lw0, w + strong_window_A)
    nu = _C_CGS / (w * 1e-8)
    for j in range(len(w)):
        i0, i1 = lo_i[j], hi_i[j]
        if i1 <= i0:
            continue
        sl = slice(i0, i1)
        near = (np.abs(lw0[sl] - w[j]) < window_A) | strong[sl]
        if not near.any():
            continue
        ii = np.arange(i0, i1)[near]
        v = (nu[j] - nu0[ii, None]) / dnuD[ii]
        a = gam[ii] / (4 * np.pi * dnuD[ii])
        H = wofz(v + 1j * a).real
        out[j] = (S[ii] * H / (np.sqrt(np.pi) * dnuD[ii])).sum(axis=0)
    return out


def default_vald_paths(root=None):
    root = Path(root) if root is not None else Path.home() / 'codex/bride_inputs/rya1282/vald_haze'
    return sorted(root.glob('vald_solar*_raw.txt'))
