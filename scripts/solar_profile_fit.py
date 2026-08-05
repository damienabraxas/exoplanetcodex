#!/usr/bin/env python3
"""
Shared in-window profile-fit machinery for the solar synthesis harnesses (RYA-643).

THE SINGLE SOURCE for the rest-frame correction and the broadening grid.

Why this module exists
----------------------
RYA-551 (Sr II) wrote a chi2 in-window profile fitter; RYA-560 (Zr II) copied it
"verbatim"; RYA-592 (Mg I) copied it again. RYA-592 then found two real defects in
its copy and fixed them THERE — which left the other two copies carrying the defect.
The RYA-643 lineage audit confirmed the executable logic of `broaden`,
`local_renorm` and `fit_line` was byte-identical across all three (the only
differences were comments, a dead assignment and statement splitting).

Three copies of one algorithm is the reason a fix could land in one place and be
missed in two. So the corrected machinery lives HERE, once, and every harness
imports it. Do not re-copy these functions into a harness.

The two defects (found RYA-592, fixed here for everyone)
--------------------------------------------------------
(a) REST FRAME. The fit compares REST-frame synthetic flux to the observed
    profile, so a residual velocity misaligns the cores, inflates chi2 and biases
    the broadening. The solar arms are NOT at rest: HARPS carries ~+0.76 km/s and
    the IAG FTS atlas ~+0.28 km/s (measured, see `measure_arm_rv`). The old
    fitter had no velocity handling at all.

    `dv` is therefore fitted as a NUISANCE parameter over `DV_GRID`, and callers
    are expected to cross-check it against `measure_arm_rv()` — an
    abundance-BLIND measurement from clean-line core centroids. That cross-check
    is the whole guard: it is what distinguishes "the spectrum has a frame
    offset" from "the fit is absorbing a profile mismatch into a shift".

    NOTE the velocity is MEASURED per arm and FITTED per line. There is no stored
    velocity constant anywhere, and none should be added — a hardcoded frame
    correction is exactly the silent-fallback this module exists to prevent.

(b) BROADENING GRID RAILING. The old grid was `np.arange(1.5, 7.0, 0.5)` and the
    fit railed at its 1.5 km/s floor, which silently trades broadening against
    abundance. Widened to `GSIG_GRID`, and `gsig_railed` is REPORTED so a railed
    nuisance parameter can never pass unnoticed again.

Scale of the effect: for the clean optical Mg I 5528/5711 both fixes moved A by
<0.005 dex. That is a result for those lines, NOT a licence to assume the same
elsewhere — RYA-643 exists to measure it for the blue/near-UV channels.
"""
import numpy as np

# Standalone-script bootstrap: these run under foreign interpreters (e.g. venv_pysme)
# and from arbitrary cwds, so put the REPO ROOT on sys.path BEFORE importing anything
# from `pipeline`. Derived from __file__, never from cwd. (RYA-313)
import os as _os_boot, sys as _sys_boot
_sys_boot.path.insert(0, _os_boot.path.dirname(_os_boot.path.dirname(
    _os_boot.path.abspath(__file__))))
from pipeline._numcompat import trapezoid as _trapezoid  # numpy>=2 removed np.trapz (RYA-313)
CLIGHT = 299792.458

# Nuisance-parameter grids. DV brackets the measured solar-arm offsets with room
# either side; GSIG spans plausible vmac+instrumental without a floor to rail on.
DV_GRID = np.round(np.arange(-1.5, 2.55, 0.10), 2)      # km/s
GSIG_GRID = np.round(np.arange(0.4, 8.01, 0.20), 2)     # km/s (vmac + instrumental)

# Clean, unblended solar lines used ONLY to measure an arm's residual velocity.
# Abundance-blind by construction: nothing but the core centroid is read. Air
# wavelengths. (Optical; the HARPS arm covers 3782-6910 A so they are available
# for the blue channels too, where they measure the arm's GLOBAL frame offset.)
RV_CHECK_LINES = [5522.446, 5525.544, 5615.644, 5679.023, 5686.530,
                  5701.104, 5701.544, 5709.378, 5731.762]


def rot_kernel(dv, vsini, eps=0.6):
    """Rotational broadening kernel on a velocity axis (linear limb darkening)."""
    x = dv / vsini
    m = np.abs(x) < 1.0
    k = np.zeros_like(dv)
    c1 = 2 * (1 - eps) / (np.pi * vsini * (1 - eps / 3.0))
    c2 = 0.5 * eps / (vsini * (1 - eps / 3.0))
    k[m] = c1 * np.sqrt(1 - x[m] ** 2) + c2 * (1 - x[m] ** 2)
    return k


def broaden(wl, fl, vsini, gsig_kms):
    """Rotational (vsini) then Gaussian (gsig_kms, absorbing vmac+instrumental)."""
    step = wl[1] - wl[0]
    cen = np.median(wl)
    out = fl
    if vsini > 0.1:
        kv = np.arange(-vsini * 1.2, vsini * 1.2 + step, step) / cen * CLIGHT
        rk = rot_kernel(kv, vsini)
        if rk.sum() > 0:
            rk /= rk.sum()
            out = np.convolve(1 - out, rk, mode='same')
            out = 1 - out
    if gsig_kms > 0.05:
        gsig_A = gsig_kms / CLIGHT * cen
        n = int(np.ceil(4 * gsig_A / step))
        gx = np.arange(-n, n + 1) * step
        gk = np.exp(-0.5 * (gx / gsig_A) ** 2)
        gk /= gk.sum()
        out = np.convolve(1 - out, gk, mode='same')
        out = 1 - out
    return out


def local_renorm(w, f, center, hw):
    """Divide observed flux by a linear continuum fit to the top-percentile points
    in the [center-hw-1.2, center+hw+1.2] window edges."""
    win = (w > center - hw - 1.2) & (w < center + hw + 1.2)
    ww, ff = w[win], f[win]
    if len(ww) < 20:
        return w, f, win
    edge = (ww < center - hw * 0.7) | (ww > center + hw * 0.7)
    xe, ye = ww[edge], ff[edge]
    if len(xe) < 6:
        return ww, ff, win
    thr = np.percentile(ye, 70)
    keep = ye >= thr
    c = np.polyfit(xe[keep], ye[keep], 1)
    cont = np.polyval(c, ww)
    return ww, ff / np.clip(cont, 1e-3, None), win


def measure_arm_rv(obs_w, obs_f, lines=None):
    """Residual velocity of an arm (km/s) from parabolic core centroids of clean,
    unblended solar lines. Returns (median_v, n_lines_used, scatter).

    Abundance-BLIND: this is the independent check that a fitted nuisance `dv` is a
    real wavelength offset rather than the fit absorbing a profile mismatch. Never
    feed an abundance-dependent quantity in here."""
    vs = []
    for lam0 in (RV_CHECK_LINES if lines is None else lines):
        m = (obs_w > lam0 - 0.12) & (obs_w < lam0 + 0.12)
        if m.sum() < 8:
            continue
        x, y = obs_w[m], obs_f[m]
        k = int(np.argmin(y))
        sl = slice(max(k - 3, 0), min(k + 4, len(x)))
        if sl.stop - sl.start < 4:
            continue
        c = np.polyfit(x[sl], y[sl], 2)
        if c[0] <= 0:                     # not a minimum — reject, never guess
            continue
        vs.append((-c[1] / (2 * c[0]) - lam0) / lam0 * CLIGHT)
    if not vs:
        return None, 0, 0.0
    return float(np.median(vs)), len(vs), float(np.std(vs))


def require_arm_rv(obs_w, obs_f, arm_label, min_lines=4):
    """Measure an arm's residual velocity and LOUD-FAIL if it cannot be sourced.

    The frame correction must always come from a measurement on the arm being fitted —
    never a hardcoded velocity, never a silent zero. If the clean check lines are not
    covered or are unusable, that is a real problem with the arm and the run must stop
    rather than quietly fit `dv` against nothing. Returns (v_kms, n_lines, scatter)."""
    v, n, sd = measure_arm_rv(obs_w, obs_f)
    if v is None or n < min_lines:
        raise SystemExit(
            f"RYA-643 FRAME NOT SOURCED [{arm_label}]: residual velocity measurable on only "
            f"{n} of {len(RV_CHECK_LINES)} clean check lines (need >= {min_lines}). The rest-frame "
            f"correction has no measured source for this arm. Refusing to fit dv unanchored or "
            f"to assume zero — fix the arm's coverage/normalisation, do not hardcode a velocity.")
    return v, n, sd


def fit_profile(center, obs_w, obs_f, synth, hw, vsini, a_lo, a_hi,
                core_hw=0.4, dv_grid=None, gsig_grid=None):
    """chi2 in-window profile fit of A over the fit window — NOT an isolated-line
    EW inversion. `synth` = {A: (wl, flux)} rest-frame syntheses.

    Free parameters: A (the measurement) + two nuisances, `gsig` (vmac +
    instrumental) and `dv` (residual velocity — see the module docstring).

    Returns None if the window has too little coverage, else a dict with:
      A, chi2/red_chi2, gsig, gsig_railed, dv, npix, obs_ew_mA, dEW_dA, railed.
    """
    dv_grid = DV_GRID if dv_grid is None else dv_grid
    gsig_grid = GSIG_GRID if gsig_grid is None else gsig_grid
    ww, ff, _ = local_renorm(obs_w, obs_f, center, hw)
    pad = (ww > center - hw - 0.05) & (ww < center + hw + 0.05)
    if pad.sum() < 10:
        return None
    wp, fp = ww[pad], ff[pad]

    best = dict(chi2=1e30)
    for a, (sw, sf) in synth.items():
        for gs in gsig_grid:
            sb = broaden(sw, sf, vsini, gs)
            for dv in dv_grid:
                xo = wp / (1 + dv / CLIGHT)        # shift the OBSERVED to rest
                sel = (xo > center - hw) & (xo < center + hw)
                if sel.sum() < 10:
                    continue
                r = fp[sel] - np.interp(xo[sel], sw, sb)
                chi2 = float(np.sum(r ** 2)) / max(int(sel.sum()) - 3, 1) / (0.01 ** 2)
                if chi2 < best['chi2']:
                    best = dict(chi2=chi2, A=float(a), gsig=float(gs), dv=float(dv),
                                npix=int(sel.sum()))
    if 'A' not in best:
        return None

    # parabolic refine in A at the best (gsig, dv)
    gs, dv = best['gsig'], best['dv']
    xo = wp / (1 + dv / CLIGHT)
    sel = (xo > center - hw) & (xo < center + hw)
    xo, yo = xo[sel], fp[sel]
    As = np.array(sorted(synth))
    chis = np.array([np.sum((yo - np.interp(xo, synth[a][0],
                                            broaden(*synth[a], vsini, gs))) ** 2)
                     for a in As])
    k = int(np.argmin(chis))
    A_ref = float(As[k])
    if 0 < k < len(As) - 1:
        d = chis[k + 1] - 2 * chis[k] + chis[k - 1]
        if d > 0:
            A_ref = float(As[k] - 0.5 * (chis[k + 1] - chis[k - 1]) / d * (As[1] - As[0]))
    best['A'] = A_ref
    best['red_chi2'] = best['chi2']
    best['gsig_railed'] = bool(gs <= gsig_grid[0] + 1e-9 or gs >= gsig_grid[-1] - 1e-9)
    best['obs_ew_mA'] = float(_trapezoid(1 - yo, xo) * 1000.0)

    # sensitivity: change in synthetic CORE EW (mA) per dex of A near the fit.
    # A weak / blend-dominated line barely responds -> low sensitivity -> unreliable.
    def _core_ew(a):
        aa = min(a_hi, max(a_lo, a))
        kk = int(np.argmin(np.abs(As - aa)))
        sw, sf = synth[float(As[kk])]
        sb = broaden(sw, sf, vsini, gs)
        m = (sw > center - core_hw) & (sw < center + core_hw)
        return float(_trapezoid(1 - sb[m], sw[m]) * 1000.0)
    best['dEW_dA'] = round(abs(_core_ew(A_ref + 0.15) - _core_ew(A_ref - 0.15)) / 0.30, 1)
    best['core_EW_mA'] = round(_core_ew(A_ref), 2)   # synthetic core EW at the fitted A
    best['railed'] = bool(A_ref <= a_lo + 0.03 or A_ref >= a_hi - 0.03)
    return best
