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

THE RELIABILITY RULE (RYA-679) — one definition, see `assess_reliability`
-------------------------------------------------------------------------
The reliability constants used to be redefined per harness, with three different
red_chi2 ceilings live at once (60.0 in RYA-564/581, 15.0 in RYA-565, 5.0 in
RYA-560's deblend path) and no ceiling at all in RYA-551/560-plain/592. RYA-679
adjudicated that spread; the ratified rule and its evidence live in
`assess_reliability` below. Harnesses import it. Do not redefine these constants.

TWO FIT ENTRYPOINTS — deliberate variants, NOT duplicates (RYA-679 §3D)
-----------------------------------------------------------------------
`fit_profile` and `fit_profile_deblend` are both first-class and both belong here.
They are not two copies of one algorithm (the situation this module was created to
end); they are two chi2 DOMAINS over the same measurement, and the choice between
them is a modelling decision a harness makes explicitly:

  * `fit_profile`        — chi2 over the whole fit window, observed renormalised by
                           `local_renorm`. Right for a clean, uncrowded window.
  * `fit_profile_deblend` — chi2 over the target's own pixels only, with a continuum
                           RATIO fitted on the blend pixels. Right for a crowded
                           window, where the full-window statistic stops being about
                           the element (see `assess_reliability` for the numbers).

Their `red_chi2` values are therefore NOT comparable with each other, and a single
threshold cannot be applied to both. That is the core RYA-679 finding.
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

# The per-pixel flux uncertainty ASSUMED by every chi2 in this module. It is a fixed
# constant, NOT a measured noise estimate — see `assess_reliability` for what that
# means for the interpretation of `red_chi2`. Named here so the two chi2 expressions
# below cannot drift apart, and so the assumption is visible rather than a literal.
SIGMA_FLUX_ASSUMED = 0.01

# ── Reliability rule (RYA-679 ratified) ──────────────────────────────────────────
# THE single definition. No harness may redefine these.
RELIABLE_DEWDA = 40.0      # mA/dex core-EW sensitivity floor (RYA-551 origin)
RCHI2_REVIEW = 5.0         # REPORTING trigger only — never gates `reliable`


def assess_reliability(fit, dewda_floor=None):
    """THE reliability rule for an in-window profile fit (RYA-679, ratified).

        reliable = (not railed) AND dEW_dA >= RELIABLE_DEWDA

    `red_chi2` is REPORTED and REVIEW-FLAGGED, but does NOT gate. Returns a dict of
    {reliable, rchi2_review, rchi2_review_reason} to merge into a per-line record.

    Why no red_chi2 ceiling — three findings, all measured
    ------------------------------------------------------
    (1) THE STATED RATIONALE FOR 60.0 IS BACKWARDS. RYA-564 set a ceiling of 60 with
        the reason "sigma_flux=0.01 floor inflates rchi2". Measured per-pixel noise in
        the actual fit windows of the actual arms (MAD of 2nd differences, which is
        blind to line/continuum structure) is sigma_pix = 0.00007-0.0051 — so the
        assumed 0.01 is 2x to 146x LARGER than the truth, everywhere. A sigma larger
        than the truth DEFLATES chi2. The floor suppresses red_chi2; it cannot inflate
        it. A perfect model on pure photon noise would score red_chi2 = 0.0000-0.26
        here, not 1. So 60.0's only written justification is wrong in sign.

    (2) red_chi2 IS NOT A CHI2 — it is a rescaled residual RMS,
        red_chi2 = (RMS_resid / SIGMA_FLUX_ASSUMED)^2. Because the assumed sigma
        swamps the real noise, anything above ~0.26 is essentially pure model and
        continuum systematic with no photon-noise content. Sr II 4077's red_chi2
        78.27 means a residual RMS of 8.9% of the continuum against a 0.42% photon
        noise — 21x the noise, i.e. ~99.8% systematic. There is no statistical
        calibration here to hang a pass/fail bar on.

    (3) FULL-WINDOW red_chi2 MEASURES THE BLEND LIST, NOT THE ELEMENT. Controlled
        experiment (RYA-560/585 Zr II, identical lines, identical spectra, only the
        chi2 domain and continuum treatment changed):

            line        fit_profile   fit_profile_deblend   ratio
            4208.98        83.12             0.39           213x
            4258.04        25.92             0.35            74x
            4442.99        15.93             1.49            11x

        while dEW_dA barely moved (36.8->36.2, 29.6->33.5, 29.7->33.5). Across
        species the full-window statistic simply tracks how crowded the window is —
        clean red windows score Eu II 6645: 0.16, Co I 5352: 0.20, Ba II 5853: 0.71;
        crowded blue/near-UV score Zr II 4208: 83.12, Sr II 4077: 78.27, Sr II 4215:
        179.97. Gating on it would gate on which part of the spectrum a line lives
        in, which is not a statement about whether A(X) is measurable.

    And it never bound anything anyway. In RYA-564's OWN data — the ticket that
    introduced 60.0 — every reliable Co line scores 0.04-4.72, and every line above
    60 is already excluded by `railed`. The term has never changed a Co disposition,
    nor Ba's (0.71), nor Eu's (fails on dEW_dA 13.9), nor Zr's (fails on dEW_dA 33.5).
    The ONLY species on which any candidate ceiling was load-bearing is Sr II — so a
    ceiling would have functioned, in practice, purely as a silent veto on the one
    live adoption candidate, on the strength of a near-UV blend model. That is the
    wrong reason to demote a measurement.

    What replaces it: RCHI2_REVIEW = 5.0 raises a LOUD, non-gating flag. It is set at
    the empirical upper edge of the clean-window population (worst reliable Co line
    4.72; Ba <=0.83; Eu 0.16; Zr deblended <=1.66) and means "the in-window model does
    not reproduce this window — inspect the blend list and continuum before adopting
    this value", which is precisely and only what a high red_chi2 licenses. Being
    non-gating, its value cannot flip any disposition by construction, so it is not a
    tunable knob (validate-don't-tune).

    NOTE the review flag is only comparable WITHIN a fit entrypoint — see the module
    docstring. A `fit_profile_deblend` red_chi2 above 5 is a much stronger signal than
    a `fit_profile` one, because the blend pixels have already been divided out.
    """
    floor = RELIABLE_DEWDA if dewda_floor is None else dewda_floor
    dewda = fit.get('dEW_dA')
    railed = bool(fit.get('railed'))
    rchi2 = fit.get('red_chi2')
    reliable = bool((not railed) and dewda is not None and dewda >= floor)
    review = bool(rchi2 is not None and rchi2 > RCHI2_REVIEW)
    reason = None
    if review:
        rms = (rchi2 ** 0.5) * SIGMA_FLUX_ASSUMED
        reason = (f"red_chi2 {rchi2:.2f} > {RCHI2_REVIEW} (residual RMS {rms * 100:.1f}% "
                  f"of continuum): the in-window model does not reproduce this window. "
                  f"Does NOT gate `reliable` (RYA-679) — inspect the blend list and "
                  f"continuum before adopting this value.")
    return dict(reliable=reliable, rchi2_review=review, rchi2_review_reason=reason)


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
                chi2 = (float(np.sum(r ** 2)) / max(int(sel.sum()) - 3, 1)
                        / (SIGMA_FLUX_ASSUMED ** 2))
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


def _cont_ratio(xc, ratio, keep0, n_clip=2):
    """Robust linear fit of the obs/synth flux ratio over a pixel mask (RYA-585).

    Closed-form least squares + `n_clip` rounds of 2-sigma clipping. `keep0` is the
    pixel mask the fit is allowed to use — callers pass the BLEND (target-free)
    pixels, which is what makes this a continuum estimate rather than a fit to the
    line being measured. Returns (slope, intercept, residual_scatter, n_used).
    """
    keep = keep0.copy()
    c1 = c0 = 0.0
    sd = 0.0
    for _ in range(n_clip + 1):
        n = int(keep.sum())
        if n < 5:                                   # too few blend pixels to fit a slope
            return 0.0, float(np.median(ratio[keep0])) if keep0.any() else 1.0, 0.0, n
        x, y = xc[keep], ratio[keep]
        sx, sy = x.sum(), y.sum()
        sxx, sxy = (x * x).sum(), (x * y).sum()
        den = n * sxx - sx * sx
        if abs(den) < 1e-30:
            c1, c0 = 0.0, sy / n
        else:
            c1 = (n * sxy - sx * sy) / den
            c0 = (sy - c1 * sx) / n
        res = ratio - (c0 + c1 * xc)
        sd = float(np.std(res[keep]))
        if sd <= 0:
            break
        keep = keep0 & (np.abs(res) < 2.0 * sd)
    return float(c1), float(c0), sd, int(keep.sum())


def fit_profile_deblend(center, obs_w, obs_f, synth, blend_only, hw, vsini,
                        a_lo, a_hi, core_hw=0.4, sens_frac=0.05):
    """In-window blend-fit variant of `fit_profile` (RYA-585, the RYA-551 pattern).

    Same measurement (A over rest-frame syntheses, `gsig`/`dv` as nuisances), same
    `dEW_dA` definition VERBATIM — the reliability floor is NOT moved by using this
    path. Two things change, both aimed at a blend/continuum systematic:

    (1) CONTINUUM. `fit_profile` calls `local_renorm`, which normalises the OBSERVED
        spectrum by a top-percentile fit over the window while the SYNTHETIC keeps
        its true continuum. That operator is asymmetric: in a crowded blue window
        there are no true continuum pixels, so the "continuum" is pinned to blend
        shoulders and the observed profile is scaled against a synthetic that was
        never scaled the same way. Here we instead fit a low-order continuum RATIO
        obs/synth over the BLEND pixels only — pixels where the target line
        contributes nothing — so the identical normalisation applies to both by
        construction, and the continuum can never eat the target line.

    (2) CHI2 DOMAIN. chi2 is restricted to the pixels the target actually touches
        (|target contribution| > `sens_frac` of its peak). The full-window red_chi2
        of `fit_profile` is dominated by VALD's rendition of the couple-hundred
        neighbouring components; it is a statistic about the blend list, not about
        the element being measured, and it is why RYA-560 saw red_chi2 41-91 on
        lines whose own profiles fit well.

    `blend_only` is a (wl, flux) synthesis of the SAME window with the target
    element suppressed; `blend_only - synth[A]` is therefore the target's own
    contribution, which defines both masks above.

    Returns None on insufficient coverage, else the `fit_profile` keys plus:
      cont0/cont1 (continuum ratio), cont_scatter, cont_npix, npix_window,
      target_EW_mA, sat_index, target_core_depth_frac.
    """
    As = np.array(sorted(synth))
    a_mid = float(As[len(As) // 2])
    sw = synth[a_mid][0]
    pad = (obs_w > center - hw - 0.05) & (obs_w < center + hw + 0.05)
    if pad.sum() < 10:
        return None
    wp, fp = obs_w[pad], obs_f[pad]
    # target-only depth in the rest frame (blend-only MINUS blend+target)
    contrib_raw = blend_only[1] - synth[a_mid][1]

    best = dict(chi2=1e30)
    cache = {}
    for gs in GSIG_GRID:
        sb_all = np.array([broaden(sw, synth[float(a)][1], vsini, gs) for a in As])
        zc = 1.0 - broaden(sw, 1.0 - contrib_raw, vsini, gs)   # broadened target depth
        cache[float(gs)] = (sb_all, zc)
        for dv in DV_GRID:
            xo = wp / (1 + dv / CLIGHT)            # shift the OBSERVED to rest
            sel = (xo > center - hw) & (xo < center + hw)
            if sel.sum() < 10:
                continue
            xs, ys = xo[sel], fp[sel]
            zci = np.interp(xs, sw, zc)
            m = np.abs(zci) > sens_frac * np.abs(zci).max()    # target pixels
            if m.sum() < 5 or (~m).sum() < 10:                 # need both populations
                continue
            xc = xs - center
            for i, a in enumerate(As):
                si = np.interp(xs, sw, sb_all[i])
                c1, c0, _, _ = _cont_ratio(xc, ys / np.clip(si, 1e-3, None), ~m)
                r = ys / np.clip(c0 + c1 * xc, 1e-3, None) - si
                chi2 = (float(np.sum(r[m] ** 2)) / max(int(m.sum()) - 3, 1)
                        / (SIGMA_FLUX_ASSUMED ** 2))
                if chi2 < best['chi2']:
                    best = dict(chi2=chi2, A=float(a), gsig=float(gs), dv=float(dv),
                                npix=int(m.sum()), npix_window=int(sel.sum()),
                                cont0=c0, cont1=c1)
    if 'A' not in best:
        return None

    # parabolic refine in A at the best (gsig, dv), continuum re-fitted per trial A
    gs, dv = best['gsig'], best['dv']
    sb_all, zc = cache[gs]
    xo = wp / (1 + dv / CLIGHT)
    sel = (xo > center - hw) & (xo < center + hw)
    xs, ys = xo[sel], fp[sel]
    zci = np.interp(xs, sw, zc)
    m = np.abs(zci) > sens_frac * np.abs(zci).max()
    xc = xs - center
    chis, sd, nu = [], 0.0, 0
    for i, a in enumerate(As):
        si = np.interp(xs, sw, sb_all[i])
        c1, c0, sd, nu = _cont_ratio(xc, ys / np.clip(si, 1e-3, None), ~m)
        chis.append(float(np.sum((ys / np.clip(c0 + c1 * xc, 1e-3, None) - si)[m] ** 2)))
    chis = np.array(chis)
    k = int(np.argmin(chis))
    A_ref = float(As[k])
    if 0 < k < len(As) - 1:
        d = chis[k + 1] - 2 * chis[k] + chis[k - 1]
        if d > 0:
            A_ref = float(As[k] - 0.5 * (chis[k + 1] - chis[k - 1]) / d * (As[1] - As[0]))
    best.update(A=A_ref, red_chi2=best['chi2'],
                gsig_railed=bool(gs <= GSIG_GRID[0] + 1e-9 or gs >= GSIG_GRID[-1] - 1e-9),
                cont_scatter=sd, cont_npix=nu)

    # sensitivity: IDENTICAL definition to fit_profile (do not diverge — the
    # reliability floor is defined against this number).
    def _core_ew(a):
        aa = min(a_hi, max(a_lo, a))
        kk = int(np.argmin(np.abs(As - aa)))
        mm = (sw > center - core_hw) & (sw < center + core_hw)
        return float(_trapezoid(1 - sb_all[kk][mm], sw[mm]) * 1000.0)
    best['dEW_dA'] = round(abs(_core_ew(A_ref + 0.15) - _core_ew(A_ref - 0.15)) / 0.30, 1)
    best['core_EW_mA'] = round(_core_ew(A_ref), 2)
    best['railed'] = bool(A_ref <= a_lo + 0.03 or A_ref >= a_hi - 0.03)

    # Diagnostics that explain a LOW dEW_dA — i.e. distinguish "the blend model was
    # wrong" (fixable here) from "the line is saturated / blend-dominated" (not).
    mm = (sw > center - core_hw) & (sw < center + core_hw)
    ew_t = float(_trapezoid(zc[mm], sw[mm]) * 1000.0)         # the TARGET's own core EW
    best['target_EW_mA'] = round(ew_t, 2)
    # sat_index = dEW/dA over the optically-thin (linear-COG) expectation ln(10)*EW.
    # ->1 is the unsaturated linear regime; <<1 means the core is on the flat part.
    best['sat_index'] = round(best['dEW_dA'] / max(np.log(10) * ew_t, 1e-9), 3)
    tot = float(np.max(1 - sb_all[int(np.argmin(np.abs(As - A_ref)))][mm]))
    # fraction of the observed core depth that is the TARGET rather than blends
    best['target_core_depth_frac'] = round(float(np.max(zc[mm])) / max(tot, 1e-9), 3)
    return best
