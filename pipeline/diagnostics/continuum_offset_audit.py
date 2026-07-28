#!/usr/bin/env python3
"""
Cross-window continuum-offset audit + gated [O I] renorm-refit proof (RYA-451 step 1).

RYA-449: the [O I] 6300 O excess is a continuum-placement artifact (obs continuum
~0.7% below synth). But "obs normalized low" (NORMALIZATION) and "synth missing CN
micro-opacity" (MICRO-OPACITY) are degenerate on one window and need OPPOSITE fixes.

PART A (audit): obs-vs-synth continuum offset across many VIS windows (405-665 nm).
Uniform offset (incl. molecule-poor windows) -> NORMALIZATION (systemic). Offset that
grows with our synth's own molecular density -> MICRO-OPACITY (line-list gap; renorm
masks it). AMENDMENT 1 also classifies the offset's WAVELENGTH SHAPE (part_a_trend):
CONSTANT / SMOOTH_DRIFT (genuine continuum shape, survives density control) /
OPACITY_DRIFT (lambda trend vanishes once density controlled) / LOCALIZED -- the
450-vs-650 question. If the normalization-vs-opacity verdict and the shape disagree,
the gate surfaces the conflict and WITHHOLDS Part B.

PART B (GATED on PART A = NORMALIZATION, no shape conflict): prove the fix -- locally
renormalize the [O I] window by the measured obs/synth continuum ratio, refit A(O),
confirm it lands ~8.69-8.74 (solar O becomes a MEASUREMENT). Withheld otherwise.

AMENDMENT 2: VIS / UV / IR are THREE distinct continuum regimes, not two. UV (e.g.
UVES-346 blue) has essentially no true continuum, so the synth-defined-continuum-pixel
estimator is MANDATORY there. The production fix must be PER-ORDER (over many pixels),
not the noisy narrow-window ratio -- [O I] 6300 amplifies continuum error (~0.26 dex
per 1%), which is why Part B's flat renorm overshot. See _print_regime_note.

Validate-don't-tune: the renorm factor is the measured continuum ratio, not a value
chosen to hit 8.69. All metrics ride the noise-robust estimator (median at synth-
defined continuum pixels), NOT the noise-biased 0.97 quantile.

RYA-452 (cross-window continuum-offset audit; RYA-451 step 1; AMENDMENTs 1 + 2).
"""
from __future__ import annotations
import argparse
import os
import numpy as np

import pipeline.cno_synthesis as cs
from pipeline.cno_synthesis import (
    REGIONS, VIS_DIAGNOSTICS, get_star_params,
    _load_atmosphere, _load_synth_resources, _load_observed_spectrum,
    _atom_codes, _fixed_ab, _synth_window, _fit_element, preflight,
)
import ispec

SOLAR_COMP = dict(C=8.491, N=7.558, O=8.69, Ni=6.20)   # converged solar (RYA-448); rest from sab
# RYA-452 AMENDMENT 1: span the full VIS arm (405-665 nm) so the 450-vs-650 offset
# TREND is actually measurable (filtered to obs coverage at runtime).
AUDIT_WINDOWS_NM = [(c - 0.8, c + 0.8) for c in
                    (405, 425, 445, 465, 485, 505, 525, 545, 565, 585,
                     605, 625, 645, 665)]
OI_FIT_WINDOW_A = (6299.5, 6301.0)
CONT_Q    = 0.70          # synth pixels above this quantile = (noiseless) continuum pixels
COARSE_NM = 0.001         # 0.01 A synth step (continuum needs no finer)
OFFSET_FLOOR = 0.003      # |offset| above this in LOW-mol windows => normalization present
OFFSET_OPAC  = 0.005      # offset only in HIGH-mol windows => micro-opacity
CONTAM_OFFSET = 0.03      # |offset| above this => NO true continuum (blue line-forest / strong-
                          # line / edge; max obs flux well below 1.0) -- AMENDMENT 2 regime, not a
                          # ~0.7% normalization scale; excluded from the wavelength-trend fit.
SUCCESS_BAND = (8.64, 8.76)   # Part B: O consistent with Asplund 8.69 / Caffau 8.73


def _robust_cont(ow, of, sw, sf):
    """Noise-robust (obs_c, syn_c) on equal footing.

    RYA-452: a bare 0.97 quantile of the OBSERVED flux is robust to weak lines but
    NOT to noise -- HARPS obs has sigma ~0.005-0.014, so the 97th percentile of a
    clean window overshoots the true continuum by ~1.9 sigma (~+0.01) while the
    NOISELESS synth's quantile does not. That positive bias masks a real negative
    offset and falsely reads AMBIGUOUS. Instead: pick continuum pixels from the
    NOISELESS synth (top CONT_Q), then take the MEDIAN of obs and of synth at those
    same pixels (noise averages out symmetrically). obs_c - syn_c is then a fair,
    noise-robust continuum offset. (Verified: overall ~ -0.0068, matches RYA-449's
    -0.0071 on the [O I] window.)
    """
    sf_i = np.interp(ow, sw, sf)
    cont = sf_i >= np.quantile(sf, CONT_Q)
    if cont.sum() < 5:
        return np.nan, np.nan
    return float(np.median(of[cont])), float(np.median(sf_i[cont]))


def _setup(star_id, region_name):
    region = REGIONS[region_name]
    rec = get_star_params(star_id)
    params = {'teff_K': float(rec['teff']), 'logg': float(rec['logg']),
              'feh': float(rec['feh_ref']), 'vturb_kms': float(rec.get('xi', 1.0))}
    broadening = preflight(region, star_id, VIS_DIAGNOSTICS)
    atm = _load_atmosphere(params['teff_K'], params['logg'], params['feh'], params['vturb_kms'])
    ll, iso, chem = _load_synth_resources()
    sab = ispec.read_solar_abundances(cs._ISPEC_SOLAR_ABUND_FILE)
    obs_w, obs_f = _load_observed_spectrum(star_id)
    codes = _atom_codes(('C', 'N', 'O', 'Ni'), chem, sab)
    return dict(params=params, broadening=broadening, atm=atm, ll=ll, iso=iso, sab=sab,
                obs_w=obs_w, obs_f=obs_f, codes=codes)


def _synth(ctx, lo_nm, hi_nm, use_mol, tmp_dir):
    sw = np.arange(lo_nm, hi_nm + COARSE_NM * 0.5, COARSE_NM)
    fa = _fixed_ab(dict(SOLAR_COMP), ctx['codes'])
    sf = _synth_window(sw, ctx['atm'], ctx['params'], ctx['ll'], ctx['iso'], ctx['sab'],
                       fa, ctx['broadening'], use_mol, tmp_dir)
    return sf


def part_a_trend(centers, offs, mols):
    """Classify the offset's wavelength behaviour, and separate a genuine continuum-shape
    drift from an opacity-density artifact via a 2-var regression offset ~ lambda + mol_dens."""
    c = np.asarray(centers, float); y = np.asarray(offs, float); md = np.asarray(mols, float)
    span = c.max() - c.min()
    b1, a1 = np.polyfit(c, y, 1)                       # 1-var: offset(lambda)
    trend_amp = abs(b1) * span
    scatter = float(np.std(y - (a1 + b1 * c)))
    off_450, off_650 = a1 + b1 * 450.0, a1 + b1 * 650.0
    # 2-var: does lambda survive once molecular/line density is controlled?
    coef, *_ = np.linalg.lstsq(np.column_stack([np.ones_like(c), c, md]), y, rcond=None)
    lam_amp_ctrl = abs(coef[1]) * span

    print("\n-- continuum-shape (offset vs wavelength) ----------------------")
    print(f"  linear slope         = {b1*100:+.5f} /100nm   amplitude across span = {trend_amp:.4f}")
    print(f"  residual scatter     = {scatter:.4f}")
    print(f"  fitted offset @450nm = {off_450:+.4f}   @650nm = {off_650:+.4f}   (the 450-vs-650 q)")
    print(f"  lambda amplitude, density-controlled (offset~lambda+mol_dens) = {lam_amp_ctrl:.4f}")

    DRIFT_AMP, SCAT = 0.002, 0.002
    if trend_amp < DRIFT_AMP and scatter < SCAT:
        print("  SHAPE = CONSTANT: offset flat vs wavelength -> a single global factor suffices.")
        return 'CONSTANT'
    if trend_amp >= DRIFT_AMP and scatter < trend_amp:
        if lam_amp_ctrl >= DRIFT_AMP:
            print("  SHAPE = SMOOTH_DRIFT (genuine continuum shape): a monotonic offset SURVIVES "
                  "controlling for molecular/line density -> VIS normalization must be wavelength-"
                  "adaptive. Fix = bounded per-order continuum (low-order spline/poly on high-"
                  "percentile points), NOT a constant, NOT a wiggly free fit.")
            return 'SMOOTH_DRIFT'
        print("  SHAPE = OPACITY-DRIVEN DRIFT: the lambda trend VANISHES once density is controlled "
              "-> the blue-ward drift is accumulating opacity, not a continuum-shape error. "
              "Fix = line-list completeness, NOT renormalization.")
        return 'OPACITY_DRIFT'
    print(f"  SHAPE = LOCALIZED: window scatter ({scatter:.4f}) dominates any trend ({trend_amp:.4f}) "
          "-> per-window local renorm; an arm-level model won't capture it.")
    return 'LOCALIZED'


def part_a(ctx, tmp_dir):
    print("== PART A: cross-window obs-vs-synth continuum-offset audit ==")
    print("   (obs_c/syn_c = noise-robust median at synth-defined continuum pixels; RYA-452)")
    print(f"{'win_nm':>13s} {'obs_c':>8s} {'syn_c':>8s} {'offset':>8s} {'mol_dens':>9s}")
    offs, mols, cens, obs_cs, syn_cs, maxfs = [], [], [], [], [], []
    for lo, hi in AUDIT_WINDOWS_NM:
        m = (ctx['obs_w'] >= lo) & (ctx['obs_w'] <= hi)
        if m.sum() < 20:
            continue
        sw = np.arange(lo, hi + COARSE_NM * 0.5, COARSE_NM)
        sf_on  = _synth(ctx, lo, hi, True,  tmp_dir)
        sf_off = _synth(ctx, lo, hi, False, tmp_dir)
        obs_c, syn_c = _robust_cont(ctx['obs_w'][m], ctx['obs_f'][m], sw, sf_on)
        if not (np.isfinite(obs_c) and np.isfinite(syn_c)):
            print(f"{lo:6.1f}-{hi:5.1f}   skipped (no synth-defined continuum pixels)")
            continue
        mol_dens = float(np.mean(sf_off - sf_on))      # >0 : molecules our synth DOES include
        offset = obs_c - syn_c
        offs.append(offset); mols.append(mol_dens); cens.append(0.5 * (lo + hi))
        obs_cs.append(obs_c); syn_cs.append(syn_c); maxfs.append(float(ctx['obs_f'][m].max()))
        print(f"{lo:6.1f}-{hi:5.1f} {obs_c:8.4f} {syn_c:8.4f} {offset:+8.4f} {mol_dens:9.4f}")

    offs, mols, cens = np.array(offs), np.array(mols), np.array(cens)
    obs_cs, syn_cs, maxfs = np.array(obs_cs), np.array(syn_cs), np.array(maxfs)

    # RYA-452 standing rule (flag/exclude badly-contaminated windows) + AMENDMENT 2:
    # windows with a GROSS offset have NO true continuum -- the blue line-forest / strong-
    # line / edge regime where max obs flux never reaches 1.0. They are not a ~0.7%
    # normalization scale and they corrupt both the molecular medians and the polyfit
    # trend. Flag them and run the OPERATIVE analysis on the true-continuum subset.
    contam = np.abs(offs) > CONTAM_OFFSET
    keep = ~contam
    if contam.any():
        print("\n-- no-true-continuum windows (gross offset; AMENDMENT 2 regime; EXCLUDED) --")
        for cc, of, oc, sc, mx in zip(cens[contam], offs[contam], obs_cs[contam],
                                      syn_cs[contam], maxfs[contam]):
            print(f"  {cc:6.1f} nm: offset {of:+.4f}  obs_c {oc:.4f} vs syn_c {sc:.4f}  "
                  f"max obs flux {mx:.3f} -> no continuum pixel reaches 1.0")

    use = keep if keep.sum() >= 3 else np.ones_like(offs, bool)
    o_u, m_u, c_u = offs[use], mols[use], cens[use]
    off_lo = float(np.median(o_u[m_u <= np.quantile(m_u, 0.34)]))
    off_hi = float(np.median(o_u[m_u >= np.quantile(m_u, 0.66)]))
    print(f"\n  (medians on the {use.sum()} true-continuum windows)")
    print(f"  median offset, LOW-molecular windows  = {off_lo:+.4f}")
    print(f"  median offset, HIGH-molecular windows = {off_hi:+.4f}")
    print(f"  overall median offset                 = {np.median(o_u):+.4f}  "
          f"(RYA-449 [O I] window ~ -0.0071)")

    print("\n-- PART A verdict (normalization vs opacity) -------------------")
    if abs(off_lo) >= OFFSET_FLOOR:
        print(f"  NORMALIZATION-SCALE: offset present even in molecule-poor windows "
              f"({off_lo:+.4f}). Global normalization scale issue, NOT missing opacity -- and "
              f"SYSTEMIC (every species' continuum). Fix = re-normalization.")
        verdict = 'NORMALIZATION'
    elif abs(off_hi) >= OFFSET_OPAC:
        print(f"  MICRO-OPACITY: offset ~0 where molecules are sparse ({off_lo:+.4f}) but grows "
              f"with molecular density ({off_hi:+.4f}) -> synth MISSING molecular (CN) opacity. "
              f"Fix = line-list completeness patch; a blind renorm would MASK it.")
        verdict = 'MICRO_OPACITY'
    else:
        print(f"  AMBIGUOUS: low-mol {off_lo:+.4f}, high-mol {off_hi:+.4f}. Post numbers for Ryan.")
        verdict = 'AMBIGUOUS'

    # RYA-452 AMENDMENT 1: classify the offset's WAVELENGTH shape. The verbatim spec runs
    # on ALL windows (printed, as specified); the OPERATIVE shape (returned/gated) is
    # re-classified on the true-continuum subset so the blue no-continuum windows can't
    # fake a drift. Both reported.
    print("\n[AMENDMENT 1, verbatim spec: part_a_trend on all windows]")
    shape_all = part_a_trend(cens, offs, mols)
    if use.sum() < len(offs) and use.sum() >= 3:
        print("\n[AMENDMENT 1, operative: part_a_trend on true-continuum subset]")
        shape = part_a_trend(c_u, o_u, m_u)
    else:
        shape = shape_all
    print(f"\n  SHAPE all-windows (verbatim) = {shape_all}   |   "
          f"SHAPE true-continuum (operative) = {shape}")
    return verdict, shape


def part_b(ctx, tmp_dir):
    print("\n== PART B: [O I] window local renorm + refit (proof: solar O -> measured) ==")
    lo_nm, hi_nm = OI_FIT_WINDOW_A[0] / 10.0, OI_FIT_WINDOW_A[1] / 10.0
    m = (ctx['obs_w'] >= lo_nm) & (ctx['obs_w'] <= hi_nm)
    sw = np.arange(lo_nm, hi_nm + COARSE_NM * 0.5, COARSE_NM)
    obs_c, syn_c = _robust_cont(ctx['obs_w'][m], ctx['obs_f'][m], sw,
                                _synth(ctx, lo_nm, hi_nm, True, tmp_dir))
    corr = syn_c / obs_c
    print(f"  obs continuum {obs_c:.4f}, synth continuum {syn_c:.4f} -> renorm factor {corr:.5f}")

    obs_corr = ctx['obs_f'].copy()
    win = (ctx['obs_w'] >= lo_nm - 0.2) & (ctx['obs_w'] <= hi_nm + 0.2)
    obs_corr[win] *= corr                              # local renorm of the [O I] window only

    diag = {d.key: d for d in VIS_DIAGNOSTICS}['OI_6300']
    fit = lambda of: _fit_element(ctx['obs_w'], of, ctx['atm'], ctx['params'], 'O',
                                  dict(SOLAR_COMP), ctx['codes'], diag.windows_A,
                                  diag.use_molecules, ctx['broadening'], 7.6, 10.0,
                                  ctx['ll'], ctx['iso'], ctx['sab'], tmp_dir)['A_X']
    a0, a1 = fit(ctx['obs_f']), fit(obs_corr)
    print(f"  A(O) before renorm = {a0}   (expect ~8.80, the RYA-449 value)")
    print(f"  A(O) after  renorm = {a1}")

    print("\n-- PART B verdict ----------------------------------------------")
    if SUCCESS_BAND[0] <= a1 <= SUCCESS_BAND[1]:
        print(f"  SUCCESS: local renorm lands A(O)={a1} in {SUCCESS_BAND} -- consistent with "
              f"Asplund 8.69 / Caffau 8.73. Solar [O I] O is now a MEASUREMENT, not a cite. "
              f"Production renorm + EW-verify layer = next brief.")
    else:
        print(f"  PARTIAL: renorm moved A(O) {a0} -> {a1}, not into {SUCCESS_BAND}. Renorm is part "
              f"of it; residual remains -- post for Ryan.")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--star', default='solar')
    ap.add_argument('--region', default='vis')
    ap.add_argument('--tmp-dir', default='/tmp/ispec_cont_audit')
    ap.add_argument('--force-part-b', action='store_true', help='debug: run Part B regardless')
    a = ap.parse_args()
    os.makedirs(a.tmp_dir, exist_ok=True)
    ctx = _setup(a.star, a.region)
    verdict, shape = part_a(ctx, a.tmp_dir)

    # RYA-452 AMENDMENT 1 reconciliation: if the normalization-vs-opacity verdict and the
    # wavelength-shape sub-analysis disagree, the two are in conflict -> surface both, do
    # NOT auto-proceed to the renorm proof.
    disagree = ((verdict == 'NORMALIZATION' and shape == 'OPACITY_DRIFT') or
                (verdict == 'MICRO_OPACITY' and shape == 'SMOOTH_DRIFT'))
    print("\n== gate ========================================================")
    if disagree and not a.force_part_b:
        print(f"  RECONCILIATION CONFLICT: primary verdict = {verdict} but wavelength shape = "
              f"{shape}. The two sub-analyses disagree -> Part B WITHHELD. Post BOTH for Ryan; "
              f"do not auto-proceed.")
    elif verdict == 'NORMALIZATION' or a.force_part_b:
        if verdict == 'NORMALIZATION':
            print(f"  verdict {verdict} + shape {shape} consistent -> Part B proceeds.")
        part_b(ctx, a.tmp_dir)
    else:
        print(f"  Part B withheld (Part A verdict = {verdict}, shape = {shape}).")

    _print_regime_note(shape)


def _print_regime_note(shape):
    """RYA-452 AMENDMENT 2 + Part-B-overshoot note: per-arm continuum regimes and the
    per-order (not narrow-window) production-fix discipline. Documentation, not a fix."""
    print("\n-- per-arm continuum regimes (AMENDMENT 2) + production-fix note -")
    print("  THREE regimes, not two -- VIS / UV / IR each need their OWN continuum model:")
    print("   - VIS: real continuum points exist -> bounded per-order adaptive continuum "
          "(low-order spline/poly on high-percentile points; physics-anchored, not wiggly).")
    print("   - IR (CRIRES+): telluric forest + thermal/blaze, entangled with RYA-373 -> a "
          "telluric-aware continuum model. Does NOT inherit the VIS approach.")
    print("   - UV (e.g. UVES-346 blue ~300-390nm: NH 3360, OH ~310, CN violet 388): the "
          "HARDEST regime -- line crowding leaves essentially NO true continuum, only a "
          "pseudo-continuum. Upper-envelope/percentile breaks down (no line-free pixels), so "
          "the synth-defined-continuum-pixel estimator is MANDATORY, not optional. Lower SNR + "
          "steeper blaze; must be model-informed from the start.")
    print("  Production fix (next RYA-451 brief): [O I] 6300 is a continuum-sensitivity "
          "amplifier (~0.26 dex per 1% continuum), so Part B's flat narrow-window renorm "
          "overshot (8.801 -> 8.598). Determine the continuum PER-ORDER over many pixels "
          f"(robust), THEN apply at 6300 -- do NOT design the fix off the noisy narrow-window "
          f"ratio. A bounded per-order model unifies the drift fix and the noise fix. "
          f"(this run's VIS shape = {shape}.)")


if __name__ == '__main__':
    main()
