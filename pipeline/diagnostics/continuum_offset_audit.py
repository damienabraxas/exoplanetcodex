#!/usr/bin/env python3
"""
Cross-window continuum-offset audit + gated [O I] renorm-refit proof (RYA-451 step 1).

RYA-449: the [O I] 6300 O excess is a continuum-placement artifact (obs continuum
~0.7% below synth). But "obs normalized low" (NORMALIZATION) and "synth missing CN
micro-opacity" (MICRO-OPACITY) are degenerate on one window and need OPPOSITE fixes.

PART A (audit): obs-vs-synth upper-continuum offset across many VIS windows. Uniform
offset (incl. molecule-poor windows) -> NORMALIZATION (systemic). Offset that grows
with our synth's own molecular density -> MICRO-OPACITY (line-list gap; renorm masks it).

PART B (GATED on PART A = NORMALIZATION): prove the fix -- locally renormalize the
[O I] window by the measured obs/synth continuum ratio, refit A(O), confirm it lands
~8.69-8.74 (solar O becomes a MEASUREMENT). Withheld on MICRO-OPACITY/AMBIGUOUS.

Validate-don't-tune: the renorm factor is the measured continuum ratio, not a value
chosen to hit 8.69.

RYA-452 (cross-window continuum-offset audit; RYA-451 step 1).
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
AUDIT_WINDOWS_NM = [(c - 0.8, c + 0.8) for c in
                    (505, 525, 545, 565, 585, 600, 615, 630, 645, 660, 675)]
OI_FIT_WINDOW_A = (6299.5, 6301.0)
CONT_Q    = 0.70          # synth pixels above this quantile = (noiseless) continuum pixels
COARSE_NM = 0.001         # 0.01 A synth step (continuum needs no finer)
OFFSET_FLOOR = 0.003      # |offset| above this in LOW-mol windows => normalization present
OFFSET_OPAC  = 0.005      # offset only in HIGH-mol windows => micro-opacity
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


def part_a(ctx, tmp_dir):
    print("== PART A: cross-window obs-vs-synth continuum-offset audit ==")
    print("   (obs_c/syn_c = noise-robust median at synth-defined continuum pixels; RYA-452)")
    print(f"{'win_nm':>13s} {'obs_c':>8s} {'syn_c':>8s} {'offset':>8s} {'mol_dens':>9s}")
    offs, mols = [], []
    for lo, hi in AUDIT_WINDOWS_NM:
        m = (ctx['obs_w'] >= lo) & (ctx['obs_w'] <= hi)
        if m.sum() < 20:
            continue
        sw = np.arange(lo, hi + COARSE_NM * 0.5, COARSE_NM)
        sf_on  = _synth(ctx, lo, hi, True,  tmp_dir)
        sf_off = _synth(ctx, lo, hi, False, tmp_dir)
        obs_c, syn_c = _robust_cont(ctx['obs_w'][m], ctx['obs_f'][m], sw, sf_on)
        mol_dens = float(np.mean(sf_off - sf_on))      # >0 : molecules our synth DOES include
        offset = obs_c - syn_c
        offs.append(offset); mols.append(mol_dens)
        print(f"{lo:6.1f}-{hi:5.1f} {obs_c:8.4f} {syn_c:8.4f} {offset:+8.4f} {mol_dens:9.4f}")

    offs, mols = np.array(offs), np.array(mols)
    off_lo = float(np.median(offs[mols <= np.quantile(mols, 0.34)]))
    off_hi = float(np.median(offs[mols >= np.quantile(mols, 0.66)]))
    print(f"\n  median offset, LOW-molecular windows  = {off_lo:+.4f}")
    print(f"  median offset, HIGH-molecular windows = {off_hi:+.4f}")
    print(f"  overall median offset                 = {np.median(offs):+.4f}  "
          f"(RYA-449 [O I] window ~ -0.0071)")

    print("\n-- PART A verdict ----------------------------------------------")
    if abs(off_lo) >= OFFSET_FLOOR:
        print(f"  NORMALIZATION-SCALE: offset present even in molecule-poor windows "
              f"({off_lo:+.4f}). Global normalization scale issue, NOT missing opacity -- and "
              f"SYSTEMIC (every species' continuum). Fix = re-normalization. Part B proceeds.")
        return 'NORMALIZATION'
    if abs(off_hi) >= OFFSET_OPAC:
        print(f"  MICRO-OPACITY: offset ~0 where molecules are sparse ({off_lo:+.4f}) but grows "
              f"with molecular density ({off_hi:+.4f}) -> synth MISSING molecular (CN) opacity. "
              f"Fix = line-list completeness patch; a blind renorm would MASK it. Part B WITHHELD.")
        return 'MICRO_OPACITY'
    print(f"  AMBIGUOUS: low-mol {off_lo:+.4f}, high-mol {off_hi:+.4f}. Part B WITHHELD; "
          f"post numbers for Ryan.")
    return 'AMBIGUOUS'


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
    verdict = part_a(ctx, a.tmp_dir)
    if verdict == 'NORMALIZATION' or a.force_part_b:
        part_b(ctx, a.tmp_dir)
    else:
        print(f"\nPart B withheld (Part A verdict = {verdict}).")


if __name__ == '__main__':
    main()
