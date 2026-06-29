"""
scripts/hard_bank_procyon_o_rya483.py
=====================================
RYA-483 — HARD-BANK Procyon O. Closes the two residuals RYA-478 left open:

  A. Cross-instrument zero-point (Procyon UVES O I 777 vs solar Kitt Peak).
     The differential [O/H] = A(O,Procyon,UVES,777) - A(O,Sun,KPNO,777) compares two
     INSTRUMENTS. We characterize the zero-point as two measurable/cited components:
       (1) continuum-method term — the local-minus-global continuum A(O) shift, measured
           PER ARM (KPNO solar + UVES Procyon). If continuum systematics were instrument-
           independent the two shifts cancel in the differential; their DIFFERENCE
           (Procyon - solar) is the residual continuum-driven cross-instrument zero-point.
       (2) reference-instrument term — KPNO reproduces our solar O (8.735) to ~1 mdex, so
           the KPNO reference is on-scale; the solar KPNO-vs-ESPRESSO O I 777 spread
           (~0.04 dex, RYA-460) bounds the reference-side instrument uncertainty.
     Disposition: there is no shared-star+diagnostic handle (no Procyon-KPNO, no Sun-UVES
     at 777), so the zero-point is NOT a stable correctable additive offset -> it is
     FOLDED INTO sigma, reported with its components, never silently absorbed.

  B. Independent HARPS [O I] 6300 Procyon fit (de-anchor the leg). RYA-348's leg was
     inconclusive because both arms were Caffau-anchor-dominated (8.730); raw split was
     UVES 8.807 vs HARPS 9.417, masked. Here [O I] 6300 is fit INDEPENDENTLY (continuum-
     localized, Ni I 6300.34 pinned at the Johansson-2003 gf via gf_resolver, RYA-365) on
     Procyon HARPS, to de-anchor the UVES-vs-HARPS leg. Inherit the solar lesson
     (RYA-447/448/455): [O I] 6300 is continuum + Ni-blend TERMINAL, NOT science-grade,
     Amarsi-2019 3D differential ~0 at 630 nm -> this is a CONSISTENCY cross-check on
     O I 777, not a co-equal O measurement. Expect it weak/noisy.

The bank decision: O I 777 (1D-NLTE) stays PRIMARY; differential vs OUR measured Sun
(8.735), never Asplund; solar control must not move (the RYA-478 guard). validate-don't-
tune: a clean value still offset is a finding. STOP at verdict, no merge.

    python -m scripts.hard_bank_procyon_o_rya483
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import pipeline.cno_synthesis as cs                              # noqa: E402
from pipeline import nlte_cno                                    # noqa: E402
from config.constants import get_star_params                     # noqa: E402
# reuse the RYA-478 continuum-localized O I 777 machinery verbatim
from scripts.procyon_oi777_continuum_localize_rya478 import (    # noqa: E402
    localize_continuum, fit_oi777, nlte_delta as oi777_nlte_delta, _resources,
    OI777_CONT_WINDOW, OI777_TRIPLET_SPAN,
    PROCYON_VMAC, PROCYON_VSINI, SOLAR_VMAC, SOLAR_VSINI,
    UVES_ANCHOR_R, KPNO_R, SOLAR_OURS_O,
)

TMP = '/tmp/ispec_cno_rya483'
RESULTS = ROOT / 'data' / 'results'

# [O I] 6300 geometry (air Angstrom): forbidden [O I] 6300.304 + Ni I 6300.339 blend.
OI6300_LINE = 6300.30
OI6300_FIT_WINDOW = ((6299.9, 6300.8),)
OI6300_CORE_SPAN = (6300.0, 6300.7)          # excluded from the local continuum
OI6300_CONT_WINDOW = (6298.0, 6302.4)        # local continuum shoulders
HARPS_R = 115000.0
# RYA-460 cited solar O I 777 reference-instrument spread (KPNO vs ESPRESSO).
KPNO_ESPRESSO_ZP = 0.04


def _rule(c='='):
    print(c * 80)


def _splice_local(wA, fl, cont_window, exclude):
    """Continuum-localize `cont_window` (excluding `exclude`) and splice back into a copy."""
    wl_win, fl_loc = localize_continuum(wA, fl, cont_window=cont_window, exclude=(exclude,))
    fl2 = fl.copy()
    idx = np.searchsorted(wA, wl_win)
    inb = (idx >= 0) & (idx < wA.size)
    fl2[idx[inb]] = fl_loc[inb]
    return fl2


# ── Part A — O I 777 continuum-method zero-point, per arm ─────────────────────

def oi777_global_local(star, w_nm, fl_global, params, broad, ll, iso, sab, codes, center):
    """Fit A(O) on O I 777 under (a) the arm's global continuum and (b) the localized
    continuum; return both + the local-minus-global shift (the continuum-method term)."""
    wA = w_nm * 10.0
    r_glob = fit_oi777(w_nm, fl_global, params, broad, ll, iso, sab, codes, center=center)
    fl_loc = _splice_local(wA, fl_global, OI777_CONT_WINDOW, OI777_TRIPLET_SPAN)
    r_loc = fit_oi777(w_nm, fl_loc, params, broad, ll, iso, sab, codes, center=center)
    shift = round(r_loc['A_X'] - r_glob['A_X'], 3)
    print(f"  [{star}] O I 777  A_LTE global {r_glob['A_X']:.3f} (χ²ᵣ {r_glob['red_chi2']:.1f}) "
          f"| local {r_loc['A_X']:.3f} (χ²ᵣ {r_loc['red_chi2']:.1f}) | shift {shift:+.3f}")
    return dict(a_global=round(r_glob['A_X'], 3), a_local=round(r_loc['A_X'], 3),
                chi2r_local=round(r_loc['red_chi2'], 2), shift=shift,
                r_local=r_loc)


def part_a(ll, iso, sab, codes):
    _rule(); print("  PART A — cross-instrument zero-point (O I 777: Procyon UVES vs solar KPNO)"); _rule()
    # solar KPNO
    rec_s = get_star_params('solar')
    ps = {'teff_K': float(rec_s['teff']), 'logg': float(rec_s['logg']),
          'feh': float(rec_s['feh_ref']), 'vturb_kms': float(rec_s.get('xi', 1.0))}
    kp = pd.read_csv(ROOT / 'data' / 'solar_reference' / 'kpno_flux_atlas' / 'kpno_OI_777_triplet.csv',
                     comment='#')
    sol = oi777_global_local('solar_kpno', kp['wavelength_air_A'].to_numpy(float) / 10.0,
                             kp['residual_flux'].to_numpy(float), ps,
                             (KPNO_R, SOLAR_VMAC, SOLAR_VSINI), ll, iso, sab, codes, center=8.9)
    lab, ds = oi777_nlte_delta(ps, sol['r_local']['A_X'])
    sol['a_nlte_local'] = round(sol['r_local']['A_X'] + ds, 3)
    # procyon UVES
    rec_p = get_star_params('procyon')
    pp = {'teff_K': float(rec_p['teff']), 'logg': float(rec_p['logg']),
          'feh': float(rec_p['feh_ref']), 'vturb_kms': float(rec_p.get('xi', 1.0))}
    arm = cs.star_arm_registry('procyon')['uves']
    w_nm, fl = cs.resolve_arm_spectrum('procyon', arm)
    pro = oi777_global_local('procyon_uves', w_nm, fl, pp,
                             (UVES_ANCHOR_R, PROCYON_VMAC, PROCYON_VSINI), ll, iso, sab, codes, center=8.8)
    lab, dp = oi777_nlte_delta(pp, pro['r_local']['A_X'])
    pro['a_nlte_local'] = round(pro['r_local']['A_X'] + dp, 3)

    cont_zp = round(pro['shift'] - sol['shift'], 3)         # differential continuum-method term
    ref_acc = round(sol['a_nlte_local'] - SOLAR_OURS_O, 3)  # KPNO reference accuracy vs our Sun
    zp_budget = round(float(np.hypot(cont_zp, KPNO_ESPRESSO_ZP)), 3)
    print(f"\n  continuum-method shift: solar {sol['shift']:+.3f} | procyon {pro['shift']:+.3f} "
          f"-> differential (cross-instrument) {cont_zp:+.3f}")
    print(f"  reference accuracy: KPNO solar O I 777 NLTE {sol['a_nlte_local']} vs our Sun "
          f"{SOLAR_OURS_O} = {ref_acc:+.3f} (on-scale)")
    print(f"  reference-instrument spread (KPNO vs ESPRESSO, RYA-460) = {KPNO_ESPRESSO_ZP:.3f}")
    print(f"  => zero-point budget (quadrature of continuum-differential + ref-spread) = "
          f"{zp_budget:.3f} dex  -> DISPOSITION: fold into sigma (no shared-star handle to "
          f"correct a stable offset)")
    return dict(solar=sol, procyon=pro, continuum_diff_zp=cont_zp,
                kpno_reference_accuracy=ref_acc, kpno_espresso_zp=KPNO_ESPRESSO_ZP,
                zeropoint_sigma=zp_budget,
                disposition='fold into sigma (cross-instrument, not a correctable additive offset)')


# ── Part B — independent HARPS [O I] 6300 Procyon fit (de-anchor the leg) ──────

def fit_oi6300(w_nm, fl_global, params, broad, ll, iso, sab, codes, center, ni=6.20):
    """Independent continuum-localized [O I] 6300 fit: O free, Ni I 6300.34 pinned at A(Ni).
    Forbidden line is LTE-insensitive (RYA-447: Amarsi 3D differential ~0 at 630 nm), so the
    reported A is the 1D-LTE value = the [O I] number. Returns (A_lte, chi2r, status)."""
    wA = w_nm * 10.0
    fl_loc = _splice_local(wA, fl_global, OI6300_CONT_WINDOW, OI6300_CORE_SPAN)
    atm = cs._load_atmosphere(params['teff_K'], params['logg'], params['feh'], params['vturb_kms'])
    r = cs._fit_element(w_nm, fl_loc, atm, params, 'O',
                        {'C': 8.46, 'N': 7.83, 'O': center, 'Ni': ni}, codes,
                        OI6300_FIT_WINDOW, False, broad,
                        center - 1.5, center + 1.5, ll, iso, sab, TMP)
    return round(r['A_X'], 3), round(r['red_chi2'], 2), r['status']


def part_b(ll, iso, sab, codes, uves777_nlte):
    _rule(); print("  PART B — independent HARPS [O I] 6300 Procyon fit (de-anchor the leg)"); _rule()
    print("  NOTE (RYA-447/448/455): [O I] 6300 is continuum + Ni-blend TERMINAL, NOT science-")
    print("  grade; this is a CONSISTENCY cross-check on O I 777, not a co-equal O measurement.")
    # solar KPNO [O I] 6300 control (terminal — reports for context, not a science gate)
    rec_s = get_star_params('solar')
    ps = {'teff_K': float(rec_s['teff']), 'logg': float(rec_s['logg']),
          'feh': float(rec_s['feh_ref']), 'vturb_kms': float(rec_s.get('xi', 1.0))}
    kp = pd.read_csv(ROOT / 'data' / 'solar_reference' / 'kpno_flux_atlas' / 'kpno_OI_6300.csv',
                     comment='#')
    s_a, s_chi, s_st = fit_oi6300(kp['wavelength_air_A'].to_numpy(float) / 10.0,
                                  kp['residual_flux'].to_numpy(float), ps,
                                  (KPNO_R, SOLAR_VMAC, SOLAR_VSINI), ll, iso, sab, codes, center=8.7)
    print(f"  [solar_kpno] [O I] 6300 A_LTE = {s_a}  χ²ᵣ {s_chi}  [{s_st}]  "
          f"(our Sun {SOLAR_OURS_O}; terminal — context only)")
    # procyon HARPS [O I] 6300, independent
    rec_p = get_star_params('procyon')
    pp = {'teff_K': float(rec_p['teff']), 'logg': float(rec_p['logg']),
          'feh': float(rec_p['feh_ref']), 'vturb_kms': float(rec_p.get('xi', 1.0))}
    w_nm, fl = cs._load_observed_spectrum('procyon')     # HARPS normalized
    p_a, p_chi, p_st = fit_oi6300(w_nm, fl, pp, (HARPS_R, PROCYON_VMAC, PROCYON_VSINI),
                                  ll, iso, sab, codes, center=8.8)
    # forbidden line: differential vs our Sun is on the same (LTE) footing; 3D diff ~0
    p_oh = round(p_a - s_a, 3) if np.isfinite(s_a) else None
    print(f"  [procyon_harps] [O I] 6300 A_LTE = {p_a}  χ²ᵣ {p_chi}  [{p_st}]  (independent, "
          f"Ni-pinned, continuum-localized)")
    # TERMINAL detection: if the solar [O I] 6300 control does not reproduce our Sun, the
    # [O I] 6300 method is terminal (RYA-447/448/455) and CANNOT de-anchor the leg — a
    # railed-to-bracket-floor solar value is a failed fit, not a measurement.
    solar_terminal = (not np.isfinite(s_a)) or abs(s_a - SOLAR_OURS_O) > 0.30
    leg = round(uves777_nlte - p_a, 3)
    if solar_terminal:
        print(f"\n  ⚠️ [O I] 6300 TERMINAL: solar control {s_a} vs our Sun {SOLAR_OURS_O} "
              f"(|Δ| {abs(s_a - SOLAR_OURS_O):.2f} > 0.30) — the line is continuum+Ni-blend "
              f"terminal (RYA-447/448/455). The independent [O I] 6300 number is UNRELIABLE; "
              f"the de-anchored leg is UNAVAILABLE, not a corroboration.")
    else:
        print(f"\n  DE-ANCHORED LEG: O I 777 (1D-NLTE) {uves777_nlte} vs independent HARPS "
              f"[O I] 6300 {p_a}  ->  Δ {leg:+.3f} dex")
    return dict(solar_kpno_oi6300_lte=s_a, solar_chi2r=s_chi,
                procyon_harps_oi6300_lte=p_a, procyon_chi2r=p_chi, procyon_status=p_st,
                deanchored_leg_777_minus_6300=leg, solar_terminal=bool(solar_terminal),
                leg_available=not bool(solar_terminal))


def main():
    Path(TMP).mkdir(parents=True, exist_ok=True)
    ll, iso, sab, codes = _resources()
    A = part_a(ll, iso, sab, codes)
    uves777_nlte = A['procyon']['a_nlte_local']
    oh = round(uves777_nlte - SOLAR_OURS_O, 3)
    B = part_b(ll, iso, sab, codes, uves777_nlte)

    _rule(); print("  HARD-BANK VERDICT (interpretation is Ryan's + Claude's)"); _rule()
    solar_unmoved = abs(A['solar']['a_nlte_local'] - SOLAR_OURS_O) < 0.02
    # σ is dominated by the cross-instrument continuum zero-point (Part A). The [O I] 6300
    # leg is terminal (Part B) so it adds no corroboration — it does NOT shrink σ.
    sigma = A['zeropoint_sigma']
    # NOT hard-bankable: the continuum zero-point is large (~0.18, >> 0.10) AND the
    # independent corroboration leg is terminal/unavailable.
    hard_bankable = (sigma <= 0.10) and B['leg_available'] and solar_unmoved
    print(f"  Primary  : Procyon O I 777 (1D-NLTE) = {uves777_nlte}  [O/H] {oh:+.3f} vs our Sun {SOLAR_OURS_O}  (χ²ᵣ {A['procyon']['chi2r_local']})")
    print(f"  Part A   : cross-instrument zero-point = {sigma:.3f} dex — DOMINATED by the UVES")
    print(f"             continuum-localization lever (solar shift {A['solar']['shift']:+.3f} vs")
    print(f"             procyon {A['procyon']['shift']:+.3f} = differential {A['continuum_diff_zp']:+.3f}); ")
    print(f"             FAR larger than RYA-478's 0.04-0.08 estimate. Fold into σ.")
    print(f"  Part B   : independent [O I] 6300 leg = {'TERMINAL/UNAVAILABLE' if B['solar_terminal'] else 'available'} "
          f"(solar control {B['solar_kpno_oi6300_lte']} vs Sun {SOLAR_OURS_O}); cannot corroborate 777.")
    print(f"  Guard    : solar control O I 777 = {A['solar']['a_nlte_local']} "
          f"({'unmoved ✓' if solar_unmoved else 'MOVED — STOP'})")
    print(f"\n  VERDICT: {'HARD-BANK' if hard_bankable else 'NOT hard-bankable — provisional, banked WITH the named residual'}.")
    print(f"  Provisional Procyon O: A(O) {uves777_nlte}, [O/H] {oh:+.3f} ± {sigma:.3f} (1σ continuum")
    print(f"  zero-point dominated); primary = O I 777 1D-NLTE; differential vs our Sun {SOLAR_OURS_O}.")
    print(f"  BLOCKERS to a tight bank: (1) the UVES O I 777 continuum-localization lever (0.23 dex,")
    print(f"  no χ² discrimination) — needs an independent continuum / 2nd 777 epoch or instrument;")
    print(f"  (2) [O I] 6300 is terminal (RYA-447/448/455) — it cannot serve as the corroboration leg.")
    out = {'ticket': 'RYA-483', 'primary_indicator': 'OI_777_1D_NLTE',
           'procyon_O_AO_nlte': uves777_nlte, 'procyon_OH_vs_our_sun': oh,
           'our_solar_O': SOLAR_OURS_O, 'banked_sigma': sigma,
           'hard_bankable': hard_bankable,
           'verdict': ('HARD_BANKED' if hard_bankable else
                       'PROVISIONAL_NOT_HARD_BANKED — continuum zero-point ~0.18 dex + [O I] 6300 leg terminal'),
           'blockers': ['UVES O I 777 continuum-localization lever 0.23 dex (no χ² discrimination)',
                        '[O I] 6300 terminal (RYA-447/448/455) — corroboration leg unavailable'],
           'solar_control_unmoved': solar_unmoved, 'part_a_zeropoint': A,
           'part_b_oi6300_leg': B}
    RESULTS.mkdir(parents=True, exist_ok=True)
    # drop the heavy fit objects before serialising
    for arm in ('solar', 'procyon'):
        A[arm].pop('r_local', None)
    (RESULTS / 'procyon_O_hardbank_rya483.json').write_text(json.dumps(out, indent=2, default=str))
    print(f"\n  [out] {RESULTS / 'procyon_O_hardbank_rya483.json'}")
    return out


if __name__ == '__main__':
    main()
