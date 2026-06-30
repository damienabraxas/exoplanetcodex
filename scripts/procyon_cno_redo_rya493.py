"""
scripts/procyon_cno_redo_rya493.py
==================================
RYA-493 — Procyon CNO redo (surgical): re-bank [O/H] regime-matched vs the closed-solar
1D denominator, re-confirm C + re-state the spread arm-matched, and ADD the new N I red
leg via the GALAH PySME NLTE grid (RYA-369). Metals untouched. validate-don't-tune;
differential vs OUR measured Sun, never Asplund; regime-match every differential (RYA-484);
all grids/atlases from Sirius.

Solar denominators (closed, RYA-485/491):  O 1D-NLTE 8.755 · C(3D) 8.449 · N (pending bank).
Procyon (6554 K > 6500 ceiling) -> the 1D-NLTE leg on every Amarsi-grid line (regime-match).
UVES 2013 RED760 (SNR~411) is the telluric-CLEAN anchor epoch (RYA-272); 2002-02 excluded.

  TRACK O — O I 777 (UVES, 1D-NLTE) vs solar 1D 8.755 = the re-banked [O/H]; [O I] 6300
    (HARPS) cross-check. Continuum lever (0.18, UVES-intrinsic, RYA-491) carried, not re-litigated.
  TRACK C — C I 5052/5380 (re-fit, atomic, Bruntt-comparable) + carry Phase-1 CH/C2 ->
    re-state the 0.253 spread per-indicator ARM-MATCHED (RYA-485 Issue 1): VIS C I internal
    spread (Bruntt comparand) vs atomic-vs-molecular spread (our novel product). Not collapsed.
  TRACK N (additive, severable) — N I 7468/8216/8683 (UVES) + GALAH NLTE (warm-Teff δ).
    FINDING to confirm: solar N I was ~-0.015; at 6554 K expect MATERIALLY larger NLTE.

    python -m scripts.procyon_cno_redo_rya493
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import pipeline.cno_synthesis as cs                              # noqa: E402
from pipeline import nlte_cno, pysme_nlte                        # noqa: E402
from config.constants import get_star_params                     # noqa: E402
from scripts.procyon_oi777_continuum_localize_rya478 import (    # noqa: E402
    localize_continuum, fit_oi777, OI777_CONT_WINDOW, OI777_TRIPLET_SPAN,
    PROCYON_VMAC, PROCYON_VSINI, UVES_ANCHOR_R)

RESULTS = ROOT / 'data' / 'results'
TMP = '/tmp/ispec_procyon_cno493'

SOLAR = {'O_1D': 8.755, 'C_3D': 8.449, 'N_banked': 8.202}   # closed-solar denominators
# Phase-1 Procyon C indicators (RYA-348/237) — re-differenced vs the new solar C
PHASE1_C = {'CH_Gband': 8.324, 'CI_5052': 8.501, 'CI_5380': 8.304, 'C2_Swan': 8.557}
PROV_OH_REGIME = 0.066                                        # the provisional regime-matched [O/H]

C_LINES = [(5052.17, 0.40), (5380.34, 0.40)]
N_LINES = [(7468.31, 0.40), (8216.34, 0.40), (8683.40, 0.40)]


def _rule(c='='):
    print(c * 80)


def _params():
    r = get_star_params('procyon')
    return {'teff_K': float(r['teff']), 'logg': float(r['logg']),
            'feh': float(r['feh_ref']), 'vturb_kms': float(r.get('xi', 1.0))}


def _resources():
    ll, iso, chem = cs._load_synth_resources()
    sab = cs.ispec.read_solar_abundances(cs._ISPEC_SOLAR_ABUND_FILE)
    return ll, iso, chem, sab


def _cno_leg(species, wl, params, a_lte):
    lab = nlte_cno.resolve_line(species, wl)
    return float(nlte_cno.cno_nlte_delta(species, lab, params['teff_K'], params['logg'],
                                         params['feh'], params['vturb_kms'], a_lte, leg='1D'))


def track_O(uves_w, uves_f, harps_w, harps_f, params, ll, iso, sab, codes):
    _rule(); print("  TRACK O — re-bank [O/H] (O I 777 UVES 1D-NLTE vs solar 1D 8.755)"); _rule()
    wA = uves_w * 10.0
    wl_win, fl_loc = localize_continuum(wA, uves_f, cont_window=OI777_CONT_WINDOW,
                                        exclude=(OI777_TRIPLET_SPAN,))
    f2 = uves_f.copy(); idx = np.searchsorted(wA, wl_win); inb = (idx >= 0) & (idx < wA.size)
    f2[idx[inb]] = fl_loc[inb]
    r = fit_oi777(uves_w, f2, params, (UVES_ANCHOR_R, PROCYON_VMAC, PROCYON_VSINI),
                  ll, iso, sab, codes, center=8.8)
    a_lte = r['A_X']; d1 = _cno_leg('OI', 7773.0, params, a_lte)
    a_1d = round(a_lte + d1, 3)
    oh = round(a_1d - SOLAR['O_1D'], 3)
    print(f"  O I 777 (UVES-red): A_LTE {a_lte:.3f}  χ²ᵣ {r['red_chi2']:.1f}  1D-NLTE δ {d1:+.3f}  "
          f"A(O) {a_1d}")
    print(f"  [O/H] = {a_1d} − solar-1D {SOLAR['O_1D']} = {oh:+.3f}  (provisional was {PROV_OH_REGIME:+.3f}, "
          f"Δ {oh-PROV_OH_REGIME:+.3f})")
    # [O I] 6300 HARPS cross-check (forbidden, LTE-insensitive)
    o6 = None
    try:
        from scripts.hard_bank_procyon_o_rya483 import fit_oi6300
        a6, chi6, st6 = fit_oi6300(harps_w, harps_f, params, (cs.HARPS_R, PROCYON_VMAC, PROCYON_VSINI),
                                   ll, iso, sab, codes, center=8.8)
        o6 = dict(a_lte=a6, chi2r=chi6, status=st6)
        print(f"  [O I] 6300 (HARPS, forbidden): A_LTE {a6}  χ²ᵣ {chi6}  [cross-check, continuum-limited]")
    except Exception as e:
        print(f"  [O I] 6300 cross-check skipped: {type(e).__name__}")
    return dict(oi777_a_lte=round(a_lte, 3), oi777_chi2r=round(r['red_chi2'], 1),
                oi777_nlte_1d=round(d1, 3), oi777_a_1d=a_1d, oh_rebanked=oh,
                oh_provisional=PROV_OH_REGIME, delta_from_provisional=round(oh - PROV_OH_REGIME, 3),
                oi6300=o6, continuum_sigma=0.186)


def track_C(uves_w, uves_f, harps_w, harps_f, params, ll, iso, sab, chem, sab2):
    _rule(); print("  TRACK C — confirm C + re-state the spread arm-matched"); _rule()
    c_codes = cs._atom_codes(('C',), chem, sab)
    atm = cs._load_atmosphere(params['teff_K'], params['logg'], params['feh'], params['vturb_kms'])
    per = {}
    # C I 5052/5380 are OPTICAL (HARPS); UVES-RED760 starts ~5680 A so it does NOT cover them
    # (the UVES-red independent atomic-C point would be a separate red C I line — not added here).
    for wl, hw in C_LINES:
        r = cs._fit_element(harps_w, harps_f, atm, params, 'C', {'C': 8.46}, c_codes,
                            ((wl - hw, wl + hw),), False, (cs.HARPS_R, PROCYON_VMAC, PROCYON_VSINI),
                            7.3, 9.2, ll, iso, sab, TMP)
        d1 = _cno_leg('CI', wl, params, r['A_X'])
        key = f"CI_{int(wl)}_HARPS"
        per[key] = dict(a_lte=round(r['A_X'], 3), a_1d=round(r['A_X'] + d1, 3),
                        chi2r=round(r['red_chi2'], 1))
        print(f"  C I {wl} (HARPS): A_LTE {per[key]['a_lte']}  1D-NLTE {per[key]['a_1d']}  χ²ᵣ {per[key]['chi2r']}")
    # re-state the spread vs new solar C (Phase-1 indicators re-differenced)
    ch = {k: round(v - SOLAR['C_3D'], 3) for k, v in PHASE1_C.items()}
    vis_ci = [PHASE1_C['CI_5052'], PHASE1_C['CI_5380']]
    vis_ci_spread = round(max(vis_ci) - min(vis_ci), 3)              # Bruntt-comparable
    atomic = [PHASE1_C['CI_5052'], PHASE1_C['CI_5380']]
    molec = [PHASE1_C['CH_Gband'], PHASE1_C['C2_Swan']]
    atomic_molec_spread = round(abs(np.mean(atomic) - np.mean(molec)), 3)   # our product
    full_spread = round(max(PHASE1_C.values()) - min(PHASE1_C.values()), 3)
    print(f"\n  spread RE-STATED arm-matched (RYA-485 Issue 1):")
    print(f"    VIS C I internal (5052/5380) = {vis_ci_spread} dex  <- Bruntt-comparable")
    print(f"    atomic(C I) vs molecular(CH/C2) = {atomic_molec_spread} dex  <- OUR product, no comparand")
    print(f"    full 4-indicator spread {full_spread} (do NOT collapse to one 'vs Bruntt' number)")
    print(f"    [C/H] per indicator vs new solar 8.449: {ch}")
    new_ci5052 = per.get('CI_5052_HARPS', {}).get('a_1d')
    holds = new_ci5052 is not None and abs(new_ci5052 - SOLAR['C_3D']) < 0.20
    print(f"    C I 5052 re-fit (HARPS, 1D-NLTE) {new_ci5052} vs new solar C 8.449 -> "
          f"{'C holds ✓' if holds else 'moved — finding'}")
    return dict(refit_per_line=per, phase1_CH=ch, vis_ci_spread=vis_ci_spread,
                atomic_vs_molecular_spread=atomic_molec_spread, full_spread=full_spread,
                ci5052_refit_harps_1d=new_ci5052, holds=holds)


def track_N(uves_w, uves_f, params, ll, iso, sab, chem):
    _rule(); print("  TRACK N — additive N I red leg (GALAH NLTE, warm-Teff finding)"); _rule()
    n_codes = cs._atom_codes(('N',), chem, sab)
    atm = cs._load_atmosphere(params['teff_K'], params['logg'], params['feh'], params['vturb_kms'])
    per = []
    for wl, hw in N_LINES:
        r = cs._fit_element(uves_w, uves_f, atm, params, 'N', {'N': 7.83}, n_codes,
                            ((wl - hw, wl + hw),), False, (UVES_ANCHOR_R, PROCYON_VMAC, PROCYON_VSINI),
                            6.5, 8.8, ll, iso, sab, TMP)
        per.append(dict(line=wl, a_lte=round(r['A_X'], 3), chi2r=round(r['red_chi2'], 1), status=r['status']))
        print(f"  N I {wl} (UVES-red): A_LTE {per[-1]['a_lte']}  χ²ᵣ {per[-1]['chi2r']}  [{r['status']}]")
    # GALAH NLTE delta at Procyon (the warm-Teff magnitude) — the new capability
    star = {'teff': params['teff_K'], 'logg': params['logg'], 'feh': params['feh'], 'vmic': params['vturb_kms']}
    nlte = None
    try:
        res = pysme_nlte.nlte_delta('N', star=star)
        nlte = round(float(res['delta_median']), 4)
        print(f"  GALAH N I NLTE δ at Procyon (6554 K) = {nlte:+.4f}  (solar was ~-0.015 -> "
              f"{'MATERIALLY LARGER ✓ (warm-Teff, strategy confirmed)' if abs(nlte) > 0.05 else 'still small'})")
    except Exception as e:
        print(f"  GALAH N NLTE unavailable ({type(e).__name__}: {e}) — grid not staged? A(N)_LTE only")
    vals = [p['a_lte'] for p in per if np.isfinite(p['a_lte'])]
    a_lte = round(float(np.median(vals)), 3) if vals else None
    a_nlte = round(a_lte + nlte, 3) if (a_lte is not None and nlte is not None) else None
    nh = round(a_nlte - SOLAR['N_banked'], 3) if a_nlte is not None else None
    print(f"  A(N)_LTE median {a_lte} (σ {round(float(np.std(vals,ddof=1)),3) if len(vals)>1 else 0.0}); "
          f"A(N)_NLTE {a_nlte}; [N/H] vs banked-solar-N {SOLAR['N_banked']} = {nh} "
          f"(FLAG: solar N I not yet banked -> Procyon-N pending-solar-N-bank)")
    return dict(per_line=per, a_lte=a_lte, nlte_delta_procyon=nlte, nlte_solar=-0.015,
                warm_teff_larger=(nlte is not None and abs(nlte) > 0.05), a_nlte=a_nlte,
                nh_provisional=nh, solar_n_bank_pending=True)


def main():
    Path(TMP).mkdir(parents=True, exist_ok=True)
    ll, iso, chem, sab = _resources()
    cs._CHEM = chem
    codes = cs._atom_codes(('C', 'N', 'O', 'Ni'), chem, sab)
    params = _params()
    print(f"\n  Procyon: Teff {params['teff_K']:.0f} logg {params['logg']:.2f} "
          f"[Fe/H] {params['feh']:+.2f} ξ {params['vturb_kms']:.1f}")

    uves_w, uves_f = cs.resolve_arm_spectrum('procyon', cs.star_arm_registry('procyon')['uves'])
    harps_w, harps_f = cs._load_observed_spectrum('procyon')

    O = track_O(uves_w, uves_f, harps_w, harps_f, params, ll, iso, sab, codes)
    C = track_C(uves_w, uves_f, harps_w, harps_f, params, ll, iso, sab, chem, sab)
    N = track_N(uves_w, uves_f, params, ll, iso, sab, chem)

    _rule(); print("  PROCYON CNO REDO — VERDICT (surgical; metals untouched; STOP at verdict)"); _rule()
    print(f"  O: [O/H] RE-BANKED {O['oh_rebanked']:+.3f} ± {O['continuum_sigma']} (777 1D-NLTE {O['oi777_a_1d']} "
          f"− solar-1D {SOLAR['O_1D']}); Δ from provisional {O['delta_from_provisional']:+.3f}")
    print(f"  C: VIS C I spread {C['vis_ci_spread']} (Bruntt-comparable) | atomic-vs-molecular {C['atomic_vs_molecular_spread']} "
          f"(our product); C I 5052 re-fit {C['ci5052_refit_harps_1d']} vs solar 8.449 ({'holds' if C['holds'] else 'moved'})")
    print(f"  N: ADDITIVE A(N)_NLTE {N['a_nlte']} (GALAH δ {N['nlte_delta_procyon']} vs solar -0.015 = "
          f"{'warm-Teff larger ✓' if N['warm_teff_larger'] else 'small'}); [N/H] {N['nh_provisional']} pending-solar-N-bank")
    out = {'ticket': 'RYA-493 Procyon CNO redo', 'validate_dont_tune': True, 'metals_untouched': True,
           'solar_denominators': SOLAR, 'params': params, 'O': O, 'C': C, 'N': N,
           'gaps': ['IR CO 2.3µm absent (CRIRES+/STAGGER)', 'IR O I 844/926 Procyon-IR absent (APOGEE H-band 1.5-1.7µm)']}
    RESULTS.mkdir(parents=True, exist_ok=True)
    (RESULTS / 'procyon_cno_redo_rya493.json').write_text(json.dumps(out, indent=2, default=str))
    print(f"\n  [out] {RESULTS / 'procyon_cno_redo_rya493.json'}")
    return out


if __name__ == '__main__':
    main()
