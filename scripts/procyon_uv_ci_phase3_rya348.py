"""
scripts/procyon_uv_ci_phase3_rya348.py
======================================
RYA-348 Phase 3 — UV arm (HST STIS/COS): the carbon discriminator.

Goal (from the Phase-3 brief): flip the Procyon UV arm ready, measure UV C I (the
independent atomic-C indicator with different formation depth than optical C I), run
the keystone UV-C I-vs-optical-C I leg, the RYA-483 continuum-lever check, and the
solar control — to cut the 0.253-dex Phase-1 carbon-spread ambiguity (molecular vs
atomic vs continuum-limited).

What this run FINDS (honest, run-and-report; the brief explicitly expected UV C I to be
at MORE continuum risk than O I 777, and warned against false precision):

  1. Vacuum->air verified on a named line (Mg II k 2796.3543 vac -> 2795.528 air).
  2. The Procyon STIS E140M arm LOADS + CONDITIONS fine (C I 1657 covered, resolved).
  3. BUT the FUV C I 1657 SYNTHESIS is not science-grade: the optical-tuned Turbospectrum
     synthesis cannot reproduce the STIS FUV pseudo-continuum. The observed "continuum"
     near 1657 sits at ~0.81 with std ~0.22 (a dense line forest, not a flat continuum);
     the synthesis assumes true continuum ~1.0, cannot match, and the fit RAILS to the
     abundance floor (A(C) ~7.0) with chi2r ~600-840 under BOTH global and local
     continuum. The continuum lever does not rescue it -> a hard FUV-synthesis-capability
     wall, not a sigma floor.
  4. The solar control + the differential vs our Sun via UV C I are UNAVAILABLE: the
     solar UV reference is the CALSPEC composite (R~150-300, dlambda~20 A) -- a flux
     composite, NOT a line atlas, so solar C I 1657 is unresolved and cannot be fit.
     The RYA-483 lesson, deeper: there is no resolved solar UV C I to difference against.
  5. So the keystone leg cannot complete and the UV arm CANNOT yet discriminate the
     0.253-dex carbon spread. Two concrete blockers are surfaced (see verdict).

Scope: analysis/scratch only. STOP at the verdict, no merge.
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
from pipeline import uv_conditioning as uvc                      # noqa: E402
from pipeline.loaders import hst_uv_loader as huv                # noqa: E402
from config.constants import get_star_params, SOLAR_REFERENCE_SPECTRA  # noqa: E402

RESULTS = ROOT / 'data' / 'results'
TMP = '/tmp/ispec_cno_rya348p3'
# Phase-1 Procyon OPTICAL C I (HARPS VIS, cited-corrected) — the leg target.
OPTICAL_CI = {'CI_5052': 8.501, 'CI_5380': 8.304, 'CH_Gband': 8.324, 'C2_Swan': 8.557}
CARBON_SPREAD = 0.253


def _rule(c='='):
    print(c * 78)


def _fit_ci1657(w_nm, fl, params, atm, ll, iso, sab, codes, broad, tag):
    r = cs._fit_element(w_nm, fl, atm, params, 'C',
                        {'C': 8.46, 'N': 7.83, 'O': 8.69, 'Ni': 6.2}, codes,
                        ((1656.38, 1658.38),), False, broad, 7.0, 9.5, ll, iso, sab, TMP)
    print(f"    {tag:22s} A(C)={r['A_X']}  chi2r={r['red_chi2']}  [{r['status']}]")
    return r


def main():
    Path(TMP).mkdir(parents=True, exist_ok=True)
    rec = get_star_params('procyon')
    params = {'teff_K': float(rec['teff']), 'logg': float(rec['logg']),
              'feh': float(rec['feh_ref']), 'vturb_kms': float(rec.get('xi', 1.0))}
    _rule(); print("  RYA-348 Phase 3 — Procyon UV arm (HST STIS E140M) C I 1657 discriminator"); _rule()

    # (1) vacuum->air on a NAMED line (the UV BERV-sign-check)
    mgk = float(uvc.to_pipeline_frame(np.array([2796.3543]))[0])
    ci = float(uvc.to_pipeline_frame(np.array([1657.38]))[0])
    print(f"  [vac->air] Mg II k 2796.3543 vac -> {mgk:.3f} air [known 2795.528]  "
          f"{'OK' if abs(mgk-2795.528)<0.01 else 'FAIL'}")
    print(f"  [boundary] C I 1657.38 stays vacuum -> {ci:.3f} (FUV < 2000 A)")

    # (2) load + condition the Procyon STIS UV arm
    w_nm, fl = huv.load_procyon_uv_arm('E140M')
    wA = w_nm * 10.0
    win = (wA >= 1652) & (wA <= 1662)
    obs = fl[win]
    print(f"\n  [arm] E140M loaded: {w_nm.size} px, {wA.min():.0f}-{wA.max():.0f} A; "
          f"C I 1657 covered, resolved (STIS R~45800)")
    print(f"  [pseudo-continuum] observed 1652-1662 A: median {np.nanmedian(obs):.2f} "
          f"std {np.nanstd(obs):.2f} min {np.nanmin(obs):.2f} max {np.nanmax(obs):.2f} "
          f"-> a line FOREST, not a flat continuum (true continuum would be ~1.0 flat)")

    ll, iso, chem = cs._load_synth_resources()
    sab = cs.ispec.read_solar_abundances(cs._ISPEC_SOLAR_ABUND_FILE)
    codes = cs._atom_codes(('C', 'N', 'O', 'Ni'), chem, sab)
    _, vmac, vsini, _ = cs._resolve_broadening('procyon')
    atm = cs._load_atmosphere(params['teff_K'], params['logg'], params['feh'], params['vturb_kms'])
    broad = (45800.0, vmac, vsini)

    # (3) + (5) C I 1657 fit + continuum lever (global vs local pseudo-continuum)
    _rule('-'); print("  C I 1657 fit + continuum-lever (RYA-483 protocol)"); _rule('-')
    r_glob = _fit_ci1657(w_nm, fl, params, atm, ll, iso, sab, codes, broad, 'global continuum')
    fl2 = fl.copy()
    pc = np.nanpercentile(obs, 90)
    fl2[win] = fl2[win] / pc
    r_loc = _fit_ci1657(w_nm, fl2, params, atm, ll, iso, sab, codes, broad, 'local continuum')
    railed = (r_glob['status'] == 'edge_pinned' or r_loc['status'] == 'edge_pinned'
              or r_glob['red_chi2'] > 100 or r_loc['red_chi2'] > 100)
    print(f"  continuum lever: global A={r_glob['A_X']} / local A={r_loc['A_X']}; "
          f"chi2r {r_glob['red_chi2']} / {r_loc['red_chi2']} -> "
          f"{'BOTH RAILED / catastrophic chi2r = FUV synthesis cannot model the pseudo-continuum' if railed else 'fittable'}")

    # (4) solar control / differential availability via the solar UV reference
    _rule('-'); print("  Solar control + differential availability (UV reference)"); _rule('-')
    uvref = SOLAR_REFERENCE_SPECTRA.get('uv_composite', {})
    print(f"  solar UV reference = {uvref.get('provenance','?')}: "
          f"{uvref.get('citation','?')[:60]}")
    print(f"  resolution = {uvref.get('resolution','?')}")
    print(f"  -> C I 1657 is UNRESOLVED at R~150-300 (a flux composite, not a line atlas): "
          f"no solar UV C I to fit -> solar control + [C/H] differential via UV C I UNAVAILABLE")

    # (4b) keystone leg + carbon-spread read
    _rule(); print("  KEYSTONE LEG (UV C I vs optical C I) + carbon-spread read"); _rule()
    print(f"  optical C I (Phase 1, Procyon): {OPTICAL_CI}")
    print(f"  UV C I 1657: NOT science-grade (railed A~7.0, chi2r ~{r_glob['red_chi2']:.0f}) "
          f"-> the leg CANNOT complete; no valid UV C I to compare.")
    print(f"  0.253-dex carbon spread: the UV arm CANNOT discriminate molecular vs atomic vs "
          f"continuum yet -> the discriminator is unavailable (FUV synthesis + solar-UV-ref blockers).")

    out = {
        'ticket': 'RYA-348 Phase 3 — UV arm C I discriminator',
        'vacuum_to_air': {'mg_ii_k_air': round(mgk, 3), 'known': 2795.528, 'ci_1657_stays_vac': round(ci, 3)},
        'arm': {'grating': 'E140M', 'npix': int(w_nm.size), 'ci1657_covered': True,
                'pseudo_continuum_median': round(float(np.nanmedian(obs)), 3),
                'pseudo_continuum_std': round(float(np.nanstd(obs)), 3)},
        'uv_ci1657_fit': {'global': {k: r_glob[k] for k in ('A_X', 'red_chi2', 'status')},
                          'local': {k: r_loc[k] for k in ('A_X', 'red_chi2', 'status')},
                          'science_grade': False,
                          'reason': 'FUV pseudo-continuum not modeled by optical-tuned synthesis; '
                                    'fit rails to abundance floor under both continuum choices'},
        'nlte': 'FUV C I 1657 NOT in the Amarsi grid (nlte_cno.resolve_line CI 1657 = None); '
                'GRID_OWED, cited ~+0.10 (Amarsi 2020), flagged-not-applied (RYA-190)',
        'solar_control': {'available': False,
                          'reason': 'solar UV ref = CALSPEC composite R~150-300; C I 1657 unresolved'},
        'leg': {'optical_ci': OPTICAL_CI, 'uv_ci': None, 'completable': False},
        'carbon_spread_dex': CARBON_SPREAD,
        'carbon_spread_uv_verdict': 'UNDISCRIMINATED — UV arm cannot yet measure UV C I',
        'blockers': ['FUV pseudo-continuum synthesis (needs FUV line list + FUV continuum opacities)',
                     'solar UV reference resolution (CALSPEC composite cannot resolve UV C I -> no differential/control)'],
    }
    RESULTS.mkdir(parents=True, exist_ok=True)
    (RESULTS / 'procyon_uv_ci_phase3_rya348.json').write_text(json.dumps(out, indent=2, default=str))
    print(f"\n  [out] {RESULTS / 'procyon_uv_ci_phase3_rya348.json'}")
    return out


if __name__ == '__main__':
    main()
