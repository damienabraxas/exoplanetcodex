"""
scripts/solar_cno_redo_rya485.py
================================
RYA-485/486 SOLAR CNO REDO — re-derive measured solar O (the 4-reference O I 777
swap-test), record both RT legs (3D-NLTE + 1D-NLTE) for regime-matched differentials,
and measure solar S I (multiplet 8). validate-don't-tune: the goal is OUR measured solar
value with honest per-arm / per-reference spreads — the Asplund numbers are comparisons,
never fit targets. All reference atlases sourced from Sirius (RYA-481/477).

Headline — O I 777 three-reference swap-test (RYA-484 Lever-3):
  Fit solar A(O) on O I 777 against each independent reference with the SAME method
  (RYA-478 continuum-localize + Amarsi-2019 3D-NLTE):
    * Kurucz-1984 KPNO   (current reference, in-repo)
    * Reiners-2016 IAG   (CDS, Sirius)
    * Baker-2020 IAG-tel (Zenodo, Sirius — telluric-corrected)
    * Wallace-2011 KPNO  (NSO, Sirius — a second KPNO variant)
  Diagnostic: IAG vs our-KPNO disagree -> Wallace breaks the tie (~Kurucz = KPNO-family
  systematic; ~IAG = Kurucz-reduction-specific). All agree -> the continuum is not
  reference-driven and the RYA-483 sigma is intrinsic.

Regime-match (RYA-484): the Sun is the 3D-NLTE leg (<6500 K), but warm stars (Procyon)
need the 1D-NLTE solar denominator. Record solar_A(O)_3D AND solar_A(O)_1D so downstream
differentials regime-match (the ~0.02 dex Procyon-O fix at its source, RYA-485 Issue 2).

Sulphur (RYA-486): S I multiplet-8 6748.68/6757.15 (VIS/HARPS), small Amarsi-2025 NLTE.
MEASURE it (the EW-pool banked 7.753 is curation-owed); report A(S)_LTE + A(S)_NLTE.

    python -m scripts.solar_cno_redo_rya485
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
from pipeline.uv_conditioning import vac_to_air                  # noqa: E402
from config.constants import get_star_params                     # noqa: E402
from scripts.procyon_oi777_continuum_localize_rya478 import (    # noqa: E402
    localize_continuum, fit_oi777, _resources,
    OI777_CONT_WINDOW, OI777_TRIPLET_SPAN, SOLAR_VMAC, SOLAR_VSINI, KPNO_R, SOLAR_OURS_O)

RESULTS = ROOT / 'data' / 'results'
SWAP = ROOT / 'data' / 'solar_reference_swaptest'
TMP = '/tmp/ispec_cno_redo485'

# 371-banked solar denominator (the values this redo re-checks).
BANKED = {'C': 8.491, 'N': 8.202, 'O': 8.735, 'S': 7.753}
ASPLUND = {'C': 8.46, 'N': 7.83, 'O': 8.69, 'S': 7.12}          # COMPARISON ONLY, not targets

# O I 777 references. Kurucz = air+residual already; the FTS atlases are wavenumber(cm-1)/vac.
REFS = [
    ('kurucz_kpno', 'csv', 'data/solar_reference/kpno_flux_atlas/kpno_OI_777_triplet.csv', KPNO_R),
    ('reiners_iag', 'wavenumber', 'data/solar_reference_swaptest/reiners_oi777_vac.txt', 1_000_000.0),
    ('baker_iag_tel', 'wavenumber', 'data/solar_reference_swaptest/baker_oi777_vac.txt', 1_000_000.0),
    ('wallace_kpno', 'wavenumber', 'data/solar_reference_swaptest/wallace_kpno_oi777_vac.txt', 700_000.0),
]
# S I multiplet-8 (VIS/HARPS), Amarsi-2025 NLTE grid (small).
S_LINES = [(6748.68, 0.35), (6757.15, 0.35)]


def _rule(c='='):
    print(c * 80)


def _load_ref(kind, path):
    """Return (wave_air_A, flux) for an O I 777 reference, sorted by wavelength."""
    p = ROOT / path
    if kind == 'csv':
        df = pd.read_csv(p, comment='#')
        wA = df['wavelength_air_A'].to_numpy(float)
        fl = df['residual_flux'].to_numpy(float)
    else:  # wavenumber cm-1 (vacuum) -> air Angstrom
        d = np.loadtxt(p)
        wA = vac_to_air(1e8 / d[:, 0])
        fl = d[:, 1]
    o = np.argsort(wA)
    return wA[o], fl[o]


def _oi777_nlte(params, a_lte, leg):
    lab = nlte_cno.resolve_line('OI', 7773.0)
    return float(nlte_cno.cno_nlte_delta('OI', lab, params['teff_K'], params['logg'],
                                         params['feh'], params['vturb_kms'], a_lte, leg=leg))


def swap_test(ll, iso, sab, codes):
    _rule(); print("  HEADLINE — O I 777 three-reference swap-test (RYA-484 Lever-3)"); _rule()
    rec = get_star_params('solar')
    params = {'teff_K': float(rec['teff']), 'logg': float(rec['logg']),
              'feh': float(rec['feh_ref']), 'vturb_kms': float(rec.get('xi', 1.0))}
    rows = []
    for name, kind, path, R in REFS:
        wA, fl = _load_ref(kind, path)
        w_nm = wA / 10.0
        # local-continuum normalize the triplet window (the RYA-478 method, identical per ref)
        wl_win, fl_loc = localize_continuum(wA, fl, cont_window=OI777_CONT_WINDOW,
                                            exclude=(OI777_TRIPLET_SPAN,))
        fl2 = fl.copy()
        idx = np.searchsorted(wA, wl_win); inb = (idx >= 0) & (idx < wA.size)
        fl2[idx[inb]] = fl_loc[inb]
        r = fit_oi777(w_nm, fl2, params, (R, SOLAR_VMAC, SOLAR_VSINI), ll, iso, sab, codes, center=8.9)
        a_lte = r['A_X']
        d3 = _oi777_nlte(params, a_lte, '3D'); d1 = _oi777_nlte(params, a_lte, '1D')
        row = dict(reference=name, a_lte=round(a_lte, 3), chi2r=round(r['red_chi2'], 1),
                   a_nlte_3d=round(a_lte + d3, 3), a_nlte_1d=round(a_lte + d1, 3))
        rows.append(row)
        print(f"  {name:14} A(O)_LTE {row['a_lte']}  χ²ᵣ {row['chi2r']:5}  "
              f"A(O)_3D-NLTE {row['a_nlte_3d']}  (1D {row['a_nlte_1d']})  vs our Sun {SOLAR_OURS_O}")
    a3 = [r['a_nlte_3d'] for r in rows]
    spread = round(float(np.max(a3) - np.min(a3)), 3)
    kpno = next(r for r in rows if r['reference'] == 'kurucz_kpno')['a_nlte_3d']
    iag = np.mean([r['a_nlte_3d'] for r in rows if 'iag' in r['reference']])
    wallace = next(r for r in rows if r['reference'] == 'wallace_kpno')['a_nlte_3d']
    print(f"\n  per-reference A(O)_3D-NLTE spread = {spread} dex (range {min(a3)}–{max(a3)})")
    print(f"  Kurucz-KPNO {kpno}  |  IAG mean {iag:.3f}  |  Wallace-KPNO {wallace}")
    iag_vs_kurucz = round(float(iag) - kpno, 3)
    wallace_side = 'KPNO-family (≈Kurucz)' if abs(wallace - kpno) < abs(wallace - iag) else 'Kurucz-reduction-specific (≈IAG)'
    verdict = ('continuum NOT reference-driven (all agree within ~0.03) — RYA-483 σ intrinsic'
               if spread < 0.05 else
               f'IAG−Kurucz {iag_vs_kurucz:+.3f}; Wallace sides with {wallace_side}')
    print(f"  DIAGNOSTIC: {verdict}")
    return dict(per_reference=rows, spread=spread, iag_minus_kurucz=iag_vs_kurucz,
                wallace_sides_with=wallace_side, verdict=verdict,
                combined_a_lte=round(float(np.mean([r['a_lte'] for r in rows])), 3),
                combined_a_nlte_3d=round(float(np.mean(a3)), 3),
                combined_a_nlte_1d=round(float(np.mean([r['a_nlte_1d'] for r in rows])), 3))


def measure_sulphur(ll, iso, sab, codes):
    _rule(); print("  SULPHUR — S I multiplet-8 6748.68/6757.15 (HARPS), Amarsi-2025 NLTE"); _rule()
    rec = get_star_params('solar')
    params = {'teff_K': float(rec['teff']), 'logg': float(rec['logg']),
              'feh': float(rec['feh_ref']), 'vturb_kms': float(rec.get('xi', 1.0))}
    w_nm, fl = cs._load_observed_spectrum('solar')               # solar HARPS
    atm = cs._load_atmosphere(params['teff_K'], params['logg'], params['feh'], params['vturb_kms'])
    _, _, chem = cs._load_synth_resources()                       # cached; for the S atom code
    codes = cs._atom_codes(('S',), chem, sab)                     # S-specific codes
    grid = pd.read_csv(ROOT / 'data' / 'nlte_grids' / 'S_Amarsi2025_PySME.csv')
    per = []
    for wl, hw in S_LINES:
        r = cs._fit_element(w_nm, fl, atm, params, 'S', {'S': 7.2}, codes,
                            ((wl - hw, wl + hw),), False, (cs.HARPS_R, SOLAR_VMAC, SOLAR_VSINI),
                            5.8, 8.2, ll, iso, sab, TMP)
        g = grid[(abs(grid['wave_A'] - wl) < 0.1) & (grid['teff_K'] == 5772)]
        d = float(g['delta_nlte'].iloc[0]) if len(g) else 0.0
        per.append(dict(line=wl, a_lte=round(r['A_X'], 3), chi2r=round(r['red_chi2'], 1),
                        nlte_delta=round(d, 4), a_nlte=round(r['A_X'] + d, 3), status=r['status']))
        print(f"  S I {wl}  A_LTE {per[-1]['a_lte']}  χ²ᵣ {per[-1]['chi2r']}  NLTE {d:+.4f}  "
              f"A_NLTE {per[-1]['a_nlte']}  [{r['status']}]")
    vals = [p['a_nlte'] for p in per if np.isfinite(p['a_nlte'])]
    a = round(float(np.median(vals)), 3) if vals else None
    sigma = round(float(np.std(vals, ddof=1)), 3) if len(vals) > 1 else 0.0
    print(f"  A(S)_NLTE median {a} (σ {sigma}, n={len(vals)})  vs Asplund {ASPLUND['S']} "
          f"({(a or 0)-ASPLUND['S']:+.3f}); banked EW-pool 7.753")
    return dict(per_line=per, a_nlte=a, sigma=sigma, n=len(vals), asplund=ASPLUND['S'])


def main():
    Path(TMP).mkdir(parents=True, exist_ok=True)
    ll, iso, sab, codes = _resources()
    cs._CHEM = None
    swap = swap_test(ll, iso, sab, codes)
    sulf = measure_sulphur(ll, iso, sab, codes)

    _rule(); print("  SOLAR CNO REDO — VERDICT (validate-don't-tune; Asplund = comparison only)"); _rule()
    print(f"  O I 777: per-reference A(O)_3D-NLTE spread {swap['spread']} dex; combined "
          f"3D-NLTE {swap['combined_a_nlte_3d']} | 1D-NLTE {swap['combined_a_nlte_1d']}")
    print(f"    new-vs-371: combined 3D {swap['combined_a_nlte_3d']} vs banked O {BANKED['O']} "
          f"(Δ {swap['combined_a_nlte_3d']-BANKED['O']:+.3f})")
    print(f"  REGIME-MATCHED denominators (RYA-484): solar_A(O)_3D {swap['combined_a_nlte_3d']} "
          f"(Sun's own) | solar_A(O)_1D {swap['combined_a_nlte_1d']} (for warm-star differentials)")
    print(f"  S I (synthesis 6748/6757): A(S)_NLTE {sulf['a_nlte']} (σ {sulf['sigma']}) "
          f"— supersedes the EW-pool 7.753 (curation-owed)")
    print(f"  swap-test diagnostic: {swap['verdict']}")
    out = {'ticket': 'RYA-485/486 solar CNO redo', 'validate_dont_tune': True,
           'banked_371': BANKED, 'asplund_comparison_only': ASPLUND,
           'oi777_swap_test': swap,
           'regime_matched_O': {'solar_A_O_3D': swap['combined_a_nlte_3d'],
                                'solar_A_O_1D': swap['combined_a_nlte_1d'],
                                'note': 'use _1D as the denominator for >6500 K stars (Procyon) to regime-match'},
           'sulphur': sulf,
           'deferred': ['full C/N 3D-1D legs across molecular arms (CH/CN/NH) — atomic C I/O I '
                        'have the Amarsi grid; molecular + N I post-hoc are a separate regime split',
                        'IR O I 844/926 nm (Lever-4) — needs telluric-gated IR solar data reaching the line']}
    RESULTS.mkdir(parents=True, exist_ok=True)
    (RESULTS / 'solar_cno_redo_rya485.json').write_text(json.dumps(out, indent=2, default=str))
    print(f"\n  [out] {RESULTS / 'solar_cno_redo_rya485.json'}")
    return out


if __name__ == '__main__':
    main()
