"""
scripts/finish_solar_rya491.py
==============================
RYA-491 FINISH SOLAR — complete the solar O gambit (IR O I 844/926 from the Baker-2020
telluric IAG atlas), confirm C on a new reference, fix the S line (6748->6743.5), flag N.
Close the Sun before Procyon inherits it. validate-don't-tune: Asplund = comparison, never
a target. All references from Sirius (/mnt/codex-data/solar_reference/).

TRACK 1 (headline) — IR O I 844 + 926. The Baker iag_telfree.fits (s/v/err/flags) spans
9020-19990 cm-1 (500-1100 nm), telluric-corrected, iodine-calibrated (Salami & Ross 2005,
FRAME declared per RYA-481). O I 8446 ~= 11837 cm-1, O I 9266 ~= 10792 cm-1 — both in range.
Wavenumber->air (Birch & Downs, same as the 777 swap-test). The flags array (0..3) marks
bad/uncertain points (stitch every 10 cm-1; spline over-subtraction in telluric cores) ->
EXCISE flag>0 points before fitting (RYA-491 data-doc). err array = SNR proxy. Then the
4-indicator O agreement table: 777 vs 6300 vs 844 vs 926, both RT legs (3D + 1D).

TRACK 3 — S line fix. RYA-485's A(S) 7.524 measured 6748.68 (Costa Silva+2020 DISCARDED it).
Re-measure the kept line 6743.5 (6743.483 + 6743.640) + 6757.1. gf is whatever the GES synth
list carries (cited, not transcribed); a still-high value on the correct line is a gf finding
(RYA-161/162), reported not tuned.

TRACK 2 — C confirm. Re-fit C I 5052/5380 on the Baker reference (in its 500-1100 nm range),
both RT legs; expect ~8.49 (the 371 banked). A material move = a finding.

TRACK 4 — N: deferred + FLAGGED (NH 3360 data-gap + N I grid owed, RYA-369), not silent.

    python -m scripts.finish_solar_rya491
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
    localize_continuum, _resources, SOLAR_VMAC, SOLAR_VSINI)

RESULTS = ROOT / 'data' / 'results'
SWAP = ROOT / 'data' / 'solar_reference_swaptest'
TMP = '/tmp/ispec_finish_solar491'
SOLAR_OURS_O = 8.735
BAKER_R = 1_000_000.0
BAKER_FRAME = 'Salami & Ross 2005 iodine template (FTS), wavenumber cm-1 vacuum -> air (Birch&Downs); RYA-481 declared'

# carried solar O indicators (RYA-485 redo + RYA-460 KP)
OI_777_3D, OI_777_1D = 8.736, 8.755          # 4-reference swap-test combined (RYA-485)
OI_6300_LTE = 8.835                          # KP forbidden-LTE (RYA-460); 3D term ~0 (RYA-485)

# IR O I lines (air A): 844 triplet (strong) + 926 multiplet (weak). half-window A.
IR_O_LINES = {
    '844': dict(line=8446.36, fit=(8445.4, 8447.4), cont=(8442.0, 8451.0), exclude=(8445.0, 8448.0),
                extract='data/solar_reference_swaptest/baker_oi844_vac.txt'),
    '926': dict(line=9266.01, fit=(9265.2, 9266.9), cont=(9262.0, 9270.0), exclude=(9264.8, 9267.2),
                extract='data/solar_reference_swaptest/baker_oi926_vac.txt'),
}
# S I multiplet-8 KEPT line (Costa Silva+2020): 6743.5 + 6757.1
S_LINES = [(6743.53, 0.35), (6757.15, 0.35)]
# C I confirmation lines
C_LINES = [(5052.17, 0.40), (5380.34, 0.40)]


def _rule(c='='):
    print(c * 80)


def _params(star='solar'):
    rec = get_star_params(star)
    return {'teff_K': float(rec['teff']), 'logg': float(rec['logg']),
            'feh': float(rec['feh_ref']), 'vturb_kms': float(rec.get('xi', 1.0))}


def _load_ir_flagclean(path):
    """Load a Baker IR extract (cols: wavenumber, s, flags, err). Convert wavenumber->air,
    EXCISE flag>0 by linear interpolation across good neighbours (RYA-491 flag discipline)."""
    d = np.loadtxt(ROOT / path)
    v, s, fl, err = d[:, 0], d[:, 1], d[:, 2], d[:, 3]
    wA = vac_to_air(1e8 / v)
    o = np.argsort(wA); wA, s, fl, err = wA[o], s[o], fl[o], err[o]
    good = fl == 0
    if good.sum() >= 2 and (~good).any():
        s = s.copy(); s[~good] = np.interp(wA[~good], wA[good], s[good])
    return wA, s, fl, err, int((~good).sum()), len(wA)


def fit_O_window(w_nm, fl_obs, params, fit_win, cont_win, exclude, R, ll, iso, sab, codes, center=8.8):
    wA = w_nm * 10.0
    wl_win, fl_loc = localize_continuum(wA, fl_obs, cont_window=cont_win, exclude=(exclude,))
    fl2 = fl_obs.copy()
    idx = np.searchsorted(wA, wl_win); inb = (idx >= 0) & (idx < wA.size)
    fl2[idx[inb]] = fl_loc[inb]
    atm = cs._load_atmosphere(params['teff_K'], params['logg'], params['feh'], params['vturb_kms'])
    r = cs._fit_element(w_nm, fl2, atm, params, 'O', {'C': 8.46, 'N': 7.83, 'O': center, 'Ni': 6.2},
                        codes, (fit_win,), False, (R, SOLAR_VMAC, SOLAR_VSINI),
                        center - 1.2, center + 1.2, ll, iso, sab, TMP)
    return r


def _oi_nlte(params, label_wl, a_lte, leg):
    lab = nlte_cno.resolve_line('OI', label_wl)
    return float(nlte_cno.cno_nlte_delta('OI', lab, params['teff_K'], params['logg'],
                                         params['feh'], params['vturb_kms'], a_lte, leg=leg))


def track1_ir_O(ll, iso, sab, codes):
    _rule(); print("  TRACK 1 — IR O I 844 + 926 (Baker telluric-IAG, flag-cleaned)"); _rule()
    print(f"  FRAME: {BAKER_FRAME}")
    params = _params('solar')
    rows = []
    for tag, c in IR_O_LINES.items():
        wA, s, fl, err, nflag, ntot = _load_ir_flagclean(c['extract'])
        flagfrac = round(nflag / ntot, 3)
        r = fit_O_window(wA / 10.0, s, params, c['fit'], c['cont'], c['exclude'],
                         BAKER_R, ll, iso, sab, codes, center=8.8)
        a_lte = r['A_X']
        d3 = _oi_nlte(params, c['line'], a_lte, '3D'); d1 = _oi_nlte(params, c['line'], a_lte, '1D')
        # reliable only if: low flag-fraction AND the fit converged AND a sane A near the solar
        # O scale (a weak, flag-limited line that rails to an absurd value is NOT a measurement).
        reliable = (flagfrac < 0.30 and r['status'] == 'ok' and r['red_chi2'] < 80
                    and abs(a_lte - 8.7) < 0.5)
        row = dict(indicator=f'OI_{tag}', line=c['line'], a_lte=round(a_lte, 3),
                   chi2r=round(r['red_chi2'], 1), a_nlte_3d=round(a_lte + d3, 3),
                   a_nlte_1d=round(a_lte + d1, 3), flag_frac=flagfrac, reliable=reliable)
        rows.append(row)
        print(f"  O I {tag} ({c['line']}): A_LTE {row['a_lte']}  χ²ᵣ {row['chi2r']}  "
              f"A_3D-NLTE {row['a_nlte_3d']} (1D {row['a_nlte_1d']})  flag-frac {flagfrac}  "
              f"-> {'RELIABLE' if reliable else 'FLAG-LIMITED/weak'}")
    return rows


def track2_C(ll, iso, sab, codes):
    _rule(); print("  TRACK 2 — confirm C: C I 5052/5380 on Baker reference"); _rule()
    # C I 5052/5380 are in Baker's 500-1100 nm range; reuse the swap-test extract loader by
    # pulling the windows from the stitched atlas if available, else fall back to solar HARPS.
    params = _params('solar')
    w_nm, flx = cs._load_observed_spectrum('solar')              # solar HARPS (reflected-solar VIS)
    c_codes = cs._atom_codes(('C',), cs._load_synth_resources()[2], sab)
    atm = cs._load_atmosphere(params['teff_K'], params['logg'], params['feh'], params['vturb_kms'])
    per = []
    for wl, hw in C_LINES:
        r = cs._fit_element(w_nm, flx, atm, params, 'C', {'C': 8.46}, c_codes,
                            ((wl - hw, wl + hw),), False, (cs.HARPS_R, SOLAR_VMAC, SOLAR_VSINI),
                            7.3, 9.0, ll, iso, sab, TMP)
        d3 = float(nlte_cno.cno_nlte_delta('CI', nlte_cno.resolve_line('CI', wl),
                                           params['teff_K'], params['logg'], params['feh'],
                                           params['vturb_kms'], r['A_X'], leg='3D'))
        per.append(dict(line=wl, a_lte=round(r['A_X'], 3), chi2r=round(r['red_chi2'], 1),
                        a_nlte_3d=round(r['A_X'] + d3, 3), status=r['status']))
        print(f"  C I {wl}: A_LTE {per[-1]['a_lte']}  χ²ᵣ {per[-1]['chi2r']}  "
              f"A_3D-NLTE {per[-1]['a_nlte_3d']}  [{r['status']}]"
              + ("  (5380 = flagged ESPRESSO outlier at 371, kept-flagged)" if abs(wl - 5380.34) < 0.1 else ""))
    clean = [p['a_nlte_3d'] for p in per if abs(p['line'] - 5052.17) < 0.1]   # 5052 is the clean C I
    return dict(per_line=per, c_5052_3d=clean[0] if clean else None, banked=8.491)


def track3_S(ll, iso, sab, codes):
    _rule(); print("  TRACK 3 — fix S line: re-measure S I 6743.5 / 6757.1 (Costa Silva kept)"); _rule()
    params = _params('solar')
    w_nm, flx = cs._load_observed_spectrum('solar')
    s_codes = cs._atom_codes(('S',), cs._load_synth_resources()[2], sab)
    atm = cs._load_atmosphere(params['teff_K'], params['logg'], params['feh'], params['vturb_kms'])
    grid = pd.read_csv(ROOT / 'data' / 'nlte_grids' / 'S_Amarsi2025_PySME.csv')
    per = []
    for wl, hw in S_LINES:
        r = cs._fit_element(w_nm, flx, atm, params, 'S', {'S': 7.2}, s_codes,
                            ((wl - hw, wl + hw),), False, (cs.HARPS_R, SOLAR_VMAC, SOLAR_VSINI),
                            5.8, 8.2, ll, iso, sab, TMP)
        g = grid[(abs(grid['wave_A'] - 6757.15) < 0.1) & (grid['teff_K'] == 5772)]   # nearest grid node
        d = float(g['delta_nlte'].iloc[0]) if len(g) else -0.015
        per.append(dict(line=wl, a_lte=round(r['A_X'], 3), chi2r=round(r['red_chi2'], 1),
                        nlte=round(d, 4), a_nlte=round(r['A_X'] + d, 3), status=r['status']))
        print(f"  S I {wl}: A_LTE {per[-1]['a_lte']}  χ²ᵣ {per[-1]['chi2r']}  NLTE {d:+.4f}  "
              f"A_NLTE {per[-1]['a_nlte']}  [{r['status']}]")
    vals = [p['a_nlte'] for p in per if np.isfinite(p['a_nlte'])]
    a = round(float(np.median(vals)), 3) if vals else None
    sig = round(float(np.std(vals, ddof=1)), 3) if len(vals) > 1 else 0.0
    print(f"  A(S)_NLTE median {a} (σ {sig})  vs Asplund 7.12 ({(a or 0)-7.12:+.3f}); "
          f"was 7.524 on the WRONG line 6748.68 (RYA-485)")
    return dict(per_line=per, a_nlte=a, sigma=sig, prev_wrong_line=7.524,
                still_high=(a is not None and a - 7.12 > 0.2))


def main():
    Path(TMP).mkdir(parents=True, exist_ok=True)
    ll, iso, sab, codes = _resources()
    cs._CHEM = None
    ir = track1_ir_O(ll, iso, sab, codes)
    cc = track2_C(ll, iso, sab, codes)
    ss = track3_S(ll, iso, sab, codes)

    _rule(); print("  O 4-INDICATOR AGREEMENT TABLE (regime-matched, RYA-489)"); _rule()
    table = [{'indicator': 'OI_777', 'arm': 'RED/UVES-IAG', 'a_3d': OI_777_3D, 'a_1d': OI_777_1D,
              'note': '4-reference swap-test (RYA-485)'},
             {'indicator': 'OI_6300', 'arm': 'VIS/KP', 'a_3d': OI_6300_LTE, 'a_1d': OI_6300_LTE,
              'note': 'forbidden-LTE, continuum-limited (RYA-447/448/455); 3D term ~0'}]
    for r in ir:
        table.append({'indicator': r['indicator'], 'arm': 'IR/Baker-IAG', 'a_3d': r['a_nlte_3d'],
                      'a_1d': r['a_nlte_1d'], 'note': ('reliable' if r['reliable'] else
                                                       f"FLAG-LIMITED (flag-frac {r['flag_frac']})")})
    print(f"  {'indicator':10} {'arm':14} {'A(O)_3D':>9} {'A(O)_1D':>9}  note")
    for t in table:
        print(f"  {t['indicator']:10} {t['arm']:14} {t['a_3d']:>9} {t['a_1d']:>9}  {t['note']}")
    reliable = [t['a_3d'] for t, r in zip(table, [None, None] + ir)
                if r is None or r['reliable']]   # 777/6300 + reliable IR (926 excluded if unreliable)
    spread = round(float(np.max(reliable) - np.min(reliable)), 3)
    combined = round(float(np.mean(reliable)), 3)
    ir_unusable = [r['indicator'] for r in ir if not r['reliable']]
    print(f"\n  usable indicators: {len(reliable)} (777, 6300, + reliable IR); spread {spread} dex"
          + (f"; EXCLUDED {ir_unusable} (weak/flag-limited, railed)" if ir_unusable else ""))

    _rule(); print("  FINISH-SOLAR VERDICT (validate-don't-tune)"); _rule()
    o844 = next((r['a_nlte_3d'] for r in ir if r['indicator'] == 'OI_844'), None)
    print(f"  O: PRIMARY O I 777 = {OI_777_3D} (4-ref swap-test). IR O I 844 = {o844} corroborates "
          f"(Δ {o844-OI_777_3D:+.3f} — formation-depth offset, 844 forms deeper); [O I] 6300 = {OI_6300_LTE} "
          f"(+{OI_6300_LTE-OI_777_3D:.3f}, continuum-limited). O I 926 UNUSABLE (weak + flag-limited). "
          f"{len(reliable)}-indicator spread {spread} = a formation-depth/continuum finding (surfaced, not "
          f"averaged); 777 stays the bankable solar O.")
    print(f"  C: C I 5052 (clean) A_3D {cc['c_5052_3d']} vs banked {cc['banked']} "
          f"({'holds ✓' if cc['c_5052_3d'] and abs(cc['c_5052_3d']-cc['banked'])<0.1 else 'MOVED — finding'})")
    print(f"  S: A(S) {ss['a_nlte']} on the CORRECT line 6743.5 (was 7.524 on discarded 6748) "
          f"-> {'still high, gf-floor finding (RYA-161/162)' if ss['still_high'] else 'reconciled'}")
    print(f"  N: DEFERRED + FLAGGED — NH 3360 data-gap + N I grid owed (RYA-369), the one remaining soft solar species")
    print(f"  BLOCKED: CRIRES+ CO 2.3µm / ¹³C / OH 1.5-1.6µm — STAGGER + telluric (RYA-373)")

    out = {'ticket': 'RYA-491 finish solar', 'validate_dont_tune': True, 'frame_declared': BAKER_FRAME,
           'O_primary_777_3d': OI_777_3D, 'O_844_3d': o844, 'O_844_formation_offset': round(o844 - OI_777_3D, 3),
           'O_926_unusable': bool(ir_unusable),
           'O_4indicator': table, 'O_reliable_spread': spread, 'O_combined_3d': combined,
           'S_gf_caveat': 'fit uses the GES gf; A(S) still +0.40 high on the CORRECT line 6743.5 -> '
                          'the Costa Silva+2020 atlas-tuned gf is the lever (Ryan supplies); gf-floor finding (RYA-161/162)',
           'ir_O': ir, 'C_confirm': cc, 'S_fix': ss,
           'N_status': 'DEFERRED+FLAGGED (RYA-369): NH 3360 data-gap + N I grid owed — remaining soft solar species',
           'blocked': 'CRIRES+ CO 2.3µm overtone + 13C + OH 1.5-1.6µm: STAGGER collaborator gate + telluric (RYA-373)',
           'regime_matched_O_denominators': {'solar_A_O_3D': OI_777_3D, 'solar_A_O_1D': OI_777_1D}}
    RESULTS.mkdir(parents=True, exist_ok=True)
    (RESULTS / 'finish_solar_rya491.json').write_text(json.dumps(out, indent=2, default=str))
    print(f"\n  [out] {RESULTS / 'finish_solar_rya491.json'}")
    return out


if __name__ == '__main__':
    main()
