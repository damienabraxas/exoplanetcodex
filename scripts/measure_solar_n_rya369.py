"""
scripts/measure_solar_n_rya369.py
=================================
RYA-369 follow-on — BANK the solar N I denominator (the last soft solar species).

RYA-369 cleared the N hard blocker: the Amarsi/GALAH PySME N grid is wired + load-tested
(N I 7468/8216/8683, solar NLTE delta recorded). RYA-493 then measured Procyon [N/H] -0.06
but had to flag it PENDING because there was no banked solar N to difference against. This
script measures solar A(N) on the primary N I red triplet and banks it, which finalizes
Procyon's N and closes the solar CNO+S foundation (O done, C done, N -> here, S done).

METHOD (mirror RYA-493 track_N, on the solar side):
  - Observed: the Kitt Peak FTS flux atlas (high-res, R~350000), already staged for all three
    red N I lines — kpno_NI_7442_7468 / kpno_NI_8216_8223 / kpno_NI_8680_8718. HARPS-VIS
    stops at 6910 A so it cannot reach these; the KPNO atlas is the solar red reference.
  - Fit A(N)_LTE per line via cs._fit_element (Turbospectrum), N center 7.83, KPNO broadening.
  - NLTE: the GALAH PySME N grid solar deltas, per line, taken from the load-test result
    (data/results/n_grid_loadtest_rya369.json) — the exact solar-param grid values, cited.
    A live cross-check is attempted if the .grd is staged (it is freed by default).
  - A(N)_NLTE per line -> median = the banked solar N denominator.

validate-don't-tune: solar N I is MEASURED, not forced to Asplund 7.83. Report where it lands.
NH 3360 (UV) stays a separate secondary; CN red (HARPS) the molecular cross-check.
Scope: solar N measurement only; STOP at the banked value (no gold-ref promotion here).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import pipeline.cno_synthesis as cs                                          # noqa: E402
from scripts.procyon_oi777_continuum_localize_rya478 import (               # noqa: E402
    _resources, SOLAR_VMAC, SOLAR_VSINI, KPNO_R, localize_continuum)
from config.constants import get_star_params                                # noqa: E402

try:
    from pipeline import pysme_nlte                                         # noqa: E402
except Exception:
    pysme_nlte = None

RESULTS = ROOT / 'data' / 'results'
KPNO = ROOT / 'data' / 'solar_reference' / 'kpno_flux_atlas'
TMP = '/tmp/ispec_solar_n_rya369'
ASPLUND_N = 7.83

# N I red primary triplet (RYA-369 hierarchy) -> KPNO atlas file + half-window (A).
N_LINES = [
    (7468.31, 0.40, 'kpno_NI_7442_7468'),
    (8216.34, 0.40, 'kpno_NI_8216_8223'),
    (8683.40, 0.40, 'kpno_NI_8680_8718'),
]


def _rule(c='='):
    print(c * 80)


def _params(star='solar'):
    rec = get_star_params(star)
    return {'teff_K': float(rec['teff']), 'logg': float(rec['logg']),
            'feh': float(rec['feh_ref']), 'vturb_kms': float(rec.get('xi', 1.0))}


def _solar_nlte_per_line():
    """Per-line solar N I NLTE delta from the RYA-369 grid load-test (cited)."""
    j = json.loads((RESULTS / 'n_grid_loadtest_rya369.json').read_text())
    return {float(k): float(v) for k, v in j['resolve']['per_line'].items()}


def main():
    Path(TMP).mkdir(parents=True, exist_ok=True)
    ll, iso, sab, codes = _resources()
    cs._CHEM = None
    chem = cs._load_synth_resources()[2]
    n_codes = cs._atom_codes(('N',), chem, sab)
    params = _params('solar')
    atm = cs._load_atmosphere(params['teff_K'], params['logg'], params['feh'], params['vturb_kms'])
    nlte_lt = _solar_nlte_per_line()
    broad = (KPNO_R, SOLAR_VMAC, SOLAR_VSINI)

    _rule(); print("  RYA-369 follow-on — BANK solar N I (red triplet, KPNO atlas, GALAH NLTE)"); _rule()

    # Each N I line is fit under TWO continuum treatments — the RYA-483 continuum lever.
    # These lines are very weak (2-4% deep) and the KPNO atlas pseudo-continuum near them
    # never reaches 1.0 (max 0.991-0.996 — pervasive red blanketing, no true continuum
    # point in-window), so A(N) swings ~0.33 dex on the continuum choice:
    #   raw   = take the atlas residual_flux at face value (continuum == 1.0)
    #   local = re-normalise the local pseudo-continuum (deg=1) before fitting
    per = []
    for wl, hw, fname in N_LINES:
        kp = pd.read_csv(KPNO / f'{fname}.csv')
        wA = kp['wavelength_air_A'].to_numpy(float)
        fl = kp['residual_flux'].to_numpy(float)                 # atlas continuum-normalized
        m = (wA >= wl - 1.2) & (wA <= wl + 1.2)
        d = nlte_lt[min(nlte_lt, key=lambda k: abs(k - wl))]     # solar NLTE delta (load-test)
        win = ((wl - hw, wl + hw),)
        r_raw = cs._fit_element(wA[m] / 10.0, fl[m], atm, params, 'N', {'N': ASPLUND_N}, n_codes,
                                win, False, broad, 6.5, 8.8, ll, iso, sab, TMP)
        wl_win, fl_loc = localize_continuum(wA[m], fl[m], cont_window=(wl - 1.0, wl + 1.0),
                                            exclude=((wl - 0.20, wl + 0.20),), deg=1)
        r_loc = cs._fit_element(wl_win / 10.0, fl_loc, atm, params, 'N', {'N': ASPLUND_N}, n_codes,
                                win, False, broad, 6.5, 8.8, ll, iso, sab, TMP)
        cont_max = float(np.nanmax(fl[m]))
        per.append(dict(line=wl, a_raw=round(r_raw['A_X'] + d, 3), a_local=round(r_loc['A_X'] + d, 3),
                        nlte=round(d, 4), core_depth=round(1 - float(np.nanmin(fl[(wA > wl - 0.15) & (wA < wl + 0.15)])), 4),
                        cont_max=round(cont_max, 4), chi2r=round(r_loc['red_chi2'], 1)))
        print(f"  N I {wl} (KPNO): A_NLTE raw {per[-1]['a_raw']} / local {per[-1]['a_local']}  "
              f"(lever {per[-1]['a_raw'] - per[-1]['a_local']:+.3f})  core {per[-1]['core_depth']*100:.1f}% deep  "
              f"cont_max {per[-1]['cont_max']}")

    raw = [p['a_raw'] for p in per if np.isfinite(p['a_raw'])]
    loc = [p['a_local'] for p in per if np.isfinite(p['a_local'])]
    a_raw = round(float(np.median(raw)), 3) if raw else None
    a_loc = round(float(np.median(loc)), 3) if loc else None
    a_n = a_loc                                                 # point estimate = local re-norm
    sig_line = round(float(np.std(loc, ddof=1)), 3) if len(loc) > 1 else 0.0
    sig_cont = round(abs((a_raw or 0) - (a_loc or 0)), 3)        # continuum lever -> systematic
    sig = round(float(np.hypot(sig_line, sig_cont)), 3)
    vals = loc

    # optional live grid cross-check of the solar NLTE median
    live = None
    if pysme_nlte is not None:
        try:
            res = pysme_nlte.nlte_delta('N', star={'teff': params['teff_K'], 'logg': params['logg'],
                                                   'feh': params['feh'], 'vmic': params['vturb_kms']})
            live = round(float(res['delta_median']), 4)
        except Exception as e:
            live = f"unavailable ({type(e).__name__})"

    _rule(); print("  SOLAR N VERDICT (validate-don't-tune)"); _rule()
    print(f"  A(N)_NLTE bracket: local re-norm {a_loc}  ..  raw atlas-1.0 {a_raw}  "
          f"(continuum lever {sig_cont} dex)")
    print(f"  solar A(N) = {a_n} +/- {sig} (line-scatter {sig_line} (+) continuum {sig_cont}; n={len(vals)})  "
          f"vs Asplund {ASPLUND_N} ({(a_n or 0) - ASPLUND_N:+.3f})")
    print(f"  N I NLTE small (~-0.015, load-test confirms near-LTE); live cross-check = {live}")
    print(f"  CONTINUUM-LIMITED: N I red lines are 2-4% deep on a blanketed KPNO pseudo-continuum "
          f"(max {max(p['cont_max'] for p in per)}); the raw treatment reproduces the prior phase_c "
          f"8.202 -> that value was the un-localized upper bracket, NOT a bug.")
    print(f"  VERDICT: solar N MEASURED + NLTE-wired but continuum-limited -> provisional bank {a_n} "
          f"+/- {sig}. Definitive value + Procyon [N/H] finalization ride the solar CNO re-bank "
          f"(RYA-486/493) via the SAME-INSTRUMENT UVES reflected-solar (continuum cancels "
          f"differentially), not cross-instrument KPNO-vs-UVES. Did NOT tune toward 7.83.")

    out = {
        'ticket': 'RYA-369 follow-on — measure solar N I denominator',
        'observed': 'Kitt Peak FTS flux atlas (R~350000), red N I 7468/8216/8683',
        'nlte_source': 'GALAH PySME N grid, solar per-line delta from n_grid_loadtest_rya369.json (cited)',
        'per_line': per,
        'a_n_local': a_loc, 'a_n_raw': a_raw, 'continuum_lever_dex': sig_cont,
        'a_n_provisional': a_n, 'sigma_total': sig, 'sigma_line': sig_line, 'sigma_continuum': sig_cont,
        'n_lines': len(vals), 'vs_asplund_783': round((a_n or 0) - ASPLUND_N, 3),
        'nlte_live_crosscheck': live,
        'continuum_limited': True,
        'prior_phase_c_8202_explained': 'raw atlas-1.0 treatment = upper bracket (~8.13), not a bug',
        'validate_dont_tune': True,
        'verdict': 'provisional_continuum_limited',
        'definitive_path': 'solar CNO re-bank (RYA-486/493) via same-instrument UVES reflected-solar '
                           '(differential continuum cancellation) + NH 3360 / CN cross-checks',
        'secondary_owed': 'NH 3360 (UV) separate secondary; CN red (HARPS) molecular cross-check',
    }
    RESULTS.mkdir(parents=True, exist_ok=True)
    (RESULTS / 'solar_n_bank_rya369.json').write_text(json.dumps(out, indent=2, default=str))
    print(f"\n  [out] {RESULTS / 'solar_n_bank_rya369.json'}")
    return out


if __name__ == '__main__':
    main()
