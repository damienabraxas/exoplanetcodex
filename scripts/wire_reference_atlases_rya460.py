#!/usr/bin/env python3
"""
scripts/wire_reference_atlases_rya460.py
========================================
RYA-460 — wire the RYA-459 Kitt Peak Solar Flux Atlas into the solar run as a
MEASURABLE source, closing the diagnostics HARPS-VIS (380-690nm) cannot reach.

Part A  register Kitt Peak as a measurable solar spectrum source + validate the
        leg by the overlap cross-check ([O I] 6300 + O I 777 measured from Kitt
        Peak must AGREE with the Phase A HARPS/ESPRESSO values).
Part B  measure the out-of-HARPS diagnostics from Kitt Peak:
        N  -- N I red (7442/7468, 8216/8223, 8680-8718) + NH 3360 + CN violet 3883
              -> solar cross-indicator agreement map. N I NLTE grid is OWED
              (RYA-369) -> measured LTE-flagged, NLTE-OWED, NOT validated.
        P/K/Co/Sc -- P I near-IR multiplet, K I, Co I, Sc II -> reclassify off
              DATA-GAP.

VALIDATE-DON'T-TUNE: 1D-LTE measured by Turbospectrum flux fit (the Phase A
_fit_element engine) on the normalized Kitt Peak windows; cited NLTE applied only
where a vendored grid exists; cross-instrument disagreement FLAGGED, never averaged;
no tuning to Asplund. Blue-edge (NH/CN/Co/Sc) carries the no-true-continuum caveat.

Out: data/audit/cno_synthesis/solar_kittpeak_rya460.json
"""
import json
import sys
import time
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from pipeline.cno_synthesis import (_load_atmosphere, _load_synth_resources,  # noqa: E402
                                    _atom_codes, _fit_element, _ISPEC_SOLAR_ABUND_FILE)
from pipeline import nlte_cno                                                  # noqa: E402
from config.constants import (SOLAR_ASPLUND2021, get_star_params,             # noqa: E402
                              NLTE_CORRECTION_ELEMENTS)
import ispec                                                                   # noqa: E402

KP_DIR = ROOT / 'data' / 'solar_reference' / 'kpno_flux_atlas'
AUDIT = ROOT / 'data' / 'audit' / 'cno_synthesis'
R_KP = 360000.0                       # Kitt Peak FTS is very high res; depth-robust A(X)

# Phase A reference (the already-validated optical legs) for the leg-validation gate.
PHASE_A_OVERLAP = {'OI_6300': {'arm': 'harps', 'raw_1dlte': 8.80, 'reconciled': 8.73},
                   'OI_777':  {'arm': 'espresso', 'raw_1dlte': 8.917, 'reconciled': 8.74}}
AGREE_TOL = 0.10                      # leg validated if KP raw agrees within this

# Diagnostic -> (segment, free element, windows_A, use_molecules, role, note)
DIAGS = [
    # --- Part A: overlap cross-check (leg validator) ---
    ('OI_6300', 'OI_6300', 'O', [(6299.8, 6300.9)], False, 'overlap',
     'forbidden [O I]; Ni I blend (RYA-365/450); cross-check vs HARPS'),
    ('OI_777', 'OI_777_triplet', 'O', [(7771.0, 7772.9), (7773.2, 7775.0), (7774.6, 7776.3)],
     False, 'overlap', 'O I 777 triplet; cross-check vs ESPRESSO (PRIMARY O)'),
    # --- Part B: N (the unblock) ---
    ('NI_7442_7468', 'NI_7442_7468', 'N', [(7441.5, 7443.2), (7467.5, 7469.2)], False, 'N_atomic',
     'N I red high-excitation; RYA-369 PRIMARY N channel'),
    ('NI_8216_8223', 'NI_8216_8223', 'N', [(8215.4, 8217.3), (8222.0, 8224.2)], False, 'N_atomic',
     'N I red; weak (~5% deep)'),
    ('NI_8680_8718', 'NI_8680_8718', 'N', [(8682.5, 8684.3), (8685.2, 8687.1), (8702.3, 8704.2)],
     False, 'N_atomic', 'N I red multiplet (8683/8686/8703)'),
    ('NH_3360', 'NH_3360', 'N', [(3358.0, 3362.0)], True, 'N_molecular',
     'NH A-X UV band head; BLUE-EDGE no-true-continuum (RYA-451/454); C/O-coupled'),
    ('CN_violet_3883', 'CN_violet_3883', 'N', [(3881.5, 3884.5)], True, 'N_molecular',
     'CN B-X violet; BLUE-EDGE; C/O-coupled (CN demoted to cross-check, RYA-369)'),
    # --- Part B: P/K/Co/Sc (were DATA-GAP) ---
    ('KI_7665_7699', 'KI_7665_7699', 'K', [(7698.4, 7699.6)], False, 'datagap',
     'K I resonance; 7665 is in the telluric O2 A-band, 7699 is the clean line'),
    ('CoI_3845', 'CoI_3845', 'Co', [(3845.0, 3846.0)], False, 'datagap',
     'Co I; BLUE-EDGE crowded'),
    ('ScII_4246', 'ScII_4246', 'Sc', [(4246.2, 4247.4)], False, 'datagap',
     'Sc II 4246; BLUE-EDGE; HFS (caveat)'),
    ('PI_10581_10596', 'PI_10581_10596', 'P', [(10580.8, 10582.4), (10596.1, 10597.7)], False,
     'datagap', 'P I near-IR multiplet -- the alternative to FUV P (HST/STIS, RYA-119)'),
]


def load_kp(segment):
    df = pd.read_csv(KP_DIR / f'kpno_{segment}.csv')
    return df['wavelength_air_A'].values / 10.0, df['residual_flux'].values   # air nm, normalized


def cont_snr(flux):
    hi = flux[flux > np.percentile(flux, 90)]
    return float(np.mean(hi) / np.std(hi)) if np.std(hi) > 0 else float('inf')


def oi_nlte(line_A, a_lte, p):
    """Cited Amarsi-2019 3D-NLTE delta for an O I line (the Phase A correction)."""
    lab = nlte_cno.resolve_line('OI', line_A)
    d = nlte_cno.cno_nlte_delta('OI', lab, p['teff_K'], p['logg'], p['feh'], p['vturb_kms'], a_lte)
    nlte_cno.assert_cno_sign('OI', lab, d)
    return float(d), lab


def main():
    rec = get_star_params('solar')
    p = {'teff_K': float(rec['teff']), 'logg': float(rec['logg']),
         'feh': float(rec['feh_ref']), 'vturb_kms': float(rec.get('xi', 1.0))}
    atm = _load_atmosphere(p['teff_K'], p['logg'], p['feh'], p['vturb_kms'])
    ll, iso, chem = _load_synth_resources()
    sab = ispec.read_solar_abundances(_ISPEC_SOLAR_ABUND_FILE)
    br = (R_KP, float(rec.get('vmac', 3.8) or 3.8), float(rec.get('vsini', 1.8) or 1.8))

    results = []
    for key, seg, el, windows, use_mol, role, note in DIAGS:
        w, f = load_kp(seg)
        snr = cont_snr(f)
        blue_edge = bool(w.min() * 10 < 4300)      # < 430 nm = blanketed blue
        # state + codes: CNO molecular equilibrium needs C,N,O; add the free element + Ni
        need = {'C', 'N', 'O', 'Ni', el}
        codes = _atom_codes(tuple(sorted(need)), chem, sab)
        state = {e: float(SOLAR_ASPLUND2021[e]) for e in ('C', 'N', 'O')}
        state['Ni'] = 6.20
        state[el] = float(SOLAR_ASPLUND2021.get(el, 5.0))
        asp = float(SOLAR_ASPLUND2021.get(el, np.nan))
        t = time.time()
        r = _fit_element(w, f, atm, p, el, state, codes, windows, use_mol, br,
                         asp - 1.2, asp + 1.2, ll, iso, sab, '/tmp/ispec_kp')
        dt = time.time() - t
        a_lte = r['A_X']
        rec_out = {
            'key': key, 'element': el, 'segment': seg, 'role': role,
            'provenance': 'kittpeak-measured', 'note': note,
            'cont_snr': round(snr, 0), 'blue_edge': blue_edge,
            'a_1dlte': round(a_lte, 3) if np.isfinite(a_lte) else None,
            'red_chi2': round(float(r['red_chi2']), 1) if np.isfinite(r['red_chi2']) else None,
            'status': r['status'], 'fit_s': round(dt, 0),
            'asplund2021': asp,
        }
        # cited NLTE: O via Amarsi grid; everything else has no wired grid here
        if el == 'O' and np.isfinite(a_lte):
            line_A = 6300.30 if key == 'OI_6300' else 7773.0
            d, lab = oi_nlte(line_A, a_lte, p)
            rec_out.update({'nlte_kind': 'amarsi2019_grid', 'nlte_delta': round(d, 3),
                            'nlte_label': lab, 'a_corr': round(a_lte + d, 3),
                            'nlte_status': 'NLTE-applied'})
        elif el == 'N':
            rec_out.update({'nlte_kind': None, 'nlte_delta': 0.0,
                            'a_corr': round(a_lte, 3) if np.isfinite(a_lte) else None,
                            'nlte_status': 'NLTE-OWED (N I grid, RYA-369)'})
        else:
            wired = el in NLTE_CORRECTION_ELEMENTS
            grid_avail = (ROOT / 'data' / 'nlte_grids' / f'{el}_Amarsi2020_PySME.csv').exists()
            rec_out.update({
                'nlte_kind': None, 'nlte_delta': 0.0,
                'a_corr': round(a_lte, 3) if np.isfinite(a_lte) else None,
                'nlte_status': ('NLTE-wired-not-applied-here' if wired else
                                ('NLTE-grid-available-not-wired' if grid_avail else
                                 'no-NLTE-grid (LTE-flagged)'))})
        delta_asp = (round(rec_out.get('a_corr') - asp, 3)
                     if rec_out.get('a_corr') is not None and np.isfinite(asp) else None)
        rec_out['delta_vs_asplund'] = delta_asp
        results.append(rec_out)
        print(f"  {key:16s} A(LTE)={rec_out['a_1dlte']}  a_corr={rec_out.get('a_corr')}  "
              f"chi2r={rec_out['red_chi2']}  SNR={rec_out['cont_snr']:.0f}  "
              f"blue={blue_edge}  {rec_out['nlte_status']}  ({dt:.0f}s)")

    # --- Part A verdict: overlap cross-check ---
    overlap = {}
    by_key = {r['key']: r for r in results}
    for k, ref in PHASE_A_OVERLAP.items():
        r = by_key.get(k)
        if r and r['a_1dlte'] is not None:
            d = round(r['a_1dlte'] - ref['raw_1dlte'], 3)
            overlap[k] = {'kittpeak_raw': r['a_1dlte'], 'phase_a_arm': ref['arm'],
                          'phase_a_raw': ref['raw_1dlte'], 'delta': d,
                          'agree': bool(abs(d) <= AGREE_TOL)}
    leg_validated = all(v['agree'] for v in overlap.values()) if overlap else False

    # --- Part B verdict: N solar cross-indicator agreement map ---
    n_inds = [r for r in results if r['element'] == 'N' and r['a_corr'] is not None]
    n_atomic = [r for r in n_inds if r['role'] == 'N_atomic' and not r['blue_edge']]
    if n_atomic:
        a_vals = [r['a_corr'] for r in n_atomic]
        n_map = {'atomic_NI_mean': round(float(np.mean(a_vals)), 3),
                 'atomic_NI_spread': round(float(np.max(a_vals) - np.min(a_vals)), 3),
                 'indicators': {r['key']: r['a_corr'] for r in n_inds},
                 'nlte_status': 'NLTE-OWED (N I grid RYA-369); NOT validated (Teff bracket owed)',
                 'molecular_blue_edge_flagged': [r['key'] for r in n_inds if r['blue_edge']]}
    else:
        n_map = {'note': 'no clean N I atomic indicator'}

    out = {
        'ticket': 'RYA-460', 'generated': date.today().isoformat(),
        'source': 'Kitt Peak Solar Flux Atlas (RYA-459, provenance=measured)',
        'leg_validation': {'overlap': overlap, 'leg_validated': leg_validated,
                           'tol': AGREE_TOL},
        'n_cross_indicator': n_map,
        'measurements': results,
    }
    AUDIT.mkdir(parents=True, exist_ok=True)
    (AUDIT / 'solar_kittpeak_rya460.json').write_text(json.dumps(out, indent=2))
    print(f"\n  LEG VALIDATION (overlap): {'PASS' if leg_validated else 'FLAG'} -> {overlap}")
    print(f"  N cross-indicator: {n_map.get('atomic_NI_mean')} (spread {n_map.get('atomic_NI_spread')})")
    print(f"  wrote {AUDIT / 'solar_kittpeak_rya460.json'}")


if __name__ == '__main__':
    main()
