#!/usr/bin/env python3
"""
RYA-545 — Ti I/II ionization-balance gate via the PRODUCTION EW->abundance path (MOOG-with-blends,
iSpec), atmosphere-matched to MARCS. Sharpens the instrument reading RYA-544's gate; does NOT widen
the pool that passes it.

Three PRE-DECLARED refinements (fixed before re-measuring — firewall: gf-provenance / blend-physics /
correct-inversion, never "adjust until it balances"):
  1. Widen the Ti I pool by QUALITY — every Lawler-2013 (LGWSC/2013ApJS) + lab-grade line, incl. the
     strong ones the RYA-544 EW<60 cut dropped; correct ABO/vdW comes from the GES linemask.
  2. Fix the 6745/6746 blend — handled by the pipeline's theoretical-EW blend filter (rejects the
     blend-contaminated line instead of mis-inverting it as a single line).
  3. Invert through the REAL pipeline: iSpec `determine_abundances` (EW_BASELINE_CODE='moog') on a
     MARCS.GES solar model — NOT the single-line PySME COG stand-in. Proper COG + damping + canonical
     single-source gf (apply_to_regions).

Ti II is the UNTOUCHED validator (no provenance filter). Re-gate ONCE: require SEM << 0.05 so the
verdict means something, then |A(Ti I)_NLTE - A(Ti II)| < FE_IONISATION_GATE. ACCEPT / STOP.
Runs on Sirius: ISPEC_DIR=/srv/codex/engines/ispec_src, venv312. NO MERGE.
"""
from __future__ import annotations
import csv
import os
import sys
from pathlib import Path

REPO = str(Path(__file__).resolve().parents[1])
if REPO not in sys.path:
    sys.path.insert(0, REPO)

import numpy as np
import pandas as pd

import pipeline.abundances_derive as ad
from config.constants import FE_IONISATION_GATE
import ispec

# solar node (same as RYA-542/544; the balance is differential so logg cancels)
STAR = {'teff_K': 5772.0, 'logg': 4.44, 'feh': 0.0, 'vturb_kms': 1.0}
DELTA_TI1 = 0.0506     # OUR Mallinson-2024 Ti I NLTE correction (RYA-544, derived not tuned)
DELTA_TI2 = -0.002     # Mallinson-2024 solar Ti II NLTE correction (~0)
MARCS_DIR = '/srv/codex/engines/ispec_src/input/atmospheres/MARCS.GES'

# RYA-545 pre-declared lab-gf provenance KEEP set (primary laboratory transition probabilities).
LAB_REFS = {'LGWSC', '2013ApJS..', '1982MNRAS.', '1983MNRAS.', '1986MNRAS.', 'NWL', 'MFW', 'BLNP', 'SK'}


def _ti1_lab_wls() -> set:
    """Ti I wavelengths whose canonical gf is a primary lab measurement (or NIST-graded)."""
    keep = set()
    with open(os.path.join(REPO, 'data/linelists/canonical_gf.csv')) as f:
        for r in csv.DictReader(f):
            if r['key_z'] == '22' and r['ion'] == '1':
                if r['loggf_reference'] in LAB_REFS or r['nist_grade'].strip():
                    keep.add(round(float(r['wavelength_air_A']), 2))
    return keep


def _build_ew_df() -> pd.DataFrame:
    """Ti I (lab-gf-graded, ALL EW) + Ti II (untouched validator) measured solar EWs."""
    lab = _ti1_lab_wls()
    rows = []
    with open(os.path.join(REPO, 'data/measured/sol_ew_results_v1.csv')) as f:
        for r in csv.DictReader(f):
            if r['element'] != 'Ti' or r['blend_flag'] != 'False':
                continue
            wl = round(float(r['wavelength_air_A']), 2)
            if r['ion'] == 'I' and wl not in lab:
                continue                       # RYA-545 pre-declared: Ti I keeps only lab-gf lines
            rows.append(dict(element='Ti', ion=r['ion'],
                             wavelength_air_A=float(r['wavelength_air_A']),
                             ew_mA=float(r['ew_mA']), ew_err_mA=float(r['ew_err_mA']),
                             blend_flag=False))
    df = pd.DataFrame(rows)
    n1 = (df.ion == 'I').sum(); n2 = (df.ion == 'II').sum()
    print(f"[pool] Ti I (lab-gf, all EW): {n1} lines   Ti II (untouched validator): {n2} lines")
    return df


def _marcs_solar():
    pack = ispec.load_modeled_layers_pack(MARCS_DIR)
    tgt = {'teff': STAR['teff_K'], 'logg': STAR['logg'], 'MH': STAR['feh'], 'alpha': 0.0}
    if not ispec.valid_atmosphere_target(pack, tgt):
        raise SystemExit(f"MARCS.GES invalid at {tgt}")
    return ispec.interpolate_atmosphere_layers(pack, tgt, code='moog')


def main():
    print("=== RYA-545 Ti ionization-balance gate — production MOOG/iSpec inversion, MARCS.GES ===")
    ew_df = _build_ew_df()
    atm = _marcs_solar()
    linemasks, A_X, _xh, _xfe = ad._ew_to_abundance(
        ew_df, STAR, atm, model_grid='MARCS.GES', code='moog', star_id='solar')
    notes = np.array([str(r['note']).strip() for r in linemasks])
    waves = np.array([float(r['wave_A']) for r in linemasks])
    A_X = np.asarray(A_X, dtype=float)

    out = {}
    for ion, note in (('I', 'Ti 1'), ('II', 'Ti 2')):
        m = (notes == note) & np.isfinite(A_X)
        a = A_X[m]; w = waves[m]
        med = float(np.median(a)) if len(a) else float('nan')
        mad = float(np.median(np.abs(a - med))) if len(a) else float('nan')
        sem = float(np.std(a, ddof=1) / np.sqrt(len(a))) if len(a) > 1 else float('nan')
        out[ion] = (med, sem, len(a))
        print(f"\n  Ti {ion}: n={len(a)}  A_LTE median = {med:.3f}  (MAD {mad:.3f}, SEM {sem:.3f})")
        for wl, ab in sorted(zip(w, a)):
            print(f"       {wl:9.2f}  A = {ab:.3f}")

    (a1, sem1, n1), (a2, sem2, n2) = out['I'], out['II']
    a1n, a2n = a1 + DELTA_TI1, a2 + DELTA_TI2
    bal_lte, bal_nlte = a1 - a2, a1n - a2n
    sem_bal = float(np.hypot(sem1, sem2))
    print(f"\n  A(Ti I)_LTE  = {a1:.3f} (n={n1}, SEM {sem1:.3f})   A(Ti II)_LTE = {a2:.3f} (n={n2}, SEM {sem2:.3f})")
    print(f"  LTE  balance = {bal_lte:+.3f}")
    print(f"  A(Ti I)_NLTE = {a1n:.3f} (+{DELTA_TI1:.4f})   A(Ti II)_NLTE = {a2n:.3f} ({DELTA_TI2:+.4f})")
    print(f"  NLTE ionization balance A(Ti I)_NLTE - A(Ti II)_NLTE = {bal_nlte:+.3f}  ± {sem_bal:.3f} (SEM)")
    print(f"  gate: SEM << {FE_IONISATION_GATE} AND |balance| < {FE_IONISATION_GATE}")

    precise = sem_bal < 0.5 * FE_IONISATION_GATE           # SEM << gate: the verdict means something
    balanced = abs(bal_nlte) < FE_IONISATION_GATE
    if not precise:
        print(f"\n  VERDICT: INCONCLUSIVE — balance SEM {sem_bal:.3f} not << gate {FE_IONISATION_GATE}; "
              f"the pool is still too noisy for a meaningful verdict.")
    elif balanced:
        print(f"\n  VERDICT: ACCEPT — with a properly-measured pool (SEM {sem_bal:.3f}), OUR Ti I(NLTE) "
              f"balances OUR Ti II within the gate. Wire the Mallinson-2024 grid.")
    else:
        print(f"\n  VERDICT: STOP — properly-measured pool (SEM {sem_bal:.3f}) still misses balance "
              f"(|{bal_nlte:+.3f}| > {FE_IONISATION_GATE}). REAL finding: Ti I is physics-limited, "
              f"not a gf-pool fix. Documented; no wire.")


if __name__ == '__main__':
    main()
