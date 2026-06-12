"""
scripts/diagnose_fe2_matching.py
==================================
RYA-247: Diagnose why abundances_derive.py loads only 3 Fe II lines.

Compare:
  1. Fe II EWs in solar_ew.csv
  2. Fe II regions in turbospectrum_synth_good_for_params_codex_extended.txt (FE file)
  3. Fe II regions in limited_but_with_missing_elements_..._all_extended.txt (ALL file)
  4. What _build_ispec_line_regions actually matches for Fe II

Prints full matching table for each Fe II EW: nearest GES wavelength, distance,
and whether it would be matched in the ALL file vs codex_extended (FE file).
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from config.constants import ISPEC_DIR, PATHS

sys.path.insert(0, str(ISPEC_DIR))
import ispec  # noqa: E402

WAVE_TOL = 0.15  # Å — default tolerance in _build_ispec_line_regions

FE_FILE  = str(ISPEC_DIR / 'input' / 'regions' / '42000_VALD' /
               'turbospectrum_synth_good_for_params_codex_extended.txt')
ALL_FILE = str(ISPEC_DIR / 'input' / 'regions' / '42000_VALD' /
               'limited_but_with_missing_elements_turbospectrum_synth_good_for_abundances_all_extended.txt')


def _get_fe2_wavs(path: str) -> np.ndarray:
    regions = ispec.read_line_regions(path)
    notes   = np.array([str(n).strip() for n in regions['note']])
    return regions[notes == 'Fe 2']['wave_A'].astype(float)


def main():
    # ── 1. Fe II EWs in solar_ew.csv ─────────────────────────────────────────
    solar_ew = pd.read_csv(str(PATHS['solar_ew']))
    solar_ew = solar_ew[(solar_ew['ew_mA'] > 0) & solar_ew['ew_mA'].notna()].copy()
    fe2_ew   = solar_ew[(solar_ew['element'] == 'Fe') & (solar_ew['ion'] == 'II')].copy()
    fe2_ew_clean = fe2_ew[
        (~fe2_ew['blend_flag'].fillna(False)) &
        (fe2_ew['ew_mA'] >= 5) & (fe2_ew['ew_mA'] <= 300)
    ]

    print(f"Fe II EWs in solar_ew.csv:                  {len(fe2_ew)} total")
    print(f"Fe II EWs passing [5,300] mÅ + no blend:    {len(fe2_ew_clean)}")

    # ── 2. Fe II GES entries in each regions file ─────────────────────────────
    fe_file_fe2  = _get_fe2_wavs(FE_FILE)
    all_file_fe2 = _get_fe2_wavs(ALL_FILE)

    print(f"Fe II entries in codex_extended (FE file):  {len(fe_file_fe2)}")
    print(f"Fe II entries in all_extended (ALL file):   {len(all_file_fe2)}")

    # ── 3. Full matching table ─────────────────────────────────────────────────
    print(f"\n{'='*80}")
    print(f"MATCHING TABLE: solar_ew.csv Fe II vs GES files (tol={WAVE_TOL} Å)")
    print(f"{'='*80}")
    print(f"{'EW wav (Å)':>12}  {'EW (mÅ)':>8}  "
          f"{'Nearest ALL (Å)':>16}  {'Δ_ALL (Å)':>10}  {'Match_ALL':>10}  "
          f"{'Nearest FE (Å)':>15}  {'Δ_FE (Å)':>9}  {'Match_FE':>9}")
    print(f"{'-'*110}")

    matched_all = 0
    matched_fe  = 0

    for _, row in fe2_ew_clean.sort_values('wavelength_air_A').iterrows():
        wav = float(row['wavelength_air_A'])
        ew  = float(row['ew_mA'])

        # Nearest in ALL file
        if len(all_file_fe2) > 0:
            d_all     = np.abs(all_file_fe2 - wav)
            near_all  = float(all_file_fe2[np.argmin(d_all)])
            delta_all = float(d_all.min())
            hit_all   = delta_all <= WAVE_TOL
        else:
            near_all, delta_all, hit_all = np.nan, np.nan, False

        # Nearest in FE (codex_extended) file
        if len(fe_file_fe2) > 0:
            d_fe     = np.abs(fe_file_fe2 - wav)
            near_fe  = float(fe_file_fe2[np.argmin(d_fe)])
            delta_fe = float(d_fe.min())
            hit_fe   = delta_fe <= WAVE_TOL
        else:
            near_fe, delta_fe, hit_fe = np.nan, np.nan, False

        if hit_all:
            matched_all += 1
        if hit_fe:
            matched_fe += 1

        print(f"{wav:12.4f}  {ew:8.2f}  "
              f"{near_all:16.4f}  {delta_all:10.4f}  {'YES' if hit_all else 'no':>10}  "
              f"{near_fe:15.4f}  {delta_fe:9.4f}  {'YES' if hit_fe else 'no':>9}")

    print(f"{'-'*110}")
    print(f"\nSummary:")
    print(f"  Fe II EWs matched by ALL file (used by abundance derivation): {matched_all}")
    print(f"  Fe II EWs matched by FE file (codex_extended, unused):        {matched_fe}")
    print(f"\n  GAP = {matched_fe - matched_all} lines matched in FE but NOT in ALL")

    # ── 4. Which lines are in FE but not ALL? ──────────────────────────────────
    print(f"\n{'='*80}")
    print("Fe II wavelengths in codex_extended (FE) but missing from all_extended (ALL):")
    print(f"{'='*80}")
    for w in sorted(fe_file_fe2):
        d = np.abs(all_file_fe2 - w).min() if len(all_file_fe2) > 0 else 999
        if d > WAVE_TOL:
            print(f"  {w:.4f} Å  ← MISSING from ALL file")
        else:
            print(f"  {w:.4f} Å    (present in ALL file)")


if __name__ == '__main__':
    main()
