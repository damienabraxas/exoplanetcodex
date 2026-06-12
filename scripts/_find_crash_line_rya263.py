"""
Find which Procyon EW line triggers SPECTRUM's expint crash (RYA-263).
Uses an anchor set of known-working Fe I lines + test candidates.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config.constants import ISPEC_DIR
sys.path.insert(0, str(ISPEC_DIR))
import numpy as np
import pandas as pd
from pipeline.abundances_derive import _load_atmosphere, _ew_to_abundance

STAR = dict(teff_K=6530, logg=3.96, feh=-0.04, vturb_kms=1.66)

# Lines we know work (rows 0-89 in wavelength-sorted order, 4807-5259 Å)
ANCHOR_WAVE_RANGE = (4807.0, 5260.0)


def crashed(ew_subset: pd.DataFrame, atm, params) -> bool:
    """True if SPECTRUM returned all-zero abundances (crash detected)."""
    try:
        lm, norm, xh, _ = _ew_to_abundance(ew_subset, params, atm)
    except Exception:
        return True
    notes = np.array([str(x) for x in lm['note']])
    fe1 = np.where(notes == 'Fe 1')[0]
    if len(fe1) == 0:
        return None  # can't tell
    return int((norm[fe1] != 0.0).sum()) == 0


if __name__ == '__main__':
    ew_df = pd.read_csv('data/processed/procyon_ew.csv')
    if 'blend_flag' in ew_df.columns:
        ew_df = ew_df[ew_df['blend_flag'] == False]
    ew_df = ew_df[(ew_df['ew_mA'] >= 5) & (ew_df['ew_mA'] <= 300)].copy()

    ew_sorted = ew_df.sort_values('wavelength_air_A').reset_index(drop=True)

    atm = _load_atmosphere(teff=6530, logg=3.96, feh=-0.04, vturb=1.66)
    params = {'teff_K': 6530, 'logg': 3.96, 'feh': -0.04, 'vturb_kms': 1.66}

    # Anchor: lines from the 4807-5260 Å range (known to NOT crash alone)
    anchor = ew_sorted[ew_sorted['wavelength_air_A'] < ANCHOR_WAVE_RANGE[1]].copy()
    print(f"Anchor set: {len(anchor)} lines (4807-5260 Å)")

    # Candidate zone: 5260-5397 Å where crash is known to occur
    candidates = ew_sorted[
        (ew_sorted['wavelength_air_A'] >= 5260.0) &
        (ew_sorted['wavelength_air_A'] <= 5400.0)
    ].reset_index(drop=True)
    print(f"Candidate zone: {len(candidates)} lines (5260-5400 Å)")
    print()

    # Test each candidate individually combined with anchor
    crash_lines = []
    for i, row in candidates.iterrows():
        test_set = pd.concat([anchor, pd.DataFrame([row])], ignore_index=True)
        result = crashed(test_set, atm, params)
        elem = row.get('element', '?')
        ion_  = row.get('ion', '?')
        wave  = row['wavelength_air_A']
        ew    = row['ew_mA']
        status = 'CRASH' if result else ('OK' if result is False else 'UNKNOWN')
        print(f"  {elem} {ion_} {wave:.3f} Å  EW={ew:.1f} mÅ  → {status}")
        if result:
            crash_lines.append(row)

    print(f"\n=== Found {len(crash_lines)} crash-triggering line(s) ===")
    for r in crash_lines:
        print(f"  {r.get('element','?')} {r.get('ion','?')} {r['wavelength_air_A']:.5f} Å  EW={r['ew_mA']:.2f} mÅ")

    if crash_lines:
        # Also verify: anchor + all candidates except crash lines = no crash?
        crash_waves = {r['wavelength_air_A'] for r in crash_lines}
        safe_cands = candidates[~candidates['wavelength_air_A'].isin(crash_waves)]
        test_no_crash = pd.concat([anchor, safe_cands], ignore_index=True)
        result_no_crash = crashed(test_no_crash, atm, params)
        print(f"\nVerification (exclude crash lines): crash={result_no_crash}")
