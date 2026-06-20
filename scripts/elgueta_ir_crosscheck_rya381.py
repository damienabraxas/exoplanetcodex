#!/usr/bin/env python3
"""
scripts/elgueta_ir_crosscheck_rya381.py
=======================================
RYA-381 step 4 — Elgueta/Jofré IR-window cross-check on the assembled solar list.

Confirms the new IR deliveries (now in linelist_solar.csv) cover the Elgueta et al.
near-IR abundance windows and the species the IR arm is meant to deliver. This is a
COVERAGE check only — the gf cross-check itself is RYA-379.

Windows (air Å) and the species Elgueta/Jofré measure there:
  Y  9800–10800   ·  J  11800–13200   ·  H  15000–17500
Key species: P I, S I, Al I, Zr I, plus n-capture (Sr, Y, Zr, Ba, La, Ce, Nd, Sm, Eu).

Flags (does NOT silently pass) any window or species with zero coverage.

Usage:  python scripts/elgueta_ir_crosscheck_rya381.py
"""
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
SOLAR = ROOT / 'data' / 'linelists' / 'linelist_solar.csv'

ELGUETA_WINDOWS = [('Y', 9800.0, 10800.0), ('J', 11800.0, 13200.0), ('H', 15000.0, 17500.0)]
KEY_SPECIES = ['P', 'S', 'Al', 'Zr']
NCAPTURE = ['Sr', 'Y', 'Zr', 'Ba', 'La', 'Ce', 'Nd', 'Sm', 'Eu']


def main():
    df = pd.read_csv(SOLAR, low_memory=False)
    df['wavelength_air_A'] = pd.to_numeric(df['wavelength_air_A'], errors='coerce')
    df = df.dropna(subset=['wavelength_air_A'])

    print("\n========= RYA-381 ELGUETA IR-WINDOW COVERAGE (assembled solar list) =========")
    gaps = []
    for wn, lo, hi in ELGUETA_WINDOWS:
        win = df[(df['wavelength_air_A'] >= lo) & (df['wavelength_air_A'] < hi)]
        elems = sorted(win['element'].dropna().unique())
        print(f"\n  {wn}-window {lo:.0f}–{hi:.0f} Å : {len(win)} lines, {len(elems)} elements")
        if win.empty:
            gaps.append(f"{wn} window EMPTY")
            continue
        # key species
        ks = {el: int((win['element'] == el).sum()) for el in KEY_SPECIES}
        print("    key species : " + ", ".join(f"{el}={n}" for el, n in ks.items()))
        for el, n in ks.items():
            if n == 0:
                gaps.append(f"{wn}: no {el} I")
        # n-capture
        nc = {el: int((win['element'] == el).sum()) for el in NCAPTURE if (win['element'] == el).any()}
        print("    n-capture   : " + (", ".join(f"{el}={n}" for el, n in nc.items()) if nc else "NONE"))
        if not nc:
            gaps.append(f"{wn}: no n-capture species")

    print("\n  ---- GAPS (flagged for RYA-379 gf work / re-extraction) ----")
    if gaps:
        for g in gaps:
            print(f"    ⚠ {g}")
    else:
        print("    none — all Elgueta windows + key species covered")
    # Whole-IR species presence (beyond the windows), for the RYA-371 IR arm context
    ir = df[df['wavelength_air_A'] >= 9500]
    print(f"\n  IR (≥9500 Å) overall: {len(ir)} lines, species present: "
          f"{', '.join(sorted(ir['element'].dropna().unique()))}")
    return 1 if gaps else 0


if __name__ == '__main__':
    sys.exit(main())
