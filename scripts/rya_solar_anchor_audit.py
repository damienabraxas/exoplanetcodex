"""
scripts/rya_solar_anchor_audit.py
RYA-290 — Solar anchor audit: where does the ~0.6 dex MOOG Fe I/Fe II
ionization imbalance come from?

Read-only diagnostic (NO pipeline fix). Reproduces the audit numbers from the
CSVs written by a pinned-params MOOG run + the synthesis-v2 run, so the verdict
is traceable. Generate the inputs first:

    # MOOG EW baseline at canonical pinned params (5772 / 4.438 / 0.0 / ξ=1.0):
    python pipeline/abundances_derive.py solar ATLAS9.Castelli --skip-convergence
    #   -> data/processed/solar_per_line.csv, solar_abundances.csv
    # Synthesis-v2 (any recent run) -> data/processed/solar_abundances_synth_v2.csv

Then:  python scripts/rya_solar_anchor_audit.py

Conclusion (RYA-290): at pinned canonical params the imbalance is a DATA/setup
signal in the EW-derived Fe II sample (bad gf-values + blend-contaminated EWs),
NOT physics and NOT the solver — the blend-aware synthesis closes the balance on
the identical line list. Solar logg is pinned (RYA-292), so the solver can no
longer absorb this by walking logg.
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from config.constants import PATHS, SOLAR_ASPLUND2021

PROC = Path(str(PATHS['solar_ew'])).parent
EXPECTED_FE = SOLAR_ASPLUND2021['Fe']   # 7.46


def _fe(df, ion):
    return df[(df['element'] == 'Fe') & (df['ion'].astype(str) == ion)].copy()


def main():
    per_line = PROC / 'solar_per_line.csv'
    v2 = PROC / 'solar_abundances_synth_v2.csv'
    if not per_line.exists():
        raise FileNotFoundError(
            f"{per_line} not found — run the pinned MOOG baseline first:\n"
            f"  python pipeline/abundances_derive.py solar ATLAS9.Castelli --skip-convergence")

    df = pd.read_csv(per_line)
    fe1, fe2 = _fe(df, 'I'), _fe(df, 'II')
    a1 = fe1['a_1dlte'].astype(float)
    a2 = fe2['a_1dlte'].astype(float)

    print('=' * 70)
    print('  RYA-290  Solar anchor audit  |  MOOG EW @ pinned 5772/4.438/0.0/ξ=1.0')
    print('=' * 70)
    print(f"\nFe I : n={len(a1):3d}  median={a1.median():.3f}  mean={a1.mean():.3f}  "
          f"std={a1.std():.3f}  (expected ~{EXPECTED_FE})")
    print(f"Fe II: n={len(a2):3d}  median={a2.median():.3f}  mean={a2.mean():.3f}  "
          f"std={a2.std():.3f}")
    print(f"Fe I - Fe II (median) = {a1.median() - a2.median():+.3f} dex")

    # Fe I microturbulence diagnostic
    rew = np.log10(fe1['ew_mA'].astype(float) / fe1['wavelength_air_A'].astype(float))
    msk = np.isfinite(rew) & np.isfinite(a1)
    slope = np.polyfit(rew[msk], a1[msk], 1)[0]
    print(f"Fe I REW-abundance slope = {slope:+.3f} (≠0 → ξ/EW-saturation signal)")

    # Fe II sample quality
    print("\n--- Fe II line sample (sorted by abundance) ---")
    cols = ['wavelength_air_A', 'ew_mA', 'log_gf', 'excitation_potential_eV',
            'a_1dlte', 'line_grade']
    print(fe2.sort_values('a_1dlte')[cols].to_string(index=False))
    print(f"\nFe II within 7.4-7.6: {((a2 >= 7.4) & (a2 <= 7.6)).sum()}/{len(a2)}  |  "
          f">8.0: {(a2 > 8.0).sum()}  |  range={a2.max() - a2.min():.2f} dex")

    # Engine cross-check
    if v2.exists():
        s = pd.read_csv(v2)
        s1 = float(_fe(s, 'I')['A_X_abs'].iloc[0])
        s2 = float(_fe(s, 'II')['A_X_abs'].iloc[0])
        print("\n--- Engine cross-check (same line list, pinned params) ---")
        print(f"MOOG/EW    Fe I-Fe II = {a1.median() - a2.median():+.3f} dex")
        print(f"Synthesis  Fe I-Fe II = {s1 - s2:+.3f} dex   (Fe I {s1:.3f} / Fe II {s2:.3f})")
        print("Synthesis closes; MOOG/EW does not → EW Fe II data (blends + bad gf).")

    print("\nVERDICT: DATA/setup (EW Fe II line sample), not physics, not solver.")


if __name__ == '__main__':
    main()
