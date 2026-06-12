"""
scripts/flag_procyon_fe2_blends_rya268.py
==========================================
Flag Fe II 6229 and 6491 Å as blend_flag=True in procyon_ew.csv.
Physically impossible a_1dlte > 9.0 dex confirms blend contamination
at F-star temperatures (Teff=6530 K). RYA-238 Phase 1 re-run finding.

Related: flag_procyon_blends_rya249.py (Fe I 5247/5762 Å),
         fix_procyon_spectrum_rya263.py (Fe II 5376 Å expint crash).
Run all three scripts to fully reproduce procyon_ew.csv blend flag state.
"""
import pandas as pd
from pathlib import Path

EW_FILE = Path(__file__).parent.parent / 'data' / 'processed' / 'procyon_ew.csv'
MATCH_TOL = 0.5   # Å — targets are rounded (6229.0, 6491.0); actual are 6229.340, 6491.254

df = pd.read_csv(EW_FILE)

targets = [6229.0, 6491.0]  # approximate — match within tolerance
flagged = []

for target_wave in targets:
    mask = (df['element'] == 'Fe') & (df['ion'] == 'II') & \
           (abs(df['wavelength_air_A'] - target_wave) < MATCH_TOL)
    matches = df[mask]
    if len(matches) == 0:
        print(f"WARNING: Fe II {target_wave} Å not found within {MATCH_TOL} Å tolerance")
        continue
    for idx in matches.index:
        prev = df.at[idx, 'blend_flag']
        df.at[idx, 'blend_flag'] = True
        wave = df.at[idx, 'wavelength_air_A']
        ew   = df.at[idx, 'ew_mA']
        print(f"Flagged Fe II {wave:.3f} Å  EW={ew:.2f} mÅ: blend_flag {prev} → True")
        flagged.append(wave)

df.to_csv(EW_FILE, index=False)
print(f"\nTotal lines flagged this run: {len(flagged)}")
print(f"All blend-flagged Fe II in procyon_ew.csv:")
fe2_flagged = df[(df['element'] == 'Fe') & (df['ion'] == 'II') & (df['blend_flag'] == True)]
print(fe2_flagged[['wavelength_air_A', 'ew_mA', 'blend_flag']].to_string())
print(f"\nTotal blend_flag=True in procyon_ew.csv: {int(df['blend_flag'].sum())}")
