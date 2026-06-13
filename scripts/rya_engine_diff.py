"""
scripts/rya_engine_diff.py
Per-element diff: full Turbospectrum-synthesis abundance run vs MOOG run.
Separates genuine engine differences from the known NLTE divergences (RYA-235).
Run from repo root: python scripts/rya_engine_diff.py <star>
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from config.constants import PATHS

# RYA-235: these divergences were diagnosed as NLTE, not engine differences.
NLTE_EXPECTED = {('Ca', 'I'), ('Ti', 'I'), ('Cr', 'I')}
DIFF_THRESHOLD = 0.10  # dex

def _load(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Abundance table not found: {path}")
    df = pd.read_csv(path)
    if 'A_X' not in df.columns:
        raise KeyError(f"No 'A_X' column in {path}; columns = {list(df.columns)}")
    return df[['element', 'ion', 'A_X']].copy()

def main(star: str):
    base = Path(str(PATHS['solar_ew'])).parent
    moog  = _load(base / f'{star}_abundances.csv').rename(columns={'A_X': 'A_X_moog'})
    synth = _load(base / f'{star}_abundances_synth.csv').rename(columns={'A_X': 'A_X_synth'})

    m = moog.merge(synth, on=['element', 'ion'], how='outer')
    m['delta_synth_minus_moog'] = (m['A_X_synth'] - m['A_X_moog']).round(3)

    def classify(row):
        if pd.isna(row['A_X_moog']) or pd.isna(row['A_X_synth']):
            return 'MISSING-ONE-ENGINE'
        if abs(row['delta_synth_minus_moog']) <= DIFF_THRESHOLD:
            return 'agree'
        if (row['element'], row['ion']) in NLTE_EXPECTED:
            return 'NLTE-EXPECTED (RYA-235)'
        return 'ENGINE-DIVERGENCE (investigate)'

    m['verdict'] = m.apply(classify, axis=1)
    m = m.sort_values(['element', 'ion']).reset_index(drop=True)

    print(f"\nEngine diff (synthesis − MOOG)  |  {star}")
    print(m.to_string(index=False))
    n_div = int((m['verdict'] == 'ENGINE-DIVERGENCE (investigate)').sum())
    print(f"\n  Unexplained engine divergences (|Δ|>{DIFF_THRESHOLD} dex, not NLTE-flagged): {n_div}")

    out = base / f'engine_diff_{star}.csv'
    m.to_csv(out, index=False)
    print(f"  Wrote {out}")

if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else 'solar')
