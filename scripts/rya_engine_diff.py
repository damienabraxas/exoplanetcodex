"""
scripts/rya_engine_diff.py
Per-element diff: full Turbospectrum-synthesis abundance run vs MOOG run.
Separates genuine engine differences from the known NLTE divergences (RYA-235).

Run from repo root:
  python scripts/rya_engine_diff.py <star>          # v1 synthesis vs MOOG
  python scripts/rya_engine_diff.py <star> --v2     # v2 flux-fit: vs MOOG AND vs v1
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

def _classify_vs_moog(row, a, b):
    if pd.isna(row[a]) or pd.isna(row[b]):
        return 'MISSING-ONE-ENGINE'
    if abs(row[b] - row[a]) <= DIFF_THRESHOLD:
        return 'agree'
    if (row['element'], row['ion']) in NLTE_EXPECTED:
        return 'NLTE-EXPECTED (RYA-235)'
    return 'ENGINE-DIVERGENCE (investigate)'

def main(star: str):
    base = Path(str(PATHS['solar_ew'])).parent
    moog  = _load(base / f'{star}_abundances.csv').rename(columns={'A_X': 'A_X_moog'})
    synth = _load(base / f'{star}_abundances_synth.csv').rename(columns={'A_X': 'A_X_synth'})

    m = moog.merge(synth, on=['element', 'ion'], how='outer')
    m['delta_synth_minus_moog'] = (m['A_X_synth'] - m['A_X_moog']).round(3)
    m['verdict'] = m.apply(lambda r: _classify_vs_moog(r, 'A_X_moog', 'A_X_synth'), axis=1)
    m = m.sort_values(['element', 'ion']).reset_index(drop=True)

    print(f"\nEngine diff (synthesis − MOOG)  |  {star}")
    print(m.to_string(index=False))
    n_div = int((m['verdict'] == 'ENGINE-DIVERGENCE (investigate)').sum())
    print(f"\n  Unexplained engine divergences (|Δ|>{DIFF_THRESHOLD} dex, not NLTE-flagged): {n_div}")

    out = base / f'engine_diff_{star}.csv'
    m.to_csv(out, index=False)
    print(f"  Wrote {out}")

def main_v2(star: str):
    """RYA-287: v2 flux-fit emits TWO comparisons — v2-vs-MOOG and v2-vs-v1."""
    base = Path(str(PATHS['solar_ew'])).parent
    moog = _load(base / f'{star}_abundances.csv').rename(columns={'A_X': 'A_X_moog'})
    v2   = _load(base / f'{star}_abundances_synth_v2.csv').rename(columns={'A_X': 'A_X_v2'})

    # ── v2 vs MOOG ───────────────────────────────────────────────────────────
    mm = moog.merge(v2, on=['element', 'ion'], how='outer')
    mm['delta_v2_minus_moog'] = (mm['A_X_v2'] - mm['A_X_moog']).round(3)
    mm['verdict'] = mm.apply(lambda r: _classify_vs_moog(r, 'A_X_moog', 'A_X_v2'), axis=1)
    mm = mm.sort_values(['element', 'ion']).reset_index(drop=True)
    print(f"\n[1] v2-synth vs MOOG  |  {star}")
    print(mm.to_string(index=False))
    n_div = int((mm['verdict'] == 'ENGINE-DIVERGENCE (investigate)').sum())
    print(f"  Unexplained divergences (|Δ|>{DIFF_THRESHOLD}, not NLTE): {n_div}")
    out_m = base / f'engine_diff_{star}_v2.csv'
    mm.to_csv(out_m, index=False)
    print(f"  Wrote {out_m}")

    # ── v2 vs v1 (proof the FLAG 1 double-count is gone) ─────────────────────
    v1_path = base / f'{star}_abundances_synth.csv'
    if not v1_path.exists():
        print(f"\n[2] v2 vs v1 SKIPPED — v1 baseline {v1_path.name} not found")
        return
    v1 = _load(v1_path).rename(columns={'A_X': 'A_X_v1'})
    vv = v1.merge(v2, on=['element', 'ion'], how='outer')
    vv['delta_v2_minus_v1'] = (vv['A_X_v2'] - vv['A_X_v1']).round(3)

    def classify_vv(row):
        if pd.isna(row['A_X_v1']) or pd.isna(row['A_X_v2']):
            return 'MISSING-ONE'
        d = row['delta_v2_minus_v1']
        if d <= -DIFF_THRESHOLD:
            return 'v2-LOWER (blend/wing fixed)'
        if d >= DIFF_THRESHOLD:
            return 'v2-HIGHER (investigate)'
        return 'unchanged'

    vv['verdict'] = vv.apply(classify_vv, axis=1)
    vv = vv.sort_values('delta_v2_minus_v1').reset_index(drop=True)
    print(f"\n[2] v2-synth vs v1-synth  |  {star}  (v2 LOWER on blended/strong lines = FLAG 1/2 fixed)")
    print(vv.to_string(index=False))
    n_lower = int((vv['verdict'] == 'v2-LOWER (blend/wing fixed)').sum())
    print(f"  Elements where v2 < v1 by >{DIFF_THRESHOLD} dex: {n_lower}")
    out_v = base / f'engine_diff_{star}_v2_vs_v1.csv'
    vv.to_csv(out_v, index=False)
    print(f"  Wrote {out_v}")

if __name__ == '__main__':
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    star = args[0] if args else 'solar'
    if '--v2' in sys.argv:
        main_v2(star)
    else:
        main(star)
