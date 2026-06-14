"""
scripts/rya_engine_diff.py
Per-element diff: full Turbospectrum-synthesis abundance run vs the MOOG EW baseline.

Run from repo root:
  python scripts/rya_engine_diff.py <star>          # v1 synthesis vs MOOG
  python scripts/rya_engine_diff.py <star> --v2     # v2 flux-fit: vs MOOG AND vs v1

RYA-289: the comparison is now strictly ABSOLUTE-vs-ABSOLUTE and scale-aware.
  * The MOOG baseline (<star>_abundances.csv) carries scale='absolute'; its A_X
    column is absolute A(X). The synthesis tables carry the absolute value in
    A_X_abs. The harness asserts both inputs are on the labelled absolute scale
    before joining and refuses a mismatched-scale join — that mismatch (MOOG
    absolute joined against the synth Asplund-2009-differential A_X) produced the
    spurious ~-7.6 dex offsets.
  * The RYA-235 NLTE carve-out is REMOVED. Synthesis (1D LTE) vs MOOG (1D LTE) is
    an LTE-vs-LTE comparison; NLTE is applied downstream of both engines, so an
    NLTE special-case here is a category error that would mask a real future
    engine divergence on Ca I / Ti I / Cr I. The verdict machinery is unchanged.
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from config.constants import PATHS

DIFF_THRESHOLD = 0.10  # dex


def _load(path: Path, out_col: str):
    """
    Load an engine abundance table and return (df[element, ion, out_col], scale).

    The returned out_col always holds the ABSOLUTE A(X). Picking that column is
    scale-aware so we never silently diff a differential column against an
    absolute one (RYA-289):
      * A_X_abs present  → synthesis table; A_X_abs is absolute by construction.
      * else scale label → MOOG baseline; require scale=='absolute', use A_X.
      * else              → ambiguous (bare A_X could be differential) → refuse.
    """
    if not path.exists():
        raise FileNotFoundError(f"Abundance table not found: {path}")
    df = pd.read_csv(path)

    if 'A_X_abs' in df.columns:
        src, scale = 'A_X_abs', 'absolute'
    elif 'scale' in df.columns:
        labels = set(df['scale'].dropna().astype(str).unique())
        if labels - {'absolute'}:
            raise ValueError(
                f"{path.name}: non-absolute scale label(s) {labels - {'absolute'}}; "
                f"refusing to diff (RYA-289 requires absolute-scale inputs).")
        scale = 'absolute'
        if 'A_X' not in df.columns:
            raise KeyError(f"{path.name}: scale labelled but no 'A_X' column "
                           f"(columns = {list(df.columns)})")
        src = 'A_X'
    else:
        raise KeyError(
            f"{path.name}: cannot determine an absolute, scale-labelled abundance "
            f"column. Need either 'A_X_abs' or a 'scale'='absolute' label. This guard "
            f"is the RYA-289 fix — a bare 'A_X' is ambiguous (absolute on the MOOG "
            f"path, Asplund-2009-differential on the synthesis path). columns = "
            f"{list(df.columns)}")

    out = df[['element', 'ion', src]].rename(columns={src: out_col})
    return out, scale


def _assert_same_scale(scale_a: str, label_a: str, scale_b: str, label_b: str):
    """Refuse a mismatched-scale join with a clear error (RYA-289)."""
    if scale_a != scale_b:
        raise ValueError(
            f"Scale mismatch: {label_a} is '{scale_a}' but {label_b} is '{scale_b}'. "
            f"Refusing to diff across scales — joining absolute against differential "
            f"is the bug RYA-289 fixes (it produced the spurious Procyon -7.6 dex).")
    if scale_a != 'absolute':
        raise ValueError(f"Both inputs are '{scale_a}', expected 'absolute' (RYA-289).")


def _classify_vs_moog(row, a, b):
    if pd.isna(row[a]) or pd.isna(row[b]):
        return 'MISSING-ONE-ENGINE'
    if abs(row[b] - row[a]) <= DIFF_THRESHOLD:
        return 'agree'
    # RYA-289: no NLTE carve-out — LTE-vs-LTE divergence is a real finding.
    return 'ENGINE-DIVERGENCE (investigate)'


def main(star: str):
    base = Path(str(PATHS['solar_ew'])).parent
    moog,  s_moog  = _load(base / f'{star}_abundances.csv', 'A_X_moog')
    synth, s_synth = _load(base / f'{star}_abundances_synth.csv', 'A_X_synth')
    _assert_same_scale(s_moog, 'MOOG baseline', s_synth, 'v1 synthesis')

    m = moog.merge(synth, on=['element', 'ion'], how='outer')
    m['delta_synth_minus_moog'] = (m['A_X_synth'] - m['A_X_moog']).round(3)
    m['verdict'] = m.apply(lambda r: _classify_vs_moog(r, 'A_X_moog', 'A_X_synth'), axis=1)
    m = m.sort_values(['element', 'ion']).reset_index(drop=True)

    print(f"\nEngine diff (synthesis - MOOG, both absolute)  |  {star}")
    print(m.to_string(index=False))
    n_div = int((m['verdict'] == 'ENGINE-DIVERGENCE (investigate)').sum())
    print(f"\n  Unexplained engine divergences (|Δ|>{DIFF_THRESHOLD} dex): {n_div}")

    out = base / f'engine_diff_{star}.csv'
    m.to_csv(out, index=False)
    print(f"  Wrote {out}")


def main_v2(star: str):
    """RYA-287: v2 flux-fit emits TWO comparisons — v2-vs-MOOG and v2-vs-v1."""
    base = Path(str(PATHS['solar_ew'])).parent
    moog, s_moog = _load(base / f'{star}_abundances.csv', 'A_X_moog')
    v2,   s_v2   = _load(base / f'{star}_abundances_synth_v2.csv', 'A_X_v2')
    _assert_same_scale(s_moog, 'MOOG baseline', s_v2, 'v2 synthesis')

    # ── v2 vs MOOG ───────────────────────────────────────────────────────────
    mm = moog.merge(v2, on=['element', 'ion'], how='outer')
    mm['delta_v2_minus_moog'] = (mm['A_X_v2'] - mm['A_X_moog']).round(3)
    mm['verdict'] = mm.apply(lambda r: _classify_vs_moog(r, 'A_X_moog', 'A_X_v2'), axis=1)
    mm = mm.sort_values(['element', 'ion']).reset_index(drop=True)
    print(f"\n[1] v2-synth vs MOOG (both absolute)  |  {star}")
    print(mm.to_string(index=False))
    n_div = int((mm['verdict'] == 'ENGINE-DIVERGENCE (investigate)').sum())
    print(f"  Unexplained divergences (|Δ|>{DIFF_THRESHOLD} dex): {n_div}")
    out_m = base / f'engine_diff_{star}_v2.csv'
    mm.to_csv(out_m, index=False)
    print(f"  Wrote {out_m}")

    # ── v2 vs v1 (proof the FLAG 1 double-count is gone) ─────────────────────
    v1_path = base / f'{star}_abundances_synth.csv'
    if not v1_path.exists():
        print(f"\n[2] v2 vs v1 SKIPPED — v1 baseline {v1_path.name} not found")
        return
    v1, s_v1 = _load(v1_path, 'A_X_v1')
    _assert_same_scale(s_v1, 'v1 synthesis', s_v2, 'v2 synthesis')
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
    print(f"\n[2] v2-synth vs v1-synth (both absolute)  |  {star}  "
          f"(v2 LOWER on blended/strong lines = FLAG 1/2 fixed)")
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
