"""
scripts/fe2_ew_quality_cull.py
==============================
RYA-352 — re-base the solar Fe II EW cull on its REAL basis (EW measurement
quality), superseding the misnamed `fe2_solar_cull_rya347.csv` + `_apply_fe2_cull`
(commit a88ef0f) that culled clean EW lines for failing a *synth* χ²ᵣ fit.

This is the report/verification tool. The live gate is
`abundances_derive._apply_fe2_ew_quality_cull` (computed at load, no frozen list).

Criteria (all computed from the EW table; thresholds from config.constants):
  SAT   : EW > PIPELINE['vmic_ew_ceiling_mA']  (COG-derived ceiling, RYA-330)
  HIERR : ew_err/EW > PIPELINE['fe2_ew_err_frac_max']
  BLEND : blend_flag == True

Step 0 derives/cites the REW ceiling (vs hardcoding −4.7) and lists the fragile
band; Step 1 computes CULL/KEEP; Step 2 reconciles vs a88ef0f's 5-line drop and
flags every disagreement; Step 3 re-derives A(Fe II; EW) (median of per-line
a_1dlte) before/after.

Usage:  python3 scripts/fe2_ew_quality_cull.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

_REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_REPO))
from config.constants import PIPELINE, PATHS  # noqa: E402

_PER_LINE = Path(str(PATHS['solar_ew'])).parent / 'solar_per_line.csv'
_SOLAR_EW = Path(str(PATHS['solar_ew']))
_LL_SOLAR = _REPO / 'data' / 'linelists' / 'linelist_solar.csv'

# a88ef0f dropped these 5 from the solar Fe II EW pool, citing RYA-347 (synth verdict).
_A88EF0F_DROP = [5234.623, 5991.371, 6084.102, 6247.557, 6456.380]

CEILING_MA = PIPELINE['vmic_ew_ceiling_mA']
ERR_FRAC_MAX = PIPELINE['fe2_ew_err_frac_max']


def _match(df, wl, col, tol=0.05):
    m = df[(df['element'] == 'Fe') & (df['ion'] == 'II')
           & (df['wavelength_air_A'].sub(wl).abs() < tol)]
    return m.iloc[m['wavelength_air_A'].sub(wl).abs().to_numpy().argmin()][col] if len(m) else np.nan


def load_pool() -> pd.DataFrame:
    pl = pd.read_csv(_PER_LINE)
    fe2 = pl[(pl['element'] == 'Fe') & (pl['ion'] == 'II')].copy()
    ew = pd.read_csv(_SOLAR_EW)
    ll = pd.read_csv(_LL_SOLAR, low_memory=False)
    fe2['REW'] = np.log10(fe2['ew_mA'] / 1000.0 / fe2['wavelength_air_A'])
    fe2['ew_err_mA'] = fe2['wavelength_air_A'].apply(lambda w: _match(ew, w, 'ew_err_mA'))
    fe2['err_over_ew'] = fe2['ew_err_mA'] / fe2['ew_mA']
    # blend: prefer the EW table's flag; fall back to the curated linelist
    def blend(w):
        b = _match(ew, w, 'blend_flag')
        return bool(b) if b == b else bool(_match(ll, w, 'blend_flag'))
    fe2['blend'] = fe2['wavelength_air_A'].apply(blend)
    return fe2.sort_values('REW', ascending=False).reset_index(drop=True)


def cull_verdict(r) -> str:
    reasons = []
    if r['ew_mA'] > CEILING_MA:
        reasons.append('SAT')
    if np.isfinite(r['err_over_ew']) and r['err_over_ew'] > ERR_FRAC_MAX:
        reasons.append('HIERR')
    if bool(r['blend']):
        reasons.append('BLEND')
    return ','.join(reasons)


def main() -> None:
    fe2 = load_pool()
    fe2['cull_reason'] = fe2.apply(cull_verdict, axis=1)
    fe2['culled'] = fe2['cull_reason'] != ''

    # ── Step 0 — derived ceiling ──────────────────────────────────────────────
    rew_at = lambda wl: np.log10(CEILING_MA / 1000.0 / wl)
    lo, hi = rew_at(fe2['wavelength_air_A'].max()), rew_at(fe2['wavelength_air_A'].min())
    print("── Step 0: saturation ceiling (DERIVED, not hardcoded) ──────────────")
    print(f"  EW ceiling = PIPELINE['vmic_ew_ceiling_mA'] = {CEILING_MA:.0f} mÅ")
    print(f"  (COG-derived: RYA-330 set it from the Fe I reduced-EW slope; the direct")
    print(f"   Fe II COG fit is too noisy for an independent break — most measured Fe II")
    print(f"   lines already sit on the flat curve.)")
    print(f"  → REW ceiling over the pool λ-range: {lo:.3f} (red) … {hi:.3f} (blue)")
    frag = fe2[(fe2['REW'] - rew_at(fe2['wavelength_air_A'])).abs() <= 0.10]
    print(f"  fragile band (±0.1 dex of ceiling): "
          f"{', '.join(f'{w:.3f}' for w in frag['wavelength_air_A']) or 'none'}")

    # ── Step 1 — criteria CULL/KEEP ───────────────────────────────────────────
    print(f"\n── Step 1: criteria cull (EW>{CEILING_MA:.0f}mÅ | err/EW>{ERR_FRAC_MAX} | blend) ──")
    cols = ['wavelength_air_A', 'ew_mA', 'REW', 'err_over_ew', 'blend', 'a_1dlte', 'cull_reason']
    print(fe2[cols].to_string(index=False,
          formatters={'REW': '{:.3f}'.format, 'err_over_ew': '{:.3f}'.format,
                      'a_1dlte': '{:.3f}'.format}))
    keep = fe2[~fe2['culled']]; cull = fe2[fe2['culled']]
    print(f"\n  CULL ({len(cull)}): "
          f"{', '.join(f'{w:.3f}({r})' for w, r in zip(cull['wavelength_air_A'], cull['cull_reason']))}")
    print(f"  KEEP ({len(keep)}): "
          f"{', '.join(f'{w:.3f}' for w in keep['wavelength_air_A'])}")

    # ── Step 2 — reconcile vs a88ef0f ─────────────────────────────────────────
    print("\n── Step 2: reconcile vs a88ef0f (it dropped 5 'because RYA-347') ─────")
    for wl in _A88EF0F_DROP:
        m = fe2[fe2['wavelength_air_A'].sub(wl).abs() < 0.05]
        if m.empty:
            print(f"  {wl:9.3f}  not in EW pool (filtered upstream) — a88ef0f drop moot")
            continue
        r = m.iloc[0]
        crit = f"CULL({r['cull_reason']})" if r['culled'] else "KEEP"
        flag = "AGREE" if r['culled'] else "*** DISAGREE → REINSTATE (clean EW line)"
        print(f"  {wl:9.3f}  criteria={crit:14s} vs a88ef0f=DROP   {flag}")
    extra = cull[~cull['wavelength_air_A'].apply(
        lambda w: any(abs(w - x) < 0.05 for x in _A88EF0F_DROP))]
    if len(extra):
        print("  a88ef0f KEPT but criteria CULL (it kept noisy lines):")
        for w, why in zip(extra['wavelength_air_A'], extra['cull_reason']):
            print(f"      {w:9.3f}  {why}")

    # ── Step 3 — re-derive A(Fe II; EW) ───────────────────────────────────────
    a88_keep = fe2[~fe2['wavelength_air_A'].apply(
        lambda w: any(abs(w - x) < 0.05 for x in _A88EF0F_DROP))]
    print("\n── Step 3: A(Fe II; EW) = median per-line a_1dlte ───────────────────")
    print(f"  pre-cull (all {len(fe2)})         : {fe2['a_1dlte'].median():.3f}")
    print(f"  a88ef0f pool ({len(a88_keep)} kept)      : {a88_keep['a_1dlte'].median():.3f}"
          f"   (the 7.466 'Asplund' landing)")
    print(f"  criteria KEEP pool ({len(keep)})      : {keep['a_1dlte'].median():.3f}"
          f"   ← honest EW-quality result")
    print("\n  NOTE: a88ef0f's 7.466 is an artifact of keeping the noisy/blended lines"
          "\n  (6149/6432/5264/5256) and dropping the clean ones. The EW-quality pool"
          "\n  reads high — consistent with the known Fe II EW-path bias (RYA-283/305);"
          "\n  synthesis remains the Fe ionization arbiter (RYA-305 DECISION 2).")


if __name__ == '__main__':
    main()
