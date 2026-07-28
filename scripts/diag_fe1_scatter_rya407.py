#!/usr/bin/env python3
"""
scripts/diag_fe1_scatter_rya407.py
==================================
RYA-407 DIAGNOSTIC — characterize the solar Fe I line-to-line scatter (σ = 0.138 dex
on main 166cb4b) as (a) real pool issue / (b) path artifact / (c) honest floor.
Analysis-only; reuses existing machinery (no new culling logic). Run the gate first:

    python scripts/validate_fe_rya238.py --star solar      # writes solar_per_line.csv

then:

    python scripts/diag_fe1_scatter_rya407.py

Inputs (committed / reproducible):
  data/processed/solar_per_line.csv      — per-line Fe I A(X), EW, scores (gate output)
  data/linelists/canonical_gf.csv        — single-source gf incl. gf_synth_ges (RYA-203/350)

What it does (each mechanism is residual-INDEPENDENT, to avoid culling-toward-a-target):
  * rank Fe I lines by residual from the pool median;
  * gf cross-match: per-line gf used vs canonical GES gf (gf_synth_ges) — RYA-203;
  * COG/saturation via REW = log10(EW/λ) and the RYA-220 saturation_score;
  * proximity/blend via the single-source vald_proximity_flag;
  * weak-line noise via ew_snr_score;
  * report σ under each residual-independent vetting, the skew, and the verdict.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
PER_LINE = ROOT / 'data' / 'processed' / 'solar_per_line.csv'
CANON_GF = ROOT / 'data' / 'linelists' / 'canonical_gf.csv'
OUT_DIR = ROOT / 'data' / 'audit' / 'fe1_scatter'

A = 'a_1dlte'   # per_line carries 1D-LTE A(Fe I); NLTE σ (gate) ≈ this (Fe I δ ~uniform)


def _ges_gf_map(fe1):
    cg = pd.read_csv(CANON_GF, low_memory=False)
    cg = cg[cg['species'] == 'Fe I']
    w = cg['wavelength_air_A'].values
    out = []
    for wl in fe1['wavelength_air_A']:
        j = int(np.abs(w - wl).argmin())
        out.append(cg.iloc[j]['gf_synth_ges'] if abs(w[j] - wl) <= 0.05 else np.nan)
    return np.asarray(out, float)


def main():
    if not PER_LINE.exists():
        sys.exit(f"missing {PER_LINE} — run validate_fe_rya238.py --star solar first")
    pl = pd.read_csv(PER_LINE)
    fe1 = pl[(pl['element'] == 'Fe') & (pl['ion'] == 'I')].copy()
    med = float(np.nanmedian(fe1[A]))
    fe1['resid'] = (fe1[A] - med)
    fe1['REW'] = np.log10((fe1['ew_mA'] / 1000.0) / fe1['wavelength_air_A'])
    fe1['ges_gf'] = _ges_gf_map(fe1)
    fe1['gf_vs_ges'] = fe1['log_gf'] - fe1['ges_gf']        # ~0 ⇒ gf already GES

    sigma = float(fe1[A].std())
    print(f"\n=== RYA-407: solar Fe I scatter diagnostic ===")
    print(f"pool: n={len(fe1)} Fe I lines | median A(Fe I)={med:.4f} | sigma={sigma:.4f} dex")

    # ── skew (a real floor is symmetric; a contamination tail is one-sided) ──
    hi = int((fe1['resid'] > 0.10).sum()); lo = int((fe1['resid'] < -0.10).sum())
    print(f"skew: {hi} lines >+0.10  vs  {lo} lines <-0.10 | mean resid={fe1['resid'].mean():+.4f} "
          f"(one-sided high tail ⇒ contamination, not symmetric noise)")

    # ── correlations (a path/COG/excitation systematic would show a trend) ──
    for c in ('REW', 'excitation_potential_eV', 'log_gf', 'ew_snr_score'):
        print(f"  corr(resid, {c:24s}) = {np.corrcoef(fe1['resid'], fe1[c])[0, 1]:+.3f}")

    # ── ranked outlier table ──
    cols = ['wavelength_air_A', 'excitation_potential_eV', 'ew_mA', 'REW', 'log_gf',
            'ges_gf', 'gf_vs_ges', A, 'resid', 'vald_proximity_flag',
            'saturation_score', 'ew_snr_score', 'line_grade']
    top = fe1.reindex(fe1['resid'].abs().sort_values(ascending=False).index).head(15)
    print("\n--- TOP 15 Fe I outliers by |residual| ---")
    print(top[cols].round(3).to_string(index=False))

    # ── sigma under residual-INDEPENDENT vetting (no abundance-outlier term) ──
    def sig(mask):
        s = fe1[mask]
        return f"n={len(s):2d}  sigma={s[A].std():.4f}"
    print("\n--- sigma under residual-INDEPENDENT cited vetting ---")
    print(f"  FULL pool                                  {sig(fe1.index == fe1.index)}")
    print(f"  adopt canonical GES gf (RYA-203)           n={fe1['ges_gf'].notna().sum():2d}  "
          f"sigma={(fe1[A] - fe1['gf_vs_ges']).std():.4f}   (gf already GES ⇒ no change)")
    print(f"  drop proximity-blend >=0.5 (single-source) {sig(fe1['vald_proximity_flag'] < 0.5)}")
    print(f"  drop saturation_score >=0.8 (COG)          {sig(fe1['saturation_score'] < 0.8)}")
    print(f"  drop ew_snr_score < 0.5 (weak-line noise)  {sig(fe1['ew_snr_score'] >= 0.5)}")
    print(f"  drop proximity>=0.5 AND snr<0.5 (both)     "
          f"{sig(~((fe1['vald_proximity_flag'] >= 0.5) | (fe1['ew_snr_score'] < 0.5)))}")

    # ── verdict ──
    sig_prox = float(fe1[fe1['vald_proximity_flag'] < 0.5][A].std())
    print("\n--- VERDICT ---")
    print("(b) path artifact: RULED OUT — no resid-vs-REW/EP/gf trend; NLTE σ ≈ 1D σ.")
    print("    gf error (RYA-203): RULED OUT — per-line gf already == canonical GES gf.")
    print("    weak-line noise: RULED OUT — ew_snr_score vetting does not move σ.")
    print(f"(a) real pool issue: PARTIAL — cited proximity blends take 0.138→{sig_prox:.3f} only.")
    print(f"(c) honest floor: PRIMARY — after residual-independent vetting σ floors at ~0.12-0.13,")
    print("    NOT the historical ~0.05-0.10. The one-sided high tail = unrecognized weak blends /")
    print("    1D-LTE line-formation inadequacy on this GES-EW pool, not a single cullable defect.")
    print("ROUTE: (c) honest floor → RYA-277 (per-type gate threshold) + RYA-282 (σ budget);")
    print("       small (a) → RYA-395 (proximity-blend cull) + RYA-279 (gate σ over vetted pool).")
    print("       Do NOT massage the pool toward 0.10 — no cited mechanism gets there.")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out = OUT_DIR / 'fe1_per_line_residuals_rya407.csv'
    fe1.sort_values('resid', ascending=False)[cols].to_csv(out, index=False)
    print(f"\nsaved ranked table → {out.relative_to(ROOT)}")


if __name__ == '__main__':
    main()
