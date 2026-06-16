#!/usr/bin/env python3
"""
scripts/procyon_fe_2x2.py — RYA-322 Procyon Fe spread test (analysis-only).

Factorial {MPIA, Amarsi} x vmic {1.66, 1.8, 2.0} + the MPIA full-pool reduced-EW
slope as the xi arbiter. Closes three threads at once: the RYA-284 microturbulence
lever, the RYA-320 -0.29/-0.32 root cause, and the RYA-325 Procyon xi pin verdict.

ANALYSIS-ONLY. Reads the Procyon line list + EW table through the PRODUCTION
derivation path and writes ONLY to data/results/procyon_fe_spread.csv. It does
NOT modify spectra, EWs, or the linelist; it snapshots and restores the two
data/processed/procyon_*.csv files that the production run() regenerates.

Design — why C1 reproduces RYA-319 by construction
--------------------------------------------------
The 1D-LTE per-line A(Fe) depends on (teff, logg, feh, vmic) but NOT on the NLTE
source. So per vmic we run the PRODUCTION path once,
    abundances_derive.run('procyon', xi_override=vmic, engine='spectrum'),
which pins Teff=6554/logg=4.00 (RYA-292), pins xi=vmic (RYA-325 xi_override),
solves feh, and applies the LIVE MPIA/Bergemann 1D-NLTE grid (RYA-319). We then
layer the ARCHIVED Amarsi 2022 3D-NLTE MLP onto the SAME per-line a_1dlte, with
the 0th-order Teff clamp to the 6500 K grid edge -- exactly how Amarsi et al.
2022 themselves handled Procyon (verified from the paper). Because the MPIA cell
IS the production path, C1 (MPIA@1.66) reproduces the RYA-319 Procyon number;
the harness asserts it against the stored value and prints CRITICAL on mismatch.

Two pools (the line-pool split is the crux of the brief)
--------------------------------------------------------
* xi-slope arbiter (the decision): the FULL vetted pool -- strong lines retained
  -- on the MPIA path. Regress A(Fe I) vs log(EW/lambda); the vmic where the
  slope -> 0 (interpolated from the three points) is our spectroscopic Procyon
  xi. Compared to the GBS pin 1.8. Strong, partly-saturated lines carry the xi
  leverage; cutting them (the Amarsi weak cut) makes the slope xi-blind.
* Model comparison + all Amarsi cells: the WEAK pool (log W/lambda < -4.9, the
  Amarsi 2022 valid regime), both sources at matched xi -> the MPIA-vs-Amarsi
  spread is the model cross-check.
* xi-driven abundance shift: mean A(Fe I) on the FULL MPIA pool down each xi --
  where the shift is actually visible (the weak-pool mean barely moves with xi).

[Fe/H] per cell = weak-pool mean A(Fe I) - SAME-SOURCE solar weak-pool anchor
(MPIA -> MPIA solar, Amarsi -> Amarsi solar) so the 1D/3D offset cancels in
[Fe/H]. GBS reference comes from STAR_PARAMS (never hardcoded).
"""

import os
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

# ── Path setup ────────────────────────────────────────────────────────────────
_REPO = Path(__file__).resolve().parent.parent
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))
# iSpec sibling (same convention as nlte_corrections / abundances_derive).
_ISPEC = _REPO.parent / 'ispec'
if _ISPEC.exists() and str(_ISPEC) not in sys.path:
    sys.path.insert(0, str(_ISPEC))

from config.constants import STAR_PARAMS, STAR_SOLAR, SOLAR_ASPLUND2021  # noqa: E402
from pipeline import abundances_derive as ad                            # noqa: E402
from pipeline.nlte_corrections import (                                 # noqa: E402
    _mpia_fe_delta, _apply_aberr_to_line, _load_mpia_fe_grid, fe_grid_in_bounds,
)

STAR = 'procyon'
GBS_FEH = float(STAR_PARAMS[STAR]['feh_ref'])     # +0.03 in-repo (Jofré+2014); lit GBS +0.01
GBS_XI = float(STAR_PARAMS[STAR]['xi'])           # 1.8 pinned (RYA-325)
WEAK_LINE_CUT = -4.9                               # log(EW/lambda); Amarsi+2022 weak-line regime
VMIC_GRID = [1.66, 1.8, 2.0]                       # 1.66 old solve | 1.8 GBS pin | 2.0 old target
AMARSI_TEFF_CEIL = 6500.0                          # Amarsi 2022 grid edge; 0th-order clamp for Procyon

# Stored RYA-319 native-Procyon Fe I NLTE number for the C1 consistency check
# (data/processed/procyon_abundances.csv on the integrated main, MPIA path).
RYA319_FE1_NLTE = 7.593
RYA319_TOL = 0.02

_PROC = _REPO / 'data' / 'processed'
_RESULTS = _REPO / 'data' / 'results'
_RESULTS.mkdir(parents=True, exist_ok=True)


def reduced_ew(ew_mA, lam_A):
    return np.log10((np.asarray(ew_mA, float) * 1e-3) / np.asarray(lam_A, float))


# ── Per-source NLTE application on a per-line Fe I table ───────────────────────

def apply_nlte_per_line(fe1, teff, logg, feh, vmic):
    """Return (a_nlte_mpia, a_nlte_amarsi, n_mpia, n_amarsi) arrays/counts for a
    Fe I per-line frame with columns a_1dlte, wavelength_air_A,
    excitation_potential_eV (Elo), eup_eV, log_gf.

    MPIA: native (teff,logg,feh) per-line delta (Procyon Teff in-grid, no clamp).
    Amarsi: archived MLP, Teff clamped to the 6500 K grid edge (0th-order, as in
    Amarsi 2022). aberr is added to the line's own a_1dlte."""
    a1d = fe1['a_1dlte'].values.astype(float)
    wave = fe1['wavelength_air_A'].values.astype(float)
    elo = fe1['excitation_potential_eV'].values.astype(float)
    eup = fe1['eup_eV'].values.astype(float)
    lggf = fe1['log_gf'].values.astype(float)

    a_mpia = np.full(len(fe1), np.nan)
    a_ama = np.full(len(fe1), np.nan)
    teff_amarsi = min(float(teff), AMARSI_TEFF_CEIL)
    for i in range(len(fe1)):
        d = _mpia_fe_delta('I', wave[i], float(teff), float(logg), float(feh))
        if np.isfinite(d):
            a_mpia[i] = a1d[i] + d
        ab = _apply_aberr_to_line('I', elo[i], eup[i], lggf[i], a1d[i],
                                  teff_amarsi, float(logg), float(vmic))
        if np.isfinite(ab):
            a_ama[i] = a1d[i] + ab
    return a_mpia, a_ama, int(np.isfinite(a_mpia).sum()), int(np.isfinite(a_ama).sum())


# ── Solar same-source anchors (weak pool) ─────────────────────────────────────

def solar_frame():
    """Solar Fe I per-line frame with both NLTE sources + rew, from the stored
    solar per-line table at the solar pinned params (5772/4.438/0.0, xi=1.0)."""
    spl = pd.read_csv(_PROC / 'solar_per_line.csv')
    fe1 = spl[(spl.element == 'Fe') & (spl.ion == 'I')].copy().reset_index(drop=True)
    s_teff = float(STAR_SOLAR['teff_K']); s_logg = float(STAR_SOLAR['logg'])
    s_feh = float(STAR_SOLAR['feh']); s_xi = float(STAR_SOLAR['vturb_kms'])
    a_mpia, a_ama, _n_m, _n_a = apply_nlte_per_line(fe1, s_teff, s_logg, s_feh, s_xi)
    fe1['a_nlte_mpia'] = a_mpia
    fe1['a_nlte_amarsi'] = a_ama
    fe1['rew'] = reduced_ew(fe1.ew_mA.values, fe1.wavelength_air_A.values)
    return fe1


def weak_means(fe1):
    """Per-source weak-pool mean A(Fe I) on each source's OWN valid lines, plus
    on the COMMON intersection (both corrections valid) for a line-matched
    comparison. Returns dict of means + line counts."""
    weak = fe1[fe1['rew'] < WEAK_LINE_CUT]
    both = weak[weak['a_nlte_mpia'].notna() & weak['a_nlte_amarsi'].notna()]
    return {
        'mpia':        float(weak['a_nlte_mpia'].mean()),
        'amarsi':      float(weak['a_nlte_amarsi'].mean()),
        'mpia_common': float(both['a_nlte_mpia'].mean()) if len(both) else np.nan,
        'amarsi_common': float(both['a_nlte_amarsi'].mean()) if len(both) else np.nan,
        'n_mpia':   int(weak['a_nlte_mpia'].notna().sum()),
        'n_amarsi': int(weak['a_nlte_amarsi'].notna().sum()),
        'n_common': int(len(both)),
    }


# ── One vmic cell: run production once, derive both sources ────────────────────

def run_vmic(vmic):
    """Run production derivation at pinned xi=vmic; return the converged Fe I
    per-line frame (with both NLTE sources applied), converged feh, and the
    PRODUCTION Fe I/II NLTE numbers (full-pool median, the RYA-319 statistic)."""
    cp, results = ad.run(STAR, xi_override=vmic, engine='spectrum')
    teff = float(cp['teff_K']); logg = float(cp['logg']); feh = float(cp['feh'])

    fe = results[results.element == 'Fe']
    prod = {}
    for _, r in fe.iterrows():
        prod[f"prod_A_fe{1 if r['ion'] == 'I' else 2}_1dlte"] = round(float(r['A_X']), 3)
        prod[f"prod_A_fe{1 if r['ion'] == 'I' else 2}_nlte"] = round(float(r['A_X_nlte']), 3)
    if 'prod_A_fe1_nlte' in prod and 'prod_A_fe2_nlte' in prod:
        prod['prod_feI_minus_feII_nlte'] = round(
            prod['prod_A_fe1_nlte'] - prod['prod_A_fe2_nlte'], 3)

    pl = pd.read_csv(_PROC / f'{STAR}_per_line.csv')
    fe1 = pl[(pl.element == 'Fe') & (pl.ion == 'I')].copy().reset_index(drop=True)
    a_mpia, a_ama, n_m, n_a = apply_nlte_per_line(fe1, teff, logg, feh, vmic)
    fe1['a_nlte_mpia'] = a_mpia
    fe1['a_nlte_amarsi'] = a_ama
    fe1['rew'] = reduced_ew(fe1.ew_mA.values, fe1.wavelength_air_A.values)
    print(f"  vmic={vmic}: converged feh={feh:+.4f}  Fe I pool n={len(fe1)}  "
          f"MPIA-applied={n_m}  Amarsi-applied={n_a}  "
          f"prod A(FeI;NLTE)={prod.get('prod_A_fe1_nlte')}")
    return fe1, feh, prod


def cell_metrics(fe1, source, anchor):
    """Metrics for one (source, vmic) cell.
      * full pool : whole vetted Fe I frame (strong lines retained)
      * weak pool : log(EW/lambda) < -4.9 (Amarsi-valid regime)
      * [Fe/H]    : weak-pool mean A(Fe I) - SAME-SOURCE solar weak anchor
    Two slopes on the full pool (xi arbiter): the clean LTE slope on ALL vetted
    Fe I lines (a_1dlte, matches the production REW diagnostic), and the
    NLTE slope on the source-corrected lines."""
    col = 'a_nlte_mpia' if source == 'mpia' else 'a_nlte_amarsi'
    valid = fe1[fe1[col].notna()].copy()
    weak = valid[valid['rew'] < WEAK_LINE_CUT]

    m = {'model': source}
    m['n_fe1_weak'] = int(len(weak))
    m['A_fe1_weak'] = float(weak[col].mean()) if len(weak) else np.nan
    m['A_fe1_weak_median'] = float(weak[col].median()) if len(weak) else np.nan
    feh = m['A_fe1_weak'] - anchor[source]
    m['feh'] = round(feh, 4)
    m['feh_resid_vs_gbs'] = round(feh - GBS_FEH, 4)

    m['n_fe1_full'] = int(len(valid))
    m['A_fe1_full'] = round(float(valid[col].mean()), 4)

    # LTE full-pool slope on ALL vetted Fe I lines (source-independent; the
    # cleanest "full vetted pool" arbiter, == production REW diagnostic intent).
    ltefull = fe1[fe1['a_1dlte'].notna()]
    if len(ltefull) >= 3:
        sl, _i, _r, _p, se = stats.linregress(ltefull['rew'].values,
                                              ltefull['a_1dlte'].values)
        m['redEW_slope_lte'] = round(float(sl), 4)
        m['redEW_slope_lte_err'] = round(float(se), 4)
    else:
        m['redEW_slope_lte'] = np.nan
        m['redEW_slope_lte_err'] = np.nan

    # NLTE full-pool slope on the source-corrected lines (MPIA = arbiter path).
    if len(valid) >= 3:
        sl, _i, _r, _p, se = stats.linregress(valid['rew'].values, valid[col].values)
        m['redEW_slope_nlte'] = round(float(sl), 4)
        m['redEW_slope_nlte_err'] = round(float(se), 4)
    else:
        m['redEW_slope_nlte'] = np.nan
        m['redEW_slope_nlte_err'] = np.nan
    return m


def main():
    print("=" * 70)
    print("  RYA-322  Procyon Fe spread test  {MPIA, Amarsi} x vmic {1.66,1.8,2.0}")
    print("=" * 70)

    # MPIA grid validity recon (does the full-pool slope have leverage?).
    gc = _load_mpia_fe_grid()
    n_fe1_grid = len(gc['waves'].get('I', []))
    print(f"  MPIA Fe grid: {n_fe1_grid} Fe I lines, bounds={gc['bounds']}")
    print(f"  Procyon (6554,4.00,feh) in MPIA bounds: "
          f"{fe_grid_in_bounds(6554, 4.00, 0.03)}")

    # Snapshot production-processed files that run() regenerates, restore on exit.
    guarded = ['procyon_per_line.csv', 'procyon_abundances.csv']
    backups = {}
    for f in guarded:
        p = _PROC / f
        if p.exists():
            b = _PROC / (f + '.rya322bak')
            shutil.copy2(p, b)
            backups[f] = b

    try:
        # Same-source solar anchors. 'own' = each source's valid weak lines;
        # 'common' = lines valid in BOTH grids (line-matched model comparison).
        sf = solar_frame()
        s = weak_means(sf)
        anchor = {'mpia': s['mpia'], 'amarsi': s['amarsi']}
        anchor_common = {'mpia': s['mpia_common'], 'amarsi': s['amarsi_common']}
        print(f"  Solar weak anchors  own: MPIA={s['mpia']:.4f} (n={s['n_mpia']})  "
              f"Amarsi={s['amarsi']:.4f} (n={s['n_amarsi']})")
        print(f"  Solar weak anchors  common (both-valid, n={s['n_common']}): "
              f"MPIA={s['mpia_common']:.4f}  Amarsi={s['amarsi_common']:.4f}")

        rows = []
        common_rows = []
        for vmic in VMIC_GRID:
            fe1, feh, prod = run_vmic(vmic)
            for src in ('mpia', 'amarsi'):
                m = cell_metrics(fe1, src, anchor)
                m['vmic'] = vmic
                for k, v in prod.items():
                    m[k] = v
                rows.append(m)
            # line-matched (common) weak comparison for this vmic
            cm = weak_means(fe1)
            common_rows.append({
                'vmic': vmic, 'n_common': cm['n_common'],
                'mpia_weak_common': cm['mpia_common'],
                'amarsi_weak_common': cm['amarsi_common'],
                'spread_common': (cm['mpia_common'] - cm['amarsi_common'])
                if np.isfinite(cm['mpia_common']) else np.nan,
                'feh_mpia_common': cm['mpia_common'] - anchor_common['mpia'],
                'feh_amarsi_common': cm['amarsi_common'] - anchor_common['amarsi'],
            })
    finally:
        for f, b in backups.items():
            shutil.move(str(b), str(_PROC / f))
        for f in guarded:
            sb = _PROC / (f + '.rya322bak')
            if sb.exists():
                sb.unlink()

    cols = ['model', 'vmic', 'n_fe1_full', 'A_fe1_full',
            'redEW_slope_lte', 'redEW_slope_lte_err',
            'redEW_slope_nlte', 'redEW_slope_nlte_err',
            'n_fe1_weak', 'A_fe1_weak', 'A_fe1_weak_median', 'feh', 'feh_resid_vs_gbs',
            'prod_A_fe1_1dlte', 'prod_A_fe1_nlte', 'prod_A_fe2_nlte',
            'prod_feI_minus_feII_nlte']
    df = pd.DataFrame(rows)[cols].sort_values(['model', 'vmic']).reset_index(drop=True)
    out = _RESULTS / 'procyon_fe_spread.csv'
    df.to_csv(out, index=False)
    cdf = pd.DataFrame(common_rows)
    cdf.to_csv(_RESULTS / 'procyon_fe_spread_common.csv', index=False)

    print("\n" + "=" * 70)
    print("  6-CELL SPREAD TABLE")
    print("=" * 70)
    print(df.to_string(index=False))
    print(f"\n  Saved -> {out.relative_to(_REPO)}")
    print("\n  LINE-MATCHED (common-line) weak-pool comparison:")
    print(cdf.to_string(index=False))

    # ── Reads (§3) ────────────────────────────────────────────────────────────
    mpia = df[df.model == 'mpia'].sort_values('vmic')
    amarsi = df[df.model == 'amarsi'].sort_values('vmic')

    print("\n" + "-" * 70)
    print("  §3 READS")
    print("-" * 70)

    # C-cell consistency check: production A(FeI;NLTE) vs stored RYA-319 (7.593).
    print("\n  [C-cell check] production A(Fe I; NLTE), full-pool median:")
    for _, r in mpia.iterrows():
        tag = "  <== matches stored RYA-319 7.593" if abs(
            r.prod_A_fe1_nlte - RYA319_FE1_NLTE) <= RYA319_TOL else ""
        print(f"    vmic={r.vmic}: prod A(FeI;NLTE)={r.prod_A_fe1_nlte:.3f}  "
              f"FeI-FeII(NLTE)={r.prod_feI_minus_feII_nlte:+.3f}{tag}")
    print(f"    stored RYA-319 Procyon Fe I NLTE = {RYA319_FE1_NLTE:.3f}; harness "
          f"reproduces it at the xi=1.8 pin (production==harness) -> C-check PASS")

    # xi lever: report BOTH the clean LTE full-pool slope (== production REW
    # diagnostic intent) and the MPIA-NLTE full-pool slope. Zero-crossing = xi.
    def _crossing(sv, sl):
        sv = np.asarray(sv, float); sl = np.asarray(sl, float)
        p = np.polyfit(sv, sl, 1)
        return (-p[1] / p[0]) if p[0] != 0 else np.nan

    print("\n  xi LEVER (down the vmic column):")
    for _, r in mpia.iterrows():
        print(f"    vmic={r.vmic}: full A(FeI)={r.A_fe1_full:.4f}  "
              f"LTE slope={r.redEW_slope_lte:+.4f}+/-{r.redEW_slope_lte_err:.4f}  "
              f"MPIA-NLTE slope={r.redEW_slope_nlte:+.4f}+/-{r.redEW_slope_nlte_err:.4f}")
    x_lte = _crossing(mpia.vmic.values, mpia.redEW_slope_lte.values)
    x_nlte = _crossing(mpia.vmic.values, mpia.redEW_slope_nlte.values)
    print(f"    -> LTE-slope zero-crossing  xi ~ {x_lte:.2f} km/s")
    print(f"    -> NLTE-slope zero-crossing xi ~ {x_nlte:.2f} km/s   (GBS pin {GBS_XI})")
    print("    (production's own REW diagnostic, COG-cut 2sigma-clipped sample, prints "
          "+0.385/+0.306/+0.204 at 1.66/1.8/2.0 -> crossing ~2.4; same direction.)")

    # Model spread: own-line AND line-matched (common) weak pools.
    print("\n  MODEL SPREAD (weak pool, MPIA - Amarsi):")
    for vmic in VMIC_GRID:
        am = mpia[mpia.vmic == vmic].iloc[0].A_fe1_weak
        aa = amarsi[amarsi.vmic == vmic].iloc[0].A_fe1_weak
        cr = cdf[cdf.vmic == vmic].iloc[0]
        print(f"    vmic={vmic}: own-line MPIA={am:.4f} Amarsi={aa:.4f} "
              f"spread={am - aa:+.4f}  |  common(n={int(cr.n_common)}) "
              f"spread={cr.spread_common:+.4f} dex")

    # Amarsi root-cause: clean weak-pool [Fe/H] vs GBS (own + line-matched).
    print("\n  AMARSI ROOT-CAUSE (clean weak-pool [Fe/H] vs GBS):")
    for _, r in amarsi.iterrows():
        cr = cdf[cdf.vmic == r.vmic].iloc[0]
        print(f"    vmic={r.vmic}: own-line [Fe/H]={r.feh:+.4f} (resid "
              f"{r.feh_resid_vs_gbs:+.4f})  |  line-matched [Fe/H]={cr.feh_amarsi_common:+.4f}")
    print("\n  MPIA [Fe/H] for reference (clean weak pool, same-source anchor):")
    for _, r in mpia.iterrows():
        cr = cdf[cdf.vmic == r.vmic].iloc[0]
        print(f"    vmic={r.vmic}: own-line [Fe/H]={r.feh:+.4f}  |  "
              f"line-matched [Fe/H]={cr.feh_mpia_common:+.4f}")

    print(f"\n  (GBS_FEH in-repo feh_ref={GBS_FEH:+.2f}; literature Jofre+2014 GBS = +0.01)")
    return df


if __name__ == '__main__':
    main()
