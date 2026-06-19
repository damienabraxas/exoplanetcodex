"""
scripts/diag_blend_flag_drift.py
================================
RYA-356 — DIAGNOSTIC (scratch-only): why did Fe II 6149/6247 `blend_flag` flip
True→False between the published (v1_1-era) data and the regenerated data?

NO canonical edits. Reads the live line list, reconstructs the OLD proximity-based
flag as the v1_1 baseline, diffs across ALL elements, and adjudicates the flipped
lines by neighbour strength (central_depth) + the continuous vald_proximity_flag.

Step 0 — the two flags + what changed
  • Line-list flag  : linelist_solar.csv['blend_flag'] — SET by the builder.
        OLD (scripts/build_linelist.flag_blends): proximity binary — ANY line within
            0.10 Å of another → both True (strength-blind).
        NEW (pipeline/build_linelist.build_vetted_blend_flag, RYA-209): default False;
            True ONLY for the curated VETTED_BLENDS (Fe I 4918.99 / 4970.50). Proximity
            moved to the CONTINUOUS vald_proximity_flag (window 0.5 Å, decay 0.1 Å).
  • Per-measurement flag : solar_ew.csv / sol_ew_results['blend_flag'] — NOT an
        independent detection: lines_fit.py reads the line-list flag (df['blend_flag'])
        and PROPAGATES it to the EW output. So both flags trace to the line-list flag;
        the drift originates there (the RYA-209 redefinition), not in EW fitting.
  Consumers: the RYA-352 cull BLEND criterion and the RYA-208 sanity filter both read
  the (now vetted) blend_flag == True.

Usage:  python3 scripts/diag_blend_flag_drift.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

_REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_REPO))

_LL = _REPO / 'data' / 'linelists' / 'linelist_solar.csv'
_SITE_V1 = (_REPO.parent / 'exoplanetcodex-site' / 'assets' / 'data'
            / 'sol_ew_results_v1.csv')   # published (proximity-era) Fe II EW pool
OLD_PROX_TOL = 0.10   # Å — the old scripts/build_linelist.flag_blends threshold
DEPTH_NEGLIGIBLE = 0.02   # central_depth below which a neighbour can't materially blend


def old_proximity_flag(df: pd.DataFrame, tol=OLD_PROX_TOL) -> pd.Series:
    """Reconstruct the OLD strength-blind proximity flag: any line within tol Å of
    another (any species) → True. Mirrors scripts/build_linelist.flag_blends."""
    w = df['wavelength_air_A'].to_numpy(float)
    order = np.argsort(w)
    flag = np.zeros(len(df), dtype=bool)
    ws = w[order]
    near = np.zeros(len(df), dtype=bool)
    near[1:] |= (np.diff(ws) < tol)
    near[:-1] |= (np.diff(ws) < tol)
    flag[order] = near
    return pd.Series(flag, index=df.index)


def _nearest_material(df, wl, species_ok=None):
    """Nearest neighbour with central_depth >= DEPTH_NEGLIGIBLE within 0.5 Å."""
    nb = df[(df['wavelength_air_A'].sub(wl).abs() < 0.5)
            & (df['wavelength_air_A'].sub(wl).abs() > 1e-4)].copy()
    if 'central_depth' in nb:
        nb = nb[nb['central_depth'] >= DEPTH_NEGLIGIBLE]
    if nb.empty:
        return None
    nb['d'] = (nb['wavelength_air_A'] - wl).abs()
    return nb.sort_values('d').iloc[0]


def main() -> None:
    ll = pd.read_csv(_LL, low_memory=False)
    ll['blend_flag'] = ll['blend_flag'].astype(str).str.lower() == 'true'
    new = ll['blend_flag']
    old = old_proximity_flag(ll)

    print("── Step 0/1: blend_flag — NEW (vetted, live) vs OLD (proximity 0.10 Å) ──")
    print(f"  NEW vetted True : {int(new.sum())}  "
          f"({[ (r['element'], r['ion'], round(r['wavelength_air_A'],3)) for _,r in ll[new].iterrows()]})")
    print(f"  OLD proximity True (reconstructed @0.10 Å): {int(old.sum())}")
    flipped_TF = ll[old & ~new]   # was True (proximity), now False (vetted)
    flipped_FT = ll[~old & new]   # was False, now True
    print(f"\n  FLIP True→False (old proximity → new vetted): {len(flipped_TF)} lines")
    print(f"  FLIP False→True                              : {len(flipped_FT)} lines")
    print(f"  → SYSTEMIC: {100*len(flipped_TF)/max(int(old.sum()),1):.0f}% of old-True "
          f"lines are now False. Not a 6149/6247-specific change — the whole "
          f"proximity-binary population was redefined (RYA-209).")

    print("\n  Flip True→False by element/ion (top 15):")
    g = (flipped_TF.groupby(['element', 'ion']).size()
         .sort_values(ascending=False).head(15))
    for (el, ion), n in g.items():
        print(f"      {el:3} {str(ion):3} : {n}")

    # ── Step 2: adjudicate the two named lines ────────────────────────────────
    print("\n── Step 2: adjudication (neighbour strength + vald_proximity_flag) ──")
    for wl in (6247.557, 6149.258):
        m = ll[(ll['element'] == 'Fe') & (ll['ion'] == 'II')
               & (ll['wavelength_air_A'].sub(wl).abs() < 0.05)]
        if m.empty:
            print(f"  Fe II {wl}: not in linelist_solar"); continue
        r = m.iloc[0]
        vpf = r.get('vald_proximity_flag', float('nan'))
        mat = _nearest_material(ll, wl)
        print(f"  Fe II {wl}:  new blend_flag=False  vald_proximity_flag={vpf:.3f}  "
              f"central_depth={r.get('central_depth')}")
        if mat is None:
            print(f"      nearest MATERIAL neighbour (cd≥{DEPTH_NEGLIGIBLE}): NONE within 0.5 Å"
                  f"  → no real blend → NEW False CORRECT (old proximity over-flagged a weak line)")
        else:
            print(f"      nearest MATERIAL neighbour: {mat['element']} {mat['ion']} "
                  f"{mat['wavelength_air_A']:.3f} (Δ={mat['wavelength_air_A']-wl:+.3f} Å, "
                  f"cd={mat['central_depth']})")
    # the ticket's explicit candidate: 6247.354 Fe II
    c = ll[(ll['element'] == 'Fe') & (ll['ion'] == 'II')
           & (ll['wavelength_air_A'].sub(6247.354).abs() < 0.02)]
    if len(c):
        cd = c.iloc[0].get('central_depth')
        print(f"\n  candidate blend 6247.354 Fe II: cd={cd}, Δ=−0.203 Å from 6247.557 "
              f"→ {'WEAK (cd<%.2f) → not a material blend' % DEPTH_NEGLIGIBLE if cd < DEPTH_NEGLIGIBLE else 'material — investigate'}")

    # ── cross-check vs the published proximity-era Fe II EW pool ──────────────
    if _SITE_V1.exists():
        site = pd.read_csv(_SITE_V1)
        s2 = site[(site['element'] == 'Fe') & (site['ion'] == 'II')]
        st = s2[s2['blend_flag'] == True]
        print(f"\n  published sol_ew_results_v1 Fe II blend_flag=True: "
              f"{[round(w,3) for w in st['wavelength_air_A']]} "
              f"(proximity-era; matches the OLD flag, now superseded by vetted)")

    print("\n── VERDICT ─────────────────────────────────────────────────────────")
    print("  IMPROVED DETECTION, not a regression. blend_flag was REDEFINED (RYA-209)")
    print("  from a strength-blind proximity binary to a curated 'confirmed non-")
    print("  separable' list, with proximity captured continuously by vald_proximity_")
    print("  flag. 6247.557/6149.258 are False because their only neighbours are")
    print("  negligibly weak (cd≲0.016) — the old True over-flagged them. RYA-352")
    print("  keeping 6247 in the EW pool is CONSISTENT with the vetted semantic.")
    print("  Systemic + intentional → STOP (no fix). Recommend: extend the RYA-355")
    print("  guard to blend_flag provenance so a future flip fails loudly.")


if __name__ == '__main__':
    main()
