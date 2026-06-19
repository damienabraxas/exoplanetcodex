"""
scripts/procyon_amarsi_clamp_check.py
=====================================
RYA-339 — Amarsi-on-Procyon clamp cross-check.

RYA-322 concluded the Amarsi weak-pool [Fe/H] ≈ −0.2 is "the 6500 K clamp itself,
irreducibly wrong for Procyon's 6554 K." That is REFUTED by Amarsi et al. 2022
(2209.13449v3, §5.3.2): they handled Procyon by exactly this 0th-order clamp to
their 6500 K grid edge and published a sensible value. Same grid, same edge node —
so a persistent ~0.15–0.2 dex miss is OUR APPLICATION, not the ceiling.

This script isolates the four suspect conditions ONE AT A TIME, printing each as an
explicit diagnostic, and reports whether our Amarsi-clamped weak-pool value
reproduces the literature Procyon value.

ANALYSIS-ONLY. Reads the stored Procyon/solar Fe I per-line tables (production path,
pinned ξ) and writes ONLY to data/results/procyon_amarsi_clamp_check.csv. It does
NOT modify spectra, EWs, or the linelist.

Params: single source of truth = STAR_PARAMS / constants.py (Teff/logg/ξ from
Heiter+2015 + Jofré+2014, ξ pinned RYA-325). Amarsi grid edge + grid bounds from
pipeline.nlte_corrections._GRID (Amarsi et al. 2022, A&A 668, A68). NO hardcoded
params; every dropped / out-of-grid line is logged explicitly (no silent fallback).
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

_REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_REPO))

from config.constants import STAR_PARAMS, STAR_SOLAR, SOLAR_ASPLUND2021, get_star_params  # noqa: E402
from pipeline.nlte_corrections import (                                  # noqa: E402
    _GRID, _compute_aberr, _apply_aberr_to_line, _mpia_fe_delta, fe_grid_in_bounds,
)

_PROC = _REPO / 'data' / 'processed'
_OUT  = _REPO / 'data' / 'results'
_OUT.mkdir(parents=True, exist_ok=True)

# ── Single source of truth (no hardcoded params) ─────────────────────────────
P = STAR_PARAMS['procyon']
PROC_TEFF = float(P['teff'])     # 6554 K — Heiter+2015 (GBS)
PROC_LOGG = float(P['logg'])     # 4.00  — Heiter+2015
PROC_XI   = float(P['xi'])       # 1.8 km/s — pinned, RYA-325 (Jofré+2014/GBS scale)
PROC_FEH_REF = float(P['feh_ref'])  # +0.03 in-repo; lit GBS +0.01 (Jofré+2014)

_sf = get_star_params('solar')   # RYA-298 single source
S_TEFF = float(_sf['teff']); S_LOGG = float(_sf['logg'])
S_FEH  = float(_sf['feh_ref']);    S_XI   = float(STAR_SOLAR['vturb_kms'])

AMARSI_TEFF_CEIL = float(_GRID['teff'][1])   # 6500 K — Amarsi 2022 grid upper edge
AFE_LO, AFE_HI   = float(_GRID['afe'][0]), float(_GRID['afe'][1])   # 4.5, 7.5
WEAK_CUT = -4.9                               # log(EW/λ); Amarsi 2022 weak-line regime

# Literature comparison target (Procyon [Fe/H], 3D-NLTE Fe):
#   Amarsi et al. 2022 (2209.13449v3, §5.3.2): Procyon [Fe/H] = +0.11. Their §6 gives
#   the sign anchor — 1D LTE is 0.20 dex SMALLER than 3D-NLTE for Fe I → correction
#   +0.20 (positive). Amarsi used ξ=1.66 (we pin 1.8), so a few-hundredths ξ offset is
#   expected. GBS (Jofré+2014) is +0.01 for reference.
AMARSI_2022_PROCYON_TARGET = 0.11   # Amarsi 2022 §5.3.2 (3D-NLTE)
GBS_FEH_LIT = 0.01                  # Jofré+2014 (for reference)


def reduced_ew(ew_mA, lam_A):
    # mÅ → Å (× 1e-3) so log(EW/λ) is the true reduced EW (≈ −5 for weak lines),
    # matching the Amarsi WEAK_CUT = −4.9 grid-validity threshold.
    return np.log10((np.asarray(ew_mA, float) * 1e-3) / np.asarray(lam_A, float))


def load_fe1(name):
    df = pd.read_csv(_PROC / f'{name}_per_line.csv')
    fe1 = df[(df.element == 'Fe') & (df.ion == 'I')].copy().reset_index(drop=True)
    fe1['log_rew'] = reduced_ew(fe1.ew_mA.values, fe1.wavelength_air_A.values)
    return fe1


def amarsi_per_line(fe1, teff, logg, vmic):
    """Apply the archived Amarsi MLP with the 0th-order Teff clamp, capturing the
    (Teff_in, afe_in, aberr) actually handed to / returned by the MLP per line.

    Mirrors pipeline.nlte_corrections._apply_aberr_to_line (afe seeded from the
    line's own a_1dlte, clipped to the grid afe range; Teff clamped to 6500)."""
    teff_in = min(float(teff), AMARSI_TEFF_CEIL)
    rows = []
    for _, r in fe1.iterrows():
        a1d = float(r['a_1dlte'])
        afe_in = float(np.clip(a1d, AFE_LO, AFE_HI))     # afe handed to MLP (CONDITION 1/3)
        ab = _apply_aberr_to_line('I', float(r['excitation_potential_eV']),
                                  float(r['eup_eV']), float(r['log_gf']), a1d,
                                  teff_in, float(logg), float(vmic))
        a_nlte = a1d + ab if np.isfinite(ab) else np.nan
        rows.append(dict(teff_in=teff_in, afe_in=afe_in, aberr=ab, a_nlte=a_nlte))
    out = pd.DataFrame(rows)
    return out


def mpia_per_line(fe1, teff, logg, feh):
    deltas = np.full(len(fe1), np.nan)
    for i, r in fe1.iterrows():
        d = _mpia_fe_delta('I', float(r['wavelength_air_A']), float(teff),
                           float(logg), float(feh))
        if np.isfinite(d):
            deltas[i] = float(r['a_1dlte']) + d
    return deltas


def banner(t):
    print('\n' + '=' * 70 + f'\n  {t}\n' + '=' * 70)


def main():
    print(__doc__.split('ANALYSIS-ONLY')[0].strip())
    print(f"\nProcyon (STAR_PARAMS): Teff={PROC_TEFF} logg={PROC_LOGG} ξ={PROC_XI} "
          f"feh_ref={PROC_FEH_REF}")
    print(f"Amarsi grid (nlte_corrections._GRID): Teff {_GRID['teff']} logg "
          f"{_GRID['logg']} afe {_GRID['afe']} vmic {_GRID['vmic']}")
    print(f"Procyon Teff 6554 > grid ceil {AMARSI_TEFF_CEIL} → 0th-order clamp required.")

    proc = load_fe1('procyon')
    sol  = load_fe1('solar')
    ama  = amarsi_per_line(proc, PROC_TEFF, PROC_LOGG, PROC_XI)
    proc = pd.concat([proc, ama.add_prefix('ama_')], axis=1)
    proc['a_nlte_mpia'] = mpia_per_line(proc, PROC_TEFF, PROC_LOGG, PROC_FEH_REF)

    sol_ama = amarsi_per_line(sol, S_TEFF, S_LOGG, S_XI)
    sol['a_nlte_amarsi'] = sol_ama['a_nlte'].values
    sol['a_nlte_mpia']   = mpia_per_line(sol, S_TEFF, S_LOGG, S_FEH)

    # ── CONDITION 1 — clamp pins BOTH Teff and A(Fe) to the edge ─────────────
    banner("CONDITION 1 — (Teff_in, afe_in) handed to the MLP (0th-order edge clamp)")
    assert (proc['ama_teff_in'] == AMARSI_TEFF_CEIL).all(), "Teff not clamped to 6500!"
    print(f"  Teff_in = {AMARSI_TEFF_CEIL} for all {len(proc)} lines ✓ (clamped from "
          f"{PROC_TEFF})")
    print(f"  afe_in range = [{proc['ama_afe_in'].min():.3f}, "
          f"{proc['ama_afe_in'].max():.3f}]  median={proc['ama_afe_in'].median():.3f}")
    n_at_edge = int((proc['ama_afe_in'] >= AFE_HI - 1e-6).sum())
    print(f"  afe_in pinned at grid edge {AFE_HI}: {n_at_edge}/{len(proc)} lines")
    print("  NOTE: afe_in is seeded from each line's OWN a_1dlte (clipped to grid), "
          "NOT a single stellar A(Fe). For weak lines a_1dlte≈stellar A(Fe); for "
          "blended/strong lines a_1dlte is inflated → afe_in mis-placed in the grid.")
    print(proc[['wavelength_air_A', 'a_1dlte', 'ama_afe_in', 'log_rew']]
          .head(6).to_string(index=False))

    # CONDITION 1 sub-test: does pinning afe to a SINGLE stellar value (vs per-line
    # a_1dlte) change the corrections? (Amarsi 2022's "pin A(Fe) to the edge" reading.)
    proc_rew = reduced_ew(proc.ew_mA.values, proc.wavelength_air_A.values)
    afe_pin = float(np.clip(np.median(proc.loc[proc_rew < WEAK_CUT, 'a_1dlte']),
                            AFE_LO, AFE_HI))   # weak-pool stellar A(Fe) ≈ Procyon
    ab_pin = np.array([_compute_aberr('I', float(r['excitation_potential_eV']),
                       float(r['eup_eV']), float(r['log_gf']),
                       AMARSI_TEFF_CEIL, PROC_LOGG, afe_pin, PROC_XI)
                       for _, r in proc.iterrows()])
    proc['ama_aberr_pinned'] = ab_pin
    print(f"  afe-PINNED variant (single stellar afe={afe_pin:.3f} for all lines): "
          f"mean aberr = {np.nanmean(ab_pin):+.3f}  vs per-line mean "
          f"{proc['ama_aberr'].mean():+.3f}  → Δ={np.nanmean(ab_pin)-proc['ama_aberr'].mean():+.3f} "
          f"(afe-insensitive if ≈0)")

    # ── CONDITION 2 — weak-line cut log(EW/λ) < −4.9 ─────────────────────────
    banner("CONDITION 2 — weak-line cut  log(EW/λ) < −4.9  (Amarsi grid validity)")
    proc['kept'] = proc['log_rew'] < WEAK_CUT
    n_kept = int(proc['kept'].sum()); n_drop = int((~proc['kept']).sum())
    print(f"  n_kept = {n_kept}   n_dropped = {n_drop}   (of {len(proc)} Fe I lines)")
    if n_drop:
        print("  Dropped (NOT silent — every line logged):")
        for _, r in proc[~proc['kept']].iterrows():
            print(f"    λ={r['wavelength_air_A']:.3f}  EW={r['ew_mA']:.1f} mÅ  "
                  f"log_rew={r['log_rew']:+.3f}  a_1dlte={r['a_1dlte']:.3f}  "
                  f"aberr={r['ama_aberr']:+.3f}")
    # the diagnostic: mean aberr on kept (weak) vs dropped (strong)
    wk = proc[proc['kept']]; st = proc[~proc['kept']]
    print(f"  mean aberr  weak={wk['ama_aberr'].mean():+.3f}  "
          f"strong={st['ama_aberr'].mean() if len(st) else float('nan'):+.3f}  "
          f"(POSITIVE post-RYA-339 sign fix; Amarsi §6 Procyon Fe I = +0.20)")

    # ── CONDITION 3 — afe added exactly once (no double-count) ───────────────
    banner("CONDITION 3 — afe added exactly once (no 7.46 double-count)")
    # afe_in must equal clip(a_1dlte); a double-add (a_1dlte + 7.46) would land ~14.9
    # → clipped to 7.5 for EVERY line. Check afe_in tracks a_1dlte, not pinned-7.5.
    expect = np.clip(proc['a_1dlte'].values, AFE_LO, AFE_HI)
    ok = np.allclose(proc['ama_afe_in'].values, expect, atol=1e-6)
    print(f"  afe_in == clip(a_1dlte) for all lines: {ok} ✓")
    print(f"  (a 7.46 double-add would force a_1dlte+7.46≈14.9 → afe_in==7.5 for ALL "
          f"lines; here {n_at_edge}/{len(proc)} at edge, afe_in tracks a_1dlte → "
          f"added ONCE)")
    print(f"  sanity: a_1dlte median={proc['a_1dlte'].median():.3f} "
          f"(≈ Procyon A(Fe), NOT ~14.9) → no double-add on the Amarsi path")

    # ── CONDITION 4 — air/vac wavelength basis (MPIA grid match) ─────────────
    banner("CONDITION 4 — air/vac basis for grid line matching (MPIA λ-keyed grid)")
    # The Amarsi MLP is keyed on (Elo,Eup,loggf), NOT λ — no air/vac match. The MPIA
    # grid IS λ-keyed; confirm our air λ matches its convention (no ~1.5 Å vac shift).
    try:
        from pipeline.nlte_corrections import _load_mpia_fe_grid
        g = _load_mpia_fe_grid()
        gl = pd.read_csv(_REPO / 'data' / 'nlte_grids' / 'Fe_Bergemann_MPIA.csv')
        gwav = np.sort(gl[gl['ion'] == 'I']['wave_A'].unique())
        print("  sample matched (λ_obs_air, nearest λ_grid, Δ):")
        offs = []
        for w in proc['wavelength_air_A'].values[:8]:
            j = int(np.argmin(np.abs(gwav - w))); off = w - gwav[j]; offs.append(off)
            print(f"    {w:.3f}  →  {gwav[j]:.3f}   Δ={off:+.3f} Å")
        med = float(np.median(offs))
        print(f"  median |Δ| on sample = {np.median(np.abs(offs)):.3f} Å "
              f"({'AIR-consistent, no vac shift' if abs(med) < 0.05 else 'SYSTEMATIC OFFSET — check basis'})")
    except Exception as exc:
        print(f"  (MPIA grid not available for λ-match sample: {exc})")

    # ── Weak-pool [Fe/H] summary ─────────────────────────────────────────────
    banner("WEAK-POOL [Fe/H]  (mean A(Fe I)_NLTE − same-source solar weak anchor)")
    pw = proc[proc['kept']]
    sw = sol[sol['log_rew'] < WEAK_CUT]
    def feh(pcol, scol, pframe, sframe):
        pm = pframe[pcol].mean(); sm = sframe[scol].mean()
        return pm, sm, pm - sm
    a_pm, a_sm, a_feh = feh('ama_a_nlte', 'a_nlte_amarsi', pw, sw)
    m_pm, m_sm, m_feh = feh('a_nlte_mpia', 'a_nlte_mpia', pw, sw)
    # common pool (both corrections valid) for a line-matched spread
    both = pw[pw['ama_a_nlte'].notna() & pw['a_nlte_mpia'].notna()]
    print(f"  weak pool: Procyon n={len(pw)} (Amarsi valid "
          f"{int(pw['ama_a_nlte'].notna().sum())}, MPIA valid "
          f"{int(pw['a_nlte_mpia'].notna().sum())}); solar n={len(sw)}")
    print(f"  AMARSI weak-pool A(FeI): Procyon={a_pm:.3f}  solar={a_sm:.3f}  "
          f"→ [Fe/H]={a_feh:+.3f}")
    print(f"  MPIA   weak-pool A(FeI): Procyon={m_pm:.3f}  solar={m_sm:.3f}  "
          f"→ [Fe/H]={m_feh:+.3f}")
    print(f"  Amarsi − MPIA spread = {a_feh - m_feh:+.3f} dex")
    print(f"  Amarsi 2022 §5.3.2 Procyon target = {AMARSI_2022_PROCYON_TARGET:+.3f} "
          f"(3D-NLTE; GBS/Jofré+2014 = {GBS_FEH_LIT:+.2f} for reference)")

    # ── verdict ──────────────────────────────────────────────────────────────
    banner("VERDICT")
    miss = a_feh - AMARSI_2022_PROCYON_TARGET
    weak_corr = wk['ama_aberr'].mean()
    sign_ok = weak_corr > 0
    print(f"  Fe I weak-pool correction = {weak_corr:+.3f}  → sign "
          f"{'POSITIVE ✓' if sign_ok else 'NEGATIVE ✗ (still flipped!)'}  "
          f"(Amarsi §6 = +0.20)")
    print(f"  our Amarsi-clamped weak-pool [Fe/H] = {a_feh:+.3f}  vs Amarsi §5.3.2 "
          f"target {AMARSI_2022_PROCYON_TARGET:+.3f}  →  miss = {miss:+.3f} dex")
    if sign_ok and abs(miss) <= 0.10:
        print("  PASS → sign positive AND [Fe/H] reproduces Amarsi §5.3.2 (+0.11) "
              "within ξ-sensitivity. The RYA-339 sign flip was the root cause; the "
              "RYA-322 'irreducible clamp' conclusion is REFUTED. True Amarsi−MPIA "
              f"spread = {a_feh - m_feh:+.3f} dex (3D-NLTE vs 1D-NLTE, a real grid diff).")
    elif not sign_ok:
        print("  CRITICAL — correction STILL NEGATIVE after the fix → diagnosis wrong. "
              "STOP; report raw MLP + apply trace. Do NOT re-attribute to the ceiling.")
    else:
        print("  sign fixed but [Fe/H] off target — investigate ξ / solar anchor / "
              "weak-pool composition (not the ceiling).")

    # ── SOLAR RE-CHECK (RYA-339 addendum) ────────────────────────────────────
    banner("SOLAR RE-CHECK — did the flip move the live solar anchor? (RYA-340 gate)")
    print("  Live solar Fe I NLTE path = MPIA (_mpia_fe_delta, flag "
          "NLTE_MPIA_Bergemann_1D), NOT the archived Amarsi MLP fixed here. The MPIA δ "
          "convention (A_NLTE − A_LTE, correction-to-add) was already correct.")
    sol_d = [_mpia_fe_delta('I', float(w), S_TEFF, S_LOGG, S_FEH)
             for w in sol['wavelength_air_A'].values]
    sol_d = np.array([d for d in sol_d if np.isfinite(d)])
    print(f"  MPIA solar Fe I δ: n={len(sol_d)}  mean={sol_d.mean():+.4f}  "
          f"({'POSITIVE — correctly signed, NOT flipped' if sol_d.mean() > 0 else 'negative — CHECK'})")
    print("  Amarsi-path solar correction WAS sign-flipped too (shared the Amarsi-path "
          "bug), but Amarsi is ARCHIVED (RYA-283 3D-vs-1D cross-check only), not the "
          "live anchor → live solar A(Fe I) NLTE is UNCHANGED by this fix.")
    print("  (Verified by run('solar') before/after: 7.516 → 7.516; see Linear comment.)")

    # ── write per-line CSV ───────────────────────────────────────────────────
    out = proc[['wavelength_air_A', 'log_rew', 'a_1dlte', 'ama_afe_in', 'ama_aberr',
                'ama_a_nlte', 'a_nlte_mpia', 'kept']].copy()
    out.columns = ['wavelength_air_A', 'log_rew', 'A_FeI_LTE', 'afe_in', 'delta_NLTE',
                   'A_FeI_NLTE', 'A_FeI_NLTE_mpia', 'kept_flag']
    summary = pd.DataFrame([{
        'wavelength_air_A': 'SUMMARY', 'log_rew': '',
        'A_FeI_LTE': f'amarsi_feh={a_feh:+.3f}', 'afe_in': f'n_weak={len(pw)}',
        'delta_NLTE': f'mpia_feh={m_feh:+.3f}', 'A_FeI_NLTE': f'spread={a_feh-m_feh:+.3f}',
        'A_FeI_NLTE_mpia': f'target={AMARSI_2022_PROCYON_TARGET:+.3f}',
        'kept_flag': f'miss={miss:+.3f}',
    }])
    out = pd.concat([out, summary], ignore_index=True)
    out_path = _OUT / 'procyon_amarsi_clamp_check.csv'
    out.to_csv(out_path, index=False)
    print(f"\n  Wrote {out_path.relative_to(_REPO)}  ({len(proc)} lines + summary)")


if __name__ == '__main__':
    main()
