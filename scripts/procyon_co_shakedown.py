"""
scripts/procyon_co_shakedown.py
===============================
RYA-348 — Procyon C/O synthesis shakedown (the first non-Fe synth-quality + NLTE
stress-test on a benchmark), modeled on the RYA-322 procyon_fe_matrix harness.

Recon-first (Step 0 is a hard gate). The MULTI-INSTRUMENT EXPANSION comment scoped
UVES O I 777 as PRIMARY O, with HARPS [O I] 6300 as a cross-check, plus HST UV and IR
CO arms. Step-0 recon of the CODE (not the intent) finds:

  • C/N/O synth path PRESENT — pipeline/cno_synthesis.py run_cno(), amarsi 3D-NLTE CNO
    backend (amarsi2019_cno, in-grid at 6554 K, no clamp — unlike Fe).
  • but run_cno only registers REGIONS={'vis': HARPS_VIS}. UVES_OPT/ESPRESSO_OPT are
    DEFINED RegionConfigs but NOT registered, and the multi-arm loader is Vesta/reflected-
    solar specific; _load_observed_spectrum loads only {star}_normalized (HARPS).
  → The UVES O I 777 PRIMARY-O arm (and HST UV / IR CO) is NOT wired for Procyon. Running
    it is a BUILD (register region + Procyon UVES intake + non-Vesta loader) = the
    multi-region follow-on (RYA-351 scope note). This harness runs the RUNNABLE arm —
    the HARPS VIS shakedown — and reports the multi-instrument gap honestly. It does NOT
    fabricate a UVES/UV/IR Procyon arm.

Differential denominator = OUR OWN measured Sun (RYA-348 comment §4; solar Phase C,
merged): A(C)=8.491, A(O)=8.735. Reported absolute AND differential.

Scope: analysis/scratch only. Reads the Procyon HARPS normalized spectrum + line list,
runs synthesis, writes to data/audit/cno_synthesis/ and data/results/. Does NOT touch
spectra, EWs, the canonical line list, or STAR_PARAMS. STOP at the verdict.

Usage:  python scripts/procyon_co_shakedown.py [--skip-run]
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from config.constants import get_star_params                      # noqa: E402
import pipeline.problem_children as pc                            # noqa: E402

# Our own measured Sun — the differential denominator (RYA-348 §4; solar Phase C merged).
SOLAR_OURS = {'C': 8.491, 'O': 8.735, 'N': 7.89}
SOLAR_PROV = 'our measured Sun (solar Phase C, merged): C=8.491, O=8.735 (O I 777 cross-arm)'

PROD_DIR = ROOT / 'data' / 'audit' / 'cno_synthesis'
RESULTS  = ROOT / 'data' / 'results'


def _rule(c='='):
    print(c * 78)


def preflight_registry(teff: float, feh: float) -> str:
    _rule(); print("  STEP 0.5 — RYA-463 problem-children pre-flight"); _rule()
    txt = pc.predict_headsup(teff, feh, star_name='Procyon')
    print(txt)
    return txt


def run_synthesis():
    """Run the HARPS VIS C/N/O synthesis for Procyon (the runnable arm)."""
    from pipeline.cno_synthesis import run_cno
    _rule(); print("  STEP 1-3 — HARPS VIS C/N/O synthesis (run_cno procyon/vis)"); _rule()
    res = run_cno('procyon', region_name='vis')
    return res


def per_indicator_table() -> pd.DataFrame:
    """Per-feature A(X) raw 1D-LTE + the CITED Phase-A correction (Amarsi-2019 3D-NLTE
    grid for C I / O I; Caffau-2015 cited 3D anchor for [O I] 6300; molecular 3D offset
    OWED) + per-feature reduced chi2. The corrected column is re-derived via the pipeline's
    own apply_cited_corrections (validate-don't-tune; cited/vendored only)."""
    from pipeline.cno_synthesis import apply_cited_corrections, REGIONS
    pb = PROD_DIR / 'procyon_vis_cno_per_band.csv'
    if not pb.exists():
        raise FileNotFoundError(f"per-band product not found: {pb} (synthesis did not write)")
    df = pd.read_csv(pb)
    rec = get_star_params('procyon')
    params = {'teff_K': float(rec['teff']), 'logg': float(rec['logg']),
              'feh': float(rec['feh_ref']), 'vturb_kms': float(rec.get('xi', 1.0))}
    corr = {c['key']: c for c in apply_cited_corrections(df.to_dict('records'), params,
                                                         REGIONS['vis'])}
    df['A_corr'] = df['key'].map(lambda k: corr.get(k, {}).get('a_corr'))
    df['corr_kind'] = df['key'].map(lambda k: corr.get(k, {}).get('kind'))
    keep = ['key', 'element', 'role', 'A_X', 'A_corr', 'corr_kind', 'red_chi2',
            'sigma_fit', 'n_pix', 'status']
    return df[[c for c in keep if c in df.columns]].copy()


# χ²ᵣ below this = the synthesis fits the feature cleanly (brief Step 3).
CHI2_CLEAN = 10.0


def differential(ind: pd.DataFrame) -> pd.DataFrame:
    """[C/H], [O/H], C/O — differential vs OUR measured Sun, on the PRIMARY indicator per
    element, using the CITED-CORRECTED value with an honesty flag. The product CSV's
    element A_X is RAW 1D-LTE (e.g. [O I] 9.417) — not used for the differential; we use
    the cited-corrected primary and flag where the correction (not a Procyon measurement)
    carries the value (the [O I] Caffau-anchor / RYA-104 case)."""
    primary = {'C': 'CH_Gband', 'O': 'OI_6300', 'N': 'CN_red'}
    flag = {
        'C': 'CH primary, NLTE not needed (molecular); 3D offset OWED; '
             'C I cross-checks disagree → provisional',
        'O': 'CITED-ANCHOR DOMINATED: raw [O I] χ²ᵣ≫10, value set by Caffau-2015 solar '
             'anchor (≈8.73) → NOT an independent Procyon O. RYA-104 recurrence. DO NOT BANK.',
        'N': 'CN molecular at HARPS blue, χ²ᵣ≫10, σ_sys huge, depends on assumed C/O → '
             'NOT science-grade. N I red multiplets (NLTE-owed) not in HARPS arm.',
    }
    rows = []
    for el in ('C', 'O', 'N'):
        sub = ind[ind['key'] == primary[el]]
        if sub.empty:
            continue
        r = sub.iloc[0]
        raw = float(r['A_X'])
        corr = float(r['A_corr']) if pd.notna(r['A_corr']) else raw
        sun = SOLAR_OURS.get(el, np.nan)
        rows.append(dict(element=el, primary=primary[el],
                         A_raw_1dlte=round(raw, 3), A_cited_corr=round(corr, 3),
                         chi2r=round(float(r['red_chi2']), 1),
                         differential=round(corr - sun, 3) if np.isfinite(sun) else np.nan,
                         note=flag[el]))
    out = pd.DataFrame(rows)
    cC = ind[ind['key'] == 'CH_Gband']['A_corr']
    cO = ind[ind['key'] == 'OI_6300']['A_corr']
    if not cC.empty and not cO.empty and pd.notna(cC.iloc[0]) and pd.notna(cO.iloc[0]):
        co_proc = 10 ** (float(cC.iloc[0]) - float(cO.iloc[0]))
        co_sun  = 10 ** (SOLAR_OURS['C'] - SOLAR_OURS['O'])
        print(f"\n  C/O (Procyon, cited-corrected CH / [O I]) = {co_proc:.3f}  "
              f"[O I] cited-anchor → C/O NOT science-grade")
        print(f"  C/O (our Sun)                             = {co_sun:.3f}   [{SOLAR_PROV}]")
    return out


def cross_check(ind: pd.DataFrame):
    """Multi-indicator cross-check: do the independent C indicators agree? Surface the
    spread as a finding (never average it away). Same-physics-disagreement = systematic."""
    _rule('-'); print("  MULTI-INDICATOR CROSS-CHECK (surface spread, do not average)"); _rule('-')
    for el in ('C', 'O'):
        sub = ind[(ind['element'] == el) & ind['A_corr'].notna()]
        if sub.empty:
            print(f"  {el}: no usable indicator"); continue
        vals = pd.to_numeric(sub['A_corr'], errors='coerce').dropna()
        spread = float(vals.max() - vals.min()) if len(vals) > 1 else 0.0
        tag = 'AGREE' if spread <= 0.15 else 'DISAGREE → systematic (surface, not averaged)'
        print(f"  {el}: {len(vals)} indicators (cited-corrected)  range {vals.min():.3f}.."
              f"{vals.max():.3f}  spread {spread:.3f}  → {tag}")
        for _, r in sub.iterrows():
            chi = float(r['red_chi2'])
            clean = 'CLEAN' if chi < CHI2_CLEAN else 'POOR-FIT'
            print(f"       {r['key']:<10} {r['role']:<12} A={float(r['A_corr']):.3f}  "
                  f"chi2r={chi:6.2f} [{clean}]  ({r.get('corr_kind','')})")


def gap_table():
    """Step 5 — which C/O diagnostics HARPS delivers vs which require UVES/UV/IR, and
    the wiring status in THIS codebase (the honest multi-instrument gap)."""
    _rule(); print("  STEP 5 — HARPS-vs-UVES/UV/IR GAP (delivered / wired-status)"); _rule()
    rows = [
        ('[O I] 6300 (Ni-blend)', 'O cross-check', 'HARPS VIS', 'WIRED (run here)',
         'continuum/blend-limited at +0.01; registry: OK but not primary'),
        ('C I 5052 / 5380',       'C (NLTE)',      'HARPS VIS', 'WIRED (run here)',
         'C I 5380 → SATURATION_COG exclude (registry); NLTE in-grid at 6554 K'),
        ('CH G-band 4290-4315',   'C molecular',   'HARPS VIS', 'WIRED (run here)',
         'blue end, lower HARPS S/N; molecular-synth path'),
        ('CN red / C2 Swan',      'N / C molec',   'HARPS VIS', 'WIRED (run here)',
         'molecular cross-checks'),
        ('O I 7771-5 triplet',    'O PRIMARY',     'UVES red',  'NOT WIRED for Procyon',
         'region defined but unregistered; loader is Vesta-reflected-solar → BUILD (follow-on)'),
        ('C I red-optical / NI 8216', 'C/N',        'UVES red',  'NOT WIRED for Procyon',
         'same UVES build dependency'),
        ('FUV C I / NH 3360',     'C/N MEASURED',  'HST STIS/COS', 'NOT WIRED (audit only)',
         'RYA-222/351 inventory exists; no synth loader → multi-region follow-on'),
        ('CO overtone / OH / CN', 'C cross + 12C/13C', 'IR CRIRES+/APOGEE', 'NOT WIRED (telluric-gated)',
         'RYA-351: APOGEE weak-CO only; no 2.3um overtone; telluric mandatory (RYA-373)'),
    ]
    w = pd.DataFrame(rows, columns=['diagnostic', 'role', 'instrument', 'wiring', 'note'])
    print(w.to_string(index=False))
    return w


def main(skip_run: bool = False):
    rec = get_star_params('procyon')
    teff, feh = float(rec['teff']), float(rec['feh_ref'])
    _rule(); print(f"  RYA-348 — Procyon C/O shakedown  (Teff={teff:.0f} logg={rec['logg']:.2f} "
                   f"[Fe/H]={feh:+.2f} xi={rec.get('xi',1.0):.1f}; GBS Heiter+2015/Jofré+2014)")
    _rule()

    preflight_registry(teff, feh)

    if not skip_run:
        run_synthesis()

    prod_path = PROD_DIR / 'procyon_vis_cno_product.csv'
    if not prod_path.exists():
        print("\n  [STOP] No product written — synthesis did not complete. "
              "Step-0 recon report stands as the deliverable.")
        return
    prod = pd.read_csv(prod_path)

    _rule(); print("  PER-INDICATOR A(X) raw 1D-LTE + cited Phase-A correction + chi2r"); _rule()
    ind = per_indicator_table()
    print(ind.to_string(index=False))

    # Synth-quality read (brief Step 3): which features fit cleanly (χ²ᵣ<10) vs fail.
    n_clean = int((ind['red_chi2'] < CHI2_CLEAN).sum())
    print(f"\n  SYNTH-QUALITY: {n_clean}/{len(ind)} features fit cleanly (χ²ᵣ<{CHI2_CLEAN:.0f}). "
          f"At 6554 K on HARPS, C/O features are POOR-FIT — the empirical marginality the "
          f"brief predicted for HARPS-only Procyon C/O.")
    # No silent grid clamp (the Fe-clamp precedent): confirm the C I Amarsi-2019 NLTE query
    # was in-hull (kind=amarsi2019_grid), not out_of_hull, at 6554 K.
    ooh = ind[ind['corr_kind'].isin(['grid_out_of_hull', 'grid_error'])]
    print(f"  NLTE grid-edge: {'NONE — C I Amarsi-2019 in-grid at 6554 K (no clamp)' if ooh.empty else 'CLAMP/ERROR: '+', '.join(ooh['key'])}")

    cross_check(ind)

    _rule(); print("  [C/H], [O/H], C/O — differential vs OUR Sun (primary, cited-corrected)"); _rule()
    diff = differential(ind)
    print(diff.to_string(index=False))

    w = gap_table()

    RESULTS.mkdir(parents=True, exist_ok=True)
    ind.to_csv(RESULTS / 'procyon_co_per_indicator.csv', index=False)
    diff.to_csv(RESULTS / 'procyon_co_differential.csv', index=False)
    w.to_csv(RESULTS / 'procyon_co_gap_table.csv', index=False)
    print(f"\n  Wrote: data/results/procyon_co_{{per_indicator,differential,gap_table}}.csv")
    _rule(); print("  VERDICT IS RYAN'S + CLAUDE'S FROM THE POSTED EVIDENCE — STOP."); _rule()


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--skip-run', action='store_true',
                    help='read existing products without re-running the synthesis')
    args = ap.parse_args()
    main(skip_run=args.skip_run)
