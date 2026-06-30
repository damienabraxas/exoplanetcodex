"""
scripts/verify_solar_o_modeling_rya485.py
=========================================
RYA-485 Issue 2 — VERIFY (do not assume) how solar O is modeled. Correctness-critical:
the solar O is the differential denominator, so a silent-LTE-where-3D-is-required error
would poison every downstream O (Procyon, alpha Cen, 55 Cnc). DIAGNOSE-ONLY pass — this
traces the radiative-transfer regime actually applied and reports; it changes nothing.

What it checks, from the live code (pipeline.nlte_cno), not assumptions:
  1. O I 777 leg selection by Teff: the Amarsi-2019 grid 3D leg (table3) reaches the
     ~6500 K STAGGER ceiling, so the SUN (5772 K) must get 3D-NLTE, Procyon (6554 K) the
     1D-NLTE leg. Confirm the Sun actually gets 3D (not silently 1D).
  2. [O I] 6300 (630.0 nm): the grid HAS a 630 nm node — confirm its 3D differential at
     solar params is ~0 (so the production 'lte_forbidden_insensitive' treatment is the
     CORRECT method, not a silent skip of a real 3D structural term — RYA-447/448).
  3. The differential REGIME-MATCH: the Procyon [O/H] differences Procyon-1D-NLTE against
     the Sun. If the Sun's denominator is 3D-NLTE while Procyon is 1D-NLTE, the differential
     mixes regimes — report the magnitude (the solar 3D-minus-1D term at 777).

Output: data/results/solar_o_modeling_verification_rya485.json + a printed report.

    python -m scripts.verify_solar_o_modeling_rya485
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline import nlte_cno as N                                # noqa: E402
from config.constants import get_star_params                     # noqa: E402

RESULTS = ROOT / 'data' / 'results'
SOLAR_OURS_O = 8.735
# solar O I 777 1D-LTE local-continuum anchor (RYA-478 solar control) — the A_lte the
# NLTE delta is applied at, for a representative regime-difference number.
SOLAR_OI777_ALTE = 8.907


def _params(star):
    r = get_star_params(star)
    return (float(r['teff']), float(r['logg']), float(r['feh_ref']), float(r.get('xi', 1.0)))


def main():
    rep = {'ticket': 'RYA-485 Issue 2', 'mode': 'diagnose-only'}

    # ── 1. O I 777 leg selection: Sun must be 3D, Procyon 1D ───────────────────
    print("=" * 78)
    print("  RYA-485 Issue 2 — solar O radiative-transfer verification (diagnose-only)")
    print("=" * 78)
    print(f"\n  [1] O I 777 leg selection (3D ceiling = {N.TEFF_3D_CEILING:.0f} K)")
    legs = {}
    for star in ('solar', 'procyon'):
        te, lg, fe, xi = _params(star)
        leg = N.select_leg(te)
        lab = N.resolve_line('OI', 7773.0)
        d_auto = N.cno_nlte_delta('OI', lab, te, lg, fe, xi, SOLAR_OI777_ALTE if star == 'solar' else 9.36)
        legs[star] = {'teff': te, 'leg': leg, 'delta_auto': round(float(d_auto), 4)}
        print(f"      {star:8} Teff={te:.0f} -> leg={leg}  O I 777 δ(auto)={d_auto:+.4f}")
    sun_gets_3d = legs['solar']['leg'] == '3D' and np.isfinite(legs['solar']['delta_auto'])
    print(f"      => Sun gets 3D-NLTE on 777: {'YES ✓ (correct, not silent-1D)' if sun_gets_3d else 'NO — FINDING'}")
    rep['oi777_leg'] = legs
    rep['sun_gets_3d_on_777'] = bool(sun_gets_3d)

    # ── 2. [O I] 6300 3D differential ~0 → forbidden-LTE is the correct method ──
    print(f"\n  [2] [O I] 6300 (630.0 nm) 3D differential at solar params")
    te, lg, fe, xi = _params('solar')
    d6300_3d = N.cno_nlte_delta('OI', '630.0nm', te, lg, fe, xi, 8.73, leg='3D')
    d6300_1d = N.cno_nlte_delta('OI', '630.0nm', te, lg, fe, xi, 8.73, leg='1D')
    forbidden_lte_ok = abs(float(d6300_3d)) < 0.01
    print(f"      630 nm 3D δ={d6300_3d:+.4f}  1D δ={d6300_1d:+.4f}  "
          f"-> forbidden-LTE treatment correct: {'YES ✓ (3D term ~0)' if forbidden_lte_ok else 'NO — 3D term non-negligible, FINDING'}")
    rep['oi6300_3d_delta'] = round(float(d6300_3d), 4)
    rep['forbidden_lte_is_correct'] = bool(forbidden_lte_ok)

    # ── 3. differential regime-match: Procyon-1D vs Sun-3D ─────────────────────
    print(f"\n  [3] differential regime-match (Procyon-1D-NLTE vs Sun-?-NLTE)")
    sun_d3 = N.cno_nlte_delta('OI', '777nm', te, lg, fe, xi, SOLAR_OI777_ALTE, leg='3D')
    sun_d1 = N.cno_nlte_delta('OI', '777nm', te, lg, fe, xi, SOLAR_OI777_ALTE, leg='1D')
    sun_3d = round(SOLAR_OI777_ALTE + float(sun_d3), 3)
    sun_1d = round(SOLAR_OI777_ALTE + float(sun_d1), 3)
    regime_term = round(sun_1d - sun_3d, 3)            # how much the denominator moves 3D->1D
    procyon_nlte = 8.82                                # RYA-483 Procyon O I 777 1D-NLTE
    oh_vs_sun3d = round(procyon_nlte - sun_3d, 3)
    oh_vs_sun1d = round(procyon_nlte - sun_1d, 3)
    print(f"      Sun O I 777: 3D-NLTE {sun_3d} | 1D-NLTE {sun_1d}  (3D−1D regime term {regime_term:+.3f})")
    print(f"      Procyon [O/H]: vs Sun-3D {oh_vs_sun3d:+.3f} (current) | vs Sun-1D {oh_vs_sun1d:+.3f} (regime-matched)")
    print(f"      => the Procyon differential mixes Procyon-1D with Sun-3D; regime-matching")
    print(f"         it (Sun-1D denominator) shifts [O/H] by {oh_vs_sun1d - oh_vs_sun3d:+.3f} dex — a FINDING to fix.")
    rep['regime_mismatch'] = {
        'sun_oi777_3d_nlte': sun_3d, 'sun_oi777_1d_nlte': sun_1d,
        'regime_term_3d_minus_1d': regime_term,
        'procyon_oh_vs_sun3d': oh_vs_sun3d, 'procyon_oh_vs_sun1d': oh_vs_sun1d,
        'shift_if_regime_matched': round(oh_vs_sun1d - oh_vs_sun3d, 3)}

    # ── 4. per-arm solar O spread (from committed RYA-460 production values) ────
    print(f"\n  [4] per-arm SOLAR O spread (production, RYA-460 Kitt Peak)")
    # RYA-460 raw (LTE) KP values; O I 777 carries the 3D-NLTE δ, [O I] 6300 is forbidden-LTE
    kp_777_raw_lte, kp_6300_raw = 8.955, 8.835
    sun_777_nlte = round(kp_777_raw_lte + float(sun_d3), 3)   # apply 3D-NLTE to the 777 LTE
    spread = round(abs(sun_777_nlte - kp_6300_raw), 3)
    print(f"      O I 777 (KP, 3D-NLTE) ≈ {sun_777_nlte}  |  [O I] 6300 (KP, forbidden-LTE) {kp_6300_raw}")
    print(f"      solar O internal per-arm spread ≈ {spread} dex  (Procyon's continuum lever was ~0.18)")
    print(f"      => solar O spread does {'NOT ' if spread < 0.12 else ''}mirror Procyon's 0.18 -> "
          f"{'Procyon continuum lever is UVES-specific, not an upstream solar bug ✓' if spread < 0.12 else 'investigate upstream'}")
    rep['per_arm_solar_O'] = {'oi777_3d_nlte': sun_777_nlte, 'oi6300_forbidden_lte': kp_6300_raw,
                              'spread': spread, 'procyon_continuum_lever': 0.18,
                              'mirrors_procyon': bool(spread >= 0.12)}

    # ── verdict ────────────────────────────────────────────────────────────────
    print("\n" + "=" * 78)
    no_silent_lte = sun_gets_3d and forbidden_lte_ok
    print(f"  VERDICT: solar O modeling is {'CORRECT — no silent-LTE' if no_silent_lte else 'has a SILENT-LTE FINDING'}.")
    print(f"  Findings to fix (named, NOT applied this pass): differential regime-mismatch "
          f"({rep['regime_mismatch']['shift_if_regime_matched']:+.3f} dex on Procyon [O/H]).")
    rep['no_silent_lte'] = bool(no_silent_lte)
    rep['findings'] = [
        'differential regime-mismatch: Procyon-1D-NLTE differenced against Sun-3D-NLTE; '
        f"regime-matched (Sun-1D denominator) shifts Procyon [O/H] by "
        f"{rep['regime_mismatch']['shift_if_regime_matched']:+.3f} dex"]
    RESULTS.mkdir(parents=True, exist_ok=True)
    (RESULTS / 'solar_o_modeling_verification_rya485.json').write_text(json.dumps(rep, indent=2, default=str))
    print(f"\n  [out] {RESULTS / 'solar_o_modeling_verification_rya485.json'}")
    return rep


if __name__ == '__main__':
    main()
