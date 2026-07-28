"""
scripts/procyon_uves_oi_phase2_rya348.py
========================================
RYA-348 Phase 2 — Procyon UVES O I 777 PRIMARY-O arm.

Phase 1 (PR #86) shook down HARPS VIS and confirmed HARPS-only Procyon C/O is not
science-grade: the only O indicator in HARPS range is [O I] 6300, which came out
cited-anchor-dominated (Caffau 2015 transplanted, χ²ᵣ≫10 — the RYA-104 recurrence).
This phase replaces that with a REAL Procyon O: the O I 7771-5 triplet on a staged,
telluric-CLEAN UVES RED760 frame (RYA-272 loader + RYA-271 audit), resolved through
the RYA-464 per-star arm registry (never silent-Vesta).

Deliverable:
  * A(O) from UVES O I 777, 1D-LTE and 1D-NLTE (Amarsi 2019 grid, self-consistent
    at the fitted A_lte — NOT a value from memory), differenced against OUR measured
    solar O (8.735; solar Phase C, RYA-469 frozen gold).
  * Cross-instrument leg: UVES O I 777 vs HARPS [O I] 6300 — agreement within error or
    a named finding; per-arm zero-point reported, NEVER averaged. (UVES also carries
    [O I] 6300, so the same line is compared across instruments as a bonus zero-point.)

Scope: analysis/scratch only (reads spectra + line list, runs synthesis, writes to
data/results/ + data/audit/). No STAR_PARAMS / line-list / EW / spectra edits. STOP at
the verdict, no merge.

Usage:  python scripts/procyon_uves_oi_phase2_rya348.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import pipeline.cno_synthesis as cs                                  # noqa: E402
from config.constants import get_star_params, SOLAR_ASPLUND2021      # noqa: E402

# Our own measured Sun — the differential denominator (solar Phase C; RYA-469 frozen gold).
SOLAR_OURS_O = 8.735
PROD_DIR = ROOT / 'data' / 'audit' / 'cno_synthesis'
RESULTS = ROOT / 'data' / 'results'
TMP = '/tmp/ispec_cno_rya348_uves'


def _rule(c='='):
    print(c * 78)


def _harps_oi6300():
    """Read the HARPS [O I] 6300 result from the Phase-1 VIS product (run on current main)
    for the leg-validation. Returns (A_raw_lte, A_cited, chi2r) or None."""
    pb = PROD_DIR / 'procyon_vis_cno_per_band.csv'
    if not pb.exists():
        return None
    df = pd.read_csv(pb)
    row = df[df['key'] == 'OI_6300']
    if row.empty:
        return None
    r = row.iloc[0]
    from pipeline.cno_synthesis import apply_cited_corrections, REGIONS
    rec = get_star_params('procyon')
    params = {'teff_K': float(rec['teff']), 'logg': float(rec['logg']),
              'feh': float(rec['feh_ref']), 'vturb_kms': float(rec.get('xi', 1.0))}
    corr = {c['key']: c for c in apply_cited_corrections(df.to_dict('records'), params, REGIONS['vis'])}
    return (float(r['A_X']), float(corr['OI_6300'].get('a_corr')), float(r['red_chi2']))


def main():
    rec = get_star_params('procyon')
    teff, logg, feh = float(rec['teff']), float(rec['logg']), float(rec['feh_ref'])
    xi = float(rec.get('xi', 1.0))
    _rule()
    print(f"  RYA-348 Phase 2 — Procyon UVES O I 777 primary-O  "
          f"(Teff={teff:.0f} logg={logg:.2f} [Fe/H]={feh:+.2f} xi={xi:.1f})")
    _rule()

    # ── resolve the UVES arm through the registry (RYA-464; never silent-Vesta) ──
    registry = cs.star_arm_registry('procyon')
    arm = registry['uves']
    print(f"  arm 'uves': ready={arm.ready}  loader={arm.loader}  provenance={arm.provenance}")
    if not arm.ready:
        raise SystemExit(f"UVES arm not ready: {arm.defer_reason}")
    obs_w, obs_f = cs.resolve_arm_spectrum('procyon', arm)   # loud-fails if it would be Vesta
    print(f"  resolved UVES spectrum: {obs_w.size:,} px, {obs_w.min()*10:.0f}-{obs_w.max()*10:.0f} A")

    # ── synthesis resources (same machinery as run_phase_a) ──
    params = {'teff_K': teff, 'logg': logg, 'feh': feh, 'vturb_kms': xi}
    atm = cs._load_atmosphere(teff, logg, feh, xi)
    ll, iso, chem = cs._load_synth_resources()
    sab = cs.ispec.read_solar_abundances(cs._ISPEC_SOLAR_ABUND_FILE)
    codes = cs._atom_codes(('C', 'N', 'O', 'Ni'), chem, sab)
    _, vmac, vsini, _ = cs._resolve_broadening('procyon')

    # pin A(C) from the Phase-1 HARPS run (for the [O I]/CO coupling); else solar anchor
    fixed = {'Ni': 6.20, 'C': SOLAR_ASPLUND2021['C'], 'N': SOLAR_ASPLUND2021['N'],
             'O': SOLAR_ASPLUND2021['O']}
    vis_prod = PROD_DIR / 'procyon_vis_cno_product.csv'
    if vis_prod.exists():
        vp = pd.read_csv(vis_prod)
        crow = vp[vp['element'] == 'C']
        if not crow.empty:
            fixed['C'] = float(crow.iloc[0]['A_X'])
            print(f"  pinned A(C)={fixed['C']:.3f} from Phase-1 HARPS product (for [O I]/CO coupling)")

    # ── fit the in-coverage O diagnostics: O I 777 (primary) + [O I] 6300 (UVES cross-check) ──
    diags = tuple(d for d in cs.PROCYON_UVES_DIAGNOSTICS if d.key in ('OI_777', 'OI_6300'))
    _rule('-')
    print(f"  fitting {[d.key for d in diags]} on UVES (RED760 covers both); "
          f"C I 5052/5380 are below RED760 5655 A — not fit here")
    _rule('-')
    Path(TMP).mkdir(parents=True, exist_ok=True)
    per_band, corrections = cs._fit_arm(cs.UVES_OPT, diags, obs_w, obs_f, params, fixed,
                                        atm, ll, iso, sab, codes,
                                        (cs.UVES_OPT.R, vmac, vsini), TMP)
    corr_by = {c['key']: c for c in corrections}
    pb_by = {r['key']: r for r in per_band}

    # ── results table ──
    _rule()
    print("  PER-INDICATOR A(O) — UVES (raw 1D-LTE + cited Amarsi-2019 1D-NLTE / Caffau anchor)")
    _rule()
    print(f"  {'key':10s} {'role':11s} {'A_lte':>7s} {'A_corr':>7s} {'delta':>7s} "
          f"{'chi2r':>7s} {'kind':>20s}")
    rows = []
    for d in diags:
        pb = pb_by[d.key]
        c = corr_by[d.key]
        a_lte = pb.get('A_X')
        a_corr = c.get('a_corr')
        delta = c.get('delta')
        chi2 = pb.get('red_chi2')
        kind = c.get('kind')
        print(f"  {d.key:10s} {d.role:11s} {a_lte:7.3f} "
              f"{(a_corr if a_corr is not None else float('nan')):7.3f} "
              f"{(delta if delta is not None else float('nan')):+7.3f} "
              f"{chi2:7.2f} {str(kind):>20s}")
        rows.append(dict(key=d.key, role=d.role, A_lte=round(float(a_lte), 3),
                         A_corr=(round(float(a_corr), 3) if a_corr is not None else None),
                         delta=(round(float(delta), 3) if delta is not None else None),
                         chi2r=round(float(chi2), 2), corr_kind=kind,
                         flag=c.get('flag'), source=c.get('source')))

    # ── the primary-O number + differential vs our Sun ──
    oi777 = next(r for r in rows if r['key'] == 'OI_777')
    a_lte = oi777['A_lte']
    a_nlte = oi777['A_corr'] if oi777['A_corr'] is not None else a_lte
    _rule()
    print("  PRIMARY O — UVES O I 777 (Amarsi 2019 1D-NLTE, self-consistent at fitted A_lte)")
    _rule()
    print(f"  A(O) 1D-LTE              = {a_lte:.3f}")
    print(f"  A(O) 1D-NLTE            = {a_nlte:.3f}   (Δ_NLTE = {a_nlte - a_lte:+.3f}, {oi777['corr_kind']})")
    print(f"  our measured solar O    = {SOLAR_OURS_O:.3f}  (solar Phase C, RYA-469 gold)")
    print(f"  [O/H] differential      = {a_nlte - SOLAR_OURS_O:+.3f}  (NLTE, vs our Sun)")
    print(f"  χ²ᵣ (O I 777 fit)        = {oi777['chi2r']:.2f}  "
          f"[{'CLEAN' if oi777['chi2r'] < 10 else 'POOR-FIT'}]")

    # ── leg validation: UVES O I 777 vs HARPS [O I] 6300; + UVES vs HARPS on [O I] 6300 ──
    _rule()
    print("  LEG VALIDATION — UVES O I 777 (primary) vs HARPS [O I] 6300 (cross-check)")
    _rule()
    leg = {'uves_oi777_nlte': a_nlte, 'uves_oi777_lte': a_lte}
    harps = _harps_oi6300()
    uves_oi6300 = next((r for r in rows if r['key'] == 'OI_6300'), None)
    if harps is not None:
        h_lte, h_cited, h_chi = harps
        leg.update(harps_oi6300_lte=round(h_lte, 3), harps_oi6300_cited=round(h_cited, 3),
                   harps_oi6300_chi2r=round(h_chi, 1))
        zp = a_nlte - h_cited
        print(f"  UVES O I 777 (NLTE)      = {a_nlte:.3f}")
        print(f"  HARPS [O I] 6300 (cited) = {h_cited:.3f}  (raw χ²ᵣ {h_chi:.1f} — Caffau-anchor "
              f"dominated, RYA-104; NOT an independent Procyon O)")
        print(f"  per-arm zero-point Δ     = {zp:+.3f} dex (UVES O I 777 minus HARPS [O I] 6300)")
        verdict = ('AGREE within ~0.1 dex' if abs(zp) <= 0.10 else
                   'DISAGREE — surfaced as a finding, NOT averaged (the HARPS leg is cited-anchor)')
        print(f"  leg verdict              = {verdict}")
        leg['zeropoint_uves777_minus_harps6300'] = round(zp, 3)
        leg['verdict'] = verdict
    else:
        print("  HARPS [O I] 6300 product absent — run the VIS shakedown first for the leg.")
    if uves_oi6300 is not None:
        same = uves_oi6300['A_corr'] if uves_oi6300['A_corr'] is not None else uves_oi6300['A_lte']
        print(f"\n  BONUS same-line zero-point — [O I] 6300 UVES vs HARPS:")
        print(f"    UVES [O I] 6300 (cited)  = {same:.3f}  χ²ᵣ {uves_oi6300['chi2r']:.1f}")
        if harps is not None:
            print(f"    HARPS [O I] 6300 (cited) = {harps[1]:.3f}  → Δ {same - harps[1]:+.3f} "
                  f"(both Caffau-anchored → not independent; sanity only)")

    # ── write results ──
    RESULTS.mkdir(parents=True, exist_ok=True)
    out = {
        'ticket': 'RYA-348 Phase 2 — Procyon UVES O I 777 primary-O',
        'star': 'procyon', 'params': params,
        'arm': {'instrument': 'UVES', 'anchor_epoch': '2013-10-08 RED760 (oi_anchor, CLEAN)',
                'resolved_via': 'resolve_arm_spectrum (RYA-464, never silent-Vesta)'},
        'indicators': rows,
        'primary_O': {'A_lte': a_lte, 'A_nlte': a_nlte, 'nlte_delta': round(a_nlte - a_lte, 3),
                      'solar_O_ours': SOLAR_OURS_O, 'OH_differential_nlte': round(a_nlte - SOLAR_OURS_O, 3),
                      'chi2r': oi777['chi2r']},
        'leg_validation': leg,
    }
    (RESULTS / 'procyon_uves_oi777_phase2.json').write_text(json.dumps(out, indent=2, default=str))
    pd.DataFrame(rows).to_csv(RESULTS / 'procyon_uves_oi777_phase2.csv', index=False)
    print(f"\n  [out] {RESULTS / 'procyon_uves_oi777_phase2.json'}")
    print(f"  [out] {RESULTS / 'procyon_uves_oi777_phase2.csv'}")
    return out


if __name__ == '__main__':
    main()
