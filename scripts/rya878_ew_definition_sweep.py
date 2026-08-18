#!/usr/bin/env python3
"""
RYA-878 — is the +0.152 dex EW offset DEFINITIONAL? Sweep the integration range.
===============================================================================
    ISPEC_DIR=... python3 scripts/rya878_ew_definition_sweep.py

The like-for-like measurement settled WHOSE the offset is: the adapter reproduces the
engine's synthetic EW to +0.0000 dex (18/18 within 10%), and the ENGINE sits +0.1523 dex
above the observed EW. So it is a MODEL-vs-OBSERVATION property, not an adapter defect,
and it SURVIVES pairing — unlike ANGLE 2, which dissolved into a line-set artifact.

WHAT IS LEFT, AFTER EIGHT REFUTATIONS
--------------------------------------
Refuted in RYA-873/875 and NOT re-tested here: same-species contamination; a uniform
continuum offset (log10 std 0.247); saturation (-0.047 vs reduced EW); line strength
(-0.100); red_chi2 (+0.059); sigma_A (-0.169).

Refuted in RYA-878 by inspection rather than by measurement, because the structure
settles them:

  * RYA-771's %.2f trial-abundance staircase — 0.01 dex is ~2.3% in EW for a weak line.
    It cannot produce a factor of 1.42. Refuted by MAGNITUDE.
  * RYA-770's residual edge pixels — `synth_ew.synthetic_ew_mA` builds its OWN wavelength
    grid and synthesises on it. It never receives the observed spectrum, so a zeroed
    observed edge pixel cannot reach it. Refuted STRUCTURALLY.

That leaves the one the numbers keep pointing at: THE TWO EWs ARE NOT THE SAME QUANTITY.

  observed EW  = ONE fitted Voigt, integrated       (profile_fit: FIT_HALF_A = 0.60 A)
  synthetic EW = ALL of the element's absorption    (synth_ew:    >= 1.0 A half-width)
                 in the window, damping wings and
                 same-species neighbours included

THE TEST
--------
Sweep the synthetic integration half-width and ask where the ratio crosses 1. If it lands
near the profile fitter's own 0.60 A extent, the offset IS the definition and the RCA is
closed with a named cause. If the ratio stays high even at the profile fitter's extent,
the definition is NOT the cause and something in the model genuinely over-absorbs — a
much more interesting finding, and the one this sweep is built to be able to return.

⚠️ NOTHING HERE CHANGES A DEFINITION. Production keeps the derived range; the override
exists so the range can be interrogated. No abundance, bar or harness residual is touched.
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
sys.path.insert(0, str(ROOT / "scripts"))

from pipeline import synth_ew                                      # noqa: E402

OUT = ROOT / "data" / "results" / "rya878"

#: Integration half-widths to sweep, in A. 0.60 is `profile_fit.FIT_HALF_A`, the extent
#: the observed EW is fitted over; 1.0 is `synth_ew`'s floor, i.e. today's definition.
HALF_WIDTHS_A = (0.30, 0.45, 0.60, 0.80, 1.00, 1.40)
PROFILE_FIT_HALF_A = 0.60


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--star", default="solar")
    ap.add_argument("--element", default="Fe")
    ap.add_argument("--ion", default="I")
    ap.add_argument("--out", type=Path, default=OUT)
    a = ap.parse_args()

    per_line = a.out / "rya878_angle1_by_line.csv"
    if not per_line.exists():
        raise SystemExit(f"run rya878_angle1_like_for_like.py first: {per_line} absent")
    m = pd.read_csv(per_line)

    from control_synthesis_handler import build_context
    from pipeline.measure import resolve_handler
    from pipeline.band_policy import resolve as resolve_band
    cat = pd.read_csv(ROOT / "data" / "catalog" / "instrument_catalog.csv")
    R = float(cat[cat.iloc[:, 0].astype(str) == "harps"].iloc[0]["resolving_power_max"])
    ctx = build_context(a.element, a.ion, R)
    handler = resolve_handler(3400.0)
    handler.prepare(resolve_band(5000.0), ctx)
    kw = handler._synth_kw(ctx, a.element)

    rows = []
    for r in m.itertuples():
        c = float(r.wave)
        for hw in HALF_WIDTHS_A:
            ew = synth_ew.synthetic_ew_mA(
                handler._synth, centre_A=c, abundance=float(r.a_engine),
                fit_half_width_A=1.0, synth_kwargs=kw, integration_half_width_A=hw)
            rows.append({"wave": c, "half_width_A": hw, "ew_synth": ew,
                         "ew_obs": float(r.ew_obs), "ratio": ew / float(r.ew_obs)})
        print(f"  {c:10.3f} done")
    d = pd.DataFrame(rows)

    print(f"\n=== synthetic/observed EW vs INTEGRATION HALF-WIDTH ===")
    print(f"{'half_A':>8}{'median':>10}{'dex':>10}{'MAD':>9}{'min':>8}{'max':>9}"
          f"{'within10%':>11}")
    summ = []
    for hw, g in d.groupby("half_width_A"):
        r = g.ratio[np.isfinite(g.ratio)]
        s = {"half_width_A": hw, "n": int(r.size),
             "median_ratio": round(float(r.median()), 4),
             "dex": round(float(np.log10(r.median())), 4),
             "mad": round(float(np.median(np.abs(r - r.median()))), 4),
             "min": round(float(r.min()), 3), "max": round(float(r.max()), 3),
             "within_10pct": int(np.sum(np.abs(r - 1.0) <= 0.10))}
        summ.append(s)
        mark = "   <- profile fitter's own extent" if hw == PROFILE_FIT_HALF_A else ""
        print(f"{hw:>8.2f}{s['median_ratio']:>10.3f}{s['dex']:>+10.4f}{s['mad']:>9.3f}"
              f"{s['min']:>8.2f}{s['max']:>9.2f}{s['within_10pct']:>8}/{s['n']}{mark}")

    at_fit = next(s for s in summ if s["half_width_A"] == PROFILE_FIT_HALF_A)
    definitional = abs(at_fit["dex"]) <= 0.05
    print(f"\n=== verdict ===")
    if definitional:
        print(f"  DEFINITIONAL. At the profile fitter's own +/-{PROFILE_FIT_HALF_A} A the")
        print(f"  ratio is {at_fit['median_ratio']:.3f} ({at_fit['dex']:+.4f} dex) — the two")
        print(f"  EWs agree once they cover the same extent, so the +0.152 is what the")
        print(f"  wider synthetic range includes and the Voigt fit does not. ANGLE 1 is")
        print(f"  comparing two different definitions; it is not a model defect.")
    else:
        print(f"  🔴 NOT DEFINITIONAL. Even at the profile fitter's own extent the model")
        print(f"  is {at_fit['dex']:+.4f} dex above the observed EW ({at_fit['median_ratio']:.3f}).")
        print(f"  The range is not the cause, so the model genuinely over-absorbs at the")
        print(f"  fitted abundance. That is a real finding about the synthesis and the")
        print(f"  remaining live root — it needs its own investigation.")

    a.out.mkdir(parents=True, exist_ok=True)
    d.to_csv(a.out / "rya878_ew_width_sweep.csv", index=False)
    (a.out / "rya878_ew_width_summary.json").write_text(json.dumps({
        "ticket": "RYA-878",
        "profile_fit_half_width_A": PROFILE_FIT_HALF_A,
        "production_half_width_A": ">= 1.0 (synth_ew.ew_half_width_A)",
        "sweep": summ,
        "definitional": bool(definitional),
        "refuted_by_measurement_873_875": [
            "same-species contamination", "uniform continuum offset (log10 std 0.247)",
            "saturation (-0.047 vs reduced EW)", "line strength (-0.100)",
            "red_chi2 (+0.059)", "sigma_A (-0.169)"],
        "refuted_by_structure_878": {
            "RYA-771 %.2f abundance quantisation":
                "0.01 dex is ~2.3% in EW; cannot produce a factor of 1.42 (MAGNITUDE)",
            "RYA-770 residual edge pixels":
                "synth_ew builds its own grid and never receives the observed spectrum, "
                "so a zeroed observed pixel cannot reach it (STRUCTURAL)"},
    }, indent=2) + "\n")
    print(f"\n  wrote {a.out.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
