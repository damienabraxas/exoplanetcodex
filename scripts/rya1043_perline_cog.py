#!/usr/bin/env python3
"""Per-line synthetic COG linearity test — the ratified `rew_max` gate (RYA-1041).

WHY PER LINE. An ENSEMBLE curve of growth for solar Fe I has no resolvable knee, and
that is physics rather than method (RYA-1041, ratified 2026-08-25): damping, blending
and core depth vary line to line, so at fixed reduced-EW the ensemble carries intrinsic
spread that smears the transition. Measured: the usable pool IS a COG (1805 lines,
corr +0.362, medians monotone -5.22 -> -4.77) and the local dREW/dX still has THREE sign
changes with no linear-then-flat structure. A single-line COG is sharp; an ensemble COG
is broad BY CONSTRUCTION.

⚠️ SO DO NOT RE-RUN THE ENSEMBLE REGRESSION. That question is answered.

CONSEQUENCE: saturation is a PER-LINE property. Every line saturates at its own reduced
EW, so `rew_max` is not a shared scalar to measure once -- the gate is "has THIS line's
own response flattened".

THE TEST. For each line, synthesise the spectrum at A-delta, A and A+delta, integrate EW
from each synthetic, and read the response directly:

    dREW/dA  =  [ log10(EW(A+d)/lam) - log10(EW(A-d)/lam) ] / (2d)

On the linear part of the curve EW is proportional to the number of absorbers and this
approaches 1. As the core saturates it collapses toward zero. The line is off the linear
COG when ITS OWN derivative has flattened -- physics per line, no regression, no quota.

The FLATTENING THRESHOLD is deliberately left as a parameter and printed with the
distribution rather than baked in: RYA-1041 Step 2 is Ryan's, and this script produces
the measurement he sets it from.

COST: three syntheses per line. Sirius-only (the grids single-source there, RYA-567) and
intended to run in the BACKGROUND while other work proceeds.
"""
from __future__ import annotations

import argparse
import glob
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

OUT = ROOT / "data" / "audit" / "rya1041"


def build_pool(lo: float, hi: float) -> pd.DataFrame:
    """in_aggregate + REW-saturation-only exclusions — clean AND untruncated.

    Neither half alone can answer this. The derived product is clean but cut at
    REW_SATURATION_CEILING = -4.9, which removes the saturated lines the test is ABOUT.
    The raw candidate set spans the knee but does not trace a COG at all (corr -0.015):
    it is full of blends where the measured EW tracks whatever sits at that wavelength.
    The stage-1 pool records ceiling-only exclusions explicitly -- "Measured and
    reported; excluded from the aggregate ONLY" -- so those lines passed every other
    quality gate and are recoverable without re-measuring anything.
    """
    files = sorted(glob.glob(str(ROOT / "data" / "measured" / "band_ew"
                                 / f"FeI_{int(lo)}_{int(hi)}_*_PROFILEFIT_ew.csv")))
    if not files:
        raise SystemExit(f"no stage-1 EW pool for {lo}-{hi} A. This test needs the "
                         f"UNTRUNCATED pool; the derived product is cut at -4.9 and "
                         f"cannot answer the question (RYA-1041).")
    d = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
    er = d.excluded_reason.astype(str)
    sat_only = er.str.startswith("REW ") & er.str.contains("saturation ceiling")
    d = d[(d["in_aggregate"] == True) | sat_only].copy()
    # NB: no leading underscore — itertuples() renames such columns positionally
    d["sat_only"] = sat_only.reindex(d.index, fill_value=False)
    # one row per line: the arms measure the same lines, and the derivative is a property
    # of the LINE and the model, not of which spectrum we happened to look at
    d["w"] = pd.to_numeric(d.wavelength_air_A, errors="coerce")
    d = d.dropna(subset=["w"]).sort_values("w")
    return d.groupby("w", as_index=False).first()


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--lo", type=float, default=4200)
    ap.add_argument("--hi", type=float, default=6910)
    ap.add_argument("--delta", type=float, default=0.10,
                    help="abundance step in dex for the two-sided derivative. 🔴 DECLARED, "
                         "NOT BORROWED (RYA-161/1043): measured on 10 lines across a 20x "
                         "range, the derivative is stable to 0.4% -- 0.02 -> +0.4224, "
                         "0.05 -> +0.4224, 0.10 -> +0.4225, 0.20 -> +0.4228, "
                         "0.40 -> +0.4241. 0.10 sits mid-plateau; drift begins only at "
                         "0.40, the expected onset of leaving the local linear regime. "
                         "The bracket is wide enough to resolve the slope and narrow "
                         "enough not to traverse the knee.")
    ap.add_argument("--half-width-A", type=float, default=0.62,
                    help="⚠️ the window integrates ALL absorption in it, so a blended "
                         "neighbour inflates EW. The DERIVATIVE is robust to that (only "
                         "Fe's abundance varies between the three syntheses) but a "
                         "blend DILUTES the response, so a blended line reads as more "
                         "saturated than it is. Reported, not corrected.")
    ap.add_argument("--limit", type=int, default=None, help="first N lines (smoke test)")
    ap.add_argument("--blend-offset", type=float, default=3.0,
                    help="dex to drive Fe DOWN by for the blend-only baseline. Large "
                         "enough that the Fe line is negligible, so what remains in the "
                         "window is everything that is not this line.")
    ap.add_argument("--a0", type=float, default=None,
                    help="the A(Fe) to evaluate dREW/dA AT. Default: the median fitted "
                         "abundance of the pool's own measurable lines.")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    from config.constants import assert_on_sirius
    assert_on_sirius("RYA-1041 per-line synthetic COG")     # grids single-source (RYA-567)

    from pipeline.abundances_derive import _synth_flux_at_abund
    from pipeline.band_products import equivalent_width
    # 🔴 THE PRODUCTION CONTEXT BUILDER, not a hand-rolled one. `build_solar_context`
    # sources every value from an existing loader (stellar params, atmosphere, broadening,
    # linelist) precisely so a synthesis here cannot silently differ from the one that
    # produced the products this test is about.
    from pipeline.nearuv_synth import build_solar_context
    from pipeline import band_policy
    from config.synth_bands import SYNTH_BANDS

    pool = build_pool(a.lo, a.hi)
    if a.limit:
        pool = pool.head(a.limit)
    print(f"[rya1041] {len(pool)} distinct lines in {a.lo:.0f}-{a.hi:.0f} A "
          f"({int(pool.sat_only.sum())} recovered from the -4.9 ceiling)")

    pol = band_policy.resolve(0.5 * (a.lo + a.hi))
    cfg = SYNTH_BANDS.get(pol.name)
    if cfg is None or not cfg.linelist.exists():
        raise SystemExit(f"no synthesis linelist for band {pol.name!r}; the near-UV and "
                         f"IR lists are BUILD STEPS, not repo files (RYA-1027).")
    print(f"[rya1041] band {pol.name}  linelist {cfg.linelist.name}")
    ctx = build_solar_context("Fe", 500000.0, linelist_file=str(cfg.linelist),
                              star="solar")

    # 🔴 ONE A0 FOR EVERY LINE, and it must not be the line's OWN fitted abundance.
    # Two reasons, the first fatal:
    #   * the lines this test is ABOUT have no fitted abundance. They are the
    #     ceiling-excluded ones, and they were excluded precisely because "the
    #     EW->abundance inversion is ill-conditioned here" -- so A is NaN for exactly
    #     the saturated lines. My first draft read `r.abundance` and returned
    #     `no_abundance` for 3 of 3 in the smoke test.
    #   * it is also the wrong question. A curve of growth is a statement at FIXED
    #     composition: "has THIS line saturated at the Sun's iron abundance". Taking the
    #     derivative at each line's own fitted A would evaluate every line at a different
    #     composition and confound saturation with the fit's own error.
    _fit = pd.to_numeric(pool.get("abundance"), errors="coerce").dropna()
    if a.a0 is not None:
        A0 = float(a.a0)
        _a0_src = "--a0 (caller-supplied)"
    elif len(_fit):
        A0 = float(_fit.median())
        _a0_src = f"median of {len(_fit)} measurable lines in this pool"
    else:
        raise SystemExit("no abundance available to evaluate the derivative at; pass --a0")
    print(f"[rya1041] evaluating dREW/dA at A(Fe) = {A0:.4f}  [{_a0_src}]")
    rows, t0 = [], time.time()
    for i, r in enumerate(pool.itertuples(), 1):
        lam = float(r.w)
        w_nm = np.arange(lam - a.half_width_A, lam + a.half_width_A, 0.001) / 10.0
        ews = {}
        try:
            # 🔴 FOURTH SYNTHESIS: THE BLEND-ONLY BASELINE. Integrating (1-f) over the
            # window captures EVERY absorber in it, and a blend contributes a roughly
            # CONSTANT EW that does not move with Fe:
            #     dREW/dA = (1/ln10) * (dEW_line/dA) / (EW_line + EW_blend)
            # so the blend inflates the denominator and DEPRESSES the derivative, worst
            # exactly where EW_line is smallest. MEASURED in the first full run: the six
            # weakest lines returned dREW/dA of 0.014 / 0.017 / 0.079 / 0.112 / -0.024 /
            # 0.040 when weak lines should be the MOST linear, with ew_mid piled at
            # 13.3-14.1 mA — the blend floor, not the line.
            # Synthesising with Fe driven far down leaves the blend alone, so subtracting
            # it recovers the LINE's own EW. One extra synthesis per line (+33%).
            for tag, A in (("blend", A0 - a.blend_offset),
                           ("lo", A0 - a.delta), ("mid", A0), ("hi", A0 + a.delta)):
                f = _synth_flux_at_abund(w_nm, ctx["atmosphere"], ctx["teff"], ctx["logg"],
                                         ctx["feh"], ctx["vturb"], ctx["linelist"],
                                         ctx["isotopes"], ctx["solar_abund"],
                                         "Fe", ctx["atom_code"], trial_A=float(A))
                # 🔴 INTEGRATE DIRECTLY; do NOT route a SYNTHETIC through
                # `band_products.equivalent_width`. That function exists for OBSERVED
                # spectra: it places (or verifies) a local continuum from side-bands, and
                # on a 1.24 A synthetic window there are none -- it refused all three
                # smoke-test lines with "too few side-band points (0)". A normalised
                # synthesis HAS no continuum problem: the continuum is unity BY
                # CONSTRUCTION, so EW = integral of (1 - f) dlambda and nothing is
                # estimated. Using the observed-spectrum machinery here would introduce a
                # continuum uncertainty that does not exist in the model.
                _wA = w_nm * 10.0
                ews[tag] = float(np.trapz(1.0 - np.asarray(f, float), _wA) * 1000.0)
        except Exception as e:                      # loud per line, never a silent drop
            rows.append(dict(wavelength_air_A=lam, status=f"synth_failed: {str(e)[:80]}"))
            continue
        blend = ews.pop("blend")
        line = {k: v - blend for k, v in ews.items()}       # the LINE's own EW
        if min(line.values()) <= 0:
            rows.append(dict(wavelength_air_A=lam, status="non_positive_line_EW",
                             ew_blend=blend, **ews)); continue
        drew_dA = (np.log10(line["hi"] / lam)
                   - np.log10(line["lo"] / lam)) / (2 * a.delta)   # units cancel here
        rows.append(dict(wavelength_air_A=lam, status="ok", A0=A0,
                         ew_window_mid=ews["mid"], ew_blend=blend,
                         ew_line_mid=line["mid"],
                         blend_fraction=float(blend / ews["mid"]) if ews["mid"] else np.nan,
                         # 🔴 REW OF THE LINE, NOT OF THE WINDOW. `rew_max` is expressed
                         # in REW, so it must be the line's own reduced EW -- the window
                         # value includes the blend and reads systematically too strong,
                         # exactly the contamination the blend baseline exists to remove.
                         # EW is in mA and lambda in A, so the 1000 is not optional:
                         # without it the smoke test reported -1.07 for a line the pool
                         # measures near -4.9.
                         rew=float(np.log10((line["mid"] / 1000.0) / lam)),
                         rew_window=float(np.log10((ews["mid"] / 1000.0) / lam)),
                         d_rew_dA=float(drew_dA),
                         recovered_from_ceiling=bool(r.sat_only)))
        if i % 25 == 0:
            el = time.time() - t0
            print(f"  {i}/{len(pool)}  {el/i:.1f}s/line  eta {(len(pool)-i)*el/i/60:.0f} min",
                  flush=True)

    out = Path(a.out) if a.out else OUT / f"rya1043_perline_cog_{int(a.lo)}_{int(a.hi)}.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(rows)
    df.to_csv(out, index=False)
    ok = df[df.status == "ok"] if "status" in df.columns else df
    print(f"\nwrote {out}   ok={len(ok)}  failed={len(df)-len(ok)}")
    if len(ok):
        d = pd.to_numeric(ok.d_rew_dA, errors="coerce").dropna()
        print("\ndREW/dA distribution (1.0 = fully linear, 0 = fully saturated):")
        for q in (0, 5, 25, 50, 75, 95, 100):
            print(f"   p{q:<3d} {np.percentile(d, q):+.3f}")
        print("\nNo flattening threshold is applied. RYA-1041 Step 2 is Ryan's.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
