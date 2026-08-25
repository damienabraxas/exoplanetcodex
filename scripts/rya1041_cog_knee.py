#!/usr/bin/env python3
"""The curve-of-growth saturation knee, per band — RYA-1041 Step 2 input for rew_max.

`rew_max` is asked for as PHYSICS: the point where the curve of growth flattens and a
line stops responding to abundance, not a gap in a histogram. Below the knee EW is
proportional to the number of absorbers and REW rises with unit slope against the line
strength; past it the core saturates, dREW/dX collapses, and an abundance recovered from
that line is barely constrained -- which is exactly why such a line cannot be gf-graded
from its own behaviour.

CONSTRUCTION. The classical empirical COG, on OUR measured EWs:

    X = log10(gf) + log10(lambda_A) - theta * chi_eV      (line strength)
    Y = REW = log10(EW / lambda)                          (measured response)
    theta = 5040 / Teff

Every Fe I line shares one solar Fe abundance, so the spread in X comes from gf and the
excitation term and the ensemble traces the curve. The knee is read from a fitted
piecewise-linear (two-segment) model: the breakpoint that minimises total residual, with
the slope on each side reported so "is it actually a knee" has numbers rather than an
eyeball.

🔴 THE gf JOIN IS THE KNOWN-DEFECTIVE STEP (RYA-1033: the measured-EW -> canonical_gf
join is keyed on ROUNDED wavelength). Joined here on a tolerance with an AMBIGUITY
REFUSAL: two canonical candidates inside the tolerance means the tolerance cannot
separate them, and picking the nearest is how a line acquires another line's gf.
Unmatched and ambiguous lines are COUNTED and REPORTED, never silently dropped.

Sets nothing. RYA-1041 Step 2 is Ryan's.
"""
from __future__ import annotations

import glob
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
BP = ROOT / "data" / "results" / "band_products"
OUT = ROOT / "data" / "audit" / "rya1041"

TEFF = 5772.0                     # solar, the node these products were derived at
THETA = 5040.0 / TEFF
WTOL = 0.02                       # A — same tolerance family as gf_grades.WAVE_TOL_A

BANDS = {"near-UV": (3000, 3780), "VIS": (4200, 6910),
         "red-optical": (6910, 9199), "IR": (9199, 12976)}


def piecewise_knee(x: np.ndarray, y: np.ndarray) -> dict:
    """Two connected straight segments; return the breakpoint that fits best.

    Reports BOTH slopes. A real knee has a high-slope segment and a low-slope segment;
    if the two slopes are similar the 'knee' is an artifact of fitting two lines to one
    line, and the caller must be able to see that.
    """
    o = np.argsort(x)
    x, y = x[o], y[o]
    best = None
    lo, hi = np.quantile(x, 0.10), np.quantile(x, 0.90)
    for xb in np.linspace(lo, hi, 120):
        left, right = x <= xb, x > xb
        if left.sum() < 8 or right.sum() < 8:
            continue
        # continuous at xb: fit y = a + s1*(x-xb) for x<=xb, a + s2*(x-xb) for x>xb
        A = np.zeros((len(x), 3))
        A[:, 0] = 1.0
        A[:, 1] = np.where(left, x - xb, 0.0)
        A[:, 2] = np.where(right, x - xb, 0.0)
        coef, res, *_ = np.linalg.lstsq(A, y, rcond=None)
        rss = float(np.sum((A @ coef - y) ** 2))
        if best is None or rss < best["rss"]:
            best = {"rss": rss, "x_break": float(xb), "intercept": float(coef[0]),
                    "slope_below": float(coef[1]), "slope_above": float(coef[2]),
                    "n_below": int(left.sum()), "n_above": int(right.sum())}
    if best is None:
        return {"verdict": "too few points"}
    # the REW at the break is what rew_max would be expressed in
    best["rew_at_break"] = best["intercept"]
    best["slope_ratio"] = (best["slope_above"] / best["slope_below"]
                           if best["slope_below"] else float("nan"))
    return best


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    cg = pd.read_csv(ROOT / "data" / "linelists" / "canonical_gf.csv", low_memory=False)
    cg = cg[cg.species == "Fe I"]
    cw = pd.to_numeric(cg.wavelength_air_A, errors="coerce").values
    cgf = pd.to_numeric(cg.log_gf, errors="coerce").values

    report = {"ticket": "RYA-1041", "quantity": "COG saturation knee (rew_max input)",
              "teff": TEFF, "theta": THETA, "sets_nothing": True, "bands": {}}

    for band, (blo, bhi) in BANDS.items():
        # 🔴 READ THE STAGE-1 MEASURED POOL, NOT THE DERIVED PRODUCT.
        # `data/results/band_products/*_lines.csv` has already had REW_SATURATION_CEILING
        # (-4.9) applied upstream: its max REW is -4.9006 with five lines piled into the
        # last 0.01 dex. Measuring a saturation knee there is circular -- the saturated
        # lines that DEFINE the knee were removed before the file was written, and the
        # two-segment fit duly came out backwards (steep side above, flat side below).
        # `data/measured/band_ew/*_ew.csv` is the pool BEFORE that cut: REW runs to -3.71.
        MEASURED = ROOT / "data" / "measured" / "band_ew"
        frames = [pd.read_csv(f) for f in
                  sorted(glob.glob(str(MEASURED / f"FeI_{blo}_{bhi}_*_PROFILEFIT_ew.csv")))]
        if not frames:
            continue
        d = pd.concat(frames, ignore_index=True)
        # 🔴 CLEAN **AND** UNTRUNCATED, and it takes both to see a curve of growth.
        # The full candidate set does NOT trace one: corr(X, REW) = -0.015 over 5335 rows,
        # with binned medians flat across 8 dex of X, because it is dominated by blends
        # and misidentifications where the measured EW tracks whatever sits at that
        # wavelength rather than the Fe line's own strength. The derived product IS clean
        # but is cut at -4.9, which removes the knee.
        # The stage-1 file records REW exclusions explicitly ("Measured and reported;
        # excluded from the aggregate ONLY"), so lines cut SOLELY by the ceiling passed
        # every other quality gate and are recoverable. That union is the right pool:
        # corr(X, REW) = +0.362 and the binned medians rise monotonically.
        _er = d.excluded_reason.astype(str)
        _sat_only = _er.str.startswith("REW ") & _er.str.contains("saturation ceiling")
        d = d[(d["in_aggregate"] == True) | _sat_only]
        d = d.assign(w=pd.to_numeric(d.wavelength_air_A, errors="coerce"),
                     rew=pd.to_numeric(d.rew, errors="coerce"),
                     ep=pd.to_numeric(d.ep_eV, errors="coerce")).dropna(subset=["w", "rew", "ep"])
        if not len(d):
            continue

        gfs, n_amb, n_miss = [], 0, 0
        for w in d.w.values:
            m = np.flatnonzero(np.abs(cw - w) <= WTOL)
            if m.size == 0:
                gfs.append(np.nan); n_miss += 1
            elif m.size > 1:
                gfs.append(np.nan); n_amb += 1        # REFUSED, not resolved by proximity
            else:
                gfs.append(cgf[m[0]])
        d = d.assign(log_gf=gfs).dropna(subset=["log_gf"])

        X = d.log_gf.values + np.log10(d.w.values) - THETA * d.ep.values
        Y = d.rew.values
        knee = piecewise_knee(X, Y)
        # The spec asks for dREW/dX falling below threshold, so measure the RESPONSE
        # PROFILE directly rather than trusting a global two-segment fit -- which put its
        # break at the 16th percentile of X, in the sparse low-X tail, not at saturation.
        prof, prev = [], None
        _o = np.argsort(X); _X, _Y = X[_o], Y[_o]
        _e = np.quantile(_X, np.linspace(0, 1, 15))
        for i in range(14):
            m = (_X >= _e[i]) & (_X <= _e[i + 1])
            if m.sum() < 20:
                continue
            xm, ym = float(np.median(_X[m])), float(np.median(_Y[m]))
            prof.append({"x": xm, "rew": ym, "n": int(m.sum()),
                         "d_rew_dx": None if prev is None else (ym - prev[1]) / (xm - prev[0])})
            prev = (xm, ym)
        sl = [r["d_rew_dx"] for r in prof if r["d_rew_dx"] is not None]
        knee["slope_profile"] = prof
        knee["profile_monotonic"] = bool(sl and all(b >= a for a, b in zip(sl, sl[1:])))
        knee["n_sign_changes"] = int(sum(1 for a, b in zip(sl, sl[1:]) if a * b < 0)) if sl else 0
        entry = {"n_used": int(len(d)), "n_gf_unmatched": n_miss,
                 "n_gf_AMBIGUOUS_refused": n_amb, "knee": knee,
                 "rew_span": [float(np.min(Y)), float(np.max(Y))]}
        report["bands"][band] = entry

        print(f"\n### {band}   n={len(d)}   (gf unmatched {n_miss}, ambiguous refused {n_amb})")
        if "verdict" in knee:
            print(f"    {knee['verdict']}"); continue
        print(f"    breakpoint X = {knee['x_break']:.3f}   REW at break = "
              f"{knee['rew_at_break']:.3f}")
        print(f"    slope below {knee['slope_below']:.3f}  (n={knee['n_below']})   "
              f"slope above {knee['slope_above']:.3f}  (n={knee['n_above']})")
        r = knee["slope_ratio"]
        # 🔴 A NEGATIVE RATIO IS NOT A KNEE. My first verdict rule was `r < 0.4`, which
        # accepted r = -0.511 as a "REAL KNEE" -- but opposite SIGNS mean a V or a peak,
        # not a response that flattens. A curve of growth needs BOTH slopes positive
        # (REW rises with line strength) and the upper one much smaller.
        sb, sa = knee["slope_below"], knee["slope_above"]
        if sb <= 0 or sa < 0:
            verdict = ("NOT A COG SHAPE — a curve of growth rises on both sides "
                       f"(slope_below {sb:+.3f}, slope_above {sa:+.3f}); this is not a "
                       "saturation knee whatever the fit residual says")
        elif r < 0.4:
            verdict = "REAL KNEE — the response collapses past the break"
        else:
            verdict = ("WEAK / NO KNEE — both segments have similar slope, so a "
                       "breakpoint here is two lines fitted to one")
        print(f"    slope ratio above/below = {r:.3f}  ->  {verdict}")
        frac = knee["n_above"] / (knee["n_above"] + knee["n_below"])
        print(f"    a cut at the knee would remove {frac*100:.1f}% of this band's pool")
        print(f"    local dREW/dX sign changes across the profile: "
              f"{knee['n_sign_changes']}  ->  "
              + ("a clean linear-then-flat response" if knee["n_sign_changes"] <= 1
                 else "OSCILLATING: no linear regime resolves, so no knee can be located"))

    (OUT / "rya1041_cog_knee.json").write_text(json.dumps(report, indent=2) + "\n")
    print(f"\nwrote {OUT / 'rya1041_cog_knee.json'}")
    print("Sets nothing. rew_max is Ryan's call (RYA-1041 Step 2).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
