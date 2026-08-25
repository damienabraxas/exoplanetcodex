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
        frames = []
        for f in sorted(glob.glob(str(BP / f"FeI_{blo}_{bhi}_*_PROFILEFIT_*1D-LTE_lines.csv"))):
            frames.append(pd.read_csv(f))
        if not frames:
            continue
        d = pd.concat(frames, ignore_index=True)
        if "in_aggregate" in d.columns:
            d = d[d["in_aggregate"] == True]
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
        verdict = ("REAL KNEE — the response collapses" if r < 0.4 else
                   "WEAK / NO KNEE — both segments have similar slope, so a breakpoint "
                   "here is two lines fitted to one")
        print(f"    slope ratio above/below = {r:.3f}  ->  {verdict}")
        frac = knee["n_above"] / (knee["n_above"] + knee["n_below"])
        print(f"    a cut at the knee would remove {frac*100:.1f}% of this band's pool")

    (OUT / "rya1041_cog_knee.json").write_text(json.dumps(report, indent=2) + "\n")
    print(f"\nwrote {OUT / 'rya1041_cog_knee.json'}")
    print("Sets nothing. rew_max is Ryan's call (RYA-1041 Step 2).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
