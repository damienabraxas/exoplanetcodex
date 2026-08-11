#!/usr/bin/env python3
"""RYA-760 — is the +0.176 dex a gf SCALE offset, or a population difference?

    python3 scripts/rya760_source_adjudication.py --lines-csv <1D-LTE per-line artifact>

THE CLAIM UNDER TEST
--------------------
RYA-760 measured that loosening the gf tier 0.04 -> 0.08 moves A(Fe) +0.176 dex away from
the optical anchor (Mann-Whitney p < 0.0001) and concluded: *"That is a gf-scale offset in
the looser sources, not new information."* The prescription follows from the diagnosis --
adjudicate FMW/GESB82c/MRW as sources, and a demonstrable scale offset becomes a finding.

The diagnosis has a confounder the ticket does not address. Comparing the GES catalogue's
Fe I entries in 6910-9199.9 A by source:

    source     n     median loggf   median EP
    BWL       24        -1.679        3.018
    FMW       67        -2.250        4.371

**FMW lines are systematically WEAKER (lower gf) and HIGHER EXCITATION than the laboratory
sources they are being compared against.** Both of those move a derived abundance on their
own, with no gf-scale error anywhere:

  * weaker lines sit on a different part of the curve of growth, and the strong-line end is
    where saturation and microturbulence bite (the RYA-523 territory);
  * higher EP lines weight a different depth of the atmosphere and carry the temperature
    sensitivity, so any Teff or structure error shows up as an EP trend in A(Fe).

So "the added 0.04-0.08 lines read +0.243 vs the current product's +0.067" is consistent
with TWO stories, and the ticket picks one:

    (a) the looser SOURCES carry a gf scale offset       <- the ticket's reading
    (b) the looser sources happen to hold weaker,
        higher-EP lines, and the offset is a curve-of-
        growth / excitation trend that any source with
        that line mix would show

They have very different consequences. Under (a) the fix is atomic data and a publishable
source offset. Under (b) the fix is ours -- an EP or saturation trend in our own analysis,
which would also be biasing the lines we DO trust, just less visibly.

WHAT THIS DOES
--------------
Separates them by controlling for the confounders instead of comparing raw medians:

  1. the raw per-source A(Fe), reproducing the ticket's numbers;
  2. the (EP, REW, loggf) distributions per source, to show the populations differ;
  3. a MATCHED comparison -- FMW lines against laboratory lines of similar EP and
     strength -- so the source term is measured at fixed line character;
  4. a multiple regression of A on EP, REW and a source indicator, so the source
     coefficient is what survives after the trends are removed.

If the source term survives the controls, the ticket is right and the offset is atomic.
If it collapses, the offset is ours and widening/narrowing the tier was never the question.
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, os.environ.get("ISPEC_DIR", "/srv/codex/engines/ispec_src"))

# Source accuracies as published, from the RYA-760 table. These are the ACCURACY of the
# source, not a quality judgement, and the tier is a cut on this number.
SOURCE_DEX = {
    "BWL": 0.02, "2014MNRAS.": 0.03, "GESHRL14": 0.03, "BK": 0.04, "BKK": 0.04,
    "FMW": 0.08, "GESB82c": 0.08, "MRW": 0.08, "K07": 0.20, "K13": 0.20, "K75": 0.20,
}
LAB_TIER = 0.04


def base_source(ref: str) -> str:
    """'BK+BWL' -> 'BK'. Composite tags take the accuracy of their first component."""
    r = str(ref).strip()
    return r.split("+")[0] if "+" in r else r


def accuracy(ref: str) -> float:
    return SOURCE_DEX.get(base_source(ref), np.nan)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--lines-csv", required=True,
                    help="the 1D-LTE per-line artifact (wavelength_air_A, abundance, ew_mA)")
    ap.add_argument("--lo", type=float, default=6910.0)
    ap.add_argument("--hi", type=float, default=9199.9)
    a = ap.parse_args()

    from pipeline.abundances_derive import _load_synth_resources
    ll, _, _ = _load_synth_resources()
    w = np.asarray(ll["wave_A"], float)
    el = np.asarray([str(x).strip() for x in ll["element"]])
    ref = np.asarray([str(x).strip() for x in ll["reference_code"]])
    gf = np.asarray(ll["loggf"], float)
    ep = np.asarray(ll["lower_state_eV"], float)
    fe1 = np.array([e.upper().startswith("FE 1") for e in el])

    d = pd.read_csv(a.lines_csv)
    d = d[d.get("in_aggregate", True).astype(bool)] if "in_aggregate" in d else d
    d = d[(d.wavelength_air_A >= a.lo) & (d.wavelength_air_A <= a.hi)]
    d = d[d.abundance.notna()].copy()

    rows = []
    for _, r in d.iterrows():
        c = float(r.wavelength_air_A)
        cand = np.where(fe1 & (np.abs(w - c) < 0.05))[0]
        if not len(cand):
            continue
        i = cand[np.argmax(gf[cand])]
        ew = float(r.ew_mA) if "ew_mA" in d.columns and np.isfinite(r.ew_mA) else np.nan
        rows.append(dict(wave=c, A=float(r.abundance), ew_mA=ew,
                         rew=(np.log10(ew / c / 1000.0) if np.isfinite(ew) and ew > 0
                              else np.nan),
                         ref=ref[i], base=base_source(ref[i]),
                         acc=accuracy(ref[i]), loggf=gf[i], ep=ep[i]))
    t = pd.DataFrame(rows)
    t = t[t.acc.notna()]
    print(f"\n{len(t)} measured Fe I lines in {a.lo:.0f}-{a.hi:.0f} A with a known source")

    # ── 1. reproduce the ticket's split ───────────────────────────────────────
    lab = t[t.acc <= LAB_TIER]
    loose = t[(t.acc > LAB_TIER) & (t.acc <= 0.08)]
    print("\n" + "=" * 78)
    print("[1] the ticket's split, reproduced")
    print("=" * 78)
    for name, s in (("tier <= 0.04 (product)", lab), ("added 0.04-0.08", loose)):
        if len(s):
            print(f"  {name:<26} n={len(s):3d}  median A={s.A.median():.3f}  "
                  f"sd={s.A.std(ddof=1) if len(s) > 1 else float('nan'):.3f}")
    if len(lab) and len(loose):
        print(f"  difference: {loose.A.median() - lab.A.median():+.3f} dex")

    # ── 2. do the populations differ? ─────────────────────────────────────────
    print("\n" + "=" * 78)
    print("[2] ARE THESE THE SAME KIND OF LINE?  (the confounder)")
    print("=" * 78)
    print(f"{'set':<26}{'n':>4}{'med EP':>9}{'med loggf':>11}{'med REW':>9}{'med EW':>9}")
    for name, s in (("tier <= 0.04", lab), ("added 0.04-0.08", loose)):
        if len(s):
            print(f"{name:<26}{len(s):>4}{s.ep.median():>9.3f}{s.loggf.median():>11.3f}"
                  f"{s.rew.median():>9.3f}{s.ew_mA.median():>9.1f}")
    print("\n  per source:")
    print(f"{'source':<12}{'acc':>6}{'n':>4}{'med A':>9}{'med EP':>9}{'med loggf':>11}{'med REW':>9}")
    for b, s in t.groupby("base"):
        print(f"{b:<12}{s.acc.iloc[0]:>6.2f}{len(s):>4}{s.A.median():>9.3f}"
              f"{s.ep.median():>9.3f}{s.loggf.median():>11.3f}{s.rew.median():>9.3f}")

    # ── 3. matched comparison ─────────────────────────────────────────────────
    print("\n" + "=" * 78)
    print("[3] MATCHED — FMW vs laboratory lines of SIMILAR EP and strength")
    print("=" * 78)
    fmw = t[t.base == "FMW"]
    if len(fmw) and len(lab):
        pairs = []
        for _, r in fmw.iterrows():
            c = lab[(np.abs(lab.ep - r.ep) < 0.30) & (np.abs(lab.rew - r.rew) < 0.20)]
            if len(c):
                pairs.append(dict(wave=r.wave, A_fmw=r.A, A_lab=c.A.median(),
                                  n_match=len(c), ep=r.ep, rew=r.rew))
        pm = pd.DataFrame(pairs)
        if len(pm):
            dlt = pm.A_fmw - pm.A_lab
            print(f"  {len(pm)} FMW lines found a laboratory match "
                  f"(|dEP|<0.30 eV, |dREW|<0.20)")
            print(f"  median A(FMW) - A(matched lab) = {dlt.median():+.3f} dex")
            print(f"  vs the RAW difference          = "
                  f"{fmw.A.median() - lab.A.median():+.3f} dex")
            print(f"\n  => matching on line character removes "
                  f"{100*(1 - abs(dlt.median())/max(abs(fmw.A.median()-lab.A.median()),1e-9)):.0f}% "
                  f"of the apparent offset")
        else:
            print("  no FMW line found a laboratory match within the window — the two")
            print("  sets do not overlap in (EP, REW) AT ALL, which is itself the answer.")

    # ── 4. regression ─────────────────────────────────────────────────────────
    print("\n" + "=" * 78)
    print("[4] REGRESSION — A ~ EP + REW + is_loose")
    print("=" * 78)
    s = t[t.rew.notna()].copy()
    s["is_loose"] = (s.acc > LAB_TIER).astype(float)
    X = np.column_stack([np.ones(len(s)), s.ep.values, s.rew.values, s.is_loose.values])
    y = s.A.values
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    dof = max(len(s) - X.shape[1], 1)
    se = np.sqrt(np.sum(resid ** 2) / dof * np.diag(np.linalg.pinv(X.T @ X)))
    names = ["intercept", "EP (eV)", "REW", "is_loose"]
    print(f"{'term':<12}{'coef':>10}{'se':>9}{'t':>8}")
    for n, b, e in zip(names, beta, se):
        print(f"{n:<12}{b:>10.4f}{e:>9.4f}{(b/e if e else np.nan):>8.2f}")
    print(f"\n  raw loose-vs-lab difference : "
          f"{loose.A.median() - lab.A.median():+.3f} dex")
    print(f"  is_loose coefficient        : {beta[3]:+.3f} dex  "
          f"(t = {beta[3]/se[3] if se[3] else float('nan'):.2f})")
    print("\n  If the coefficient is much smaller than the raw difference, the offset was")
    print("  the LINE MIX (EP/strength), not the source's gf scale.")


if __name__ == "__main__":
    main()
