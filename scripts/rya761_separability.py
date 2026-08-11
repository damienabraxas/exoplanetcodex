#!/usr/bin/env python3
"""RYA-761 — MEASURE the two-component separability floor. Do not assume one.

    python3 scripts/rya761_separability.py [--out data/audit/rya761]

WHY THIS EXISTS
---------------
RYA-761 sets three preconditions for recovering a quarantined line, and the second is
*"the two centres are separated enough for the fit to be identifiable — **needs a
measured floor, not an assumed one**."* Picking 0.15 A because it looks reasonable is
exactly the move the ticket forbids, and it is how a recovery campaign quietly
manufactures 290 bad measurements. So this measures it by injection and recovery.

WHAT IS INJECTED
----------------
A two-component Voigt blend on a unit continuum, over a grid of the four things that
actually decide identifiability:

  * **separation**      0.02 - 0.60 A, the axis the floor lives on
  * **depth ratio**     interloper/target 0.1 - 3.0; a neighbour DEEPER than the target
                        is the BLEND-DOMINATED case (105 of the 290) and is the hard one
  * **noise**           0.0005 - 0.005 in continuum units, bracketing the Kitt Peak FTS
  * **catalogue error** 0.00 - 0.05 A applied to the interloper centre the fitter is
                        GIVEN, while the injected line stays put

That last axis is the one an idealised study would omit and the one that decides
whether this works on real data. The fit holds the interloper centre FIXED (RYA-761
forbids freeing it), so a catalogue that is wrong by 0.03 A — the RYA-704 disagreement
between two catalogues quoting the same transition — is pinning the second component in
the wrong place. If recovery collapses there, the floor is set by catalogue quality,
not by the optimiser.

WHAT IS MEASURED
----------------
The target's equivalent width, integrated from its OWN fitted component, against the
injected truth. Reported as the fractional EW error, alongside the single-component fit
on the identical spectrum so the comparison is like-for-like. The floor is the smallest
separation at which the median error stays inside tolerance AND beats the
single-component fit.

It writes a CSV of every trial: this is evidence, not a claim.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline.lines_fit import _blend_model, _fit_blend, _fit_profile, _voigt_abs  # noqa: E402

# Solar metal lines at Kitt Peak in the 6910-9199 A band: the observed width is star
# (+) instrument in quadrature and the star dominates at R ~ 500 000 (RYA-713).
SIGMA_TRUE = 0.045
GAMMA_TRUE = 0.008
SIGMA_INIT, SIGMA_MIN = 0.045, 0.005

SEPARATIONS = np.array([0.02, 0.04, 0.06, 0.08, 0.10, 0.13, 0.16, 0.20,
                        0.25, 0.30, 0.40, 0.50, 0.60])
DEPTH_RATIOS = [0.1, 0.3, 1.0, 3.0]
NOISES = [0.0005, 0.002, 0.005]
CATALOGUE_ERRORS = [0.00, 0.01, 0.03, 0.05]
N_REALISATIONS = 12
TARGET_DEPTH = 0.30
CENTRE = 8000.0
TOL = 0.05                      # 5 % EW error — the pass bar


def component_ew(depth: float, sigma: float, gamma: float,
                 half: float = 3.0, n: int = 4001) -> float:
    """EW (in A) of ONE fitted component, integrated from its own parameters."""
    x = np.linspace(-half, half, n)
    return float(np.trapezoid(1.0 - _voigt_abs(x, 0.0, depth, sigma, gamma), x))


def run_grid(rng: np.random.Generator) -> pd.DataFrame:
    rows = []
    true_ew = component_ew(TARGET_DEPTH, SIGMA_TRUE, GAMMA_TRUE)

    for sep in SEPARATIONS:
        for ratio in DEPTH_RATIOS:
            d_int = TARGET_DEPTH * ratio
            for noise in NOISES:
                for cat_err in CATALOGUE_ERRORS:
                    for k in range(N_REALISATIONS):
                        half = max(0.60, sep + 0.45)
                        wav = np.linspace(CENTRE - half, CENTRE + half, 400)
                        flux = np.ones_like(wav)
                        for c, d in ((CENTRE, TARGET_DEPTH),
                                     (CENTRE + sep, d_int)):
                            flux -= 1.0 - _voigt_abs(wav, c, d, SIGMA_TRUE, GAMMA_TRUE)
                        flux = flux + rng.normal(0.0, noise, wav.size)

                        kw = dict(sigma_init=SIGMA_INIT, sigma_min=SIGMA_MIN)

                        # what the fitter is TOLD the neighbour's centre is
                        given = CENTRE + sep + cat_err
                        p2, _, kind2, chi2 = _fit_blend(
                            wav, flux, [CENTRE, given], **kw)
                        p1, _, kind1, chi1 = _fit_profile(wav, flux, CENTRE, **kw)

                        def ew_of(p, kind):
                            if p is None:
                                return np.nan
                            g = p[3] if kind == 'voigt' else 0.0
                            return component_ew(p[1], p[2], g)

                        ew2, ew1 = ew_of(p2, kind2), ew_of(p1, kind1)
                        rows.append(dict(
                            separation_A=sep, depth_ratio=ratio, noise=noise,
                            catalogue_error_A=cat_err, realisation=k,
                            true_ew_A=true_ew,
                            ew_two_A=ew2, ew_one_A=ew1,
                            err_two=abs(ew2 - true_ew) / true_ew,
                            err_one=abs(ew1 - true_ew) / true_ew,
                            interloper_depth=(p2[-1] if p2 is not None else np.nan),
                            chi2_two=chi2, chi2_one=chi1))
    return pd.DataFrame(rows)


def report(df: pd.DataFrame) -> None:
    print("\n" + "=" * 78)
    print("MEASURED SEPARABILITY — median fractional EW error on the TARGET")
    print("=" * 78)

    clean = df[df.catalogue_error_A == 0.0]
    print("\n[A] perfect catalogue wavelength (cat_err = 0.00 A)")
    print(f"{'sep (A)':>8} {'2-comp':>9} {'1-comp':>9} {'better?':>8}   by depth ratio")
    for sep, g in clean.groupby('separation_A'):
        e2, e1 = g.err_two.median(), g.err_one.median()
        per = "  ".join(f"{r}:{gg.err_two.median():.3f}"
                        for r, gg in g.groupby('depth_ratio'))
        print(f"{sep:8.2f} {e2:9.3f} {e1:9.3f} {'YES' if e2 < e1 else 'no':>8}   {per}")

    ok = [s for s, g in clean.groupby('separation_A')
          if g.err_two.median() < TOL and g.err_two.median() < g.err_one.median()]
    floor = min(ok) if ok else None
    print(f"\n  => FLOOR (median err < {TOL:.0%} and beats 1-component): "
          f"{floor:.2f} A" if floor else "\n  => no separation met the bar")
    if floor:
        print(f"     = {floor / SIGMA_TRUE:.2f} x sigma  (sigma = {SIGMA_TRUE} A)")

    print("\n[B] the hard case — interloper DEEPER than target (ratio 3.0), cat_err 0")
    hard = clean[clean.depth_ratio == 3.0]
    for sep, g in hard.groupby('separation_A'):
        print(f"{sep:8.2f} 2-comp {g.err_two.median():6.3f}   "
              f"1-comp {g.err_one.median():6.3f}")

    print("\n[C] CATALOGUE WAVELENGTH ERROR — the axis that decides real-data use")
    print(f"{'cat_err':>9}  " + "  ".join(f"sep{s:.2f}" for s in
                                          [0.10, 0.20, 0.30, 0.40, 0.60]))
    for ce, g in df.groupby('catalogue_error_A'):
        cells = []
        for s in [0.10, 0.20, 0.30, 0.40, 0.60]:
            sub = g[g.separation_A == s]
            cells.append(f"{sub.err_two.median():7.3f}" if len(sub) else "      -")
        print(f"{ce:9.2f}  " + "  ".join(cells))

    print("\n[D] does it INVENT a neighbour? median fitted interloper depth when the")
    print("    true interloper is shallow (ratio 0.1, true depth 0.030)")
    inv = clean[clean.depth_ratio == 0.1]
    for sep, g in inv.groupby('separation_A'):
        print(f"{sep:8.2f}  fitted {g.interloper_depth.median():.4f}  (true 0.0300)")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default=str(ROOT / "data" / "audit" / "rya761"))
    a = ap.parse_args()

    rng = np.random.default_rng(20260811)
    df = run_grid(rng)
    report(df)

    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    f = out / "separability_injection_recovery.csv"
    df.to_csv(f, index=False)
    print(f"\n  wrote {f.relative_to(ROOT)}  ({len(df)} trials)")


if __name__ == "__main__":
    main()
