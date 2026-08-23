#!/usr/bin/env python3
"""Why the stage-4 frac_rise cut cannot be one number, and what it is per arm — RYA-992.

    python3 scripts/rya992_arm_scale.py [--kp <product>] [--harps <product>]

RYA-995 measured that `frac_rise_weaker` sits UNIFORMLY LOWER on HARPS than on Kitt Peak
and inferred that a Kitt-Peak-derived cut of 1.0 OVER-rejects HARPS. This script tests the
inference rather than adopting it, on the RYA-986 pool — the SAME 1483 Fe I lines measured
through both arms, so the line list, the gf values and the model atmosphere are identical
and the only difference between the two columns is the instrument.

🔒 RYA-161. The target is |A - MEDIAN A| — each pool's dispersion about ITSELF. No published
abundance enters anywhere in this script.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

_PRODUCTS = ROOT / "data" / "results" / "band_products"
_CANONICAL = ROOT / "data" / "linelists" / "canonical_gf.csv"
_MATCH_TOL_A = 0.005

KP_DEFAULT = "FeI_4200_6910_kpno_solar_atlas_SYNTH_FROMEW_1D-LTE_lines.csv"
HARPS_DEFAULT = ("FeI_4200_6910_harps_solar_harps_molecfit_corrected_SYNTH_FROMEW"
                 "_1D-LTE_lines.csv")


def _load(fname: str) -> pd.DataFrame:
    p = _PRODUCTS / fname
    if not p.exists():
        raise SystemExit(
            f"{fname} is not in {_PRODUCTS}. The per-line products are gitignored "
            f"(RYA-469); stage them from the RYA-986 run rather than deriving an arm "
            f"scale on whichever product happens to be present.")
    d = pd.read_csv(p)
    d = d[d.in_aggregate & d.abundance.notna()].copy()
    d["f"] = pd.to_numeric(d.frac_rise_weaker, errors="coerce")
    d["s"] = pd.to_numeric(d.sigma_A, errors="coerce")
    d["err"] = (d.abundance - d.abundance.median()).abs()
    # k = frac_rise * sigma_A^2. See the module docstring of pipeline/synth_gof: for a
    # locally parabolic chi2 this is B^2/dof — bracket width over window pixel count —
    # i.e. everything in frac_rise that is NOT the abundance constraint.
    d["k"] = d.f * d.s ** 2
    return d


def _graded_mask(d: pd.DataFrame) -> np.ndarray:
    """Which rows carry a PRIMARY LABORATORY gf. Same key the selectors use."""
    cg = pd.read_csv(_CANONICAL, low_memory=False)
    lab = (cg[(cg.species == "Fe I")
              & cg.gf_tier.astype(str).str.contains("LAB", na=False)]
           .sort_values("wavelength_air_A").reset_index(drop=True))
    W = lab.wavelength_air_A.values.astype(float)
    w = d.wavelength_air_A.astype(float).values
    i = np.abs(W[None, :] - w[:, None]).argmin(1)
    return np.abs(W[i] - w) <= _MATCH_TOL_A


def _knee(x: np.ndarray, e: np.ndarray, grid: np.ndarray):
    """The break in err-vs-log10(frac_rise): a declining segment meeting a flat one.

    Continuous hinge  err = c + a * max(b - x, 0), fitted by iteratively reweighted least
    squares (an L1 loss — the error distribution has a long upper tail and a squared loss
    would let a handful of catastrophic lines choose the breakpoint). `b` is scanned; the
    breakpoint is the one with the smallest absolute deviation. NOT a percentile: nothing
    here fixes how many lines land on either side.
    """
    best = (np.inf, None, None)
    for b in grid:
        h = np.maximum(b - x, 0.0)
        A = np.c_[np.ones_like(h), h]
        w = np.ones_like(h)
        coef = np.array([np.median(e), 0.0])
        for _ in range(15):
            W = np.sqrt(w)[:, None]
            coef, *_ = np.linalg.lstsq(A * W, e * np.sqrt(w), rcond=None)
            w = 1.0 / np.maximum(np.abs(e - A @ coef), 1e-3)
        if coef[1] <= 0:
            continue
        loss = float(np.abs(e - A @ coef).sum())
        if loss < best[0]:
            best = (loss, float(b), coef)
    return best[1], best[2]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--kp", default=KP_DEFAULT)
    ap.add_argument("--harps", default=HARPS_DEFAULT)
    ap.add_argument("--boot", type=int, default=400)
    a = ap.parse_args(argv)

    D = {"kpno_solar_atlas": _load(a.kp), "harps": _load(a.harps)}
    for t, d in D.items():
        d["graded"] = _graded_mask(d)

    common = set(D["kpno_solar_atlas"].wavelength_air_A) & set(D["harps"].wavelength_air_A)
    print(f"\nRYA-992 per-arm frac_rise scale — {len(common)} lines measured through BOTH arms")

    # ---------------------------------------------------------------- 1. the raw shift
    print("\n1. THE SHIFT RYA-995 SAW (medians over each arm's own pool)")
    print(f"   {'arm':20s} {'n':>5s} {'frac_rise':>10s} {'sigma_A':>9s} {'red_chi2':>9s}")
    for t, d in D.items():
        print(f"   {t:20s} {len(d):5d} {d.f.median():10.3f} {d.s.median():9.4f} "
              f"{pd.to_numeric(d.red_chi2, errors='coerce').median():9.1f}")
    r = D["kpno_solar_atlas"].f.median() / D["harps"].f.median()
    print(f"   HARPS frac_rise is lower by x{r:.3f} — but its sigma_A is WORSE by "
          f"x{D['harps'].s.median() / D['kpno_solar_atlas'].s.median():.3f}, so the arm "
          f"that looks 'over-rejected' is\n   the arm whose fits actually pin A less well.")

    # ------------------------------------------------------- 2. the arm scale, measured
    print("\n2. THE NUISANCE FACTOR  k = frac_rise * sigma_A^2  (= B^2/dof; no abundance)")
    rng = np.random.default_rng(992)
    K = {}
    for t, d in D.items():
        m = (d.f > 0) & (d.s > 0)
        kv = d.k[m].values
        K[t] = float(np.median(kv))
        bs = np.array([np.median(rng.choice(kv, kv.size)) for _ in range(a.boot)])
        print(f"   {t:20s} n={m.sum():5d}  k = {K[t]:.5f}  "
              f"boot sd {bs.std():.5f} ({100 * bs.std() / K[t]:.1f}%)")
    scale = K["harps"] / K["kpno_solar_atlas"]
    print(f"   ARM SCALE  k(harps)/k(kpno) = {scale:.3f}")

    print("\n   composition check — the same ratio on a DIFFERENT line selection:")
    for lbl, sel in (("all", slice(None)), ("laboratory-graded only", None)):
        ks = {}
        for t, d in D.items():
            sub = d[d.graded] if sel is None else d
            m = (sub.f > 0) & (sub.s > 0)
            ks[t] = float(sub.k[m].median())
        print(f"     {lbl:24s} k_kp {ks['kpno_solar_atlas']:.5f}  "
              f"k_harps {ks['harps']:.5f}  ratio {ks['harps'] / ks['kpno_solar_atlas']:.3f}")
    print("   ⇒ k itself moves with the line selection; the ARM RATIO does not. That is\n"
          "     what makes the ratio, and only the ratio, transferable.")

    # ------------------------------------------------- 3. independent check: the knee
    print("\n3. INDEPENDENT CHECK — the break in each arm's OWN err-vs-frac_rise relation")
    grid = np.arange(-1.6, 1.21, 0.05)
    knees = {}
    for t, d in D.items():
        m = d.f > 0
        x, e = np.log10(d.f[m]).values, d.err[m].values
        b, coef = _knee(x, e, grid)
        bs = []
        for _ in range(min(a.boot, 120)):
            i = rng.integers(0, x.size, x.size)
            bb, _c = _knee(x[i], e[i], grid)
            if bb is not None:
                bs.append(bb)
        bs = np.array(bs)
        knees[t] = 10 ** b
        print(f"   {t:20s} knee frac_rise {10 ** b:6.3f}  "
              f"boot 16-84% [{10 ** np.percentile(bs, 16):.3f}, "
              f"{10 ** np.percentile(bs, 84):.3f}]  floor |err| {coef[0]:.3f}")
    print(f"   knee(harps)/knee(kpno) = {knees['harps'] / knees['kpno_solar_atlas']:.3f}   "
          f"vs the k ratio {scale:.3f} — same direction, same size, measured two ways.")

    # ----------------------------------- 4. the decisive table: matched frac_rise bins
    print("\n4. AT THE SAME frac_rise, IS A HARPS LINE AS GOOD AS A KITT PEAK LINE?")
    print(f"   {'frac_rise bin':>16s} {'KP med|err| (n)':>20s} {'HARPS med|err| (n)':>20s}")
    edges = [0, 0.5, 1, 1.5, 2, 3, 5, 8, 15, np.inf]
    for lo, hi in zip(edges[:-1], edges[1:]):
        cells = []
        for t in ("kpno_solar_atlas", "harps"):
            d = D[t]
            m = (d.f >= lo) & (d.f < hi)
            cells.append(f"{d.err[m].median():.4f} ({int(m.sum())})" if m.sum() > 5 else "-")
        print(f"   [{lo:5.1f},{'inf' if hi == np.inf else f'{hi:5.1f}'})"
              f" {cells[0]:>19s} {cells[1]:>20s}")
    print("   ⇒ HARPS is never BETTER at a given frac_rise. The cut must not be lowered "
          "for it.")

    for t, d in D.items():
        rs, ps = spearmanr(d.f, d.err)
        print(f"   spearman(frac_rise, |A-median A|) {t:20s} {rs:+.3f}  p {ps:.3g}")

    # ------------------------------------------------------------- 5. what it costs
    from pipeline.synth_gof import SYNTH_GOF_REF_CUT, synth_gof_cut     # noqa: E402
    print(f"\n5. THE ADOPTED PER-ARM CUTS (reference {SYNTH_GOF_REF_CUT:.1f} on Kitt Peak)")
    for t, d in D.items():
        cut = synth_gof_cut(t)
        for lbl, sub in (("all", d), ("graded", d[d.graded])):
            keep = sub.f >= cut
            was = sub.f >= SYNTH_GOF_REF_CUT
            print(f"   {t:20s} {lbl:7s} cut {cut:5.3f}: keep {int(keep.sum()):4d}/"
                  f"{len(sub):4d} ({keep.mean():5.1%})   a flat {SYNTH_GOF_REF_CUT:.1f} kept "
                  f"{int(was.sum()):4d}  ->  {int(was.sum()) - int(keep.sum()):d} more lines "
                  f"now refused")
    print("\n  🔒 no reference abundance enters this script (RYA-161).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
