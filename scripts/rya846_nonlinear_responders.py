#!/usr/bin/env python3
"""
RYA-846 — the 12 near-UV lines whose continuum response is not a derivative
==========================================================================
RYA-841 perturbed the continuum under all 40 near-UV lines and found 12 whose response is
NOT linear (r^2 <= 0.98), several of them running in the physically impossible direction —
a continuum placed HIGHER giving a LOWER abundance. Two return slopes above +36 at r^2 of
0.53 and 0.82, which are not derivatives of anything.

The ticket asks whether those are a continuum problem or a line problem. They are a line
problem, and the evidence is that they differ from the other 28 in a way the continuum
cannot explain:

    property                 linear (28)   non-linear (12)   permutation p
    median A(Fe I)               7.413          7.880           0.0029
    median reduced chi2           81.6          136.6           0.059
    median theoretical depth     0.898          0.899             --

They sit **+0.467 dex high**, significantly, at IDENTICAL line strength. A continuum
misplacement cannot do that selectively to 12 of 40 lines while leaving their depths alone,
so whatever they are, they are not a continuum artifact.

⚠️ WHAT THE MECHANISM IS, IS NOT ESTABLISHED. The natural story is unmodelled opacity in
the window: it would explain the high abundance, the worse fits, and the absent derivative
at once, because if the residual is dominated by something else then moving iron barely
moves it. But the blend metrics available do NOT confirm it — core dominance runs
0.774 -> 0.529 at p = 0.26, and depth, excitation potential and window occupancy are flat.
So the finding is that these lines are a distinct population that the continuum cannot
explain; naming the cause needs per-line work this ticket does not do.

WHAT THIS DOES AND DOES NOT DECIDE
  * The continuum term is computed from the lines where a derivative EXISTS. Averaging a
    slope over lines that have no slope is not a measurement (the ticket: "do not average
    them in"). The linear-only lever is reported here.
  * It does NOT drop them from the near-UV product. That is an RYA-844 adjudication —
    non-physical lines are excluded with a stated PER-LINE physical reason and the appendix
    carries the burden of proof — and "its continuum response is not a derivative" is
    evidence for that case, not the case itself. RYA-161/777 forbid tightening a bar by
    removing lines without one.
  * Every line is listed with its evidence either way (RYA-711: quarantine, never cull).

⚠️ The tempting move is visible and is NOT taken: dropping the 12 moves the near-UV product
7.488 -> 7.413 and its scatter 0.413 -> 0.347. That is exactly the shape RYA-777 warns
about, so it is REPORTED as what would happen, not applied.

Usage:
    python3 scripts/rya846_nonlinear_responders.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

SENS = ROOT / "data" / "results" / "rya841" / "rya841_continuum_sensitivity_per_line.csv"
BLEND = ROOT / "data" / "results" / "rya841" / "rya841_scatter_driver_per_line.csv"
OUT = ROOT / "data" / "results" / "rya846"
PER_LINE = OUT / "rya846_responder_diagnosis.csv"
SUMMARY = OUT / "rya846_responders.json"

#: RYA-841's linearity threshold, reused unchanged so the two tickets partition the same
#: 40 lines the same way.
R2_LINEAR = 0.98

#: A physically impossible response: placing the continuum HIGHER must make lines deeper
#: and the abundance larger. A negative slope is not a small error, it is the wrong sign.
IMPOSSIBLE_SLOPE = 0.0


def per_line_slopes() -> pd.DataFrame:
    d = pd.read_csv(SENS)
    rows = []
    for wv, g in d.groupby("wave_A"):
        g = g[np.isfinite(g.a_synth)].sort_values("delta")
        if len(g) < 5:
            continue
        x = g.delta.to_numpy(float)
        y = g.a_synth.to_numpy(float)
        slope, icept = np.polyfit(x, y, 1)
        pred = slope * x + icept
        ss_tot = float(((y - y.mean()) ** 2).sum())
        r2 = 1.0 - float(((y - pred) ** 2).sum()) / ss_tot if ss_tot > 0 else np.nan
        base = g[g.delta == 0.0]
        rows.append({
            "wave_A": float(wv),
            "slope": float(slope),
            "r2": float(r2),
            "red_chi2": float(g.red_chi2.median()),
            "a_unperturbed": float(base.a_synth.iloc[0]) if len(base) else np.nan,
            "theo_depth": float(g.theo_depth.iloc[0]),
            "ep_eV": float(g.ep_eV.iloc[0]),
            "loggf": float(g.loggf.iloc[0]),
        })
    t = pd.DataFrame(rows).sort_values("wave_A").reset_index(drop=True)
    t["linear"] = t.r2 > R2_LINEAR
    t["impossible_sign"] = t.slope < IMPOSSIBLE_SLOPE
    return t


def permutation_p(a: np.ndarray, b: np.ndarray, *, seed: int = 0,
                  n: int = 20000) -> tuple[float, float]:
    """Median difference and its two-sided permutation p. Never quote one without the
    other — RYA-841 found a +0.117 dex 'explanation' that was p = 0.28."""
    rng = np.random.default_rng(seed)
    obs = float(np.median(a) - np.median(b))
    comb = np.concatenate([a, b])
    na = len(a)
    perm = np.empty(n)
    for i in range(n):
        rng.shuffle(comb)
        perm[i] = np.median(comb[:na]) - np.median(comb[na:])
    return obs, float((np.abs(perm) >= abs(obs)).mean())


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    if not SENS.exists():
        raise SystemExit(f"RYA-841 sensitivity artifact absent: {SENS}")

    t = per_line_slopes()
    if BLEND.exists():
        b = pd.read_csv(BLEND)
        b = b[b.pool == "rya832"][["wave_A", "depth_frac_core", "n_in_core",
                                   "blend_depth", "nn_sep_A"]]
        t = t.merge(b, on="wave_A", how="left")

    lin, non = t[t.linear], t[~t.linear]
    print(f"=== {len(non)} non-linear of {len(t)} lines (r2 <= {R2_LINEAR}) ===")
    print(f"    of those, {int(non.impossible_sign.sum())} respond in the physically "
          f"IMPOSSIBLE direction (slope < 0)")

    comparisons = {}
    print("\n=== how the two populations differ ===")
    for col, label in (("a_unperturbed", "median A(Fe I)"),
                       ("red_chi2", "median reduced chi2"),
                       ("theo_depth", "median theoretical depth"),
                       ("ep_eV", "median excitation potential"),
                       ("depth_frac_core", "median core dominance")):
        if col not in t.columns or t[col].isna().all():
            continue
        a = non[col].dropna().to_numpy(float)
        c = lin[col].dropna().to_numpy(float)
        if len(a) < 3 or len(c) < 3:
            continue
        obs, p = permutation_p(a, c)
        comparisons[col] = {"non_linear_median": float(np.median(a)),
                            "linear_median": float(np.median(c)),
                            "difference": obs, "permutation_p": p}
        flag = "  <-- significant" if p < 0.05 else ""
        print(f"  {label:<28} {np.median(c):8.3f} -> {np.median(a):8.3f}   "
              f"diff {obs:+.3f}  p={p:.4f}{flag}")

    # The lever, restricted to lines that HAVE one.
    lever_all = float(np.median(np.abs(t.slope)))
    lever_lin = float(np.median(np.abs(lin.slope)))
    print(f"\n=== the lever, restricted to lines where a derivative exists ===")
    print(f"  all {len(t)} lines            : {lever_all:.2f}")
    print(f"  the {len(lin)} linear responders : {lever_lin:.2f}   <- use this")
    print(f"  averaging a slope over lines that have no slope is not a measurement")

    # What dropping them WOULD do — reported, never applied.
    drop = {"median_all": float(t.a_unperturbed.median()),
            "median_linear_only": float(lin.a_unperturbed.median()),
            "scatter_all": float(t.a_unperturbed.std(ddof=1)),
            "scatter_linear_only": float(lin.a_unperturbed.std(ddof=1))}
    print(f"\n⚠️  IF the 12 were dropped (NOT done here — RYA-844 adjudication):")
    print(f"    product  {drop['median_all']:.3f} -> {drop['median_linear_only']:.3f}"
          f"   ({drop['median_linear_only']-drop['median_all']:+.3f} dex)")
    print(f"    scatter  {drop['scatter_all']:.3f} -> {drop['scatter_linear_only']:.3f}")
    print(f"    That is the RYA-777 shape: a tighter bar bought by removing lines. It "
          f"needs a per-line physical reason (RYA-844), not a slope.")

    print(f"\n=== the 12, carried with their evidence (RYA-711) ===")
    cols = [c for c in ("wave_A", "slope", "r2", "red_chi2", "a_unperturbed",
                        "theo_depth", "depth_frac_core") if c in non.columns]
    print(non.sort_values("red_chi2", ascending=False)[cols].round(3).to_string(index=False))

    t.to_csv(PER_LINE, index=False)
    SUMMARY.write_text(json.dumps({
        "ticket": "RYA-846",
        "measurement": "diagnosis of the non-linear continuum responders",
        "r2_threshold": R2_LINEAR,
        "n_lines": int(len(t)), "n_non_linear": int(len(non)),
        "n_impossible_sign": int(non.impossible_sign.sum()),
        "comparisons": comparisons,
        "lever_all_lines": lever_all,
        "lever_linear_only": lever_lin,
        "if_dropped": drop,
        "verdict": ("NOT a continuum problem: the non-linear responders sit +0.467 dex "
                    "high at p=0.003 with identical line strength, which a continuum "
                    "misplacement cannot do selectively to 12 of 40 lines. The MECHANISM "
                    "is not established — the blend metrics do not discriminate "
                    "(core dominance p=0.26) — so these are candidate problem-children "
                    "for RYA-844 adjudication, NOT dropped here, and NOT averaged into "
                    "the continuum term."),
        "non_linear_lines": non.wave_A.round(3).tolist(),
    }, indent=2))
    print(f"\n[out] {PER_LINE}\n[out] {SUMMARY}")


if __name__ == "__main__":
    main()
