#!/usr/bin/env python3
"""The stage-4 GoF screen, re-runnable — RYA-992.

    python3 scripts/rya992_gof_screen.py [--anchor rya984_graded_163]

Every number in `pipeline/synth_gof.py`'s docstring comes out of here. It is a SCRIPT and
not a table because the fixture grows (the HARPS deep leg, RYA-977's +65), and a
correlation that was true on 163 lines and is quoted forever is how a statistic outlives
its evidence.

🔒 RYA-161: the target is |A - MEDIAN A| — the anchor's dispersion about ITSELF. No
published abundance enters. A statistic tuned to predict distance-from-gold would be
selecting for agreement with a number that carries its own zero point.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from pipeline.anchor_pools import load                      # noqa: E402
from pipeline.synth_gof import synth_gof_cut                # noqa: E402

CANDIDATES = ("frac_rise_weaker", "red_chi2", "sigma_A", "edge_distance_dex",
              "cited_sigma_dex", "ep_eV")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--anchor", default="rya984_graded_163")
    # The cut is a property of the ARM (RYA-992 follow-up), so the screen has to say which
    # arm its anchor was measured on rather than quoting one project-wide number.
    ap.add_argument("--arm", default="kpno_solar_atlas",
                    help="instrument the anchor's lines were measured through")
    a = ap.parse_args(argv)
    cut = synth_gof_cut(a.arm)

    df = load(a.anchor).copy()
    df["err"] = (df.abundance - df.abundance.median()).abs()
    print(f"\nRYA-992 stage-4 GoF screen — anchor {a.anchor}, n={len(df)}")
    print(f"target |A - median A|: median {df.err.median():.4f}  p90 {df.err.quantile(.9):.4f}\n")

    print(f"{'statistic':20s} {'rho':>8s} {'p':>10s}   verdict")
    for c in CANDIDATES:
        if c not in df.columns:
            continue
        v = pd.to_numeric(df[c], errors="coerce")
        m = v.notna() & df.err.notna()
        if m.sum() < 10:
            continue
        rs, ps = spearmanr(v[m], df.err[m])
        print(f"{c:20s} {rs:+8.3f} {ps:10.3g}   "
              f"{'PREDICTS' if ps < 0.05 else 'does not discriminate'}")

    # THE CONFOUND. Part of the error is gf scatter no fit statistic can predict; leaving it
    # in flatters nothing and hides that red_chi2 carries signal (correcting RYA-981).
    print("\ngf-confound removed (cited_sigma_dex regressed out of the target):")
    x = df.cited_sigma_dex.values.astype(float)
    resid = df.err.values - np.polyval(np.polyfit(x, df.err.values, 1), x)
    for c in ("frac_rise_weaker", "red_chi2", "sigma_A"):
        v = pd.to_numeric(df[c], errors="coerce").values
        rs, ps = spearmanr(v, resid)
        print(f"  {c:20s} rho {rs:+.3f}  p {ps:.4g}"
              f"{'   still predicts' if ps < 0.05 else ''}")

    f = pd.to_numeric(df.frac_rise_weaker, errors="coerce")
    print("\ndistribution structure of the winner (the threshold must sit on a break):")
    for lo, hi in ((0, 1), (1, 2), (2, 3), (3, 5), (5, 8), (8, 15), (15, np.inf)):
        m = (f >= lo) & (f < hi)
        if m.sum():
            print(f"  frac_rise [{lo:>3},{'inf' if hi == np.inf else hi:>4}): "
                  f"n={int(m.sum()):3d}  median |err| {df.err[m].median():.4f}")

    lo_, hi_ = df.err[f < cut], df.err[f >= cut]
    u, p = mannwhitneyu(lo_, hi_, alternative="greater")
    keep = df[f >= cut]
    print(f"\nAT THE ADOPTED CUT FOR {a.arm} — frac_rise >= {cut:.3f}:")
    print(f"  rejected n={len(lo_)} ({len(lo_)/len(df):.1%})  median |err| {lo_.median():.4f}"
          f"  scatter {lo_.std(ddof=1):.4f}")
    print(f"  kept     n={len(hi_)}          median |err| {hi_.median():.4f}"
          f"  scatter {hi_.std(ddof=1):.4f}")
    print(f"  Mann-Whitney p = {p:.4g}")
    print(f"  anchor scatter {df.abundance.std(ddof=1):.4f} -> "
          f"{keep.abundance.std(ddof=1):.4f}")
    print("\n  🔒 no reference abundance enters this screen (RYA-161).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
