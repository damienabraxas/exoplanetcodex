#!/usr/bin/env python3
"""The RYA-847 sweep report: which metric separates, and where does the cut fall?

    python3 scripts/rya847_constraint_sweep.py --sweep data/results/rya847_sweep

READS ONLY. Every number here comes from per-line CSVs the sweep already wrote, so this
can be re-run, corrected, and re-run again without touching Turbospectrum — the RYA-836
split between the expensive half and the reporting half, which is what let RYA-836
recover from a reporting crash that landed after 122 fits.

WHAT IT IS FOR
--------------
RYA-847 item 2 asks which lines are unconstrained, per band. Item 3 forbids setting the
threshold from one product (RYA-161) and, as refined, hands the SWEEP the choice of the
metric's SHAPE as well as its value. So this script ranks the three candidates on one
question:

    does the metric put a GAP in the same place in every band?

A metric that separates cleanly in one band and not in another cannot be a pipeline gate,
however good it looks where it was invented. That is precisely how the absolute
`red_chi2` threshold fails: RYA-342 ratified 10.0 against a bimodal solar Fe II
distribution, and by measurement a perfectly good NIR line sits at red_chi2 = 72 while
every one of the 40 near-UV lines sits between 27.7 and 999.5.

WHAT IT REFUSES TO DO
---------------------
Pick a number on its own. It prints the candidate cuts each metric's gap would support
and says which bands agree; the choice is a decision, recorded in the ticket, not a
side effect of running a script.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

#: The metrics, and which direction is WORSE. Ordered so the report reads consistently.
METRICS = {
    "sigma_A": ("higher is worse", "dex"),
    "frac_rise_weaker": ("LOWER is worse", "ratio"),
    "edge_distance_dex": ("LOWER is worse", "dex"),
    "red_chi2": ("higher is worse", "ratio"),
}


def load(sweep: Path) -> pd.DataFrame:
    """Every per-line row the sweep produced, tagged with its cell."""
    rows = []
    for f in sorted(sweep.glob("*_lines.csv")):
        d = pd.read_csv(f)
        # The cell identity is in the filename: element/ion, band bounds, treatment.
        d["cell"] = f.name.replace("_lines.csv", "")
        d["source_file"] = f.name
        rows.append(d)
    if not rows:
        raise SystemExit(f"no *_lines.csv under {sweep} — has the sweep run?")
    d = pd.concat(rows, ignore_index=True)
    for c in METRICS:
        if c not in d.columns:
            d[c] = np.nan
        d[c] = pd.to_numeric(d[c], errors="coerce")
    return d


def synthesis_rows(d: pd.DataFrame) -> pd.DataFrame:
    """Only lines with a chi2 surface behind them.

    An EW-route line carries None for every metric, and that is NOT a missing measurement
    — there is no chi2 to probe behind a curve-of-growth inversion, so the constraint
    question does not apply. Counting those as unconstrained would manufacture a defect
    (RYA-833: an absence is a hypothesis, never a conclusion).
    """
    return d[d.red_chi2.notna()].copy()


def gap_report(v: pd.Series, higher_is_worse: bool) -> dict:
    """The largest RELATIVE gap in the sorted values, and where it sits.

    A gate wants a threshold in a region the data does not populate. The widest
    multiplicative gap is the honest way to find one without choosing it: it makes no
    assumption about how many lines should fail, which is what a quantile would smuggle in.
    """
    s = v.dropna().sort_values().to_numpy()
    if s.size < 4:
        return {}
    s = s[s > 0]
    if s.size < 4:
        return {}
    ratios = s[1:] / s[:-1]
    i = int(np.argmax(ratios))
    return {"n": int(s.size), "min": float(s[0]), "max": float(s[-1]),
            "median": float(np.median(s)),
            "gap_lo": float(s[i]), "gap_hi": float(s[i + 1]),
            "gap_factor": float(ratios[i]),
            "cut_would_be": float(np.sqrt(s[i] * s[i + 1])),   # geometric midpoint
            "n_worse_than_cut": int((s > np.sqrt(s[i] * s[i + 1])).sum()
                                    if higher_is_worse
                                    else (s < np.sqrt(s[i] * s[i + 1])).sum())}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--sweep", default="data/results/rya847_sweep")
    ap.add_argument("--out", default=None, help="write the per-line union here")
    a = ap.parse_args(argv)
    sweep = Path(a.sweep)

    d = load(sweep)
    syn = synthesis_rows(d)
    print(f"{len(d)} per-line rows across {d.cell.nunique()} cells; "
          f"{len(syn)} carry a chi2 surface (synthesis-fit), "
          f"{len(d) - len(syn)} are EW-route (metrics not applicable)\n")

    print("PER-CELL COVERAGE")
    print(f"  {'cell':58s} {'n':>4s} {'synth':>6s} {'in_agg':>7s} {'sigma_A NaN':>12s}")
    for cell, g in d.groupby("cell"):
        gs = synthesis_rows(g)
        print(f"  {cell[:58]:58s} {len(g):4d} {len(gs):6d} "
              f"{int(g.in_aggregate.astype(str).str.lower().eq('true').sum()):7d} "
              f"{int(gs.sigma_A.isna().sum()):12d}")

    if syn.empty:
        print("\nno synthesis-fit rows yet — nothing to rank")
        return 0

    print("\n" + "=" * 78)
    print("DOES THE METRIC PUT A GAP IN THE SAME PLACE IN EVERY BAND? (item 3)")
    print("=" * 78)
    for m, (direction, unit) in METRICS.items():
        hw = "LOWER" not in direction
        print(f"\n{m}  ({direction}, {unit})")
        cuts = []
        for cell, g in syn.groupby("cell"):
            r = gap_report(g[m], hw)
            if not r:
                print(f"    {cell[:52]:52s} too few finite values to rank")
                continue
            cuts.append(r["cut_would_be"])
            print(f"    {cell[:52]:52s} range {r['min']:.4g}-{r['max']:.4g}  "
                  f"widest gap x{r['gap_factor']:.2f} at {r['cut_would_be']:.4g} "
                  f"({r['n_worse_than_cut']}/{r['n']} worse)")
        if len(cuts) >= 2:
            spread = max(cuts) / min(cuts)
            verdict = ("TRANSFERS — the per-band cuts agree within a factor of 2"
                       if spread <= 2.0 else
                       f"DOES NOT TRANSFER — per-band cuts span x{spread:.1f}")
            print(f"    -> {verdict}")

    print("\n" + "=" * 78)
    print("THE QUIET FAILURES: accepted into an aggregate, worst-constrained first")
    print("=" * 78)
    agg = syn[syn.in_aggregate.astype(str).str.lower().eq("true")]
    cols = ["cell", "wavelength_air_A", "abundance", "sigma_A",
            "frac_rise_weaker", "red_chi2", "edge_distance_dex"]
    worst = agg.sort_values("sigma_A", ascending=False).head(25)
    print(worst[cols].to_string(index=False, float_format=lambda x: f"{x:.4g}"))
    print(f"\n  sigma_A NaN (not measurable) AND in an aggregate: "
          f"{int(agg.sigma_A.isna().sum())} — these are railed or flat fits carrying a "
          f"number they did not earn")

    if a.out:
        p = Path(a.out)
        p.parent.mkdir(parents=True, exist_ok=True)
        d.to_csv(p, index=False)
        print(f"\nwrote {p}")
    print("\nNO THRESHOLD IS SET HERE. The cut is a decision recorded in RYA-847, chosen "
          "from\nthe agreement above and not from whichever band happened to look "
          "cleanest (RYA-161).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
