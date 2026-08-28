#!/usr/bin/env python3
"""
scripts/rya1095_amarsi_error_budget.py
======================================
RYA-1095 — build the systematic budget the Amarsi 3D-NLTE Fe I leg has never had, and
put its statistical bar on the same footing as every other product's.

    python3 scripts/rya1095_amarsi_error_budget.py --run-dir <dir with rya817_3dnlte_*>

🔴 TWO SEPARATE DEFECTS, AND THEY LOOK LIKE ONE.

1. `sigma_syst` is NULL. The EW-3D route calls `band_products.build_product`, which
   computes a value and a spread and stops; it never assembles an `error_budget`. So the
   product reaches the feed with no systematic term at all, and `A +/- sigma_stat` states
   a precision nobody computed (RYA-968).
2. `sigma_stat` IS THE WRONG STATISTIC. `build_product` publishes `np.std(vals, ddof=1)`
   -- the raw per-line scatter -- while the band route publishes a STANDARD ERROR
   (`error_budget.py:609`, scatter/sqrt(n) RMS'd with the other averaging-down terms).
   Nothing in the feed records which, so the Amarsi bar renders beside the others and
   reads about sqrt(n) times worse than it is. RYA-1092 first mistook that for anomalous
   scatter; it is a schema defect, and this is where it gets fixed for this leg.

WHAT THIS SCRIPT DOES NOT DO
----------------------------
⚠️ It does not copy a sibling's `sigma_syst`. The Amarsi leg uses a SUBSET of its 1D base
pool (50 of 67 lines survive the training-domain check), so the gf rung is re-decided on
the lines this product actually used -- `gf_rung.for_lines`, the same call the band route
makes. A pool's bar is priced by the pool's own lines (RYA-855), and a quarantined line
does not price the measurement it was quarantined out of.

⚠️ IT DOES NOT INVENT THE NETWORK'S OWN ACCURACY. The vendored package says, in its own
README, that the authors "cannot supply errors as such". So the MLP's intrinsic accuracy
is recorded as an UNMEASURED term -- `error_budget` has a first-class third state for
exactly this -- and the published `sigma_syst` is a documented FLOOR rather than a total.
Writing a plausible number there would be the RYA-907 failure this codebase keeps
catching: "measured and negligible" published for a quantity nobody determined.

WHAT IT DOES MEASURE, from our own artifacts
--------------------------------------------
* the pool's line-to-line scatter, and the standard error that follows from it;
* the gf rung, re-decided on the surviving lines;
* the harness residual for the handler that produced them;
* the 3D-NLTE AXIS term. `aberr` is the correction at the pinned abundance axis and
  `aberr_axis_line` the same correction at the network's own per-line axis -- two
  defensible choices, so the paired spread between them on the lines that have both is a
  measured component of the model uncertainty. It is small (~0.004 dex), and being small
  is a result rather than a reason to leave it out.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline import gf_rung, harness_residual                    # noqa: E402
from pipeline.band_products import LineMeasurement                # noqa: E402
from pipeline.error_budget import Term, build as build_budget     # noqa: E402
from pipeline.error_budget import assert_stat_publishable, round_dex   # noqa: E402

#: The vendor's own words. Quoted rather than paraphrased because it is the SOURCE of an
#: unmeasured term, and a term's source is the thing that makes it inspectable.
_VENDOR_NO_ERRORS = (
    "vendor/1L-3NErrors/README.md: the authors state they \"cannot supply errors as "
    "such\" for the network's own abundance corrections. Amarsi, Liljegren & Nissen 2022 "
    "(A&A 668 A68) publish no per-line accuracy we can read off, and nothing we hold "
    "measures the MLP against the 3D non-LTE calculations it was trained on -- so this "
    "term is UNMEASURED, not zero (RYA-907).")


def load_run(run_dir: Path) -> tuple:
    per_line = pd.read_csv(run_dir / "rya817_3dnlte_per_line.csv")
    products = pd.read_csv(run_dir / "rya817_3dnlte_products.csv")
    summary = json.loads((run_dir / "rya817_run_summary.json").read_text())
    return per_line, products, summary


def pool_of(per_line: pd.DataFrame, element: str, ion: str, band: str) -> pd.DataFrame:
    """The lines that ACTUALLY entered the product -- `run_3dnlte_bands`' own rule:
    in the training domain, in the 1D base's aggregate, and with a finite 3D value."""
    g = per_line[(per_line.element == element) & (per_line.ion == ion)
                 & (per_line.band == band)]
    return g[g.in_domain.map(bool) & g.in_aggregate_1dlte.map(bool)
             & np.isfinite(pd.to_numeric(g.a_3dnlte, errors="coerce"))].copy()


def axis_term(pool: pd.DataFrame) -> Term:
    """The 3D-NLTE ABUNDANCE-AXIS term, measured as a paired difference.

    The correction depends on the A(Fe) axis it is evaluated at. `aberr` uses the pinned,
    self-consistently converged axis (the product); `aberr_axis_line` uses the network's
    own per-line axis. Both are defensible, so the spread between them on the lines that
    carry both IS the size of that choice. Paired, because the two differ only in the axis
    -- everything else about the line is held still.
    """
    both = pool.dropna(subset=["aberr", "aberr_axis_line"])
    if len(both) < 3:
        return Term("3D-NLTE abundance axis", None, False,
                    f"only {len(both)} line(s) carry both the pinned-axis and per-line-"
                    f"axis correction, which is too few to measure the spread between "
                    f"them. UNMEASURED rather than assumed small (RYA-907).")
    d = (both.aberr - both.aberr_axis_line).to_numpy(float)
    rms = float(np.sqrt((d ** 2).mean()))
    return Term("3D-NLTE abundance axis", rms, False,
                f"paired |aberr(pinned axis) - aberr(per-line axis)| over the {len(both)} "
                f"lines that carry both; RMS {rms:.5f} dex. A systematic: every line in "
                f"the product is evaluated at the SAME pinned axis, so the choice shifts "
                f"them together and adding lines cannot average it away.")


def network_term() -> Term:
    return Term("3D-NLTE network accuracy", None, False, _VENDOR_NO_ERRORS)


def budget_from_pool(pool: pd.DataFrame, *, element: str, ion: str, instrument: str,
                     handler: str, scatter_dex: float):
    """The budget for an Amarsi pool held IN MEMORY. The one assembly, two callers.

    `scripts/rya817_run_3dnlte_bands.py` calls this while the run still has the pool, so
    the product it writes carries `stat_dex`/`syst_dex` like every band product; the audit
    entry point below calls it from files. Two copies of "what the Amarsi budget is" would
    be free to drift, which is the RYA-845 shape.
    """
    n = len(pool)
    measurements = [
        LineMeasurement(element=element, ion=ion,
                        wavelength_air_A=float(r.wavelength_air_A),
                        instrument=str(instrument), ew_mA=float(r.ew_mA),
                        ew_method="amarsi 3D-NLTE on the 1D-LTE base",
                        abundance=float(r.a_3dnlte), treatment="ENGINE-A-3DNLTE",
                        in_aggregate=True, ep_eV=float(r.elo_eV))
        for r in pool.itertuples()]
    linelist = pd.DataFrame({
        "element": [f"{element} {1 if ion == 'I' else 2}"] * n,
        "wave_A": pool.wavelength_air_A.to_numpy(float),
        "lower_state_eV": pool.elo_eV.to_numpy(float),
        "loggf": pool.loggf.to_numpy(float)})
    rung = gf_rung.for_lines(element, ion, measurements, linelist=linelist)
    hres = harness_residual.for_handler(str(handler))
    lo, hi = float(pool.wavelength_air_A.min()), float(pool.wavelength_air_A.max())
    b = build_budget(element, 0.5 * (lo + hi), n, scatter_dex=scatter_dex,
                     **rung.budget_kwargs(), **hres.budget_kwargs())
    b.add(axis_term(pool))
    b.add(network_term())
    return b, rung


def budget_for(run_dir: Path, *, element="Fe", ion="I", band="VIS") -> dict:
    per_line, products, summary = load_run(run_dir)
    row = products[(products.element == element) & (products.ion == ion)
                   & (products.band == band)]
    if row.empty:
        raise SystemExit(f"no {element} {ion} {band} product in {run_dir}")
    row = row.iloc[0]
    pool = pool_of(per_line, element, ion, band)
    n = len(pool)
    if n != int(row.n_lines):
        raise SystemExit(
            f"reconstructed pool is {n} lines but the product says {int(row.n_lines)}. "
            f"Refusing to price a pool I cannot reproduce -- the membership rule here "
            f"must be the one the run used, not a similar one.")

    a3 = pool.a_3dnlte.to_numpy(float)
    scatter = float(np.std(a3, ddof=1))
    published = float(row.stat_dex) if "stat_dex" in products.columns else float(row.sigma)
    b, rung = budget_from_pool(pool, element=element, ion=ion,
                               instrument=str(row.instrument),
                               handler=str(row.handler), scatter_dex=scatter)
    stat, syst = b.total()
    assert_stat_publishable(stat, cell=f"{element} {ion} {band} {row.treatment} "
                                       f"({row.instrument}, n={n})")
    return {
        "run_dir": str(run_dir), "instrument": str(row.instrument),
        "holding_hint": summary.get("band_products_dir"),
        "element": element, "ion": ion, "band": band,
        "treatment": str(row.treatment), "handler": str(row.handler),
        "A": float(row.A if "A" in products.columns else row.value), "n_lines": n,
        "published_statistical_field": published,
        "line_scatter_dex": round_dex(scatter),
        "sigma_stat_STANDARD_ERROR": round_dex(stat),
        "sigma_syst_FLOOR": round_dex(syst),
        "gf_rung": rung.describe(),
        "dominant_term": (b.dominant().name if b.dominant() else None),
        "stat_basis": b.stat_basis(),
        "unmeasured_terms": [{"name": t.name, "source": t.source}
                             for t in b.unmeasured_terms()],
        "terms": [{"name": t.name, "dex": t.dex, "averages_down": t.averages_down,
                   "source": t.source} for t in b.terms],
        "axis_railed_at_grid_ceiling": summary.get("per_band", {})
            .get(f"{element} {ion} {band}", {}).get("afe_axis_railed_at_grid_ceiling"),
        "axis_sensitivity_dex": summary.get("per_band", {})
            .get(f"{element} {ion} {band}", {}).get("afe_axis_sensitivity_dex"),
        "sklearn_version_skew": summary.get("run_environment", {})
            .get("sklearn_version_skew"),
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run-dir", action="append", required=True, type=Path)
    ap.add_argument("--element", default="Fe")
    ap.add_argument("--ion", default="I")
    ap.add_argument("--band", default="VIS")
    ap.add_argument("--out", type=Path,
                    default=ROOT / "data" / "results" / "rya1095"
                    / "rya1095_amarsi_error_budget.json")
    a = ap.parse_args()

    out = [budget_for(d, element=a.element, ion=a.ion, band=a.band) for d in a.run_dir]
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps({"ticket": "RYA-1095", "budgets": out}, indent=2) + "\n")

    print(f"\n{'instrument':22s} {'n':>3s} {'A':>8s} {'scatter':>8s} "
          f"{'stat(SE)':>9s} {'syst':>8s}  published")
    for r in out:
        print(f"{r['instrument']:22s} {r['n_lines']:3d} {r['A']:8.4f} "
              f"{r['line_scatter_dex']:8.4f} {r['sigma_stat_STANDARD_ERROR']:9.4f} "
              f"{r['sigma_syst_FLOOR']:8.4f}  {r['published_statistical_field']:.4f}")
    print(f"\n  gf rung        : {out[0]['gf_rung']}")
    print(f"  dominant term  : {out[0]['dominant_term']}")
    print(f"  UNMEASURED     : {[t['name'] for t in out[0]['unmeasured_terms']]}")
    print(f"  -> sigma_syst is a FLOOR, not a total.")
    try:
        shown = a.out.relative_to(ROOT)
    except ValueError:
        shown = a.out          # an out-of-repo scratch path is legitimate for a dry look
    print(f"\n  -> {shown}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
