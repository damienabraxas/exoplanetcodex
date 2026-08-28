#!/usr/bin/env python3
"""
scripts/rya311_xi_ceiling_sensitivity.py
========================================
RYA-311 — how much of the solar xi answer is the SATURATION CEILING rather than the data.

    python3 scripts/rya311_xi_ceiling_sensitivity.py --pool all

🔴 WHY THIS RUN EXISTS. The FITEXY solve returns a formal Delta-chi2 = 1 error on xi, and
that error answers exactly one question: how well the slope is pinned GIVEN the stated EW
errors. It says nothing about the line selection. The reduced-EW solve has one selection
knob that is not a matter of taste -- the equivalent-width ceiling above which solar Fe I
stops responding to microturbulence and starts responding to van der Waals damping -- and
if xi moves further under a defensible change of that knob than its own formal error, then
quoting the formal error alone UNDERSTATES delta_xi. That is a measurement, so it is made
rather than argued.

The ceilings scanned are not invented for this script: 100 mA is the production constant
(RYA-330), 150 mA is the value RYA-196 originally proposed, and 80/120 bracket them.

Output: one artifact carrying every leg's xi, its formal error, the OLS root beside it and
the spread across ceilings, so the comparison "formal error vs selection sensitivity" can
be read off rather than reconstructed.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

# Sibling import: `python3 scripts/<name>.py` puts `scripts/` on sys.path[0], which is the
# same mechanism `derive_band_products` uses for `line_accounting_rya709`. Importing the
# solve rather than reimplementing it is the point -- two copies of "the xi solve" would be
# free to drift, and the sensitivity leg must be the SAME measurement at a different cut.
from rya311_solar_xi_fitexy import (EW_CEILING_MA, OUT_DIR,             # noqa: E402
                                    XI_MAX, XI_MIN, XI_STEP, solve)

#: RYA-330 production, RYA-196's original proposal, and two brackets.
CEILINGS_MA = (80.0, 100.0, 120.0, 150.0)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--star", default="solar")
    ap.add_argument("--element", default="Fe")
    ap.add_argument("--ion", default="I")
    ap.add_argument("--pool", default="all", choices=("graded", "all"))
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    grid = np.round(np.arange(XI_MIN, XI_MAX + 0.5 * XI_STEP, XI_STEP), 3)
    legs = []
    for c in CEILINGS_MA:
        res, _ = solve(a.star, a.element, a.ion, a.pool, c, grid)
        legs.append({"ew_ceiling_mA": c,
                     "root_found": bool(res.get("root_found", False)),
                     "n_fitted": res.get("n_fitted"),
                     "xi_measured": res.get("xi_measured"),
                     "delta_xi_formal": res.get("delta_xi_formal"),
                     "delta_xi_chi2_scaled": res.get("delta_xi_chi2_scaled"),
                     "xi_ols_root": res.get("xi_ols_root"),
                     "A_at_measured_xi": res.get("A_at_measured_xi"),
                     "verdict": res.get("verdict", "solved")})

    solved = [l for l in legs if l["root_found"]]
    xis = [l["xi_measured"] for l in solved]
    spread = (max(xis) - min(xis)) if len(xis) > 1 else None
    prod = next((l for l in solved if l["ew_ceiling_mA"] == EW_CEILING_MA), None)
    out = {
        "ticket": "RYA-311",
        "what": "xi vs the EW saturation ceiling -- selection sensitivity beside the "
                "formal FITEXY error",
        "star": a.star, "species": f"{a.element} {a.ion}", "pool": a.pool,
        "production_ceiling_mA": EW_CEILING_MA,
        "legs": legs,
        "xi_spread_across_ceilings": spread,
        "delta_xi_formal_at_production_ceiling": (prod or {}).get("delta_xi_formal"),
        "selection_over_formal": (
            spread / prod["delta_xi_formal"]
            if spread and prod and prod.get("delta_xi_formal") else None),
    }
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    path = Path(a.out) if a.out else OUT_DIR / f"rya311_xi_ceiling_sensitivity_{a.pool}.json"
    path.write_text(json.dumps(out, indent=2))
    print(f"\n  ceiling  n    xi        delta_xi(formal)   xi_OLS")
    for l in legs:
        if l["root_found"]:
            print(f"  {l['ew_ceiling_mA']:6.0f}  {l['n_fitted']:3d}  "
                  f"{l['xi_measured']:.4f}   {l['delta_xi_formal']:.4f}         "
                  f"{l['xi_ols_root']}")
        else:
            print(f"  {l['ew_ceiling_mA']:6.0f}   --  NO SOLVE")
    if spread is not None:
        print(f"\n  xi spread across ceilings: {spread:.4f} km/s vs a formal error of "
              f"{out['delta_xi_formal_at_production_ceiling']:.4f} "
              f"({out['selection_over_formal']:.1f}x)")
    print(f"  -> {path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
