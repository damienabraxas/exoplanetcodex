#!/usr/bin/env python3
"""RYA-1043 — the decision surface for `min_d_rew_dA`. SETS NOTHING.

RYA-1043 delivers a per-line saturation MEASUREMENT; RYA-1041 Step 2 sets the threshold
from it, and that step is Ryan's. This script is the bridge: it turns 617 per-line
derivatives into the keep/drop consequences of every candidate cut, so the number is
chosen against its cost rather than picked off a distribution.

🔴 IT DOES NOT PICK. There is no "recommended" row, no argmin, no default. A script that
proposed a value would be setting the gate for an entire grid column under the cover of
reporting one, and the whole reason `Thresholds` has no defaults is that a number nobody
consciously chose is indistinguishable, downstream, from one somebody measured (RYA-161).

WHAT dREW/dA MEANS. On the linear part of a line's own curve of growth, EW is proportional
to the number of absorbers and the derivative approaches 1. As the core saturates it
collapses toward 0. So the threshold answers: how much flattening still counts as linear?
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

AUDIT = ROOT / "data" / "audit" / "rya1041"
DEFAULT_COG = AUDIT / "rya1041_perline_cog_4200_6910.csv"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cog", default=str(DEFAULT_COG))
    ap.add_argument("--out", default=str(AUDIT / "rya1043_threshold_surface.json"))
    a = ap.parse_args()

    cog = Path(a.cog)
    df = pd.read_csv(cog)
    ok = df[df.status == "ok"].copy()
    d = pd.to_numeric(ok.d_rew_dA, errors="coerce")
    ok = ok[d.notna()].copy()
    ok["d_rew_dA"] = d[d.notna()]

    prov_path = cog.with_suffix(".provenance.json")
    prov = json.loads(prov_path.read_text()) if prov_path.exists() else None

    print(f"RYA-1043 threshold surface — {cog.name}\n")
    if prov is None:
        print("  ⚠️  NO PROVENANCE SIDECAR. This CSV predates RYA-1043's provenance fix, so")
        print("      the `delta` its derivatives were computed over is NOT RECORDED. The")
        print("      script's default is 0.10 dex and the run most likely used it, but that")
        print("      is an inference, not a record — and dREW/dA is a two-sided derivative")
        print("      OVER delta. Any threshold set from this file inherits that uncertainty.")
        print("      Re-run to certify it; the script now always writes the sidecar.\n")
    else:
        print(f"  delta = {prov['delta_dex']} dex (recorded)   band {prov['band_A']}\n")

    n = len(ok)
    print(f"  n = {n} lines with a measured derivative")
    q = ok.d_rew_dA.quantile([0, .05, .25, .5, .75, .95, 1]).round(4)
    print("  dREW/dA quantiles:", {f"{int(k*100)}%": v for k, v in q.items()})

    # 🔴 The physical anchors, stated so the cut is read against physics not percentiles.
    neg = int((ok.d_rew_dA < 0).sum())
    print(f"\n  UNPHYSICAL: {neg} lines have dREW/dA < 0 — EW FALLING as abundance rises.")
    print(f"    A line cannot absorb less when given more absorbers; these are blend")
    print(f"    artifacts, and any threshold above 0 removes them as a side effect.")

    rows = []
    for t in [0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.60, 0.70]:
        keep = ok[ok.d_rew_dA >= t]
        rows.append(dict(
            threshold=t, kept=len(keep), dropped=n - len(keep),
            kept_pct=round(100 * len(keep) / n, 1),
            median_blend_frac_kept=round(float(keep.blend_fraction.median()), 3) if len(keep) else None,
            median_rew_kept=round(float(keep.rew.median()), 3) if len(keep) else None,
        ))
    surf = pd.DataFrame(rows)
    print("\n  KEEP/DROP SURFACE (a line is admitted when its own dREW/dA >= threshold):")
    print(surf.to_string(index=False))

    # ── 🔴 THE COUPLING THAT MUST BE SEEN BEFORE THE NUMBER IS CHOSEN ────────────────
    # A stricter saturation gate makes the surviving pool MORE blend-dominated, not less.
    c_blend = float(ok.d_rew_dA.corr(ok.blend_fraction))
    c_rew = float(ok.d_rew_dA.corr(ok.rew))
    c_br = float(ok.blend_fraction.corr(ok.rew))
    print("\n  🔴 THE GATE IS ANTI-CORRELATED WITH PURITY:")
    print(f"    corr(dREW/dA, blend_fraction) = {c_blend:+.3f}   more blended -> looks MORE linear")
    print(f"    corr(dREW/dA, rew)            = {c_rew:+.3f}   the gate is largely a STRENGTH gate")
    print(f"    corr(blend_fraction, rew)     = {c_br:+.3f}   weak lines are blend-dominated")
    print("    The chain is physical, not a bug: the gate admits weak lines because weak")
    print("    lines ARE linear, and weak lines are disproportionately blend-dominated. But")
    print("    it means dREW/dA cannot be set in isolation -- raising it buys linearity by")
    print("    spending purity, and the tier would fill with lines whose measured EW is")
    print("    mostly somebody else's.")
    print("\n  BLEND CONTEXT — the derivative is the LINE's only if the blend is small:")
    for lo, hi in [(0.0, 0.25), (0.25, 0.5), (0.5, 0.8), (0.8, 1.01)]:
        b = ok[(ok.blend_fraction >= lo) & (ok.blend_fraction < hi)]
        if len(b):
            print(f"    blend_fraction {lo:.2f}-{hi:.2f}: n={len(b):4d}  "
                  f"median dREW/dA {b.d_rew_dA.median():+.3f}")
    print("\n  PURITY COST OF EACH CUT (share of ADMITTED lines that are mostly blend):")
    for t in (0.30, 0.40, 0.50, 0.60):
        k = ok[ok.d_rew_dA >= t]
        if len(k):
            print(f"    threshold {t:.2f}: kept {len(k):3d}   blend>0.5: {int((k.blend_fraction > 0.5).sum()):3d} "
                  f"({100 * (k.blend_fraction > 0.5).mean():.0f}%)   blend>0.8: {int((k.blend_fraction > 0.8).sum())}")
    print("\n  ⚠️ `blend_fraction` is MEASURED here and consumed NOWHERE: stage3_purity keys")
    print("     on `problem_class`, not on the blend the COG actually quantifies. So the")
    print("     confound this section documents currently has no gate to answer it.")

    out = Path(a.out)
    out.write_text(json.dumps({
        "ticket": "RYA-1043",
        "sets_no_threshold": True,
        "source": cog.name,
        "source_provenance": prov,
        "n": n,
        "quantiles": {f"p{int(k*100)}": float(v) for k, v in q.items()},
        "n_negative_derivative": neg,
        "coupling": {"corr_d_rew_dA_blend_fraction": c_blend,
                     "corr_d_rew_dA_rew": c_rew,
                     "corr_blend_fraction_rew": c_br,
                     "note": "a stricter saturation gate yields a MORE blend-dominated pool; "
                             "blend_fraction is measured but consumed by no stage"},
        "purity_cost": [
            {"threshold": t,
             "kept": int((ok.d_rew_dA >= t).sum()),
             "kept_blend_gt_0.5": int(((ok.d_rew_dA >= t) & (ok.blend_fraction > 0.5)).sum()),
             "kept_blend_gt_0.8": int(((ok.d_rew_dA >= t) & (ok.blend_fraction > 0.8)).sum())}
            for t in (0.30, 0.40, 0.50, 0.60)],
        "surface": rows,
    }, indent=2) + "\n", encoding="utf-8")
    print(f"\n  wrote {out}")
    print("\n  🔴 SETS NOTHING. `min_d_rew_dA` stays undeclared until Ryan declares it;")
    print("     `Thresholds.require()` keeps the CONSISTENT tier closed until then.")


if __name__ == "__main__":
    main()
