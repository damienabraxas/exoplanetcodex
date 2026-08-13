#!/usr/bin/env python3
"""RYA-764 scope 3 — the Engine-A ∩ Engine-B overlap, re-derived.

    python3 scripts/rya764_engine_overlap.py --products <dir>

WHY THIS NUMBER MATTERS. Engines are separate products and are NEVER combined (RYA-712),
so the ONLY place they can be compared is line by line, on lines that carry both. RYA-764
measured that overlap at **6 lines** and made growing it the point of the lever: "the
cross-engine diagnostic therefore rests on 6 lines."

Scope 2-4 were reopened because they "need the Fe measurement, which is not on main".
RYA-783 landed it, so this recomputes the overlap from the per-line artifacts rather than
from an assertion.

⚠️ IT REPORTS A DIAGNOSTIC, NOT A COMBINATION. The per-line A-vs-B difference is an
RYA-525 model-atom diagnostic. It is never averaged into a value, never used to arbitrate,
and neither engine is "primary" (ratified: *"All Engines, LTE and NLTE, are products that
get presented."*).
"""
from __future__ import annotations

import argparse
import glob
import os
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "data" / "results" / "rya764"
WCOL = "wavelength_air_A"


def load_treatment(products: Path, stem: str, treatment: str) -> pd.DataFrame:
    hits = sorted(glob.glob(str(products / "**" / f"{stem}_{treatment}_lines.csv"),
                            recursive=True))
    if not hits:
        return pd.DataFrame()
    d = pd.read_csv(hits[0])
    d["_file"] = os.path.relpath(hits[0], products)
    return d


def used(d: pd.DataFrame) -> pd.DataFrame:
    """Lines that actually entered the product — in_aggregate AND an abundance."""
    if not len(d):
        return d
    m = d.get("in_aggregate", pd.Series([True] * len(d))).fillna(False).astype(bool)
    if "abundance" in d:
        m &= d.abundance.notna()
    return d[m]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--products", default=str(ROOT / "data" / "results" / "band_products"))
    ap.add_argument("--tol", type=float, default=0.05,
                    help="wavelength match tolerance in A (both sides are OUR own "
                         "measured wavelengths, so this is bookkeeping, not line ID)")
    a = ap.parse_args()
    products = Path(a.products)

    stems = sorted({os.path.basename(f).rsplit("_", 2)[0]
                    for f in glob.glob(str(products / "**" / "*_lines.csv"),
                                       recursive=True)
                    if os.path.basename(f).startswith("Fe")})
    if not stems:
        raise SystemExit(f"no per-line artifacts under {products}")

    rows = []
    print("=" * 96)
    print("RYA-764 scope 3 — ENGINE-A ∩ ENGINE-B overlap, per band")
    print("The only place two engines may be compared is line by line (RYA-712).")
    print("=" * 96)

    for stem in stems:
        a_lines = used(load_treatment(products, stem, "ENGINE-A"))
        b_lte = used(load_treatment(products, stem, "ENGINE-B"))
        b_nlte = used(load_treatment(products, stem, "ENGINE-B-NLTE"))
        if not len(a_lines):
            continue
        print(f"\n### {stem}")
        print(f"  ENGINE-A used       : {len(a_lines)}")
        for name, b in (("ENGINE-B", b_lte), ("ENGINE-B-NLTE", b_nlte)):
            if not len(b):
                print(f"  {name:<20}: (absent)")
                continue
            wa = a_lines[WCOL].to_numpy(float)
            wb = b[WCOL].to_numpy(float)
            pairs = []
            for i, w in enumerate(wa):
                j = int(np.argmin(np.abs(wb - w)))
                if abs(wb[j] - w) <= a.tol:
                    pairs.append((w, float(a_lines.abundance.iloc[i]),
                                  float(b.abundance.iloc[j])))
            print(f"  {name:<20}: used {len(b):3d}   OVERLAP with Engine-A = "
                  f"{len(pairs)}")
            if pairs:
                d = np.array([p[2] - p[1] for p in pairs])
                print(f"      per-line B-minus-A (RYA-525 DIAGNOSTIC, never a bar): "
                      f"median {np.median(d):+.4f}  mean {np.mean(d):+.4f}  "
                      f"sd {d.std(ddof=1) if len(d) > 1 else float('nan'):.4f}")
                print(f"      range {d.min():+.4f}..{d.max():+.4f}")
                for w, av, bv in pairs:
                    rows.append(dict(stem=stem, engine_b=name, wavelength_air_A=w,
                                     engine_a=av, engine_b_value=bv, delta=bv - av))

    if rows:
        OUT_DIR.mkdir(parents=True, exist_ok=True)
        out = OUT_DIR / "engine_a_b_overlap.csv"
        pd.DataFrame(rows).sort_values(["stem", "engine_b", WCOL]).to_csv(out, index=False)
        print(f"\nwrote {out.relative_to(ROOT)}  ({len(rows)} overlapping line-pairs)")
    print("\n  RYA-764 recorded the overlap as SIX lines and made growing it the point of")
    print("  the lever. Compare the counts above against that baseline.")


if __name__ == "__main__":
    main()
