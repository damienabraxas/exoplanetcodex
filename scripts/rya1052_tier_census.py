#!/usr/bin/env python3
"""RYA-1052 — the tier census: which (species x band) cells CAN hold a graded product.

Ryan, 2026-08-26: *"I want a double check there are no graded lines in those bands, and
then we check deep grade next. We still go through the motions even if we say hit a
non blocker like this."*

That is the whole design. An empty tier is a RESULT and has to be MEASURED, published and
reason-coded like any other — not inferred from a run that failed to start. This script
produces that measurement for every band we define and every band CRIRES+ reaches.

🔴 IT USES THE PRODUCT PATH'S OWN DEPTH, NOT A SECOND OPINION. The GRADED/DEEPGRADED split
is `_feature_depth` from `derive_band_products` — the measured `central_depth` of the
blended FEATURE a line belongs to, grouped by `GROUP_A`, exactly as
`line_accounting_rya709.features()` does it. A census that computed depth any other way
would report on a different population than the one the harness actually gates, and the
whole point is to know what a run WOULD find before spending a night on it.

⚠️ `predicted_depth` in the accounting table is NOT this quantity. An earlier pass of mine
used it and got 0 deep NIR lines; it is a theoretical column and the gate is measured.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from derive_band_products import _feature_depth          # noqa: E402
from line_accounting_rya709 import DEPTH_HI, DEPTH_LO    # noqa: E402
from pipeline.gf_grades import canonical_species         # noqa: E402

OUT = ROOT / "data" / "audit" / "rya1052_tier_census"

#: Bands we define, plus the slices CRIRES+ reaches that have no band yet. The CRIRES+
#: span is 9479.3-24855.0 A, MEASURED from real pixels (QUAL==0, finite non-zero flux) in
#: RYA-805 — not from WAVELMIN/MAX, because coverage is a comb of settings with real
#: inter-band gaps. ⚠️ CRIRES+ does NOT reach red-optical: its floor is 280.3 A above that
#: band's red edge (6910-9199). Ryan 2026-08-26: keep CRIRES NIR, no red-optical product.
CRIRES_LO, CRIRES_HI = 9479.3, 24855.0
BANDS = [
    ("near-UV",      3000.0,  3780.0, True),
    ("VIS",          3780.0,  6910.0, True),
    ("red-optical",  6910.0,  9199.0, True),
    ("NIR",          9199.0, 13000.0, True),
    ("H (proposed)",15007.0, 17500.0, False),
    ("K (proposed)",17500.0, 21400.0, False),
]


def census(species: str) -> list[dict]:
    cg = canonical_species(species)
    lab = cg[cg.gf_tier.astype(str).str.contains("LAB", na=False)]
    rows = []
    for name, lo, hi, defined in BANDS:
        inband = lab[lab.wavelength_air_A.between(lo, hi)]
        n = len(inband)
        if n:
            d = _feature_depth(inband.wavelength_air_A.values.astype(float))
            shallow = int(((d >= DEPTH_LO) & (d <= DEPTH_HI)).sum())
            deep = int((d > DEPTH_HI).sum())
            faint = int((d < DEPTH_LO).sum())
            dmin, dmed, dmax = float(d.min()), float(np.median(d)), float(d.max())
        else:
            shallow = deep = faint = 0
            dmin = dmed = dmax = float("nan")
        rows.append(dict(
            species=species, band=name, lo_A=lo, hi_A=hi, band_defined=defined,
            crires_reaches=bool(hi > CRIRES_LO and lo < CRIRES_HI),
            lab_lines=n, graded_shallow=shallow, deepgraded=deep, below_floor=faint,
            depth_min=round(dmin, 3) if n else None,
            depth_median=round(dmed, 3) if n else None,
            depth_max=round(dmax, 3) if n else None,
        ))
    return rows


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rows = census("Fe I") + census("Fe II")
    df = pd.DataFrame(rows)

    print(f"RYA-1052 tier census — depth gate [{DEPTH_LO}, {DEPTH_HI}], "
          f"measured central_depth per FEATURE\n")
    hdr = (f"{'species':7} {'band':14} {'lab':>5} {'GRADED':>7} {'DEEP':>6} {'<floor':>7} "
           f"{'depth min/med/max':>22}  {'CRIRES+':>8}")
    print(hdr); print("-" * len(hdr))
    for _, r in df.iterrows():
        dm = ("—" if r.lab_lines == 0
              else f"{r.depth_min:.2f}/{r.depth_median:.2f}/{r.depth_max:.2f}")
        flag = "reaches" if r.crires_reaches else "—"
        star = "" if r.band_defined else "  *undefined band*"
        print(f"{r.species:7} {r.band:14} {r.lab_lines:5d} {r.graded_shallow:7d} "
              f"{r.deepgraded:6d} {r.below_floor:7d} {dm:>22}  {flag:>8}{star}")

    print("\nEMPTY TIERS — these are RESULTS and belong in the store's `gaps`:")
    empty = df[(df.lab_lines > 0) & ((df.graded_shallow == 0) | (df.deepgraded == 0))]
    if empty.empty:
        print("  none — every band with lab lines can fill both tiers")
    for _, r in empty.iterrows():
        which = "GRADED" if r.graded_shallow == 0 else "DEEPGRADED"
        print(f"  {r.species} {r.band}: {which} is EMPTY — {r.lab_lines} lab lines, "
              f"{r.graded_shallow} shallow / {r.deepgraded} deep "
              f"(depth median {r.depth_median})")
    nolab = df[df.lab_lines == 0]
    for _, r in nolab.iterrows():
        print(f"  {r.species} {r.band}: NO LAB GF AT ALL — neither tier is reachable")

    df.to_csv(OUT / "tier_census.csv", index=False)
    (OUT / "tier_census.json").write_text(json.dumps({
        "ticket": "RYA-1052", "depth_gate": [DEPTH_LO, DEPTH_HI],
        "depth_source": "measured central_depth per feature (_feature_depth, GROUP_A max)",
        "crires_span_A": [CRIRES_LO, CRIRES_HI],
        "crires_reaches_red_optical": False,
        "rows": rows,
    }, indent=2) + "\n", encoding="utf-8")
    print(f"\nwrote {OUT}/tier_census.csv + .json")


if __name__ == "__main__":
    main()
