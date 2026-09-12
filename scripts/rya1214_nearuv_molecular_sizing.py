#!/usr/bin/env python3
"""RYA-1214 — sizing the near-UV CNO MOLECULAR region Ryan named as a priority.

Ryan, 2026-09-12 (comment): *"PRIORITY for CNO molecular: near-UV and VIS are where the
molecular bands are strongest (C2 Swan peaks 473-516 nm, CH A-X 387-432 nm, CN violet
~388 nm, OH A-X 306-330 nm)."*

The VIS molecular route is wired and validated (`cno_synthesis.VIS_DIAGNOSTICS`, gates in
`SOLAR_VIS_GATES`) and is what this ticket RAN. The near-UV region does not exist:
`cno_synthesis.REGIONS` is `{'vis'}`. This script measures everything needed to build it,
so the decision to build it rests on numbers rather than on an estimate — and so that
"not run" is a sized item and not an omission.

WHAT IS IN BAND, AND WHAT IS NOT
--------------------------------
The near-UV band is 3000-3780 A. Of the four bands Ryan lists:

    OH A-X       3060-3300 A   IN BAND
    NH A-X       ~3360 A       IN BAND (not in Ryan's list; it is the nitrogen counterpart)
    CN violet    ~3883 A       OUT (redward of 3780)
    CH A-X       3870-4320 A   OUT (it is the VIS G-band this ticket already ran)

So near-UV molecular CNO is specifically an OXYGEN and NITROGEN opportunity, and it is the
only route we have to a molecular A(O) and a molecular A(N) that does not need a region
above 1 um.

THE WINDOWS ARE CHOSEN BY THE TARGET MOLECULE'S DOMINANCE, NOT BY EYE
---------------------------------------------------------------------
A molecular band fit is only a measurement of X if X's lines carry the window. Measured on
the RYA-1207 near-UV lists, per 3 A window (the width `VIS_DIAGNOSTICS` uses):

    3063-3066 A    OH  30   NH  10      OH-dominated 3:1   <- OH A-X (0,0) band head
    3122-3125 A    OH  26   NH   7      OH-dominated 4:1
    3358-3361 A    NH 109   OH  11      NH-dominated 10:1  <- NH A-X (0,0) band head
    3370-3373 A    NH  75   OH   8      NH-dominated 9:1

Those are real band heads with real dominance ratios, so the windows are defensible.

🔴 AND THERE IS A MEASURED, ONE-SIDED OBSTACLE THAT THE VIS BAND DOES NOT HAVE.
The near-UV synthesis is known to carry LESS opacity than the observation, and the size of
the deficit is measured, not suspected:

  * RYA-1204/1207 wired molecular opacity into the near-UV and measured the lever at
    -0.050 dex (paired) on IRON. That ticket's own conclusion is that it STILL
    under-corrects afterwards: "synth/obs above 1 on 24 of 40 windows after", which it
    states is "a floor on the missing opacity, not a correction stopped at a target".
  * RYA-1189 measured the band as blend-dominated at 59 lines/A with 0 of 10 clean
    side-bands and no isolated Fe I line existing in it at all.

A fit compensates missing opacity by RAISING the abundance. So a near-UV OH or NH band fit
returns an UPPER BOUND on A(O) / A(N), not a measurement, until the deficit is closed —
and the honest product would have to be labelled that way. That is why this script sizes
the region and does not run it: emitting a one-sidedly biased number as an abundance is
the failure this repo's gates exist to prevent, and the bias direction is known in advance.

WHAT BUILDING IT COSTS
----------------------
Four pieces, enumerated from the VIS region's own structure:

  1. a `RegionConfig` for 3000-3780 A on `kpno_solar_atlas` / `solar_kpno_kurucz2005_corrected`
  2. `NEARUV_DIAGNOSTICS` — OH A-X and NH A-X, windows as measured above
  3. a loader dispatch for the Kitt Peak atlas (the VIS arms use `harps_normalized` /
     `reflected_solar`; `measure_band_ew.load_window_ex` already reads this holding)
  4. an acceptance gate in the shape of `SOLAR_VIS_GATES`, and the opacity-deficit label
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "audit" / "rya1214_cno_products"

#: Band edges from `pipeline.band_policy` (near-UV).
BAND_LO_A, BAND_HI_A = 3000.0, 3780.0
#: The window width `cno_synthesis.VIS_DIAGNOSTICS` uses for a molecular band head.
WINDOW_A = 3.0
#: Candidate band heads. OH A-X and NH A-X vibrational sequences; NOT tuned to an answer,
#: and every candidate is reported with its dominance ratio whether it wins or loses.
CANDIDATES = {
    "OH": (3063, 3078, 3090, 3122, 3144, 3167, 3180),
    "NH": (3358, 3370, 3380, 3400, 3430),
}
#: Ryan's four named bands, and whether the near-UV band contains them.
RYAN_BANDS = (
    ("OH A-X", 3060, 3300, "O"),
    ("NH A-X", 3350, 3450, "N"),
    ("CN violet", 3860, 3900, "N"),
    ("CH A-X", 3870, 4320, "C"),
)


def _load(path: Path) -> np.ndarray:
    w = []
    with path.open() as fh:
        for ln in fh:
            parts = ln.split()
            if not parts:
                continue
            try:
                w.append(float(parts[0]))
            except ValueError:
                continue
    return np.asarray(w, float)


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    ispec = os.environ.get("ISPEC_DIR")
    if not ispec:
        raise SystemExit("ISPEC_DIR unset — the near-UV .bsyn lists are RYA-1207's, "
                         "written into the iSpec molecules dir. Run on Sirius.")
    d = Path(ispec) / "input" / "linelists" / "turbospectrum" / "molecules"
    lists = {}
    for mol in ("OH", "NH", "CH", "CN"):
        p = d / f"{mol}_300-378.bsyn"
        if p.exists():
            lists[mol] = _load(p)

    in_band = [{"band": n, "lo_A": lo, "hi_A": hi, "element": el,
                "in_near_uv": bool(lo >= BAND_LO_A and hi <= BAND_HI_A),
                "why": ("in band" if lo >= BAND_LO_A and hi <= BAND_HI_A
                        else f"redward of the {BAND_HI_A:.0f} A band edge")}
               for n, lo, hi, el in RYAN_BANDS]

    windows = []
    for mol, heads in CANDIDATES.items():
        if mol not in lists:
            continue
        for lo in heads:
            hi = lo + WINDOW_A
            counts = {m: int(((w >= lo) & (w <= hi)).sum()) for m, w in lists.items()}
            other = sum(v for k, v in counts.items() if k != mol)
            windows.append({
                "target": mol, "lo_A": lo, "hi_A": hi,
                "counts": counts,
                "dominance": (round(counts[mol] / other, 2) if other else None),
                "usable": bool(counts[mol] >= 10 and (other == 0 or counts[mol] / other >= 2.0)),
            })

    doc = {
        "ticket": "RYA-1214",
        "what": "sizing the near-UV CNO molecular region (Ryan's priority comment)",
        "status": "SIZED, NOT RUN — and the reason is measured, not a preference",
        "band_A": [BAND_LO_A, BAND_HI_A],
        "ryan_named_bands": in_band,
        "near_uv_is_an_O_and_N_opportunity": (
            "OH A-X and NH A-X are in band; CN violet (~3883 A) and CH A-X (3870-4320 A) "
            "are redward of the 3780 A edge — CH A-X IS the VIS G-band this ticket ran."),
        "lists_present": {m: {"lines": int(len(w)),
                              "span_A": [round(float(w.min()), 1), round(float(w.max()), 1)]}
                          for m, w in lists.items()},
        "candidate_windows": windows,
        "recommended": [w for w in windows if w["usable"]][:4],
        "blocker": {
            "kind": "MEASURED ONE-SIDED BIAS, not a missing input",
            "detail": "the near-UV synthesis carries LESS opacity than the observation. "
                      "RYA-1204/1207 measured the molecular lever at -0.050 dex (paired) "
                      "on iron and concluded it STILL under-corrects after: 'synth/obs "
                      "above 1 on 24 of 40 windows', 'a floor on the missing opacity, not "
                      "a correction stopped at a target'. RYA-1189 measured the band as "
                      "blend-dominated at 59 lines/A with 0 of 10 clean side-bands.",
            "consequence": "a fit compensates missing opacity by RAISING the abundance, so "
                           "a near-UV OH or NH band fit returns an UPPER BOUND on A(O) / "
                           "A(N), not a measurement. The bias direction is known in "
                           "advance, which is exactly why the number must not be emitted "
                           "as an abundance without that label.",
        },
        "to_build": [
            "RegionConfig for 3000-3780 A on kpno_solar_atlas / solar_kpno_kurucz2005_corrected",
            "NEARUV_DIAGNOSTICS: OH A-X and NH A-X at the windows above",
            "a loader dispatch for the Kitt Peak atlas (measure_band_ew.load_window_ex "
            "already reads this holding; cno_synthesis has no dispatch for it)",
            "an acceptance gate shaped like SOLAR_VIS_GATES, plus the opacity-deficit label",
        ],
    }
    (OUT / "nearuv_molecular_sizing.json").write_text(json.dumps(doc, indent=2) + "\n")

    print("=== RYA-1214 — near-UV CNO molecular region, SIZED ===")
    for b in in_band:
        print(f"  {b['band']:12s} {b['lo_A']}-{b['hi_A']} A  ({b['element']})  "
              f"{'IN BAND' if b['in_near_uv'] else 'OUT: ' + b['why']}")
    print("\n  candidate windows (3 A, the width VIS_DIAGNOSTICS uses):")
    for w in windows:
        mark = "USABLE" if w["usable"] else "      "
        print(f"    {mark} {w['target']} {w['lo_A']}-{w['hi_A']:.0f} A  "
              f"{w['counts']}  dominance={w['dominance']}")
    print(f"\n  🔴 {doc['blocker']['kind']}: {doc['blocker']['consequence'][:120]}...")
    print(f"  wrote {(OUT / 'nearuv_molecular_sizing.json').relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
