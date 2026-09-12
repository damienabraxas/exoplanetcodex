#!/usr/bin/env python3
"""RYA-1214 — the AGSS21 Table 3 ATOMIC indicator set, per element, as a runnable line set.

Step 4 names the lines: "[C I] 8727 + the permitted C I set", "the 5 N I lines",
"[O I] 6300 + O I 777 triplet + O I 8446". This writes exactly those, from the RYA-1136
census rather than from the ticket's prose, so the set is traceable to an artifact and a
re-run cannot quietly drift from what was measured.

🔴 WHY A NAMED SET AND NOT THE PRODUCTION SELECTOR. `select_lines` ranks candidates by
theoretical central depth and refuses anything below 0.15. Measured on the GES v6 list
this band actually synthesises:

    O I  VIS          117 lines   0 above the floor   max depth 0.040 ([O I] 6300)
    N I  VIS           97 lines   0                   max depth 0.000
    N I  red-optical   68 lines   0                   max depth 0.040
    C I  VIS          768 lines   7                   max depth 0.260
    C I  red-optical  284 lines   7                   max depth 0.420 ([C I] 8727 = 0.030)
    O I  red-optical   46 lines   4                   the 777 triplet + 8446.36

So the production rule reaches 11 of the 23 atomic indicators and reaches NEITHER
forbidden line and NONE of nitrogen. The floor is not wrong — a line that shallow is not
visibly above the crowding — but the light elements' indicators are weak BECAUSE they are
unsaturated, which is exactly what makes them good diagnostics. A strength-ranked rule
selects against the lines the campaign exists to measure, and no amount of re-running it
will produce an [O I] 6300 number.

⚠️ WHAT THIS SET IS AND IS NOT. It is AGSS21's ADOPTED indicator list as the RYA-1136
census records it, restricted to rows whose `use_status` says the source ADOPTED them.
The Amarsi-2019 grid-input rows are carried too, flagged, because the census itself says
"grid input is not by itself proof of final AGSS21 adopted-line use" — an ambiguity the
census refused to collapse and this script does not collapse either. The two are written
to SEPARATE files so a product measured on one cannot be reported as the other.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
CENSUS = ROOT / "data" / "audit" / "rya1136_cno_intake" / "atomic_source_census.csv"
OUT = ROOT / "data" / "linelists" / "reference_sets"

#: The census's own `use_status` for rows AGSS21 states it adopted.
ADOPTED = "AGSS21_ADOPTED_FIVE_LINE_SET"
#: Rows that entered the SOURCE ANALYSIS grid. The census flags these as not proof of
#: final adopted-line use; kept, separated, never merged with the line above.
GRID = "SOURCE_ANALYSIS_GRID_SET"


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    cen = pd.read_csv(CENSUS)
    cen = cen[cen.element.isin(["C", "N", "O"]) &
              (cen.use_status != "BLEND_COMPONENT_NOT_A_SOURCE_LINE")].copy()

    written = []
    for el, g in cen.groupby("element"):
        for status, tag in ((ADOPTED, "adopted"), (GRID, "grid")):
            s = g[g.use_status == status]
            if s.empty:
                continue
            out = OUT / f"agss21_cno_{el}_{tag}_rya1214.csv"
            s[["element", "species", "line_label", "wavelength_air_A", "lower_EP_eV",
               "published_loggf", "source_band", "use_status",
               "reference_line_set"]].sort_values("wavelength_air_A").to_csv(out, index=False)
            written.append({"element": el, "use_status": status, "n_lines": int(len(s)),
                            "file": str(out.relative_to(ROOT)),
                            "span_A": [float(s.wavelength_air_A.min()),
                                       float(s.wavelength_air_A.max())]})
            print(f"  {el} {tag:8s} {len(s):3d} lines -> {out.relative_to(ROOT)}")

    prov = {
        "ticket": "RYA-1214",
        "what": "AGSS21 Table 3 ATOMIC indicator line sets, per element, for "
                "`derive_band_products --lines-from-set`",
        "source": str(CENSUS.relative_to(ROOT)),
        "source_ticket": "RYA-1136 / RYA-1183 (atomic_source_census.csv)",
        "why_not_select_lines": "select_lines ranks by theoretical central depth and "
                                "refuses below 0.15. That excludes BOTH forbidden "
                                "indicators ([C I] 8727 depth 0.030, [O I] 6300 depth "
                                "0.040) and ALL FIVE N I lines (max 0.040) — 12 of the 23 "
                                "atomic indicators. The light elements' indicators are "
                                "weak because they are unsaturated.",
        "adopted_vs_grid": "kept in SEPARATE files. The census states that a row's "
                           "presence in the Amarsi-2019 analysis grid is NOT proof AGSS21 "
                           "adopted it; merging the two would publish a set the paper "
                           "does not claim.",
        "files": written,
    }
    (OUT / "agss21_cno_lineset_rya1214.prov.json").write_text(json.dumps(prov, indent=2) + "\n")
    print(f"\n  wrote {(OUT / 'agss21_cno_lineset_rya1214.prov.json').relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
