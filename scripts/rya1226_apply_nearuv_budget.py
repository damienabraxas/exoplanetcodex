#!/usr/bin/env python3
"""
RYA-1226 part D -- attach the migrated budgets to the feed and re-emit.

Only rows whose budget VALIDATES are attached. A row that still holds a component keeps
its legacy shape and its legacy bar: the point of the migration is that an incomplete
budget cannot publish, and that has to stay true for the rows that are still incomplete.

🔴 THE BAR MOVES, AND IT MOVES A LOT. sigma_reported goes ~0.122 -> ~0.321 on the Fe II
1D-LTE rows -- 2.6x -- because the legacy budget charged four terms and the contract
charges sixteen. 85% of the new variance comes from BOUNDS rather than measurements, and
profile_ew alone is 45%. The abundance does not move, on any row.
"""

from __future__ import annotations

import json
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
    sys.path.insert(0, str(ROOT / "scripts"))

from pipeline.uncertainty_contract import (  # noqa: E402
    UncertaintyError, assert_publication_feed, validate)
import rya1226_migrate_nearuv_budget as M                             # noqa: E402

FEED = ROOT / "data/products/solar/Fe.json"
REPORT = ROOT / "data/results/rya1226/nearuv_budget_applied.json"


def main() -> int:
    feed = json.loads(FEED.read_text())
    part_c = json.loads((ROOT / "data/results/rya1226/nearuv_fe1_xi_dadxi.json").read_text())
    fe1 = {(r["holding"], r["treatment"]): r["dA_dxi"] for r in part_c["pools"]
           if r["xi_state"] == "MEASURED"}

    applied, skipped = [], []
    for p in feed["products"]:
        if p.get("band") != "near-UV" or p.get("tier") != "DEEPGRADED":
            continue
        if p.get("treatment") not in ("1D-LTE", "ENGINE-A"):
            continue
        xi = -0.1100 if p["ion"] == "II" else fe1.get((p["holding"], p["treatment"]))
        built = M.build(p, xi_slope=xi)
        try:
            validate(built["budget"], scope=built["scope"])
        except UncertaintyError as exc:
            skipped.append({"ion": p["ion"], "holding": p["holding"],
                            "treatment": p["treatment"], "reason": str(exc),
                            "sigma_reported_kept": p.get("sigma_reported")})
            continue
        before = p.get("sigma_reported")
        p["star"] = "solar"
        p["uncertainty_indicator_ids"] = built["evidence"]["indicator_ids"]
        p["uncertainty"] = built["budget"]
        p["sigma_reported"] = built["budget"]["sigma_reported"]
        p["sigma_reported_basis"] = (
            "RYA-587 canonical total over all 16 contract components (RYA-1226). Supersedes "
            "the legacy quadrature, which charged four terms. See `uncertainty` for the "
            "component table; terms whose evidence is a bound carry evidence.bound = true.")
        applied.append({"ion": p["ion"], "holding": p["holding"],
                        "treatment": p["treatment"], "A": p["A"],
                        "sigma_reported_before": before,
                        "sigma_reported_after": p["sigma_reported"],
                        "ratio": round(p["sigma_reported"] / before, 3) if before else None})

    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(json.dumps(
        {"ticket": "RYA-1226", "part": "D", "n_applied": len(applied),
         "n_skipped": len(skipped), "applied": applied, "skipped": skipped}, indent=2) + "\n")

    #: 🔴 THROUGH THE GATE, NEVER AROUND IT. Writing the feed directly would skip the very
    #: contract this migration exists to satisfy -- the point is that RYA-587 ACCEPTS these
    #: rows, and that is only demonstrated by asking it.
    assert_publication_feed(feed, previous=json.loads(FEED.read_text()))
    FEED.write_text(json.dumps(feed, indent=2) + "\n")
    print(f"applied {len(applied)}, skipped {len(skipped)}")
    for a in applied:
        print(f"  Fe {a['ion']:3} {a['holding'][:31]:31} {a['treatment']:9} "
              f"A={a['A']} sigma {a['sigma_reported_before']} -> "
              f"{a['sigma_reported_after']:.6f}  ({a['ratio']}x)")
    for sk in skipped:
        print(f"  SKIP Fe {sk['ion']:3} {sk['holding'][:31]:31} {sk['treatment']:9} "
              f"{sk['reason'][:60]}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
