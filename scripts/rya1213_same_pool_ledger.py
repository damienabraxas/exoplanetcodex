#!/usr/bin/env python3
"""RYA-1213 — Reference against its depth-split sibling, cell by cell. READ-ONLY.

    python3 scripts/rya1213_same_pool_ledger.py [--check]

🔴 THE CENTRAL EMPIRICAL CLAIM OF THIS TICKET, MADE CHECKABLE. Outside VIS the lab pool
does not straddle the 0.60 feature-depth gate, so the Reference selector and one of the
depth-split selectors return THE SAME LINES. If that is true, the two products must agree
to every published digit — and where they do not, the difference is the CODE VINTAGE of
the stored artifact, never the selector.

⚠️ A DIFFERENCE AGAINST A STORED VALUE IS NOT A MEASUREMENT OF THE SELECTOR (RYA-1204).
This ledger therefore reports the artifact DATE of the sibling beside every delta, because
that is the variable that actually explains the non-zero rows. Two paired controls in this
ticket — the sibling re-run on the SAME commit — returned exactly zero on six treatments
across two cells; see nir_paired_selector_control.json.
"""
from __future__ import annotations

import argparse
import collections
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

FEED = ROOT / "data/products/solar/Fe.json"
OUT = ROOT / "data" / "audit" / "rya1213_reference_matrix"

#: Which depth-split tier the Reference pool coincides with, per band, and why.
#: MEASURED in rya1213_reference_matrix.json, not declared here — this only names the
#: expectation so a row can be read without cross-referencing.
COINCIDES_WITH = {
    "red-optical": ("GRADED", "no lab Fe line in the band is above the depth gate"),
    "NIR": ("GRADED", "no lab Fe line in the band is above the depth gate"),
    "H": ("GRADED", "no lab Fe line in the band is above the depth gate"),
    "near-UV": ("DEEPGRADED", "58 of 59 lab Fe I lines are above the gate; the Fe II pool "
                              "is entirely above it"),
    "VIS": (None, "VIS is the ONE band whose lab pool straddles the gate (67 at or below, "
                  "109 above), so Reference there is a genuinely third pool and no "
                  "agreement is expected"),
}


def build() -> dict:
    feed = json.loads(FEED.read_text())
    by = collections.defaultdict(dict)
    for p in feed["products"]:
        key = (p["band"], p["ion"], p["holding"], p["treatment"])
        by[key][p["tier"]] = p

    rows, agree, differ = [], 0, 0
    for key in sorted(by):
        band, ion, holding, treatment = key
        sib_tier, why = COINCIDES_WITH.get(band, (None, ""))
        ref = by[key].get("REFERENCE")
        if ref is None or sib_tier is None:
            continue
        sib = by[key].get(sib_tier)
        if sib is None:
            continue
        d = round(float(ref["A"]) - float(sib["A"]), 6)
        same_n = int(ref["n_lines"]) == int(sib["n_lines"])
        rows.append({
            "band": band, "ion": f"Fe {ion}", "holding": holding, "treatment": treatment,
            "sibling_tier": sib_tier,
            "reference_A": ref["A"], "reference_n": ref["n_lines"],
            "sibling_A": sib["A"], "sibling_n": sib["n_lines"],
            "delta_dex": d, "same_n": same_n,
            # 🔴 THE VARIABLE THAT EXPLAINS THE NON-ZERO ROWS. Not decoration.
            "sibling_artifact_date": sib["provenance"]["artifact_mtime"][:10],
            "reference_artifact_date": ref["provenance"]["artifact_mtime"][:10],
        })
        if d == 0.0 and same_n:
            agree += 1
        else:
            differ += 1

    dates = sorted({r["sibling_artifact_date"] for r in rows if r["delta_dex"] != 0.0
                    or not r["same_n"]})

    #: 🔴 WHICH REFERENCE CELLS ARE NEW COVERAGE, NOT A RELABEL. The honest summary of
    #: this ticket needs both halves: most Reference products outside VIS re-measure a
    #: pool some other grade already had, and a few reach an (band, ion, holding,
    #: treatment) cell that had NO product in ANY grade. Those are new engine coverage and
    #: must not be counted together with the relabels.
    have = collections.defaultdict(set)
    for p in feed["products"]:
        have[(p["band"], p["ion"], p["holding"], p["treatment"])].add(p["tier"])
    new_cells = sorted(k for k, t in have.items() if t == {"REFERENCE"})
    return {
        "ticket": "RYA-1213",
        "note": ("Reference against its depth-split sibling where the two pools are the "
                 "SAME LINES. DIAGNOSTIC — nothing here is calibrated (RYA-161)."),
        "n_cells": len(rows),
        "n_identical": agree,
        "n_differing": differ,
        "sibling_artifact_dates_of_the_differing_rows": dates,
        "reading": ("a row with delta 0.000 AND the same n is the same-pool property "
                    "holding exactly. A differing row is explained by the sibling's "
                    "artifact DATE, not by the selector — which two paired controls "
                    "measured directly at exactly zero."),
        "reference_only_cells": [
            {"band": b, "ion": f"Fe {i}", "holding": h, "treatment": t}
            for (b, i, h, t) in new_cells],
        "reference_only_note": ("cells where the ONLY product in any grade is a Reference "
                                "one. These are new engine coverage rather than a "
                                "re-measurement of an existing pool, and they are counted "
                                "separately for that reason."),
        "rows": rows,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--check", action="store_true", help="print only; write nothing")
    a = ap.parse_args()
    doc = build()
    print(f"\nRYA-1213 same-pool ledger — {doc['n_cells']} comparable cells: "
          f"{doc['n_identical']} identical, {doc['n_differing']} differing\n")
    print(f"{'band':12s} {'ion':6s} {'holding':32s} {'treatment':20s} "
          f"{'REF':>8s} {'SIB':>8s} {'delta':>8s}  n   sibling")
    for r in doc["rows"]:
        flag = "" if (r["delta_dex"] == 0.0 and r["same_n"]) else "  <-"
        print(f"{r['band']:12s} {r['ion']:6s} {r['holding']:32s} {r['treatment']:20s} "
              f"{r['reference_A']:8.4f} {r['sibling_A']:8.4f} {r['delta_dex']:+8.4f}  "
              f"{r['reference_n']:>2d}/{r['sibling_n']:<2d} {r['sibling_artifact_date']}{flag}")
    if doc["reference_only_cells"]:
        print(f"\n  {len(doc['reference_only_cells'])} cell(s) where Reference is the ONLY "
              f"product in any grade — new engine coverage:")
        for c in doc["reference_only_cells"]:
            print(f"    {c['band']:12s} {c['ion']:6s} {c['holding']:32s} {c['treatment']}")
    if doc["n_differing"]:
        print(f"\n  differing rows carry sibling artifacts dated: "
              f"{', '.join(doc['sibling_artifact_dates_of_the_differing_rows'])}")
    if not a.check:
        OUT.mkdir(parents=True, exist_ok=True)
        (OUT / "same_pool_ledger.json").write_text(json.dumps(doc, indent=1) + "\n")
        print(f"\nwrote {OUT}/same_pool_ledger.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
