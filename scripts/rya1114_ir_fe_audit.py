"""
scripts/rya1114_ir_fe_audit.py — RYA-1114
=========================================
The IR Fe xi-aware uncertainty re-verification, as an executable audit.

READ-ONLY. It moves no value and writes no product; it reads the live feed and the
modules that BUILD a budget, and reports what they actually contain.

WHY IT RECONSTRUCTS THE BUDGET RATHER THAN READING syst_dex
-----------------------------------------------------------
`sigma_syst` in the feed is one float. The question this audit exists to answer is
WHICH TERMS ARE IN IT, and a float cannot be asked that. So the budget is rebuilt
through `pipeline.error_budget.build()` — the same call `derive_band_products.py`
makes — and every term is printed with its measured/unmeasured state. A term whose
`dex` is None contributes NOTHING to `total()` and leaves no trace in the published
number: that silence is the finding, and it is only visible from the term list.

THE THREE THINGS THIS AUDIT CANNOT DO, AND WHY
-----------------------------------------------
🔴 Per-line work (residual-telluric per line, model-domain per line) needs a
`_lines.csv`. NONE of the six live IR products has one. The only per-line NIR data in
the repo (`data/products/solar/Fe_perline.csv`, 34 rows, kpno only) was generated
2026-08-18 from `rya847/gated` + `rya877`, a week BEFORE the live NIR artifacts
(mtime 2026-08-25), and its 34 rows reconcile with no live product's line count. So
per-line results below are labelled INDICATIVE: they describe the KP NIR region, not
the published pool.
"""
from __future__ import annotations

import csv
import json
import math
import pathlib
import sys

REPO = pathlib.Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

FEED = REPO / "data" / "products" / "solar" / "Fe.json"
PERLINE = REPO / "data" / "products" / "solar" / "Fe_perline.csv"
XI_HARNESS = REPO / "data" / "results" / "rya1120" / "xi_sigma_reported.json"
OUT = REPO / "data" / "audit" / "rya1114_ir_fe_audit"

#: The bound a holding-to-holding disagreement may not exceed (RYA-1133).
CONTINUUM_BOUND_DEX = 0.10
#: The dig-in threshold this ticket sets.
DIG_IN_DEX = 0.10
#: NIR synthesis half-width, config/synth_bands.yaml. A line is measured over
#: +/- this, so telluric-band membership of the CENTRE is not the whole question.
NIR_HALF_WIDTH_A = 1.40


def ir_products() -> list[dict]:
    d = json.loads(FEED.read_text())
    return [x for x in d["products"] if x.get("band") == "NIR"]


def rebuild_budget(p: dict) -> list[dict]:
    """The same `build()` call the deriver makes, with every term exposed."""
    import dataclasses
    from pipeline.error_budget import build
    from pipeline import harness_residual
    raw_scatter = p["sigma_stat"] * math.sqrt(p["n_lines"])
    b = build("Fe", 10141.0, p["n_lines"], scatter_dex=raw_scatter, gf_graded=True,
              **harness_residual.for_handler("SynthesisHandler").budget_kwargs())
    return [{f.name: getattr(t, f.name) for f in dataclasses.fields(t)} for t in b.terms]


def xi_coverage(p: dict) -> str:
    """Is this product in the RYA-1120 per-pool dA/dxi harness the ticket requires?"""
    if not XI_HARNESS.exists():
        return "HARNESS-ABSENT"
    h = json.loads(XI_HARNESS.read_text())
    key = (round(p["A"], 3), p["n_lines"], p.get("holding"))
    for q in h["products"]:
        if (round(q["A"], 3), q["n_lines"], q.get("holding")) == key:
            return f"MEASURED sigma_xi={q['sigma_xi']:.4f}"
    # ⚠️ ABSENT does NOT mean this arm was skipped. The xi campaign is BAND-PINNED:
    # run_campaign.py:132 hardcodes `--lo 4200 --hi 6910` on every unit, so all 15
    # pools measured VIS only. IAG WAS in the campaign (3 of the 15 units). What is
    # missing is the BAND axis, not the arm.
    return "ABSENT — campaign is band-pinned to VIS (run_campaign.py:132); arm not skipped"


def perline_indicative() -> dict:
    """Part E-a / E-c on the ONLY per-line NIR data there is. INDICATIVE, not the pool."""
    from pipeline.telluric_policy import TELLURIC_BANDS
    rows = [r for r in csv.DictReader(
        l for l in PERLINE.read_text().splitlines() if not l.startswith("#"))]
    nir = [r for r in rows if r["arm"] == "NIR"]
    inside, window_overlap = [], []
    for r in nir:
        lam = float(r["wavelength_air_A"])
        for lo, hi, name in TELLURIC_BANDS:
            if lo <= lam <= hi:
                inside.append((lam, name, r["status"]))
            elif lam + NIR_HALF_WIDTH_A > lo and lam - NIR_HALF_WIDTH_A < hi:
                window_overlap.append((lam, name, r["status"],
                                       round(min(abs(lam - lo), abs(lam - hi)), 3)))
    lab = [r for r in nir if "LAB" in r["gf_source"].upper()]
    return {
        "n_rows": len(nir),
        "centre_inside_telluric_band": inside,
        "window_overlaps_telluric_band": window_overlap,
        "lab_lines": len(lab),
        "lab_fraction_pct": round(100.0 * len(lab) / len(nir), 1) if nir else None,
        "ungraded": sum(1 for r in nir if r["gf_grade"] == "ungraded"),
    }


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    prods = ir_products()

    rows = []
    for p in prods:
        terms = rebuild_budget(p)
        unmeasured = [t["name"] for t in terms if t["dex"] is None]
        syst_rebuilt = math.sqrt(sum(
            t["dex"] ** 2 for t in terms
            if t["dex"] is not None and not t["averages_down"]))
        rows.append({
            "instrument": p["instrument"],
            "holding": p.get("holding"),
            "treatment": p["treatment"],
            "display": p.get("display"),
            "A": p["A"],
            "n_lines": p["n_lines"],
            "n_excluded": p.get("n_excluded"),
            "sigma_stat_SE": p["sigma_stat"],
            "raw_scatter_dex": round(p["sigma_stat"] * math.sqrt(p["n_lines"]), 4),
            "sigma_syst_published": p["sigma_syst"],
            "sigma_syst_rebuilt_graded_bound": round(syst_rebuilt, 4),
            "budget_terms_unmeasured": "|".join(unmeasured) or "none",
            "xi_dA_dxi": xi_coverage(p),
            "dig_in_gt_0.1": "YES" if p["sigma_stat"] > DIG_IN_DEX else "no",
            "per_line_artifact": "ABSENT",
        })

    tbl = OUT / "ir_fe_audit_table.csv"
    with tbl.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    spreads = {}
    for treat in sorted({r["treatment"] for r in rows}):
        vals = {r["instrument"]: r["A"] for r in rows if r["treatment"] == treat}
        spread = max(vals.values()) - min(vals.values())
        spreads[treat] = {
            "values": vals, "spread_dex": round(spread, 4),
            "verdict": ("BREACH" if spread > CONTINUUM_BOUND_DEX else "within")
                       + f" the {CONTINUUM_BOUND_DEX} bound (RYA-1133)"}

    summary = {
        "n_live_ir_products": len(rows),
        "holding_spread_by_treatment": spreads,
        "per_line_indicative_kpno_only": perline_indicative(),
    }
    (OUT / "ir_fe_audit_summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    print(f"wrote {tbl.relative_to(REPO)}  ({len(rows)} live IR products)")
    for k, v in spreads.items():
        print(f"  {k:9s} spread {v['spread_dex']:.4f} — {v['verdict']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
