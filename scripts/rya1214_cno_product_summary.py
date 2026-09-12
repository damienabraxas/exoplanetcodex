#!/usr/bin/env python3
"""RYA-1214 Step 4 — the CNO products, and WHICH AGSS21 indicators each one measured.

A(X) alone does not say what was measured. `select_lines` takes the strongest N lines by
theoretical depth inside the band, which is OUR selection rule — it is not AGSS21's
adopted set, and line selection is the known dominant lever on an abundance (RYA-842). So
a C I number that happens to sit near 8.48 has to say whether it got there on the lines
Asplund used or on a different pool that agrees by coincidence.

This joins each emitted per-line artifact to the AGSS21 Table 3 atomic indicator list and
reports the overlap per product, so the comparison is stated rather than implied.

⚠️ A PRODUCT HERE IS NOT A REPLICATION. A replication is measured ON the published line
set and published with an explicit `line_set` (RYA-1127/1185) — that is what the ratified
vocabulary calls Reference Grade, and it is a separate run, not a relabelling of this one.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
BP = ROOT / "data" / "results" / "band_products"
OUT = ROOT / "data" / "audit" / "rya1214_cno_products"
CENSUS = ROOT / "data" / "audit" / "rya1136_cno_intake" / "atomic_source_census.csv"

#: The census prints 0.01 nm; the window IS its quantisation (RYA-1109).
TOL_A = 0.1


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    cen = pd.read_csv(CENSUS)
    cen = cen[cen.element.isin(["C", "N", "O"])]

    rows = []
    for prod in sorted(BP.glob("[CNO]I_*_SYNTH_products.csv")):
        stem = prod.name[: -len("_products.csv")]
        p = pd.read_csv(prod)
        for _, r in p.iterrows():
            lines_f = BP / f"{stem}_{r.treatment}_lines.csv"
            hit, used = [], []
            if lines_f.exists():
                lf = pd.read_csv(lines_f)
                inagg = lf[lf.in_aggregate.astype(str).str.lower() == "true"]
                used = [round(float(w), 3) for w in inagg.wavelength_air_A]
                ce = cen[cen.element == r.element]
                for _, c in ce.iterrows():
                    w = float(c.wavelength_air_A)
                    if any(abs(w - u) <= TOL_A for u in used):
                        hit.append(c.line_label)
            budget = BP / f"{stem}_budgets.txt"
            rung = ""
            if budget.exists():
                for ln in budget.read_text().splitlines():
                    if ln.strip().startswith("gf rung:"):
                        rung = ln.split("gf rung:", 1)[1].strip()
            rows.append({
                "element": r.element, "ion": r.ion, "band": r.band,
                "instrument": r.instrument,
                "holding": stem.split(f"_{r.instrument}_", 1)[-1].replace("_SYNTH", ""),
                "treatment": r.treatment,
                "A": r.A, "sigma_stat": r.stat_dex, "sigma_syst": r.syst_dex,
                "n_lines": r.n_lines, "n_excluded": r.n_excluded,
                "dominant_term": r.dominant,
                "gf_rung": rung,
                "n_agss21_indicators_in_pool": len(set(hit)),
                "agss21_indicators_in_pool": "|".join(sorted(set(hit))),
                "pool_wavelengths_A": "|".join(f"{w:.3f}" for w in used),
            })

    d = pd.DataFrame(rows).sort_values(["element", "band", "holding", "treatment"])
    d.to_csv(OUT / "cno_products_summary.csv", index=False)

    pd.set_option("display.width", 250)
    print("=== RYA-1214 Step 4 — CNO products ===")
    if d.empty:
        print("  no product artifact found in data/results/band_products/")
        return 0
    print(d[["element", "ion", "band", "holding", "treatment", "A", "sigma_stat",
             "sigma_syst", "n_lines", "n_excluded",
             "n_agss21_indicators_in_pool"]].to_string(index=False))
    print("\n  gf rung per product:")
    for _, r in d.iterrows():
        print(f"    {r.element} {r.ion} {r.band:12s} {r.holding:34s} "
              f"{r.gf_rung[:110]}")
    (OUT / "cno_products_summary.json").write_text(json.dumps({
        "ticket": "RYA-1214", "step": "4 — Solar CNO atomic-indicator products",
        "n_products": int(len(d)),
        "selection_rule": "strongest N by theoretical depth inside the band, min "
                          "separation 4.0 A (select_lines) — OUR selection, not AGSS21's "
                          "adopted set. Line selection is the dominant lever (RYA-842), "
                          "so the indicator overlap is reported per product.",
        "products": d.to_dict("records"),
    }, indent=2, default=str) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
