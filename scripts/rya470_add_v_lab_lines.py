#!/usr/bin/env python3
"""RYA-470 -- add the V I lab lines that canonical_gf does not hold, through the shared add route
(pipeline.published_line_add) and its clash guard, as RYA-728 did for Cr.

  Lawler+2014  7944.6472  (red-optical)
  Holmes+2016  19998.913, 20230.568  (K band; the only lab V I gf there)

Lawler 5845.2695 stays out: no solar hfs-on VALD extract covers 3780-6910 A, so the route has
no components for it (atomic_gap_known.csv records it). Run on Sirius, where the extracts live.
Without --write it only prints the plan.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from pipeline import published_line_add as pla  # noqa: E402

WAVES = {"lawler2014_v1": (7944.6472,), "holmes2016_v1": (19998.913, 20230.568)}
AUDIT = ROOT / "data" / "audit" / "orchestrator_prepare" / "solar_V_added_lab_lines.json"
CITE = {"lawler2014_v1": {"doi": "10.1088/0067-0049/215/2/20", "citation": "Lawler et al. 2014 (ApJS 215, 20)"},
        "holmes2016_v1": {"doi": "10.3847/0067-0049/224/2/35", "citation": "Holmes et al. 2016 (ApJS 224, 35)"}}


def main() -> None:
    lab = pd.read_csv(ROOT / "data" / "reference" / "lab_gf_tables" / "lab_gf_lines.csv", low_memory=False)
    plans = []
    for paper, waves in WAVES.items():
        for w in waves:
            r = lab[(lab.paper == paper) & ((lab.wavelength_air_A - w).abs() < 1e-4)]
            if len(r) != 1:
                raise SystemExit(f"{paper} {w}: {len(r)} lab rows")
            r = r.iloc[0]
            plans.append(pla.plan({"species": r.species, "wavelength_A": float(r.wavelength_air_A),
                                   "ep_eV": float(r.ep_eV), "set": paper, "authors_loggf": float(r.loggf),
                                   "authors_sigma_dex": float(r.loggf_sigma_dex), "authors_gf_source": paper}))
    for p in plans:
        print(p["species"], p["wavelength_A"], "extract", p["extract"], "components", p["n_components"],
              "VALD", p["vald_total_loggf"], "lab", p["authors_loggf"], "clashes", p["clashes"][:4])
    if "--write" not in sys.argv:
        return
    clean = [p for p in plans if not p["clashes"]]
    for p in plans:
        if p["clashes"]:
            print("NOT ADDED (clash guard; needs a ruling):", p["wavelength_A"], p["clashes"])
    written = pla.apply(clean, citation=CITE, tier="LAB", status="adopted_rya470_lab_table_added",
                        element_scoped=True, at_centroid=True)
    AUDIT.write_text(json.dumps({"added": written,
                                 "refused_clash": [{"wavelength_A": p["wavelength_A"], "clashes": p["clashes"]}
                                                   for p in plans if p["clashes"]]},
                                indent=1, default=str) + "\n")
    print("written", len(written))


if __name__ == "__main__":
    main()
