#!/usr/bin/env python3
"""RYA-728 -- add the Sobeck+2007 Cr I lab lines that neither canonical_gf nor the IR synthesis list
holds (Ryan 2026-11-08: "yes add all lines"), through the shared add route and its clash guard.

The components come from our solar VALD extract, keyed on wavelength and lower level; the canonical
row carries the lab gf and sigma at tier LAB, so the lines join the Reference pool. Cr I 9294.23 sits
0.3 A from a finished Fe I pool line, so the clash guard names it; Ryan's ruling is that a line leaves
only when it cannot be resolved, so it is added, its clash recorded for the run report. Cr II 4627.654
(Gurell+2010) stays out: its list is iSpec's GES list, outside the repo. Run on Sirius.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from pipeline import published_line_add as pla  # noqa: E402

WAVES = (9294.23, 9446.81, 9571.75, 9667.2)
ACCEPT_CLASH = (9294.23,)    # Ryan 2026-11-08
ALREADY_ADDED = (9446.81, 9571.75, 9667.2)
UNRESOLVED = [{"species": "Cr II", "wavelength_A": 4627.654, "ep_eV": 5.6699, "gurell_log_gf": -1.269,
               "theory_log_gf": [-4.682, -4.74],
               "reason": "Gurell 2010 Table 3 sits 3.4 dex above both theory values in the same row. The "
                         "transition is not in the GES or solar VALD lists. The Delbouille atlas shows only "
                         "a shoulder of about 3% on the red wing of Fe I 4627.547. Held under the 0.5 dex rule."}]
AUDIT = ROOT / "data" / "audit" / "orchestrator_prepare" / "solar_Cr_added_lab_lines.json"
CITE = {"sobeck2007": {"doi": "10.1086/519987", "citation": "Sobeck et al. 2007 (ApJ 667, 1267)"}}


def main() -> None:
    lab = pd.read_csv(ROOT / "data" / "reference" / "lab_gf_tables" / "lab_gf_lines.csv", low_memory=False)
    lab = lab[(lab.paper == "sobeck2007") & lab.wavelength_air_A.round(2).isin(WAVES)]
    if len(lab) != len(WAVES):
        raise SystemExit(f"expected {len(WAVES)} Sobeck rows, found {len(lab)}")
    plans = [pla.plan({"species": r.species, "wavelength_A": r.wavelength_air_A, "ep_eV": r.ep_eV,
                       "set": "sobeck2007", "authors_loggf": r.loggf, "authors_sigma_dex": r.loggf_sigma_dex,
                       "authors_gf_source": "Sobeck+2007"}) for r in lab.itertuples()]
    for p in plans:
        print(p["species"], p["wavelength_A"], "extract", p["extract"], "components", p["n_components"],
              "VALD", p["vald_total_loggf"], "lab", p["authors_loggf"], "clashes", p["clashes"][:4])
    if "--write" not in sys.argv:
        return
    todo = [p for p in plans if p["wavelength_A"] not in ALREADY_ADDED]
    written = pla.apply(todo, citation=CITE, tier="LAB", status="adopted_rya728_lab_table_added",
                        element_scoped=True, at_centroid=True, accept_clash=ACCEPT_CLASH)
    done = json.loads(AUDIT.read_text()) if AUDIT.exists() else {"added": []}
    done["added"] += written
    done["accepted_clash"] = [{"wavelength_A": p["wavelength_A"], "clashes": p["clashes"]}
                              for p in plans if p["wavelength_A"] in ACCEPT_CLASH]
    done.pop("refused_clash", None)
    done["unresolved"] = UNRESOLVED
    AUDIT.write_text(json.dumps(done, indent=1, default=str) + "\n")
    print("written", len(written))


if __name__ == "__main__":
    main()
