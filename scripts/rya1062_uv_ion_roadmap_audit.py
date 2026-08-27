#!/usr/bin/env python3
"""Build the read-only RYA-1062 UV ion-roadmap/source reconciliation."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
TRACKER = ROOT / "data/audit/element_status_tracker.csv"
GOLD = ROOT / "data/reference/solar/solar_abundances_v3.csv"
UV = ROOT / "data/audit/rya1061_uv_completeness/per_species_verdict.csv"
OUT = ROOT / "data/audit/rya1062_uv_ion_roadmap"

# Source-level claims only. Counts are publication-wide, never silently treated as the
# number that overlaps Codex. Blank uncertainty/count means transition-table verification
# is still owed.
SOURCES = [
    ("Sc II", "Lawler et al.", 2019, "UV-visible", "PRIMARY_LAB", "BF+LIF", 259, "paper-reported", "yes", "yes", "10.3847/1538-4365/ab08ef", "RYA-1061", "NEW_LAB_SOURCE"),
    ("Ti II", "Wood et al.", 2013, "UV-visible", "PRIMARY_LAB", "BF+LIF", None, "paper-reported", "no", "yes", "", "none found", "NEW_LAB_SOURCE"),
    ("V II", "Wood et al.", 2014, "UV-visible", "PRIMARY_LAB", "BF+LIF", None, "paper-reported", "yes", "yes", "", "RYA-470 seed family", "SOURCE_KNOWN_NOT_INGESTED"),
    ("V II", "Kramida & Saloman", 2017, "UV-visible", "EVALUATED", "critical compilation", None, "accuracy codes", "mixed/trace", "yes", "", "RYA-1061", "NEW_EVALUATED_SOURCE"),
    ("Cr II", "Lawler et al.", 2017, "UV-visible", "PRIMARY_LAB", "BF+LIF", 183, "paper-reported", "no", "verify", "", "none found", "NEW_LAB_SOURCE"),
    ("Cr II", "Ward et al.", 2023, "2080-4140", "PRIMARY_LAB_MIXED_LIFETIME", "BF + experimental/theory lifetimes", 268, "10-26%", "no", "yes", "10.3847/1538-4357/acfafe", "RYA-1061", "SOURCE_KNOWN_NOT_INGESTED"),
    ("Mn II", "Kling & Griesmann", 2000, "vacuum UV", "PRIMARY_LAB", "experimental f-values", None, "paper-reported", "preserve separately", "verify", "", "RYA-1061", "SOURCE_KNOWN_NOT_INGESTED"),
    ("Mn II", "Den Hartog et al.", 2011, "UV-visible", "PRIMARY_LAB", "BF+LIF", None, "paper-reported", "yes", "yes", "", "Mn I cited; ion coverage not reconciled", "SOURCE_KNOWN_NOT_INGESTED"),
    ("Fe II", "Lawler/Den Hartog source family", None, "UV-blue", "UNRESOLVED_PRIMARY_FAMILY", "BF+LIF candidates", None, "unknown", "no", "unknown", "", "RYA-715/909", "SOURCE_KNOWN_NOT_INGESTED"),
    ("Co II", "Mullman, Cooper & Lawler", 1998, "UV resonance", "PRIMARY_LAB", "LIF+UV BF", None, "paper-reported", "no", "verify", "", "none found", "NEW_LAB_SOURCE"),
    ("Co II", "Lawler et al.", 2018, "UV-visible", "PRIMARY_LAB", "echelle BF+LIF", 84, "paper-reported", "yes (28 levels)", "yes", "10.3847/1538-4365/aac773", "RYA-1061", "SOURCE_KNOWN_NOT_INGESTED"),
    ("Ni II", "modern experimental/evaluated search", None, "FUV-NUV", "UNVERIFIED_SEED", "unknown", None, "unknown", "unknown", "unknown", "", "none found", "NEW_LAB_SOURCE"),
    ("Cu II", "UV resonance/subordinate search", None, "FUV-NUV", "UNVERIFIED_SEED", "unknown", None, "unknown", "yes", "unknown", "", "none found", "ROADMAP_ION_MISSING"),
    ("Zn II", "UV resonance/subordinate search", None, "FUV-NUV", "UNVERIFIED_SEED", "unknown", None, "unknown", "isotope check", "unknown", "", "none found", "ROADMAP_ION_MISSING"),
    ("Sr II", "primary/evaluated UV source search", None, "NUV-blue", "UNVERIFIED_SEED", "unknown", None, "unknown", "isotope check", "unknown", "", "RYA-421/551", "SOURCE_KNOWN_NOT_INGESTED"),
    ("Y II", "Hannaford et al.", 1982, "UV-visible", "PRIMARY_LAB", "BF+lifetimes", None, "paper-reported", "no", "verify", "", "RYA-475 seed", "SOURCE_KNOWN_NOT_INGESTED"),
    ("Zr II", "modern laboratory source search", None, "NUV-blue", "UNVERIFIED_SEED", "unknown", None, "unknown", "isotope check", "unknown", "", "RYA-560", "NEW_LAB_SOURCE"),
    ("Ba II", "laboratory/evaluated source family", None, "NUV-blue", "UNVERIFIED_SEED", "unknown", None, "unknown", "yes", "unknown", "", "RYA-559/581", "SOURCE_KNOWN_NOT_INGESTED"),
    ("Ce II", "Lawler et al.", 2009, "UV-visible", "PRIMARY_LAB", "BF+LIF", None, "paper-reported", "no", "yes", "", "none found", "NEW_LAB_SOURCE"),
    ("Nd II", "Den Hartog et al.", 2003, "UV-visible", "PRIMARY_LAB", "BF+LIF", None, "paper-reported", "no", "yes", "", "none found", "NEW_LAB_SOURCE"),
    ("Eu II", "Lawler et al.", 2001, "UV-visible", "PRIMARY_LAB", "BF+LIF", None, "paper-reported", "yes + isotope", "yes", "", "RYA-565 source family needs check", "SOURCE_KNOWN_NOT_INGESTED"),
    ("P II", "UV laboratory/evaluated search", None, "FUV-NUV", "UNVERIFIED_SEED", "unknown", None, "unknown", "no", "unknown", "", "none found", "ROADMAP_ION_MISSING"),
    ("S II", "UV laboratory/evaluated search", None, "FUV-NUV", "UNVERIFIED_SEED", "unknown", None, "unknown", "no", "unknown", "", "none found", "ROADMAP_ION_MISSING"),
]

ROLE = {"Fe II": "arbiter", "V II": "future ionization anchor"}
EXPLICIT_TICKETS = {
    "Fe II": "RYA-715;RYA-909", "Sc II": "RYA-460;RYA-683", "V II": "RYA-470;RYA-733",
    "Cr II": "RYA-715-family;RYA-1061", "Sr II": "RYA-421;RYA-551;RYA-683",
    "Y II": "RYA-458;RYA-475;RYA-683", "Zr II": "RYA-279;RYA-560;RYA-683",
    "Ba II": "RYA-559;RYA-581;RYA-683", "Eu II": "RYA-565;RYA-683",
}
MISSING = [
    ("P", "II", "strong FUV/NUV ion diagnostics plausible", "ROADMAP_ION_MISSING", "Literature/source census, then assess STIS/COS photospheric reach."),
    ("S", "II", "strong FUV/NUV ion diagnostics plausible", "ROADMAP_ION_MISSING", "Separate photospheric candidates from ISM and transition-region lines."),
    ("Cu", "II", "neutral VIS sources do not cover the ion", "ROADMAP_ION_MISSING", "Search resonance/subordinate laboratory data; preserve HFS."),
    ("Zn", "II", "neutral VIS sources do not cover the ion", "ROADMAP_ION_MISSING", "Search permitted photospheric lines; exclude forbidden plasma diagnostics."),
    ("Ni", "II", "modern UV source family unresolved", "ROADMAP_ION_MISSING", "Separate permitted abundance lines from forbidden low-density-plasma lines."),
]


def read_tracker() -> pd.DataFrame:
    return pd.read_csv(TRACKER, comment="#")


def build(out: Path) -> dict:
    out.mkdir(parents=True, exist_ok=True)
    tracker = read_tracker()
    gold = pd.read_csv(GOLD, comment="#")
    uv = pd.read_csv(UV).set_index("species")
    rows = []
    for _, r in tracker.iterrows():
        species = f"{r.element} {r.ion.split('+')[0]}"
        if r.element == "Fe" and r.ion == "I+II":
            species = "Fe I"
        current = uv.loc[species] if species in uv.index else None
        role = ROLE.get(species, "headline" if r.verdict in ("PASS", "CURATION-OWED") else "future")
        tickets = EXPLICIT_TICKETS.get(species, str(r.source_tickets))
        rows.append((r.element, r.ion, tickets, role, str(r.engine_reach),
                     int(current.canonical_uv_rows) if current is not None else 0,
                     int(current.lab_rows) if current is not None else 0,
                     "yes" if species in set(s[0] for s in SOURCES) else "no",
                     str(r.classification), str(r.action_needed)))
    inventory = pd.DataFrame(rows, columns=["element", "ion", "existing_tickets", "current_role", "wavelength_regions_considered", "canonical_uv_rows", "canonical_primary_lab_rows", "uv_source_considered", "status", "current_action"])
    inventory.to_csv(out / "ion_ticket_inventory.csv", index=False)

    cols = ["ion", "source", "year", "wavelength_coverage", "evidence_type", "method", "published_transition_count", "uncertainty", "hfs_isotope_support", "machine_readable", "doi", "already_cited_in_ion_ticket", "gap_class"]
    rec = pd.DataFrame(SOURCES, columns=cols)
    rec["already_represented_in_canonical_gf"] = rec.ion.map(lambda x: "partial" if x in uv.index and int(uv.loc[x].lab_rows) else "no")
    rec["uv_overlap_with_holdings"] = rec.wavelength_coverage.map(lambda x: "potential: STIS/COS or near-UV" if x else "unresolved")
    rec["action"] = rec.gap_class.map({
        "SOURCE_KNOWN_NOT_INGESTED": "Acquire table; reconcile physical transitions; propose separate ingest ticket.",
        "NEW_LAB_SOURCE": "Verify bibliography/table and route to RYA-1061/946 plus ion ticket.",
        "NEW_EVALUATED_SOURCE": "Trace every adopted value to underlying primary work and retain accuracy code.",
        "ROADMAP_ION_MISSING": "Open roadmap discovery ticket only after photospheric diagnostic/reach validation.",
    })
    rec.to_csv(out / "uv_source_reconciliation.csv", index=False)

    pd.DataFrame(MISSING, columns=["element", "ion", "evidence", "gap_class", "action"]).to_csv(out / "roadmap_missing_ions.csv", index=False)
    absent = rec[rec.gap_class.isin(["NEW_LAB_SOURCE", "NEW_EVALUATED_SOURCE"])].copy()
    absent.to_csv(out / "modern_sources_absent_from_ion_tickets.csv", index=False)

    mismatches = []
    g = gold.set_index("element")
    for _, r in tracker.iterrows():
        # Fe II is an intentionally separate arbiter row beside the Fe I headline. It is
        # not the species label of the frozen headline abundance and therefore cannot be
        # compared as though it replaced that row.
        if r.element not in g.index or r.ion == "I+II" or (r.element == "Fe" and r.ion == "II"):
            continue
        gold_ion = str(g.loc[r.element].ion)
        if gold_ion != str(r.ion):
            mismatches.append((r.element, gold_ion, r.ion, "ION_LABEL_MISMATCH", "RYA-683",
                               "Record only; frozen gold remains immutable until the next freeze."))
    pd.DataFrame(mismatches, columns=["element", "gold_ion", "measured_ion", "gap_class", "route_to", "action"]).to_csv(out / "ion_label_mismatches.csv", index=False)

    verdict_rows = []
    for ion in sorted(set(rec.ion)):
        rr = rec[rec.ion == ion]
        if (rr.gap_class == "ROADMAP_ION_MISSING").any(): verdict = "material source gap"
        elif (rr.evidence_type == "UNRESOLVED_PRIMARY_FAMILY").any(): verdict = "material source gap"
        elif ion in uv.index and int(uv.loc[ion].lab_rows) > 0: verdict = "likely complete"
        elif rr.evidence_type.str.startswith("UNVERIFIED").any(): verdict = "no current UV path"
        else: verdict = "material source gap"
        verdict_rows.append((ion, verdict, ";".join(sorted(set(rr.gap_class))), "RYA-1061;RYA-946;ion-specific ticket"))
    pd.DataFrame(verdict_rows, columns=["ion", "verdict", "basis", "route_findings_to"]).to_csv(out / "per_ion_verdict.csv", index=False)

    follow = [
        (1, "Cr II Ward 2023 transition-table reconciliation", "RYA-1061;RYA-946;Cr ion ticket", "current near-UV plus future NUV", "largest verified modern UV table; mixed lifetime provenance must remain per row"),
        (2, "Sc II/Ti II/V II/Co II machine-readable table acquisition", "RYA-1061;RYA-946;ion tickets", "near-UV/STIS", "majority-ion and HFS-sensitive anchors"),
        (3, "Resolve exact Fe II primary UV papers", "RYA-715;RYA-909;RYA-1061", "solar/proxy UV path", "do not let a NIST compilation stand in for the experiment"),
        (4, "P II/S II/Ni II/Cu II/Zn II diagnostic feasibility", "new follow-ups;RYA-683 if adopted", "FUV/NUV STIS/COS", "roadmap stages missing; formation/ISM/plasma-lane screen first"),
        (5, "Heavy-ion UV branch census", "RYA-946;per-element 700-series", "near-UV/STIS", "Ce II has 2082 canonical UV rows but zero LAB-tier rows"),
    ]
    pd.DataFrame(follow, columns=["priority", "followup", "route_to", "data_path", "reason"]).to_csv(out / "followup_recommendations.csv", index=False)

    summary = {"ticket": "RYA-1062", "roadmap_rows": len(inventory), "source_pairs": len(rec),
               "missing_ion_stages": len(MISSING), "ion_label_mismatches": len(mismatches),
               "new_source_leads": len(absent), "mutation_policy": "READ_ONLY_NO_STATE_CHANGES"}
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    (out / "README.md").write_text(f"""# RYA-1062 — UV ion-roadmap cross-scan

Read-only reconciliation of the committed ion tracker, frozen solar gold, RYA-1061 UV
inventory, and source-level laboratory leads. No canonical gf, line pool, ion label,
product, abundance, or roadmap state is changed.

## Verdict

The live tracker has {len(inventory)} ion/element rows. Frozen gold disagrees with the
measured ion for **{len(mismatches)}** elements; these are routed to RYA-683 and remain
unmodified here. Five scientifically plausible ion stages (P II, S II, Ni II, Cu II,
Zn II) are absent from the roadmap and require a photospheric-line/data-reach screen
before tickets can promote them.

No surveyed ion is declared `complete`. The modern primary-lab tables are generally not
represented as LAB-tier canonical UV rows. Cr II is the clearest material opportunity:
Ward et al. (2023) publishes 268 transitions over 208–414 nm, while Codex has 826 Cr II
near-UV rows and zero primary-lab rows. Sc II, Ti II, V II, Mn II, Co II, and the
rare-earth majority ions also retain source-level reconciliation work.

`UNVERIFIED_SEED` and `UNRESOLVED_PRIMARY_FAMILY` are deliberately not evidence grades.
They are stop signs: exact bibliography, tables, wavelength frame, physical transition,
uncertainty, and HFS/isotope provenance must be resolved before adoption.

Reproduce with `python3 scripts/rya1062_uv_ion_roadmap_audit.py`.
""")
    return summary


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=OUT)
    args = ap.parse_args()
    print(json.dumps(build(args.out), indent=2))


if __name__ == "__main__":
    main()
