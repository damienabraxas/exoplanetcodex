#!/usr/bin/env python3
"""RYA-728 -> RYA-1378: list every Cr I / Cr II line without a primary laboratory gf.

Each row is a candidate for the future theoretical grade (RYA-1351 THEORETICAL_CODEX)
or for a laboratory measurement. Nothing here changes a gf, a grade or a product.

Flags use RYA-1378's vocabulary and may coexist:
  UNRATED_GF        gf present, no primary laboratory source attached
  MISSING_GF        a laboratory line with no row in canonical_gf at all
  PROVENANCE_GAP    gf taken from a compilation (MFW, CSE, PhFl) or an opaque VALD3 tag:
                    the value may trace to a lab paper that has not been attached
  GF_DISCREPANCY    a lab gf and the list gf disagree by more than 0.5 dex (RYA-728 held lines)
  SOLAR_DISCREPANCY a Codex fit lands more than 0.5 dex from the solar Cr abundance, or is
                    non-physical (from the RYA-728 runs and band diagnostics)

Priority (higher = more useful to a theoretical or laboratory effort):
  1  the line's band has NO laboratory route at all (every Cr cell there is HELD)
  2  a GF_DISCREPANCY or SOLAR_DISCREPANCY line in a measured band
  3  an ungraded line in the synthesis list in a band that is measured
  4  not in the synthesis list (not fitted by any engine today)
Observability (depth, blending) is NOT established here and is not inferred.

Usage: python scripts/rya728_cr_gap_candidates.py [--out PATH]
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
CANON = ROOT / "data/linelists/canonical_gf.csv"
DIAG = ROOT / "data/audit/orchestrator_prepare/solar_Cr_band_diagnostics.csv"
OUT = ROOT / "data/audit/orchestrator_prepare/solar_Cr_atomic_gap_candidates.csv"

BANDS = [("near-UV", 3000, 3780), ("VIS", 3780, 6910), ("red-optical", 6910, 9200),
         ("NIR", 9200, 13000), ("H", 15000, 17500), ("K", 19400, 25000)]
# (species, band) with no laboratory gf source at all: every Cr cell there is HELD (RYA-728
# literature sweep, data/reference/litscan/Cr.yaml).
NO_LAB_ROUTE = {("Cr I", "H"), ("Cr I", "K"), ("Cr II", "red-optical"), ("Cr II", "NIR"),
                ("Cr II", "H"), ("Cr II", "K")}
COMPILATIONS = {"MFW", "CSE", "1968PhFl.."}
SOLAR_CR = 5.62  # AGSS21 (litscan Cr.yaml); only used to flag a >0.5 dex fit
# RYA-728 held / conflicting lines: (species, wavelength, evidence)
KNOWN = [
    ("Cr II", 3665.979, "GF_DISCREPANCY", "Gurell vs list >0.5 dex; no fit row returned"),
    ("Cr II", 4627.654, "GF_DISCREPANCY", "Gurell gf 3.4 dex above both theories; unresolved"),
    ("Cr II", 4891.4948, "GF_DISCREPANCY;SOLAR_DISCREPANCY",
     "lab vs list >0.5 dex; fit non-physical with either gf (KP 7.32-7.33, HARPS 7.42)"),
    ("Cr II", 6421.7063, "GF_DISCREPANCY;SOLAR_DISCREPANCY",
     "list gf non-physical (7.99 / 8.11); Lawler gf gives 2.64 KP, 4.40 HARPS"),
    ("Cr II", 3179.449, "GF_DISCREPANCY", "Gurell -2.119 vs Ward -1.263 (0.86 dex); held"),
    ("Cr II", 3183.33, "GF_DISCREPANCY", "Gurell -0.965 vs Ward -0.214 (0.75 dex); held"),
    ("Cr I", 4026.21, "GF_DISCREPANCY", "Sobeck -1.92 vs list; not adopted, not yet tested"),
]
# Lab values found but not yet attached (none after the 11-09 Tozzi / Ward adoption)
LAB_PENDING: dict = {}
MATCH_A = 0.02


def band_of(w: float) -> str:
    for name, lo, hi in BANDS:
        if lo <= w < hi:
            return name
    return "outside"


def build() -> pd.DataFrame:
    c = pd.read_csv(CANON, low_memory=False)
    c = c[c.species.isin(["Cr I", "Cr II"])].copy()
    c["band"] = c.wavelength_air_A.map(band_of)
    lab = c.lab_source_tag.notna() & (c.gf_tier == "LAB")
    every = c
    c = c[~lab].copy()

    flags = {i: ["UNRATED_GF"] for i in c.index}
    evidence = {i: [] for i in c.index}
    ref = c.loggf_reference.astype(str)
    for i in c.index[ref.isin(COMPILATIONS) | (ref == "VALD3")]:
        flags[i].append("PROVENANCE_GAP")
        evidence[i].append(f"gf via {ref[i]}: source paper not attached")

    def near(sp, w):
        m = (c.species == sp) & ((c.wavelength_air_A - w).abs() <= MATCH_A)
        return list(c.index[m])

    absent = []
    for sp, w, fl, ev in KNOWN:
        hits = near(sp, w)
        for i in hits:
            flags[i] += [f for f in fl.split(";") if f not in flags[i]]
            evidence[i].append(ev)
        in_table = ((every.species == sp) & ((every.wavelength_air_A - w).abs() <= MATCH_A)).any()
        if not in_table:  # a lab line with no row in canonical_gf: still a candidate
            absent.append({"species": sp, "wavelength_air_A": w, "band": band_of(w),
                           "in_synth": False, "flags": "MISSING_GF;" + fl,
                           "evidence": ev + " | no row in canonical_gf"})
    for (sp, w), ev in LAB_PENDING.items():
        for i in near(sp, w):
            evidence[i].append("lab gf found, not attached: " + ev)

    d = pd.read_csv(DIAG)
    for r in d.itertuples():
        sp = f"Cr {r.ion}"
        bad = (not bool(r.in_aggregate)) or abs(float(r.A_1d_lte) - SOLAR_CR) > 0.5
        if not bad:
            continue
        for i in near(sp, float(r.wavelength_air_A)):
            if "SOLAR_DISCREPANCY" not in flags[i]:
                flags[i].append("SOLAR_DISCREPANCY")
            evidence[i].append(f"ungraded diagnostic {r.holding}: A={r.A_1d_lte}"
                               + ("" if bool(r.in_aggregate) else " (excluded)"))

    c["flags"] = [";".join(flags[i]) for i in c.index]
    c["evidence"] = [" | ".join(evidence[i]) for i in c.index]
    if absent:
        c = pd.concat([c, pd.DataFrame(absent)], ignore_index=True)
    c["no_lab_route_in_band"] = [(s, b) in NO_LAB_ROUTE for s, b in zip(c.species, c.band)]
    disc = c["flags"].str.contains("DISCREPANCY")
    c["priority"] = 3
    c.loc[disc, "priority"] = 2
    c.loc[~c.in_synth.astype(bool) & ~disc, "priority"] = 4
    # IR lines are fitted from their own band lists, not the optical synthesis list, so
    # in_synth does not measure reach there: every line in a no-lab-route band is priority 1.
    c.loc[c.no_lab_route_in_band, "priority"] = 1
    c["candidate_for"] = "theoretical grade (RYA-1351) or laboratory measurement"
    c["review_state"] = "candidate"
    cols = ["line_id", "physical_id", "species", "wavelength_air_A", "excitation_potential_eV",
            "log_gf", "loggf_reference", "gf_tier", "band", "in_synth",
            "no_lab_route_in_band", "flags", "priority", "evidence", "candidate_for",
            "review_state"]
    return c[cols].sort_values(["priority", "species", "wavelength_air_A"])


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=Path, default=OUT)
    a = ap.parse_args()
    t = build()
    t.to_csv(a.out, index=False)
    print(f"{len(t)} Cr rows without a primary laboratory gf -> {a.out}")
    print(t.groupby(["priority", "species", "band"]).size().to_string())
    for f in ["PROVENANCE_GAP", "GF_DISCREPANCY", "SOLAR_DISCREPANCY"]:
        print(f, int(t["flags"].str.contains(f).sum()))


if __name__ == "__main__":
    main()
