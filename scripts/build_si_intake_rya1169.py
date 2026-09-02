#!/usr/bin/env python3
"""Build the evidence-first RYA-1169 silicon intake artifacts.

This builder deliberately separates reference membership, laboratory gf
provenance, and Solar-depth routing.  It does not promote a line to an active
grade until an approved Solar spectrum supplies a measured feature depth.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "audit" / "rya1169_si_intake"
CANONICAL = ROOT / "data" / "linelists" / "canonical_gf.csv"

REFERENCE = [
    # Amarsi & Asplund (2017), table 1. The first two are explicitly excluded
    # from the final Solar abundance; the remaining Si I rows enter the mean.
    ("Si I", 3905.53, 1.9087, -0.99, "explicitly_rejected", "violet test line; excluded from final Solar mean"),
    ("Si I", 4102.93, 1.9087, -3.04, "explicitly_rejected", "violet test line; excluded from final Solar mean"),
    ("Si I", 5645.61, 4.9296, -2.04, "used", ""),
    ("Si I", 5684.48, 4.9538, -1.55, "used", ""),
    ("Si I", 5690.43, 4.9296, -1.77, "used", ""),
    ("Si I", 5701.11, 4.9296, -1.95, "used", ""),
    ("Si I", 5772.15, 5.0823, -1.65, "used", ""),
    ("Si I", 5793.07, 4.9296, -1.96, "used", ""),
    ("Si I", 6741.64, 5.9840, -1.65, "used", ""),
    ("Si I", 7034.90, 5.8708, -0.78, "used", ""),
    ("Si I", 7226.21, 5.6135, -1.41, "used", ""),
    # Scott et al. (2015), table 2; retained in the weighted Solar mean.
    ("Si II", 6371.370, 8.121, -0.044, "used", "3D non-LTE correction assumed zero"),
]

# Den Hartog et al. (2023), tables 3 and 4. Percent A-value uncertainty is
# converted to dex by sigma_loggf = log10(1 + p/100).
DH23_SI1 = [
    (2207.978, 0.0000, -1.248, 5), (2210.892, 0.00956, -0.895, 5),
    (2211.745, 0.00956, -1.388, 5), (2218.057, 0.02767, -1.402, 5),
    (2218.916, 0.02767, -2.603, 6), (2438.768, 0.0000, -2.723, 9),
    (2443.365, 0.00956, -2.828, 9), (2452.118, 0.02767, -2.868, 9),
    (2506.897, 0.00956, -0.595, 5), (2514.316, 0.0000, -0.672, 5),
    (2516.112, 0.02767, -0.098, 5), (2519.202, 0.00956, -0.810, 5),
    (2528.508, 0.02767, -0.583, 5), (2563.679, 0.78085, -3.922, 20),
    (2564.825, 0.78085, -4.228, 26), (2881.578, 0.78085, -0.088, 5),
    (2970.353, 0.78085, -3.531, 11), (2987.643, 0.78085, -2.035, 8),
    (3905.523, 1.9087, -1.077, 10), (4102.936, 1.9087, -3.026, 18),
]
DH23_SI2 = [(2334.407, 0.03634, -5.088, 16), (2350.172, 0.03634, -5.116, 16)]

REF_DOI = "10.1093/mnras/stw2445"
SCOTT_DOI = "10.1051/0004-6361/201424109"
DH23_DOI = "10.3847/1538-4365/acb642"


def band(w: float) -> str:
    if w < 3000: return "FUV"
    if w < 3800: return "NUV"
    if w < 6910: return "VIS"
    if w < 10000: return "red-optical"
    return "IR"


def load_canonical() -> list[dict[str, str]]:
    with CANONICAL.open(newline="") as handle:
        return [r for r in csv.DictReader(handle) if r["species"] in {"Si I", "Si II"}]


def match(rows: list[dict[str, str]], species: str, wave: float, ep: float) -> tuple[dict[str, str] | None, str]:
    # Physical-transition join: species + air wavelength + EP. Tolerances only
    # absorb published rounding and never resolve multiple candidates silently.
    hits = [r for r in rows if r["species"] == species
            and abs(float(r["wavelength_air_A"]) - wave) <= 0.03
            and abs(float(r["excitation_potential_eV"]) - ep) <= 0.01]
    return (hits[0], "matched") if len(hits) == 1 else (None, "missing" if not hits else "ambiguous")


def write_csv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("w", newline="") as handle:
        out = csv.DictWriter(handle, fieldnames=fields)
        out.writeheader(); out.writerows(rows)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    canonical = load_canonical()
    reference_rows = []
    for species, wave, ep, loggf, status, note in REFERENCE:
        c, join = match(canonical, species, wave, ep)
        reference_rows.append({
            "reference_line_set": "AGSS21_Si_Amarsi2017_Scott2015",
            "adopted_solar_value_lineage": "AGSS21 7.51 -> Amarsi2017 7.51 -> Scott2015 lines",
            "source_paper_table": "Amarsi2017 Table 1" if species == "Si I" else "Scott2015 Table 2",
            "species": species, "wavelength_air_A": f"{wave:.3f}", "ep_eV": f"{ep:.4f}",
            "published_loggf": f"{loggf:.3f}", "published_gf_source": "Garz1973+0.097dex (O'Brian&Lawler lifetimes)",
            "reference_status": status, "band": band(wave), "canonical_line_id": c["line_id"] if c else "",
            "canonical_loggf": c["log_gf"] if c else "", "canonical_gf_tier": c["gf_tier"] if c else "",
            "delta_reference_vs_codex": f"{loggf-float(c['log_gf']):.3f}" if c else "",
            "join_status": join, "notes": note,
        })
    ref_fields = list(reference_rows[0])
    write_csv(OUT / "si_agss21_reference_lines.csv", reference_rows, ref_fields)

    lab_rows = []
    for species, source in (("Si I", DH23_SI1), ("Si II", DH23_SI2)):
        for wave, ep, loggf, pct in source:
            c, join = match(canonical, species, wave, ep)
            lab_rows.append({
                "species": species, "wavelength_air_A": f"{wave:.3f}", "ep_eV": f"{ep:.5f}",
                "experimental_loggf": f"{loggf:.3f}", "uncertainty_percent": pct,
                "lab_sigma_dex": f"{math.log10(1+pct/100):.5f}",
                "primary_lab_source": "Den Hartog et al. 2023", "doi": DH23_DOI,
                "method": "branching fractions + published radiative lifetimes",
                "canonical_line_id": c["line_id"] if c else "", "canonical_loggf": c["log_gf"] if c else "",
                "delta_lab_vs_codex": f"{loggf-float(c['log_gf']):.3f}" if c else "", "join_status": join,
                "band": band(wave), "solar_feature_depth": "", "depth_status": "NOT_MEASURED",
                "grade": "ungraded", "grade_reason": "approved Solar feature depth not yet measured",
            })
    lab_fields = list(lab_rows[0])
    write_csv(OUT / "si_primary_lab_gf_census.csv", lab_rows, lab_fields)

    joined = []
    for r in lab_rows:
        refs = [x for x in reference_rows if x["species"] == r["species"] and abs(float(x["wavelength_air_A"])-float(r["wavelength_air_A"])) <= .03]
        joined.append({**r, "AGSS21_reference_used": "true" if refs and refs[0]["reference_status"] == "used" else "false",
                       "reference_status": refs[0]["reference_status"] if refs else "not_member",
                       "reference_loggf": refs[0]["published_loggf"] if refs else ""})
    write_csv(OUT / "si_joined_line_product.csv", joined, list(joined[0]))

    coverage = []
    for b in ("FUV", "NUV", "VIS", "red-optical", "IR"):
        rr = [x for x in reference_rows if x["band"] == b]; ll = [x for x in lab_rows if x["band"] == b]
        coverage.append({"band": b, "reference_used": sum(x["reference_status"] == "used" for x in rr),
                         "reference_rejected": sum(x["reference_status"] != "used" for x in rr),
                         "codex_matched": sum(x["join_status"] == "matched" for x in ll),
                         "primary_lab_available": len(ll), "codex_grade": 0, "deep_grade": 0,
                         "below_floor": 0, "ungraded_depth_pending": len(ll),
                         "missing_from_codex": sum(x["join_status"] == "missing" for x in ll),
                         "ambiguous_crossmatches": sum(x["join_status"] == "ambiguous" for x in ll)})
    write_csv(OUT / "coverage_matrix.csv", coverage, list(coverage[0]))

    sources = [
        {"role":"AGSS21 adopted abundance", "citation":"Asplund, Amarsi & Grevesse 2021", "doi":"10.1051/0004-6361/202140445", "article_url":"https://doi.org/10.1051/0004-6361/202140445", "download_url":"https://arxiv.org/abs/2105.01661", "used":"yes"},
        {"role":"3D non-LTE Si abundance and final line membership", "citation":"Amarsi & Asplund 2017", "doi":REF_DOI, "article_url":"https://doi.org/"+REF_DOI, "download_url":"https://arxiv.org/abs/1609.07283", "used":"yes"},
        {"role":"underlying Solar Si line table", "citation":"Scott et al. 2015", "doi":SCOTT_DOI, "article_url":"https://doi.org/"+SCOTT_DOI, "download_url":"https://arxiv.org/abs/1405.0279", "used":"yes"},
        {"role":"primary laboratory UV/blue gf", "citation":"Den Hartog et al. 2023", "doi":DH23_DOI, "article_url":"https://doi.org/"+DH23_DOI, "download_url":"https://arxiv.org/abs/2301.11391", "used":"yes"},
        {"role":"Scott/Amarsi relative gf scale", "citation":"Garz 1973, A&A 26, 471", "doi":"", "article_url":"https://ui.adsabs.harvard.edu/abs/1973A%26A....26..471G", "download_url":"", "used":"lineage; row values obtained through Scott/Amarsi"},
        {"role":"absolute lifetime normalization", "citation":"O'Brian & Lawler 1991a", "doi":"10.1016/0375-9601(91)90990-E", "article_url":"https://doi.org/10.1016/0375-9601(91)90990-E", "download_url":"", "used":"lineage"},
        {"role":"absolute lifetime normalization", "citation":"O'Brian & Lawler 1991b", "doi":"10.1103/PhysRevA.44.7134", "article_url":"https://doi.org/10.1103/PhysRevA.44.7134", "download_url":"", "used":"lineage"},
    ]
    write_csv(OUT / "literature_sources.csv", sources, list(sources[0]))
    summary = {"reference_rows": len(reference_rows), "reference_used": sum(x["reference_status"] == "used" for x in reference_rows),
               "primary_lab_rows": len(lab_rows), "active_graded_rows": 0,
               "freeze_status": "BLOCKED", "blockers": ["approved Solar-depth measurement absent", "historical primary-lab census beyond DH23/Garz lineage not yet exhaustive", "NIST row-level cross-check absent"]}
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__": main()
