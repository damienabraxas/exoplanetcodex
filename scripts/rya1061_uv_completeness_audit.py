#!/usr/bin/env python3
"""Build the read-only RYA-1061 UV completeness audit (1150 <= lambda < 3780 A)."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
CANONICAL = ROOT / "data/linelists/canonical_gf.csv"
HOLDINGS = ROOT / "data/catalog/holdings_manifest_registry.csv"
DEFAULT_OUT = ROOT / "data/audit/rya1061_uv_completeness"

REGIMES = (("FUV", 1150.0, 2000.0), ("NUV", 2000.0, 3000.0), ("NEAR_UV", 3000.0, 3780.0))

LAB_SOURCES = [
    ("Cr II", "Ward et al. 2023", "2080-3780", 268, "10-26%", "10.3847/1538-4357/acfafe", "PRIMARY_LAB", "ACQUIRE_TABLE", "Highest-volume non-Fe candidate; published table extends to 4140 A, audit window clips at 3780 A."),
    ("Mn II", "Kling & Griesmann 2000", "UV", None, "paper-reported", "", "PRIMARY_LAB", "VERIFY_AND_ACQUIRE", "Seed named by ticket; transition-level table must be acquired and wavelength-frame checked."),
    ("Co II", "Lawler et al. 2018", "UV/visible", 84, "paper-reported", "10.3847/1538-4365/aac773", "PRIMARY_LAB", "ACQUIRE_TABLE", "82 absolute transition probabilities from branching fractions plus lifetimes; HFS constants for 28 levels."),
    ("Sc I; Sc II", "Lawler et al. 2019", "UV/visible", None, "paper-reported", "", "PRIMARY_LAB", "VERIFY_AND_ACQUIRE", "Retain ion and HFS metadata; do not collapse Sc I and Sc II."),
    ("V II", "Kramida & Saloman 2017", "UV/visible", None, "evaluated", "", "EVALUATED", "VERIFY_AND_ACQUIRE", "Evaluated compilation is a separate evidence tier, not primary laboratory gf."),
    ("Fe I", "Belmonte et al. 2017", "<3780", 70, "0.02-0.10 dex in Codex table", "", "PRIMARY_LAB", "ALREADY_PARTIAL", "Known RYA-836 source; 26 lines lie in 3000-3780 in the source census."),
    ("Fe I", "Den Hartog et al. 2014", "<3780", 31, "0.02-0.12 dex in Codex table", "", "PRIMARY_LAB", "ALREADY_PRESENT", "Known RYA-836 source; all 31 lines lie in 3000-3780."),
    ("Fe I", "Ruffoni et al. 2014", "<3780", 4, "0.02-0.08 dex in Codex table", "", "PRIMARY_LAB", "ALREADY_PRESENT", "Known RYA-836 source; all four lines lie in 3000-3780."),
    ("Ti I; Ti II", "laboratory/evaluated source family", "1150-3780", None, "unknown", "", "UNVERIFIED_SEED", "LITERATURE_SEARCH", "Large canonical VALD-only population makes this a high-priority source search."),
    ("Ni I; Ni II", "laboratory/evaluated source family", "1150-3780", None, "unknown", "", "UNVERIFIED_SEED", "LITERATURE_SEARCH", "Prioritize transitions overlapping actual HST holdings."),
    ("Cu I; Cu II; Zn I; Zn II", "laboratory/evaluated source family", "1150-3780", None, "unknown", "", "UNVERIFIED_SEED", "LITERATURE_SEARCH", "Preserve HFS/isotope structure where available."),
    ("P I; P II; S I; S II", "laboratory/evaluated source family", "1150-3000", None, "unknown", "", "UNVERIFIED_SEED", "LITERATURE_SEARCH", "FUV/NUV species cannot overlap the current canonical list because it has no rows below 3000 A."),
]

PAPERS = [
    ("ASTRAL spectral library", "Ayres et al.; MAST HLSP", "Procyon; Alpha Cen A", "1150-3100", "STIS E140H/E140M/E230H/E230M", "PHOTOSPHERIC_LINE_IDENTIFICATION_OR_VALIDATION", "SUPPORTING", "High-resolution calibrated spectra; archive cautions and artifacts remain applicable.", "https://archive.stsci.edu/prepds/astral/"),
    ("Alpha Cen A ultraviolet spectrum", "Pagano et al. 2004", "Alpha Cen A", "1133-3150", "STIS E140H/E230H", "CHROMOSPHERE_TRANSITION_REGION_CORONA", "SUPPORTING", "Atmospheric-structure and emission-line analysis, not a photospheric abundance authority.", "https://arxiv.org/abs/astro-ph/0310901"),
    ("Alpha Cen A/B Be II 3130", "HST/GHRS abundance literature", "Alpha Cen A; Alpha Cen B", "near 3130", "GHRS", "PHOTOSPHERIC_ABUNDANCE", "CANDIDATE_ABUNDANCE", "Verify exact bibliographic record, component attribution, adopted gf, and blend treatment before extraction.", ""),
    ("Procyon GHRS ultraviolet atlas", "HST atlas product", "Procyon", "UV", "GHRS", "CHROMOSPHERE_TRANSITION_REGION_CORONA", "SUPPORTING", "Useful line-identification context; not an abundance publication.", "https://archive.stsci.edu/hst/atlasprocyon/"),
]

GAPS = [
    (1, "FUV canonical vacuum", "all", "1150-2000", "0 canonical rows", "Acquire primary/evaluated transition tables, then overlap with STIS/COS wavelength support."),
    (2, "NUV canonical vacuum", "all", "2000-3000", "0 canonical rows", "Start with Cr II, Mn II, Co II, Sc II, V II and target-held intervals."),
    (3, "Fe II primary-lab deficit", "Fe II", "3000-3780", "canonical population exists but primary-lab tier is absent", "Source census and transition-level reconciliation; do not infer Fe I authority transfers to Fe II."),
    (4, "Cr II high-volume replacement", "Cr II", "2080-3780", "large VALD-only near-UV population", "Ingest Ward table to quarantine, normalize air/vacuum wavelengths, match with ambiguity reporting."),
    (5, "HFS-sensitive odd-Z ions", "Mn II; Co II; Sc II; V II; Cu I/II", "1150-3780", "HFS metadata incomplete or unverified", "Require component/constant provenance; never substitute unsplit gf silently."),
    (6, "non-Fe near-UV VALD dependence", "Cr I/II; Co I; Ti I/II; Mn I/II; Ni I; V I/II", "3000-3780", "dominant tier VALD3", "Rank laboratory searches by canonical count and actual target coverage."),
    (7, "photospheric-abundance literature gap", "all", "1150-3780", "spectra exist without a verified abundance-paper census", "Bibliographic follow-up per target; keep chromospheric and planetary/ISM lanes out of abundance authority."),
]


def regime(wavelength: float) -> str:
    for name, lo, hi in REGIMES:
        if lo <= wavelength < hi:
            return name
    raise ValueError(wavelength)


def build(out: Path) -> dict:
    out.mkdir(parents=True, exist_ok=True)
    gf = pd.read_csv(CANONICAL, low_memory=False)
    uv = gf[(gf.wavelength_air_A >= 1150.0) & (gf.wavelength_air_A < 3780.0)].copy()
    uv["regime"] = uv.wavelength_air_A.map(regime)

    grid = pd.MultiIndex.from_product(
        [[r[0] for r in REGIMES], sorted(set(uv.species.dropna()))], names=["regime", "species"]
    ).to_frame(index=False)
    counts = uv.pivot_table(index=["regime", "species"], columns="gf_tier", values="line_id", aggfunc="count", fill_value=0).reset_index()
    inventory = grid.merge(counts, how="left", on=["regime", "species"]).fillna(0)
    for col in ["LAB", "NIST-C+", "OTHER", "VALD3"]:
        if col not in inventory:
            inventory[col] = 0
        inventory[col] = inventory[col].astype(int)
    inventory["total"] = inventory[["LAB", "NIST-C+", "OTHER", "VALD3"]].sum(axis=1)
    inventory = inventory[inventory.total > 0].sort_values(["regime", "total", "species"], ascending=[True, False, True])
    inventory.to_csv(out / "current_uv_inventory.csv", index=False)

    pd.DataFrame(LAB_SOURCES, columns=["species", "source", "published_coverage_A", "published_transition_count", "published_uncertainty", "doi", "evidence_tier", "codex_status", "audit_note"]).to_csv(out / "laboratory_source_census.csv", index=False)
    pd.DataFrame(PAPERS, columns=["work", "citation", "targets", "coverage_A", "instrument", "science_lane", "authority", "audit_note", "url"]).to_csv(out / "stellar_uv_paper_census.csv", index=False)
    pd.DataFrame(GAPS, columns=["priority", "gap_class", "species", "coverage_A", "evidence", "next_action"]).to_csv(out / "candidate_gap_analysis.csv", index=False)

    holdings = pd.read_csv(HOLDINGS)
    hu = holdings[holdings.instrument_id.isin(["hst_stis", "hst_cos"])].copy()
    hu[["holding_id", "system_id", "instrument_id", "manifest_path", "evidence_state", "source_issue_ids"]].to_csv(out / "target_uv_holdings.csv", index=False)
    targets = sorted(hu.system_id.unique())
    overlap = []
    display = {"alpha_cen_a": "Alpha Cen A", "alpha_cen_b": "Alpha Cen B", "procyon": "Procyon"}
    for target in targets:
        name = display.get(target, target)
        mentions = [p for p in PAPERS if name.lower() in p[2].lower()]
        abundance = [p for p in mentions if p[5] == "PHOTOSPHERIC_ABUNDANCE"]
        overlap.append((target, name, ";".join(sorted(hu[hu.system_id == target].holding_id)), bool(mentions),
                        "CANDIDATE_PRESENT_NEEDS_VERIFICATION" if abundance else "NO_PHOTOSPHERIC_ABUNDANCE_PAPER_IN_CENSUS",
                        "Paper mention is not adoption: verify line list, gf authority, photospheric lane, and target component."))
    pd.DataFrame(overlap, columns=["system_id", "target", "uv_holdings", "appears_in_census", "abundance_literature_status", "routing_note"]).to_csv(out / "target_literature_overlap.csv", index=False)

    verdicts = []
    for species, group in uv.groupby("species"):
        tiers = group.gf_tier.value_counts()
        lab = int(tiers.get("LAB", 0)); total = len(group)
        verdict = "PARTIAL_PRIMARY_LAB" if lab else "NO_PRIMARY_LAB_IN_CANONICAL_UV"
        route = "KEEP_CURRENT_AND_AUDIT_GAPS" if lab else "SOURCE_CENSUS_BEFORE_ANY_ADOPTION"
        verdicts.append((species, total, lab, int(tiers.get("NIST-C+", 0)), int(tiers.get("VALD3", 0)), verdict, route))
    pd.DataFrame(verdicts, columns=["species", "canonical_uv_rows", "lab_rows", "nist_cplus_rows", "vald3_rows", "verdict", "route"]).sort_values(["canonical_uv_rows", "species"], ascending=[False, True]).to_csv(out / "per_species_verdict.csv", index=False)

    summary = {
        "ticket": "RYA-1061", "window_A": {"lower_inclusive": 1150.0, "upper_exclusive": 3780.0},
        "canonical_rows": int(len(uv)),
        "canonical_by_regime": {name: int(((uv.wavelength_air_A >= lo) & (uv.wavelength_air_A < hi)).sum()) for name, lo, hi in REGIMES},
        "canonical_by_tier": {str(k): int(v) for k, v in uv.gf_tier.value_counts().items()},
        "uv_holdings": int(len(hu)), "targets_with_uv_holdings": targets,
        "mutation_policy": "READ_ONLY_AUDIT_NO_CANONICAL_OR_HOLDING_MUTATIONS",
        "important_limits": [
            "Published source counts are paper-wide unless explicitly described as in-window.",
            "No transition is recommended for canonical adoption until the source table is acquired and reconciled line by line.",
            "Air/vacuum wavelength normalization is required below any match or overlap claim.",
            "Chromosphere/TR/corona and planetary/ISM/CSM studies cannot serve as photospheric abundance authority."
        ],
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    readme = f"""# RYA-1061 — UV completeness audit

This is a read-only evidence audit for **1150 <= wavelength < 3780 A**. It changes no
canonical gf row, line pool, holding, product, abundance, or astrophysical-gf authority.

## Verdict

The committed canonical list contains **{len(uv):,}** rows in the audit window, all in
3000–3780 A. FUV and NUV each contain **zero** canonical rows. The near-UV list is
{int((uv.gf_tier == 'VALD3').sum()):,}/{len(uv):,} VALD3 rows; only
{int((uv.gf_tier == 'LAB').sum())} rows carry the primary-laboratory tier. This is a
coverage discontinuity, not evidence that no transitions exist below 3000 A.

RYA-836 remains the verified Fe I anchor: 105 primary-lab Fe I transitions below 3780 A
in its source census, of which 61 lie in 3000–3780 A. The current canonical inventory is
not equivalent to that source census (it carries 71 total LAB rows across species), so
source availability, canonical adoption, and target reachability must remain separate.

The highest-leverage acquisition is Ward et al. (2023) Cr II (268 published transitions
over 208–414 nm, clipped here at 378 nm), followed by transition-table verification for
Mn II, Co II, Sc I/II, and V II. These are **source candidates**, not approved canonical
rows. No line-level adoption verdict is possible until tables are acquired, wavelength
frames normalized, ions matched, HFS/isotope structure retained, and ambiguous matches
reported.

## Target and literature routing

The holdings registry has four UV holdings: Alpha Cen A/STIS, Alpha Cen B/STIS,
Procyon/STIS, and Procyon/COS. ASTRAL is supporting spectral/identification evidence.
Pagano et al. is routed to the chromosphere/TR/corona lane. The Alpha Cen Be II 3130 work
is a candidate photospheric-abundance source requiring exact bibliographic and gf/blend
verification. The Procyon GHRS atlas is supporting chromospheric/identification evidence;
this census found no Procyon photospheric-abundance paper to pair with its UV spectra.

## Files

- `current_uv_inventory.csv`: per-species, per-regime canonical counts by gf tier.
- `laboratory_source_census.csv`: source-level evidence and acquisition routing.
- `stellar_uv_paper_census.csv`: strict science-lane and authority classification.
- `target_uv_holdings.csv`: UV rows derived from the holdings registry.
- `target_literature_overlap.csv`: target-level spectra/literature gap.
- `candidate_gap_analysis.csv`: prioritized gaps and next actions.
- `per_species_verdict.csv`: current-state canonical verdicts only.
- `summary.json`: machine-readable headline and audit limits.

Reproduce with `python3 scripts/rya1061_uv_completeness_audit.py`.
"""
    (out / "README.md").write_text(readme)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    print(json.dumps(build(args.out), indent=2))


if __name__ == "__main__":
    main()
