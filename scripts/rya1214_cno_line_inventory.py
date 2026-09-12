#!/usr/bin/env python3
"""RYA-1214 Step 2 — the measurable CNO line inventory, per element x band x line type.

This sizes the product campaign. It answers four questions per cell and refuses to
average over the difference between them:

  1. How many CNO lines does `canonical_gf` hold in this band?
  2. How many of those carry a GRADED gf after Step 1?
  3. Which solar holdings actually cover the band -- i.e. is there a spectrum at all?
  4. What molecular indicators does the RYA-1136 crossmatch place in the band?

And separately, the set that decides what gets RUN: the AGSS21 Table 3 atomic indicators
from `atomic_source_census.csv`, re-joined to `canonical_gf`.

🔴 THE CENSUS'S OWN JOIN IS RE-DERIVED HERE, BECAUSE ITS TOLERANCE HID LINES THAT EXIST.

`atomic_source_census.csv` marks 24 of its 48 CNO indicator rows `ABSENT` from
`canonical_gf`. Four of them are not absent, and two of those four are the headline
forbidden indicators of the whole campaign:

    [C I] 872.7 nm   census 8727.120 A   store gf_002448 at 8727.139 A   (0.019 A away)
    [O I] 630.0 nm   census 6300.300 A   store gf_108773 at 6300.304 A   (0.004 A away)
    C I   658.8 nm   census 6587.610 A   store gf_002178 at 6587.610 A   (exact)
    C I   505.2 nm   census 5052.160 A   store gf_001680 at 5052.167 A   (0.007 A away)

That is the RYA-1109 shape exactly -- a match tolerance tighter than the source's own
printed precision reports a coverage hole that is not there. The census quotes its lines
to 0.01 nm (0.1 A), so a 0.1 A window is the source's OWN quantisation and anything
tighter is measuring rounding. This script therefore joins at 0.1 A with an EP check,
reports the per-row separation so the join can be audited, and prints BOTH verdicts --
the census's and this one -- rather than quietly replacing it.

⚠️ A wider window buys ambiguity, so uniqueness is required and a tie is refused. The
one real tie is C I 5052.16, where the store holds two VALD components 0.023 A apart
(gf_001679 / gf_001680) for a census row quoted at 0.01 nm: that is a REPORTED tie, not
a pick. Genuinely absent after all this: N I 10108.90 and the nine O I 926 nm lines.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from pipeline.band_policy import POLICIES  # noqa: E402

CANON = ROOT / "data" / "linelists" / "canonical_gf.csv"
CENSUS = ROOT / "data" / "audit" / "rya1136_cno_intake" / "atomic_source_census.csv"
MOLEC = ROOT / "data" / "audit" / "rya1136_cno_intake" / "molecular_physical_crossmatch.csv"
INSTRUMENTS = ROOT / "data" / "catalog" / "instrument_catalog.csv"
HOLDINGS = ROOT / "data" / "catalog" / "holdings_manifest_registry.csv"
OUT = ROOT / "data" / "audit" / "rya1214_cno_products"

#: The census prints wavelengths to 0.01 nm. 0.1 A IS that quantisation -- derived from
#: the source's precision, never chosen to make a count come out (RYA-1109).
CENSUS_TOL_A = 0.1
#: Second key. Wide enough to admit the census's own rounding of EP, tight enough that a
#: different multiplet cannot pass.
CENSUS_EP_TOL_EV = 0.02

GRADED_TIERS = ("LAB", "NIST-C+")


def band_of(w: float) -> str:
    for p in POLICIES:                      # first-match, H before NIR (RYA-1094)
        if p.lo_A <= w < p.hi_A:
            return p.name
    return "out-of-band"


def solar_coverage() -> dict[str, list[str]]:
    """Which solar holdings' instruments span each band. An UPPER BOUND, and says so.

    The instrument catalogue's span is the instrument's capability, not proof that a
    given holding's file reaches the edge of it. It is the right question for SIZING
    (is there any spectrum at all in this band?) and the wrong one for a product claim,
    which is what `rya1187_applicability_matrix` resolves per holding.
    """
    inst = pd.read_csv(INSTRUMENTS).set_index("instrument_id")
    hold = pd.read_csv(HOLDINGS)
    hold = hold[hold.system_id == "solar"]
    out: dict[str, list[str]] = {p.name: [] for p in POLICIES}
    for _, h in hold.iterrows():
        if h.instrument_id not in inst.index:
            continue
        lo = float(inst.loc[h.instrument_id, "wavelength_min_nm"]) * 10.0
        hi = float(inst.loc[h.instrument_id, "wavelength_max_nm"]) * 10.0
        for p in POLICIES:
            if lo < p.hi_A and hi > p.lo_A:
                out[p.name].append(h.holding_id)
    return out


def rejoin_census(canon: pd.DataFrame) -> pd.DataFrame:
    """Re-derive the AGSS21 indicator -> canonical_gf join at the source's own precision."""
    cen = pd.read_csv(CENSUS)
    cen = cen[cen.element.isin(["C", "N", "O"])].copy()
    rows = []
    for _, r in cen.iterrows():
        w, ep, sp = float(r.wavelength_air_A), r.lower_EP_eV, str(r.species)
        c = canon[canon.species == sp]
        near = c[(c.wavelength_air_A - w).abs() <= CENSUS_TOL_A]
        if pd.notna(ep):
            near = near[(near.excitation_potential_eV - float(ep)).abs() <= CENSUS_EP_TOL_EV]
        rec = {
            "element": r.element, "species": sp, "line_label": r.line_label,
            "wavelength_air_A": w, "lower_EP_eV": ep,
            "band": band_of(w),
            "use_status": r.use_status,
            "census_join_status": r.join_status,
            "census_line_id": r.canonical_line_id,
            "n_store_candidates": int(len(near)),
        }
        if len(near) == 1:
            s = near.iloc[0]
            rec.update({"store_line_id": s.line_id,
                        "store_wavelength_A": float(s.wavelength_air_A),
                        "separation_A": round(float(s.wavelength_air_A) - w, 4),
                        "store_log_gf": float(s.log_gf),
                        "store_gf_tier": s.gf_tier,
                        "store_nist_grade": s.nist_grade,
                        "store_gf_sigma_dex": s.gf_sigma_dex,
                        "verdict": "PRESENT"})
        elif len(near) > 1:
            rec.update({"store_line_id": "|".join(near.line_id),
                        "separation_A": np.nan,
                        "verdict": "AMBIGUOUS_TIE_REPORTED_NOT_PICKED"})
        else:
            rec["verdict"] = "ABSENT_FROM_CANONICAL_GF"
        rows.append(rec)
    return pd.DataFrame(rows)


#: Per-species NIST pull files, for the indicator cross-check below.
_PULLS = ROOT / "data" / "linelists" / "primary_gf"
_PULL_TAG = {"C I": "CI", "C II": "CII", "N I": "NI", "N II": "NII",
             "O I": "OI", "O II": "OII"}


def indicator_gf_check(ind: pd.DataFrame, canon: pd.DataFrame) -> pd.DataFrame:
    """Every AGSS21 atomic indicator, store value vs the NIST ASD pull, LINE BY LINE.

    🔴 THIS EXISTS BECAUSE THE BULK RULE BURIED THE MOST IMPORTANT LINE IN THE CAMPAIGN.

    `rya1214_adjudicate_cno_gf` admits a NIST row only when the multiplicities agree,
    which is right in bulk and wrong for [O I] 6300.304. NIST serves that wavelength
    TWICE — the M1 transition (log gf -9.776, class B+) and an E2 partner (-12.20, class
    C+) — and those are two MULTIPOLE CHANNELS of one level pair, not two fine-structure
    components. VALD lists the M1 alone, so `hfs_n_components` is 1, the rule sees 2
    candidates against 1, and refuses. The line then falls out of every downstream
    comparison, including the one that checks protected rows against the pull.

    What it hides is a real disagreement on the dominant forbidden oxygen indicator:

        store   -9.717   class A    "NIST ASD v5.11 grade A (Storey & Zeippen 2000)",
                                    adjudicated RYA-367
        pull    -9.776   class B+   TP T4539,T5081  (the E2 partner adds 10^-12.2,
                                    i.e. nothing, so the sum is the M1 to 4 decimals)

    0.059 dex apart, and the store claims one accuracy class BETTER than the source it
    cites now publishes. That is RYA-1171's finding (the 777 triplet carried A+ where ASD
    publishes A) recurring one line over. The sibling [O I] 6363.776 reproduces the pull
    EXACTLY (-10.2580 vs -10.2581), which makes this line-specific rather than a scale
    offset and removes the comfortable explanation.

    Reported, not resolved — same referral as the 777 Tachiev offset.
    """
    rows = []
    for _, r in ind[ind.verdict == "PRESENT"].iterrows():
        sp = str(r.species)
        tag = _PULL_TAG.get(sp)
        if tag is None:
            continue
        n = pd.concat([pd.read_csv(_PULLS / f"nist_asd_{tag}_{b}.tsv", sep="\t")
                       for b in ("900_3000", "3000_25000")], ignore_index=True)
        n = n[n.log_gf.notna() & n.ei_eV.notna()]
        # Ladenburg gate, as in the adjudicator -- a non-physical row must not enter a
        # comparison either.
        gk = n.gi * n.fik / (1.4992e-16 * n.wavelength_A ** 2 * n["aki_s-1"])
        n = n[gk.between(0.5, 50.0)]
        store_w = float(r.store_wavelength_A)
        c = n[(n.wavelength_A - store_w).abs() <= 0.05]
        if pd.notna(r.lower_EP_eV):
            c = c[(c.ei_eV - float(r.lower_EP_eV)).abs() <= 0.02]
        if c.empty:
            rows.append({"element": r.element, "species": sp, "line_label": r.line_label,
                         "store_line_id": r.store_line_id, "store_wavelength_A": store_w,
                         "store_log_gf": r.store_log_gf, "store_gf_tier": r.store_gf_tier,
                         "store_nist_grade": r.store_nist_grade,
                         "n_pull_rows": 0, "verdict": "NO NIST ROW IN PULL"})
            continue
        gf = 10.0 ** c.log_gf
        summed = float(np.log10(gf.sum()))
        # ⚠️ THE WORST GRADE AMONG THE COMPONENTS THAT ACTUALLY CARRY THE LINE.
        # Worst-of-all is the safe rule for a fine-structure blend, where the components
        # are comparable. It is the WRONG rule here: [O I] 6300's E2 channel is 2.4 dex
        # weaker than its M1 partner and contributes 0.4% of the summed gf, so letting
        # its C+ set the grade would describe the line by a term that is not in it. The
        # 1% floor is two orders of magnitude below "comparable" and its effect is
        # stated, not assumed -- it changes exactly one line in this set.
        share = gf / gf.sum()
        carry = c[share >= 0.01]
        worst = max(carry.nist_grade, key=lambda g: {"AAA": .3, "AA": 1, "A+": 2, "A": 3,
                                                     "B+": 7, "B": 10, "C+": 18, "C": 25,
                                                     "D+": 40, "D": 50, "E": 100}.get(g, 999))
        dropped = sorted(set(c.nist_grade) - set(carry.nist_grade))
        d = summed - float(r.store_log_gf)
        rows.append({
            "element": r.element, "species": sp, "line_label": r.line_label,
            "store_line_id": r.store_line_id, "store_wavelength_A": store_w,
            "store_log_gf": r.store_log_gf, "store_gf_tier": r.store_gf_tier,
            "store_nist_grade": r.store_nist_grade,
            "n_pull_rows": int(len(c)),
            "n_pull_rows_carrying_1pct": int(len(carry)),
            "pull_log_gf_summed": round(summed, 4),
            "pull_grade_worst": worst,
            "pull_grades_below_1pct_ignored": "|".join(dropped),
            "pull_tp_codes": "|".join(sorted(set(c.ref_transition_probability.astype(str)))),
            "delta_pull_minus_store": round(d, 4),
            "verdict": ("AGREES" if abs(d) <= 0.02 else
                        "🔴 DISAGREES — reported to Ryan, not resolved"),
            "grade_verdict": ("" if str(r.store_nist_grade) in ("", "nan") or
                              str(r.store_nist_grade) == worst
                              else f"🔴 store claims {r.store_nist_grade}, pull publishes {worst}"),
        })
    return pd.DataFrame(rows)


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    canon = pd.read_csv(CANON, low_memory=False)
    cno = canon[canon.key_z.isin([6, 7, 8])].copy()
    cno["band"] = cno.wavelength_air_A.map(band_of)
    cno["graded"] = cno.gf_tier.isin(GRADED_TIERS)

    cov = solar_coverage()

    # ------------------------------------------------------- atomic inventory per cell
    atomic = (cno.groupby(["species", "band"])
              .agg(lines=("line_id", "size"), graded=("graded", "sum"),
                   lo_A=("wavelength_air_A", "min"), hi_A=("wavelength_air_A", "max"))
              .reset_index())
    atomic["element"] = atomic.species.str.split().str[0]
    atomic["line_type"] = "atomic"
    atomic["solar_holdings_covering_band"] = atomic.band.map(lambda b: len(cov.get(b, [])))
    atomic["holdings"] = atomic.band.map(lambda b: "|".join(cov.get(b, [])))

    # ---------------------------------------------------- molecular indicators per cell
    mol = pd.read_csv(MOLEC)
    # The crossmatch is in VACUUM nm. Band edges are AIR Angstroms, and at 1 um the two
    # frames differ by ~2.7 A -- small against a band edge, but converted rather than
    # assumed equal, because "nm vs A" is exactly how RYA-1190 read the wrong file.
    lam_vac_A = mol.wavelength_vac_nm.to_numpy(float) * 10.0
    s2 = (1e4 / lam_vac_A) ** 2                     # Ciddor/IAU air-vacuum, as in iSpec
    nref = 1.0 + 0.0000834254 + 0.02406147 / (130.0 - s2) + 0.00015998 / (38.9 - s2)
    mol["wavelength_air_A"] = lam_vac_A / nref
    mol["band"] = mol.wavelength_air_A.map(band_of)
    mol["element"] = mol.species.map(lambda s: "C" if s in ("C2", "CH", "12C16O", "CN")
                                     else ("N" if s in ("NH",) else "O"))
    # CN carries BOTH C and N -- it is AGSS21's largest nitrogen indicator (522 lines)
    # and a carbon one. Duplicated deliberately: collapsing it onto one element would
    # delete the nitrogen campaign's biggest cell.
    cn = mol[mol.species == "CN"].copy()
    cn["element"] = "N"
    mol = pd.concat([mol, cn], ignore_index=True)
    molecular = (mol.groupby(["element", "species", "band"])
                 .agg(lines=("source_row", "size"))
                 .reset_index())
    molecular["line_type"] = "molecular"
    molecular["graded"] = 0
    molecular["solar_holdings_covering_band"] = molecular.band.map(lambda b: len(cov.get(b, [])))

    inv = pd.concat([
        atomic[["element", "species", "band", "line_type", "lines", "graded",
                "solar_holdings_covering_band", "holdings"]],
        molecular[["element", "species", "band", "line_type", "lines", "graded",
                   "solar_holdings_covering_band"]],
    ], ignore_index=True).sort_values(["element", "line_type", "species", "band"])
    inv.to_csv(OUT / "cno_line_inventory.csv", index=False)

    # ---------------------------------------------------------- AGSS21 indicator re-join
    ind = rejoin_census(canon)
    ind.to_csv(OUT / "agss21_indicator_join.csv", index=False)
    chk = indicator_gf_check(ind, canon)
    chk.to_csv(OUT / "agss21_indicator_gf_check.csv", index=False)

    # ------------------------------------------------------------------------- report
    print("=== RYA-1214 Step 2 — measurable CNO line inventory ===\n")
    print("ATOMIC (canonical_gf, graded = gf_tier in %s)" % (GRADED_TIERS,))
    piv = atomic.pivot_table(index="species", columns="band",
                             values=["lines", "graded"], aggfunc="sum", fill_value=0)
    print(piv.to_string())
    print("\nMOLECULAR indicators (RYA-1136 crossmatch, 408 rows; CN counted under BOTH C and N)")
    print(molecular.pivot_table(index=["element", "species"], columns="band",
                                values="lines", aggfunc="sum", fill_value=0).to_string())
    print("\nSOLAR HOLDINGS COVERING EACH BAND (instrument-catalogue span — an UPPER BOUND)")
    for b, hs in cov.items():
        print(f"  {b:12s} {len(hs):2d}  {', '.join(hs) if hs else '(none)'}")

    print("\n=== AGSS21 Table 3 atomic indicators, re-joined at the census's own 0.1 A ===")
    print(ind.verdict.value_counts().to_string())
    moved = ind[(ind.census_join_status == "ABSENT") & (ind.verdict != "ABSENT_FROM_CANONICAL_GF")]
    print(f"\n🔴 rows the census called ABSENT that ARE in canonical_gf: {len(moved)}")
    if len(moved):
        print(moved[["element", "species", "line_label", "wavelength_air_A", "band",
                     "store_line_id", "separation_A", "store_gf_tier",
                     "verdict"]].to_string(index=False))
    still = ind[ind.verdict == "ABSENT_FROM_CANONICAL_GF"]
    print(f"\ngenuinely absent after the re-join: {len(still)}")
    if len(still):
        print(still[["element", "species", "line_label", "wavelength_air_A",
                     "band"]].to_string(index=False))

    bad = chk[chk.verdict.str.startswith("🔴")] if len(chk) else chk
    gbad = chk[chk.grade_verdict.astype(str).str.startswith("🔴")] if len(chk) else chk
    print(f"\n=== AGSS21 indicator gf: store vs the NIST ASD pull, line by line "
          f"({len(chk)} lines) ===")
    print(f"  value disagreements (> 0.02 dex): {len(bad)}")
    if len(bad):
        print(bad[["element", "line_label", "store_line_id", "store_log_gf",
                   "pull_log_gf_summed", "delta_pull_minus_store", "store_nist_grade",
                   "pull_grade_worst"]].to_string(index=False))
    print(f"  grade over-claims               : {len(gbad)}")
    if len(gbad):
        print(gbad[["element", "line_label", "store_line_id",
                    "grade_verdict"]].to_string(index=False))

    prov = {
        "ticket": "RYA-1214", "step": "2 — measurable CNO line inventory",
        "agss21_indicator_gf_check": {
            "n_lines": int(len(chk)),
            "value_disagreements_gt_0p02_dex": int(len(bad)),
            "grade_over_claims": int(len(gbad)),
            "why_separate_from_the_bulk_adjudication":
                "the multiplicity rule refuses [O I] 6300.304 — NIST serves it twice as "
                "the M1 and E2 MULTIPOLE channels of one level pair, not as two "
                "fine-structure components — so the line falls out of every bulk "
                "comparison. This check is per-line and unconditional.",
        },
        "canonical_cno_rows": int(len(cno)),
        "canonical_cno_graded_rows": int(cno.graded.sum()),
        "molecular_crossmatch_rows": int(len(pd.read_csv(MOLEC))),
        "census_rejoin": {
            "tolerance_A": CENSUS_TOL_A,
            "tolerance_basis": "the census prints 0.01 nm = 0.1 A; the window IS the "
                               "source's own quantisation (RYA-1109)",
            "ep_tolerance_eV": CENSUS_EP_TOL_EV,
            "verdicts": {k: int(v) for k, v in ind.verdict.value_counts().items()},
            "census_ABSENT_rows_that_are_present": int(len(moved)),
            "genuinely_absent": int(len(still)),
        },
        "solar_band_coverage_is_an_upper_bound": "instrument-catalogue span, not a "
            "per-holding file span; rya1187_applicability_matrix resolves the latter",
    }
    (OUT / "cno_line_inventory.prov.json").write_text(json.dumps(prov, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
