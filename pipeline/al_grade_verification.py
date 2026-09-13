"""RYA-1134: source-backed Al dispositions, separate from abundance measurement.

The frozen intake is the denominator, never the authority for its own grades.
Source records, candidates, and decisions are retained separately. No source input
or canonical gf is rewritten. Missing uncertainty is not a numeric fallback.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from pipeline import line_match
from pipeline.wavelength_util import air_to_vac

ROOT = Path(__file__).resolve().parents[1]
INTAKE = Path("data/audit/rya1132_al_intake")
SOURCE = Path("data/reference/al_grade_verification_rya1134")
VERIFIED = Path("data/results/rya1134/verified_v1")


def raw_component_records(root: Path) -> pd.DataFrame:
    """Recover term/J identities from the preserved HFS-on VALD deliveries."""
    names = ("vald_solar_fuv_1150_2000_hfson_raw.txt", "vald_solar_nearuv_2000_3780_hfson_raw.txt",
             "vald_solar_raw.txt", "vald_solar_redopt_6910_9500_hfson_raw.txt",
             "vald_solar_ir_9500_17000_hfson_raw.txt", "vald_solar_ir_17000_25000_hfson_raw.txt")
    records = []
    for name in names:
        lines = (root / "data/linelists" / name).read_text().splitlines()
        for i, line in enumerate(lines):
            if not re.match(r"'Al [12]'", line):
                continue
            fields = [v.strip() for v in line.split(",")]
            ion = "I" if fields[0] == "'Al 1'" else "II"
            records.append(dict(species="Al " + ion, wavelength_air=float(fields[1]),
                loggf=float(fields[2]), lower_EP=float(fields[3]), lower_J=float(fields[4]),
                upper_EP=float(fields[5]), upper_J=float(fields[6]),
                lower_term=lines[i+1].strip("' "), upper_term=lines[i+2].strip("' "),
                source_annotation=lines[i+3], raw_path="data/linelists/" + name, raw_line_number=i+1))
    return pd.DataFrame(records)


def component_audit(ledger: pd.DataFrame, root: Path):
    """Re-sum the actual VALD components and preserve every contributing row.

    Term/J assignments are recovered from the raw deliveries, never inferred
    from the merged component count. Multi-transition groups remain explicit.
    """
    sys.path.insert(0, str(root / "scripts"))
    from rya1001_al_census import HFS_SPAN_CM1, HFS_EP_TOL_EV
    from line_accounting_rya709 import DEPTH_LO, DEPTH_HI
    ll = pd.read_csv(root / "data/linelists/linelist_solar.csv", low_memory=False)
    ll["merged_row_number"] = np.arange(len(ll)) + 2
    ll = ll[ll.element.eq("Al")].sort_values("wavelength_air_A", kind="stable").copy()
    separation = 1e8 * ll.wavelength_air_A.diff().abs() / (ll.wavelength_air_A * ll.wavelength_air_A.shift())
    same = separation.le(HFS_SPAN_CM1) & ll.excitation_potential_eV.diff().abs().le(HFS_EP_TOL_EV) & ll.ion.eq(ll.ion.shift())
    ll["group"] = (~same).cumsum()
    raw = raw_component_records(root)
    groups, components = [], []
    for group_id, part in ll.groupby("group", sort=False):
        gf = np.power(10., part.log_gf.to_numpy(float))
        groups.append(dict(species="Al " + str(part.ion.iloc[0]),
            wavelength_air=float(np.average(part.wavelength_air_A, weights=gf)),
            lower_EP=float(part.excitation_potential_eV.min()), group=group_id,
            sum_loggf=float(np.log10(gf.sum())), n_components=len(part),
            feature_depth=float(part.central_depth.max()), blend=bool(part.blend_flag.any())))
    groups = pd.DataFrame(groups)
    out = ledger.copy()
    for i, r in out.iterrows():
        found = unique_match(r.to_dict(), groups, .005, .0001)
        if len(found) != 1:
            out.loc[i, "hfs_component_status"] = "NO_UNIQUE_RAW_COMPONENT_GROUP"
            out.loc[i, "feature_depth"] = np.nan
            out.loc[i, "codex_membership"] = out.loc[i, "deep_membership"] = "HOLD" if r.verified_source_class == "PRIMARY_LABORATORY" else "NONMEMBER"
            continue
        group = found.iloc[0]
        part = ll[ll.group.eq(group.group)]
        out.loc[i, "feature_depth"] = group.feature_depth
        out.loc[i, "component_count"] = group.n_components
        out.loc[i, "raw_component_sum_loggf"] = group.sum_loggf
        out.loc[i, "blend_flag"] = group.blend
        proof = []
        for component in part.to_dict("records"):
            target = dict(species=r.species, wavelength_air=component["wavelength_air_A"], lower_EP=component["excitation_potential_eV"])
            found_raw = unique_match(target, raw, .00001, .000001)
            found_raw = found_raw[(found_raw.loggf - component["log_gf"]).abs().le(.000001)]
            proof.append(found_raw.iloc[0].to_dict() if len(found_raw) == 1 else None)
        identities = {(p["lower_J"], p["upper_J"], p["lower_term"], p["upper_term"]) for p in proof if p}
        intact = all(p is not None for p in proof) and len(identities) == 1
        out.loc[i, "hfs_component_status"] = ("RAW_PHYSICAL_COMPONENTS_RECOVERED" if intact else
            "MULTIPLE_FINE_STRUCTURE_TRANSITIONS" if len(identities) > 1 else "RAW_COMPONENT_IDENTITY_UNRESOLVED")
        lab = r.verified_source_class == "PRIMARY_LABORATORY"
        out.loc[i, "codex_membership"] = "MEMBER" if lab and DEPTH_LO <= group.feature_depth <= DEPTH_HI else "NONMEMBER"
        out.loc[i, "deep_membership"] = "MEMBER" if lab and group.feature_depth > DEPTH_HI else "NONMEMBER"
        # Preserve normalized weights as proof only; deployment also requires
        # the source-identity and exact-holding gates below.
        for component, p in zip(part.to_dict("records"), proof):
            weight = 10**(component["log_gf"] - group.sum_loggf)
            components.append(dict(canonical_line_id=r.canonical_line_id,
                merged_row_number=component["merged_row_number"], species=r.species,
                wavelength_air_A=component["wavelength_air_A"], lower_EP=component["excitation_potential_eV"],
                component_loggf=component["log_gf"], source=component["loggf_source"],
                gf_fraction=weight, parent_sum_loggf=group.sum_loggf,
                source_total_loggf=r.verified_loggf,
                rescaled_component_loggf=(r.verified_loggf + math.log10(weight)) if intact and math.isfinite(r.verified_loggf) else np.nan,
                raw_path=p["raw_path"] if p else "", raw_line_number=p["raw_line_number"] if p else "",
                lower_J=p["lower_J"] if p else np.nan, upper_J=p["upper_J"] if p else np.nan,
                raw_lower_EP=p["lower_EP"] if p else np.nan, raw_upper_EP=p["upper_EP"] if p else np.nan,
                lower_term=p["lower_term"] if p else "", upper_term=p["upper_term"] if p else "",
                source_annotation=p["source_annotation"] if p else "",
                component_semantics="RAW_PHYSICAL_COMPONENT" if intact else "UNRESOLVED_GROUP_MEMBER"))
    return out, pd.DataFrame(components)


def loggf_from_a(wavelength_vac_A: float, upper_j: float, a_s: float) -> float:
    """Einstein A to absorption gf; wavelength is explicitly VACUUM Angstrom."""
    if not all(math.isfinite(x) for x in (wavelength_vac_A, upper_j, a_s)) or a_s <= 0:
        raise ValueError("finite positive transition probability required")
    return math.log10(1.49919e-16 * (2 * upper_j + 1) * wavelength_vac_A**2 * a_s)


def fractional_bound_dex(percent: float) -> float:
    """Positive logarithmic bound, NOT an implicit Gaussian 1-sigma conversion."""
    return math.log10(1 + percent / 100) if math.isfinite(percent) and percent > 0 else np.nan


def unique_match(row: dict, sources: pd.DataFrame, tol_A: float, ep_tol: float):
    """Shared strict physical matcher, with no nearest-choice among coincidences."""
    same = sources[sources.species.eq(row["species"])].copy()
    if same.empty:
        return same
    result = line_match.match([row["wavelength_air"]], same.wavelength_air,
                             want_ep=[row["lower_EP"]], src_ep=same.lower_EP,
                             require_ep=True, tol_A=tol_A, ep_tol_eV=ep_tol)
    candidates = same[(same.wavelength_air - row["wavelength_air"]).abs().le(tol_A)
                      & (same.lower_EP - row["lower_EP"]).abs().le(ep_tol)]
    if len(candidates) == 1 and result.n_resolved != 1:
        raise RuntimeError("strict matcher disagrees with candidate enumeration")
    return candidates


def _j(level: str) -> float:
    # Source strings use either LaTeX or CDS subscript notation.
    values = re.findall(r"(\d+)\s*/\s*(\d+)", level)
    return float(values[-1][0]) / float(values[-1][1]) if values else np.nan


def burheim_sources(root: Path) -> pd.DataFrame:
    lab = pd.read_csv(root / "data/reference/al_gf_lab/al1_lab_loggf.csv")
    quantities = pd.read_csv(root / SOURCE / "burheim_table1_quantities.csv")
    rows = []
    for q, r in zip(quantities.to_dict("records"), lab.to_dict("records")):
        if abs(q["wavelength_vac_A"] - r["lam_vac_A"]) > .00001:
            raise ValueError("Burheim table transcription/order mismatch")
        derived = loggf_from_a(r["lam_vac_A"], _j(r["upper_level"]), q["Aul_s"])
        # Table A is rounded to 3 significant digits; loggf to 2-4 decimals.
        if abs(derived - r["loggf"]) > .006:
            raise ValueError(f"Burheim A->gf does not reproduce {q['source_record_id']}")
        rows.append({**q, "species": "Al I", "wavelength_air": r["wavelength_air_A"],
                     "lower_EP": r["elo_eV"], "upper_EP": r["eup_eV"],
                     "lower_level": r["lower_level"], "upper_level": r["upper_level"],
                     "lower_J": _j(r["lower_level"]), "upper_J": _j(r["upper_level"]),
                     "loggf": r["loggf"], "recomputed_loggf": derived,
                     "gf_bound_dex": fractional_bound_dex(r["unc_pct"]),
                     "uncertainty_basis": "published fractional uncertainty; confidence not specified",
                     "source_class": "PRIMARY_LABORATORY" if q["lifetime_kind"] == "experimental"
                     else "MIXED_LAB_THEORY", "source_key": "Burheim2023",
                     "source_path": "data/reference/al_gf_lab/al1_lab_loggf.csv",
                     "source_semantics": "FINE_STRUCTURE_TOTAL", "priority": 0 if q["lifetime_kind"] == "experimental" else 1})
    if len(rows) != 12:
        raise ValueError("Burheim Table 1 must contain exactly twelve derived gf records")
    return pd.DataFrame(rows)


def vujnovic_sources(root: Path) -> pd.DataFrame:
    """Re-read CDS bytes, including flags the intake dropped; keep all 106 rows.

    Lower-level energies below come from the ground fine-structure separation
    (112.061 cm-1, independently checked RYA-1141 A4). Other lower levels remain
    unresolved here; their source values are retained without wavelength-only use.
    """
    rows = []
    for table, count in ((2, 29), (3, 22), (4, 24), (5, 31)):
        lines = (root / "data/reference/vujnovic2002_al/raw" / f"table{table}.dat").read_text().splitlines()
        if len(lines) != count:
            raise ValueError(f"CDS table {table} count changed")
        for n, s in enumerate(lines, 1):
            def number(a, b):
                return float(s[a:b]) if s[a:b].strip() else np.nan
            if table == 2:
                upper, lower, wave = s[:20].strip(), s[21:36].strip(), number(37, 46)
                a, pct = number(70, 75) * 1e8, number(77, 79)
                flags = dict(aki_limit=s[69:70].strip(), aki_note=s[75:76].strip(),
                             uncertainty_limit=s[76:77].strip(), intensity_limit=s[47:48].strip(),
                             theory_flag="", branching_limit="")
                other_a, other_ref = number(80, 85) * 1e8, s[86:].strip()
            elif table == 5:
                upper, lower, wave = s[:13].strip(), s[14:27].strip(), number(28, 36)
                a, pct = number(48, 55) * 1e8, number(58, 60)
                flags = dict(aki_limit=s[47:48].strip(), aki_note="", uncertainty_limit=s[57:58].strip(),
                             intensity_limit="", theory_flag=s[36:37].strip(), branching_limit=s[37:38].strip())
                other_a, other_ref = number(61, 66) * 1e8, "ChangWang1987"
            else:
                upper, lower, wave = s[:19].strip(), s[20:35].strip(), number(36, 44)
                a = pct = other_a = np.nan
                other_ref = ""
                flags = dict.fromkeys(("aki_limit", "aki_note", "uncertainty_limit", "intensity_limit", "theory_flag", "branching_limit"), "")
            # Table 2 rows 5/6 use theoretical LS ratios; row 10 was NOT observed.
            theory = bool(flags["theory_flag"]) or (table == 2 and n in (5, 6, 10))
            ep = (0.0 if "1/2" in lower else 112.061 / 8065.544005) if table == 2 and lower.startswith("3p ") else np.nan
            finite = math.isfinite(a) and a > 0 and not flags["aki_limit"]
            bound = fractional_bound_dex(pct) if finite and not flags["uncertainty_limit"] else np.nan
            vac = float(air_to_vac(np.array([wave]))[0]) if wave >= 2000 else wave
            derived = loggf_from_a(vac, _j(upper), a) if finite and math.isfinite(_j(upper)) else np.nan
            rows.append({"source_record_id": f"vuj2002_t{table}_{n:03d}", "source_key": "Vujnovic2002",
                         "source_path": f"data/reference/vujnovic2002_al/raw/table{table}.dat",
                         "species": "Al I" if table in (2, 3) else "Al II", "wavelength_air": wave,
                         "wavelength_vac_A": vac, "lower_EP": ep, "lower_level": lower, "upper_level": upper,
                         "lower_J": _j(lower), "upper_J": _j(upper),
                         "Aul_s": a, "unc_pct": pct, "loggf": derived, "recomputed_loggf": derived,
                         "gf_bound_dex": bound, "source_class": "MIXED_LAB_THEORY" if theory else "PRIMARY_LABORATORY",
                         "source_semantics": "RATIO_ONLY" if not finite else "FINE_STRUCTURE_TOTAL",
                         "uncertainty_basis": "published bound; Vujnovic sums relative errors, not Gaussian quadrature",
                         "priority": 2, "other_Aul_s": other_a, "other_reference": other_ref, **flags})
    return pd.DataFrame(rows)


def evaluated_sources(root: Path) -> pd.DataFrame:
    rows = []
    # This is the same fixed snapshot consumed by the frozen census, not a fresh search.
    path = "data/linelists/primary_gf/nist_asd_AlI_6600_42000.tsv"
    for n, r in enumerate(pd.read_csv(root / path, sep="\t").to_dict("records"), 1):
        if not math.isfinite(r["log_gf"]):
            continue
        raw = str(r["ei_ek_raw"]).split("-")
        eup = float(re.sub(r"[^0-9.]", "", raw[-1])) if len(raw) == 2 else np.nan
        ref = str(r["ref_transition_probability"])
        # u36/LS is the Opacity Project chain; keep unknown source codes unresolved.
        cls = "EVALUATED_THEORY" if "u36" in ref else "EVALUATED_SOURCE_UNRESOLVED"
        pct = r["nist_acc_pct"]
        bound = fractional_bound_dex(pct) if str(r["nist_grade"]) != "E" else np.nan
        vac = float(air_to_vac(np.array([r["wavelength_A"]]))[0])
        g_upper = r["gi"] * r["fik"] / (1.49919e-16 * vac**2 * r["aki_s-1"])
        # A and f are rounded in the frozen snapshot. A unique integer weight
        # within 3% reconstructs J; otherwise leave it unestablished. The raw
        # VALD J is an independent referee, checked downstream.
        g_int = round(g_upper)
        upper_j = (g_int - 1) / 2 if g_int > 0 and abs(g_upper / g_int - 1) < .03 else np.nan
        rows.append(dict(source_record_id=f"nist_snapshot_{n:04d}", species="Al I",
                         wavelength_air=r["wavelength_A"], lower_EP=r["ei_eV"], upper_EP=eup,
                         lower_J=(r["gi"]-1)/2, upper_J=upper_j,
                         upper_statistical_weight_from_A_f=g_upper,
                         loggf=math.log10(r["gi"] * r["fik"]), recomputed_loggf=math.log10(r["gi"] * r["fik"]),
                         gf_bound_dex=bound, source_key="NIST_ASD_Al_snapshot", source_path=path,
                         source_class=cls, nist_grade=r["nist_grade"], source_reference=ref,
                         uncertainty_basis="NIST accuracy bound; E is >50%, no finite bound assigned",
                         source_semantics="FINE_STRUCTURE_TOTAL", priority=3))
    # NL2017 provides source-resolved level identities outside the limited ASD snapshot.
    path = "data/reference/asplund2021_al/nordlander_lind_2017_analysis_lines.csv"
    for n, r in enumerate(pd.read_csv(root / path).to_dict("records"), 1):
        authority = str(r["gf_source_per_line"])
        if "TOPbase" not in authority:
            continue  # Vujnovic is judged from the original source, never double-counted.
        rows.append(dict(source_record_id=f"nl2017_a1_{n:03d}", species="Al I",
                         wavelength_air=r["wavelength_air_A"], lower_EP=r["elo_eV"], upper_EP=r["eup_eV"],
                         lower_level=r["lower_level"], upper_level=r["upper_level"],
                         lower_J=_j(r["lower_level"]), upper_J=_j(r["upper_level"]),
                         loggf=r["loggf"], recomputed_loggf=r["loggf"], gf_bound_dex=r["loggf_sigma_dex"],
                         source_key="NordlanderLind2017_TOPbase", source_path=path,
                         source_class="EVALUATED_THEORY", source_reference=authority,
                         uncertainty_basis="NL2017 Table A.1, inherited from NIST; not independent laboratory evidence",
                         source_semantics="FINE_STRUCTURE_TOTAL", priority=4))
    return pd.DataFrame(rows)


def johnson_source() -> pd.DataFrame:
    vac = float(air_to_vac(np.array([2669.157]))[0])
    return pd.DataFrame([dict(source_record_id="johnson1986_2669", source_key="Johnson1986",
        source_path="data/reference/al_grade_verification_rya1134/source_review.md",
        species="Al II", wavelength_air=2669.157, lower_EP=0.0,
        lower_level="3s2 1S0", upper_level="3s3p 3P1o", lower_J=0., upper_J=1., upper_EP=12398.419843 / vac,
        loggf=loggf_from_a(vac, 1., 3330.), recomputed_loggf=loggf_from_a(vac, 1., 3330.),
        gf_bound_dex=fractional_bound_dex(100 * 230 / 3330),
        uncertainty_basis="90% confidence bound, NOT a one-sigma uncertainty",
        source_class="PRIMARY_LABORATORY", source_semantics="FINE_STRUCTURE_TOTAL", priority=0)])


def source_candidates(row: dict, sources: pd.DataFrame) -> pd.DataFrame:
    candidates = []
    for key, subset in sources.groupby("source_key", sort=False):
        # Precision from source: NL2017 EP is 0.001 eV; other sources are finer.
        tol, ep_tol = ((.008, .0006) if key == "NordlanderLind2017_TOPbase" else
                       (.01, .0001) if key == "Vujnovic2002" else (.008, .0001))
        found = unique_match(row, subset, tol, ep_tol)
        if len(found):
            found = found.copy()
            found["identity_status"] = "UNIQUE_SPECIES_WAVELENGTH_EP" if len(found) == 1 else "AMBIGUOUS_FINE_STRUCTURE"
            candidates.append(found)
    return pd.concat(candidates, ignore_index=True) if candidates else sources.iloc[:0].copy()


def build(root: Path = ROOT) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    manifest = pd.read_csv(root / INTAKE / "al_line_manifest.csv").fillna("")
    if len(manifest) != 505 or not manifest.canonical_line_id.is_unique:
        raise ValueError("frozen 505-row denominator changed")
    sources = pd.concat([burheim_sources(root), vujnovic_sources(root), evaluated_sources(root), johnson_source()], ignore_index=True)
    bibliography = pd.read_csv(root / "data/refs/bibliography.csv").fillna("").set_index("key")
    keys = {"Burheim2023": "burheim2023", "Vujnovic2002": "vujnovic2002",
            "NIST_ASD_Al_snapshot": "kelleher2008_al", "NordlanderLind2017_TOPbase": "nordlander_lind2017",
            "Johnson1986": "johnson_smith_parkinson1986"}
    sources["bibliography_key"] = sources.source_key.map(keys)
    chains = {"Burheim2023": "burheim2023", "Vujnovic2002": "vujnovic2002",
              "NIST_ASD_Al_snapshot": "kelleher2008_al",
              "NordlanderLind2017_TOPbase": "nordlander_lind2017;kelleher2008_al;mendoza1995_al",
              "Johnson1986": "johnson_smith_parkinson1986;trabert1999_al"}
    sources["bibliography_chain"] = sources.source_key.map(chains)
    op = sources.source_key.eq("NIST_ASD_Al_snapshot") & sources.source_class.eq("EVALUATED_THEORY")
    sources.loc[op, "bibliography_chain"] += ";mendoza1995_al"
    lifetime_keys = {"Buurman1986": "buurman1986_al", "BuurmanDonszelmann1990": "buurman_donszelmann1990_al",
                     "Papoulia2019": "papoulia2019"}
    for i, record in sources[sources.source_key.eq("Burheim2023")].iterrows():
        sources.loc[i, "bibliography_chain"] += ";" + lifetime_keys[record.lifetime_source]
    for chain in sources.bibliography_chain:
        if any(key not in bibliography.index for key in chain.split(";")):
            raise ValueError(f"unresolved bibliography chain: {chain}")
    sources["source_doi"] = sources.bibliography_key.map(bibliography.doi)
    sources["citation_role"] = sources.source_key.map(lambda key: "source-lineage compilation; numeric values from ASD snapshot"
        if key == "NIST_ASD_Al_snapshot" else "published source of the retained quantities")
    if sources.bibliography_key.isna().any() or sources.source_doi.isna().any():
        raise ValueError("source missing from authoritative bibliography")
    decisions, links = [], []
    for row in manifest.to_dict("records"):
        row["wavelength_air"] = float(row["wavelength_air"])
        row["lower_EP"] = float(row["lower_EP"]) if row["lower_EP"] != "" else np.nan
        candidates = source_candidates(row, sources)
        for candidate in candidates.to_dict("records"):
            links.append({"canonical_line_id": row["canonical_line_id"], **candidate})
        usable = candidates[candidates.identity_status.eq("UNIQUE_SPECIES_WAVELENGTH_EP")
                            & candidates.gf_bound_dex.notna()
                            & candidates.loggf.notna()
                            & candidates.source_semantics.eq("FINE_STRUCTURE_TOTAL")
                            & ~candidates.source_class.eq("EVALUATED_SOURCE_UNRESOLVED")] if len(candidates) else candidates
        chosen = usable.sort_values("priority", kind="stable").iloc[0].to_dict() if len(usable) else None
        # A transition total cannot silently stand in for the 11254.9 unresolved feature.
        blend_total = row["canonical_line_id"] == "alphys_I_11254.9239_0407"
        if blend_total:
            chosen = None
        why = "BEST_AVAILABLE_SOURCE_WITH_IDENTITY_AND_BOUND" if chosen else (
            "COMPONENT_NOT_BLEND_TOTAL" if blend_total else
            "SOURCE_CANDIDATES_UNRESOLVED_OR_NO_FINITE_UNCERTAINTY" if len(candidates) else
            "NO_QUALIFIED_MATCH_IN_FROZEN_SOURCE_ASSETS")
        decisions.append(dict(canonical_line_id=row["canonical_line_id"], species=row["species"],
            wavelength_air=row["wavelength_air"], lower_EP=row["lower_EP"], historical_band=row["band"],
            historical_gf_source=row["gf_source"], historical_gf_grade=row["gf_grade"],
            historical_source_class=row["gf_source_type"], historical_loggf=row["loggf_adopted"],
            historical_gf_sigma_dex=row["gf_sigma_dex"],
            verified_source_record_id=chosen["source_record_id"] if chosen else "",
            bibliography_chain=chosen["bibliography_chain"] if chosen else "",
            verified_source_class=chosen["source_class"] if chosen else "UNRESOLVED",
            verified_loggf=chosen["loggf"] if chosen else np.nan,
            gf_uncertainty_bound_dex=chosen["gf_bound_dex"] if chosen else np.nan,
            uncertainty_basis=chosen["uncertainty_basis"] if chosen else "NO_DEFENSIBLE_BOUND",
            reference_membership="MEMBER" if chosen else "HOLD", disposition_reason=why,
            physical_identity_status=chosen["identity_status"] if chosen else "UNRESOLVED",
            component_semantics="UNRESOLVED_BLEND_TOTAL" if blend_total else "FINE_STRUCTURE_TOTAL" if chosen else "UNVERIFIED",
            hfs_component_status="PENDING_COMPONENT_LEDGER", codex_membership="PENDING_DEPTH_AND_COMPONENT_AUDIT",
            deep_membership="PENDING_DEPTH_AND_COMPONENT_AUDIT",
            external_membership=row["literature_line_set_membership"],
            source_candidate_count=len(candidates)))
    return pd.DataFrame(decisions), sources, pd.DataFrame(links)


def finalize(ledger: pd.DataFrame, sources: pd.DataFrame, components: pd.DataFrame):
    """Require independent raw upper/lower identity before a measurement handoff."""
    ledger = ledger.copy()
    src = sources.set_index("source_record_id")
    for i, r in ledger.iterrows():
        if not r.verified_source_record_id:
            continue
        evidence = src.loc[r.verified_source_record_id]
        part = components[components.canonical_line_id.eq(r.canonical_line_id)]
        status = "SOURCE_ONLY_NO_RAW_HOLDING_COMPONENTS"
        if len(part) and r.hfs_component_status == "RAW_PHYSICAL_COMPONENTS_RECOVERED":
            p = part.iloc[0]
            lower = evidence.get("lower_J", np.nan)
            upper = evidence.get("upper_J", np.nan)
            if r.species == "Al II" and evidence.source_key == "Johnson1986":
                lower, upper = 0., 1.
            if math.isfinite(lower) and math.isfinite(upper):
                status = "VERIFIED_RAW_J_AND_SOURCE_IDENTITY" if (part.lower_J.eq(lower).all() and part.upper_J.eq(upper).all()) else "SOURCE_RAW_J_CONFLICT"
                source_upper = evidence.get("upper_EP", np.nan)
                # NL2017 energies are printed to 0.001 eV; raw VALD to 0.0001.
                # Preserve this independent energy check when the source has it.
                if math.isfinite(source_upper) and not (part.raw_upper_EP-source_upper).abs().le(.0006).all():
                    status = "SOURCE_RAW_UPPER_ENERGY_CONFLICT"
            else:
                # NIST snapshot carries energies but no upper J. Do not invent it.
                status = "SOURCE_J_NOT_IN_SNAPSHOT"
        ledger.loc[i, "physical_identity_status"] = status
        if status in ("SOURCE_RAW_J_CONFLICT", "SOURCE_RAW_UPPER_ENERGY_CONFLICT"):
            ledger.loc[i, "reference_membership"] = "HOLD"
            ledger.loc[i, "disposition_reason"] = status
        ledger.loc[i, "atomic_handoff_status"] = "READY" if status == "VERIFIED_RAW_J_AND_SOURCE_IDENTITY" else "HOLD"
        # The Träbert successor is cited in the accessible source chain, but its
        # article is unavailable for direct review. Retain Johnson's measured
        # value and record that limitation instead of blocking the handoff.
        if evidence.source_key == "Johnson1986":
            ledger.loc[i, "disposition_reason"] = "JOHNSON1986_RETAINED_TRABERT1999_UNAVAILABLE"
    ledger["atomic_handoff_status"] = ledger.atomic_handoff_status.fillna("HOLD")
    ledger["replication_membership"] = ledger.external_membership.map(
        lambda value: "EXCLUDED_BY_SOURCE" if "AGSS21" in value and "EXCLUDED" in value else
        "MEMBER" if "AGSS21" in value else "NONMEMBER")
    return ledger


def load_pool_dispositions(line_set: str, root: Path = ROOT) -> pd.DataFrame:
    """Return the complete pool denominator, including HOLD/NONMEMBER rows.

    This is the RYA-1217 intake surface. It is not permission to run an engine;
    eligibility and exact-holding conditioning must also pass. External replication
    gf remain owned by reference_lineset.load('asplund-al'), not our chosen gf.
    """
    from pipeline.al_eligibility import POOL_AXES
    if line_set not in POOL_AXES:
        raise ValueError(f"unknown Al line set: {line_set}")
    ledger = pd.read_csv(root / VERIFIED / "line_dispositions.csv").fillna({"verified_source_record_id": ""})
    memberships = pd.read_csv(root / VERIFIED / "pool_memberships.csv")
    membership = memberships[memberships.line_set.eq(line_set)]
    if len(ledger) != 505 or not ledger.canonical_line_id.is_unique or len(membership) != len(ledger):
        raise ValueError("incomplete or duplicate Al disposition denominator")
    joined = ledger.merge(membership[["canonical_line_id", "line_set", "membership"]], on="canonical_line_id", validate="one_to_one")
    if len(joined) != len(ledger):
        raise ValueError("Al pool membership IDs do not reconcile")
    return joined


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    ledger, sources, links = build()
    ledger, components = component_audit(ledger, ROOT)
    ledger = finalize(ledger, sources, components)
    from pipeline.al_eligibility import matrices, POOL_AXES
    line_matrix, engine_matrix = matrices(ledger, ROOT)
    memberships = []
    for r in ledger.to_dict("records"):
        for pool, axis in POOL_AXES.items():
            memberships.append(dict(canonical_line_id=r["canonical_line_id"], line_set=pool,
                membership=r[axis], atomic_handoff_status=r["atomic_handoff_status"],
                verified_source_record_id=r["verified_source_record_id"],
                gf_semantics="EXTERNAL_REFERENCE_GF_REQUIRED" if pool == "asplund-al" else "VERIFIED_SELECTED_GF",
                reason=r["disposition_reason"]))
    grade_delta = ledger[["canonical_line_id", "historical_source_class", "historical_gf_grade",
        "verified_source_class", "reference_membership", "codex_membership", "deep_membership",
        "replication_membership", "disposition_reason"]]
    summary = dict(ticket="RYA-1134", denominator=len(ledger),
        historical_source_classes=ledger.historical_source_class.value_counts().to_dict(),
        verified_source_classes=ledger.verified_source_class.value_counts().to_dict(),
        memberships={pool: ledger[axis].value_counts().to_dict() for pool, axis in POOL_AXES.items()},
        atomic_handoff=ledger.atomic_handoff_status.value_counts().to_dict(),
        n_line_holding_cells=len(line_matrix), n_solar_holdings=line_matrix.holding_id.nunique(),
        n_engine_cells=len(engine_matrix), engine_dispositions=engine_matrix.status.value_counts().to_dict(),
        abundance_runs=0, measurement_gate="HOLD_RYA1176_AND_EXACT_HOLDING_INTAKE",
        scientific_freeze_status="HOLD_SOURCE_HIERARCHY_BIBLIOGRAPHY_AND_HOLDING_GATES",
        uncertainty_contract="source bounds retain their confidence semantics; no blanket/sentinel sigma",
        candidate_scope="frozen source assets and RYA1173 lineage; absence is bounded to inspected assets")
    frames = dict(line_dispositions=ledger, source_quantities=sources, source_candidates=links,
                  component_proof=components, pool_memberships=pd.DataFrame(memberships),
                  grade_delta=grade_delta, line_holding_matrix=line_matrix, engine_matrix=engine_matrix)
    frames["band_grade_eligibility"] = line_matrix.groupby(
        ["band", "verified_source_class", "status"], dropna=False).agg(
        n_line_holding_cells=("canonical_line_id", "size"),
        n_unique_lines=("canonical_line_id", "nunique")).reset_index()
    files = {name + ".csv": frame.round(10).to_csv(index=False) for name, frame in frames.items()}
    files["summary.json"] = json.dumps(summary, indent=2, sort_keys=True) + "\n"
    inputs = {str(INTAKE / "al_line_manifest.csv"), str(SOURCE / "burheim_table1_quantities.csv"),
              "data/linelists/linelist_solar.csv", "pipeline/al_grade_verification.py", "pipeline/al_eligibility.py",
              "data/catalog/holdings_manifest_registry.csv", "data/catalog/instrument_catalog.csv",
              "data/catalog/model_registry.csv", "scripts/measure_band_ew.py", "pipeline/gerber_nlte.py",
              "data/refs/bibliography.csv", "pipeline/model_registry.py", "pipeline/band_policy.py",
              "pipeline/line_match.py", "pipeline/wavelength_util.py",
              "config/constants.py", "data/nlte_grids/Al_Amarsi2020_PySME.csv",
              "data/nlte_grids/Al_Amarsi2020_PySME.prov.json", "scripts/derive_band_products.py",
              "data/reference/vujnovic2002_al/raw/ReadMe",
              "pipeline/telluric_policy.py", "scripts/rya1001_al_census.py", "scripts/line_accounting_rya709.py"}
    inputs.update(sources.source_path.dropna())
    inputs.update(str(p.relative_to(ROOT)) for p in (ROOT / SOURCE).iterdir() if p.is_file())
    from pipeline.telluric_policy import EVIDENCE
    inputs.add(str(EVIDENCE.relative_to(ROOT)))
    inputs.update(components.raw_path[components.raw_path.ne("")])
    files["input_hashes.json"] = json.dumps({p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(inputs)}, indent=2) + "\n"
    if args.check:
        changed = [name for name, value in files.items() if not (args.out/name).exists() or (args.out/name).read_text() != value]
        if changed:
            raise SystemExit("stale RYA-1134 outputs: " + ", ".join(changed))
        print("RYA-1134 generated artifacts reproduce exactly")
    else:
        args.out.mkdir(parents=True, exist_ok=False)
        for name, value in files.items():
            (args.out/name).write_text(value)
    print(json.dumps(ledger.verified_source_class.value_counts().to_dict(), indent=2))


if __name__ == "__main__":
    main()
