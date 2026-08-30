#!/usr/bin/env python3
"""Build the RYA-1132 frozen Al atomic-data intake assets.

Inventory only: this script does not mutate canonical_gf and does not derive an
abundance.  The RYA-1001 physical-feature census is the denominator; later
canonical, Burheim, IGRINS, and CRIRES+ evidence is overlaid conservatively.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import re

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data/audit/rya1132_al_intake"
CENSUS = ROOT / "data/results/rya1001/rya1001_al_line_census.csv"
CANONICAL = ROOT / "data/linelists/canonical_gf.csv"
BURHEIM = ROOT / "data/reference/al_gf_lab/al1_lab_loggf.csv"
IGRINS = ROOT / "data/literature/igrins_nandakumar_2024/igrins_al_completeness_audit.csv"
CHIAPPINO = ROOT / "data/audit/rya1059_chiappino/al_completeness_delta.csv"
VUJ_RAW = ROOT / "data/reference/vujnovic2002_al/raw"

# Manual physical adjudication after reading Vujnovic Table 2.  These are finite
# A-values with stated uncertainties and unique fine-structure identities.  Limits,
# ratio-only rows, 3092.839 (no independent uncertainty), and the two lines for which
# Burheim 2023 is the stronger source are deliberately absent.
VUJ_PROMOTE_A = {2652.484, 2660.393, 3082.153, 3092.710, 3944.006, 3961.520}

WEB_FOLLOWUP = [
    # species, wavelength, source, DOI, evidence, uncertainty, disposition, note
    ("Al I", 3944.006, "Vujnovic et al. 2002", "10.1051/0004-6361:20020560", "PRIMARY_LAB_COMPOSITE", "Aki 10%", "INGESTED_GF_LAB", "Finite Aki and unique fine-structure identity."),
    ("Al I", 3961.520, "Vujnovic et al. 2002", "10.1051/0004-6361:20020560", "PRIMARY_LAB_COMPOSITE", "Aki 6%", "INGESTED_GF_LAB", "Finite Aki and unique fine-structure identity."),
    ("Al I", 2652.484, "Vujnovic et al. 2002", "10.1051/0004-6361:20020560", "PRIMARY_LAB_COMPOSITE", "Aki 12%", "INGESTED_GF_LAB", "Measured intensity ratio; absolute Aki from selected published lifetime."),
    ("Al I", 2660.393, "Vujnovic et al. 2002", "10.1051/0004-6361:20020560", "PRIMARY_LAB_COMPOSITE", "Aki 11%", "INGESTED_GF_LAB", "Measured intensity ratio; absolute Aki from selected published lifetime."),
    ("Al I", 3082.153, "Vujnovic et al. 2002", "10.1051/0004-6361:20020560", "PRIMARY_LAB_COMPOSITE", "Aki 12%", "INGESTED_GF_LAB", "Measured branching/intensity ratio; absolute Aki from lifetime."),
    ("Al I", 3092.710, "Vujnovic et al. 2002", "10.1051/0004-6361:20020560", "PRIMARY_LAB_COMPOSITE", "Aki 2%", "INGESTED_GF_LAB", "Level-resolved separately from 3092.839."),
    ("Al I", 3092.839, "Vujnovic et al. 2002", "10.1051/0004-6361:20020560", "PRIMARY_LAB_COMPOSITE", "Aki not independently stated", "PHYSICAL_CROSSMATCH_REQUIRED", "Intensity-ratio component; do not inherit the 3092.710 uncertainty."),
    ("Al II", 2669.157, "Johnson, Smith & Parkinson 1986", "10.1086/164569", "PRIMARY_LAB", "Aki=(3.33+/-0.23)e3 s-1 at 90% confidence", "INGESTED_GF_LAB", "Direct time-resolved ion-storage lifetime measurement; one decay channel."),
    ("Al II", 1670.78861, "Murphy & Berengut 2014 / Griesmann & Kling 2000", "10.1093/mnras/stt2120", "MIXED", "wavelength 20 m/s; f=theory", "WAVELENGTH_ONLY_NOT_GF_LAB", "Excellent laboratory wavelength, but quoted oscillator strength is theoretical."),
    ("Al I;Al II", np.nan, "Kelleher & Podobedova 2008 NIST critical compilation", "10.1063/1.2734564", "CRITICALLY_EVALUATED", "per-line accuracy grade", "INGEST_EVALUATED_TIER_SEPARATELY", "Reference values only; evaluated is not Codex/Deep primary-lab grade under RYA-946."),
    ("Al I;Al II", np.nan, "Vujnovic et al. 2002 CDS J/A+A/388/704", "10.1051/0004-6361:20020560", "SOURCE_TABLE", "106 rows across tables 2-5", "ACQUIRED_AND_NORMALIZED", "Raw fixed-width tables preserved; limits and lifetime provenance retained."),
    ("Al I", np.nan, "IR literature follow-up", "", "NEGATIVE_RESULT", "", "NO_NEW_GRADING_SOURCE", "Buurman 1986, Davidson 1990, and Buurman & Donszelmann 1990 are already represented through Burheim/NIST; no new source found for held 11254.9, 7835/7836, 8772/8773, or 21208 A."),
]


def text(v: object) -> str:
    return "" if pd.isna(v) else str(v).strip()


def band(w: float) -> str:
    if w < 2000: return "FUV"
    if w < 3000: return "NUV"
    if w < 3780: return "near-UV"
    if w < 6910: return "VIS"
    if w < 9199: return "red-optical"
    if w < 13000: return "NIR"
    if 13195.23 <= w < 15007.11: return "J"
    if 15007.11 <= w < 17493.69: return "H"
    if 19510.4 <= w < 24857.7: return "K"
    return "OUTSIDE_CURRENT_INSTRUMENT_REACH"


def source_type(tier: str, source: str) -> str:
    t, s = tier.upper(), source.upper()
    if t == "LAB" or "BURHEIM" in s: return "PRIMARY_LABORATORY"
    if t.startswith("NIST") or "NIST" in s: return "CRITICALLY_EVALUATED"
    if "THEORY" in s or "P19" in s or "OP95" in s: return "THEORETICAL"
    return "FALLBACK"


def nearest(frame: pd.DataFrame, w: float, ep: float, wcol: str, epcol: str | None,
            wtol: float = .06, eptol: float = .02) -> pd.Series | None:
    if frame.empty: return None
    ok = (pd.to_numeric(frame[wcol], errors="coerce") - w).abs() <= wtol
    if epcol and np.isfinite(ep):
        ok &= (pd.to_numeric(frame[epcol], errors="coerce") - ep).abs() <= eptol
    c = frame[ok]
    if len(c) != 1: return None
    return c.iloc[0]


def _number(raw: str) -> float:
    return np.nan if not raw.strip() else float(raw.strip())


def _upper_j(level: str) -> float:
    matches = re.findall(r"_(\d+)(?:/(\d+))?_", level)
    if not matches:
        return np.nan
    n, d = matches[-1]
    return float(n) / float(d or 1)


#: 🔴 ARTIFACTS MUST NOT CARRY RAW FLOAT REPR — RYA-1084's lesson, applied at the write.
#: `loggf_adopted` and the Vujnovic derived columns come out of `math.log10`, and a
#: transcendental differs by ONE ULP between numpy builds. Written unrounded at full
#: 17-digit repr, that ULP becomes a visible byte change, so the reproducibility test
#: compared bit-exact floats across environments and failed on CI (numpy 2.5.1) against
#: artifacts committed from the Mac (numpy 2.2.6) -- 14 lines of 506, max drift 8.9e-16,
#: no value, row, column or ordering actually different.
#:
#: 10 decimals is five orders of magnitude coarser than the drift and far finer than
#: anything physical here (gf uncertainties run 0.01-0.1 dex), so it removes the noise
#: without touching a single meaningful digit.
FLOAT_DECIMALS = 10


def _stable(df: "pd.DataFrame") -> "pd.DataFrame":
    """The frame as written: float columns rounded so the bytes are portable."""
    return df.round(FLOAT_DECIMALS)


def load_vujnovic() -> pd.DataFrame:
    """Parse all four CDS J/A+A/388/704 fixed-width source tables."""
    rows: list[dict] = []
    specs = {
        "table2.dat": ("Al I", 2), "table3.dat": ("Al I", 3),
        "table4.dat": ("Al II", 4), "table5.dat": ("Al II", 5),
    }
    for filename, (species, table) in specs.items():
        for source_row, line in enumerate((VUJ_RAW / filename).read_text().splitlines(), 1):
            if not line.strip():
                continue
            if table == 2:
                upper, lower, wave = line[0:20].strip(), line[21:36].strip(), _number(line[37:46])
                intensity, intensity_unc = _number(line[48:52]), _number(line[53:55])
                branching, branching_unc = _number(line[61:66]), _number(line[66:68])
                aki_limit, aki, aki_unc = line[69:70].strip(), _number(line[70:75]), _number(line[77:79])
            elif table in (3, 4):
                upper, lower, wave = line[0:19].strip(), line[20:35].strip(), _number(line[36:44])
                intensity, intensity_unc = _number(line[45:49]), _number(line[50:52])
                branching = branching_unc = aki = aki_unc = np.nan
                aki_limit = ""
            else:
                upper, lower, wave = line[0:13].strip(), line[14:27].strip(), _number(line[28:36])
                intensity = intensity_unc = np.nan
                branching, branching_unc = _number(line[38:42]), _number(line[43:46])
                aki_limit, aki, aki_unc = line[47:48].strip(), _number(line[48:55]), _number(line[58:60])
            j_upper = _upper_j(upper)
            finite_aki = np.isfinite(aki) and not aki_limit
            loggf = (math.log10(1.49919e-16 * (2*j_upper + 1) * wave**2 * aki * 1e8)
                     if finite_aki and np.isfinite(j_upper) else np.nan)
            sigma = math.log10(1 + aki_unc/100) if finite_aki and np.isfinite(aki_unc) else np.nan
            rows.append({"source_row_id":f"vuj2002_t{table}_{source_row:03d}", "table":table,
                "species":species, "upper_level":upper, "lower_level":lower,
                "wavelength_A":wave, "intensity_ratio":intensity,
                "intensity_unc_pct":intensity_unc, "branching_ratio":branching,
                "branching_unc_pct":branching_unc, "aki_limit":aki_limit,
                "aki_1e8_s-1":aki, "aki_unc_pct":aki_unc, "upper_J":j_upper,
                "derived_loggf":loggf, "derived_sigma_dex":sigma,
                "doi":"10.1051/0004-6361:20020560"})
    return pd.DataFrame(rows)


def ingest_new_lab_sources(m: pd.DataFrame, out: Path) -> pd.DataFrame:
    """Overlay only physically adjudicated finite laboratory measurements."""
    vuj = load_vujnovic()
    cross_rows = []
    for _, src in vuj.iterrows():
        candidates = m[m.species.eq(src.species)].copy()
        candidates["delta_A"] = (candidates.wavelength_air - src.wavelength_A).abs()
        near = candidates[candidates.delta_A <= .08].sort_values("delta_A")
        matched = len(near) == 1
        target = near.iloc[0] if matched else None
        promotable = (matched and src.species == "Al I" and
                      any(abs(src.wavelength_A-w) < .0005 for w in VUJ_PROMOTE_A) and
                      np.isfinite(src.derived_loggf) and np.isfinite(src.derived_sigma_dex))
        disposition = "GF_LAB_PROMOTED" if promotable else (
            "MATCHED_NOT_PROMOTED" if matched else "NO_UNIQUE_MANIFEST_MATCH")
        cross_rows.append({**src.to_dict(), "canonical_line_id":
            (target.canonical_line_id if matched else ""), "wavelength_delta_A":
            (target.delta_A if matched else np.nan), "disposition":disposition})
        if not promotable:
            continue
        idx = target.name
        m.loc[idx, ["loggf_adopted","gf_source","gf_source_type","gf_grade",
                    "gf_sigma_dex","gf_source_doi","upper_lower_level_identity",
                    "intake_status","source_ticket"]] = [
            src.derived_loggf, "EXP-VUJNOVIC2002", "PRIMARY_LABORATORY", "GF-LAB",
            src.derived_sigma_dex, src.doi,
            f"{src.lower_level} - {src.upper_level}", "FROZEN",
            "RYA-1132;RYA-1001;Vujnovic2002"]
        m.loc[idx, "notes"] = ("Vujnovic 2002 finite laboratory Aki ingested; loggf "
            "derived reproducibly from wavelength, upper J, and Aki. No abundance adoption.")

    # Johnson et al. 1986: direct one-channel Al II 2669.157 measurement.
    q = m[m.species.eq("Al II")].copy()
    q["delta_A"] = (q.wavelength_air - 2669.157).abs()
    q = q[q.delta_A <= .08]
    if len(q) != 1:
        raise AssertionError("Al II 2669.157 must have one physical manifest match")
    idx = q.index[0]
    johnson_loggf = math.log10(1.49919e-16 * 3 * 2669.157**2 * 3.33e3)
    johnson_sigma = math.log10(1 + .23/3.33)  # published 90%-confidence bound, conservative
    m.loc[idx, ["loggf_adopted","gf_source","gf_source_type","gf_grade",
                "gf_sigma_dex","gf_source_doi","upper_lower_level_identity",
                "intake_status","source_ticket"]] = [
        johnson_loggf, "EXP-JOHNSON1986", "PRIMARY_LABORATORY", "GF-LAB",
        johnson_sigma, "10.1086/164569", "3s2 1S0 - 3s3p 3P1o", "FROZEN",
        "RYA-1132;RYA-1001;Johnson1986"]
    m.loc[idx, "notes"] = ("Johnson 1986 direct ion-storage Aki ingested; uncertainty "
        "stored conservatively as the published 90%-confidence logarithmic bound. "
        "Träbert 1999/NIST remains the higher-precision comparison.")

    _stable(vuj).to_csv(out / "vujnovic2002_normalized.csv", index=False)
    _stable(pd.DataFrame(cross_rows)).to_csv(out / "vujnovic2002_crossmatch.csv", index=False)
    return m


def build(out: Path = OUT) -> dict:
    out.mkdir(parents=True, exist_ok=True)
    c = pd.read_csv(CENSUS, low_memory=False)
    c = c[c.ion.isin(["I", "II"])].copy()  # Al III is outside this ticket's atomic scope.
    can = pd.read_csv(CANONICAL, low_memory=False)
    can = can[can.species.isin(["Al I", "Al II"])].copy()
    bur = pd.read_csv(BURHEIM)
    igr = pd.read_csv(IGRINS)
    chi = pd.read_csv(CHIAPPINO)

    rows: list[dict] = []
    for row_number, (_, r) in enumerate(c.iterrows()):
        w, ep = float(r.wave_air_A), float(r.ep_eV)
        species = f"Al {r.ion}"
        cm = nearest(can[can.species == species], w, ep, "wavelength_air_A", "excitation_potential_eV")
        canonical_id = text(cm.line_id) if cm is not None else ""
        # The manifest freezes the strongest defensible source, while retaining the
        # current canonical value separately so scale mismatches cannot disappear.
        adopted = float(r.best_log_gf)
        source = text(r.best_gf_source)
        tier = ("GF-LAB" if "BURHEIM" in source.upper() else
                (text(r.nist_grade) or text(r.tier)))
        sigma = r.best_sigma_dex
        doi = ("10.1051/0004-6361/202245394" if "BURHEIM" in source.upper()
               else (text(cm.gf_source_doi) if cm is not None else ""))
        if abs(w - 11254.924) <= .08 and int(r.hfs_n_components) > 1:
            # Burheim measured the strong component, while this census row is the
            # unresolved feature total.  It is evidence, not an adoptable total gf.
            adopted = float(cm.log_gf) if cm is not None else float(r.log_gf_linelist_sum)
            source = text(cm.loggf_reference) if cm is not None else "VALD3_BLEND_TOTAL"
            tier = text(cm.gf_tier) if cm is not None else "VALD3"
            sigma = cm.gf_sigma_dex if cm is not None else np.nan
            doi = text(cm.gf_source_doi) if cm is not None else ""
        memberships = ["RYA1001_PHASE0"]
        if ((igr.rya1001_wavelength_candidate_A - w).abs() <= .08).any(): memberships.append("NANDAKUMAR2024_IGRINS")
        if ((chi.wavelength_air_A - w).abs() <= .08).any(): memberships.append("CHIAPPINO2026_CRIRES+")
        problem = []
        if text(r.scale_mismatch).lower() == "true": problem.append("SCALE_MISMATCH")
        if float(r.central_depth) < .05: problem.append("BELOW_OBSERVABILITY")
        if float(r.central_depth) > .60: problem.append("SATURATED_DEEP")
        if bool(r.blend_flag): problem.append("BLEND_FLAG")
        if not canonical_id: problem.append("ABSENT_CANONICAL")
        b = band(w)
        context = "AVAILABLE" if text(r.band_methods) else "NO_DECLARED_BAND_POLICY"
        if b in {"J", "H", "K"} and not text(r.band_methods): context = "LINELIST_OR_ROUTE_NOT_WIRED"
        intake = "FROZEN" if canonical_id else "CROSSMATCH_REVIEW"
        suitability = "CANDIDATE_NOT_SELECTED"
        if .05 <= float(r.central_depth) <= .60 and canonical_id:
            suitability = "ELIGIBLE_WITH_STATED_GF_TIER"
        rows.append({
            "canonical_line_id": f"alphys_{r.ion}_{w:.4f}_{row_number:04d}",
            "canonical_source_line_id": canonical_id,
            "species": species, "wavelength_air": w, "wavelength_vac": r.wave_vac_A,
            "lower_EP": ep, "upper_lower_level_identity": text(r.burheim_transition),
            "band": b, "instrument_reach": text(r.instruments_coverage_module),
            "transition_source": "RYA-1001 physical-feature census",
            "loggf_adopted": adopted, "gf_source": source,
            "gf_source_type": source_type(tier, source), "gf_grade": tier,
            "gf_sigma_dex": sigma, "gf_source_doi": doi,
            "current_canonical_loggf": (cm.log_gf if cm is not None else np.nan),
            "current_canonical_source": (text(cm.loggf_reference) if cm is not None else ""),
            "competing_gf_summary": (f"Burheim={r.burheim_log_gf}; canonical={r.canonical_log_gf}; "
                                      f"NIST={r.nist_log_gf}"),
            "HFS_status": "COMPONENT_SUM_VERIFIED" if int(r.hfs_n_components) > 1 else "NO_SPLIT_COMPONENTS_IN_CENSUS",
            "component_or_total": "TOTAL_TRANSITION_GF",
            "literature_line_set_membership": "|".join(memberships),
            "telluric_risk": "VERIFICATION_REQUIRED" if bool(r.telluric_required_band) else "NOT_FLAGGED",
            "synthesis_context_status": context, "intake_status": intake,
            "measurement_suitability_status": suitability,
            "rejection_problem_code": "|".join(problem), "source_ticket": "RYA-1132;RYA-1001",
            "notes": "No abundance adoption; provenance ranking fixed before measurement.",
        })

    # Preserve new empirical IGRINS candidates absent from the older physical census.
    for _, r in igr[igr.rya1001_wavelength_candidate_A.isna()].iterrows():
        w = float(r.wavelength_air_A)
        rows.append({"canonical_line_id": f"alphys_I_{w:.4f}_igrins", "canonical_source_line_id":"", "species":"Al I",
            "wavelength_air":w, "wavelength_vac":np.nan, "lower_EP":np.nan,
            "upper_lower_level_identity":"", "band":band(w), "instrument_reach":"IGRINS",
            "transition_source":"Nandakumar et al. 2024 empirical line list", "loggf_adopted":np.nan,
            "gf_source":"UNRESOLVED", "gf_source_type":"FALLBACK", "gf_grade":"UNRESOLVED",
            "gf_sigma_dex":np.nan, "gf_source_doi":"", "current_canonical_loggf":np.nan,
            "current_canonical_source":"", "competing_gf_summary":"",
            "HFS_status":"HFS_IDENTITY_REQUIRED", "component_or_total":"UNKNOWN",
            "literature_line_set_membership":"NANDAKUMAR2024_IGRINS", "telluric_risk":"VERIFICATION_REQUIRED",
            "synthesis_context_status":"IDENTITY_AND_GF_NOT_RESOLVED", "intake_status":"CROSSMATCH_REVIEW",
            "measurement_suitability_status":"HOLD", "rejection_problem_code":"MISSING_PHYSICAL_IDENTITY|MISSING_GF_PROVENANCE",
            "source_ticket":"RYA-1132;RYA-1056", "notes":"Wavelength evidence only; promotion forbidden."})

    # Burheim's full 12-row table is retained, including four mid-IR controls outside
    # current reach.  These prove completeness and prevent a range from becoming a list.
    held_waves = np.array([float(x["wavelength_air"]) for x in rows])
    for j, r in bur.iterrows():
        w=float(r.wavelength_air_A)
        if len(held_waves) and np.min(np.abs(held_waves-w)) <= .08: continue
        rows.append({"canonical_line_id":f"alphys_I_{w:.4f}_burheim{j:02d}","canonical_source_line_id":"",
            "species":"Al I","wavelength_air":w,"wavelength_vac":r.lam_vac_A,"lower_EP":r.elo_eV,
            "upper_lower_level_identity":f"{r.lower_level} - {r.upper_level}","band":band(w),
            "instrument_reach":"OUTSIDE_CURRENT_REACH","transition_source":"Burheim2023 Table 3",
            "loggf_adopted":r.loggf,"gf_source":"EXP-BURHEIM23","gf_source_type":"PRIMARY_LABORATORY",
            "gf_grade":"GF-LAB","gf_sigma_dex":r.e_loggf_dex,"gf_source_doi":"10.1051/0004-6361/202245394",
            "current_canonical_loggf":np.nan,"current_canonical_source":"","competing_gf_summary":f"P19={r.loggf_papoulia19}; K95={r.loggf_kurucz95}; TOPbase={r.loggf_topbase00}",
            "HFS_status":"SOURCE_TOTAL_TRANSITION","component_or_total":"TOTAL_TRANSITION_GF",
            "literature_line_set_membership":"BURHEIM2023_TABLE3_COMPLETE_CONTROL","telluric_risk":"OUTSIDE_REACH",
            "synthesis_context_status":"OUTSIDE_CURRENT_INSTRUMENT_REACH","intake_status":"FROZEN_SOURCE_CONTROL",
            "measurement_suitability_status":"OUTSIDE_CURRENT_REACH","rejection_problem_code":"OUTSIDE_CURRENT_INSTRUMENT_REACH",
            "source_ticket":"RYA-1132;RYA-1002","notes":"Full-table completeness control; not a measurement candidate."})

    m = pd.DataFrame(rows).sort_values(["wavelength_air", "species"]).reset_index(drop=True)
    if m.canonical_line_id.duplicated().any():
        raise AssertionError("manifest IDs must be unique")
    m = ingest_new_lab_sources(m, out)
    _stable(m).to_csv(out / "al_line_manifest.csv", index=False)

    grade = m.groupby(["band","gf_source_type"]).size().unstack(fill_value=0).reset_index()
    grade.to_csv(out / "gf_grade_matrix.csv", index=False)
    cov = m.groupby(["band","instrument_reach","gf_source_type"]).size().unstack(fill_value=0).reset_index()
    cov.to_csv(out / "coverage_matrix.csv", index=False)

    conflicts = m[m.rejection_problem_code.str.contains("SCALE_MISMATCH|MISSING_PHYSICAL_IDENTITY", na=False)].copy()
    special = m[np.isclose(m.wavelength_air, 11254.925, atol=.08)].copy()
    special["rejection_problem_code"] = "BURHEIM_STRONG_COMPONENT_VS_CANONICAL_BLEND_TOTAL"
    special["notes"] = "Burheim +0.327 is the strong component; observed feature total is +0.354. Never substitute one for the other."
    conflicts = pd.concat([conflicts, special]).drop_duplicates("canonical_line_id", keep="last")
    _stable(conflicts).to_csv(out / "conflict_ledger.csv", index=False)
    _stable(m[m.rejection_problem_code.ne("")]).to_csv(out / "problem_child_ledger.csv", index=False)

    bib = pd.DataFrame([
        ["Burheim2023","Burheim, Hartman & Nilsson 2023, A&A 672 A197","10.1051/0004-6361/202245394","2023A&A...672A.197B","Table 3","PRIMARY_LAB_GF","https://www.aanda.org/articles/aa/full_html/2023/04/aa45394-22/aa45394-22.html","https://www.aanda.org/articles/aa/pdf/2023/04/aa45394-22.pdf"],
        ["Vujnovic2002","Vujnovic et al. 2002, A&A 388, 704-711","10.1051/0004-6361:20020560","2002A&A...388..704V","CDS J/A+A/388/704 tables 2-5","PRIMARY_LAB_COMPOSITE","https://www.aanda.org/articles/aa/abs/2002/23/aa7151/aa7151.html","https://www.aanda.org/articles/aa/pdf/2002/23/aa7151.pdf"],
        ["Trabert1999","Trabert et al. 1999, J. Phys. B 32, 537-552","10.1088/0953-4075/32/2/031","1999JPhB...32..537T","Al II 2669 intercombination lifetime","PRIMARY_LAB_GF","https://doi.org/10.1088/0953-4075/32/2/031",""],
        ["Johnson1986","Johnson, Smith & Parkinson 1986, ApJ 308, 1013","10.1086/164569","1986ApJ...308.1013J","Al II 2669 direct measurement","PRIMARY_LAB_GF","https://ntrs.nasa.gov/citations/19870032227","https://adsabs.harvard.edu/pdf/1986ApJ...308.1013J"],
        ["KelleherPodobedova2008","Kelleher & Podobedova 2008, J. Phys. Chem. Ref. Data 37, 709","10.1063/1.2734564","2008JPCRD..37..709K","Al I-Al XII compilation","CRITICALLY_EVALUATED","https://www.nist.gov/publications/atomic-transition-probabilities-aluminuma-critical-compilation","https://www.nist.gov/system/files/documents/srd/jpcrd372008911p.pdf"],
        ["Papoulia2019","Papoulia, Ekman & Jonsson 2019, A&A 621, A16","10.1051/0004-6361/201833764","2019A&A...621A..16P","CDS calculated transition tables","THEORETICAL_COMPARATOR","https://www.aanda.org/articles/aa/full_html/2019/01/aa33764-18/aa33764-18.html","https://arxiv.org/pdf/1808.09478"],
        ["GriesmannKling2000","Griesmann & Kling 2000, ApJ 536, L113-L115","10.1086/312738","2000ApJ...536L.113G","Al II 1670 laboratory wavelength","WAVELENGTH_ONLY_NOT_GF_LAB","https://www.nist.gov/publications/interferometric-measurement-resonance-transition-wavelengths-civ-siiv-aliii-al-ii-and","https://arxiv.org/pdf/astro-ph/0004190"],
        ["RoedererLawler2021","Roederer & Lawler 2021, ApJ 912, 119","10.3847/1538-4357/abf142","2021ApJ...912..119R","Al II 2669 stellar use and source chain","EMPIRICAL_AND_PROVENANCE_GUIDE","https://iopscience.iop.org/article/10.3847/1538-4357/abf142","https://arxiv.org/pdf/2103.12764"],
        ["Lind2022","Lind et al. 2022, A&A 665, A33","10.1051/0004-6361/202142195","2022A&A...665A..33L","Al line list and non-LTE context","MODELING_CONTEXT","https://www.aanda.org/articles/aa/full_html/2022/09/aa42195-21/aa42195-21.html","https://openresearch-repository.anu.edu.au/bitstreams/f1d056b4-19fb-466d-9ae7-b800f7b408fd/download"],
        ["JonssonLundberg1983","Jonsson & Lundberg 1983, Z. Phys. A 313, 151-154","10.1007/BF01417221","","Al I 2S and 2D lifetime sequences","PRIMARY_LAB_LIFETIME_PROVENANCE","https://portal.research.lu.se/en/publications/natural-radiative-lifetimes-in-the-2s12-and-2d5232-sequences-of-a/",""],
        ["Davidson1990","Davidson, Volten & Donszelmann 1990, A&A 238, 452-454","","1990A&A...238..452D","Al I nd 2D lifetimes and oscillator strengths","PRIMARY_LAB_PROVENANCE_ALREADY_PROPAGATED","https://www.researchgate.net/publication/234442756_Lifetimes_and_oscillator_strengths_of_the_3s2nd_2D_series_in_neutral_aluminum",""],
        ["NIST_ASD","Kramida et al., NIST Atomic Spectra Database","","NIST_ASD","Al transition export","CRITICALLY_EVALUATED","https://physics.nist.gov/asd",""],
        ["Nandakumar2024","Nandakumar et al. 2024 IGRINS abundance lines","10.3847/1538-4357/ad4451","2024ApJ...964...96N","CDS A15","EMPIRICAL_MEMBERSHIP_ONLY","https://doi.org/10.3847/1538-4357/ad4451",""],
        ["Chiappino2026","Chiappino et al. 2026 CRIRES+ J/H/K line set","10.3847/1538-4357/ae7de8","2026ApJ...","publisher tables","EMPIRICAL_MEMBERSHIP_ONLY","https://doi.org/10.3847/1538-4357/ae7de8",""],
        ["RYA1001","Codex Al Phase-0 VALD/manual physical-feature census","","","rya1001_al_line_census.csv","CANDIDATE_DENOMINATOR","",""] ,
    ], columns=["source_id","citation","doi","ads_bibcode","table_catalog","role","article_url","download_url"])
    bib.to_csv(out / "source_bibliography.csv", index=False)

    follow = pd.DataFrame(WEB_FOLLOWUP, columns=["species","wavelength_A","source","doi","evidence_class","published_uncertainty","disposition","note"])
    follow.to_csv(out / "web_source_followup.csv", index=False)

    verdicts = {
        "UV":"PARTIAL_GF_LAB_INGESTED_POLICY_BLOCKED", "VIS":"FROZEN_WITH_DOCUMENTED_FALLBACKS",
        "IR":"BLOCKED_PIPELINE_COVERAGE", "overall":"BLOCKED_PIPELINE_COVERAGE",
        "measurement_unblocked":False,
        "reasons":{
            "UV":"Six finite uncertainty-bearing Vujnovic Al I transitions and the direct Johnson Al II 2669 measurement are physically crossmatched and ingested as GF-LAB. Limits, ratio-only rows, and 3092.839 without an independent uncertainty remain non-promoted. FUV/NUV still lack a declared measurement policy.",
            "VIS":"Physical identities and adopted sources are explicit; primary-lab Burheim lines remain distinct from fallbacks.",
            "IR":"RYA-1003 telluric verification and RYA-1004 red-edge/context coverage remain open; wavelength-only empirical candidates stay HOLD.",
        }}
    (out / "intake_verdict.json").write_text(json.dumps(verdicts, indent=2)+"\n")
    summary = {"ticket":"RYA-1132", "candidate_rows":len(m),
        "by_ion":m.species.value_counts().sort_index().to_dict(),
        "by_band":m.band.value_counts().sort_index().to_dict(),
        "by_gf_source_type":m.gf_source_type.value_counts().to_dict(),
        "crossmatch_review":int((m.intake_status=="CROSSMATCH_REVIEW").sum()),
        "abundances_generated":False, "canonical_gf_mutated":False, "verdicts":verdicts}
    (out / "summary.json").write_text(json.dumps(summary, indent=2)+"\n")
    report = f"""# RYA-1132 - Al atomic-data intake closure

This is an inventory and provenance freeze, not an abundance result. It preserves
**{len(m)}** physical candidates ({int((m.species == 'Al I').sum())} Al I,
{int((m.species == 'Al II').sum())} Al II) from UV through IR. No canonical gf
or Solar abundance was changed.

## Verdict

- UV: `PARTIAL_GF_LAB_INGESTED_POLICY_BLOCKED` - all 106 Vujnovic CDS rows are
  preserved and normalized. Six finite, uncertainty-bearing Al I transitions plus
  Johnson et al. 1986 Al II 2669 are physically crossmatched and ingested as GF-LAB.
  Limits, ratio-only rows, and 3092.839 without an independent uncertainty remain
  non-promoted; FUV/NUV still lack a declared measurement policy.
- VIS: `FROZEN_WITH_DOCUMENTED_FALLBACKS` - physical identities and evidence
  ceilings are explicit; Burheim laboratory lines are not blurred with fallback gf.
- IR: `BLOCKED_PIPELINE_COVERAGE` - RYA-1003 telluric verification and RYA-1004
  synthesis red-edge/context coverage remain open.
- Overall: `BLOCKED_PIPELINE_COVERAGE`; no new Solar Al measurement is unblocked.

The 6696.015 Burheim transition remains physically distinct from 6696.185. The
11254.9 conflict ledger explicitly distinguishes Burheim's strong-component
`log gf=+0.327` from the observed blended-feature total near `+0.354`. IGRINS and
CRIRES+ wavelength-only evidence remains HOLD unless wavelength plus EP/levels
establish a unique physical transition.

The ingestion changes seven manifest rows from fallback to GF-LAB without mutating
`canonical_gf` or deriving an abundance. `vujnovic2002_normalized.csv` and
`vujnovic2002_crossmatch.csv` retain the complete source accounting.

The web follow-up found no new independent IR source beyond the older Buurman,
Davidson, and Buurman-Donszelmann measurements already propagated through
Burheim/NIST. In particular it did not resolve 11254.9, 7835/7836, 8772/8773,
or the 21208 A IGRINS candidate. See `web_source_followup.csv`.

## Reproduce

`python3 scripts/build_al_intake_rya1132.py`
"""
    (out / "README.md").write_text(report)

    if out.resolve() == OUT.resolve():
        ledger_path = ROOT / "data/audit/rya1129_atomic_intake/intake_status_ledger.csv"
        ledger = pd.read_csv(ledger_path)
        mask = ledger.element.eq("Al")
        if mask.sum() != 1: raise AssertionError("RYA-1129 ledger must contain exactly one Al row")
        ledger.loc[mask, "intake_status"] = verdicts["overall"]
        ledger.loc[mask, "primary_source_status"] = (
            f"PRIMARY_LAB={int((m.gf_source_type=='PRIMARY_LABORATORY').sum())};"
            f"EVALUATED_NIST={int((m.gf_source_type=='CRITICALLY_EVALUATED').sum())};"
            f"FALLBACK={int((m.gf_source_type=='FALLBACK').sum())}")
        ledger.loc[mask, "canonical_gf_status"] = f"{len(m)} Al I/II physical candidates frozen in RYA-1132 manifest"
        ledger.loc[mask, "HFS_status"] = "COMPONENT SUMS PRESERVED; 11254.9 COMPONENT/BLEND CONFLICT EXPLICIT"
        ledger.loc[mask, "band_coverage_status"] = "UV PARTIAL_GF_LAB_INGESTED_POLICY_BLOCKED; VIS FROZEN_WITH_DOCUMENTED_FALLBACKS; IR BLOCKED_PIPELINE_COVERAGE"
        ledger.loc[mask, "source_ticket"] = "RYA-1132"
        ledger.loc[mask, "last_verified"] = "2026-08-30"
        ledger.loc[mask, "notes"] = "Vujnovic CDS tables normalized; six Al I lines plus Johnson Al II 2669 ingested GF-LAB. Limits/ratio-only rows remain held. No new IR source or abundance."
        ledger.to_csv(ledger_path, index=False)
    return summary


def main() -> None:
    p=argparse.ArgumentParser(); p.add_argument("--out",type=Path,default=OUT); a=p.parse_args()
    print(json.dumps(build(a.out),indent=2))


if __name__ == "__main__": main()
