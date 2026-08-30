#!/usr/bin/env python3
"""Build the RYA-1132 frozen Al atomic-data intake assets.

Inventory only: this script does not mutate canonical_gf and does not derive an
abundance.  The RYA-1001 physical-feature census is the denominator; later
canonical, Burheim, IGRINS, and CRIRES+ evidence is overlaid conservatively.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data/audit/rya1132_al_intake"
CENSUS = ROOT / "data/results/rya1001/rya1001_al_line_census.csv"
CANONICAL = ROOT / "data/linelists/canonical_gf.csv"
BURHEIM = ROOT / "data/reference/al_gf_lab/al1_lab_loggf.csv"
IGRINS = ROOT / "data/literature/igrins_nandakumar_2024/igrins_al_completeness_audit.csv"
CHIAPPINO = ROOT / "data/audit/rya1059_chiappino/al_completeness_delta.csv"

WEB_FOLLOWUP = [
    # species, wavelength, source, DOI, evidence, uncertainty, disposition, note
    ("Al I", 2652.484, "Vujnovic et al. 2002", "10.1051/0004-6361:20020554", "PRIMARY_LAB_COMPOSITE", "Aki 12%", "PHYSICAL_CROSSMATCH_REQUIRED", "Measured intensity ratio; absolute Aki from selected published lifetime."),
    ("Al I", 2660.393, "Vujnovic et al. 2002", "10.1051/0004-6361:20020554", "PRIMARY_LAB_COMPOSITE", "Aki 11%", "PHYSICAL_CROSSMATCH_REQUIRED", "Measured intensity ratio; absolute Aki from selected published lifetime."),
    ("Al I", 3082.153, "Vujnovic et al. 2002", "10.1051/0004-6361:20020554", "PRIMARY_LAB_COMPOSITE", "Aki 12%", "PHYSICAL_CROSSMATCH_REQUIRED", "Measured branching/intensity ratio; absolute Aki from lifetime."),
    ("Al I", 3092.710, "Vujnovic et al. 2002", "10.1051/0004-6361:20020554", "PRIMARY_LAB_COMPOSITE", "Aki 2%", "PHYSICAL_CROSSMATCH_REQUIRED", "The 3092.710/3092.839 pair must be level-resolved, never proximity-selected."),
    ("Al I", 3092.839, "Vujnovic et al. 2002", "10.1051/0004-6361:20020554", "PRIMARY_LAB_COMPOSITE", "Aki not independently stated", "PHYSICAL_CROSSMATCH_REQUIRED", "Intensity-ratio component; do not inherit the 3092.710 uncertainty."),
    ("Al II", 2669.157, "Johnson, Smith & Parkinson 1986", "10.1086/164569", "PRIMARY_LAB", "Aki=(3.33+/-0.23)e3 s-1 at 90% confidence", "PHYSICAL_CROSSMATCH_REQUIRED", "Direct time-resolved ion-storage lifetime measurement; one decay channel."),
    ("Al II", 1670.78861, "Murphy & Berengut 2014 / Griesmann & Kling 2000", "10.1093/mnras/stt2120", "MIXED", "wavelength 20 m/s; f=theory", "WAVELENGTH_ONLY_NOT_GF_LAB", "Excellent laboratory wavelength, but quoted oscillator strength is theoretical."),
    ("Al I;Al II", np.nan, "Kelleher & Podobedova 2008 NIST critical compilation", "10.1063/1.2734564", "CRITICALLY_EVALUATED", "per-line accuracy grade", "INGEST_EVALUATED_TIER_SEPARATELY", "Reference values only; evaluated is not Codex/Deep primary-lab grade under RYA-946."),
    ("Al I;Al II", np.nan, "Vujnovic et al. 2002 CDS J/A+A/388/704", "10.1051/0004-6361:20020554", "SOURCE_TABLE", "29 Al I + 31 Al II absolute-probability rows", "ACQUIRE_AND_NORMALIZE", "Machine-readable tables 2-5; preserve limits and lifetime provenance."),
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
    m.to_csv(out / "al_line_manifest.csv", index=False)

    grade = m.groupby(["band","gf_source_type"]).size().unstack(fill_value=0).reset_index()
    grade.to_csv(out / "gf_grade_matrix.csv", index=False)
    cov = m.groupby(["band","instrument_reach","gf_source_type"]).size().unstack(fill_value=0).reset_index()
    cov.to_csv(out / "coverage_matrix.csv", index=False)

    conflicts = m[m.rejection_problem_code.str.contains("SCALE_MISMATCH|MISSING_PHYSICAL_IDENTITY", na=False)].copy()
    special = m[np.isclose(m.wavelength_air, 11254.925, atol=.08)].copy()
    special["rejection_problem_code"] = "BURHEIM_STRONG_COMPONENT_VS_CANONICAL_BLEND_TOTAL"
    special["notes"] = "Burheim +0.327 is the strong component; observed feature total is +0.354. Never substitute one for the other."
    conflicts = pd.concat([conflicts, special]).drop_duplicates("canonical_line_id", keep="last")
    conflicts.to_csv(out / "conflict_ledger.csv", index=False)
    m[m.rejection_problem_code.ne("")].to_csv(out / "problem_child_ledger.csv", index=False)

    bib = pd.DataFrame([
        ["Burheim2023","Burheim, Hartman & Nilsson 2023, A&A 672 A197","10.1051/0004-6361/202245394","2023A&A...672A.197B","Table 3","PRIMARY_LAB_GF","https://www.aanda.org/articles/aa/full_html/2023/04/aa45394-22/aa45394-22.html","https://www.aanda.org/articles/aa/pdf/2023/04/aa45394-22.pdf"],
        ["Vujnovic2002","Vujnovic et al. 2002, A&A 388, 704-711","10.1051/0004-6361:20020554","2002A&A...388..704V","CDS J/A+A/388/704 tables 2-5","PRIMARY_LAB_COMPOSITE","https://www.aanda.org/articles/aa/abs/2002/23/aa7151/aa7151.html","https://www.aanda.org/articles/aa/pdf/2002/23/aa7151.pdf"],
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
        "UV":"CROSSMATCH_REVIEW", "VIS":"FROZEN_WITH_DOCUMENTED_FALLBACKS",
        "IR":"BLOCKED_PIPELINE_COVERAGE", "overall":"BLOCKED_PIPELINE_COVERAGE",
        "measurement_unblocked":False,
        "reasons":{
            "UV":"Web follow-up found Vujnovic 2002 laboratory transition probabilities and the direct Johnson 1986 Al II 2669 measurement; tables require physical crossmatch, uncertainty normalization, and ingestion. FUV/NUV still lack a declared band policy.",
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

- UV: `CROSSMATCH_REVIEW` - Vujnovic et al. 2002 and Johnson et al. 1986 provide
  missed laboratory evidence for several manifest lines. Their tables still need
  level-resolved crossmatching and uncertainty normalization; FUV/NUV also lack
  a declared measurement policy.
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
        ledger.loc[mask, "band_coverage_status"] = "UV CROSSMATCH_REVIEW; VIS FROZEN_WITH_DOCUMENTED_FALLBACKS; IR BLOCKED_PIPELINE_COVERAGE"
        ledger.loc[mask, "source_ticket"] = "RYA-1132"
        ledger.loc[mask, "last_verified"] = "2026-08-30"
        ledger.loc[mask, "notes"] = "Web follow-up found UV laboratory candidates (Vujnovic 2002; Johnson 1986); physical crossmatch/ingest owed. No new IR source found. No abundance generated."
        ledger.to_csv(ledger_path, index=False)
    return summary


def main() -> None:
    p=argparse.ArgumentParser(); p.add_argument("--out",type=Path,default=OUT); a=p.parse_args()
    print(json.dumps(build(a.out),indent=2))


if __name__ == "__main__": main()
