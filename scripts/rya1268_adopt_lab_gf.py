#!/usr/bin/env python3
"""RYA-1268 -- adopt the lab-referenced gf of every remaining element into canonical_gf: the
"any lab-referenced line" half of Reference Grade (Ryan 2026-10-21).

Order of precedence per canonical row (species, wavelength <= WTOL_A, lower level <= 0.02 eV,
unique):
  1. LAB      -- a lab paper's own machine-readable table (lab_gf_lines.csv): its gf and its
                 A +/- dA error. gf_tier LAB, i.e. the orchestrator's `reference` pool.
  2. LAB      -- a VALD line whose gf code is a lab paper we cannot hold (LHSNWC Sc, LGWSC Ti I,
                 LWHS Eu II) AND that NIST ASD grades A..C on the same gf (|d| <= SAME_GF_DEX):
                 the lab gf, priced by NIST's class.
  3. NIST-C+  -- any other NIST ASD line of class A..C: NIST's gf and class (the RYA-1250 Mg
                 adjudication convention). Graded (Codex/Deep pools), not `reference`.
Guards (each refusal is counted in the summary):
  * UPPER LEVEL -- the candidate's upper energy must equal (<= 0.02 eV) that of a VALD line at
    the canonical row's wavelength and EP. Wavelength + lower EP alone matched unrelated weak
    Kurucz lines that share a lower level (V I 7936.845: -6.07 -> -1.83).
  * NIST A-CONSISTENCY -- NIST's log gf must equal log10(1.4992e-16 lambda^2 g_k A) from its own
    A-value (<= 0.05 dex, g_k from VALD's J_up): astroquery returned fik=16 for Zn I 11054.25
    (log gf +1.68; its A gives -0.36, VALD -0.30). A corrupt row is refused, never repaired.
  * no VALD counterpart -> the gf may not move more than NO_VALD_MAX_DEX.
A row that is already graded, or that another ticket adjudicated, is never touched. Every change is logged with the old values in
data/audit/rya1268_line_coverage/canonical_adoption.csv. Culls are not this script's business.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from pipeline.gf_grades import GRADED_NIST_CLASSES, is_gf_graded, nist_sigma_dex  # noqa: E402

CANON = ROOT / "data" / "linelists" / "canonical_gf.csv"
LAB = ROOT / "data" / "reference" / "lab_gf_tables" / "lab_gf_lines.csv"
AUD = ROOT / "data" / "audit" / "rya1268_line_coverage"
PRIMARY = ROOT / "data" / "linelists" / "primary_gf"
WTOL_A = 0.03           # lab / NIST / VALD wavelengths share the same energy levels
SAME_GF_DEX = 0.03
UP_TOL_EV = 0.02
A_TOL_DEX = 0.05
NO_VALD_MAX_DEX = 0.5
# Ryan 2026-10-23: "super bad dex should not be included". A gf that would move a line by more
# than this is HELD -- not adopted -- and listed for review in held_large_gf_moves.csv. Blends
# are the synthesis' job; a gross gf disagreement is not.
HOLD_MOVE_DEX = 0.5
# A row another ticket already adjudicated (nist / adjudicated_kr1968_rya466 / isotope_rya1075 /
# adjudicated_rya468 / flagged_manual_pull_rya354 ...) is finished work: never overwritten here.
# Only these workflow states are open to adoption.
OPEN_STATUS = {"", "single_source", "pending", "agreed"}
VALD_LAB = {"LHSNWC": ("lawler2019_sc", "Lawler et al. 2019 (ApJS 241, 21)", "10.3847/1538-4365/ab08ef"),
            "LGWSC": ("lawler2013", "Lawler et al. 2013 (ApJS 205, 11)", "10.1088/0067-0049/205/2/11"),
            "LWHS": ("lawler2001_eu2", "Lawler et al. 2001 (ApJ 563, 1075)", "10.1086/323407")}
CITE = {"wood2013_ti2": "Wood et al. 2013 (ApJS 208, 27)", "lawler2014_v1": "Lawler et al. 2014 (ApJS 215, 20)",
        "wood2014_v2": "Wood et al. 2014 (ApJS 214, 18)", "lawler2017_cr2": "Lawler et al. 2017 (ApJS 228, 10)",
        "lawler2015": "Lawler et al. 2015 (ApJS 220, 13)", "wood2014": "Wood et al. 2014 (ApJS 211, 20)",
        "denhartog2011": "Den Hartog et al. 2011 (ApJS 194, 35)", "sobeck2007": "Sobeck et al. 2007 (ApJ 667, 1267)",
        "lawler2001_eu2": "Lawler et al. 2001 (ApJ 563, 1075)",
        # parsed from the paper's own PDF (rya1268_parse_lab_pdfs.py); Ca I keys carry ":<ExptRef>"
        "bw2006_ti1": "Blackwell-Whitehead et al. 2006 (MNRAS 373, 1603)",
        "bw2011_mn1": "Blackwell-Whitehead et al. 2011 (A&A 525, A44)",
        "denhartog2021_ca1": "Den Hartog et al. 2021 (ApJS 255, 27)",
        "lawler2018_co2": "Lawler et al. 2018 (ApJS 238, 7)"}


def nist() -> pd.DataFrame:
    parts = []
    for f in sorted(PRIMARY.glob("nist_asd_*_*_3000_25000.tsv")):
        d = pd.read_csv(f, sep="\t")
        el, ion = f.name.split("_")[2:4]
        parts.append(d.assign(species=f"{el} {ion}"))
    n = pd.concat(parts, ignore_index=True)
    n = n.rename(columns={"aki_s-1": "aki"})
    return n[n.log_gf.notna() & n.wavelength_A.notna() & n.nist_grade.isin(GRADED_NIST_CLASSES)]


def match(cg: pd.DataFrame, by_sp: dict, sp: str, w: float, ep: float):
    c = by_sp.get(sp)
    if c is None:
        return None
    h = c[((c.wavelength_air_A - w).abs() <= WTOL_A) & ((c.excitation_potential_eV - ep).abs() <= 0.02)]
    return h.index[0] if len(h) == 1 else None


def nist_upper(raw: str) -> float:
    try:
        return float(str(raw).split("-")[-1].strip().strip("[]()?+x "))
    except ValueError:
        return np.nan


def main() -> None:
    cg = pd.read_csv(CANON, low_memory=False)
    species = sorted(set(pd.read_csv(LAB).species) | set(nist().species))
    graded0 = is_gf_graded(cg) | ~cg.adjudication_status.fillna("").isin(OPEN_STATUS)
    by_sp = {sp: cg[(cg.species == sp) & ~graded0] for sp in species}
    log, taken = [], set()
    vall = pd.read_csv(AUD / "vald_solar_gf_refs.csv")
    vall = vall[vall.species.isin(species)]
    vby = {sp: g for sp, g in vall.groupby("species")}
    refused = {"upper_level": 0, "nist_A_inconsistent": 0, "no_vald_big_move": 0}

    def vald_levels(i):
        g = vby.get(cg.at[i, "species"])
        if g is None:
            return g
        return g[((g.wavelength_air_A - cg.at[i, "wavelength_air_A"]).abs() <= 0.01)
                 & ((g.ep_eV - cg.at[i, "excitation_potential_eV"]).abs() <= 0.02)]

    held = []

    def ok(i, gf, e_up, a_hz=None, w=None, src=""):
        if abs(gf - float(cg.at[i, "log_gf"])) > HOLD_MOVE_DEX:
            held.append({"line_id": cg.at[i, "line_id"], "species": cg.at[i, "species"],
                         "wavelength_air_A": cg.at[i, "wavelength_air_A"],
                         "ep_eV": cg.at[i, "excitation_potential_eV"],
                         "canonical_loggf": cg.at[i, "log_gf"], "canonical_tier": cg.at[i, "gf_tier"],
                         "candidate_loggf": gf, "candidate_source": src,
                         "d_dex": round(gf - float(cg.at[i, "log_gf"]), 3)})
            return False
        lv = vald_levels(i)
        if lv is None or lv.empty:
            if abs(gf - float(cg.at[i, "log_gf"])) > NO_VALD_MAX_DEX:
                refused["no_vald_big_move"] += 1
                return False
            return True
        hit = lv[(lv.ep_upper_eV - e_up).abs() <= UP_TOL_EV]
        if hit.empty:
            refused["upper_level"] += 1
            return False
        if a_hz is not None:
            gk = 2 * float(hit.j_upper.iloc[0]) + 1
            if abs(np.log10(1.4992e-16 * w * w * gk * a_hz) - gf) > A_TOL_DEX:
                refused["nist_A_inconsistent"] += 1
                return False
        return True

    def write(i, gf, sig, tier, ref, doi, tag, status, grade=None):
        log.append({"line_id": cg.at[i, "line_id"], "species": cg.at[i, "species"],
                    "wavelength_air_A": cg.at[i, "wavelength_air_A"], "old_loggf": cg.at[i, "log_gf"],
                    "old_tier": cg.at[i, "gf_tier"], "old_reference": cg.at[i, "loggf_reference"],
                    "new_loggf": gf, "new_sigma_dex": sig, "new_tier": tier, "new_reference": ref,
                    "nist_grade": grade or "", "status": status})
        cg.at[i, "log_gf"] = gf
        cg.at[i, "gf_sigma_dex"] = sig
        cg.at[i, "gf_tier"] = tier
        cg.at[i, "loggf_reference"] = ref
        cg.at[i, "gf_source_doi"] = doi
        cg.at[i, "lab_source_tag"] = tag
        cg.at[i, "adjudication_status"] = status
        if grade:
            cg.at[i, "nist_grade"] = grade
        taken.add(i)

    # 1. lab tables
    # Scott+2015's Ti I precedence is Lawler+2013 > Nitz+1998 > Blackwell-Whitehead+2006: BW2006
    # yields wherever VALD carries one of the first two (LGWSC / NWL) -- at 5252.100 they differ
    # by 0.27 dex and AGSS21 adopted Lawler's.
    vr = pd.read_csv(AUD / "vald_solar_gf_refs.csv")
    outranked = {"bw2006_ti1": vr[(vr.species == "Ti I") & vr.gf_ref_code.isin({"LGWSC", "NWL"})]}
    for r in pd.read_csv(LAB).itertuples():
        if r.paper in outranked and (
                (outranked[r.paper].wavelength_air_A - r.wavelength_air_A).abs() <= 0.06).any():
            continue
        i = match(cg, by_sp, r.species, r.wavelength_air_A, r.ep_eV)
        if i is None or i in taken or not ok(i, r.loggf, r.ep_upper_eV, src=r.paper):
            continue
        key, _, sub = r.paper.partition(":")
        where = (f"machine-readable table {r.vizier}" if isinstance(r.vizier, str) and r.vizier
                 else "table parsed from the paper's PDF") + (f", experiment {sub}" if sub else "")
        write(i, r.loggf, r.loggf_sigma_dex, "LAB", f"{CITE[key]}, {where} (stated uncertainty)",
              r.doi, r.paper,
              "adopted_rya1268_lab_table")
    n = nist()
    # 2. VALD lab code + NIST class on the same gf
    v = pd.read_csv(AUD / "vald_solar_gf_refs.csv")
    v = v[v.gf_ref_code.isin(VALD_LAB)]
    for r in v.itertuples():
        i = match(cg, by_sp, r.species, r.wavelength_air_A, r.ep_eV)
        if i is None or i in taken:
            continue
        c = n[(n.species == r.species) & ((n.wavelength_A - r.wavelength_air_A).abs() <= WTOL_A)
              & ((n.ei_eV - r.ep_eV).abs() <= 0.02) & ((n.log_gf - r.loggf).abs() <= SAME_GF_DEX)]
        if len(c) != 1 or not ok(i, r.loggf, r.ep_upper_eV, src=f"VALD {r.gf_ref_code}"):
            continue
        g = str(c.nist_grade.iloc[0])
        bib, cite, doi = VALD_LAB[r.gf_ref_code]
        write(i, r.loggf, round(nist_sigma_dex(g), 4), "LAB",
              f"{cite} via VALD ({r.gf_ref_code}); sigma: NIST ASD class {g} on the same gf",
              doi, bib, "adopted_rya1268_vald_lab_nist_sigma", grade=g)
    # 2b. an Asplund-cited line whose adopting paper STATES its gf error (Gallagher+2020 Table 1,
    #     Storm+2024 Eu II 12% / Y II 10%, Grevesse+2015 Zr II 5-10%, Scott+2015 [Ca II]): the
    #     published gf at that sigma. Zr II (Ljung+2006) and Eu II (Lawler+2001) have no free
    #     machine-readable table, so this is the only route that prices them.
    ref = pd.read_csv(ROOT / "data" / "reference" / "agss21_cited_sets" / "agss21_reference_lines.csv")
    ref = ref[ref.gf_sigma_basis.fillna("").str.startswith("stated") & ref.loggf_sigma_dex.notna()
              & ref.canonical_line_id.notna()]
    # a line two sets share takes the AGSS21 set's own route (Y II: Grevesse over Storm+2024)
    claim = pd.read_csv(ROOT / "data" / "reference" / "agss21_cited_sets" / "agss21_reference_lines.csv")
    claim = claim.groupby("canonical_line_id").set_name.nunique()
    ref = ref[ref.canonical_line_id.map(claim).eq(1)]
    lid = pd.Series(cg.index, index=cg.line_id)
    for r in ref.itertuples():
        i = lid.get(r.canonical_line_id)
        if i is None or i in taken or graded0.iloc[i]:
            continue
        # the audit already joined this line on lambda + EP; the published table prints no upper
        # level, so only the hold guard applies (e_up is VALD's own, the closest-gf component)
        lv = vald_levels(i)
        e_up = (float(lv.ep_upper_eV.iloc[(lv.loggf - r.published_loggf).abs().argmin()])
                if lv is not None and len(lv) else np.nan)
        if not ok(i, float(r.published_loggf), e_up, src=f"{r.set_name} stated sigma"):
            continue
        write(i, float(r.published_loggf), float(r.loggf_sigma_dex), "LAB",
              f"{r.citation}; gf {r.gf_source_per_line}; sigma {r.gf_sigma_basis}", r.doi, r.paper,
              "adopted_rya1268_published_stated_sigma")
    # 3. NIST A..C
    for r in n.itertuples():
        i = match(cg, by_sp, r.species, r.wavelength_A, r.ei_eV)
        if (i is None or i in taken or pd.isna(r.aki)
                or not ok(i, float(r.log_gf), nist_upper(r.ei_ek_raw), float(r.aki), r.wavelength_A,
                                       src=f"NIST {r.nist_grade}")):
            continue
        write(i, round(float(r.log_gf), 3), round(nist_sigma_dex(r.nist_grade), 4), "NIST-C+",
              f"NIST ASD (accuracy {r.nist_grade}; TP ref {r.ref_transition_probability})", None,
              "NIST_ASD", "adopted_rya1268_nist_asd", grade=r.nist_grade)
    cg.to_csv(CANON, index=False)
    a = pd.DataFrame(log)
    a.to_csv(AUD / "canonical_adoption.csv", index=False)
    print(a.groupby(["species", "new_tier"]).size().unstack(fill_value=0).to_string())
    pd.DataFrame(held).drop_duplicates("line_id").to_csv(AUD / "held_large_gf_moves.csv", index=False)
    print("refused:", refused, "| held (>", HOLD_MOVE_DEX, "dex):", len({h["line_id"] for h in held}))
    print(len(a), "rows adopted;", int((a.old_loggf - a.new_loggf).abs().gt(0.1).sum()), "moved > 0.1 dex")


if __name__ == "__main__":
    main()
