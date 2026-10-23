#!/usr/bin/env python3
"""RYA-1268 -- Reference Grade line coverage for every remaining element, per band.

Reference Grade = the Asplund (AGSS21) line set PLUS any lab-referenced line (Ryan
2026-10-21). This script builds both halves for every TARGET_ELEMENTS element not yet done
(Fe, C, N, O, Si, Al done; Mg is RYA-1250's) and says, per element x species x band, what
exists and what is missing:

  1. agss21_reference_lines.csv -- the Asplund set: AGSS21 Table A.1 for Na/K/Ca, otherwise
     the Scott+2015a/b / Grevesse+2015 table AGSS21 adopts (Sect. 2.4-2.6), with the
     AGSS21 changes to it recorded per element (AGSS21_NOTES). Each line carries a per-line
     gf sigma ONLY where the NIST ASD class belongs to the same gf: NIST's log gf within
     SAME_GF_DEX of the published one. A NIST class for a different gf is not this gf's error.
  2. coverage_matrix.csv -- per band: Asplund lines / with sigma, the NIST-graded lab pool
     (class A..C, in canonical_gf, measurable in the Sun), Elgueta 2026 Sun-robust IR lines,
     and a verdict: REFERENCE (>= MIN_LINES lines), THIN (1..MIN_LINES-1), or EMPTY.

Nothing here touches canonical_gf, runs a measurement, or writes SOURCES.yaml.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from config.constants import TARGET_ELEMENTS  # noqa: E402
from pipeline.gf_grades import GRADED_NIST_CLASSES, nist_sigma_dex  # noqa: E402

REF = ROOT / "data" / "reference"
SETS = REF / "agss21_cited_sets"
OUT = ROOT / "data" / "audit" / "rya1268_line_coverage"
PRIMARY = ROOT / "data" / "linelists" / "primary_gf"
DONE = {"Fe", "C", "N", "O", "Si", "Al", "Mg"}
REMAINING = [e for e in TARGET_ELEMENTS if e not in DONE] + ["Zn"]  # Zn: 28th canonical
BANDS = [("near-UV", 3000.0, 3780.0), ("VIS", 3780.0, 6910.0), ("red-optical", 6910.0, 9199.0),
         ("NIR", 9199.0, 13000.0), ("H", 15007.11, 17493.69), ("K", 19452.42, 24845.62)]
MIN_LINES = 3          # a band needs >= 3 Reference lines to stand without Codex/Deep help
SAME_GF_DEX = 0.03     # NIST's class applies to the published gf only if it IS that gf
DEPTH = (0.01, 0.75)   # measurable in the Sun: detectable, not saturated (linelist_solar)
JOIN_TOL_A = 0.06      # printed solar wavelength vs VALD; EP within 0.02 eV decides
TABLE_A1 = {"Na", "K", "Ca"}
SOURCE_PAPER = {"P": "scott2015_NatoCa", "S": "scott2015_NatoCa",
                **{e: "scott2015_FePeak" for e in ("Sc", "Ti", "V", "Cr", "Mn", "Co", "Ni")},
                **{e: "grevesse2015_Heavies" for e in ("Cu", "Zn", "Sr", "Y", "Zr")}}
LAB = ROOT / "data" / "reference" / "lab_gf_tables" / "lab_gf_lines.csv"
VALD_REFS = OUT / "vald_solar_gf_refs.csv"
# AGSS21's replacement gf for Sc (Lawler+2019) is not on VizieR; VALD carries it as LHSNWC.
VALD_LAB_CODE = {"Sc": "LHSNWC"}
# A gf error the adopting paper STATES in its text (no per-line table): species -> (sigma, quote).
# Upper end of a stated range; a species-wide statement applies to that species' set lines only.
STATED_IN_TEXT = {
    "Zr II": (0.0414, "Grevesse+2015 Sect. Zr: 'uncertainties of the mean gf-values we adopt are "
                      "of order 5-10%' (10% -> log10(1.10))"),
    ("Ca II", 7323.89): (0.03, "Scott+2015 Sect. Ca, forbidden line: 'the uncertainty should be of "
                               "order 0.03 dex'"),
}
# What AGSS21 (Sect. 2.2-2.6) changed relative to the table it adopts. gf_status
# UPDATE_PENDING = AGSS21 replaced the printed gf with a newer lab source; attach_sigma swaps in
# that source's gf where we hold it (lab_gf_lines.csv, or VALD's code for Sc).
AGSS21_NOTES = {
    "Li": ("AS_PRINTED", "Wang et al. 2021: Li I 670.8 nm doublet, Yan+1998 gf (NIST AAA) -- extra_solar_sets.csv"),
    "Na": ("AS_PRINTED", "Table A.1: same five lines and gf as Scott+2015b"),
    "K": ("AS_PRINTED", "Table A.1: 404.4/580.2 dropped, NIST (Kramida 2019) gf adopted"),
    "Ca": ("AS_PRINTED", "Table A.1: 616.13 dropped; 616.90+616.95 fitted as one feature"),
    "P": ("AS_PRINTED", "Scott+2015b 3D LTE adopted unchanged (Berzinsh+1997 gf)"),
    "S": ("AS_PRINTED", "Scott+2015b adopted unchanged"),
    "Sc": ("UPDATE_PENDING", "gf -> Lawler et al. 2019 (ApJS 241, 21) for all lines"),
    "Ti": ("AS_PRINTED", "value from Ti II ONLY (Wood+2013 gf); Ti I (Lawler+2013) kept as lab lines"),
    "V": ("UPDATE_PENDING", "gf -> Lawler et al. 2014 (V I) + Wood et al. 2014a (V II)"),
    "Cr": ("PARTIAL_UPDATE", "Cr II gf -> Lawler et al. 2017; Cr I as printed"),
    "Mn": ("AS_PRINTED", "Scott+2015a advocated (not Bergemann+2019); sigma inflated to 0.06"),
    "Co": ("UPDATE_PENDING", "gf -> Lawler et al. 2015 for lines in common"),
    "Ni": ("AS_PRINTED", "Scott+2015a 3D LTE, Wood+2014b gf"),
    "Cu": ("AS_PRINTED", "Grevesse+2015 retained (old Kock & Richter gf)"),
    "Zn": ("AS_PRINTED", "Grevesse+2015 adopted"),
    "Sr": ("AS_PRINTED", "Grevesse+2015 adopted"),
    "Y": ("AS_PRINTED", "Grevesse+2015 adopted"),
    "Zr": ("AS_PRINTED", "Grevesse+2015 adopted"),
    "Ba": ("AS_PRINTED", "Gallagher+2020 Table 1, 4 Ba II with per-line gf error (AGSS21's Ba source)"),
    "Eu": ("AS_PRINTED", "AGSS21/Grevesse print no Eu table; Lodders 2025 -> Storm+2024: Eu II 6645.104"),
}
UPDATED_IONS = {"Cr": {"Cr II"}}


def band_of(w: float) -> str:
    for name, lo, hi in BANDS:
        if lo <= w < hi:
            return name
    return "outside"


def load_nist() -> pd.DataFrame:
    parts = []
    for f in sorted(PRIMARY.glob("nist_asd_*_*_3000_25000.tsv")):  # RYA-1268's El_Ion pulls
        d = pd.read_csv(f, sep="\t")
        el, ion = f.name.split("_")[2:4]
        d["species"] = f"{el} {ion}"
        parts.append(d)
    n = pd.concat(parts, ignore_index=True)
    n = n[n.log_gf.notna() & n.wavelength_A.notna()].copy()
    n["graded"] = n.nist_grade.isin(GRADED_NIST_CLASSES)
    return n


def asplund_set() -> pd.DataFrame:
    cited = pd.read_csv(SETS / "agss21_cited_lines.csv")
    a1 = pd.read_csv(SETS / "agss21_tableA1_NaKCa.csv")
    rf = pd.read_csv(SETS / "agss21_cited_gf_refs.csv")
    # "\\citet{Lawler89}" -> "Lawler89": the paper's own bib key, readable and greppable
    refs = {(r.paper, r.gf_ref): re.sub(r"\\cite[tp]?\{([^}]*)\}", r"\1",
                                        re.sub(r"\\(mbox|nopagebreak)", "", r.reference))
            for r in rf.itertuples()}
    rows = []
    for _, r in a1.iterrows():
        rows.append({"species": r.species, "wavelength_air_A": r.wavelength_air_A,
                     "ep_eV": r.ep_eV, "published_loggf": r.published_loggf, "ew_pm": r.ew_pm,
                     "weight": np.nan, "paper": "agss21_tableA1", "doi": r.doi,
                     "citation": r.citation, "gf_source_per_line": "AGSS21 Table A.1",
                     "note": r.note if isinstance(r.note, str) else ""})
    for el, paper in SOURCE_PAPER.items():
        c = cited[(cited.paper == paper) & (cited.species.str.split().str[0] == el)]
        for _, r in c.iterrows():
            rows.append({"species": r.species, "wavelength_air_A": r.wavelength_air_A,
                         "ep_eV": r.ep_eV, "published_loggf": r.published_loggf,
                         "ew_pm": r.ew_pm, "weight": r.weight, "paper": paper, "doi": r.doi,
                         "citation": r.citation,
                         "gf_source_per_line": refs.get((paper, r.gf_ref), ""),
                         "note": ("gf weight reduced (*)" if r.gf_weight_reduced else "")
                         + (r.blend_note if isinstance(r.blend_note, str) else "")})
    x = pd.read_csv(SETS / "extra_solar_sets.csv", comment="#")
    for _, r in x.iterrows():
        rows.append({"set_name": r.set_name, "species": r.species,
                     "wavelength_air_A": r.wavelength_air_A, "ep_eV": r.ep_eV,
                     "published_loggf": r.published_loggf, "ew_pm": np.nan, "weight": np.nan,
                     "paper": r.paper, "doi": r.doi, "citation": r.citation,
                     "gf_source_per_line": r.gf_source_per_line, "stated_sigma_dex": r.loggf_sigma_dex,
                     "note": r.source_table})
    s = pd.DataFrame(rows)
    s["element"] = s.species.str.split().str[0]
    s["set_name"] = s.get("set_name").fillna(s.element.str.upper() + "_AGSS21")
    s["band"] = s.wavelength_air_A.map(band_of)
    st = {e: v[0] for e, v in AGSS21_NOTES.items()}
    s["agss21_gf_status"] = [
        ("UPDATE_PENDING" if sp in UPDATED_IONS.get(e, set()) else "AS_PRINTED")
        if st.get(e) == "PARTIAL_UPDATE" else st.get(e, "AS_PRINTED")
        for e, sp in zip(s.element, s.species)]
    s["agss21_note"] = s.element.map(lambda e: AGSS21_NOTES.get(e, ("", ""))[1])
    return s


def _nearest(df: pd.DataFrame, sp: str, w: float, ep: float, gf: float, wtol: float,
             wcol: str = "wavelength_air_A", ecol: str = "ep_eV", gcol: str = "loggf"):
    c = df[(df.species == sp) & ((df[wcol] - w).abs() <= wtol) & ((df[ecol] - ep).abs() <= 0.02)]
    # the closest in gf: a multiplet's components share EP and sit within the window
    return None if c.empty else c.loc[(c[gcol] - gf).abs().idxmin()]


def attach_sigma(s: pd.DataFrame, nist: pd.DataFrame) -> pd.DataFrame:
    """Per-line gf sigma, best source first: (1) the set's own stated error; (2) the lab paper's
    A +/- dA for the same transition (lab_gf_lines.csv) -- which is ALSO where AGSS21's
    replacement gf comes from for UPDATE_PENDING lines; (3) NIST's class on the same gf."""
    lab = pd.read_csv(LAB)
    vald = pd.read_csv(VALD_REFS)
    vald = vald[vald.gf_ref_code.isin(set(VALD_LAB_CODE.values()))]
    out = []
    for _, r in s.iterrows():
        tol = 0.1 if r.paper == "agss21_tableA1" else JOIN_TOL_A
        rec = {"printed_loggf": r.published_loggf, "published_loggf": r.published_loggf,
               "agss21_gf_status": r.agss21_gf_status, "loggf_sigma_dex": np.nan,
               "gf_sigma_basis": "", "lab_paper": "", "lab_loggf": np.nan}
        gf = r.published_loggf
        h = _nearest(lab, r.species, r.wavelength_air_A, r.ep_eV, gf, tol)
        if h is not None:
            rec.update(lab_paper=h.paper, lab_loggf=h.loggf)
        if r.agss21_gf_status == "UPDATE_PENDING":
            if h is not None:
                gf = float(h.loggf)
                rec.update(published_loggf=gf, agss21_gf_status=f"UPDATED_FROM_LAB ({h.paper})",
                           gf_source_per_line=f"{h.paper} (as AGSS21 adopted)")
            elif r.element in VALD_LAB_CODE:
                v = _nearest(vald, r.species, r.wavelength_air_A, r.ep_eV, gf, tol)
                if v is not None:
                    gf = float(v.loggf)
                    rec.update(published_loggf=gf,
                               agss21_gf_status=f"UPDATED_FROM_VALD ({VALD_LAB_CODE[r.element]})",
                               gf_source_per_line=f"{VALD_LAB_CODE[r.element]} via VALD (Lawler+2019, as AGSS21 adopted)")
            else:
                rec["agss21_gf_status"] = "PRINTED (no lab row for this line)"
        if pd.notna(r.get("stated_sigma_dex", np.nan)):
            rec.update(loggf_sigma_dex=r.stated_sigma_dex, gf_sigma_basis="stated by the source table")
        elif h is not None and abs(float(h.loggf) - gf) <= SAME_GF_DEX:
            rec.update(loggf_sigma_dex=h.loggf_sigma_dex,
                       gf_sigma_basis=f"lab A +/- dA, {h.paper} (same gf, |d|={abs(float(h.loggf) - gf):.3f})")
        n = _nearest(nist, r.species, r.wavelength_air_A, r.ep_eV, gf, tol, "wavelength_A", "ei_eV", "log_gf")
        rec.update(nist_loggf=np.nan if n is None else round(n.log_gf, 3),
                   nist_grade="" if n is None or pd.isna(n.nist_grade) else str(n.nist_grade),
                   d_loggf_nist=np.nan if n is None else round(float(n.log_gf) - gf, 3))
        stated = STATED_IN_TEXT.get((r.species, r.wavelength_air_A)) or STATED_IN_TEXT.get(r.species)
        if rec["gf_sigma_basis"]:
            pass
        elif stated:
            rec.update(loggf_sigma_dex=stated[0], gf_sigma_basis=f"stated in text: {stated[1]}")
        elif n is None:
            rec["gf_sigma_basis"] = "no lab row; transition not in NIST ASD"
        elif n.nist_grade in GRADED_NIST_CLASSES and abs(n.log_gf - gf) <= SAME_GF_DEX:
            rec.update(loggf_sigma_dex=round(nist_sigma_dex(n.nist_grade), 4),
                       gf_sigma_basis=f"NIST ASD class {n.nist_grade} on the same gf (|d|={abs(n.log_gf - gf):.3f})")
        elif n.nist_grade in GRADED_NIST_CLASSES:
            rec["gf_sigma_basis"] = f"NIST {n.nist_grade} is for a DIFFERENT gf (d={n.log_gf - gf:+.3f})"
        else:
            rec["gf_sigma_basis"] = f"no lab row; NIST class '{rec['nist_grade'] or 'none'}' is not a grade"
        out.append(rec)
    o = pd.DataFrame(out)
    # a swapped-in gf brings its own source; every other line keeps the paper's
    src = o.pop("gf_source_per_line") if "gf_source_per_line" in o else pd.Series(np.nan, index=o.index)
    base = s.drop(columns=[c for c in o.columns if c in s.columns] + ["stated_sigma_dex"],
                  errors="ignore").reset_index(drop=True)
    base["gf_source_per_line"] = src.fillna(base["gf_source_per_line"])
    return pd.concat([base, o], axis=1)


def lab_pool(nist: pd.DataFrame) -> pd.DataFrame:
    sol = pd.read_csv(ROOT / "data" / "linelists" / "linelist_solar.csv", low_memory=False)
    sol["species"] = sol.element.astype(str) + " " + sol.ion.astype(str)
    g = nist[nist.graded].copy()
    keep = []
    for sp, grp in g.groupby("species"):
        ss = sol[sol.species == sp].sort_values("wavelength_air_A")
        if ss.empty:
            continue
        w = ss.wavelength_air_A.to_numpy()
        i = np.clip(np.searchsorted(w, grp.wavelength_A.to_numpy()), 1, len(w) - 1)
        left, right = w[i - 1], w[i]
        j = np.where(np.abs(grp.wavelength_A.to_numpy() - left)
                     <= np.abs(grp.wavelength_A.to_numpy() - right), i - 1, i)
        dist = np.abs(grp.wavelength_A.to_numpy() - w[j])
        depth = ss.central_depth.to_numpy()[j]
        grp = grp.assign(solar_depth=np.where(dist <= 0.05, depth, np.nan))
        keep.append(grp)
    p = pd.concat(keep, ignore_index=True)
    p["measurable"] = p.solar_depth.between(*DEPTH)
    p["band"] = p.wavelength_A.map(band_of)
    p["element"] = p.species.str.split().str[0]
    return p


def elgueta() -> pd.DataFrame:
    """Elgueta 2026 Sun/G-dwarf Rob=Y from RYA-1058's normalised atomic tables -- the same
    selection SI_ELGUETA2026 registers (42 Si I). The paper's table5.dat 'Robust (YYYY)' list
    is a stricter subset (30 Si I); both counts are reported."""
    n = pd.read_csv(ROOT / "data" / "audit" / "rya1058_elgueta" / "normalized_lines.csv")
    n = n[n.sun_Gd_rob == "Y"]
    t5 = set()
    for line in (REF / "elgueta2026_vizier" / "table5.dat").read_text().splitlines():
        t5.add((line[2:6].strip(), round(float(line[7:16]), 1)))  # Lambda in 0.1 nm = A
    sp = n.species.str.extract(r"([A-Z][a-z]?)([IV]+)")
    e = pd.DataFrame({"species": sp[0] + " " + sp[1], "wavelength_A": n.wavelength_A,
                      "band": n.wavelength_A.map(band_of),
                      "in_table5": [(s, round(w, 1)) in t5 for s, w in zip(n.species, n.wavelength_A)]})
    e["element"] = e.species.str.split().str[0]
    return e


def canonical_join(s: pd.DataFrame) -> pd.DataFrame:
    """The set's OWN join to canonical_gf: species, wavelength within JOIN_TOL_A (0.1 A for a
    Table A.1 0.01-nm print) AND lower level within 0.02 eV -- the FE_AGSS21 rule (0.05 A +
    0.02 eV). Scott/Grevesse print solar-atlas wavelengths up to 0.055 A from VALD's (Cr I
    5785.024 vs 5784.969, same EP), so 0.012 A lost 51 real joins. The EP is the
    discriminator: a wavelength-only window finds two rows for V I 5727.655 and Cr I 5300.743
    (an unrelated high-EP Kurucz line). No unique hit -> blank, never a guess."""
    canon = pd.read_csv(ROOT / "data" / "linelists" / "canonical_gf.csv", low_memory=False)
    ids, why = [], []
    for _, r in s.iterrows():
        tol = 0.1 if (r.paper == "agss21_tableA1" and r.wavelength_air_A == round(r.wavelength_air_A, 1)) else JOIN_TOL_A
        h = canon[(canon.species == r.species)
                  & ((canon.wavelength_air_A - r.wavelength_air_A).abs() <= tol)
                  & ((canon.excitation_potential_eV - r.ep_eV).abs() <= 0.02)]
        ids.append(h.line_id.iloc[0] if len(h) == 1 else "")
        why.append("" if len(h) == 1 else ("MISSING from canonical_gf" if len(h) == 0
                                          else f"AMBIGUOUS ({len(h)} rows)"))
    return s.assign(canonical_line_id=ids, canonical_join_note=why)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    nist = load_nist()
    s = canonical_join(attach_sigma(asplund_set(), nist))
    s.to_csv(SETS / "agss21_reference_lines.csv", index=False)
    pool = lab_pool(nist)
    pool[pool.measurable & pool.element.isin(REMAINING)].to_csv(OUT / "nist_lab_pool_measurable.csv", index=False)
    el = elgueta()
    have_nist = set(nist.species)
    rows = []
    for e in REMAINING:
        sps = sorted(set(s[s.element == e].species) | set(pool[pool.element == e].species)
                     | set(el[el.element == e].species), key=lambda x: (len(x), x))
        for sp in sps:
            for band, _, _ in BANDS:
                a = s[(s.species == sp) & (s.band == band)]
                p = pool[(pool.species == sp) & (pool.band == band) & pool.measurable]
                g = el[(el.species == sp) & (el.band == band)]
                n_ref = len(a) + len(p[~p.wavelength_A.round(0).isin(a.wavelength_air_A.round(0))])
                verdict = ("REFERENCE" if n_ref >= MIN_LINES else "THIN" if n_ref else "EMPTY")
                # Elgueta's IR lines carry VALD3 lab gf but no per-line sigma: candidates that
                # need a graded gf before they count, as SI_ELGUETA2026's did (RYA-1233)
                n_ir = len(g[~g.wavelength_A.round(0).isin(a.wavelength_air_A.round(0))])
                rows.append({"element": e, "species": sp, "band": band,
                             "asplund_lines": len(a), "asplund_with_sigma": int(a.loggf_sigma_dex.notna().sum()),
                             "asplund_gf_status": ";".join(sorted(set(a.agss21_gf_status))) if len(a) else "",
                             "nist_lab_measurable": len(p),
                             "nist_lab_measurable_AB": int(p.nist_grade.isin(["AAA", "AA", "A+", "A", "B+", "B"]).sum()),
                             "elgueta_sun_robust": len(g),
                             "elgueta_in_table5": int(g.in_table5.sum()),
                             "reference_lines_total": n_ref, "verdict": verdict,
                             "with_elgueta_candidates": n_ref + n_ir,
                             "nist_pulled": sp in have_nist})
    m = pd.DataFrame(rows)
    m.to_csv(OUT / "coverage_matrix.csv", index=False)
    print(m.pivot_table(index=["element", "species"], columns="band", values="reference_lines_total",
                        aggfunc="sum").reindex(columns=[b[0] for b in BANDS]).fillna(0).astype(int).to_string())


if __name__ == "__main__":
    main()
