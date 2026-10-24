#!/usr/bin/env python3
"""RYA-1268 -- minimum litscan anchors (runbook step 4) for every remaining element.

Values are transcribed from the papers' recommended-abundance tables, never from memory:
  AGSS21   Asplund, Amarsi & Grevesse 2021, A&A 653, A141, Table 2 (local aa40445-21.pdf)
  LBP25    Lodders, Bergemann & Palme 2025, Space Sci. Rev., Table 2 (local AbundancesLodders.pdf)
  AG26     Amarsi & Grevesse 2026, "Solar chemical composition" (arXiv 2608.23155), Table 3
           (local Amarsi_Grevesse_2026_solar_composition_2608.23155.pdf) -- the AGSS21 successor.
ANCHOR (Ryan 2026-10-24, RYA-1251): AG26 where its value or sigma differs from AGSS21, AGSS21
everywhere else; the other is kept as corroborating. Published elements (Fe, Si, Al, C, N, O,
Mg) are not written by this script and keep their anchors.
The band is anchor +/- sigma, DEVIATE at 2 sigma -- the Si/Al rule. Each file is a precursor:
the element's own dossier ticket re-projects it when its literature sweep runs.
Existing files are never overwritten, except with --force and only a file this script wrote
(source_of_truth: RYA-1268).
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LIT = ROOT / "data" / "reference" / "litscan"
# element: (ions, AGSS21 value, sigma, AGSS21 basis, LBP25 value, sigma, flag)
T = {
    "Li": ("I", 0.96, 0.06, "Wang et al. 2021 3D-NLTE, Li I 670.8 nm", 1.04, 0.09, "B"),
    "Na": ("I", 6.22, 0.03, "AGSS21 3D-NLTE, Table A.1 (5 Na I)", 6.29, 0.05, "A"),
    "P": ("I", 5.41, 0.03, "Scott et al. 2015 (A&A 573, A25) 3D-LTE", 5.44, 0.12, "C+"),
    "S": ("I", 7.12, 0.03, "Scott et al. 2015 (A&A 573, A25) 3D-LTE + 1D-NLTE", 7.16, 0.11, "B-"),
    "K": ("I", 5.07, 0.03, "AGSS21 3D-NLTE, Table A.1 (3 K I)", 5.09, 0.09, "B"),
    "Ca": ("I, II", 6.30, 0.03, "AGSS21 3D-NLTE, Table A.1 (Ca I + Ca II)", 6.35, 0.06, "A-"),
    "Sc": ("I, II", 3.14, 0.04, "Scott et al. 2015 (A&A 573, A26) + Lawler et al. 2019 gf", 3.13, 0.11, "B-"),
    "Ti": ("II", 4.97, 0.05, "Scott et al. 2015 (A&A 573, A26), Ti II only", 4.97, 0.11, "B-"),
    "V": ("I", 3.90, 0.08, "Scott et al. 2015 (A&A 573, A26) + Lawler et al. 2014 gf", 3.89, 0.16, "C-"),
    "Cr": ("I, II", 5.62, 0.04, "Scott et al. 2015 (A&A 573, A26) + Lawler et al. 2017 Cr II gf", 5.74, 0.11, "B-"),
    "Mn": ("I", 5.42, 0.06, "Scott et al. 2015 (A&A 573, A26), sigma inflated", 5.52, 0.05, "A"),
    "Co": ("I", 4.94, 0.05, "Scott et al. 2015 (A&A 573, A26) + Lawler et al. 2015 gf", 4.95, 0.11, "B-"),
    "Ni": ("I", 6.20, 0.04, "Scott et al. 2015 (A&A 573, A26) 3D-LTE", 6.24, 0.06, "A-"),
    "Cu": ("I", 4.18, 0.05, "Grevesse et al. 2015 (A&A 573, A27) 3D-LTE", 4.24, 0.11, "B-"),
    "Zn": ("I", 4.56, 0.05, "Grevesse et al. 2015 (A&A 573, A27)", 4.55, 0.11, "B-"),
    "Sr": ("I, II", 2.83, 0.06, "Grevesse et al. 2015 (A&A 573, A27)", 2.93, 0.11, "B-"),
    "Y": ("II", 2.21, 0.05, "Grevesse et al. 2015 (A&A 573, A27)", 2.30, 0.06, "A-"),
    "Zr": ("I, II", 2.59, 0.04, "Grevesse et al. 2015 (A&A 573, A27)", 2.68, 0.11, "B-"),
    "Ba": ("II", 2.27, 0.05, "Gallagher et al. 2020 (A&A 634, A55) 3D-NLTE", 2.27, 0.06, "A-"),
    "Eu": ("II", 0.52, 0.04, "Grevesse et al. 2015 update of Lawler et al. 2001 (HM synthesis)", 0.57, 0.06, "A-"),
}
# Amarsi & Grevesse 2026 Table 3: (value, sigma, modelling note as printed)
AG26 = {"Li": (0.96, 0.06, "3N"), "Na": (6.22, 0.03, "3N"), "P": (5.35, 0.04, "3L+1N"),
        "S": (7.06, 0.04, "3N"), "K": (5.07, 0.03, "3N"), "Ca": (6.30, 0.03, "3N"),
        "Sc": (3.13, 0.04, "3L+1N"), "Ti": (4.97, 0.05, "3L, full weight to Ti II"),
        "V": (3.90, 0.08, "3L+1N"), "Cr": (5.65, 0.04, "3L, full weight to Cr II"),
        "Mn": (5.47, 0.08, "3N, mean of two studies"), "Co": (4.94, 0.07, "3L+1N"),
        "Ni": (6.20, 0.05, "3L"), "Cu": (4.18, 0.05, "3L+1N"), "Zn": (4.56, 0.05, "3L+1N"),
        "Sr": (2.84, 0.07, "3L+1N"), "Y": (2.29, 0.07, "3N"), "Zr": (2.59, 0.05, "3L, full weight to Zr II"),
        "Ba": (2.27, 0.05, "3N"), "Eu": (0.49, 0.05, "3N")}
# further per-element studies named by AG26 / Grevesse+2015 / Storm+2024 (transcribed from each paper)
EXTRA = {"Eu": [
    ("Lawler et al. 2001", "ApJ 563, 1075, Table 6 (doi:10.1086/323407)", 0.52, 0.04,
     "1D LTE (HM74) synthesis of disc-centre intensity, 14 Eu II lines weighted 1-3; AGSS21's base value"),
    ("Mucciarelli et al. 2008", "A&A 484, 841 (doi:10.1051/0004-6361:20079327)", 0.52, 0.02,
     "3D LTE (CO5BOLD) EWs, 5 Eu II lines (4129.7, 6049.5, 6645.1, 7077.1, 7217.6), flux + intensity"),
    ("Storm et al. 2024", "A&A 683, A200 (doi:10.1051/0004-6361/202348971)", 0.57, 0.06,
     "3D NLTE, Eu II 6645.1 only (blended by weak unknown lines); Lodders 2025's Eu source")]}


def main() -> None:
    force = "--force" in sys.argv
    for el, (ions, v, s, basis, lv, ls, flag) in T.items():
        f = LIT / f"{el}.yaml"
        if f.exists() and not (force and "source_of_truth: RYA-1268" in f.read_text()):
            print(f"keep {f.name} (exists)")
            continue
        av, asg, anote = AG26[el]
        extra = "".join(f"""  - source: {src}
    citation: {cite}
    value: {xv:.2f}
    sigma: {xs:.2f}
    role: corroborating
    method: {m}
""" for src, cite, xv, xs, m in EXTRA.get(el, []))
        ion_list = [i.strip() for i in ions.split(",")]
        agss = ("Asplund, Amarsi & Grevesse 2021",
                "A&A 653, A141, Table 2 (doi:10.1051/0004-6361/202140445)", v, s, basis,
                "Asplund et al. 2021 (A&A 653, A141)")
        ag26 = ("Amarsi & Grevesse 2026", "Solar chemical composition, arXiv 2608.23155, Table 3",
                av, asg, anote, "Amarsi & Grevesse 2026 (arXiv 2608.23155)")
        moved = (av, asg) != (v, s)
        anc, other = (ag26, agss) if moved else (agss, ag26)
        other_role = "corroborating (superseded by AG26)" if moved else "corroborating (AGSS21 successor; unchanged)"
        why = ("AG26 changed this element (AGSS21 %.2f +/- %.2f -> AG26 %.2f +/- %.2f)"
               % (v, s, av, asg)) if moved else "AG26 unchanged from AGSS21; AGSS21 kept"
        ion_list = [i.strip() for i in ions.split(",")]
        cv, cs = anc[2], anc[3]
        f.write_text(f"""# LITSCAN -- {el} (RYA-1268 minimum anchor; transcribed from AGSS21 Table 2, Lodders,
# Bergemann & Palme 2025 Table 2 and Amarsi & Grevesse 2026 Table 3). Precursor: the
# element's own literature sweep re-projects this file.
# Band = anchor +/- sigma, DEVIATE at 2 sigma (Si/Al rule). Anchor: {why} (RYA-1251).
element: {el}
ion: {ion_list[0]}
ions: [{", ".join(ion_list)}]
source_of_truth: RYA-1268
scale: as adopted by {"AG26" if moved else "AGSS21"} ({anc[4]})
entries:
  - source: {anc[0]}
    citation: {anc[1]}
    value: {cv:.2f}
    sigma: {cs:.2f}
    role: best external (anchor)
    method: {anc[4]}
  - source: Lodders, Bergemann & Palme 2025
    citation: Space Science Reviews (accepted 2025-02-07), Table 2; local AbundancesLodders.pdf
    value: {lv:.2f}
    sigma: {ls:.2f}
    role: corroborating
    note: Lodders quality flag {flag}.
  - source: {other[0]}
    citation: {other[1]}
    value: {other[2]:.2f}
    sigma: {other[3]:.2f}
    role: {other_role}
    method: {other[4]}
{extra}range:
  central: {cv:.2f}
  min: {cv - cs:.2f}
  max: {cv + cs:.2f}
  best_external: {anc[5]}
  sigma_external: {cs:.2f}
  deviate_beyond: {2 * cs:.2f}
  derivation: best-external-plus-sigma
  is_statistical_spread: false
verification_owed:
  - Minimum anchor only; the element's literature sweep adds per-band references.
""")
        print(f"wrote {f.name}")


if __name__ == "__main__":
    main()
