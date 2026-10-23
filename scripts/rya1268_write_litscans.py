#!/usr/bin/env python3
"""RYA-1268 -- minimum litscan anchors (runbook step 4) for every remaining element.

Values are transcribed from the papers' recommended-abundance tables, never from memory:
  AGSS21   Asplund, Amarsi & Grevesse 2021, A&A 653, A141, Table 2 (local aa40445-21.pdf)
  LBP25    Lodders, Bergemann & Palme 2025, Space Sci. Rev., Table 2 (local AbundancesLodders.pdf)
The band is AGSS21 +/- sigma, DEVIATE at 2 sigma -- the Si/Al rule. Each file is a precursor:
the element's own dossier ticket re-projects it when its literature sweep runs.
Existing files are never overwritten.
"""
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


def main() -> None:
    for el, (ions, v, s, basis, lv, ls, flag) in T.items():
        f = LIT / f"{el}.yaml"
        if f.exists():
            print(f"keep {f.name} (exists)")
            continue
        ion_list = [i.strip() for i in ions.split(",")]
        f.write_text(f"""# LITSCAN -- {el} (RYA-1268 minimum anchor; transcribed from AGSS21 Table 2 and
# Lodders, Bergemann & Palme 2025 Table 2). Precursor: the element's own literature sweep
# re-projects this file. Band = AGSS21 +/- sigma, DEVIATE at 2 sigma (Si/Al rule).
element: {el}
ion: {ion_list[0]}
ions: [{", ".join(ion_list)}]
source_of_truth: RYA-1268
scale: as adopted by AGSS21 ({basis})
entries:
  - source: Asplund, Amarsi & Grevesse 2021
    citation: A&A 653, A141, Table 2 (doi:10.1051/0004-6361/202140445)
    value: {v:.2f}
    sigma: {s:.2f}
    role: best external (anchor)
    method: {basis}
  - source: Lodders, Bergemann & Palme 2025
    citation: Space Science Reviews (accepted 2025-02-07), Table 2; local AbundancesLodders.pdf
    value: {lv:.2f}
    sigma: {ls:.2f}
    role: corroborating
    note: Lodders quality flag {flag}.
range:
  central: {v:.2f}
  min: {v - s:.2f}
  max: {v + s:.2f}
  best_external: Asplund et al. 2021 (A&A 653, A141)
  sigma_external: {s:.2f}
  deviate_beyond: {2 * s:.2f}
  derivation: best-external-plus-sigma
  is_statistical_spread: false
verification_owed:
  - Minimum anchor only; the element's literature sweep adds per-band references.
""")
        print(f"wrote {f.name}")


if __name__ == "__main__":
    main()
