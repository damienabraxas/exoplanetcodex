#!/usr/bin/env python3
"""Vendor the primary-laboratory Fe II log gf table — RYA-953.

🔴 WHY Fe II WAS UNGRADED WITH A LAB TABLE SITTING RIGHT THERE. `gf_rung.decide` gates
on `LAB_GRADED_SPECIES`, which RYA-1002 derived from `gf_grades.LAB_TABLES` so the two
could never disagree. Fe II had no entry, so every Fe II pool returned rung 1 with the
message "no primary-laboratory gf table exists for Fe II" -- while 22 Den Hartog 2019
LAB rows sat in `canonical_gf`, ingested by RYA-945 with their DOI and a per-line sigma.
The message was true about the TABLE and false about the DATA. Same shape as the Al
case RYA-1002 fixed: a species gained a lab table and the ladder did not notice.

WHERE THESE VALUES COME FROM, precisely. This does NOT re-fetch from the journal. It
extracts the rows RYA-945 already ingested and adjudicated (`adjudication_status ==
'lab_rya945'`, `gf_tier == 'LAB'`, `lab_source_tag == 'DH19'`), each carrying
DOI 10.3847/1538-4365/ab322e. The ingest is the primary-source step; this is the
projection of it into the shape `gf_grades` reads. Deriving it rather than re-keying it
is deliberate -- a second hand-entry of the same 22 numbers is a second place to make a
typo, and RYA-353's single-sourcing rule says the copy that is not the source drifts.

⚠️ WHAT THIS DOES NOT MAKE GRADEABLE. 22 lines, 3002.6-4583.8 A. Only 9 fall in VIS
(4200-6910) and every one of those sits below 4584 A -- the blue edge. There is no Fe II
laboratory gf anywhere in the red half of VIS, none in red-optical, and none in the IR.
A Fe II graded VIS product is therefore a BLUE product wearing a band name, and that has
to be said out loud rather than discovered from the wavelength column (RYA-489).
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

CANONICAL = ROOT / "data" / "linelists" / "canonical_gf.csv"
OUT = ROOT / "data" / "reference" / "fe2_gf_lab" / "fe2_lab_loggf.csv"

#: The one source in this table. A second tag appearing here without a CITATIONS entry
#: would produce an unattributable pedigree, so it is asserted rather than assumed.
EXPECTED_TAG = "DH19"
EXPECTED_DOI = "10.3847/1538-4365/ab322e"


def build() -> pd.DataFrame:
    cg = pd.read_csv(CANONICAL, comment="#", low_memory=False)
    fe2 = cg[cg.species.astype(str).str.strip() == "Fe II"]
    lab = fe2[fe2.gf_tier.astype(str).str.strip() == "LAB"].copy()
    if lab.empty:
        raise SystemExit(
            "no LAB-tier Fe II rows in canonical_gf. RYA-945 is the ingest that puts "
            "them there; this script only projects them. Run that first.")

    tags = sorted(set(lab.lab_source_tag.astype(str)))
    if tags != [EXPECTED_TAG]:
        raise SystemExit(
            f"Fe II LAB rows carry sources {tags}, expected only ['{EXPECTED_TAG}']. "
            f"Every source in a lab table needs a `gf_grades.CITATIONS` entry -- an "
            f"uncited lab value is an unattributable pedigree. Add the citation and "
            f"widen this check deliberately, do not drop the rows.")
    dois = sorted(set(lab.gf_source_doi.astype(str)))
    if dois != [EXPECTED_DOI]:
        raise SystemExit(f"unexpected DOI(s) on the Fe II LAB rows: {dois}")

    miss = [c for c in ("wavelength_air_A", "excitation_potential_eV", "log_gf",
                        "gf_sigma_dex") if lab[c].isna().any()]
    if miss:
        raise SystemExit(
            f"Fe II LAB rows are missing values in {miss}. A lab table without a "
            f"per-line sigma cannot carry a MEASURED systematic, and falling back to "
            f"the generic bound would manufacture precision (RYA-968).")

    out = pd.DataFrame({
        "source": "DenHartog2019",
        "wavelength_air_A": lab.wavelength_air_A.astype(float),
        "elo_eV": lab.excitation_potential_eV.astype(float),
        "loggf": lab.log_gf.astype(float),
        "e_loggf_dex": lab.gf_sigma_dex.astype(float),
    }).sort_values("wavelength_air_A").reset_index(drop=True)
    return out


def main() -> None:
    out = build()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT, index=False)
    vis = out.wavelength_air_A.between(4200, 6910).sum()
    print(f"wrote {OUT.relative_to(ROOT)}  ({len(out)} lines)")
    print(f"  span            : {out.wavelength_air_A.min():.1f} - "
          f"{out.wavelength_air_A.max():.1f} A")
    print(f"  in VIS 4200-6910: {vis}   (all below "
          f"{out[out.wavelength_air_A.between(4200,6910)].wavelength_air_A.max():.1f} A "
          f"-- a BLUE subset, not a VIS-spanning pool)")
    print(f"  cited sigma     : {out.e_loggf_dex.min():.3f} - "
          f"{out.e_loggf_dex.max():.3f} dex, RMS "
          f"{(out.e_loggf_dex**2).mean()**0.5:.4f}")


if __name__ == "__main__":
    main()
