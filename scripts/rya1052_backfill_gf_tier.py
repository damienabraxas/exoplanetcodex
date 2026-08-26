#!/usr/bin/env python3
"""RYA-1052 — backfill `gf_tier` on the RYA-1047 H-band rows. THEY WERE INVISIBLE.

🔴 THE BUG. RYA-1047 appended 2800 Fe rows (12976-21400 A) and adjudicated 27 of them
against Ruffoni 2013's laboratory gf. It set `lab_source_tag`, `gf_sigma_dex`,
`gf_source_doi` and `adjudication_status` — and left `gf_tier` as NaN on all 2800.

`_cand_graded` and `_cand_deep_graded` in `derive_band_products` select on

    cg.gf_tier.astype(str).str.contains("LAB")

so those 27 laboratory lines COULD NOT BE SELECTED FOR ANY PRODUCT. They graded correctly
through `grade_line()` — which reads `lab_lines()`, a different path — so every check I ran
said rung 3 while the product path could not see them at all. I reported "CRIRES+ H-band
Fe I is off rung 1" on the strength of the grader. That was premature: the value was right
and the SELECTOR could not reach it, which is this project's most repeated failure shape
and the reason RYA-1052's census exists.

Found by the tier census, which reported 0 lab lines in H when Ruffoni put 25 there.

THE MAPPING IS DERIVED FROM PRECEDENT, NOT INVENTED. Each reference code is tiered the way
canonical_gf ALREADY tiers that same code elsewhere:

    K14   -> KURUCZ   (315 precedent rows)
    VALD3 -> VALD3    (40597)
    RU    -> OTHER    (8074)   [Raassen & Uylings — NOT Ruffoni, RYA-945]
    BRW, KCN -> OTHER (no precedent; OTHER is the existing unclassified bucket)
    PRIMARY LAB Ruffoni2013 -> LAB
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
CANON = ROOT / "data" / "linelists" / "canonical_gf.csv"

BY_PRECEDENT = {"K14": "KURUCZ", "VALD3": "VALD3", "RU": "OTHER"}
NO_PRECEDENT = {"BRW": "OTHER", "KCN": "OTHER"}


def main() -> None:
    c = pd.read_csv(CANON, comment="#", low_memory=False)
    mask = c.seed_source.astype(str).str.contains("RYA-1047", na=False)
    n = int(mask.sum())
    if not n:
        raise SystemExit("no RYA-1047 rows in canonical_gf — nothing to backfill")
    already = int(c.loc[mask, "gf_tier"].notna().sum())
    print(f"RYA-1047 rows: {n}   already tiered: {already}")

    lab = mask & c.adjudication_status.astype(str).eq("lab_rya1047")
    c.loc[lab, "gf_tier"] = "LAB"
    # 🔴 `lab_source_tag` uses MACHINE-SHORT codes (DH14/RU14/BEL17/DH19), not the long
    # source names the lab CSVs carry. RYA-1047 wrote the long name and nothing caught it
    # until RYA-945's citation test. Import the map rather than restating it — a second
    # copy of this vocabulary is how the two drift.
    sys.path.insert(0, str(ROOT / "scripts"))
    from rya945_ingest_fe_lab_gf import LAB_TAG
    long_names = set(c.loc[lab, "lab_source_tag"].dropna().astype(str))
    for nm in long_names:
        if nm not in LAB_TAG:
            raise SystemExit(
                f"lab_source_tag {nm!r} has no short code in rya945_ingest_fe_lab_gf."
                f"LAB_TAG. Add it there — that dict is the single source (RYA-353).")
        c.loc[lab & c.lab_source_tag.astype(str).eq(nm), "lab_source_tag"] = LAB_TAG[nm]
        print(f"  lab_source_tag {nm} -> {LAB_TAG[nm]}")
    print(f"  LAB (adjudicated): {int(lab.sum())}")

    for ref, tier in {**BY_PRECEDENT, **NO_PRECEDENT}.items():
        m = mask & ~lab & c.loggf_reference.astype(str).eq(ref)
        if int(m.sum()):
            src = "precedent" if ref in BY_PRECEDENT else "no precedent -> OTHER bucket"
            c.loc[m, "gf_tier"] = tier
            print(f"  {ref:6} -> {tier:7} {int(m.sum()):5d}   ({src})")

    left = mask & c.gf_tier.isna()
    if int(left.sum()):
        # 🔴 REFUSE rather than sweeping the remainder into OTHER. An unrecognised
        # reference is a fact about our coverage, not a row to tidy away (RYA-833).
        raise SystemExit(
            f"{int(left.sum())} RYA-1047 row(s) carry a reference with no mapping: "
            f"{sorted(set(c.loc[left, 'loggf_reference'].astype(str)))}. Add it to the "
            f"table above WITH its precedent rather than defaulting it.")

    c.to_csv(CANON, index=False)
    print(f"\n  wrote {CANON.name}")
    lab_rows = c[c.gf_tier.astype(str).str.contains("LAB", na=False)]
    w = pd.to_numeric(lab_rows.wavelength_air_A, errors="coerce")
    print(f"  gf_tier=LAB rows now: {len(lab_rows)}   reach {w.min():.1f}-{w.max():.1f} A")


if __name__ == "__main__":
    main()
