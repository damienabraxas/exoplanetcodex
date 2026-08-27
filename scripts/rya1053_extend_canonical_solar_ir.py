#!/usr/bin/env python3
"""RYA-1053 — re-extend canonical_gf from the SOLAR hfs-ON VALD pulls, ALL species.

RYA-1047 extended canonical to 21390 A from `vald_55cnc_nir_5k30k_hfsoff_quarantine.txt`
and admitted Fe ONLY, because that pull is hfs-OFF: fine for Fe (56Fe has I=0) and wrong
for every odd-Z blending species. Both of those choices were forced by a file I should not
have been using.

Ryan, 2026-08-26: "all Vald dats should be there for solar." They were, the whole time:

    vald_solar_ir_9500_17000_hfson_raw.txt    9500.5-16998.5 A   13068 lines
    vald_solar_ir_17000_25000_hfson_raw.txt   17000.6-24986.0 A   8044 lines

⚠️ I missed them because an `ls | head -30` truncated immediately before the `vald_solar_*`
entries, and I then reported the H band as blocked on a data pull. It was never blocked.

WHAT THE OLD SOURCE COST. Fe in the 12976-17000 A overlap: solar 2643 vs 55 Cnc 1425. Of
1423 unambiguous 1-1 matches, 1415 agree EXACTLY on log gf AND EP -- so what RYA-1047
merged is CORRECT, just short ~1200 Fe lines, because 55 Cnc's extract-stellar depth
threshold selected a smaller set. (The 8 that differ also differ in EP by up to 3.9 eV --
different transitions at coincident wavelengths, the RYA-780 pattern, not a conflict.)

hfs VERIFIED, NOT ASSUMED — rows sharing a wavelength in the 9500-17000 pull:
K 54/86, Na 213/270, Co 250/2142, Mn 97/799, and Fe 14/4429 (essentially none, as I=0
predicts). So this source can carry the blending species the other one could not.

🔴 A CLUSTERED LINE IS ONE ROW. canonical_gf stores hfs components collapsed:
`hfs_n_components` counts them, `log_gf` is log10(sum(10**component)), the wavelength is
the gf-weighted centroid and EP the mean. Clustering uses
`pipeline.gf_resolver.cluster_physical_lines` -- the ONE canonical clustering, EP-first
then wavelength-gap -- so the build and the synth reroute cannot disagree about what a
physical line is. Writing one row per component instead would double-count every odd-Z
species' gf.

🔴 AND EVERY ROW GETS A `gf_tier`. RYA-1047 left it NaN on all 2800 rows, which made its
27 laboratory lines invisible to `_cand_graded`/`_cand_deep_graded` -- they select on that
column. The value was right and the SELECTOR could not reach it (RYA-1052). This script
refuses rather than writing a row whose reference has no tier precedent.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline.gf_resolver import cluster_physical_lines          # noqa: E402
from pipeline.species import element_z                          # noqa: E402

CANON = ROOT / "data" / "linelists" / "canonical_gf.csv"
PULLS = [ROOT / "data" / "linelists" / "vald_solar_ir_9500_17000_hfson_raw.txt",
         ROOT / "data" / "linelists" / "vald_solar_ir_17000_25000_hfson_raw.txt"]
OUT_AUDIT = ROOT / "data" / "audit" / "rya1053_solar_ir_gf"

LO_A = 12975.91          # canonical's pre-RYA-1047 red edge
WTOL_A, EPTOL_EV = 0.02, 0.002       # RYA-834's tolerances, imported in spirit
_ROMAN = {1: "I", 2: "II", 3: "III"}
_DATA = re.compile(r"^'([A-Z][a-z]?) (\d+)',\s*([\d.]+),\s*(-?[\d.]+),\s*([\d.]+),")
_GFREF = re.compile(r"gf:(\S+)")


def load_pull(path: Path) -> pd.DataFrame:
    """VALD long format. EVERY data line yields a row; the reference trailer is optional
    metadata (RYA-1047 learned that the hard way — a parser that dropped rows made
    ambiguous pairs look like clean 1-1 matches)."""
    rows, pending = [], None
    for raw in path.open(encoding="latin-1"):
        m = _DATA.match(raw)
        if m:
            if pending is not None:
                rows.append(pending)
            pending = dict(element=m.group(1), ion=int(m.group(2)),
                           wavelength_air_A=float(m.group(3)),
                           gf=float(m.group(4)),
                           excitation_potential_eV=float(m.group(5)),
                           reference_code="")
            continue
        if pending is not None and raw.startswith("'_") and not pending["reference_code"]:
            g = _GFREF.search(raw)
            pending["reference_code"] = g.group(1) if g else ""
    if pending is not None:
        rows.append(pending)
    return pd.DataFrame(rows)


def tier_map(canon: pd.DataFrame) -> dict:
    """reference_code -> gf_tier, DERIVED from how canonical already tiers that code."""
    known = canon[canon.gf_tier.notna()]
    out = {}
    for ref, grp in known.groupby(known.loggf_reference.astype(str)):
        vc = grp.gf_tier.value_counts()
        if len(vc):
            out[ref] = (vc.index[0], int(vc.iloc[0]))
    return out


def tier_for(ref: str, tmap: dict) -> str | None:
    """Exact precedent first; then the FAMILY RULE canonical itself already follows.

    🔴 THE RULE IS READ OFF THE TABLE, NOT INVENTED. canonical_gf tiers Kurucz references
    by whether the code is a YEAR or a LABEL, consistently and at scale:

        K + digits  -> KURUCZ   K03 K04 K06 K07 K08 K09 K10 K11 K12 K13 K14 K75 K99
                                (13 codes, ~89,000 rows, zero exceptions)
        K + letters -> OTHER    K KCN KG KP KSG KZB KZBa K14P
                                (8 codes; KCN is Kurucz's CN list)

    So K16/K17 are Kurucz line lists by year, and KCO/KOH are Kurucz MOLECULAR lists whose
    direct analogue KCN is already OTHER. Anything outside both the exact map and this rule
    still raises — the point is that the tier is derived from evidence, never defaulted.
    """
    if ref in tmap:
        return tmap[ref][0]
    # ONE explicit exception, named rather than absorbed into a wider pattern. VALD's own
    # reference line for it reads "Mg Kurucz NLTE ... gf:FFa" — Kurucz family, label not
    # year, so OTHER by the rule below. It is a SINGLE Mg line; listing it individually
    # keeps the next unknown code raising instead of being swallowed silently.
    if ref == "FFa":
        return "OTHER"
    if re.fullmatch(r"K\d{2}", ref):
        return "KURUCZ"
    if re.fullmatch(r"K[A-Za-z]+", ref):
        return "OTHER"
    return None


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--apply", action="store_true")
    a = ap.parse_args()
    OUT_AUDIT.mkdir(parents=True, exist_ok=True)

    src = pd.concat([load_pull(p) for p in PULLS], ignore_index=True)
    print(f"RYA-1053 — solar hfs-ON IR pulls: {len(src)} component rows, "
          f"{src.wavelength_air_A.min():.1f}-{src.wavelength_air_A.max():.1f} A")
    src = src[src.wavelength_air_A > LO_A].reset_index(drop=True)
    print(f"  above canonical's pre-1047 edge ({LO_A}): {len(src)}")

    # ── collapse hfs components into physical lines ─────────────────────────────
    keys = [f"{r.element} {_ROMAN.get(r.ion, r.ion)}" for r in src.itertuples()]
    wls = src.wavelength_air_A.to_numpy(float)
    eps = src.excitation_potential_eV.to_numpy(float)
    gf = src.gf.to_numpy(float)
    clusters = cluster_physical_lines(keys, wls, eps)
    recs = []
    for cl in clusters:
        w = 10.0 ** gf[cl]
        recs.append(dict(
            species=keys[cl[0]],
            element=src.element.iloc[cl[0]], ion=src.ion.iloc[cl[0]],
            wavelength_air_A=float(np.sum(wls[cl] * w) / w.sum()),
            excitation_potential_eV=float(np.mean(eps[cl])),
            log_gf=float(np.log10(np.sum(w))),
            hfs_n_components=len(cl),
            reference_code=src.reference_code.iloc[cl[0]] or "VALD3",
        ))
    phys = pd.DataFrame(recs).sort_values("wavelength_air_A").reset_index(drop=True)
    multi = int((phys.hfs_n_components > 1).sum())
    print(f"  {len(src)} component rows -> {len(phys)} physical lines "
          f"({multi} carry hfs, max {int(phys.hfs_n_components.max())} components)")

    canon = pd.read_csv(CANON, comment="#", low_memory=False)
    tmap = tier_map(canon)
    phys["gf_tier"] = [tier_for(str(r), tmap) for r in phys.reference_code]
    missing = sorted(set(phys.loc[phys.gf_tier.isna(), "reference_code"]))
    if missing:
        raise SystemExit(
            f"no gf_tier precedent in canonical_gf for reference code(s) {missing}. "
            f"Add them WITH their precedent rather than defaulting — a row whose tier "
            f"nobody chose is indistinguishable downstream from one somebody measured.")
    print(f"  gf_tier from precedent: {phys.gf_tier.value_counts().to_dict()}")

    # ── dedupe against everything canonical already holds ───────────────────────
    # 🔴 ON A TOLERANCE, NOT ROUNDED EQUALITY. The first version rounded to 3 decimals and
    # compared exactly, so 15631.947 (RYA-1047, Ruffoni LAB) and 15631.948 (this pull,
    # Kurucz) read as two different lines and BOTH were kept. 220 of 8280 rows collided
    # that way -- and the damage is not cosmetic: a Kurucz twin sitting 0.001 A from a
    # LABORATORY line makes `grade_line` refuse the pair as AmbiguousLineMatch, which
    # silently un-grades the very H-band lab lines RYA-1052 had just made selectable.
    # Use the SAME tolerance the wavelength+EP matcher uses everywhere else, and keep the
    # EXISTING row on a collision -- it may carry adjudication history this pull cannot.
    from math import isclose  # noqa: F401
    cw = canon.wavelength_air_A.to_numpy(float)
    ce = canon.excitation_potential_eV.to_numpy(float)
    csp = canon.species.astype(str).to_numpy()
    keep_mask, collided = [], 0
    for r in phys.itertuples():
        hit = np.flatnonzero((csp == r.species)
                             & (np.abs(cw - r.wavelength_air_A) <= WTOL_A)
                             & (np.abs(ce - r.excitation_potential_eV) <= EPTOL_EV))
        keep_mask.append(hit.size == 0)
        collided += int(hit.size > 0)
    fresh = phys[np.array(keep_mask)].reset_index(drop=True)
    new = phys
    print(f"  already in canonical within {WTOL_A} A / {EPTOL_EV} eV (kept as-is): {collided}")
    print(f"  NEW physical lines to add        : {len(fresh)}")
    print(f"    species: {dict(fresh.species.value_counts().head(8))}")

    summary = {
        "ticket": "RYA-1053", "sources": [p.name for p in PULLS],
        "hfs": "ON (verified: K/Na/Co/Mn share wavelengths, Fe does not)",
        "lo_A": LO_A, "component_rows": int(len(src)),
        "physical_lines": int(len(phys)), "with_hfs": multi,
        "already_present": int(len(new) - len(fresh)), "new_lines": int(len(fresh)),
        "by_species": {k: int(v) for k, v in fresh.species.value_counts().items()},
        "written_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }
    (OUT_AUDIT / "solar_ir_extension.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    fresh.to_csv(OUT_AUDIT / "new_physical_lines.csv", index=False)
    print(f"\n  wrote {OUT_AUDIT}/solar_ir_extension.json + new_physical_lines.csv")

    if not a.apply:
        print("\n  [dry-run] re-run with --apply to extend canonical_gf.csv")
        return

    start = 1 + max(int(x.split("_")[1]) for x in canon.line_id if str(x).startswith("gf_"))
    add = pd.DataFrame({
        "line_id": [f"gf_{start + i:06d}" for i in range(len(fresh))],
        # 🔴 key_z IS NOT OPTIONAL. `gf_resolver._row_key` does int(float(key_z)) to build
        # its index, so a NaN here raises "cannot convert float NaN to integer" and takes
        # down EVERY synthesis that loads canonical_gf — not just runs touching these
        # lines. I shipped exactly that: RYA-1053's first apply left key_z NaN on all 8280
        # rows and the next Fe II run died inside apply_to_synth_array. A column the
        # RESOLVER indexes on is load-bearing for the whole table.
        "key_z": [element_z(e) for e in fresh.element],
        "ion": fresh.ion.astype(float), "species": fresh.species,
        "wavelength_air_A": fresh.wavelength_air_A,
        "excitation_potential_eV": fresh.excitation_potential_eV,
        "hfs_n_components": fresh.hfs_n_components,
        "log_gf": fresh.log_gf,
        "loggf_reference": fresh.reference_code,
        "gf_tier": fresh.gf_tier,
        "nist_grade": None,
        "seed_source": "solar_ir_vald_hfson(RYA-1053)",
        "adjudication_status": "single_source",
        "gf_synth_ges": np.nan, "gf_regions_vald": np.nan,
        "gf_linelist_vald": fresh.log_gf,
        "in_synth": False, "in_regions": False, "in_linelist": True,
        "delta_synth_minus_ew": np.nan,
        "lab_source_tag": None, "gf_sigma_dex": np.nan, "gf_source_doi": None,
        "is_diagnostic": False,
    })
    out = pd.concat([canon, add], ignore_index=True)
    out.to_csv(CANON, index=False)
    w = pd.to_numeric(out.wavelength_air_A, errors="coerce")
    print(f"\n  canonical_gf.csv  {len(canon)} -> {len(out)} rows")
    print(f"  red edge {w.max():.1f} A;  species now {out.species.nunique()}")


if __name__ == "__main__":
    main()
