#!/usr/bin/env python3
"""RYA-1232 -- put the O I 926 nm multiplet into canonical_gf, from our own NIST ASD pull.

    python3 scripts/rya1232_canonical_gf_ir_cno.py [--apply]

WHY. The budget joins every accepted line to canonical_gf on lambda AND EP (RYA-1037) and
held every O I 926 nm product: "9260.806: 0 canonical_gf rows". The store was seeded from
the optical synthesis list (GESv6 420-920 nm) and the 926 nm lines were only ever added to
the NIR synthesis list (rya1214_augment_nir_linelist_cno), never to the store -- so
rya1214_adjudicate_cno_gf, which grades EXISTING rows, could not reach them.

SOURCE. data/linelists/primary_gf/nist_asd_OI_3000_25000.tsv (RYA-1160 pull): all nine
components grade A (3%), TP code T5803 (Wiese, Fuhr & Deters 1996, Opacity Project), and
identical to AGSS21's adopted values (Amarsi+2019 Table 1: -0.242, +0.111, +0.224 ...).

ALSO C I 16419.33 (Elgueta+2026 H set): the store carries VALD3 -0.713 with no sigma, which
held both H-band C I products. NIST grades the SAME transition (levels 9.33047 -> 10.08538 eV)
at -0.710, grade D -- but at air 16419.217 A, 0.11 A from VALD's 16419.33, outside the RYA-1214
adjudicator's wavelength tolerance. Graded here by LEVEL identity (lower EP within 0.002 eV,
log gf within 0.1, lambda within 0.3 A, exactly one NIST match), the store's wavelength kept.
Tier NIST-C+ (critically evaluated theory, never LAB -- RYA-1172). sigma from
pipeline.gf_grades.nist_sigma_dex, the store's own grade -> sigma mapping. Air wavelength =
NIST ritz (NIST quotes air above 2000 A).
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline.gf_grades import nist_sigma_dex  # noqa: E402
from pipeline.physical_line_id import physical_id  # noqa: E402

STORE = ROOT / "data/linelists/canonical_gf.csv"
NIST = ROOT / "data/linelists/primary_gf/nist_asd_OI_3000_25000.tsv"
LO, HI = 9260.0, 9267.0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--apply", action="store_true")
    a = ap.parse_args()
    store = pd.read_csv(STORE, low_memory=False)
    nist = pd.read_csv(NIST, sep="\t")
    n = nist[nist.wavelength_ritz_A.between(LO, HI)].copy()
    have = store[(store.species == "O I") & store.wavelength_air_A.between(LO, HI)]
    if len(have) == len(nist[nist.wavelength_ritz_A.between(LO, HI)]):
        print(f"O I 926 nm: {len(have)} rows already present -- nothing to add")
        n = n.iloc[0:0]
    elif len(have):
        raise SystemExit(f"{len(have)} O I 926 nm row(s) partly present -- refusing to guess")
    nxt = max(int(str(x).split("_")[1]) for x in store.line_id if str(x).startswith("gf_")) + 1
    rows = []
    for i, r in enumerate(n.itertuples()):
        w, ep = round(float(r.wavelength_ritz_A), 3), round(float(r.ei_eV), 7)
        rows.append({
            "line_id": f"gf_{nxt + i}", "physical_id": physical_id("O I", w, ep),
            "key_z": 8, "ion": 1.0, "species": "O I", "wavelength_air_A": w,
            "excitation_potential_eV": ep, "hfs_n_components": 1,
            "log_gf": round(float(r.log_gf), 6),
            "loggf_reference": f"NIST ASD grade {r.nist_grade} ({r.ref_transition_probability})",
            "nist_grade": r.nist_grade, "seed_source": "nist_asd_OI (RYA-1160 pull)",
            "adjudication_status": "nist_rya1232_oi926", "in_synth": False, "in_regions": False,
            "in_linelist": True, "gf_sigma_dex": round(nist_sigma_dex(r.nist_grade), 6),
            "is_diagnostic": True, "gf_tier": "NIST-C+", "diagnostic_excluded_by_registry": False})
    add = pd.DataFrame(rows, columns=list(store.columns))
    if len(add):
        print(add[["line_id", "wavelength_air_A", "excitation_potential_eV", "log_gf",
               "nist_grade", "gf_sigma_dex"]].to_string(index=False))
    # IR C/N/O I rows with NO sigma -- grade by LEVEL identity against our NIST pulls.
    # (First case: C I 16419.33, 0.11 A from NIST's air wavelength; then C I 17274.94, 0.06 A.)
    # A row is graded only when EXACTLY ONE NIST transition shares its lower level (EP within
    # 0.002 eV), lies within 0.3 A and agrees in log gf within 0.1 dex; every other case is
    # reported and left as it is.
    pulls = {sp: pd.read_csv(ROOT / f"data/linelists/primary_gf/nist_asd_{sp.replace(' ', '')}_3000_25000.tsv",
                             sep="\t") for sp in ("C I", "N I", "O I")}
    todo = store.index[store.species.isin(list(pulls)) & (store.wavelength_air_A > 9200)
                       & store.gf_sigma_dex.isna()]
    graded, skipped = 0, 0
    for i in todo:
        r0 = store.loc[i]
        nc = pulls[r0.species]
        k = nc[((nc.ei_eV - r0.excitation_potential_eV).abs() < 0.002)
               & ((nc.wavelength_ritz_A - r0.wavelength_air_A).abs() < 0.3)
               & ((nc.log_gf - r0.log_gf).abs() < 0.1)]
        if len(k) != 1 or not isinstance(k.iloc[0].nist_grade, str):
            skipped += 1
            continue
        k = k.iloc[0]
        store.loc[i, ["log_gf", "loggf_reference", "nist_grade", "adjudication_status",
                      "gf_sigma_dex", "gf_tier"]] = [
            round(float(k.log_gf), 6),
            f"NIST ASD grade {k.nist_grade} ({k.ref_transition_probability}); level-identity match "
            f"at air {k.wavelength_ritz_A:.3f} A",
            k.nist_grade, "nist_rya1232_level_identity",
            round(nist_sigma_dex(k.nist_grade), 6), "NIST-C+"]
        graded += 1
        print(f"  {r0.species} {r0.wavelength_air_A}: {r0.loggf_reference} {r0.log_gf} -> NIST "
              f"{k.log_gf:.3f} grade {k.nist_grade} (NIST air {k.wavelength_ritz_A:.3f})")
    print(f"level identity: {graded} graded, {skipped} left (no unique graded NIST match)")
    # Fe I 8446.575 -- the blend on O I 844.6 nm's red component (Ryan's call, option 1).
    # Ruffoni+2014 (lab, listed at 8446.5685) give -1.44 +/- 0.13; VALD3 -1.871. The solar
    # profile 8446.0-8447.0 with A(O) free, on IAG, KP molecfit and Kurucz 2005 independently,
    # is best fitted at -1.7 to -2.0 (rms 0.0195 -> 0.009); at -1.44 the synthesis is 4% too
    # deep at 8446.50-8446.65 even with no oxygen, and the fit pulled A(O) down. VALD's value
    # is adopted for the BLEND; the lab-vs-adopted difference is carried as sigma and priced
    # on 844.6 by a CODEX_GF_OVERRIDE leg. No Fe product measures this line.
    j = store.index[(store.species == "Fe I") & ((store.wavelength_air_A - 8446.575).abs() < 0.002)]
    if len(j) != 1:
        raise SystemExit(f"Fe I 8446.575: {len(j)} store rows")
    store.loc[j[0], ["log_gf", "loggf_reference", "adjudication_status", "gf_sigma_dex"]] = [
        -1.871,
        ("VALD3 -1.871 adopted for a BLEND (RYA-1232): solar-profile fit of O I 844.6 prefers "
         "-1.7..-2.0 on 3 spectra; Ruffoni+2014 lab -1.44 +/- 0.13 (at 8446.5685) is the alternative"),
        "rya1232_blend_solar_profile", round(abs(-1.44 - -1.871), 3)]
    print(f"Fe I 8446.575: -> -1.871 (VALD3), sigma {abs(-1.44 - -1.871):.3f} = lab - adopted")
    if a.apply:
        out = pd.concat([store, add], ignore_index=True)[store.columns]
        out.to_csv(STORE, index=False)
        print(f"wrote {len(add)} rows -> {STORE.relative_to(ROOT)} ({len(out)} rows)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
