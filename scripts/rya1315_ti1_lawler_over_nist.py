#!/usr/bin/env python3
"""RYA-1315 -- Ti I: Lawler et al. 2013 gf over the older NIST ASD compilation.

    python scripts/rya1315_ti1_lawler_over_nist.py            # report only
    python scripts/rya1315_ti1_lawler_over_nist.py --apply    # write canonical_gf + audit log

RYA-1268 rule 2 took a VALD lab gf (code LGWSC = Lawler 2013) only when NIST ASD agreed within
0.03 dex; when NIST disagreed the row fell to rule 3 and took NIST's gf. For Ti I that is
backwards: NIST's Ti I values predate Lawler 2013, the lab measurement Scott et al. 2015 (and so
AGSS21) used. Verified against the local Lawler 2013 PDF (Table 4 sample) and Scott 2015 Table 1:
the VALD LGWSC gf reproduces both (4281.363 -1.26, 5022.866 -0.33, 3199.914 +0.31, 3717.391
-1.19), the NIST row does not (-1.359, -0.434, +0.21, -1.222).

Scope: Ti I rows with adjudication_status `adopted_rya1268_nist_asd` whose transition VALD
attributes to LGWSC. The VALD extracts are HFS/isotope split, so the transition gf is the
log of the summed component gf (one upper level, one J pair, one extract file).

sigma: Lawler's per-line errors are in the machine-readable Table 3 only (VizieR J/ApJS/205/11,
unreachable from both hosts on 2026-10-29). Each row takes the paper's stated upper bound,
0.10 dex ("0.02 dex for dominant branches to ~0.10 dex for weak branches"), never NIST's class
(which prices a different gf). Tighten when the table is in hand.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
CANON = ROOT / "data" / "linelists" / "canonical_gf.csv"
VALD = ROOT / "data" / "audit" / "rya1268_line_coverage" / "vald_solar_gf_refs.csv"
AUD = ROOT / "data" / "audit" / "rya1315_ti"
FROM_STATUS = "adopted_rya1268_nist_asd"
NEW_STATUS = "adopted_rya1315_lawler2013_over_nist"
WTOL_A = 0.03            # RYA-1268's own wavelength tolerance
EP_TOL = 0.02
HOLD_MOVE_DEX = 0.5      # RYA-1268 "super bad dex" hold
SIGMA_DEX = 0.10         # Lawler 2013's stated upper bound
DOI = "10.1088/0067-0049/205/2/11"


def lawler_gf(v: pd.DataFrame, w: float, ep: float):
    """(log gf summed over the transition's components, n components) or None."""
    g = v[((v.wavelength_air_A - w).abs() <= 0.15) & ((v.ep_eV - ep).abs() <= EP_TOL)]
    if not len(g):
        return None
    k = g.iloc[(g.wavelength_air_A - w).abs().argmin()]
    if abs(k.wavelength_air_A - w) > WTOL_A or "LGWSC" not in str(k.gf_ref_code):
        return None
    g = g[((g.ep_upper_eV - k.ep_upper_eV).abs() <= 0.002) & (g.j_lower == k.j_lower)
          & (g.j_upper == k.j_upper) & (g.file == k.file)]
    if not g.gf_ref_code.astype(str).str.contains("LGWSC").all():
        return None
    return float(np.log10((10.0 ** g.loggf).sum())), len(g)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--apply", action="store_true")
    a = ap.parse_args(argv)
    cg = pd.read_csv(CANON, low_memory=False, dtype=str, keep_default_na=False)
    v = pd.read_csv(VALD)
    v = v[v.species == "Ti I"].drop_duplicates(
        ["wavelength_air_A", "loggf", "ep_eV", "ep_upper_eV", "file"])
    rows = cg[(cg.species == "Ti I") & (cg.adjudication_status == FROM_STATUS)]
    log, held = [], []
    for i, r in rows.iterrows():
        hit = lawler_gf(v, float(r.wavelength_air_A), float(r.excitation_potential_eV))
        if hit is None:
            continue
        gf, n = hit
        old = float(r.log_gf)
        entry = {"line_id": r.line_id, "wavelength_air_A": r.wavelength_air_A,
                 "ep_eV": r.excitation_potential_eV, "old_loggf": r.log_gf,
                 "old_tier": r.gf_tier, "old_sigma_dex": r.gf_sigma_dex,
                 "old_reference": r.loggf_reference, "old_status": r.adjudication_status,
                 "nist_grade": r.nist_grade, "new_loggf": round(gf, 3),
                 "move_dex": round(gf - old, 3), "vald_components": n}
        if abs(gf - old) > HOLD_MOVE_DEX:
            held.append(entry)
            continue
        log.append(entry)
        if a.apply:
            cg.at[i, "log_gf"] = f"{gf:.3f}"
            cg.at[i, "gf_tier"] = "LAB"
            cg.at[i, "lab_source_tag"] = "lawler2013"
            cg.at[i, "gf_sigma_dex"] = f"{SIGMA_DEX}"
            cg.at[i, "gf_source_doi"] = DOI
            cg.at[i, "loggf_reference"] = (
                f"Lawler et al. 2013 (ApJS 205, 11) via VALD (LGWSC); sigma: Lawler's stated "
                f"upper bound {SIGMA_DEX} dex (per-line errors in the MRT only); supersedes "
                f"NIST ASD {r.log_gf} (class {r.nist_grade or '?'}), RYA-1315")
            cg.at[i, "adjudication_status"] = NEW_STATUS
    d = pd.DataFrame(log)
    print(f"{len(rows)} Ti I rows at {FROM_STATUS}; {len(d)} carry a VALD LGWSC gf; "
          f"{(d.move_dex.abs() > 0.03).sum() if len(d) else 0} move > 0.03 dex; "
          f"{len(held)} held (> {HOLD_MOVE_DEX} dex)")
    if len(d):
        print(d.move_dex.describe().round(3).to_string())
    if a.apply:
        AUD.mkdir(parents=True, exist_ok=True)
        d.to_csv(AUD / "lawler2013_over_nist.csv", index=False)
        pd.DataFrame(held).to_csv(AUD / "lawler2013_over_nist_held.csv", index=False)
        cg.to_csv(CANON, index=False)
        print(f"wrote {CANON.relative_to(ROOT)} and {AUD.relative_to(ROOT)}/lawler2013_over_nist*.csv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
