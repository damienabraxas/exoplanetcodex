#!/usr/bin/env python3
"""
scripts/rya1009_al_gf_pool_rebuild.py — RYA-1009

Re-source the Burheim-2023-covered Al I lines in `canonical_gf.csv` from EXPERIMENTAL
oscillator strengths, as a POOL REBUILD.

🔴 THE FIREWALL IS THE WHOLE TICKET (RYA-161). The swap is made because Burheim, Hartman
& Nilsson 2023 is an EXPERIMENTAL measurement (FTS branching fractions x radiative
lifetimes, 2-11% per-line accuracy) and the incumbent is Kurucz-1995 semi-empirical.
Better provenance, full stop. The justification is fixed BEFORE the abundance effect is
looked at, and it would be identical if the swap moved A(Al) the other way. "It agrees
better with 6.43" is never the reason and never appears in this file.

WHAT IS SWAPPED, AND WHAT IS DELIBERATELY NOT
---------------------------------------------
Of the 8 census lines Burheim covers, only 4 exist in `canonical_gf` at all (its Al span
ends at 11254.9), and of those:

  6696.023  SWAP   -1.347 -> -1.460   Kurucz-1995 -> experiment
  6698.673  SWAP   -1.647 -> -1.760   Kurucz-1995 -> experiment
  11253.189 RE-SOURCE ONLY, value unchanged: canonical ALREADY carries +0.167, identical
            to Burheim to 3 dp. The gf does not move; what it gains is a cited source, a
            per-line sigma and a DOI in place of the bare tag "VALD3". Worth doing
            precisely because a correct number with unstated provenance is the thing this
            project keeps getting caught by.
  11254.924 ⚠️ HELD BACK — NOT SWAPPED, and this is a physics decision, not caution.
            Burheim's +0.327 is the STRONG COMPONENT ALONE: his Table 1 omits the weak
            3d 2D5/2 - 4f 2F5/2 partner as "more than an order of magnitude weaker".
            canonical's +0.3538 carries hfs_n_components=2, i.e. the BLENDED feature the
            spectrograph actually sees. Substituting a single-component value for a
            two-component sum would UNDER-COUNT the feature's opacity by 0.027 dex — a
            real error introduced while "improving" provenance. It needs the weak partner
            restored before it can be swapped, which is a line-list job, not this one.

The remaining 4 (12749.9, 12757.3, 13123.4, 13150.8) are absent from canonical_gf
entirely and are recorded NOT-SWAPPED with that reason.
"""
from __future__ import annotations
import json, sys
from pathlib import Path
import numpy as np, pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path: sys.path.insert(0, str(ROOT))

CANON = ROOT/"data/linelists/canonical_gf.csv"
LAB   = ROOT/"data/reference/al_gf_lab/al1_lab_loggf.csv"
OUT   = ROOT/"data/results/rya1009"

CITE = "Burheim, Hartman & Nilsson 2023, A&A 672, A197 (experimental Al I log gf)"
DOI  = "10.1051/0004-6361/202245394"
WTOL, EPTOL = 0.06, 0.02

#: Lines Burheim covers that must NOT be swapped, with the reason carried in the artifact.
HOLD = {
    11254.924: ("Burheim's +0.327 is the STRONG COMPONENT ALONE (his Table 1 omits the "
                "weak 3d 2D5/2 - 4f 2F5/2 partner as >1 order of magnitude weaker); "
                "canonical carries hfs_n_components=2, the BLENDED feature. Swapping a "
                "single-component value onto a two-component sum would under-count the "
                "opacity by 0.027 dex. Restore the weak partner first."),
}


def main() -> int:
    lab = pd.read_csv(LAB)
    can = pd.read_csv(CANON, low_memory=False)
    rows, n_swap = [], 0

    for _, b in lab[lab.wavelength_air_A < 25000].iterrows():
        lam, ep = float(b.wavelength_air_A), float(b.elo_eV)
        m = can[(can.species == "Al I")
                & ((can.wavelength_air_A - lam).abs() <= WTOL)
                & ((can.excitation_potential_eV - ep).abs() <= EPTOL)]
        rec = dict(burheim_lam_air_A=round(lam, 4), burheim_lam_vac_A=float(b.lam_vac_A),
                   elo_eV=round(ep, 4), burheim_log_gf=float(b.loggf),
                   burheim_sigma_dex=round(float(b.e_loggf_dex), 4),
                   burheim_unc_pct=float(b.unc_pct),
                   transition=f"{b.lower_level} - {b.upper_level}")
        if m.empty:
            rec.update(action="NOT-SWAPPED", old_log_gf=None, new_log_gf=None,
                       delta_log_gf=None,
                       reason="absent from canonical_gf (its Al span ends at 11254.9 A)")
            rows.append(rec); continue

        i = m.index[0]
        old = float(can.at[i, "log_gf"])
        held = next((v for k, v in HOLD.items() if abs(k - float(can.at[i, "wavelength_air_A"])) < 0.01), None)
        if held:
            rec.update(action="HELD", old_log_gf=old, new_log_gf=None, delta_log_gf=None,
                       canonical_lam_A=float(can.at[i, "wavelength_air_A"]), reason=held)
            rows.append(rec); continue

        can.at[i, "log_gf"] = float(b.loggf)
        can.at[i, "loggf_reference"] = CITE
        can.at[i, "gf_sigma_dex"] = float(b.e_loggf_dex)
        can.at[i, "gf_source_doi"] = DOI
        can.at[i, "lab_source_tag"] = "Burheim2023"
        can.at[i, "gf_tier"] = "LAB"
        can.at[i, "adjudication_status"] = "adjudicated_rya1009"
        n_swap += 1
        rec.update(action=("SWAP" if abs(float(b.loggf) - old) > 1e-4 else "RE-SOURCE"),
                   old_log_gf=old, new_log_gf=float(b.loggf),
                   delta_log_gf=round(float(b.loggf) - old, 4),
                   canonical_lam_A=float(can.at[i, "wavelength_air_A"]),
                   reason=("experimental gf supersedes the semi-empirical incumbent "
                           "(provenance, RYA-161)"))
        rows.append(rec)

    can.to_csv(CANON, index=False)
    OUT.mkdir(parents=True, exist_ok=True)
    d = pd.DataFrame(rows)
    d.to_csv(OUT/"rya1009_al_gf_swap_per_line.csv", index=False)
    (OUT/"rya1009_al_gf_swap_summary.json").write_text(json.dumps(dict(
        ticket="RYA-1009", source=CITE, doi=DOI,
        n_covered=int(len(d)), n_rows_written=n_swap,
        n_swapped=int((d.action == "SWAP").sum()),
        n_resourced_value_unchanged=int((d.action == "RE-SOURCE").sum()),
        n_held=int((d.action == "HELD").sum()),
        n_not_in_canonical=int((d.action == "NOT-SWAPPED").sum()),
        firewall=("RYA-161 — swapped because Burheim is EXPERIMENTAL and the incumbent is "
                  "Kurucz-1995 semi-empirical. Provenance only. The justification was "
                  "fixed before the abundance effect was examined and is unchanged by "
                  "its direction."),
    ), indent=2) + "\n")

    pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 34)
    print(d[["burheim_lam_air_A", "elo_eV", "action", "old_log_gf", "new_log_gf",
             "delta_log_gf", "burheim_sigma_dex"]].to_string(index=False))
    print(f"\n{n_swap} canonical_gf row(s) rewritten "
          f"({int((d.action=='SWAP').sum())} value changes, "
          f"{int((d.action=='RE-SOURCE').sum())} provenance-only)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
