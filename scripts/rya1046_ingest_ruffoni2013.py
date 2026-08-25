#!/usr/bin/env python3
"""Ingest Ruffoni et al. 2013 H-band Fe I laboratory log gf — RYA-1046.

WHY THIS TABLE. `canonical_gf`'s laboratory pool stops at 11316.1 A (Belmonte 2017
reaches 1033 nm), so Elgueta's CRIRES+ J (11796-13195 A) and H (15007-17494 A) arms hold
ZERO cited-lab Fe I lines between them and would grade rung 1 — the ungraded 0.17 blanket
— no matter how well measured. Ruffoni, Pickering, Allende Prieto & Nave 2013 (ApJ 779,
17) measured 28 Fe I transitions across 1.4-1.7 um precisely because no experimental
values existed there. 25 of them land inside Elgueta's H arm.

🔴 THE LADENBURG COLUMN, NOT "RECOMMENDED". Table 6 carries three log gf columns and the
choice is a FIREWALL decision, not a preference:

  * "BF & Effective tau" — branching fractions times upper-level lifetimes that were
    REFINED BY FITTING THE SOLAR SPECTRUM (their Table 3, a `delta log(gf)` column), and
    whose uncertainty explicitly includes 12% from the 0.05 dex solar Fe abundance of
    Asplund et al. (2005). Using those to MEASURE a solar Fe abundance is circular.
  * "Recommended" — tracks BF & Effective tau, so it inherits the same dependence.
  * "Ladenburg" — relative line intensities in laboratory emission and absorption
    (their Tables 4-5). No solar input.

That is the same test that keeps `melendez1999` out of the lab tier despite covering J
and H: solar-tuned gf cannot referee a solar abundance (RYA-161). The two methods agree
to a median 0.010 dex (max 0.050), which is reassuring but is NOT the reason to choose —
"it barely matters" would still admit a solar-anchored number into a rung-3 tier.

⚠️ STATED LIMITATION, not resolved here. The Ladenburg network is placed on an absolute
scale via a reference transition chosen BECAUSE its log gf was unchanged by the
refinement. That reads as independence, but the paper does not state that the anchor's
absolute value never touched the solar fit. Recorded rather than claimed away.

🔴 UPPER LIMITS ARE EXCLUDED, and the data excludes them for us. Seven transitions are
non-detections, flagged `<` in the BF column (their footnote c: Kurucz predicts BFs of at
least 1% "yet they are unobserved in our spectra"). NOT ONE has a Ladenburg value, so
taking that column drops them by construction. An upper limit is a bound, not a
measurement, and 7 of them entering the lab tier would manufacture rung-3 lines out of
non-detections.

WAVELENGTHS ARE AIR, PROVEN not assumed. Converting sigma (cm^-1) to vacuum and dividing
by n_air reproduces the tabulated values to 0.001 A, against 4.4 A if they were vacuum —
a ~4000x discriminator. At 1.6 um a 4.4 A error would place every line on the wrong
feature and still fit something (RYA-944's lesson: the ticket said vacuum, the data was
air).

LEVEL ENERGIES: Table 6 gives only the transition wavenumber, so E_lower = E_upper minus
sigma_ul, with E_upper from their Table 2 (itself from Nave et al. 1994).
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "reference" / "fe_gf_lab" / "fe1_lab_loggf_ruffoni2013.csv"
SOURCE_TAG = "Ruffoni2013"
CM1_TO_EV = 1.239841984e-4


def n_air(lam_um: float) -> float:
    s2 = (1.0 / lam_um) ** 2
    return 1 + 1e-8 * (8342.13 + 2406030 / (130 - s2) + 15997 / (38.9 - s2))


def parse_table2(p: Path) -> dict[str, float]:
    """upper-level label -> energy (cm^-1). Their Table 2, from Nave et al. (1994)."""
    out, cfg = {}, None
    for ln in p.read_text(errors="replace").splitlines():
        f = [x.strip() for x in ln.split("\t")]
        if len(f) < 4:
            continue
        term, J, E = f[1], f[2], f[3]
        try:
            energy = float(E)
            j = int(J)
        except ValueError:
            continue
        out[f"{term}_{j}"] = energy
    return out


def parse_table6(p: Path) -> list[dict]:
    rows, cur = [], None
    for ln in p.read_text(errors="replace").splitlines():
        f = [x.strip() for x in ln.rstrip("\n").split("\t")]
        if len(f) < 12:
            continue
        if f[0]:
            cur = f[0]
        try:
            lam_nm = float(f[2].replace("*", "").strip())
            sigma = float(f[3].replace("^b", "").strip())
        except ValueError:
            continue
        is_limit = f[4].strip().startswith("<")

        def num(i):
            try:
                return float(f[i])
            except (ValueError, IndexError):
                return None

        rows.append(dict(upper=cur, lower=f[1], lam_nm=lam_nm, sigma_cm1=sigma,
                         is_upper_limit=is_limit,
                         lad=num(9), lad_e=num(10), eff=num(7), rec=num(11)))
    return rows


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tables", required=True,
                    help="directory holding apj485713t2_ascii.txt and t6")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    d = Path(a.tables)
    t2, t6 = d / "apj485713t2_ascii.txt", d / "apj485713t6_ascii.txt"
    for f in (t2, t6):
        if not f.exists():
            raise SystemExit(f"missing {f}")

    energies = parse_table2(t2)
    rows = parse_table6(t6)
    print(f"Table 2: {len(energies)} upper levels;  Table 6: {len(rows)} transitions "
          f"({sum(r['is_upper_limit'] for r in rows)} upper limits)")

    # the Ladenburg column drops the limits by construction — assert it, never assume it
    bad = [r for r in rows if r["is_upper_limit"] and r["lad"] is not None]
    if bad:
        raise SystemExit(f"{len(bad)} upper-limit rows carry a Ladenburg value; the "
                         f"column can no longer be trusted to exclude non-detections")

    out, no_energy, not_air = [], [], []
    for r in rows:
        if r["lad"] is None:
            continue
        lam_A = r["lam_nm"] * 10.0
        # AIR CHECK, per line: sigma -> vacuum -> /n_air must reproduce the tabulation
        vac_nm = 1e8 / r["sigma_cm1"] / 10.0
        air_nm = vac_nm / n_air(vac_nm / 1000.0)
        if abs(air_nm - r["lam_nm"]) * 10.0 > 0.05:
            not_air.append((r["lam_nm"], abs(air_nm - r["lam_nm"]) * 10.0))
            continue
        eup = energies.get(str(r["upper"]).strip())
        if eup is None:
            no_energy.append(r["upper"]); continue
        elo = eup - r["sigma_cm1"]
        out.append(dict(source=SOURCE_TAG, wavelength_air_A=round(lam_A, 4),
                        elo_cm1=round(elo, 3), eup_cm1=round(eup, 3),
                        elo_eV=round(elo * CM1_TO_EV, 6),
                        eup_eV=round(eup * CM1_TO_EV, 6),
                        loggf=r["lad"], e_loggf_dex=r["lad_e"]))
    if not_air:
        print(f"  🔴 {len(not_air)} rows failed the AIR check and were dropped: {not_air[:3]}")
    if no_energy:
        print(f"  ⚠️  {len(no_energy)} rows had no Table-2 upper level: "
              f"{sorted(set(no_energy))[:4]}")
    df = pd.DataFrame(out).sort_values("wavelength_air_A").reset_index(drop=True)
    print(f"\n  ingestible Ladenburg lines: {len(df)}")
    if len(df):
        print(f"  span {df.wavelength_air_A.min():.1f}-{df.wavelength_air_A.max():.1f} A")
        print(f"  in Elgueta H (15007-17494): "
              f"{int(df.wavelength_air_A.between(15007,17494).sum())}")
        print(f"  per-line sigma: median {df.e_loggf_dex.median():.3f} dex "
              f"(min {df.e_loggf_dex.min():.3f} max {df.e_loggf_dex.max():.3f})")
    if a.dry_run:
        print("\n[dry-run] nothing written"); return 0
    OUT.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT, index=False)
    print(f"\nwrote {OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
