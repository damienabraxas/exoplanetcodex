"""
RYA-263: Flag Fe II 5376.466 Å as SPECTRUM-incompatible for Procyon.

Root cause: SPECTRUM's expint() function (Taylor series for E_2(x)) fails to
converge in 100 iterations for Fe II 5376.466 Å when computed with Procyon's
atmosphere (Teff=6530 K, logg=3.96). The failure calls nrerror() → exit(1),
killing the child process before it can enqueue results. All 87 other lines then
return their default 0.0 values — the pipeline "passes" gates only because the
RYA-261 gate conversion adds 7.46 to the zero, producing A(Fe I) = 7.46.

The line is physically real (Fe II, EW=35.53 mÅ, no blend contamination).
The failure is a SPECTRUM numerical limitation at F-star temperatures and logg.
Fix: exclude via blend_flag so _ew_to_abundance() skips it.

Since data/processed/procyon_ew.csv is gitignored, this script provides the
reproducible record of the fix.

Related: flag_procyon_blends_rya249.py (applies two impossible-EW blend flags).
Run both scripts to fully reproduce the procyon_ew.csv state.
"""

import pandas as pd
from pathlib import Path

EW_CSV = Path(__file__).parent.parent / "data" / "processed" / "procyon_ew.csv"

FIX_TARGET = {
    "wavelength_air_A": 5376.46583,
    "ion": "II",
    "reason": (
        "Fe II 5376.466 Å (EW=35.53 mÅ) — SPECTRUM expint series fails to "
        "converge for Procyon Teff=6530K/logg=3.96 atmosphere. The expint(x,2) "
        "call in flux.c exits via nrerror(), killing the abundance subprocess "
        "and causing all lines to return 0.0. Identified by binary search "
        "in scripts/_find_crash_line_rya263.py. Not a physical blend; "
        "SPECTRUM numerical limitation at F-star parameters."
    ),
}

TOL = 0.002  # Å


def main() -> None:
    df = pd.read_csv(EW_CSV)
    n_before = int(df["blend_flag"].sum())

    mask = (
        (df["wavelength_air_A"] - FIX_TARGET["wavelength_air_A"]).abs() < TOL
    ) & (df["ion"] == FIX_TARGET["ion"])

    hits = mask.sum()
    if hits == 0:
        print(f"WARNING: {FIX_TARGET['wavelength_air_A']} Å not found — skip")
        return
    df.loc[mask, "blend_flag"] = True
    print(f"  Flagged Fe {FIX_TARGET['ion']} {FIX_TARGET['wavelength_air_A']} Å")
    print(f"  Reason: {FIX_TARGET['reason'][:80]}...")

    n_after = int(df["blend_flag"].sum())
    df.to_csv(EW_CSV, index=False)
    print(f"\nDone. blend_flag=True: {n_before} → {n_after}  (saved {EW_CSV.name})")


if __name__ == "__main__":
    main()
