"""
RYA-249: Flag physically impossible Procyon Fe I lines as blend_flag=True.

Two Fe I lines show EWs that are physically impossible for Procyon's stellar
parameters (Teff=6530K, logg=3.96, [Fe/H]=-0.04) based on SPECTRUM synthesis:
  - Fe I 5247.050 Å  (EW=127.69 mÅ) — predicted <30 mÅ at F5 subgiant
  - Fe I 5762.413 Å  (EW=123.64 mÅ) — predicted <20 mÅ; confirmed blend

Both are set to blend_flag=True so _ew_to_abundance() excludes them (RYA-209).

data/processed/procyon_ew.csv is gitignored, so this script provides the
reproducible record of the flagging decision.
"""

import pandas as pd
from pathlib import Path

EW_CSV = Path(__file__).parent.parent / "data" / "processed" / "procyon_ew.csv"

BLEND_TARGETS = [
    {"wavelength_air_A": 5247.0498, "ion": "I",  "reason": "EW=127.69 mÅ — impossible at Teff=6530K; SPECTRUM predicts <30 mÅ"},
    {"wavelength_air_A": 5762.41276, "ion": "I", "reason": "EW=123.64 mÅ — impossible at Teff=6530K; confirmed atmospheric blend"},
]

TOL = 0.002  # Å matching tolerance


def main() -> None:
    df = pd.read_csv(EW_CSV)
    n_before = int(df["blend_flag"].sum())

    for target in BLEND_TARGETS:
        mask = (
            (df["wavelength_air_A"] - target["wavelength_air_A"]).abs() < TOL
        ) & (df["ion"] == target["ion"])
        hits = mask.sum()
        if hits == 0:
            print(f"WARNING: {target['wavelength_air_A']} Å not found in EW CSV — skip")
            continue
        if hits > 1:
            print(f"WARNING: {target['wavelength_air_A']} Å matched {hits} rows — flagging all")
        df.loc[mask, "blend_flag"] = True
        print(f"  Flagged Fe {target['ion']} {target['wavelength_air_A']} Å  [{target['reason']}]")

    n_after = int(df["blend_flag"].sum())
    df.to_csv(EW_CSV, index=False)
    print(f"\nDone. blend_flag=True: {n_before} → {n_after}  (saved {EW_CSV.name})")


if __name__ == "__main__":
    main()
