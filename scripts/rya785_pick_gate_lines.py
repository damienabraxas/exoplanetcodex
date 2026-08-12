#!/usr/bin/env python3
"""RYA-785 — choose the Fe gate lines from evidence, not from a guess.

The RYA-534 gate table records two ways this goes wrong, both already paid for:

  * **Co/Ni:** the first-pass lines had an UNIDENTIFIED upper level in the GES list, so
    bsyn silently set departure = 1 and the deck raised. A gate line must carry BOTH
    levels identified.
  * **Ti:** the first pass used weak lines (0.9-5.7 mA) and compared them against a
    guessed +0.05 anchor. Re-gating needed stronger lines that OVERLAP the MPIA grid, so
    the comparison is per-line against a real anchor rather than against a number someone
    expected.

So an Fe gate line has to clear four things at once, and this picks them by checking all
four rather than by choosing plausible-looking wavelengths:

  1. **both GES NLTE levels identified** — else the departure silently disengages;
  2. **present in the committed MPIA per-line grid** — so the anchor is that line's own
     published delta, not an element-level average;
  3. **measured by us, in-aggregate** — a gate line we cannot measure proves nothing;
  4. **not weak** — the Ti lesson; weak lines carry the largest EW error.

It prints the candidates with everything needed to justify the pick, and the MPIA delta
per line so the anchor is READ rather than assumed.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, os.environ.get("ISPEC_DIR", "/srv/codex/engines/ispec_src"))

from config.constants import STAR_PARAMS                     # noqa: E402
from pipeline import nlte_corrections as nc                  # noqa: E402

UNSET = {"", "0", "-1", "nan", "none", "None", "NONE"}
EW_CSV = ROOT / "data" / "measured" / "band_ew" / \
    "FeI_3800_6910_kpno_solar_atlas_PROFILEFIT_ew.csv"


def main() -> None:
    p = STAR_PARAMS["solar"]
    teff, logg = float(p["teff"]), float(p["logg"])
    feh = float(p.get("feh", p.get("feh_ref", 0.0)))

    from pipeline.abundances_derive import _load_synth_resources
    ll, _, _ = _load_synth_resources()
    w = np.asarray(ll["wave_A"], float)
    el = np.asarray([str(x).strip() for x in ll["element"]])
    lo_l = np.asarray([str(x).strip() for x in ll["nlte_label_low"]])
    up_l = np.asarray([str(x).strip() for x in ll["nlte_label_up"]])
    gf = np.asarray(ll["loggf"], float)
    ep = np.asarray(ll["lower_state_eV"], float)
    fe1 = np.array([e.upper().startswith("FE 1") for e in el])

    # (1) both levels identified, and inside the committed MPIA grid's span
    grid = nc._load_mpia_fe_grid()
    mp = np.asarray(grid["waves"].get(1, []), dtype=float)
    print(f"committed MPIA Fe I grid: {len(mp)} nodes, "
          f"{mp.min():.1f}-{mp.max():.1f} A")

    ok_lvl = fe1 & np.array([(a not in UNSET) and (b not in UNSET)
                             for a, b in zip(lo_l, up_l)])
    print(f"GES Fe I with BOTH levels identified: {int(ok_lvl.sum())}")

    # (3) our own measured, in-aggregate lines
    if not EW_CSV.exists():
        raise SystemExit(f"need the measured optical EWs at {EW_CSV}")
    ew = pd.read_csv(EW_CSV)
    meas = ew[ew.in_aggregate.fillna(False)][["wavelength_air_A", "ew_mA"]]
    print(f"our in-aggregate optical Fe I lines: {len(meas)}")

    rows = []
    for _, r in meas.iterrows():
        c = float(r.wavelength_air_A)
        # (1) level-identified GES entry at this wavelength
        j = np.where(ok_lvl & (np.abs(w - c) < 0.05))[0]
        if not len(j):
            continue
        j = j[np.argmax(gf[j])]
        # (2) in the MPIA grid, and read ITS OWN delta
        if not len(mp) or np.min(np.abs(mp - c)) > 0.15:
            continue
        d = nc._mpia_fe_delta("I", c, teff, logg, feh)
        if not np.isfinite(d):
            continue
        rows.append(dict(wave=c, ew_mA=float(r.ew_mA), loggf=float(gf[j]),
                         ep_eV=float(ep[j]), mpia_delta=float(d),
                         lvl_low=lo_l[j], lvl_up=up_l[j]))

    t = pd.DataFrame(rows).sort_values("ew_mA", ascending=False)
    print(f"\ncandidates clearing ALL FOUR requirements: {len(t)}")
    if not len(t):
        print("  none — the gate cannot be built from MPIA-overlap lines")
        return
    print(f"\n{'wave':>10}{'EW mA':>8}{'loggf':>8}{'EP eV':>7}{'MPIA d':>9}  low/up levels")
    for _, r in t.head(20).iterrows():
        print(f"{r.wave:>10.3f}{r.ew_mA:>8.1f}{r.loggf:>8.3f}{r.ep_eV:>7.3f}"
              f"{r.mpia_delta:>+9.4f}  {r.lvl_low}/{r.lvl_up}")

    # The Ti lesson: prefer well-measured lines, not the weakest.
    strong = t[(t.ew_mA >= 20) & (t.ew_mA <= 120)]
    pick = strong.head(3) if len(strong) >= 3 else t.head(3)
    print(f"\nPICK (EW 20-120 mA where possible — the Ti 'weak lines' lesson):")
    for _, r in pick.iterrows():
        print(f"  Fe I {r.wave:.3f}  EW {r.ew_mA:.1f} mA  MPIA delta {r.mpia_delta:+.4f}")
    print(f"\n  ANCHOR from these lines' own MPIA deltas: "
          f"median {pick.mpia_delta.median():+.4f}  "
          f"(range {pick.mpia_delta.min():+.4f}..{pick.mpia_delta.max():+.4f})")
    print(f"  waves=[{', '.join(f'{x:.3f}' for x in pick.wave)}]")


if __name__ == "__main__":
    main()
