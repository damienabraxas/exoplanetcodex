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
  2. **present in the committed MPIA per-line grid** — so a per-line MPIA delta exists to
     report as the CROSS-ENGINE DIAGNOSTIC (RYA-525). It is no longer the anchor: see below;
  3. **measured by us, in-aggregate** — a gate line we cannot measure proves nothing;
  4. **not weak** — the Ti lesson; weak lines carry the largest EW error;
  5. **ISOLATED from other Fe I** within the gate's own EW integration half-width.

⚠️ THE FIFTH REQUIREMENT IS NEW, AND THE FIRST PASS DID NOT HAVE IT. `ts_gerber_gate.ew()`
integrates `1 - flux` over +/- `EW_HW` = 1.2 A, which is wide. A blend from a DIFFERENT
species is harmless — it adds the same amount to the NLTE spectrum and to every LTE
curve-of-growth point, so it cancels exactly in the EW -> A* inversion. **A second Fe I line
does not cancel**: it responds to the Fe abundance in both spectra, so it dilutes the
line's sensitivity and drags the recovered delta toward the neighbour's. Selecting on
strength alone can therefore pick a pair and report it as one line.

⚠️ THE ANCHOR IS NO LONGER MPIA. The first pass set it to the median of these lines' own
Bergemann-MPIA deltas (+0.0124) and returned CHECK at |median - anchor| = 0.051 vs tol
0.05 — a miss by 0.001 dex. But RYA-785 says to validate against the deck's OWN published
anchor and explicitly forbids MPIA as the target, and RYA-525/712 classify a Gerber-vs-MPIA
difference as a cross-engine DIAGNOSTIC, never a validation criterion. Gerber+2023 (A&A
669 A43) Table 3/4 — the deck's own paper — publishes the solar correction directly:
[Fe/H] -0.04 (1D LTE) -> +0.02 (1D NLTE) = **+0.06 dex**. That is the anchor.

It prints the candidates with everything needed to justify the pick, the MPIA delta per
line as the diagnostic, and the isolation distance actually measured.
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

import argparse

UNSET = {"", "0", "-1", "nan", "none", "None", "NONE"}
EW_HW = 1.2          # MUST track ts_gerber_gate.EW_HW — the integration half-width
ISO_DGF = 1.5        # a neighbour this many dex weaker than the target is not a threat
EW_CSV = ROOT / "data" / "measured" / "band_ew" / \
    "FeI_3800_6910_kpno_solar_atlas_PROFILEFIT_ew.csv"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n", type=int, default=15,
                    help="how many gate lines to pick (default 15; the first pass used 3)")
    ap.add_argument("--ew-min", type=float, default=20.0)
    ap.add_argument("--ew-max", type=float, default=120.0)
    ap.add_argument("--no-isolation", action="store_true",
                    help="reproduce the first pass's four-requirement selection")
    a = ap.parse_args()

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
        # (5) isolation: the nearest OTHER Fe I line inside the EW window. Only a
        # comparably strong one matters -- a neighbour >= ISO_DGF dex weaker contributes
        # negligibly to an integrated EW.
        near = np.where(fe1 & (np.abs(w - c) < EW_HW) & (np.abs(w - c) > 0.02))[0]
        near = [k for k in near if gf[k] > gf[j] - ISO_DGF]
        if near:
            k = near[int(np.argmin(np.abs(w[np.array(near)] - c)))]
            iso_d, iso_gf = float(abs(w[k] - c)), float(gf[k])
        else:
            iso_d, iso_gf = float("inf"), float("nan")
        rows.append(dict(wave=c, ew_mA=float(r.ew_mA), loggf=float(gf[j]),
                         ep_eV=float(ep[j]), mpia_delta=float(d),
                         iso_dA=iso_d, iso_loggf=iso_gf,
                         lvl_low=lo_l[j], lvl_up=up_l[j]))

    t = pd.DataFrame(rows).sort_values("ew_mA", ascending=False)
    print(f"\ncandidates clearing ALL FOUR requirements: {len(t)}")
    if not len(t):
        print("  none — the gate cannot be built from MPIA-overlap lines")
        return
    clean = t[np.isinf(t.iso_dA)]
    print(f"  of those, ISOLATED (no Fe I within {EW_HW} A at > loggf-{ISO_DGF}): "
          f"{len(clean)}")
    print(f"\n{'wave':>10}{'EW mA':>8}{'loggf':>8}{'EP eV':>7}{'MPIA d':>9}"
          f"{'nearest FeI':>12}  low/up levels")
    for _, r in t.head(24).iterrows():
        iso = "clean" if np.isinf(r.iso_dA) else f"{r.iso_dA:.3f}A"
        print(f"{r.wave:>10.3f}{r.ew_mA:>8.1f}{r.loggf:>8.3f}{r.ep_eV:>7.3f}"
              f"{r.mpia_delta:>+9.4f}{iso:>12}  {r.lvl_low}/{r.lvl_up}")

    pool = t if a.no_isolation else clean
    strong = pool[(pool.ew_mA >= a.ew_min) & (pool.ew_mA <= a.ew_max)]
    if len(strong) < a.n:
        print(f"\n  NOTE: only {len(strong)} lines clear every requirement; taking all.")
    # Spread over excitation potential rather than taking the strongest n. Gerber+2023's
    # own Fig. 6 plots the Fe I NLTE correction AGAINST excitation potential, so a gate
    # bunched at one EP cannot see the behaviour the paper reports.
    strong = strong.sort_values("ep_eV")
    idx = np.unique(np.linspace(0, len(strong) - 1, min(a.n, len(strong))).round().astype(int))
    pick = strong.iloc[idx].sort_values("wave")
    print(f"\nPICK ({len(pick)} lines; EW {a.ew_min:.0f}-{a.ew_max:.0f} mA, "
          f"{'isolation enforced' if not a.no_isolation else 'ISOLATION OFF'}, "
          f"spread over EP):")
    for _, r in pick.iterrows():
        print(f"  Fe I {r.wave:9.3f}  EW {r.ew_mA:5.1f} mA  EP {r.ep_eV:5.3f}  "
              f"MPIA delta {r.mpia_delta:+.4f}")
    print(f"\n  EP span {pick.ep_eV.min():.2f}-{pick.ep_eV.max():.2f} eV; "
          f"wavelength span {pick.wave.min():.0f}-{pick.wave.max():.0f} A")
    print(f"  CROSS-ENGINE DIAGNOSTIC (RYA-525), not the anchor: these lines' own MPIA "
          f"deltas median {pick.mpia_delta.median():+.4f} "
          f"({pick.mpia_delta.min():+.4f}..{pick.mpia_delta.max():+.4f})")
    print(f"  ANCHOR is Gerber+2023 A&A 669 A43 Table 3/4 solar 1D NLTE - 1D LTE = +0.06")
    print(f"  waves=[{', '.join(f'{x:.3f}' for x in pick.wave)}]")


if __name__ == "__main__":
    main()
