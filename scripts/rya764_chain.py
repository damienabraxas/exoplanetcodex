#!/usr/bin/env python3
"""RYA-764 — walk the 447 -> 14 chain with real counts and find what actually cuts it.

    python3 scripts/rya764_chain.py --ew-csv <measured EWs> [--lines-csv <1D-LTE lines>]

WHY
---
RYA-764 states the chain as

    measured+verified 103  ->  inverts 101  ->  laboratory gf tier 29
                            ->  "in GES below 9199.9 A"  ->  14

and concludes that step 4 is where Engine B loses its lines. `rya764_engine_b_coverage.py`
measured step 4 and found it is not a filter at all: 99 of the 103 verified lines are
Engine-B-capable. So the 29 -> 14 drop has an unidentified cause, and quoting "14" as an
Engine-B number rests on it.

This walks every step against the artifacts and splits the last one by ENGINE, because
the two engines have completely different constraints:

  * Engine A is limited by MPIA's per-line departure table — a SUPPLIER limit (RYA-763).
  * Engine B is limited by GES level identification — which we have now measured at ~96%.

If the 14 is an Engine-A number that was carried across to Engine B by symmetry, that is
the whole answer, and Engine B's real ceiling is much higher than the ticket says.

⚠️ WHICH ENGINE-A SOURCE. There are TWO, and they do not have the same reach:

  * the COMMITTED per-line table `Fe_Bergemann_MPIA.csv`, which this script reads via
    `nlte_corrections._mpia_fe_delta`. Measured span: Fe I 3791.7-6843.7 A, Fe II
    3805.5-6586.7 A -> **zero nodes in 6910-9199.9 A**. RYA-710 says the same thing
    ("stops at 6844 A").
  * the LIVE MPIA query used by `derive_band_products.py`, which served 44 of 101 lines
    in the RYA-713 run and which RYA-763 counts at 28 Fe I lines past 6910 A.

So a "0" from this script is NOT the model's limit — it is the committed extract's limit,
which is exactly the "a derived extract is not the model's limit" trap RYA-710 flags as a
standing lesson at its third occurrence. The number is reported here labelled as the
extract, never as Engine A's capability.
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, os.environ.get("ISPEC_DIR", "/srv/codex/engines/ispec_src"))

from config.constants import STAR_PARAMS                                  # noqa: E402
from pipeline.curate_nonfe_pools import ACCEPTED_GF_TIERS, _gf_tier       # noqa: E402
from pipeline import nlte_corrections as nc                               # noqa: E402
from scripts.rya764_engine_b_coverage import classify                     # noqa: E402
from pipeline.abundances_derive import _load_synth_resources              # noqa: E402

CANON = ROOT / "data" / "linelists" / "canonical_gf.csv"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--element", default="Fe")
    ap.add_argument("--ion", default="I")
    ap.add_argument("--lo", type=float, default=6910.0)
    ap.add_argument("--hi", type=float, default=9199.9)
    ap.add_argument("--ew-csv", required=True)
    ap.add_argument("--lines-csv", help="the 1D-LTE per-line artifact (for the 'inverts' step)")
    a = ap.parse_args()

    p = STAR_PARAMS["solar"]
    teff, logg = float(p["teff"]), float(p["logg"])
    feh = float(p.get("feh", p.get("feh_ref", 0.0)))

    ew = pd.read_csv(a.ew_csv)
    ew = ew[(ew.wavelength_air_A >= a.lo) & (ew.wavelength_air_A <= a.hi)]
    verified = ew[ew.in_aggregate.fillna(False)].copy()

    if a.lines_csv and Path(a.lines_csv).exists():
        ln = pd.read_csv(a.lines_csv)
        inverts = ln[ln.in_aggregate.fillna(False)].copy()
        inv_waves = set(np.round(inverts.wavelength_air_A.astype(float), 3))
    else:
        inv_waves = set(np.round(verified.wavelength_air_A.astype(float), 3))

    # ── gf tier, from the project's own function ──────────────────────────────
    canon = pd.read_csv(CANON, low_memory=False)
    canon["element"] = canon["species"].astype(str).str.split().str[0]
    ci = canon[(canon.element == a.element) &
               (canon.ion == (1 if a.ion.upper() == "I" else 2))]
    cw = ci.wavelength_air_A.astype(float).values

    ll, _, _ = _load_synth_resources()

    rows = []
    for _, r in verified.iterrows():
        c = float(r.wavelength_air_A)
        inverted = round(c, 3) in inv_waves

        j = int(np.abs(cw - c).argmin()) if len(cw) else -1
        if j >= 0 and abs(cw[j] - c) <= 0.05:
            g, ref = ci.iloc[j]["nist_grade"], ci.iloc[j]["loggf_reference"]
            tier = _gf_tier(g, ref)
        else:
            tier, ref = "NOT-IN-CANONICAL", ""

        # COMMITTED extract only — see the docstring. Not Engine A's capability.
        d_a = nc._mpia_fe_delta(a.ion, c, teff, logg, feh)
        served_a = bool(np.isfinite(d_a))

        vb, _, _, _ = classify(ll, c, a.element, a.ion)

        rows.append(dict(wavelength_air_A=c, inverted=inverted, gf_tier=tier,
                         gf_reference=str(ref), engine_a_served=served_a,
                         engine_b_verdict=vb))

    df = pd.DataFrame(rows)
    lab = df[df.gf_tier.isin(ACCEPTED_GF_TIERS)]

    print("\n" + "=" * 78)
    print(f"THE CHAIN, RECOUNTED — {a.element} {a.ion} {a.lo:.0f}-{a.hi:.0f} A")
    print("=" * 78)
    print(f"  measured                                  {len(ew):5d}")
    print(f"  verified (in_aggregate)                   {len(df):5d}")
    print(f"  inverts to an abundance                   {int(df.inverted.sum()):5d}")
    print(f"  laboratory gf tier {sorted(ACCEPTED_GF_TIERS)}        {len(lab):5d}")
    print(f"\n  gf tier breakdown of the verified set:")
    for k, v in df.gf_tier.value_counts().items():
        print(f"    {k:<18} {v:5d}")

    print(f"\n  --- the last step, SPLIT BY ENGINE (the ticket merges them) ---")
    print(f"  of the {len(lab)} laboratory-gf lines:")
    a_served = int(lab.engine_a_served.sum())
    b_capable = int((lab.engine_b_verdict == "ENGINE-B-CAPABLE").sum())
    print(f"    Engine A  committed MPIA extract        {a_served:5d}"
          f"   <- the EXTRACT stops at 6843.7 A")
    print(f"    Engine B  level-identified in GES       {b_capable:5d}"
          f"   <- what Engine B can actually correct")
    both = int(((lab.engine_b_verdict == "ENGINE-B-CAPABLE") & lab.engine_a_served).sum())
    print(f"    BOTH (the cross-engine diagnostic pool) {both:5d}")

    print(f"\n  same split over ALL {len(df)} verified lines (not just laboratory gf):")
    print(f"    Engine A committed extract              {int(df.engine_a_served.sum()):5d}")
    print(f"    Engine B capable                        "
          f"{int((df.engine_b_verdict == 'ENGINE-B-CAPABLE').sum()):5d}")

    # The committed Fe extract's actual span, printed so the 0 is never read as
    # "Engine A cannot do the IR".
    g = nc._load_mpia_fe_grid()
    spans = {i: (float(np.min(w)), float(np.max(w))) for i, w in g['waves'].items() if len(w)}
    print("\n  committed Fe_Bergemann_MPIA.csv span, per ion:")
    for i, (w0, w1) in sorted(spans.items()):
        print(f"    ion {i}: {w0:.1f} - {w1:.1f} A"
              f"   in-band nodes: "
              f"{int(np.sum((np.asarray(g['waves'][i]) >= a.lo) & (np.asarray(g['waves'][i]) <= a.hi)))}")

    print("\n" + "-" * 78)
    print("  VERDICT")
    print(f"  Engine B is NOT limited by GES here: {b_capable} of the {len(lab)} accepted-gf")
    print(f"  lines are level-identified, and {int((df.engine_b_verdict == 'ENGINE-B-CAPABLE').sum())}"
          f" of all {len(df)} verified lines are.")
    print("  So RYA-764's step 4 ('in GES below 9199.9 A') is not the filter it is stated")
    print("  to be, and a '14' attributed to it does not describe Engine B.")
    print(f"\n  The committed MPIA extract contributes {a_served} in-band — but that is the")
    print("  EXTRACT's ceiling (6843.7 A), not the model's; the live MPIA query served 44")
    print("  of 101 in the RYA-713 run. Any Engine-A count must say which source it used.")

    out = ROOT / "data" / "audit" / "rya764"
    out.mkdir(parents=True, exist_ok=True)
    f = out / f"{a.element}{a.ion}_chain_recount.csv"
    df.to_csv(f, index=False)
    print(f"\n  wrote {f.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
