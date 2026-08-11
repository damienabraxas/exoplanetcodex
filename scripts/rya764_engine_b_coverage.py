#!/usr/bin/env python3
"""RYA-764 — what Engine B can actually reach, counted from the deck rather than quoted.

    python3 scripts/rya764_engine_b_coverage.py --element Fe --ion I --lo 6910 --hi 9199.9
    python3 scripts/rya764_engine_b_coverage.py --element Fe --ion I --lo 6910 --hi 9199.9 \
        --ew-csv <measured EWs>          # adds the per-line classification

THE QUESTION
------------
RYA-764 asks, as its first scope item, what fraction of our lines GES actually carries
and, for those it does not, *"genuinely absent vs. present-but-unmatched by our lookup."*
That is the right question, and answering it needs one distinction the ticket does not
draw: **being in GES is not the same as being usable by Engine B.**

Engine B is Turbospectrum's TS-native NLTE path. A departure coefficient is indexed by
LEVEL, so a line is only NLTE-capable if BOTH of its levels are identified in the line
list. GES marks every Fe I row `nlte = 'T'` and always populates `nlte_label_low`, but
`nlte_label_up` is frequently the literal string **`none`** — the upper level is not
mapped. Such a line looks NLTE-ready by every surface check and silently synthesises in
LTE, which is exactly the failure mode RYA-764 warns about ("wavelength matching silently
falls back to LTE") arriving through a different door.

So this reports THREE outcomes per line, not two:

  ENGINE-B-CAPABLE   in GES, both levels identified -> a departure coefficient applies
  NO-UPPER-LEVEL     in GES, `nlte_label_up` is 'none' -> WOULD SILENTLY RUN IN LTE
  ABSENT-FROM-GES    no GES entry within the match tolerance at all

Matching is by wavelength ONLY to FIND the row; the verdict is then taken from that row's
level identification. RYA-764 forbids matching *into the deck* by wavelength, and this
does not: it never treats a wavelength match as evidence of NLTE capability.
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

from pipeline.abundances_derive import _load_synth_resources  # noqa: E402

OUT_DIR = ROOT / "data" / "audit" / "rya764"
# GES carries no level identification past this, so Engine B does not exist beyond it.
GES_NLTE_CEILING_A = 9199.9
UNSET = {"", "0", "-1", "nan", "none", "None", "NONE"}


def _s(col) -> np.ndarray:
    return np.asarray([str(x).strip() for x in col])


def deck_inventory(ll, element: str, ion: str, lo: float, hi: float) -> dict:
    """Count what the DECK holds — the numerator RYA-764 quotes, recounted."""
    w = np.asarray(ll["wave_A"], dtype=float)
    el = _s(ll["element"])
    want = f"{element} {1 if ion.upper() == 'I' else 2}".upper()
    is_sp = np.array([e.upper().startswith(want) for e in el])

    m = is_sp & (w >= lo) & (w <= hi)
    lo_lab, up_lab = _s(ll["nlte_label_low"])[m], _s(ll["nlte_label_up"])[m]
    both = np.array([(a not in UNSET) and (b not in UNSET)
                     for a, b in zip(lo_lab, up_lab)])
    return dict(entries=int(m.sum()), level_identified=int(both.sum()),
                no_upper_level=int((~both).sum()), mask=m, both=both,
                waves=w[m])


def classify(ll, wave_A: float, element: str, ion: str, tol: float = 0.05) -> tuple:
    """(verdict, matched_wave, low_label, up_label) for one measured line."""
    w = np.asarray(ll["wave_A"], dtype=float)
    el = _s(ll["element"])
    want = f"{element} {1 if ion.upper() == 'I' else 2}".upper()
    idx = [i for i in np.where(np.abs(w - wave_A) <= tol)[0]
           if el[i].upper().startswith(want)]
    if not idx:
        return ("ABSENT-FROM-GES", np.nan, "", "")
    i = min(idx, key=lambda k: abs(w[k] - wave_A))
    a = str(ll["nlte_label_low"][i]).strip()
    b = str(ll["nlte_label_up"][i]).strip()
    if a not in UNSET and b not in UNSET:
        return ("ENGINE-B-CAPABLE", float(w[i]), a, b)
    return ("NO-UPPER-LEVEL" if a not in UNSET else "NO-LEVEL-ID", float(w[i]), a, b)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--element", default="Fe")
    ap.add_argument("--ion", default="I")
    ap.add_argument("--lo", type=float, default=6910.0)
    ap.add_argument("--hi", type=float, default=GES_NLTE_CEILING_A)
    ap.add_argument("--ew-csv", help="measured EWs, to classify our own lines")
    a = ap.parse_args()

    if a.hi > GES_NLTE_CEILING_A:
        print(f"  NOTE: clipping {a.hi} -> {GES_NLTE_CEILING_A} A; GES carries no level "
              f"identification past it, so Engine B does not exist there (RYA-762).")
        a.hi = GES_NLTE_CEILING_A

    ll, _, _ = _load_synth_resources()

    print("\n" + "=" * 78)
    print(f"ENGINE-B DECK INVENTORY — {a.element} {a.ion}")
    print("=" * 78)
    for lo, hi, label in ((4200.0, 9200.0, "GES full span"),
                          (a.lo, a.hi, f"the band {a.lo:.0f}-{a.hi:.0f}")):
        inv = deck_inventory(ll, a.element, a.ion, lo, hi)
        pct = 100.0 * inv["level_identified"] / max(inv["entries"], 1)
        print(f"\n  {label}")
        print(f"    entries in GES              {inv['entries']:6d}")
        print(f"    BOTH levels identified      {inv['level_identified']:6d}   ({pct:.1f}%)"
              f"   <- what Engine B can actually correct")
        print(f"    upper level 'none'          {inv['no_upper_level']:6d}   "
              f"<- would SILENTLY run in LTE")

    if not a.ew_csv:
        print("\n  (pass --ew-csv to classify our own measured lines)")
        return

    ew = pd.read_csv(a.ew_csv)
    print("\n" + "=" * 78)
    print(f"OUR LINES — {Path(a.ew_csv).name}")
    print("=" * 78)
    rows = []
    for _, r in ew.iterrows():
        c = float(r.wavelength_air_A)
        if not (a.lo <= c <= a.hi):
            continue
        v, mw, lo_l, up_l = classify(ll, c, a.element, a.ion)
        rows.append(dict(wavelength_air_A=c,
                         in_aggregate=bool(r.get("in_aggregate", False)),
                         verdict=v, ges_wave_A=mw,
                         nlte_label_low=lo_l, nlte_label_up=up_l,
                         excluded_reason=str(r.get("excluded_reason", ""))[:80]))
    df = pd.DataFrame(rows)

    for scope, sub in (("ALL measured", df),
                       ("in-aggregate (verified)", df[df.in_aggregate])):
        print(f"\n  {scope}  (n = {len(sub)})")
        for v, k in sub.verdict.value_counts().items():
            print(f"    {v:<20} {k:5d}   ({100.0*k/max(len(sub),1):.1f}%)")

    cap = df[df.in_aggregate & (df.verdict == "ENGINE-B-CAPABLE")]
    silent = df[df.in_aggregate & (df.verdict == "NO-UPPER-LEVEL")]
    print(f"\n  => Engine-B capable among verified lines: {len(cap)}")
    print(f"  => verified lines that would SILENTLY run in LTE: {len(silent)}")
    if len(silent):
        print(f"     {', '.join(f'{x:.3f}' for x in silent.wavelength_air_A.head(10))}"
              f"{' ...' if len(silent) > 10 else ''}")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    f = OUT_DIR / f"{a.element}{a.ion}_{int(a.lo)}_{int(a.hi)}_engine_b_coverage.csv"
    df.to_csv(f, index=False)
    print(f"\n  wrote {f.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
