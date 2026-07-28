#!/usr/bin/env python3
"""
scripts/smoke_gdas_vesta.py
===========================
RYA-380 Phase 1 — GDAS feasibility smoke test on the Vesta CRIRES+ nights
(2022-11-21…25, Paranal).

Confirms (without running molecfit) that the standing per-night GDAS recipe is sound:
  (a) a REAL per-night GDAS profile is retrievable for every Vesta night via the
      reusable pipeline.telluric.gdas_fetch (cache → ESO tarball → NOAA ARL),
  (b) the loud-fail guard works — an unavailable profile RAISES GDASUnavailable rather
      than silently falling back to a standard atmosphere (the RYA-373 CRITICAL bug),
  (c) the nearest-3-hourly slot resolution lands on a real GDAS hour.

The molecfit before/after RYA-390 residual-telluric metric is Phase 2
(`python -m pipeline.crires_telluric --condition` → `python -m
pipeline.co_validation_rya390`), because regenerating mtrans needs the molecfit run.
NOTE: at branch tip the driver ALREADY feeds real GDAS, so the Phase-1 "before"
(standard-atm) product no longer exists on disk — the relevant comparison is the
RYA-390 verdict on the real-GDAS product, recorded in Phase 2.

Usage:  python -m scripts.smoke_gdas_vesta
"""
from __future__ import annotations

import datetime as dt

from config.constants import PATHS, get_site
from pipeline.telluric.gdas_fetch import (GDASUnavailable, fetch_gdas,
                                          gdas_cache_path, nearest_3hourly)

VESTA_NIGHTS = [dt.date(2022, 11, d) for d in (21, 22, 23, 24, 25)]
SITE = 'paranal'


def main() -> int:
    rec = get_site(SITE)
    print("=" * 78)
    print("  RYA-380 Phase 1 — GDAS feasibility smoke test (Vesta CRIRES+, Paranal)")
    print("=" * 78)
    print(f"  site={SITE} lat={rec['lat']} lon={rec['lon']} gdas_loc={rec['gdas_loc']}")
    print(f"  cache={PATHS['gdas_cache']}\n")

    ok = []
    for n in VESTA_NIGHTS:
        slot = nearest_3hourly(n)
        try:
            p = fetch_gdas(SITE, night=n)
            cached = gdas_cache_path(PATHS['gdas_cache'], rec['gdas_loc'], slot)
            print(f"  GDAS {n} (slot {slot:%Y-%m-%dT%H}UT): OK -> {p.name} "
                  f"({p.stat().st_size} B){'  [cached]' if p == cached else ''}")
            ok.append(n)
        except GDASUnavailable as e:
            print(f"  GDAS {n}: UNAVAILABLE -> {e}")

    print(f"\n  {len(ok)}/{len(VESTA_NIGHTS)} nights have a REAL per-night GDAS profile")

    # Loud-fail guard: an unknown site / unstaged manual pull must RAISE, never silently
    # return a standard atmosphere. Prove the guard fires.
    print("\n  loud-fail guard (no silent standard-atm fallback):")
    try:
        fetch_gdas('nowhere-no-tarball', night=VESTA_NIGHTS[0],
                   lat=0.0, lon=0.0, gdas_loc='C-999.9-99.9')
        print("    FAIL — expected GDASUnavailable, got a profile")
        guard_ok = False
    except GDASUnavailable as e:
        print(f"    OK — GDASUnavailable raised: {str(e)[:72]}…")
        guard_ok = True

    passed = (len(ok) == len(VESTA_NIGHTS)) and guard_ok
    print("\n" + "-" * 78)
    print(f"  Phase 1 {'PASS' if passed else 'INCOMPLETE'}: "
          f"{'all Vesta nights retrievable + loud-fail confirmed' if passed else 'see above'}")
    print("  → Phase 2: python -m pipeline.crires_telluric --condition  (re-condition +")
    print("     persist mtrans) then python -m pipeline.co_validation_rya390  (validate)")
    print("=" * 78)
    return 0 if passed else 1


if __name__ == '__main__':
    raise SystemExit(main())
