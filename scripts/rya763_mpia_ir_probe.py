#!/usr/bin/env python3
"""RYA-763 — ask MPIA directly how far Engine A reaches for Fe. Bounded, read-only.

The committed extract stops at 6843.7 A and the live query has only ever been quoted
second-hand ("28 Fe I past 6910 A of 530"). This asks the service itself, with real
in-band Fe I wavelengths from our own accounting catalogue, so the answer is measured
rather than inherited.

Deliberately small: a few batches of 40 (the size the existing scraper uses because
"MPIA dislikes very long batches"). This is a reach probe, not a harvest.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.build_nlte_grids_mpia import _submit_batch  # noqa: E402
from config.constants import STAR_PARAMS                 # noqa: E402

ACCOUNTING = ROOT / "data" / "audit" / "line_accounting" / "per_line.csv"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ion", default="I", choices=["I", "II"])
    ap.add_argument("--batches", type=int, default=3)
    ap.add_argument("--per-batch", type=int, default=40)
    ap.add_argument("--bands", default="4200:6800,6910:9199,9199:13000")
    a = ap.parse_args()

    p = STAR_PARAMS["solar"]
    node = [dict(name="sun", teff_K=int(round(float(p["teff"]))),
                 logg=float(p["logg"]), feh=0.0)]
    code = {"I": "26.01", "II": "26.02"}[a.ion]

    acc = pd.read_csv(ACCOUNTING)
    # per_line.csv stores ion as a ROMAN STRING (I/II/III), not an integer.
    fe = acc[(acc.element.astype(str) == "Fe") & (acc.ion.astype(str) == a.ion)]
    print(f"Fe {a.ion}: {len(fe)} lines in the accounting catalogue")

    out = []
    for band in a.bands.split(","):
        lo, hi = (float(x) for x in band.split(":"))
        sel = fe[(fe.wave_air_A >= lo) & (fe.wave_air_A <= hi)]
        waves = np.sort(sel.wave_air_A.astype(float).values)
        if not len(waves):
            print(f"\n  {lo:.0f}-{hi:.0f} A: no catalogued Fe {a.ion} lines")
            continue
        # even sample across the band so one crowded region cannot dominate
        n = min(a.batches * a.per_batch, len(waves))
        idx = np.unique(np.linspace(0, len(waves) - 1, n).astype(int))
        sample = waves[idx]

        served = 0
        for i in range(0, len(sample), a.per_batch):
            chunk = [float(x) for x in sample[i:i + a.per_batch]]
            try:
                r = _submit_batch(chunk, code, node).get("sun", {})
            except Exception as e:
                print(f"      batch {i//a.per_batch} failed: {type(e).__name__}: "
                      f"{str(e)[:70]}")
                continue
            for w in chunk:
                if not r:
                    continue
                k = min(r, key=lambda x: abs(x - w))
                v = r[k]
                if abs(k - w) < 0.6 and v is not None and np.isfinite(v) and v != 0.0:
                    served += 1
        rate = served / max(len(sample), 1)
        print(f"\n  {lo:7.0f}-{hi:7.0f} A: catalogued {len(waves):5d}   "
              f"probed {len(sample):4d}   SERVED {served:4d}   rate {rate:.3f}")
        out.append(dict(lo=lo, hi=hi, catalogued=len(waves),
                        probed=len(sample), served=served, rate=rate))

    if out:
        d = pd.DataFrame(out)
        print("\n" + "=" * 70)
        print("  MPIA Engine-A reach for Fe " + a.ion + ", measured live:")
        print(d.to_string(index=False))
        opt = d[d.lo < 6900]
        ir = d[d.lo >= 6900]
        if len(opt) and len(ir):
            print(f"\n  optical rate {opt.rate.mean():.3f}  vs  IR rate {ir.rate.mean():.3f}")


if __name__ == "__main__":
    main()
