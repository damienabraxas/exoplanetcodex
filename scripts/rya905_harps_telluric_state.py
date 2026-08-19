#!/usr/bin/env python3
"""RYA-905 — determine `telluric_applied` for the direct-solar HARPS holding, from FLUX.

RYA-806 registered `solar_harps` with `telluric_applied=unknown` and said why: *"product
files were not reachable at backfill time and no prior ticket measured this holding's
telluric state."* It also said how to resolve it — determine it, then record it. Never
default it: assuming `applied` fabricates a correction, assuming `not-applied` risks
correcting an already-corrected product twice (RYA-786).

`telluric_intake.from_headers()` cannot answer here — the committed HARPS product is a CSV
with no FITS headers. So this uses the other instrument RYA-805 established: measure the
DEPTH in a window where the terrestrial atmosphere absorbs and the Sun does not, against a
control window where neither does. Same metric as RYA-805 (`pct_below_0.5`), so the numbers
are comparable to the CRIRES+, Kitt Peak and IAG figures already on the record.

🔴 THE CONTROL IS THE POINT. A band that looks absorbed proves nothing on its own — deep
solar lines are also absorbed. What discriminates is a clean window measured the same way:
if the O2 band bottoms out near zero while a window 60 A away never drops below ~0.57, the
absorber is terrestrial.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from config.constants import PATHS  # noqa: E402

#: O2 A-band. Terrestrial only — there is no solar O2 counterpart, which is what makes it
#: a clean discriminator. Core/wing split per pipeline.telluric_policy's own band edges.
WINDOWS = [
    ("O2 A-band core", 6867.0, 6884.0, "telluric"),
    ("O2 A-band wing", 6884.0, 6910.0, "telluric"),
    ("control 6800-6860", 6800.0, 6860.0, "control"),
    ("control 5000-5100", 5000.0, 5100.0, "control"),
]

#: Already on the record, same metric, for scale (RYA-783/786/794/805).
REFERENCE_SCALE = {
    "kpno_solar_atlas (Kurucz 1984), O2 A-band": 51.3,
    "iag_fts_solar_atlas (Baker 2020), O2 A-band": 0.1,
    "elgueta2026 Sun_Y_rv (corrected), O2 A-band": 0.10,
}


def stats(w, f, lo, hi):
    m = (w >= lo) & (w <= hi)
    n = int(m.sum())
    if not n:
        return None
    y = f[m]
    return {"n_px": n,
            "median": round(float(np.median(y)), 4),
            "min": round(float(np.min(y)), 4),
            "pct_below_0.7": round(float((y < 0.7).mean() * 100), 2),
            "pct_below_0.5": round(float((y < 0.5).mean() * 100), 2)}


def main() -> int:
    path = Path(str(PATHS["solar_normalized"]))
    d = pd.read_csv(path)
    w = d["wavelength_air_A"].to_numpy(float)
    f = d["flux_normalized"].to_numpy(float)
    ok = np.isfinite(w) & np.isfinite(f)
    w, f = w[ok], f[ok]

    rows = {}
    for label, lo, hi, kind in WINDOWS:
        s = stats(w, f, lo, hi)
        if s:
            s["kind"] = kind
            rows[label] = s

    tell = max(r["pct_below_0.5"] for r in rows.values() if r["kind"] == "telluric")
    ctrl_min = min(r["min"] for r in rows.values() if r["kind"] == "control")
    tell_min = min(r["min"] for r in rows.values() if r["kind"] == "telluric")

    # A corrected product measures ~0.1% here. An uncorrected one measures percent-level
    # or worse. The verdict is the CONTRAST, not the absolute.
    verdict = "not-applied" if (tell >= 1.0 and tell_min < 0.1 <= ctrl_min) else "INCONCLUSIVE"

    out = {
        "ticket": "RYA-905",
        "holding_id": "solar_harps",
        "instrument": "harps",
        "source": str(path),
        "method": ("flux depth in a terrestrial-only band against a clean control, "
                   "same metric as RYA-805 (pct_below_0.5)"),
        "windows": rows,
        "discriminator": {
            "worst_telluric_pct_below_0.5": tell,
            "min_flux_in_telluric_window": tell_min,
            "min_flux_in_control_window": ctrl_min,
            "reads": (f"flux reaches {tell_min} inside the O2 A-band while never dropping "
                      f"below {ctrl_min} in a clean window 60 A away"),
        },
        "reference_scale_same_metric": REFERENCE_SCALE,
        "verdict_telluric_applied": verdict,
        "consequence": ("None for the HARPS band products. `harps` is registered "
                        "telluric_required=no, so gate_holding() permits it and tellurics "
                        "are handled by per-line clean-line selection over the enumerated "
                        "bands (RYA-460/786) — which is what quarantined the 5 O2 B-band "
                        "lines in the RYA-897 run, with the reason stated. This records "
                        "the state; it does not change the method."),
    }
    outdir = ROOT / "data" / "audit" / "telluric_state"
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "rya905_harps_telluric_state.json").write_text(json.dumps(out, indent=2))

    print("=== RYA-905 — solar_harps telluric state, measured from flux")
    for k, v in rows.items():
        print(f"  {k:20} n={v['n_px']:6}  median {v['median']:.4f}  min {v['min']:.4f}  "
              f"%<0.5 {v['pct_below_0.5']:6.2f}  [{v['kind']}]")
    print(f"\n  discriminator: {out['discriminator']['reads']}")
    print(f"  VERDICT telluric_applied = {verdict}")
    print(f"\n[out] {outdir / 'rya905_harps_telluric_state.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
