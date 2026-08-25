#!/usr/bin/env python3
"""RYA-1041 Step 1 — the two diagnostics the five thresholds must be SET from.

Produces evidence, sets nothing. The numbers are Ryan's (RYA-1041 gate), and the
firewall (RYA-161) is that a threshold is derived from our own data or physics, never
from proximity to a reference and never borrowed.

D1. REW and red_chi2 distributions PER BAND on the current corrected pool.
    Tests for a NATURAL BREAK rather than assuming one: the largest ratio-gap between
    consecutive order statistics in the central 90%, reported against the spread. A
    percentile is a QUOTA, not a derivation, so no percentile is proposed as a cut.

D2. The Fe anchor's DIFFERENTIAL line-to-line scatter — the SETTING input for
    admission_tol_dex (RYA-898/1041), and NOT the solo A(Fe) scatter.

    🔴 WHY SOLO SCATTER IS THE WRONG INPUT. Every published value is
    [X/H] = A(X)star - A(X)sun, and it works because the systematics CANCEL: same
    instrument, resolution, continuum treatment, EW method, line list on both sides.
    Solo scatter measures what does NOT cancel plus what does. Sizing an admission
    tolerance from it would admit lines on a spread the differential never actually
    pays. Every pool here is solar, so the differential available is SAME LINE, TWO
    ARMS: d_i = A_i(arm1) - A_i(arm2), paired per line.

    Two flavours, and the distinction is the whole point:
      * SAME-INSTRUMENT  (KP1984-molecfit vs KP2005) — one telescope, one atlas
        lineage, differing only in reduction. Closest thing we hold to a true
        cancelling chain: this is the floor.
      * CROSS-INSTRUMENT (HARPS vs KP, IAG vs KP) — what does NOT cancel when the
        chains differ, which is exactly the case RYA-898 says must not be used as a
        differential baseline.
"""
from __future__ import annotations

import glob
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
BP = ROOT / "data" / "results" / "band_products"
OUT = ROOT / "data" / "audit" / "rya1041"

BANDS = {"near-UV": (3000, 3780), "VIS": (4200, 6910),
         "red-optical": (6910, 9199), "IR": (9199, 12976)}


def _load(pattern: str) -> pd.DataFrame | None:
    frames = []
    for f in sorted(glob.glob(str(BP / pattern))):
        d = pd.read_csv(f)
        d["_src"] = Path(f).name
        # holding is the product identity, and it is in the stem (RYA-933/934)
        for h in ("solar_harps_molecfit_corrected", "solar_kpno_molecfit_corrected",
                  "solar_kpno_kurucz2005_corrected", "solar_iag"):
            if h in Path(f).name:
                d["_holding"] = h
                break
        else:
            d["_holding"] = "unknown"
        frames.append(d)
    return pd.concat(frames, ignore_index=True) if frames else None


def _keep(d: pd.DataFrame) -> pd.DataFrame:
    return d[d["in_aggregate"] == True] if "in_aggregate" in d.columns else d


def largest_gap(v: np.ndarray, *, log: bool = False) -> dict:
    """The biggest jump between consecutive sorted values, in the central 90%.

    Reported as a RATIO to the median step, so "is there a break" has a number rather
    than an eyeball. Edges are trimmed because the largest gap in any sample is almost
    always at the tail, which says nothing about structure.

    🔴 TEST A HEAVY-TAILED POSITIVE QUANTITY IN LOG SPACE. red_chi2 spans 2.5 to 3404
    here, and on raw values the step size grows WITH the value, so the largest gap sits
    near the top by construction and every band reports "structure" that is an artifact
    of scale. My first version did exactly that and flagged 5 of 6 distributions as
    structured. Log-spacing makes the steps comparable, so a gap means a gap.
    """
    v = np.asarray(v, float)
    if log:
        v = v[v > 0]
        v = np.log10(v)
    v = np.sort(v)
    v = v[(v >= np.quantile(v, 0.05)) & (v <= np.quantile(v, 0.95))]
    if len(v) < 8:
        return {"n": int(len(v)), "verdict": "too few points to test"}
    steps = np.diff(v)
    med = float(np.median(steps)) or 1e-12
    i = int(np.argmax(steps))
    return {"n": int(len(v)), "largest_gap": float(steps[i]),
            "gap_at": float(v[i]), "median_step": med,
            "gap_over_median_step": float(steps[i] / med),
            "range_p5_p95": [float(v[0]), float(v[-1])]}


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    report: dict = {"ticket": "RYA-1041", "sets_nothing": True, "d1": {}, "d2": {}}

    # ── D1 ────────────────────────────────────────────────────────────────────────
    print("=" * 78)
    print("D1  REW and red_chi2 per band, on the CURRENT corrected pool")
    print("=" * 78)
    for band, (lo, hi) in BANDS.items():
        for route, pat in (("SYNTH", f"FeI_{lo}_{hi}_*_SYNTH_*1D-LTE_lines.csv"),
                           ("PROFILEFIT", f"FeI_{lo}_{hi}_*_PROFILEFIT_*1D-LTE_lines.csv")):
            d = _load(pat)
            if d is None:
                continue
            k = _keep(d)
            if not len(k):
                continue
            entry = {"n": int(len(k)), "arms": sorted(k["_holding"].unique().tolist())}
            for col in ("rew", "red_chi2"):
                if col not in k.columns:
                    continue
                v = pd.to_numeric(k[col], errors="coerce").dropna().values
                if len(v) < 8:
                    continue
                q = {p: float(np.quantile(v, p / 100)) for p in (0, 5, 50, 95, 100)}
                use_log = col == "red_chi2"          # positive, spans 3+ decades
                bt = largest_gap(v, log=use_log)
                bt["tested_in"] = "log10" if use_log else "linear"
                entry[col] = {"quantiles": q, "break_test": bt}
                gt = entry[col]["break_test"].get("gap_over_median_step")
                print(f"\n  {band:<12} {route:<11} {col:<9} n={len(v):4d}  "
                      f"min={q[0]:.4g}  med={q[50]:.4g}  max={q[100]:.4g}")
                if gt is not None:
                    _b = entry[col]["break_test"]
                    _at = 10 ** _b["gap_at"] if _b["tested_in"] == "log10" else _b["gap_at"]
                    print(f"      [{_b['tested_in']}] largest central gap "
                          f"at {_at:.4g}  "
                          f"= {gt:.1f}x the median step  ->  "
                          f"{'STRUCTURE?' if gt > 10 else 'NO NATURAL BREAK (continuous)'}")
            report["d1"][f"{band}|{route}"] = entry

    # ── D2 ────────────────────────────────────────────────────────────────────────
    print("\n" + "=" * 78)
    print("D2  Fe anchor DIFFERENTIAL line-to-line scatter  (admission_tol_dex input)")
    print("=" * 78)
    d = _load("FeI_4200_6910_*_SYNTH_GRADED_1D-LTE_lines.csv")
    k = _keep(d)
    k = k.assign(w=pd.to_numeric(k["wavelength_air_A"], errors="coerce"),
                 A=pd.to_numeric(k["abundance"], errors="coerce")).dropna(subset=["w", "A"])
    wide = k.pivot_table(index="w", columns="_holding", values="A", aggfunc="first")
    PAIRS = [("solar_kpno_molecfit_corrected", "solar_kpno_kurucz2005_corrected", "SAME-INSTRUMENT"),
             ("solar_harps_molecfit_corrected", "solar_kpno_molecfit_corrected", "CROSS-INSTRUMENT"),
             ("solar_iag", "solar_kpno_molecfit_corrected", "CROSS-INSTRUMENT"),
             ("solar_harps_molecfit_corrected", "solar_iag", "CROSS-INSTRUMENT")]
    for a, b, kind in PAIRS:
        if a not in wide.columns or b not in wide.columns:
            continue
        pair = wide[[a, b]].dropna()
        if len(pair) < 5:
            continue
        dif = (pair[a] - pair[b]).values
        mad = float(np.median(np.abs(dif - np.median(dif))))
        rec = {"kind": kind, "n_lines_paired": int(len(pair)),
               "median_offset_dex": float(np.median(dif)),
               "std_dex": float(np.std(dif, ddof=1)),
               "MAD_dex": mad, "robust_sigma_dex": float(1.4826 * mad)}
        report["d2"][f"{a} - {b}"] = rec
        print(f"\n  {kind:<18} {a.replace('solar_','')} - {b.replace('solar_','')}")
        print(f"      paired lines {rec['n_lines_paired']:3d}   median offset "
              f"{rec['median_offset_dex']:+.4f}   std {rec['std_dex']:.4f}   "
              f"robust sigma (1.4826*MAD) {rec['robust_sigma_dex']:.4f} dex")

    # solo scatter, reported ONLY to show it is the larger, wrong number
    solo = float(1.4826 * np.median(np.abs(k["A"] - k["A"].median())))
    report["d2"]["_solo_scatter_dex_WRONG_INPUT"] = solo
    print(f"\n  For contrast, SOLO A(Fe) scatter (NOT the setting input): {solo:.4f} dex")

    (OUT / "rya1041_threshold_diagnostics.json").write_text(json.dumps(report, indent=2) + "\n")
    print(f"\nwrote {OUT / 'rya1041_threshold_diagnostics.json'}")
    print("\nNo threshold is set by this script. RYA-1041 Step 2 is Ryan's.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
