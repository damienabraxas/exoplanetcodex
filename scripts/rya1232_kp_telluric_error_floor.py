#!/usr/bin/env python3
"""RYA-1232 -- measure the KP 1984 molecfit correction's per-pixel error against an
independent telluric-free reference, and apply RYA-940's own 5% budget with it.

    python3 scripts/rya1232_kp_telluric_error_floor.py            # measure + report only
    python3 scripts/rya1232_kp_telluric_error_floor.py --apply    # also quarantine in place

WHY. RYA-940 quarantines a corrected pixel when its relative error from the transmission
model would exceed CORRECTED_RELATIVE_ERROR_BUDGET = 0.05, estimating that error as
dT/T with dT the fit residual MAD over the band. That MAD is dominated by telluric-free
pixels (~0.001), so T_min came out ~0.2. Measured against the IAG telluric-free atlas the
correction is UNBIASED (median within 0.007 in every bin) but its error is not constant:
it is sigma(T) = k (1 - T) / T -- a fixed fraction k of the telluric line's depth -- with
k = 0.23-0.34 in the optical and 0.06-0.09 in the NIR. At T = 0.5 that is 20-40%, not 0.2%.

METHOD (per band, per 20 A chunk): both spectra normalised by their own p95 on the
chunk's telluric-free pixels (T > 0.97) -- which also cancels IAG's slow normalisation
defect beyond 10000 A -- then KP_corrected - reference on pixels the reference shows as
solar CONTINUUM (> 0.97), so solar-line registration and resolution never enter. k is the
median over T bins (0.50-0.95) of rms / ((1 - T) / T); the floor is the T where
sigma(T) = budget, i.e. T_floor = k / (k + budget).

REFERENCE: IAG (Reiners+2016 / Baker+2020) to 11075 A; beyond it the CRIRES+ J rest-frame
holding (another site, epoch and Doppler offset, so its telluric pixels do not coincide).
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from rya940_kp1984_correct import CORRECTED_RELATIVE_ERROR_BUDGET  # noqa: E402

KP_DIR = ROOT / "data/processed/kp1984_telluric_corrected"
OUT = ROOT / "data/results/rya1232_kp_telluric_error/error_law.json"
IAG_MAX_A = 11075.0
BIN_EDGES = (0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.95)
MIN_PER_BIN = 40
CHUNK_A = 20.0


def _reference(band_lo: float):
    if band_lo < IAG_MAX_A - 100:
        return "iag_fts_solar_atlas", "solar_iag"
    return "crires_plus", "solar_crires_plus_j_rya1219"


def measure(fn: Path) -> dict:
    import measure_band_ew as H
    d = np.loadtxt(fn)
    w, f, t = d[:, 0], d[:, 1], d[:, 2]
    lo, hi = map(float, fn.stem.split("_")[-2:])
    inst, hold = _reference(lo)
    top = min(hi, IAG_MAX_A) if inst.startswith("iag") else hi
    X, T = [], []
    for c in np.arange(lo + CHUNK_A / 2, top - CHUNK_A / 2, CHUNK_A):
        try:
            a = H.load_window_ex(inst, c, CHUNK_A / 2, holding=hold)
        except Exception:
            continue
        aw, af = np.asarray(a.wave, float), np.asarray(a.flux, float)
        g = np.isfinite(af)
        if g.sum() < 50:
            continue
        aw, af = aw[g], af[g]
        sel = ((w > max(c - CHUNK_A / 2 + 0.5, aw.min())) & (w < min(c + CHUNK_A / 2 - 0.5, aw.max()))
               & np.isfinite(f))
        clean = sel & (t > 0.97)
        if clean.sum() < 20:
            continue
        ref = np.interp(w, aw, af)
        ref = ref / np.percentile(ref[clean], 95)
        m = sel & (t >= BIN_EDGES[0]) & (t < BIN_EDGES[-1]) & (ref > 0.97)
        X.append(f[m] / np.percentile(f[clean], 95) - ref[m])
        T.append(t[m])
    X, T = np.concatenate(X), np.concatenate(T)
    bins, ks = [], []
    for a_, b_ in zip(BIN_EDGES[:-1], BIN_EDGES[1:]):
        m = (T >= a_) & (T < b_)
        if m.sum() < MIN_PER_BIN:
            continue
        x = X[m]
        rms = float(1.4826 * np.median(np.abs(x - np.median(x))))
        tm = float(np.median(T[m]))
        k = rms / ((1 - tm) / tm)
        ks.append(k)
        bins.append({"T_lo": a_, "T_hi": b_, "n": int(m.sum()), "median": round(float(np.median(x)), 4),
                     "rms": round(rms, 4), "k": round(k, 3)})
    k = float(np.median(ks))
    floor = k / (k + CORRECTED_RELATIVE_ERROR_BUDGET)
    return {"band": fn.stem.replace("kp1984_corrected_", ""), "reference_holding": hold,
            "k": round(k, 4), "T_floor": round(floor, 4), "bins": bins}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--apply", action="store_true")
    a = ap.parse_args()
    report = {"ticket": "RYA-1232", "budget": CORRECTED_RELATIVE_ERROR_BUDGET,
              "error_law": "sigma(T) = k (1 - T) / T", "bands": []}
    for fn in sorted(KP_DIR.glob("kp1984_corrected_*.txt")):
        if a.apply and "RYA-1232: pixels with T <" in fn.read_text()[:2000]:
            raise SystemExit(f"{fn.name} already carries the RYA-1232 floor; regenerate it with "
                             "rya940_kp1984_correct.py before re-applying (k must be measured "
                             "on the unquarantined correction)")
        lo = float(fn.stem.split("_")[-2])
        if lo < 5800:          # 3000-5800: 51 pixels below T 0.95 -- nothing to measure
            continue
        r = measure(fn)
        hdr = [ln for ln in fn.read_text().splitlines()[:5] if ln.startswith("#")]
        d = np.loadtxt(fn)
        before = float(np.mean(~np.isfinite(d[:, 1])))
        newq = np.isfinite(d[:, 1]) & (d[:, 2] < r["T_floor"])
        r["quarantined_frac_before"] = round(before, 4)
        r["quarantined_frac_after"] = round(before + float(newq.mean()), 4)
        report["bands"].append(r)
        print(f"{r['band']:12s} k={r['k']:.3f} floor T>={r['T_floor']:.3f}  "
              f"quarantined {before:.3f} -> {r['quarantined_frac_after']:.3f}")
        if a.apply and newq.any():
            d[newq, 1] = np.nan
            note = (f"# RYA-1232: pixels with T < {r['T_floor']:.4f} quarantined -- measured error "
                    f"k(1-T)/T, k={r['k']:.4f} vs {r['reference_holding']}, exceeds the 5% budget.")
            np.savetxt(fn, d, fmt=["%12.5f", "%14.6f", "%10.6f"],
                       header="\n".join(h.lstrip("# ") for h in hdr + [note]))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(report, indent=2) + "\n")
    print(f"wrote {OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
