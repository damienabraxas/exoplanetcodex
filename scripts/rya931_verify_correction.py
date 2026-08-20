#!/usr/bin/env python3
"""RYA-931 verification: O2 B gate, clean-line preservation, IAG discriminator.

Every number here is measured from the products; nothing is tuned to pass.  The
three checks answer three different questions, and a product must clear all of
them:

1. **Did the terrestrial absorption come out?**  O2 B 6867-6884 A residual
   statistics before and after, against a clean control window that the
   correction never touches.
2. **Did the Sun survive?**  Correction must not flatten or manufacture solar
   structure.  Outside the fitted region this is an exact-equality invariant;
   inside it, clean-window solar line depths must be preserved.
3. **Does it agree with an independent instrument?**  The Baker+2020 IAG
   telluric-free atlas is a different telescope and a different correction
   method, so agreement is not self-confirmation.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits
from scipy.ndimage import gaussian_filter1d

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from pipeline.wavelength_util import vac_to_air  # noqa: E402

O2_B = (6867.0, 6884.0)
CONTROL = (6800.0, 6860.0)
FITTED = (6857.0, 6911.0)
# Solar-line preservation tolerance.  Derived, not chosen: the correction can
# only perturb a line by the transmission model's own accuracy, so a clean-window
# line (transmission >= 0.98) may move by at most dT/T ~ 1% in relative depth.
LINE_DEPTH_TOL = 0.01
# The continuum spline may move outside the fitted region only by the amount the
# recovered O2 B flux moves its own knot anchor.  1% bounds that; the measured
# value is reported alongside so the bound is never mistaken for the result.
CONTINUUM_SHIFT_TOL = 0.01
# The clean control window is renormalised like everything else, so it is checked
# for a bounded shift, not for bit equality.
CONTROL_SHIFT_TOL = 0.01


def load_csv(path: Path):
    df = pd.read_csv(path, comment="#")
    return (df["wavelength_air_A"].to_numpy(float), df["flux_normalized"].to_numpy(float),
            df["flux_raw"].to_numpy(float), df["continuum"].to_numpy(float))


def band_stats(w: np.ndarray, f: np.ndarray, lo: float, hi: float) -> dict:
    k = (w >= lo) & (w <= hi)
    v = f[k]
    finite = v[np.isfinite(v)]
    return {"window_A": [lo, hi], "n_pixels": int(k.sum()),
            "n_finite": int(finite.size),
            "min": float(finite.min()) if finite.size else None,
            "median": float(np.median(finite)) if finite.size else None,
            "pct_below_0.5": float(100 * np.mean(finite < 0.5)) if finite.size else None,
            "pct_below_0.95": float(100 * np.mean(finite < 0.95)) if finite.size else None}


def iag_reference(atlas: Path, w: np.ndarray) -> np.ndarray:
    with fits.open(atlas) as hdul:
        d = hdul[1].data
        aw = vac_to_air(1.0e8 / np.asarray(d["v"], float))
        af = np.asarray(d["s"], float)
    o = np.argsort(aw)
    aw, af = aw[o], af[o]
    # Degrade the R~1e6 atlas to the HARPS resolution before comparing.
    step = float(np.median(np.diff(aw)))
    fwhm_a = 6875.0 / 101_000.0
    af = gaussian_filter1d(af, (fwhm_a / 2.3548) / step)
    return np.interp(w, aw, af, left=np.nan, right=np.nan)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("uncorrected_csv", type=Path)
    ap.add_argument("corrected_csv", type=Path)
    ap.add_argument("atlas", type=Path)
    ap.add_argument("out_json", type=Path)
    ap.add_argument("--plot", type=Path, default=None)
    args = ap.parse_args()

    wu, fu, ru, cu = load_csv(args.uncorrected_csv)
    wc, fc, rc, cc = load_csv(args.corrected_csv)
    if wu.shape != wc.shape or not np.allclose(wu, wc, atol=1e-9):
        raise SystemExit("corrected and uncorrected spectra are not on the same grid")

    # ── 1. O2 B before/after, with an untouched control ──────────────────────
    gate = {"o2b_before": band_stats(wu, fu, *O2_B), "o2b_after": band_stats(wc, fc, *O2_B),
            "control_before": band_stats(wu, fu, *CONTROL),
            "control_after": band_stats(wc, fc, *CONTROL)}

    # ── 2a. containment: the correction touches only the fitted region ───────
    # The exact-equality invariant lives on the RAW co-add.  Normalised flux
    # cannot carry it, because continuum fitting is global: removing the O2 band
    # raises the 95th-percentile anchor of the knot window that contains it, so
    # the spline -- and therefore every normalised pixel -- shifts slightly.
    # That shift is a real, wanted consequence (the old anchor was depressed by
    # telluric absorption), so it is measured and reported, not asserted away.
    outside = (wu < FITTED[0]) | (wu > FITTED[1])
    raw_identical = bool(np.array_equal(ru[outside], rc[outside]))
    finite_out = outside & np.isfinite(fu) & np.isfinite(fc)
    max_norm_shift = float(np.max(np.abs(fc[finite_out] - fu[finite_out])))
    max_continuum_shift = float(np.max(np.abs(cc[finite_out] - cu[finite_out])
                                       / cu[finite_out]))

    # ── 2b. clean-window solar-line preservation inside the fitted region ────
    inside = (wu >= FITTED[0]) & (wu <= FITTED[1]) & np.isfinite(fc) & np.isfinite(fu)
    solar = iag_reference(args.atlas, wu)
    # "clean" = the IAG atlas says there IS a solar line here and the uncorrected
    # HARPS flux is close to it, i.e. no meaningful telluric absorption.
    clean = inside & np.isfinite(solar) & (np.abs(fu - solar) < 0.02)
    line = clean & (solar < 0.9)
    depth_shift = float(np.max(np.abs(fc[line] - fu[line]))) if line.any() else 0.0

    # ── 3. independent-instrument agreement in the O2 B band ─────────────────
    k = (wu >= O2_B[0]) & (wu <= O2_B[1]) & np.isfinite(solar)
    ku, kc = k & np.isfinite(fu), k & np.isfinite(fc)
    iag = {"n_compared_before": int(ku.sum()), "n_compared_after": int(kc.sum()),
           "rms_vs_iag_before": float(np.sqrt(np.mean((fu[ku] - solar[ku]) ** 2))),
           "rms_vs_iag_after": float(np.sqrt(np.mean((fc[kc] - solar[kc]) ** 2))),
           "median_offset_before": float(np.median(fu[ku] - solar[ku])),
           "median_offset_after": float(np.median(fc[kc] - solar[kc]))}

    report = {
        "ticket": "RYA-931",
        "uncorrected": str(args.uncorrected_csv), "corrected": str(args.corrected_csv),
        "o2b_gate": gate,
        "solar_preservation": {
            "raw_coadd_bit_identical_outside_fitted_region": raw_identical,
            "max_normalised_flux_shift_outside": max_norm_shift,
            "max_fractional_continuum_shift_outside": max_continuum_shift,
            "continuum_shift_tolerance": CONTINUUM_SHIFT_TOL,
            "n_clean_solar_line_pixels": int(line.sum()),
            "max_abs_depth_shift": depth_shift,
            "tolerance": LINE_DEPTH_TOL,
            "pass": bool(raw_identical and line.sum() > 0
                         and depth_shift <= LINE_DEPTH_TOL
                         and max_continuum_shift <= CONTINUUM_SHIFT_TOL),
        },
        "iag_discriminator": iag,
    }
    report["verdict"] = {
        "o2b_min_improved": gate["o2b_after"]["min"] > gate["o2b_before"]["min"],
        "o2b_saturation_cleared": gate["o2b_after"]["pct_below_0.5"] == 0.0,
        "control_shift_bounded": bool(
            abs(gate["control_after"]["median"] - gate["control_before"]["median"])
            <= CONTROL_SHIFT_TOL
            and gate["control_after"]["pct_below_0.5"]
            == gate["control_before"]["pct_below_0.5"]
            and gate["control_after"]["n_finite"] == gate["control_before"]["n_finite"]),
        "solar_preserved": report["solar_preservation"]["pass"],
        "iag_agreement_improved": iag["rms_vs_iag_after"] < iag["rms_vs_iag_before"],
    }
    report["verdict"]["pass"] = all(report["verdict"].values())
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report["verdict"], indent=2))
    print(json.dumps(gate, indent=2))
    print(json.dumps(iag, indent=2))

    if args.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(2, 1, figsize=(12, 7), sharex=True)
        m = (wu >= 6855) & (wu <= 6900)
        axes[0].plot(wu[m], fu[m], lw=0.7, color="#b03030", label="uncorrected solar_harps")
        axes[0].plot(wu[m], solar[m], lw=0.7, color="#7a7a7a", label="IAG telluric-free (Baker 2020)")
        axes[1].plot(wc[m], fc[m], lw=0.7, color="#20609a", label="RYA-931 corrected")
        axes[1].plot(wu[m], solar[m], lw=0.7, color="#7a7a7a", label="IAG telluric-free (Baker 2020)")
        for ax in axes:
            ax.axvspan(*O2_B, color="#d8b400", alpha=0.15, label="O2 B 6867-6884")
            ax.set_ylim(0, 1.15); ax.legend(fontsize=8, loc="lower left"); ax.grid(alpha=0.2)
        axes[1].set_xlabel("air wavelength (A)")
        fig.suptitle("RYA-931 HARPS solar telluric correction — O2 B band")
        fig.tight_layout()
        args.plot.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(args.plot, dpi=140)
        print(f"plot -> {args.plot}")


if __name__ == "__main__":
    main()
