#!/usr/bin/env python3
"""Kitt Peak 1984 vs Kurucz 2005: product-level cross-check (RYA-938).

RYA-929 compared these two at window scale and got a clean GO, then compared them
line by line and got nonsense -- Kurucz showing almost no absorption at Al I 6696
and K I 7665/7699, where the 1984 atlas and the telluric-free IAG atlas agree on
deep lines.  A discrepancy that is invisible over a 300 A window and fatal over a
1 A window is the signature of a WAVELENGTH REGISTRATION error, not of a
correction.  RYA-929 said so and deferred it; this resolves it.

Three products, three different conventions, and the point of this script is that
none of them is assumed:

* **1984 Kitt Peak** -- NSO Atlas No. 1.  Its README (recovered under RYA-938;
  the staged copy was a saved HTTP 500 page) states the columns outright:
  wavelength in AIR (nm), pseudo-residual flux, irradiance in uW/cm^2/nm.
  Unity IS the continuum, so nothing is re-fitted (RYA-911).
* **Kurucz 2005** -- `irradthu.dat`, absolute irradiance in W/m^2/nm on a nm grid
  whose MEDIUM THE HEADER NEVER STATES.  It is measured here, not assumed.
* **IAG/Baker 2020** -- telluric-free, vacuum wavenumber, converted to air.  It is
  the independent third opinion: a different telescope and a different correction
  method, so agreement with it is not self-confirmation.

The medium test is deliberately blunt: cross-correlate each product against the
1984 atlas over telluric-free windows and read off the lag.  Air-vs-vacuum at
6700 A is ~1.85 A, roughly 200 sampled pixels -- if that is the defect it cannot
hide, and if it is not, the lag comes back at zero.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits
from scipy.ndimage import maximum_filter1d, uniform_filter1d

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from pipeline.wavelength_util import vac_to_air, air_to_vac  # noqa: E402
from pipeline.telluric_policy import TELLURIC_BANDS  # noqa: E402

#: Telluric-free windows with strong, well-separated solar lines. Chosen for the
#: registration test only: they must contain structure to correlate, and no
#: terrestrial absorption that could differ between a corrected and an
#: uncorrected product and bias the lag.
REGISTRATION_WINDOWS = ((5160.0, 5200.0), (6100.0, 6160.0), (6600.0, 6650.0))
#: Clean controls that no registered telluric band touches.
CLEAN_CONTROLS = ((6600.0, 6650.0), (6690.0, 6702.0), (5160.0, 5200.0))
#: The lines RYA-929 flagged and could not explain, plus its own clean controls.
DIAGNOSTIC_LINES = (
    ("Al I 6631", 6631.218), ("Al I 6696", 6696.185),
    ("Fe I 6872", 6872.1617), ("Fe I 6875", 6875.44503),
    ("Si I 6876", 6876.35917), ("Fe I 6881", 6881.44175),
    ("K I 7665", 7664.899), ("K I 7699", 7698.964), ("N I 8216", 8216.336),
)
#: Continuum envelope width for products that ship no continuum. Much wider than
#: any photospheric line, narrow enough to track the instrument response.
ENVELOPE_A = 8.0
CONTINUUM_PERCENTILE = 95.0


# ── readers: one per product, each stating its own conventions ───────────────

def read_kp1984(directory: Path, lo: float, hi: float):
    """(air_A, residual_flux). Unity IS the continuum -- nothing is re-fitted."""
    w, f = [], []
    for path in sorted(Path(directory).glob("lm[0-9]*")):
        stem = path.name[2:]
        if not stem.isdigit():
            continue
        start_a = float(stem) * 10.0
        if start_a > hi + 20.0 or start_a + 41.0 < lo - 20.0:
            continue
        try:
            d = np.atleast_2d(np.loadtxt(path))
        except Exception:
            continue                      # integrity is kp_atlas_integrity's job
        wa = d[:, 0] * 10.0
        m = (wa >= lo) & (wa <= hi)
        w.append(wa[m]); f.append(d[m, 1])
    if not w:
        return np.array([]), np.array([])
    w = np.concatenate(w); f = np.concatenate(f)
    o = np.argsort(w)
    w, f = w[o], f[o]
    keep = np.r_[True, np.diff(w) > 0]     # segments overlap 0.1 nm at the red end
    return w[keep], f[keep]


def read_kurucz2005(path: Path, lo: float, hi: float, assume_vacuum: bool):
    """(air_A, irradiance). `assume_vacuum` is the hypothesis under test."""
    w, f = [], []
    with Path(path).open(errors="replace") as fh:
        for line in fh:
            q = line.split()
            if len(q) < 2:
                continue
            try:
                x, y = float(q[0]) * 10.0, float(q[1])
            except ValueError:
                continue                   # header lines
            if lo - 5.0 <= x <= hi + 5.0:
                w.append(x); f.append(y)
    w, f = np.asarray(w), np.asarray(f)
    if not w.size:
        return w, f
    if assume_vacuum:
        w = np.asarray(vac_to_air(w))
    m = (w >= lo) & (w <= hi)
    return w[m], f[m]


def read_iag(path: Path, lo: float, hi: float):
    """(air_A, normalised flux) from vacuum wavenumbers."""
    with fits.open(path, memmap=True) as hdul:
        d = hdul[1].data
        wn = np.asarray(d["v"], float)
        fl = np.asarray(d["s"], float)
    w = np.asarray(vac_to_air(1.0e8 / wn))
    o = np.argsort(w)
    w, fl = w[o], fl[o]
    m = (w >= lo) & (w <= hi) & np.isfinite(fl)
    return w[m], fl[m]


def envelope_normalise(w, f):
    """Continuum for a product that ships none: a high running envelope.

    Stated plainly because it is a DIAGNOSTIC normalisation, not a production
    continuum, and the difference matters (RYA-911). It is only ever applied to
    Kurucz irradiance; the 1984 atlas and IAG already carry theirs.
    """
    if w.size < 8:
        return f
    step = float(np.median(np.diff(w)))
    width = max(3, int(round(ENVELOPE_A / step)) | 1)
    env = maximum_filter1d(f, size=width, mode="nearest")
    env = uniform_filter1d(env, size=width, mode="nearest")
    env = np.where(env > 0, env, np.nan)
    return f / env


# ── the registration measurement ─────────────────────────────────────────────

def cross_correlate_shift(w_ref, f_ref, w_test, f_test, max_shift_A=3.0,
                          step_A=0.005):
    """Lag in Angstrom that best aligns `test` onto `ref`, by absorption pattern.

    Both are reduced to (1 - flux) so only line structure drives the match, which
    keeps a difference in continuum treatment or flux units out of the answer.
    """
    if w_ref.size < 50 or w_test.size < 50:
        return float("nan"), float("nan")
    lo = max(w_ref.min(), w_test.min()) + max_shift_A
    hi = min(w_ref.max(), w_test.max()) - max_shift_A
    if hi - lo < 5.0:
        return float("nan"), float("nan")
    grid = np.arange(lo, hi, step_A)
    a = 1.0 - np.interp(grid, w_ref, f_ref)
    a = a - np.nanmean(a)
    shifts = np.arange(-max_shift_A, max_shift_A + step_A, step_A)
    best, best_lag = -np.inf, float("nan")
    scores = []
    for s in shifts:
        b = 1.0 - np.interp(grid + s, w_test, f_test)
        b = b - np.nanmean(b)
        denom = np.sqrt(np.nansum(a * a) * np.nansum(b * b))
        c = float(np.nansum(a * b) / denom) if denom > 0 else np.nan
        scores.append(c)
        if np.isfinite(c) and c > best:
            best, best_lag = c, float(s)
    return best_lag, best


def line_depth(w, f, centre, half=0.35):
    """1 - min(flux) in a tight window. The window is deliberately tight: a
    registration error must show up as a MISSING line, not be papered over."""
    m = (w >= centre - half) & (w <= centre + half) & np.isfinite(f)
    if m.sum() < 5:
        return {"n": int(m.sum()), "depth": float("nan"), "min_wave_A": float("nan")}
    x, y = w[m], f[m]
    j = int(np.nanargmin(y))
    return {"n": int(m.sum()), "depth": float(1.0 - y[j]),
            "min_wave_A": float(x[j])}


def band_metrics(w, f, lo, hi):
    """The RYA-805/905 depth instrument, so numbers compare to HARPS/IAG/Elgueta."""
    m = (w >= lo) & (w <= hi) & np.isfinite(f)
    if m.sum() < 10:
        return {"n": int(m.sum()), "min": None, "median": None,
                "pct_below_0.5": None, "pct_below_0.95": None}
    y = f[m]
    return {"n": int(m.sum()), "min": float(np.nanmin(y)),
            "median": float(np.nanmedian(y)),
            "pct_below_0.5": float(100 * np.mean(y < 0.5)),
            "pct_below_0.95": float(100 * np.mean(y < 0.95))}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--kp1984", type=Path, required=True)
    ap.add_argument("--kurucz", type=Path, required=True)
    ap.add_argument("--iag", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    report = {"ticket": "RYA-938", "products": {
        "kp1984": {"path": str(args.kp1984),
                   "wavelength_medium": "air",
                   "medium_basis": ("NSO Atlas No. 1 README, recovered under RYA-938: "
                                    "'the first column contains the wavelength in air'"),
                   "flux": "pseudo-residual flux; unity IS the continuum",
                   "continuum_treatment": "none applied -- the product ships it"},
        "kurucz2005": {"path": str(args.kurucz),
                       "wavelength_medium": "UNDER TEST",
                       "medium_basis": "header states only 'nm'; measured below",
                       "flux": "absolute irradiance W/m^2/nm",
                       "continuum_treatment": (f"diagnostic {CONTINUUM_PERCENTILE:.0f}th "
                                               f"running envelope, {ENVELOPE_A} A")},
        "iag_baker2020": {"path": str(args.iag),
                          "wavelength_medium": "vacuum wavenumber -> air",
                          "flux": "normalised, telluric-free",
                          "continuum_treatment": "none applied"}}}

    # ── 1. wavelength medium, measured ───────────────────────────────────────
    registration = []
    for lo, hi in REGISTRATION_WINDOWS:
        wk, fk = read_kp1984(args.kp1984, lo, hi)
        wi, fi = read_iag(args.iag, lo, hi)
        row = {"window_A": [lo, hi], "kp1984_n": int(wk.size)}
        lag_i, c_i = cross_correlate_shift(wk, fk, wi, fi)
        row["iag_vs_kp1984_lag_A"] = lag_i
        row["iag_vs_kp1984_corr"] = c_i
        for label, vac in (("as_read_nm", False), ("vac_to_air", True)):
            wu, fu = read_kurucz2005(args.kurucz, lo, hi, assume_vacuum=vac)
            fu = envelope_normalise(wu, fu)
            lag, c = cross_correlate_shift(wk, fk, wu, fu)
            row[f"kurucz_{label}_lag_A"] = lag
            row[f"kurucz_{label}_corr"] = c
            row[f"kurucz_{label}_n"] = int(wu.size)
        mid = 0.5 * (lo + hi)
        row["air_vac_offset_at_mid_A"] = float(air_to_vac(np.array([mid]))[0] - mid)
        registration.append(row)
    report["registration"] = registration

    lags = [r["kurucz_as_read_nm_lag_A"] for r in registration
            if np.isfinite(r["kurucz_as_read_nm_lag_A"])]
    offsets = [r["air_vac_offset_at_mid_A"] for r in registration]
    median_lag = float(np.median(lags)) if lags else float("nan")
    median_offset = float(np.median(offsets))
    # Declared BEFORE measuring: the lag either sits at zero (air) or at the
    # air-vacuum offset (vacuum). Half the offset separates the two hypotheses.
    decision_threshold = 0.5 * median_offset
    verdict = ("vacuum" if np.isfinite(median_lag) and median_lag > decision_threshold
               else "air" if np.isfinite(median_lag) and abs(median_lag) <= decision_threshold
               else "UNDETERMINED")
    report["medium_verdict"] = {
        "kurucz2005": verdict,
        "median_lag_as_read_A": median_lag,
        "median_air_vac_offset_A": median_offset,
        "decision_threshold_A": decision_threshold,
        "rule": ("lag ~ 0 => the column is already air; lag ~ +air-vac offset => the "
                 "column is vacuum and must be converted. Threshold declared before "
                 "measurement as half the offset."),
    }
    use_vacuum = verdict == "vacuum"
    report["products"]["kurucz2005"]["wavelength_medium"] = verdict

    # ── 2. telluric state, on the same instrument as HARPS/IAG/Elgueta ───────
    bands = [(lo, hi, name) for lo, hi, name in TELLURIC_BANDS]
    rows = []
    for lo, hi, name in bands + [(a, b, "clean control") for a, b in CLEAN_CONTROLS]:
        wk, fk = read_kp1984(args.kp1984, lo, hi)
        wu, fu = read_kurucz2005(args.kurucz, lo, hi, assume_vacuum=use_vacuum)
        fu = envelope_normalise(wu, fu)
        wi, fi = read_iag(args.iag, lo, hi)
        rows.append({"band": name, "lo_A": lo, "hi_A": hi,
                     "kp1984": band_metrics(wk, fk, lo, hi),
                     "kurucz2005": band_metrics(wu, fu, lo, hi),
                     "iag": band_metrics(wi, fi, lo, hi)})
    report["telluric_bands"] = rows

    # ── 3. the lines RYA-929 could not explain ───────────────────────────────
    lines = []
    for name, centre in DIAGNOSTIC_LINES:
        lo, hi = centre - 6.0, centre + 6.0
        wk, fk = read_kp1984(args.kp1984, lo, hi)
        wi, fi = read_iag(args.iag, lo, hi)
        entry = {"line": name, "air_A": centre,
                 "kp1984": line_depth(wk, fk, centre),
                 "iag": line_depth(wi, fi, centre)}
        for label, vac in (("as_read_nm", False), ("vac_to_air", True)):
            wu, fu = read_kurucz2005(args.kurucz, lo, hi, assume_vacuum=vac)
            fu = envelope_normalise(wu, fu)
            entry[f"kurucz_{label}"] = line_depth(wu, fu, centre)
        lines.append(entry)
    report["diagnostic_lines"] = lines

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report["medium_verdict"], indent=2))


if __name__ == "__main__":
    main()
