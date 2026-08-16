#!/usr/bin/env python3
"""
RYA-846 — sigma_delta, the near-UV continuum-placement uncertainty
==================================================================
RYA-841 measured the LEVER: dA/d(delta) = 2.42 dex per unit fractional continuum error,
against the 1.0 the budget silently assumes. The honest near-UV continuum systematic is

    term = lever x sigma_delta

and sigma_delta -- how well the atlas continuum is actually placed in this band -- has
never been measured. This measures it.

🔴 THE OBVIOUS ROUTE IS DEAD, AND IT WOULD HAVE ANSWERED ZERO. The Kitt Peak atlas ships as
250 `lm####` segments and every seam OVERLAPS its neighbour by ~1.0 A -- 20 overlapping
seams inside 3000-3780 A. Two independent scans of the same wavelengths would measure the
atlas's own normalisation repeatability directly, with no model involved. They are not
independent: the overlapping pixels are on an IDENTICAL wavelength grid with a flux
difference of EXACTLY 0.000e+00. The atlas is one continuous spectrum split into files with
1 A of padding. Reporting sigma_delta from those seams would have produced a spuriously
perfect ~0, and hence a continuum term of ~0. Checked before use, not after.

WHAT IS MEASURED INSTEAD. A second, independent reduction of the same band exists on this
machine: Wallace, Hinkle, Livingston & Davis 2011 (NSO), a telluric-removed KPNO solar flux
atlas, fetched under RYA-485 for exactly this purpose -- its own PROVENANCE says it is the
"tie-breaker for the continuum lever: distinguishes a KPNO-family systematic from a
Kurucz-reduction-specific one".

    region 6   2899.5-3525.2 A   167,936 px   clean (no overflow markers)
    region 5   2995.4-4106.4 A   512,000 px   2,227 overflow markers in band, flux to +9.98

Region 6 is used; region 5 is reported as a cross-check where it is physically sane and NOT
used for the headline, because a normalised flux of 9.98 is not a measurement.

The difference between two independent normalisations of the same spectrum is an EMPIRICAL
LOWER BOUND on sigma_delta: it captures reduction-to-reduction disagreement, not distance
from an unobservable truth. That distinction is stated in the output, not buried here.

TWO CONTROLS, BOTH OF WHICH CAN FAIL LOUDLY
  1. WAVELENGTH CONVENTION. Wallace ships WAVENUMBERS, so 1e8/nu is a VACUUM wavelength,
     while the Kitt Peak `lm` files are AIR. At 3300 A that is a 0.94 A error -- far wider
     than a line -- and it would masquerade as a normalisation difference. The conversion
     is verified by cross-correlation, and the run ABORTS if the residual shift is large.
  2. THE SPECTRA MUST AGREE ON THE LINES. If the two atlases did not correlate, a ratio
     between them would be meaningless. Measured and asserted.

WHAT THIS DOES NOT DO
  * It does not extrapolate RYA-841's lever off the +/-4% grid it was measured on.
  * It does not tune sigma_delta toward any target (RYA-161).

Usage (Sirius -- needs the Kitt Peak atlas and the Wallace atlas):
    python3 scripts/rya846_nearuv_sigma_delta.py
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from config.constants import codex_path  # noqa: E402
from rya759_nearuv_synth import _kp_segments, _load_kp_window  # noqa: E402

LO_A, HI_A = 3000.0, 3780.0

#: RYA-841, measured over 200 fits on a +/-4% grid. Quoted, never re-derived here.
LEVER_DEX_PER_UNIT_DELTA = 2.42
LEVER_P16, LEVER_P84 = 0.9076, 4.6776

#: The value under test. Never used in a computation, only in the comparison.
ASSUMED_TERM_DEX = 0.100

#: Bin width for the normalisation comparison. Wide enough to contain many pixels and
#: several lines, narrow enough that a real continuum drift is not averaged away.
BIN_A = 10.0
BIN_HALF_A = 5.0

#: The cross-correlation control must land within this, or the wavelength convention is
#: wrong and every number downstream is meaningless.
MAX_RESIDUAL_SHIFT_A = 0.05
#: Below this the two atlases are not measuring the same spectrum and a ratio is nonsense.
MIN_SPECTRAL_CORR = 0.90

#: A normalised solar flux outside this is an artifact, not a measurement. Region 5 carries
#: values to +9.98 and -1.00 in this band; they are excluded and COUNTED.
SANE_FLUX = (-0.05, 1.20)

OUT = ROOT / "data" / "results" / "rya846"
SUMMARY = OUT / "rya846_sigma_delta.json"
PER_BIN = OUT / "rya846_normalisation_per_bin.csv"


def wallace_path(region: int) -> Path:
    """The Wallace 2011 atlas lives beside the other solar references (RYA-485)."""
    return Path(codex_path("data.solar_wallace2011_kpno")) / f"sptr.reg{region}"


def vac_to_air(lam_vac: np.ndarray) -> np.ndarray:
    """Edlen/IAU standard. Wallace ships wavenumbers => 1e8/nu is VACUUM; Kitt Peak is AIR.

    Skipping this puts the two atlases 0.94 A apart at 3300 A, which is wider than a line
    and would be read as a normalisation difference rather than a unit error.
    """
    s2 = (1e4 / lam_vac) ** 2
    n = 1.0 + 8.34254e-5 + 2.406147e-2 / (130.0 - s2) + 1.5998e-4 / (38.9 - s2)
    return lam_vac / n


def load_wallace(region: int) -> tuple[np.ndarray, np.ndarray, int]:
    """Air wavelength and flux, with Fortran field overflows made explicit.

    ⚠️ The files contain `*******` where a value overflowed its format width. Those must
    become NaN and be COUNTED -- `genfromtxt(invalid_raise=False)` would drop the rows
    silently, which quietly deletes exactly the pixels where something went wrong.
    """
    d = pd.read_csv(wallace_path(region), sep=r"\s+", header=None,
                    na_values=["*", "**", "***", "****", "*****", "******",
                               "*******", "********"])
    lam = vac_to_air(1e8 / d[0].to_numpy(dtype=float))
    flux = d[1].to_numpy(dtype=float)
    n_overflow = int(np.isnan(flux).sum())
    order = np.argsort(lam)
    return lam[order], flux[order], n_overflow


def verify_alignment(segs, lam_w: np.ndarray, flux_w: np.ndarray,
                     centre: float = 3300.0) -> dict:
    """CONTROL 1+2: the conversion is right and the two atlases see the same spectrum."""
    kw, kf = _load_kp_window(segs, centre, 2.0)
    m = (lam_w >= kw.min() - 2.0) & (lam_w <= kw.max() + 2.0) & np.isfinite(flux_w)
    if m.sum() < 100:
        raise SystemExit(f"Wallace does not cover {centre} A — cannot control the alignment")

    best_shift, best_corr = None, -2.0
    for sh in np.arange(-1.5, 1.501, 0.01):
        fi = np.interp(kw, lam_w[m] + sh, flux_w[m])
        r = float(np.corrcoef(kf, fi)[0, 1])
        if r > best_corr:
            best_shift, best_corr = float(sh), r

    if abs(best_shift) > MAX_RESIDUAL_SHIFT_A:
        raise SystemExit(
            f"WAVELENGTH CONVENTION WRONG: residual shift {best_shift:+.3f} A at "
            f"{centre} A exceeds {MAX_RESIDUAL_SHIFT_A} A. A vacuum/air mix-up is 0.94 A "
            f"here and would be read as a normalisation difference.")
    if best_corr < MIN_SPECTRAL_CORR:
        raise SystemExit(
            f"the two atlases correlate at only {best_corr:.3f} — they are not measuring "
            f"the same spectrum and a ratio between them means nothing")

    print(f"[control] alignment at {centre:.0f} A: residual shift {best_shift:+.3f} A, "
          f"correlation {best_corr:+.4f}  (vacuum-uncorrected would be ~-0.94 A)")
    return {"centre_A": centre, "residual_shift_A": best_shift,
            "spectral_correlation": best_corr}


def compare_normalisation(segs, lam_w, flux_w, *, lo: float, hi: float) -> pd.DataFrame:
    """Per-bin ratio of the two independent normalisations of the same band."""
    rows = []
    for c in np.arange(lo + BIN_HALF_A, hi - BIN_HALF_A + 1e-9, BIN_A):
        try:
            kw, kf = _load_kp_window(segs, float(c), BIN_HALF_A)
        except LookupError:
            continue
        m = (lam_w >= kw.min()) & (lam_w <= kw.max()) & np.isfinite(flux_w)
        if m.sum() < 50:
            continue
        fi = np.interp(kw, lam_w[m], flux_w[m])

        sane = ((fi >= SANE_FLUX[0]) & (fi <= SANE_FLUX[1])
                & (kf >= SANE_FLUX[0]) & (kf <= SANE_FLUX[1]))
        n_insane = int((~sane).sum())
        # Deep cores carry no continuum information and divide badly; the ratio is taken
        # where there is flux to compare.
        ok = sane & (kf > 0.05) & (fi > 0.05)
        if ok.sum() < 50:
            continue

        rows.append(dict(
            centre_A=float(c), n_px=int(ok.sum()), n_insane=n_insane,
            # the pixel-wise ratio: 1.0 if both reductions put the continuum in the same place
            median_ratio=float(np.median(fi[ok] / kf[ok])),
            # the high-percentile ratio traces each reduction's own pseudo-continuum LEVEL,
            # which is the quantity the normalisation actually sets
            p90_ratio=float(np.percentile(fi[ok], 90) / np.percentile(kf[ok], 90)),
            kp_p90=float(np.percentile(kf[ok], 90)),
            wallace_p90=float(np.percentile(fi[ok], 90)),
            # NOT named `corr`: that shadows DataFrame.corr and reads as a method.
            spectral_corr=float(np.corrcoef(kf[ok], fi[ok])[0, 1]),
        ))
    return pd.DataFrame(rows)


def report(d: pd.DataFrame, controls: dict, meta: dict) -> dict:
    print(f"\n=== normalisation comparison, {len(d)} bins of {BIN_A:.0f} A ===")
    print(f"  median per-bin spectral correlation: {d.spectral_corr.median():.4f}")

    out = {}
    for key, label in (("median_ratio", "pixel-wise median ratio"),
                       ("p90_ratio", "90th-percentile (pseudo-continuum) ratio")):
        v = d[key]
        med = float(v.median())
        mad = float(np.median(np.abs(v - med)))
        std = float(v.std(ddof=1))
        out[key] = {"median": med, "MAD": mad, "std": std,
                    "p16": float(np.percentile(v, 16)),
                    "p84": float(np.percentile(v, 84)),
                    "min": float(v.min()), "max": float(v.max())}
        print(f"  {label:<44} median {med:.4f}  MAD {mad:.4f}  std {std:.4f}")

    # sigma_delta is the SPREAD of the disagreement, not its mean: a constant offset
    # between two reductions is a scale both share with the truth or neither does, while
    # the scatter is how much the placement wanders within the band.
    sd_mad = out["p90_ratio"]["MAD"]
    sd_std = out["p90_ratio"]["std"]
    offset = abs(out["p90_ratio"]["median"] - 1.0)

    print(f"\n=== sigma_delta (empirical LOWER BOUND — two reductions, not truth) ===")
    print(f"  systematic offset between reductions : {offset:.4f}  ({offset*100:.2f}%)")
    print(f"  scatter across the band, MAD         : {sd_mad:.4f}  ({sd_mad*100:.2f}%)")
    print(f"  scatter across the band, std         : {sd_std:.4f}  ({sd_std*100:.2f}%)")

    print(f"\n=== the honest continuum term = lever x sigma_delta ===")
    print(f"  lever (RYA-841, measured on a +/-4% grid) = {LEVER_DEX_PER_UNIT_DELTA:.2f}")
    terms = {}
    for nm, sd in (("offset", offset), ("MAD", sd_mad), ("std", sd_std)):
        t = LEVER_DEX_PER_UNIT_DELTA * sd
        lo = LEVER_P16 * sd
        hi = LEVER_P84 * sd
        terms[nm] = {"sigma_delta": sd, "term_dex": t,
                     "term_p16": lo, "term_p84": hi}
        print(f"  sigma_delta = {sd:.4f} ({nm:6s}) -> term {t:.3f} dex "
              f"[{lo:.3f}-{hi:.3f} across the lever's 16-84 pct]")
    print(f"\n  assumed by the budget today: {ASSUMED_TERM_DEX:.3f} dex")

    if max(sd_mad, sd_std, offset) > 0.04:
        print("  ⚠️ sigma_delta exceeds the ~4.1% the 0.100 dex term implicitly asserts.")

    summary = {
        "ticket": "RYA-846",
        "measurement": "sigma_delta, near-UV continuum-placement uncertainty",
        "band_A": [LO_A, HI_A],
        "method": ("difference between two INDEPENDENT reductions of the same band — "
                   "Kitt Peak (Kurucz/Brault) vs Wallace+2011 (NSO). This is an empirical "
                   "LOWER BOUND on sigma_delta: it measures reduction-to-reduction "
                   "disagreement, not distance from an unobservable true continuum."),
        "dead_route": ("the Kitt Peak segment overlaps (20 seams, ~1.0 A each, in band) are "
                       "BYTE-IDENTICAL duplicates — same grid, flux difference exactly 0.0 "
                       "— so they carry no independent information and would have returned "
                       "sigma_delta ~ 0 and a continuum term of ~0."),
        "controls": controls,
        "ratios": out,
        "sigma_delta": {"offset": offset, "MAD": sd_mad, "std": sd_std},
        "lever": {"value": LEVER_DEX_PER_UNIT_DELTA, "p16": LEVER_P16, "p84": LEVER_P84,
                  "source": "RYA-841, 200 fits, measured on a +/-4% grid ONLY"},
        "terms_dex": terms,
        "assumed_term_dex": ASSUMED_TERM_DEX,
        **meta,
    }
    return summary


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--region", type=int, default=6,
                    help="Wallace region: 6 is clean over 2899-3525 A (default); "
                         "5 spans 2995-4106 A but carries overflow markers")
    ap.add_argument("--lo", type=float, default=None)
    ap.add_argument("--hi", type=float, default=None)
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    lam_w, flux_w, n_overflow = load_wallace(a.region)
    lo = a.lo if a.lo is not None else max(LO_A, float(np.nanmin(lam_w)))
    hi = a.hi if a.hi is not None else min(HI_A, float(np.nanmax(lam_w)))
    print(f"[wallace] region {a.region}: {lam_w.min():.1f}-{lam_w.max():.1f} A (air), "
          f"{len(lam_w)} px, {n_overflow} overflow markers -> NaN")
    print(f"[band]    comparing {lo:.1f}-{hi:.1f} A "
          f"({hi-lo:.0f} A of the {HI_A-LO_A:.0f} A near-UV band)")

    segs = _kp_segments()
    print(f"[kp]      {len(segs)} atlas segments")
    controls = verify_alignment(segs, lam_w, flux_w)

    d = compare_normalisation(segs, lam_w, flux_w, lo=lo, hi=hi)
    if d.empty:
        raise SystemExit("no comparable bins — check coverage")
    d.to_csv(PER_BIN, index=False)

    summary = report(d, controls, {
        "wallace_region": a.region,
        "wallace_overflow_markers": n_overflow,
        "compared_A": [lo, hi],
        "band_fraction_compared": (hi - lo) / (HI_A - LO_A),
        "n_bins": int(len(d)),
        "bin_width_A": BIN_A,
    })
    SUMMARY.write_text(json.dumps(summary, indent=2))
    print(f"\n[out] {PER_BIN}\n[out] {SUMMARY}")


if __name__ == "__main__":
    main()
