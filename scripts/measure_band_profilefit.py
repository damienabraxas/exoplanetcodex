#!/usr/bin/env python3
"""Measure EWs by PROFILE FIT for one element/band/instrument — RYA-713 / RYA-429.

    python3 scripts/measure_band_profilefit.py --element Fe --ion I --lo 3924 --hi 6905

Replaces the interval-integration harness, which the optical control falsified: against
the HARPS pool over 146 lines it returned a median EW ratio of 0.773 (−0.112 dex) with a
5× spread, because a separation-derived window clips the wings of crowded strong lines and
over-reaches on isolated weak ones. You cannot both exclude a neighbour and include the
wings by choosing an interval — which is why the reference method fits a profile and
integrates the MODEL.

REUSES the RYA-429 fitter rather than reimplementing it. `pipeline/lines_fit` owns the
continuum renormalisation, the Voigt/Gaussian fit and the model integration; this driver
supplies the spectrum, the line list, and the instrument's resolution. Copying those
functions here would have produced a second EW definition that drifts from the first —
the RYA-701 failure mode.

INSTRUMENT RESOLUTION IS AN ARGUMENT, NOT A CONSTANT
----------------------------------------------------
The instrumental profile width belongs to the instrument. HARPS at R ≈ 115 000 has
σ ≈ 0.020 Å; the Kitt Peak FTS at R ≈ 500 000 has σ ≈ 0.005 Å — which is exactly the old
hardcoded lower bound, so a naive reuse would pin the fit against its own constraint and
report a width it never measured. Resolution is resolved from
`data/catalog/instrument_catalog.csv`, never typed in here.
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

from pipeline.lines_fit import (  # noqa: E402
    _local_renorm, _fit_profile, _integrate_profile, _ew_error)
from pipeline.band_products import LineMeasurement, assert_single_element  # noqa: E402
from pipeline.band_policy import check_intake, resolve as resolve_band  # noqa: E402
from scripts.measure_band_ew import (  # noqa: E402
    kp_segments, load_kp_window, telluric_reason, PRE_NORMALISED, verify_feature,
    attribute_root_cause)

ACCOUNTING = ROOT / "data" / "audit" / "line_accounting" / "per_line.csv"
CATALOG = ROOT / "data" / "catalog" / "instrument_catalog.csv"
OUT = ROOT / "data" / "measured" / "band_ew"

FWHM_TO_SIGMA = 1.0 / 2.35482
# Fit window. Wide enough that the model's wings are integrated rather than truncated,
# and the continuum anchors sit outside the line. This is NOT the blend-avoidance window
# that broke the interval method -- a fit models the neighbour, it does not exclude it.
FIT_HALF_A = 0.60
CONT_HALF_A = 1.20


# Stellar broadening scale for solar-type Fe lines: thermal + micro + macroturbulence,
# expressed as a velocity sigma. ~1.7 km/s reproduces the observed solar Fe I widths
# (sigma ≈ 0.031 Å at 5500 Å). It is a starting GUESS for the fit, not an imposed width —
# the fit is free to move away from it between the bounds below.
STELLAR_SIGMA_KMS = 1.7
C_KMS = 299792.458


def instrument_sigma(instrument: str, wavelength_A: float) -> tuple[float, float, float]:
    """(sigma_init, sigma_min, sigma_max) in Å for the fit at this wavelength.

    THE LINE WIDTH IS NOT THE INSTRUMENT'S (RYA-713). An observed line is broadened by the
    star AND the instrument, added in quadrature, and at high resolving power the star
    dominates completely:

        at 5500 Å   HARPS sigma_inst 0.0203   Kitt Peak sigma_inst 0.0047
                    stellar sigma    0.0312   -> observed on Kitt Peak 0.0315

    Seeding the fit from the instrumental width alone put the start four times below the
    truth on Kitt Peak; the optimiser walked to the lower bound and stayed there, and
    Fe I 4995 came back 20.6 mÅ against a pool value of 138. So:

      * `sigma_init` = quadrature sum of instrumental and stellar — the physical expectation.
      * `sigma_min`  = the INSTRUMENTAL sigma, which is a hard floor: no observed feature
                       can be narrower than the instrument's own profile. A fit that lands
                       on this bound is reporting a width it never measured, and the driver
                       flags it rather than trusting the EW.
      * `sigma_max`  = generous; rotation and unresolved blends legitimately broaden lines.
    """
    cat = pd.read_csv(CATALOG)
    idcol = cat.columns[0]
    row = cat[cat[idcol].astype(str) == instrument]
    if not len(row):
        raise SystemExit(f"instrument {instrument!r} not in {CATALOG.name}; add it there "
                         f"rather than hardcoding a width here")
    R = float(row.iloc[0]["resolving_power_max"])
    sigma_inst = (wavelength_A / R) * FWHM_TO_SIGMA
    sigma_star = wavelength_A * STELLAR_SIGMA_KMS / C_KMS
    return float(np.hypot(sigma_inst, sigma_star)), float(sigma_inst), 0.40


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--element", required=True)
    ap.add_argument("--ion", default="I")
    ap.add_argument("--lo", type=float, required=True)
    ap.add_argument("--hi", type=float, required=True)
    ap.add_argument("--instrument", default="kpno_solar_atlas")
    ap.add_argument("--max-lines", type=int, default=None)
    ap.add_argument("--depth-min", type=float, default=0.05)
    ap.add_argument("--depth-max", type=float, default=0.60)
    ap.add_argument("--out", default=str(OUT))
    a = ap.parse_args()

    # Band policy decides whether a profile fit is valid here (RYA-713). Where the band
    # bans it on its MEDIAN gap, individual lines may still escape -- isolation is a
    # per-line property, and ~306 of 1612 near-UV Fe lines are isolated enough to fit.
    from pipeline.band_policy import permits_profile_fit_for_line, BandPolicyError as _BPE
    pol = resolve_band((a.lo + a.hi) / 2.0)
    band_permits = "profile-fit" in pol.permitted_methods
    if not band_permits:
        print(f"  {pol.name} bans profile fitting by default (median gap "
              f"{pol.median_gap_A:.3f} A). Running PER-LINE: only lines whose nearest "
              f"neighbour is >= {__import__('pipeline.band_policy', fromlist=['x']).PROFILE_FIT_MIN_GAP_A} A "
              f"are measured; the rest are quarantined with that reason.")
    print(f"{a.element} {a.ion}  {a.lo:.0f}-{a.hi:.0f} A  {a.instrument}  band={pol.name}")
    if pol.telluric_required:
        print("  NOTE: this band requires telluric correction; telluric windows are skipped, "
              "not corrected — lines outside them are unaffected.")

    acc = pd.read_csv(ACCOUNTING)
    sel = acc[(acc.element == a.element) & (acc.ion == a.ion) &
              (acc.wave_air_A >= a.lo) & (acc.wave_air_A <= a.hi) &
              acc.predicted_depth.between(a.depth_min, a.depth_max) &
              acc.instruments.notna()].sort_values("wave_air_A").reset_index(drop=True)
    if a.max_lines:
        idx = np.unique(np.linspace(0, len(sel) - 1, a.max_lines).astype(int))
        sel = sel.iloc[idx].reset_index(drop=True)
    print(f"  candidates {len(sel)}")

    segs = kp_segments()
    allw = acc.wave_air_A.values
    pre = PRE_NORMALISED[a.instrument]

    rows, skipped, causes = [], [], []
    for _, r in sel.iterrows():
        c = float(r.wave_air_A)
        why = telluric_reason(c)
        if why:
            skipped.append(dict(wave=c, reason=why)); continue

        # Per-line isolation: nearest catalogued neighbour, either side.
        _d = np.abs(allw - c); _d = _d[_d > 1e-4]
        gap = float(_d.min()) if len(_d) else 99.0
        allowed, gap_reason = permits_profile_fit_for_line(c, gap)
        if not allowed:
            lm = LineMeasurement(element=a.element, ion=a.ion, wavelength_air_A=c,
                                 instrument=a.instrument, ew_mA=float("nan"),
                                 ew_method=f"profile fit not attempted in {pol.name}")
            lm.in_aggregate = False
            lm.excluded_reason = f"BAND-POLICY: {gap_reason}"
            rows.append(lm)
            continue
        try:
            w, f, src = load_kp_window(segs, c, pad=CONT_HALF_A * 1.3)
            m = np.abs(w - c) <= CONT_HALF_A
            wf, ff = w[m], f[m]
            if wf.size < 40:
                raise ValueError(f"only {wf.size} points in the fit window")

            # CONTINUUM. The same policy the interval method needed, which I failed to
            # apply here on the first pass -- and it cost -0.088 dex.
            #
            # Kitt Peak ships RESIDUAL FLUX: unity IS the continuum, by construction.
            # Re-fitting a local continuum through percentile-filtered edge strips lands
            # BELOW unity in crowded spectrum (those strips contain lines), so dividing by
            # it makes every line shallower and every EW smaller -- a coherent deficit, not
            # scatter. That is exactly the residual the control still showed after the
            # profile fit fixed the scatter.
            #
            # So on pre-normalised data the atlas continuum is used as-is. On
            # un-normalised data (HARPS) the local renorm is required and is applied.
            if pre:
                fn, cont = ff, np.ones_like(ff)
            else:
                fn, cont = _local_renorm(wf, ff, c, window=CONT_HALF_A)

            fit_m = np.abs(wf - c) <= FIT_HALF_A
            s_init, s_min, s_max = instrument_sigma(a.instrument, c)
            popt, pcov, ptype, chi2 = _fit_profile(
                wf[fit_m], fn[fit_m], c,
                sigma_init=s_init, sigma_min=s_min, sigma_max=s_max,
                core_half_A=max(3 * s_init, 0.03))
            if popt is None or ptype == "failed":
                raise ValueError("profile fit did not converge")
            ew = _integrate_profile(wf[fit_m], popt, ptype)
            err = _ew_error(fn[~fit_m], ew, popt, pcov, ptype)
            ok, vwhy = verify_feature(wf, fn, c, FIT_HALF_A, float(r.predicted_depth))
        except Exception as e:
            skipped.append(dict(wave=c, reason=f"{type(e).__name__}: {e}"))
            continue

        lm = LineMeasurement(
            element=a.element, ion=a.ion, wavelength_air_A=c, instrument=a.instrument,
            ew_mA=float(ew),
            ew_method=(f"PROFILE-FIT ({ptype}); gap {gap:.3f} A"
                       + ("" if band_permits else f" [PER-LINE ESCAPE from {pol.name} ban]")
                       + f"; model integrated; chi2_red={chi2:.4g}; "
                       f"sigma_fit={popt[2]:.4f} A (init {s_init:.4f}, floor {s_min:.4f}); "
                       f"ew_err={err:.2f} mA; continuum linear anchors at +/-{CONT_HALF_A} A; "
                       f"segment(s) {src}"))
        if abs(float(popt[2]) - s_min) < 1e-4:
            lm.in_aggregate = False
            lm.excluded_reason = (
                f"FIT-PINNED: fitted sigma {popt[2]:.4f} A sits on the instrumental floor "
                f"{s_min:.4f} A. No observed feature can be narrower than the instrument "
                f"profile, so the optimiser hit a bound rather than a minimum and the "
                f"integrated EW is not trustworthy.")
        elif not ok:
            lm.in_aggregate = False
            lm.excluded_reason = f"FEATURE-VERIFICATION: {vwhy}"
        if not lm.in_aggregate:
            causes.append(dict(wave=c, symptom=lm.excluded_reason[:90],
                               **attribute_root_cause(wf, fn, c, FIT_HALF_A,
                                                      float(r.predicted_depth or 0.0),
                                                      lm.excluded_reason, allw)))
        rows.append(lm)

    assert_single_element(rows, a.element)
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    stem = f"{a.element}{a.ion}_{int(a.lo)}_{int(a.hi)}_{a.instrument}_PROFILEFIT"
    if a.max_lines:
        # A SUBSET run never writes over a full one. A --max-lines smoke test
        # silently clobbered a documented 445-line result; the fixed output path
        # made a partial run indistinguishable from a complete one on disk.
        stem += f"_SUBSET{a.max_lines}"
    pd.DataFrame([vars(l) for l in rows]).to_csv(out / f"{stem}_ew.csv", index=False)
    (out / f"{stem}_skipped.json").write_text(json.dumps(skipped, indent=2))
    pd.DataFrame(causes).to_csv(out / f"{stem}_root_causes.csv", index=False)

    df = pd.DataFrame([vars(l) for l in rows])
    print(f"\n  measured {len(rows)}, skipped {len(skipped)}")
    if len(df):
        print(f"  EW range {df.ew_mA.min():.1f} - {df.ew_mA.max():.1f} mA")
        print(f"  in aggregate {int(df.in_aggregate.sum())}, quarantined "
              f"{int((~df.in_aggregate).sum())}")
    print(f"  wrote {out / (stem + '_ew.csv')}")


if __name__ == "__main__":
    main()
