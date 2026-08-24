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
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline.lines_fit import (  # noqa: E402
    _local_renorm, _fit_profile, _integrate_profile, _ew_error)
from config.constants import PIPELINE as _PIPE, STAR_PARAMS  # noqa: E402  (RYA-911: quote the params)
from pipeline.band_products import carried_ep, LineMeasurement, assert_single_element  # noqa: E402
from pipeline.band_policy import check_intake, resolve as resolve_band  # noqa: E402
from scripts.measure_band_ew import (  # noqa: E402
    # RYA-913: load_kp_window removed -- it was imported and never called, a
    # Kitt-Peak-specific flux reader sitting in the import line of a route that
    # correctly dispatches. An unused import of the thing the rule forbids is how
    # the next instance starts.
    # RYA-911 rebase: UNION, not swap -- 913 drops load_kp_window, 911 adds
    # reference_continuum. Both changes are kept.
    kp_segments, load_window, load_window_ex, telluric_reason,
    reference_continuum, verify_feature,
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

#: 🔴 RYA-911 — how far above its own continuum the observed flux may sit before the
#: continuum is disqualified. NOT a tuned threshold: 1.0 is the physical bound and this
#: is the bound plus a noise allowance. HARPS's local re-fit reached 1.204 and the
#: product's own column never exceeds 1.007, so the decision is nowhere near this number
#: -- which is the property a threshold needs when no sweep can be run on it (RYA-847).
CONTINUUM_MAX_OVER = 1.02


# Stellar broadening scale for solar-type Fe lines: thermal + micro + macroturbulence,
# expressed as a velocity sigma. ~1.7 km/s reproduces the observed solar Fe I widths
# (sigma ≈ 0.031 Å at 5500 Å). It is a starting GUESS for the fit, not an imposed width —
# the fit is free to move away from it between the bounds below.
STELLAR_SIGMA_KMS = 1.7
C_KMS = 299792.458


#: RYA-906/911 — the width physics moved to `pipeline/line_width.py` so the OTHER
#: profile fitter (`pipeline/measure/profile_fit.py`) could import it. It could not:
#: this module's import chain loads the Kitt Peak atlas, so the corrected guard was
#: unreachable from the pipeline and that site kept the refuted sigma-only test.
#: Re-exported under the old private names so this file's call sites read unchanged.
from pipeline.line_width import (  # noqa: E402
    voigt_fwhm,
    irreducible_sigma_kms as _irreducible_sigma_kms,
    physical_floor_fwhm as _physical_floor_fwhm,
    total_width_below_physical_floor as _total_width_below_physical_floor,
    under_physical_width_reason as _under_physical_width_reason,
    # RYA-959 — the ceiling. The floor above had no mirror, so a fit that failed escaped
    # UPWARD into `sigma_max` and returned a 400-1000 mA "Fe I line" that looked converged.
    gaussian_sigma_above_physical_ceiling as _sigma_above_physical_ceiling,
    over_physical_width_reason as _over_physical_width_reason,
    implied_width_A as _implied_width_A,
    implied_width_exceeds_ceiling as _implied_width_exceeds_ceiling,
    over_implied_width_reason as _over_implied_width_reason,
)


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
    ap.add_argument("--holding", default=None, help="Name ONE holding explicitly instead of taking the instrument's first covering candidate (RYA-933/934). An instrument that serves several products otherwise answers with whichever is listed first, which is how two telluric-corrected Kitt Peak holdings sat unmeasurable while `--instrument kpno_solar_atlas` looked like it covered them.")
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

    segs = kp_segments() if a.instrument == "kpno_solar_atlas" else None
    allw = acc.wave_air_A.values
    # 🔴 RYA-904 — `pre` USED TO BE READ ONCE, FROM AN INSTRUMENT-KEYED DICT. It is a
    # property of the HOLDING, and one instrument can serve holdings that disagree:
    # CRIRES+'s raw Vesta IDPs are un-normalised adu while its RYA-794 corrected Y arm
    # arrives normalised. A single per-run answer is therefore wrong for one of them, and
    # being wrong in the "already normalised, renormalise anyway" direction is the
    # RYA-713 double-normalisation defect (EWs low by a median 11.7 %, up to 71.4 %). It
    # now comes off the window that was actually loaded, per line.

    rows, skipped, causes = [], [], []
    for _, r in sel.iterrows():
        c = float(r.wave_air_A)
        why = telluric_reason(c, a.instrument)
        if why:
            skipped.append(dict(wave=c, reason=why)); continue

        # Per-line isolation: nearest catalogued neighbour, either side.
        _d = np.abs(allw - c); _d = _d[_d > 1e-4]
        gap = float(_d.min()) if len(_d) else 99.0
        allowed, gap_reason = permits_profile_fit_for_line(c, gap)
        if not allowed:
            lm = LineMeasurement(element=a.element, ion=a.ion, wavelength_air_A=c,
                                 instrument=a.instrument, ew_mA=float("nan"),
                                 ep_eV=carried_ep(r, wavelength_A=c, element=a.element,
                                                  ion=a.ion),
                                 ew_method=f"profile fit not attempted in {pol.name}")
            lm.in_aggregate = False
            lm.excluded_reason = f"BAND-POLICY: {gap_reason}"
            rows.append(lm)
            continue
        try:
            # RYA-783: dispatched by instrument so the IAG arm runs through the
            # same harness rather than a parallel copy of it.
            win = load_window_ex(a.instrument, c, CONT_HALF_A * 1.3, segs,
                                 holding=a.holding)
            w, f, src, pre = win.wave, win.flux, win.provenance, win.pre_normalised
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
                cont_method = (f"pre-normalised: {win.holding.holding_id} ships its own "
                               f"continuum and unity IS it, by construction — no local "
                               f"fit is made")
            else:
                fn, cont = _local_renorm(wf, ff, c, window=CONT_HALF_A)
                cont_method = (
                    f"local-linear-refit: polyfit deg 1 through the top "
                    f"{100 - _PIPE['cont_anchor_percentile']:.0f}% of pixels in the outer "
                    f"{_PIPE['cont_edge_frac']:.0%} of each +/-{CONT_HALF_A} A edge strip")
            # 🔴 RYA-911 — RECORD IT. `_local_renorm` has always returned the continuum
            # it placed and this driver has always bound it and thrown it away, which is
            # why the HARPS Fe II -0.34 dex deficit had to be diagnosed by REFITTING an
            # approximation of this code rather than by reading it. Taken AT THE LINE
            # CENTRE, because that is the divisor that sets this line's depth.
            _j = int(np.argmin(np.abs(wf - c)))
            cont_level = float(cont[_j])
            cont_ref = reference_continuum(win.holding, c)
            # 🔴 RYA-911 — THE ONE-SIDED PHYSICAL TEST, APPLIED TO EVERY ARM.
            #
            # A continuum-normalised stellar spectrum cannot exceed 1.0 by any real
            # margin: flux above the continuum is not a small error, it is impossible.
            # So `max(F/C) > 1` CONVICTS a continuum outright, with no reference, no
            # reconstruction and no second opinion needed.
            #
            # This is what was missing. The HARPS arm re-fit its continuum through 0.3 A
            # anchor strips of the crowded solar optical, landed BELOW the flux on half
            # the pool, and shrank every EW on it -- and nothing objected, because a
            # coherent deficit looks exactly like a smaller line. RYA-897 could only see
            # it as a 0.34 dex hole at the far end of the pipeline.
            #
            # ⚠️ ONE-SIDED ON PURPOSE. The converse is NOT a defect: in the blanketed
            # blue there may be no unabsorbed pixel in the window at all, so a correct
            # continuum SHOULD sit above every observed point. Testing that direction too
            # would quarantine correct continua in exactly the bands that need them most.
            _over = float(np.max(fn)) if fn.size else float("nan")
            cont_unphysical = (np.isfinite(_over) and _over > CONTINUUM_MAX_OVER)

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
            # RYA-959 — THE OBSERVED CORE, taken where THIS fit put its centre.
            # The implied-width test pairs the integral with the core it came from, so
            # the depth has to be the one the same fit was describing; reading it at the
            # catalogued position instead would convict every line the fit legitimately
            # re-centred within its +/-0.30 A freedom.
            _core = np.abs(wf - float(popt[0])) <= max(3 * s_init, 0.03)
            obs_depth = float(1.0 - np.min(fn[_core])) if _core.any() else float("nan")
        except Exception as e:
            skipped.append(dict(wave=c, reason=f"{type(e).__name__}: {e}"))
            continue

        lm = LineMeasurement(
            element=a.element, ion=a.ion, wavelength_air_A=c, instrument=a.instrument,
            ew_mA=float(ew),
            # RYA-871 — carried from the SAME accounting row `c` came from, so the
            # wavelength and the EP describe one transition by construction.
            ep_eV=carried_ep(r, wavelength_A=c, element=a.element, ion=a.ion),
            ew_method=(f"PROFILE-FIT ({ptype}); gap {gap:.3f} A"
                       + ("" if band_permits else f" [PER-LINE ESCAPE from {pol.name} ban]")
                       + f"; model integrated; chi2_red={chi2:.4g}; "
                       f"sigma_fit={popt[2]:.4f} A (init {s_init:.4f}, floor {s_min:.4f}); "
                       f"ew_err={err:.2f} mA; continuum linear anchors at +/-{CONT_HALF_A} A; "
                       f"source {src}"),   # RYA-904: names the HOLDING
            # RYA-911 — the continuum is a first-class column, not prose. A number a
            # downstream RCA has to parse out of a sentence is a number it will
            # reconstruct instead.
            continuum_level=cont_level, continuum_method=cont_method,
            continuum_ref=cont_ref,
            # RYA-911 — the fit's own numbers, as COLUMNS. `chi2` and `popt[2]` were
            # computed here, used to decide the width verdict, and then written only into
            # the `ew_method` sentence. `sigma_A` is deliberately NOT used: it means one
            # sigma on A in DEX, and this is a width in ANGSTROM.
            profile_sigma_A=float(popt[2]), profile_sigma_floor_A=float(s_min),
            # RYA-906 — the Lorentzian half-width, WITHOUT which `profile_sigma_A` cannot
            # be interpreted: the two are degenerate in a Voigt fit.
            profile_gamma_A=(float(popt[3]) if ptype == "voigt" and len(popt) > 3 else None),
            red_chi2=float(chi2),
            # RYA-959 — recorded on EVERY row, not only the rejected ones. The quantity
            # that convicts a bad EW is the quantity an auditor needs on the good ones too.
            observed_depth=(obs_depth if np.isfinite(obs_depth) else None),
            implied_width_A=(_implied_width_A(float(ew), obs_depth)
                             if np.isfinite(obs_depth) else None))
        if cont_unphysical:
            # Quarantined, never culled (RYA-711): measured, kept, reported with its
            # number. The EW is still written down -- it is the CONTINUUM that is
            # disqualified, and a reader needs to see how far off it was.
            lm.in_aggregate = False
            lm.excluded_reason = (
                f"CONTINUUM-UNPHYSICAL: {_over:.4f} of the continuum is exceeded by the "
                f"observed flux inside the fit window (limit {CONTINUUM_MAX_OVER:.3f}). "
                f"Flux above the continuum is impossible, so this continuum sits BELOW "
                f"the spectrum and every EW measured on it is too small. Method was: "
                f"{cont_method}. This is the RYA-911 defect that cost HARPS Fe II a "
                f"median 23.8% of its EW.")
        elif _total_width_below_physical_floor(popt, ptype, c, s_min):
            # 🔴 RYA-906 — THIS TEST WAS WRONG AND IS CORRECTED HERE.
            #
            # It used to read `abs(popt[2] - s_min) < 1e-4` and reject on the Gaussian
            # sigma alone, arguing "no observed feature can be narrower than the
            # instrument profile". Measured: the guard fired on 25 HARPS Fe II lines and
            # EVERY ONE of them was a VOIGT fit — it never once fired on a Gaussian.
            #
            # In a Voigt profile sigma and gamma are DEGENERATE. A railed sigma means the
            # optimiser put the width into the LORENTZIAN, which is a statement about
            # where a degeneracy landed, not about how wide the line is. The total width
            # can be entirely physical while sigma sits on its bound, and the EW — which
            # is the integral of the WHOLE profile — is then fine.
            #
            # So the test now asks the question it always meant to ask: is the TOTAL
            # width below what physics permits? The floor is instrumental (+) thermal (+)
            # microturbulent, from the star's ratified parameters, with macroturbulence
            # deliberately excluded so it stays a true lower bound.
            lm.in_aggregate = False
            lm.excluded_reason = _under_physical_width_reason(popt, ptype, c, s_min)
        elif _sigma_above_physical_ceiling(popt, c, s_min):
            # 🔴 RYA-959 — THE MIRROR THAT WAS MISSING, and the fault that built
            # RYA-958's contaminated pool. `instrument_sigma` allows sigma out to 0.40 A,
            # ~5x the widest Gaussian solar Doppler physics permits, so a fit with nothing
            # good to converge on ran to the bound and integrated a pedestal instead of a
            # line. Quarantined, never culled (RYA-711): the EW stays on the row so a
            # reader can see how far off it was, and `attribute_root_cause` below still
            # records WHY the fit had nothing to hold on to.
            lm.in_aggregate = False
            lm.excluded_reason = _over_physical_width_reason(popt, ptype, c, s_min, s_max)
        elif _implied_width_exceeds_ceiling(float(ew), obs_depth, 2 * FIT_HALF_A):
            # RYA-959 — the INTEGRAL bound. A runaway Lorentzian gamma inflates the EW
            # without deepening the core, so it clears the sigma ceiling above and is
            # caught only here. This is the check RYA-958's 251-of-444 diagnosis ran by
            # hand, now wired into the measurement that produces the number.
            lm.in_aggregate = False
            lm.excluded_reason = _over_implied_width_reason(
                float(ew), obs_depth, c, s_min, 2 * FIT_HALF_A, float(r.predicted_depth))
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
    # RYA-933/934 -- the HOLDING goes in the stem whenever one was named. The stem
    # carried only the INSTRUMENT, so `solar_harps` and `solar_harps_molecfit_corrected`
    # both wrote FeII_3800_6910_harps_PROFILEFIT_ew.csv and the second silently
    # overwrote the first. Two products of the same instrument differing precisely by
    # whether tellurics were removed is the one pair this repo must never collapse
    # (RYA-911, RYA-927 identity canary).
    stem = f"{a.element}{a.ion}_{int(a.lo)}_{int(a.hi)}_{a.instrument}"
    if a.holding:
        stem += f"_{a.holding}"
    stem += "_PROFILEFIT"
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
