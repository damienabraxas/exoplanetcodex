"""The standing per-band continuum rule -- a local polynomial on the upper envelope of
each fit window. RYA-1230, ratified by Ryan 2026-09-26.

WHY IT EXISTS
-------------
Ryan, 2026-09-26: *"There should be a standing rule, even though it has little effect on
the abundance measurement, where we do a continuum for each band."* Chosen shape: a LOCAL
polynomial per fit window, on every holding, including the pre-normalised atlases.

The measurement that forced it (RYA-1230): the "pre-normalised" IAG and Kitt Peak atlases
sit 0.5-1.0% BELOW unity around the solar N I lines (p95 0.990-0.995 within +/-1.5 A,
stable out to +/-5 A). Those lines are only 0.3-1% deep, so the synthesis route -- which
compares observed flux against a synthesis normalised to exactly 1 -- spent the offset
on A(N): median -0.142 dex over 12 line x holding cells, up to -0.44 on one line, at
0.02-0.04 dex per 0.1% of continuum. For a 100 mA Fe line the same 0.5% is a few
thousandths of a dex, which is why nobody saw it on iron.

WHAT THIS IS NOT (the RYA-1026 boundary)
----------------------------------------
RYA-1026 forbids re-normalising a pre-normalised product, after two SILENT corruptions:
a free BAND-WIDE polynomial that followed a saturated telluric band down (RYA-940), and a
wrong `pre_normalised` flag that re-normalised KP2005 for months (RYA-929). This rule is
bounded against both:

  * LOCAL, never band-wide: a straight line over +/-ENV_HALF_WIDTH_A of one window. It
    cannot tilt a band, because it never sees one.
  * UPPER ENVELOPE, not a free fit: the p95 of each 1 A bin (the RYA-1000 /
    verify_feature estimator, reused rather than invented), with the bins that sit below
    the envelope -- the ones a line or telluric band owns -- clipped before the fit.
  * BOUNDED OUTPUT: a local level more than MAX_SHIFT from unity is NOT applied. A 3%
    "continuum" inside +/-5 A of a solar line is a telluric band or a line forest, not a
    normalisation error, and dividing by it is the RYA-940 failure. The line keeps the
    shipped continuum and the record says why (CONTINUUM_UNPLACEABLE).
  * RECORDED ON EVERY LINE: level at the line centre, slope, bins used, applied or not.
    Nothing is divided silently.

`prenormalised_guard.assert_not_renormalising(..., local_window_envelope=True)` is the one
operation the guard now admits, and only this module passes it.

NEAR-UV IS THE STATED EXCEPTION. There the true continuum is never observed (median flux
0.283-0.805, RYA-759/1189), so an upper envelope is a PSEUDO-continuum and dividing by it
would remove real opacity (RYA-1189: "a pseudo-continuum is not a continuum"). The band
config says `apply: false` there; the envelope is still MEASURED and recorded.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

#: Half-width of the envelope window, A. MEASURED on the 12 N I cells (RYA-1230): the p95
#: level at +/-1.5 A and +/-5 A agrees to <= 0.002 on 11 of 12; at +/-15 A it drifts (IAG
#: 7442 0.990 -> 0.993) and on one holding diverges outright (KP-mf 8216: 1.001 at 1.5 A,
#: 1.045 at 15 A -- a broken composite seam). Wide enough to find continuum bins in a
#: 0.3-1%-deep line's neighbourhood, narrow enough to stay local.
ENV_HALF_WIDTH_A = 5.0
#: RYA-1230 RYA-587 `continuum` leg: the SAME envelope, same bins, same window, with only
#: the per-bin estimator varied (p90 / p99 around the nominal p95), so the envelope stays
#: constrained and the leg prices PLACEMENT, not "rule off". Env-scoped and loud.
#: (An earlier +/-1.5 A leg left < 5 bins, went CONTINUUM_UNCONSTRAINED and measured the
#: whole correction instead -- withdrawn.)
import os as _os
ENV_PERCENTILE = float(_os.environ.get("CODEX_CONT_PCT") or 95.0)
if ENV_PERCENTILE != 95.0:
    print(f"  \u26a0\ufe0f  CONTINUUM ENVELOPE PERCENTILE OVERRIDE (RYA-1230 budget leg): p{ENV_PERCENTILE:g}")
#: MODEL-GUIDED selection quantile: the pixels where the SYNTHESIS is in its top decile
#: inside the window. The budget legs vary it (CODEX_CONT_Q = 80 / 97).
MODEL_Q = float(_os.environ.get("CODEX_CONT_Q") or 90.0)
if MODEL_Q != 90.0:
    print(f"  \u26a0\ufe0f  CONTINUUM MODEL-QUANTILE OVERRIDE (RYA-1230 budget leg): q{MODEL_Q:g}")
#: Minimum selected pixels for a model-guided fit.
MIN_MODEL_PIX = 20
#: Model-guided anchors whose obs/reference residual exceeds this many MADs are clipped
#: (iterated); see fit_model_guided.
CLIP_MAD_MODEL = 3.0
#: Envelope bin width, A. Each bin contributes its p95.
BIN_A = 1.0
#: Polynomial degree. Linear: over 10 A the only shape a normalisation error can have that a
#: line cannot is a level and a slope. A quadratic starts to fit line wings.
DEGREE = 1
#: Bins more than this many MADs BELOW the first-pass envelope are line/telluric bins.
CLIP_MAD = 2.0
#: Minimum surviving bins for a fit. Below this the envelope is not constrained.
MIN_BINS = 5
#: Refuse to apply a local level further than this from unity (see module docstring).
MAX_SHIFT = 0.03


@dataclass
class LocalContinuum:
    level_at_centre: float
    slope_per_A: float
    n_bins: int
    applied: bool
    reason: str
    method: str = (f"RYA-1230 local envelope: p{ENV_PERCENTILE:g} per {BIN_A:g} A bin over +/-"
                   f"{ENV_HALF_WIDTH_A:g} A, {CLIP_MAD:g}-MAD low clip, degree {DEGREE}")


def fit(wave_A, flux, centre_A: float, *, env_half_width_A: float = ENV_HALF_WIDTH_A,
        apply: bool = True) -> tuple[LocalContinuum, np.ndarray | None]:
    """Fit the local upper envelope around `centre_A`.

    Returns (record, continuum array on `wave_A`) -- the array is None when the level is
    not applied, so a caller cannot divide by a refused continuum by accident.
    """
    w = np.asarray(wave_A, float)
    f = np.asarray(flux, float)
    ok = np.isfinite(w) & np.isfinite(f)
    m = ok & (np.abs(w - centre_A) <= env_half_width_A)
    edges = np.arange(centre_A - env_half_width_A, centre_A + env_half_width_A + 1e-9, BIN_A)
    bx, by = [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        s = m & (w >= lo) & (w < hi)
        if s.sum() >= 5:
            bx.append(0.5 * (lo + hi)); by.append(float(np.percentile(f[s], ENV_PERCENTILE)))
    bx, by = np.asarray(bx), np.asarray(by)
    if bx.size < MIN_BINS:
        return LocalContinuum(float("nan"), float("nan"), int(bx.size), False,
                              f"CONTINUUM_UNCONSTRAINED: {bx.size} populated bins "
                              f"< {MIN_BINS}"), None
    keep = np.ones(bx.size, bool)
    for _ in range(3):
        c = np.polyfit(bx[keep] - centre_A, by[keep], DEGREE)
        r = by - np.polyval(c, bx - centre_A)
        mad = float(np.median(np.abs(r[keep] - np.median(r[keep])))) or 1e-6
        new = r > -CLIP_MAD * 1.4826 * mad
        if new.sum() < MIN_BINS or np.array_equal(new, keep):
            break
        keep = new
    c = np.polyfit(bx[keep] - centre_A, by[keep], DEGREE)
    level = float(np.polyval(c, 0.0))
    slope = float(c[0]) if DEGREE >= 1 else 0.0
    if abs(level - 1.0) > MAX_SHIFT:
        return LocalContinuum(level, slope, int(keep.sum()), False,
                              f"CONTINUUM_UNPLACEABLE: local level {level:.4f} is more "
                              f"than {MAX_SHIFT:.0%} from unity -- a telluric band or line "
                              f"forest, not a normalisation error; shipped continuum kept"), None
    if not apply:
        return LocalContinuum(level, slope, int(keep.sum()), False,
                              "MEASURED, NOT APPLIED: band declares no observable true "
                              "continuum (pseudo-continuum regime)"), None
    return (LocalContinuum(level, slope, int(keep.sum()), True, "applied"),
            np.polyval(c, w - centre_A))


def fit_model_guided(wave_A, flux, centre_A: float, model_wave_A, model_flux, *,
                     exclude_half_width_A: float,
                     env_half_width_A: float = ENV_HALF_WIDTH_A,
                     apply: bool = True,
                     reference: str = "synthesis") -> tuple[LocalContinuum, np.ndarray | None]:
    """THE STANDING RULE's estimator where a synthesis is in hand (RYA-1230, second pass).

    WHY NOT THE ENVELOPE. An absolute upper envelope is only a continuum where the true
    continuum is reached. Measured on our own atlases: at C I 5052 the p95 envelope reads
    0.005-0.012 BELOW the model-guided level, and at C I 6587 (the H-alpha wing) and N I
    8216 (the CN forest) the synthesis has NO pixel at 0.998 inside +/-5 A -- there the
    envelope divides out absorption the synthesis ALSO models, i.e. counts it twice. It
    moved A(C) VIS by about -0.14 dex before this was caught.

    THE ESTIMATOR compares like with like: observed / synthetic on the pixels where the
    SYNTHESIS is highest in the window (its top decile, `MODEL_Q`), the fitted line's own
    window (+/- `exclude_half_width_A`) excluded, straight line in lambda. Wherever the
    synthesis carries the local pseudo-continuum (line wings, molecular forests, H wings),
    the ratio cancels it and only the NORMALISATION error remains. Same MAX_SHIFT bound,
    same record, same refusal to divide by a refused level.
    """
    w = np.asarray(wave_A, float)
    f = np.asarray(flux, float)
    mw = np.asarray(model_wave_A, float)
    mf = np.asarray(model_flux, float)
    inside = (np.isfinite(w) & np.isfinite(f) & (np.abs(w - centre_A) <= env_half_width_A)
              & (np.abs(w - centre_A) > exclude_half_width_A)
              & (w >= mw.min()) & (w <= mw.max()))
    mod = np.interp(w, mw, mf)
    inside &= np.isfinite(mod) & (mod > 0)
    method = (f"RYA-1230 model-guided: obs/{reference} on the {reference}'s top {100 - MODEL_Q:g}% "
              f"pixels within +/-{env_half_width_A:g} A (line window +/-"
              f"{exclude_half_width_A:g} A excluded), degree {DEGREE}, {CLIP_MAD_MODEL:g}-MAD clipped")
    if inside.sum() < MIN_MODEL_PIX:
        return LocalContinuum(float("nan"), float("nan"), int(inside.sum()), False,
                              f"CONTINUUM_UNCONSTRAINED: {int(inside.sum())} pixels", method), None
    thr = float(np.percentile(mod[inside], MODEL_Q))
    sel = inside & (mod >= thr)
    if sel.sum() < MIN_MODEL_PIX:
        return LocalContinuum(float("nan"), float("nan"), int(sel.sum()), False,
                              f"CONTINUUM_UNCONSTRAINED: {int(sel.sum())} selected pixels", method), None
    # RYA-1232: robust fit. Anchor pixels carrying a residual (a telluric leftover after
    # correction, a cosmic, an unmodelled feature) steered the straight line: at [O I] 6300 on
    # corrected HARPS the level read 1.022 / 1.001 / 0.996 for envelopes of +4/+5/+6 A;
    # clipped at CLIP_MAD_MODEL it reads 0.992 / 1.001 / 0.999, and clean windows move < 0.1%.
    x, y = w[sel] - centre_A, f[sel] / mod[sel]
    for _ in range(5):
        c = np.polyfit(x, y, DEGREE)
        r = y - np.polyval(c, x)
        mad = 1.4826 * float(np.median(np.abs(r - np.median(r))))
        keep = np.abs(r) <= CLIP_MAD_MODEL * max(mad, 1.0e-4)
        if keep.all() or keep.sum() < MIN_MODEL_PIX:
            break
        x, y = x[keep], y[keep]
    level = float(np.polyval(c, 0.0))
    slope = float(c[0]) if DEGREE >= 1 else 0.0
    if abs(level - 1.0) > MAX_SHIFT:
        return LocalContinuum(level, slope, int(sel.sum()), False,
                              f"CONTINUUM_UNPLACEABLE: local level {level:.4f} is more than "
                              f"{MAX_SHIFT:.0%} from unity; shipped continuum kept", method), None
    if not apply:
        return LocalContinuum(level, slope, int(sel.sum()), False,
                              "MEASURED, NOT APPLIED: pseudo-continuum regime", method), None
    return (LocalContinuum(level, slope, int(len(x)), True, "applied", method),
            np.polyval(c, w - centre_A))


def apply_to_windows(wave_A, flux, windows_A, *, apply: bool = True):
    """The same rule for a route that fits several windows on one loaded band (RYA-1230:
    `pipeline.cno_synthesis`). Each window gets its OWN envelope, centred on it and
    widened by its half-width, and only the pixels inside that window (+/-1 A) are
    divided -- a window never inherits a neighbour's level.

    Returns (flux copy, [record per window]).
    """
    w = np.asarray(wave_A, float)
    out = np.asarray(flux, float).copy()
    records = []
    for lo, hi in windows_A:
        c, half = 0.5 * (lo + hi), 0.5 * (hi - lo)
        rec, cont = fit(w, out, c, env_half_width_A=half + ENV_HALF_WIDTH_A, apply=apply)
        records.append(dict(window_A=(float(lo), float(hi)), level=rec.level_at_centre,
                            slope_per_A=rec.slope_per_A, n_bins=rec.n_bins,
                            applied=rec.applied, reason=rec.reason, method=rec.method))
        if cont is not None:
            sel = (w >= lo - 1.0) & (w <= hi + 1.0)
            out[sel] = out[sel] / cont[sel]
    return out, records


#: RYA-1232 — THE IAG ATLAS AS A CONTINUUM REFERENCE (now the `contref` LEG, see CONT_REF;
#: the first pass made it nominal, which the N I A/B overturned). Measured 2026-10-02: at its own top-decile pixels our synthesis sits
#: 1-3% ABOVE the real Sun (obs/synth 0.979-0.994 where obs/IAG is 1.001-1.013 on the same
#: corrected HARPS: [O I] 6300, C I 5052, C2 5163, CN 6127) -- weak absorption missing
#: from the line list. Dividing by that "continuum" made every line shallower; on the 4%
#: [O I] 6300 blend it cost -0.75 dex. The solar literature anchors on the OBSERVED atlas
#: (Amarsi+2021: local maxima, two independent atlases; Asplund+2021: own local placement on
#: the atlases), and Jofre+2014 note the synthesis-guided placement fails where the list is
#: incomplete. The IAG atlas carries the real H-alpha wings and molecular forests, so the
#: ratio still cancels them -- the reason the rule became model-guided -- without the
#: missing-opacity bias.
#:
#: 🔴 THE REFERENCE IS USED ONLY WHERE IT PASSES THE TWO-ATLAS CHECK (Amarsi+2021 practice).
#: Measured 25 A windows, level = holding / IAG on IAG's top decile:
#:   5000-10000 A: HARPS 0.988-1.007, KP-molecfit 0.997-1.008, KP-Kurucz2005 0.993-1.006
#:                 (band medians) -- three independent reductions agree with Baker+2020 to ~1%.
#:   4000-4500 A:  HARPS 0.900, KP 0.962, K2005 0.958 against Reiners+2016 -- and HARPS vs KP
#:                 disagree by 6% with EACH OTHER. Baker+2020 (and Reiners) normalise on local
#:                 maxima; where lines never let the true continuum through (the blue, like
#:                 the near-UV of RYA-1189) that is a pseudo-continuum. No observed reference
#:                 is verified there, so it is not used.
#: So the span is Baker+2020's processed range, 5001.1-10000 A. 🔴 RYA-1232: the file runs to
#: 11086 A but Baker+2020 processed 500-1000 nm and the agreement above was only ever verified
#: to 10000 A. Beyond it IAG sits 2-7% LOW against both KP and CRIRES+ at 10490-10870 A and its
#: C I Y lines read A(C) 8.89-9.27 where KP/CRIRES+ give 8.56-8.62 on the same lines. Outside the
#: span the `contref` leg has no atlas and equals the nominal (reference spread 0, recorded).
SOLAR_ATLAS_SPAN_A = (5001.1, 10000.0)


#: 🔴 RYA-1232 FINAL RULE (A/B measured 2026-10-03): the NOMINAL reference is the
#: SYNTHESIS (3-MAD clipped); the IAG atlas is the RYA-587 `contref` LEG (CODEX_CONT_REF=iag).
#: A fit compares the observation to OUR synthesis, so the continuum must be placed
#: relative to that synthesis -- the fitting analogue of Amarsi+2020's local-maxima EWs,
#: which measure obs and model each against its own local continuum. Measured: N I
#: (AGSS21, IAG holding) 1D-LTE 7.950 on the synthesis reference vs 8.12 on the IAG
#: reference; the literature-method EW inversion gives 7.83-7.99, Lodders+2025 7.94,
#: Mashonkina+2024 7.92 (1D-NLTE). The IAG reference fixed the absolute level but let the
#: synthesis's missing weak opacity be spent on 0.3-1%-deep lines. [O I] 6300's -0.75 was
#: unclipped anchors on O2-gamma residuals; clipped it reads 8.93 (synthesis) / 8.87 (IAG).
#: The reference spread is priced, not hidden: |A(iag) - A(nominal)| enters `continuum`.
CONT_REF = (_os.environ.get("CODEX_CONT_REF") or "synthesis").strip().lower()
if CONT_REF != "synthesis":
    print(f"  \u26a0\ufe0f  CONTINUUM REFERENCE LEG (RYA-587 contref): {CONT_REF}")


def solar_atlas_reference(lo_A: float, hi_A: float, resolving_power: float):
    """IAG solar atlas over [lo_A, hi_A], broadened to `resolving_power` -> (wave_A, flux),
    or None outside SOLAR_ATLAS_SPAN_A."""
    if CONT_REF != "iag":
        return None
    if lo_A < SOLAR_ATLAS_SPAN_A[0] or hi_A > SOLAR_ATLAS_SPAN_A[1]:
        return None
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
    import measure_band_ew as _mbe                                     # noqa: E402
    from scipy.ndimage import gaussian_filter1d
    w, f = _mbe.iag_atlas()
    k = (w >= lo_A - 1.0) & (w <= hi_A + 1.0)
    w, f = w[k], f[k]
    if w.size < 50:
        return None
    step = float(np.median(np.diff(w)))
    fwhm = 0.5 * (lo_A + hi_A) / float(resolving_power)
    return w, gaussian_filter1d(f, fwhm / 2.3548 / step)


def apply_to_windows_model_guided(wave_A, flux, windows_A, model_fn, *, apply: bool = True,
                                  solar_R: float | None = None):
    """`apply_to_windows` with the model-guided estimator (`fit_model_guided`): for each
    fit window, `model_fn(lo_A, hi_A) -> (wave_A, flux)` synthesises the window +/-
    ENV_HALF_WIDTH_A at the route's current composition, and the window itself is excluded
    from the continuum pixels. Used by `pipeline.cno_synthesis`, whose CN/CH/OH windows sit
    in molecular forests where an absolute envelope double-counts the forest.
    """
    w = np.asarray(wave_A, float)
    out = np.asarray(flux, float).copy()
    records = []
    for lo, hi in windows_A:
        c, half = 0.5 * (lo + hi), 0.5 * (hi - lo)
        env = half + ENV_HALF_WIDTH_A
        ref = solar_atlas_reference(c - env - 0.5, c + env + 0.5, solar_R) if solar_R else None
        win_apply = apply
        if ref is not None:
            mw, mf, label = ref[0], ref[1], "IAG solar atlas"
        else:
            mw, mf = model_fn(c - env - 0.5, c + env + 0.5)
            label = "synthesis"
            pass                                     # no atlas here: the leg equals nominal
        mw, mf = np.asarray(mw, float), np.asarray(mf, float)
        keep = np.abs(mw - c) <= env                 # iSpec zeroes synthesis edges
        rec, cont = fit_model_guided(w, out, c, mw[keep], mf[keep],
                                     exclude_half_width_A=half, env_half_width_A=env,
                                     apply=win_apply, reference=label)

        records.append(dict(window_A=(float(lo), float(hi)), level=rec.level_at_centre,
                            slope_per_A=rec.slope_per_A, n_pix=rec.n_bins,
                            applied=rec.applied, reason=rec.reason, method=rec.method,
                            reference=label))
        if cont is not None:
            sel = (w >= lo - 1.0) & (w <= hi + 1.0)
            out[sel] = out[sel] / cont[sel]
    return out, records


def place_for_synthesis(wave_A, flux, centre_A: float, ctx: dict, element: str, *,
                        band_half_width_A: float, use_molecules: bool, apply: bool = True,
                        star: str):
    """THE standing rule for a synthesis route, in one place (RYA-1230 / RYA-1232).

    Synthesises +/- ENV_HALF_WIDTH_A around `centre_A` from the route's OWN context
    (atmosphere, line list, molecules, the target element at the context's solar value),
    then `fit_model_guided`, excluding the BAND's fit window (never a leg's, so the
    core-window leg cannot move the continuum). Returns (flux, record): the flux divided by
    the placed continuum when applied, else unchanged.
    """
    env = ENV_HALF_WIDTH_A
    if star == "solar":
        ref = solar_atlas_reference(centre_A - env - 0.5, centre_A + env + 0.5,
                                    float(ctx["resolving_power"]))
        if ref is not None:
            k = np.abs(ref[0] - centre_A) <= env
            rec, cont = fit_model_guided(wave_A, flux, centre_A, ref[0][k], ref[1][k],
                                         exclude_half_width_A=float(band_half_width_A),
                                         apply=apply, reference="IAG solar atlas")
            out = np.asarray(flux, float) / cont if cont is not None else flux
            return out, rec
        pass                                           # no atlas here: the leg equals nominal
    import os
    from pathlib import Path
    from pipeline.abundances_derive import _synth_flux_at_abund
    tmp = f"/tmp/ispec_cont_{os.getpid()}"
    Path(tmp).mkdir(parents=True, exist_ok=True)
    env = ENV_HALF_WIDTH_A
    mw = np.arange((centre_A - env - 0.5) / 10.0, (centre_A + env + 0.5) / 10.0, 0.0005)
    mf = _synth_flux_at_abund(
        mw, ctx["atmosphere"], ctx["teff"], ctx["logg"], ctx["feh"], ctx["vturb"],
        ctx["linelist"], ctx["isotopes"], ctx["solar_abund"], element,
        int(ctx["atom_code"]), float(ctx["solar_A"]),
        R=float(ctx["resolving_power"]), macroturbulence=float(ctx["macroturbulence"]),
        vsini=float(ctx["vsini"]), use_molecules=bool(use_molecules), tmp_dir=tmp)
    edge = np.abs(mw * 10.0 - centre_A) <= env          # iSpec zeroes synthesis edges
    rec, cont = fit_model_guided(wave_A, flux, centre_A, mw[edge] * 10.0, np.asarray(mf)[edge],
                                 exclude_half_width_A=float(band_half_width_A), apply=apply,
                                 reference="synthesis")
    out = np.asarray(flux, float) / cont if cont is not None else flux
    return out, rec

