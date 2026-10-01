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
                     apply: bool = True) -> tuple[LocalContinuum, np.ndarray | None]:
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
    method = (f"RYA-1230 model-guided: obs/synth on the synthesis's top {100 - MODEL_Q:g}% "
              f"pixels within +/-{env_half_width_A:g} A (line window +/-"
              f"{exclude_half_width_A:g} A excluded), degree {DEGREE}")
    if inside.sum() < MIN_MODEL_PIX:
        return LocalContinuum(float("nan"), float("nan"), int(inside.sum()), False,
                              f"CONTINUUM_UNCONSTRAINED: {int(inside.sum())} pixels", method), None
    thr = float(np.percentile(mod[inside], MODEL_Q))
    sel = inside & (mod >= thr)
    if sel.sum() < MIN_MODEL_PIX:
        return LocalContinuum(float("nan"), float("nan"), int(sel.sum()), False,
                              f"CONTINUUM_UNCONSTRAINED: {int(sel.sum())} selected pixels", method), None
    c = np.polyfit(w[sel] - centre_A, f[sel] / mod[sel], DEGREE)
    level = float(np.polyval(c, 0.0))
    slope = float(c[0]) if DEGREE >= 1 else 0.0
    if abs(level - 1.0) > MAX_SHIFT:
        return LocalContinuum(level, slope, int(sel.sum()), False,
                              f"CONTINUUM_UNPLACEABLE: local level {level:.4f} is more than "
                              f"{MAX_SHIFT:.0%} from unity; shipped continuum kept", method), None
    if not apply:
        return LocalContinuum(level, slope, int(sel.sum()), False,
                              "MEASURED, NOT APPLIED: pseudo-continuum regime", method), None
    return (LocalContinuum(level, slope, int(sel.sum()), True, "applied", method),
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
