"""Did the objective actually pin the abundance? — RYA-847.

THE DEFECT THIS EXISTS TO CLOSE
-------------------------------
RYA-843 found that the band-product synthesis route reimplements accept/reject as
`status != 'ok'` and throws `red_chi2` away, so the RYA-342 merit gate never reaches it.
The consequence is not the loud one:

    UNCONSTRAINED FITS ARE ACCEPTED, AND ONLY THE ONES THAT LAND SOMEWHERE
    IMPLAUSIBLE ARE VISIBLE.

Two NIR lines entered the published aggregate at 7.833 and 7.979 — perfectly plausible
numbers — whose chi2 moves 2.2% and 1.4% across EIGHT DEX of iron. They measured nothing
and landed somewhere believable by luck. A railed fit at A_solar+5 is the same defect
wearing a costume loud enough to notice.

`edge_pinned` cannot express this. It asks "is the answer near a bound", which is a
question about the bracket, not about whether the data determined anything. RYA-843
measured how sharp that is: `a_hi` is 12.496 exactly, and two lines both PRINTING 12.486
fell on opposite sides of `< 1e-2` — the verdict was decided in the 4th decimal of a
rounded number.

WHAT THIS MODULE DOES, AND DELIBERATELY DOES NOT DO
---------------------------------------------------
It MEASURES three candidate constraint metrics per line and takes no view on which is
right or where the cut belongs. That is not indecision, it is the ticket's instruction
(RYA-847 item 1 as refined, and RYA-161): the threshold and the metric SHAPE are chosen
from a sweep across every synthesis band, never from one product.

    sigma_A          curvature at the minimum, chi2 rescaled so red_chi2 == 1, in DEX.
                     BRACKET-FREE and A_solar-FREE. Delta-chi2 is a DIFFERENCE, so the
                     depressed-baseline offset that makes an ABSOLUTE red_chi2
                     untransferable cancels exactly -- that is the property that lets one
                     number cross bands. Measures PRECISION, not accuracy: a sharp
                     minimum in the wrong place still scores well.

    frac_rise        fractional chi2 rise from the minimum to each end of the caller's
                     bracket. Separated cleanly on RYA-843's nine lines, but it is
                     CALLER-DEPENDENT (handlers bracket differently, so it does not
                     transfer) and it partly measures PROXIMITY TO THE BRACKET EDGE --
                     with a bracket of [A_solar-3, A_solar+5] that is a plausibility
                     prior tied to A_solar entering through the back door (RYA-161).
                     Carried for comparison, with its leak stated.

    edge_distance    what `edge_pinned` used, in dex rather than as a boolean. Carried
                     for completeness so the sweep can show what the old criterion would
                     have decided.

COST
----
Four extra chi2 evaluations per line: two curvature probes and the two bracket ends. On
the measured Fe windows that is roughly a 15% overhead on a fit, which is why the metrics
are computed inside the fitter rather than by re-probing lines afterwards.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np

#: Curvature-probe step, in dex. Small enough to stay in the parabolic neighbourhood of
#: the minimum, large enough that Delta-chi2 clears numerical noise. Shared by every
#: caller so two paths cannot probe at different scales and be compared as if they had
#: not (RYA-845: declare it once).
CURVATURE_PROBE_STEP_DEX = 0.05


def as_float_or_none(v) -> float | None:
    """float, or None for anything unmeasurable. Shared so callers cannot disagree.

    None and NaN are NOT interchangeable in an artifact. None reads as "not applicable" —
    an EW inversion has no chi2 surface, so the constraint question does not apply to it.
    NaN reads as "a measurement that failed", which would put every EW-route line under
    permanent suspicion. Both of this module's consumers had started writing their own
    version of this, which is how two columns come to mean slightly different things.
    """
    try:
        f = float(v)
    except (TypeError, ValueError):
        return None
    return f if np.isfinite(f) else None


@dataclass(frozen=True)
class ConstraintMetrics:
    """Everything measured about whether the fit determined the abundance.

    No verdict. `sigma_A` NaN means the curvature was not measurable -- a railed fit, or
    an objective with no curvature to invert -- which is itself the strongest statement
    available and must never be rendered as a number (RYA-848: the old code returned
    0.000 for a railed fit, the tightest possible bar for the least constrained possible
    fit, and clipped a flat objective to a plausible-looking 1.000).
    """
    sigma_A: float
    frac_rise_lo: float
    frac_rise_hi: float
    frac_rise_weaker: float
    edge_distance_dex: float
    red_chi2: float
    n_probe_evals: int

    def as_dict(self) -> dict:
        return asdict(self)


def curvature_sigma(chi2_fn, *, a_best: float, chi2_min: float, red_chi2: float,
                    a_lo: float, a_hi: float, edge: bool,
                    step: float = CURVATURE_PROBE_STEP_DEX) -> float:
    """1 sigma on the abundance from the local chi2 curvature. NaN when not measurable.

    Hoisted here by RYA-847 from `cno_synthesis`, where RYA-848 fixed it; it is the
    published sigma_stat for CNO's N and O, and it is now also the constraint metric for
    every synthesis band. Two paths computing this separately is exactly how the three
    accept/reject implementations drifted apart in the first place.

    (a) DELTA-CHI2 = 1 MEANS ONE SIGMA ONLY IF CHI2 IS CALIBRATED, AND IT IS NOT. The
        per-pixel sigma on both synthesis paths is a 0.01 model-adequacy PLACEHOLDER, so
        the raw Delta-chi2 scale is arbitrary: measured red_chi2 runs 2.6-1226 on Fe
        windows and 17-67 on CNO bands. A chi2 N times too large makes Delta-chi2 N times
        too large and sigma sqrt(N) times too SMALL. Rescaling to red_chi2 == 1 is the
        standard correction, and because Delta-chi2 is a DIFFERENCE the constant
        model-inadequacy offset cancels rather than propagating.

    (b) PROBE BOTH SIDES, AND NEVER CLAMP. The pre-RYA-848 probe was
        `min(a_best + step, a_hi)`: one-sided and clamped, so a fit sitting ON the bound
        got a zero-length step and reported sigma = 0.000. Out-of-bracket probes are
        SKIPPED here; a railed fit returns NaN, because on the rail side the abundance is
        not bounded at all and there is no curvature to invert.

    The worse-constrained side is taken rather than the average: for a true parabola they
    agree, and where they do not, the honest sigma is the larger one.
    """
    if edge:
        return float('nan')
    if not np.isfinite(red_chi2) or red_chi2 <= 0.0:
        return float('nan')
    sig = []
    for probe in (a_best - step, a_best + step):
        if not (a_lo <= probe <= a_hi):
            continue
        d_raw = float(chi2_fn(probe)) - float(chi2_min)
        if d_raw <= 0.0:
            continue
        sig.append(abs(probe - a_best) / np.sqrt(d_raw / red_chi2))
    return float(max(sig)) if sig else float('nan')


def measure_constraint(chi2_fn, *, a_best: float, chi2_min: float, red_chi2: float,
                       a_lo: float, a_hi: float,
                       step: float = CURVATURE_PROBE_STEP_DEX,
                       edge_tol_dex: float = 1e-2) -> ConstraintMetrics:
    """All three candidate metrics, measured on the caller's own objective.

    `chi2_fn` must be THE FITTER'S OWN objective, not a reconstruction of it. RYA-843
    captured it by spying on `minimize_scalar` for exactly this reason: a rebuilt chi2
    agrees with itself while the product does something else (RYA-845).
    """
    n0 = [0]

    def _chi2(a):
        n0[0] += 1
        return float(chi2_fn(a))

    edge_distance = float(min(abs(a_best - a_lo), abs(a_best - a_hi)))
    edge = edge_distance < edge_tol_dex

    sigma_A = curvature_sigma(_chi2, a_best=a_best, chi2_min=chi2_min,
                              red_chi2=red_chi2, a_lo=a_lo, a_hi=a_hi,
                              edge=edge, step=step)

    # Fractional rise to the bracket ends. Reported as a RATIO to chi2_min so it is
    # dimensionless, which is also precisely why it cannot be compared between callers
    # whose brackets differ -- stated in the module docstring, carried anyway so the
    # sweep can show what it would have decided.
    def _frac(a_end):
        if not np.isfinite(chi2_min) or chi2_min <= 0.0:
            return float('nan')
        return (_chi2(a_end) - chi2_min) / chi2_min

    lo = _frac(a_lo)
    hi = _frac(a_hi)
    weaker = float(np.nanmin([lo, hi])) if np.isfinite([lo, hi]).any() else float('nan')

    return ConstraintMetrics(
        sigma_A=float(sigma_A), frac_rise_lo=float(lo), frac_rise_hi=float(hi),
        frac_rise_weaker=weaker, edge_distance_dex=edge_distance,
        red_chi2=float(red_chi2), n_probe_evals=int(n0[0]))
