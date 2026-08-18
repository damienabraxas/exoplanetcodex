"""The SYNTHETIC equivalent width of one line — one definition, one home. RYA-878.

`SynthesisHandler.measure_line` computed this inline, and the production synth-v2 path
(`abundances_derive._run_synthesis_v2_mode`) computed it not at all: it banks `A_X`,
`red_chi2` and `status` per line and no EW. That asymmetry is why the SynthesisHandler
control's ANGLE 1 was never a measured second angle —

    ANGLE 1 = handler's SYNTHETIC EW  vs  the pool's MEASURED profile-fit EW

— and the claim that the production engine shares its +0.1562 dex offset could only be
INFERRED (RYA-875), because the engine banked nothing to compare against. A comparison is
an angle only if both sides are measured the same way; that is the lesson RYA-875 closed
on ANGLE 2, where RYA-770 had held the lines fixed on one side and left the reference a
scalar from a different sample.

So the definition lives here, once, and every side of the comparison calls it. Writing a
second copy into the production path would have made the two sides free to disagree about
what "the synthetic EW" means, which is the RYA-845 defect shape and the exact thing this
ticket exists to remove.

WHY IT IS A DIFFERENCE OF TWO SYNTHESES
---------------------------------------
The window is synthesised twice, once at the fitted abundance and once with the element
depleted, and the difference is integrated. Everything that is not this element absorbs
identically in both and subtracts away, so other species cancel exactly and no blend
model is needed.

⚠️ ANOTHER LINE OF THE SAME ELEMENT DOES NOT CANCEL — it vanishes from the depleted
synthesis too, so its absorption lands in the difference. RYA-873 tested whether that
explains ANGLE 1 and REFUTED it (15 of 18 lines carry no predicted contamination and
still span observed ratios 0.50–2.07), but the property is real and belongs in the
docstring of the thing that has it.

⚠️ THE INTEGRATION RANGE IS DELIBERATELY WIDER THAN THE FIT WINDOW. The fit window is
tight so the line dominates chi2; integrating the EW over that same window truncates the
damping wings and drove the EW ratio to 0.250 by construction. Since the difference is
zero wherever the element does not absorb, widening adds wings and nothing else.
"""
from __future__ import annotations

import numpy as np

#: Integration half-width, in Å: three fit half-widths, floored so a narrow line still
#: gets its wings. The floor is what makes a weak line's range wide enough to close the
#: profile rather than clipping it at the fit bound.
EW_HALF_WIDTH_FLOOR_A = 1.0
EW_HALF_WIDTH_FACTOR = 3.0

#: Sampling of the integration grid, in Å, and the minimum number of points. 0.004 Å is
#: well inside the narrowest solar line core at optical resolutions, so the trapezoid is
#: not the limiting error.
EW_STEP_A = 0.004
EW_MIN_POINTS = 200

#: How far the element is depleted to build the "without" spectrum. Four dex removes it
#: entirely for any line that matters; it is not a fit parameter and nothing is tuned to
#: it — at this depth the residual absorption of the element is far below the trapezoid
#: error, so the difference is the line.
EW_DEPLETION_DEX = 4.0


def ew_half_width_A(fit_half_width_A: float) -> float:
    """The EW integration half-width implied by a fit half-width."""
    return max(EW_HALF_WIDTH_FACTOR * float(fit_half_width_A), EW_HALF_WIDTH_FLOOR_A)


def synthetic_ew_mA(synth, *, centre_A: float, abundance: float,
                    fit_half_width_A: float, synth_kwargs: dict) -> float:
    """Synthetic EW of the target element's absorption around `centre_A`, in mÅ.

    `synth(wave_nm, trial_A=..., **synth_kwargs)` is the project's ONE flux generator
    (`abundances_derive._synth_flux_at_abund`); it is passed in rather than imported so
    this module stays free of the import cycle between the handler and the deriver, and
    so a caller cannot accidentally use a different generator than the one that fitted
    the line.

    Returns NaN if the synthesis raises — the caller decides what an unmeasurable EW
    means, and for the control it means the line reports no EW rather than a zero that
    would look like a measurement.
    """
    from pipeline._numcompat import trapezoid as _trapz
    hw = ew_half_width_A(fit_half_width_A)
    n = max(int(2 * hw / EW_STEP_A), EW_MIN_POINTS)
    w_A = np.linspace(centre_A - hw, centre_A + hw, n)
    try:
        with_el = synth(w_A / 10.0, trial_A=float(abundance), **synth_kwargs)
        without = synth(w_A / 10.0, trial_A=float(abundance) - EW_DEPLETION_DEX,
                        **synth_kwargs)
    except Exception:
        return float("nan")
    return float(_trapz(np.clip(without - with_el, 0.0, None), w_A)) * 1000.0
