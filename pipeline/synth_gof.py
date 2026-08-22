"""Stage-4 goodness-of-fit for the SYNTHESIS route — RYA-992.

RYA-968's grading tree gates a line at stage 4 on whether the fit that produced its
abundance was any good. On the EW route that is `red_chi2`. On the synthesis route it could
not be, and RYA-981 measured why: corr(red_chi2, |A - mean A|) = 0.025. A synthesis chi2 is
dominated by how well the whole window matches — continuum, blanketing, broadening — and
almost none of that is about the target line's abundance. Gating on it would grade on noise.

🔴 WHAT THE SCREEN FOUND, ON THE RYA-984 163-LINE GRADED FIXTURE.

    statistic            spearman vs |A - median A|      p
    red_chi2                   +0.106                 0.178     no
    sigma_A                    +0.148                 0.059     no
    edge_distance_dex          -0.030                 0.707     no
    frac_rise_weaker           -0.200                 0.010     YES
    cited_sigma_dex            +0.245                 0.002     YES, but see below

`frac_rise_weaker` wins. It is how far chi2 climbs as the abundance is driven WEAKER —
i.e. how hard the objective fights a change in A. That is the quantity stage 4 is actually
asking about: not "does the model match the spectrum" but "did this fit pin A".

🔴 THE CONFOUND THAT ALMOST MADE THIS WRONG, AND WHICH CORRECTS RYA-981.

`cited_sigma_dex` — the laboratory gf uncertainty — predicts the error too, and more
strongly. That is not a fit statistic: it says a line whose gf is poorly known yields a
poorer abundance, which is true and irreducible. Part of |A - median A| on this fixture is
GF SCATTER THAT NO GOODNESS-OF-FIT STATISTIC CAN PREDICT, and it caps how well any
candidate can score.

Regressing the gf term out first changes the ranking:

    frac_rise_weaker   vs gf-corrected error   rho -0.158   p 0.044   still predicts
    red_chi2           vs gf-corrected error   rho +0.176   p 0.024   ALSO predicts
    sigma_A            vs gf-corrected error   rho +0.144   p 0.067   no

So **RYA-981's "red_chi2 fails on synth" is too strong**: it does carry signal, which the
gf scatter was masking in the raw correlation. `frac_rise_weaker` is still preferred, on two
grounds that are not about the correlation coefficient:

  * IT IS OFFSET-INVARIANT. It is a property of the chi2 surface's SHAPE under a change in
    A, so shifting every abundance by a constant leaves it identical. `red_chi2` is not, and
    RYA-992 requires the non-invariant stage to owe a per-line audit list — the one
    non-invariant stage must not become the unaudited backdoor.
  * IT IS ALREADY ON EVERY SYNTH ROW, including ungraded ones. `cited_sigma_dex` exists only
    for graded lines BY DEFINITION, so it cannot grade the RYA-986 ungraded pool at all —
    the population stage 4 exists to classify.

⚠️ `cited_sigma_dex` is therefore NOT adopted as the GoF statistic despite the strongest
correlation. It answers a different question, and on the pool this is for it is undefined.
"""
from __future__ import annotations

import math

#: 🔴 DERIVED FROM THE OBJECTIVE, NOT FROM A PERCENTILE (the RYA-981 lesson).
#:
#: `frac_rise_weaker` is a FRACTIONAL rise in chi2, so 1.0 is the point at which the
#: objective merely DOUBLES as the abundance is driven weaker. Below that the fit barely
#: resists a change in A at all, and "the abundance this fit reports" stops meaning much.
#: The cut is that statement about the objective — not a quota, and not tuned to reject a
#: chosen fraction.
#:
#: MEASURED on the 163-line graded fixture, it separates decisively:
#:     frac_rise <  1.0 : n= 12  median |err| 0.2220  scatter 0.2513
#:     frac_rise >= 1.0 : n=151  median |err| 0.0940  scatter 0.1225
#:     Mann-Whitney p = 1.2e-4
#: and it costs 7.4% of the anchor, taking its scatter 0.1992 -> 0.1720.
#:
#: ⚠️ THE BREAK IS NOT RAZOR-SHARP, AND SAYING SO IS PART OF THE RESULT. Binned, the error
#: falls 0.222 / 0.132 / 0.064 across [0,1) / [1,2) / [2,3) and is then flat and noisy
#: (0.114 / 0.095 / 0.078 / 0.093). The data supports "at least 1"; it does not sharply
#: distinguish 1 from 2. 1.0 is chosen because it is the value with a meaning, and a cut at
#: 2.0 would reject 23% of the anchor on a distinction the fixture cannot resolve.
SYNTH_GOF_MIN_FRAC_RISE = 1.0


def synth_gof_ok(frac_rise_weaker: float | None) -> bool:
    """Did this synthesis fit pin the abundance well enough to grade on?

    `None` is NOT a pass. A row without the statistic is a row whose fit was never scored,
    and admitting it would let the un-scored population through the one stage that exists
    to score them (RYA-833: an absence is a hypothesis, never a conclusion).
    """
    if frac_rise_weaker is None:
        return False
    v = float(frac_rise_weaker)
    if not math.isfinite(v):
        return False
    return v >= SYNTH_GOF_MIN_FRAC_RISE


def synth_gof_reason(frac_rise_weaker: float | None) -> str:
    """The coded refusal, so a rejected line carries WHY (RYA-711 quarantine-not-cull)."""
    if frac_rise_weaker is None or not math.isfinite(float(frac_rise_weaker)):
        return ("SYNTH-GOF-UNSCORED: this row carries no frac_rise_weaker, so the fit was "
                "never scored. Not graded — an unscored fit is not a passing one.")
    v = float(frac_rise_weaker)
    return (f"SYNTH-GOF: chi2 rises only {v:.2f}x toward weaker abundances, below the "
            f"{SYNTH_GOF_MIN_FRAC_RISE:.1f}x at which the objective merely doubles. The "
            f"fit does not pin A(X) tightly enough for its abundance to be graded "
            f"(RYA-992; red_chi2 cannot make this call on the synthesis route — RYA-981).")
