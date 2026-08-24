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

🔴 THE CUT IS PER ARM, AND NOT FOR THE REASON RYA-995 GAVE.

RYA-995 measured `frac_rise_weaker` sitting uniformly LOWER on HARPS than on Kitt Peak
(median 3.26 vs 4.81 deep, 1.82 vs 2.46 shallow) and read that as a flat 1.0 OVER-rejecting
HARPS. `scripts/rya992_arm_scale.py` tests the inference on the RYA-986 pool — the SAME
1483 Fe I lines through both arms, so the line list, the gf values and the atmosphere are
identical and the only difference between the two columns is the instrument. It comes out
the other way round:

    at the SAME frac_rise, a HARPS line is never BETTER than a Kitt Peak line

      frac_rise bin   KP med |A-median A|   HARPS med |A-median A|
        [0.0, 0.5)         0.2840                 0.3110
        [0.5, 1.0)         0.1545                 0.1820
        [1.0, 1.5)         0.1490                 0.1595
        [1.5, 2.0)         0.1000                 0.1490
        [2.0, 3.0)         0.1150                 0.1090
        [3.0, 5.0)         0.0990                 0.1630

HARPS's frac_rise is lower by x1.27 — and its `sigma_A` is WORSE by x1.29. The arm that
looks over-rejected is the arm whose fits actually pin A less well. A flat 1.0 was UNDER-
rejecting HARPS, not over-rejecting it, and lowering the HARPS cut to equalise the rejected
FRACTION would have been the RYA-981 error in its purest form: a quota dressed as a cut.

🔴 WHAT IS ACTUALLY ARM-DEPENDENT: A WINDOW-GEOMETRY FACTOR, AND IT IS MEASURABLE.

`fit_constraint` builds the two statistics from one chi2 surface:

    sigma_A   = step / sqrt(dchi2_raw(step) / red_chi2)      calibrated, RYA-848
    frac_rise = dchi2_raw(bracket end) / chi2_min            NOT calibrated

`chi2_min = red_chi2 * dof`, so red_chi2 cancels out of frac_rise entirely and, for a
locally parabolic surface, frac_rise = B^2 / (sigma_A^2 * dof) with B the bracket half-
width. Everything in frac_rise that is not the abundance constraint collects into

    k = frac_rise * sigma_A^2  =  B^2 / dof

— bracket width over the number of fitted pixels. That is pure window geometry: it says how
much fractional chi2 rise an arm hands you per unit of abundance constraint, and it is why
`fit_constraint` already warns that frac_rise "cannot be compared between callers whose
brackets differ". MEASURED on the RYA-986 pools:

    kpno_solar_atlas   k = 0.00415   (bootstrap sd 3.7%)
    harps              k = 0.00564   (bootstrap sd 3.6%)
    ARM SCALE  k(harps)/k(kpno) = 1.357

k itself moves with the line selection — it is 1.90x larger on the laboratory-graded subset
of either arm, because those windows are narrower — but that factor is IDENTICAL on both
arms, so the RATIO is 1.357 on the graded subset too. The ratio, and only the ratio, is a
property of the instrument. That is what `ARM_SCALE` stores.

🔴 CONFIRMED A SECOND WAY, WITH NO PARABOLA ASSUMED. Fit a declining segment meeting a flat
one to each arm's own err-vs-log10(frac_rise) relation and take the breakpoint:

    kpno_solar_atlas   knee frac_rise 0.126   bootstrap 16-84% ~[0.06, 0.16]
    harps              knee frac_rise 0.158   bootstrap 16-84% ~[0.11, 0.22]
    knee(harps)/knee(kpno) = 1.26   against the geometry ratio 1.357

(the breakpoint is scanned on a 0.05-dex grid, so the bootstrap interval moves by a bin
between runs; the point estimates and their ratio do not)

Same direction, same size, measured two ways. ⚠️ THE KNEE'S ABSOLUTE POSITION IS NOT USED,
and saying why is part of the result: a knee sits where fit-driven error drops below the
population's irreducible gf scatter, so it lands at ~0.13 on the full pool and at ~2-3 on
the graded anchor. That is a property of the POOL, not of the arm. Only the ratio between
two arms measured on the SAME lines divides that floor out.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

#: 🔴 THE REFERENCE CUT, DERIVED FROM THE OBJECTIVE AND NOT FROM A PERCENTILE (RYA-981).
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
#:
#: 🔴 IT IS THE CUT **ON KITT PEAK**, because that is the arm the fixture was measured on.
#: It is not a project-wide constant; every other arm derives its own from `ARM_SCALE`.
SYNTH_GOF_REF_CUT = 1.0


@dataclass(frozen=True)
class ArmScale:
    """How much fractional chi2 rise one arm hands you per unit of abundance constraint.

    `scale` multiplies `SYNTH_GOF_REF_CUT`. It is k(arm)/k(`relative_to`) with
    k = frac_rise_weaker * sigma_A**2, measured on lines COMMON to both arms — see the
    module docstring, and `scripts/rya992_arm_scale.py` for the derivation.

    🔴 `relative_to` IS PART OF THE MEASUREMENT, WHICH IS WHY IT LIVES ON THE ROW AND NOT
    IN A MODULE CONSTANT. A ratio without its denominator is not a number, so every row
    names the arm it was measured against; the reference arm is then DERIVED from the
    registry (`reference_arm()`) rather than declared once and inherited. That also keeps
    RYA-913/922 satisfied — no module-level instrument constant exists here, and
    `synth_gof_cut` takes its arm from the caller.
    """
    instrument: str
    scale: float
    relative_to: str
    measured_on: str
    note: str


#: 🔴 A REGISTRY, NOT A FORMULA. Each row is a MEASUREMENT against the reference arm on a
#: named pair of products, dated by the ticket that made them. Add a row; do not edit one.
#:
#: ⚠️ AN UNKNOWN INSTRUMENT RAISES. It does not fall back to 1.0. Defaulting an unmeasured
#: arm to the Kitt Peak scale is exactly the transfer this ticket exists to stop, and it
#: would do it SILENTLY (RYA-833: an absence is a hypothesis, never a conclusion).
ARM_SCALE: dict[str, ArmScale] = {
    "kpno_solar_atlas": ArmScale(
        instrument="kpno_solar_atlas", scale=1.0, relative_to="kpno_solar_atlas",
        measured_on="RYA-992 — the arm the reference fixture was measured on",
        note=("Scale 1.0 BY DEFINITION, not by measurement: this is the arm the others "
              "are measured against, so measuring it against itself would report only "
              "the noise in its own k. `relative_to` naming itself is what marks it as "
              "the reference — see `reference_arm()`."),
    ),
    "harps": ArmScale(
        instrument="harps", scale=1.357, relative_to="kpno_solar_atlas",
        measured_on=("RYA-986 VIS pools, 1483 Fe I lines common to "
                     "FeI_4200_6910_kpno_solar_atlas_SYNTH_FROMEW_1D-LTE_lines.csv and "
                     "FeI_4200_6910_harps_solar_harps_molecfit_corrected_SYNTH_FROMEW"
                     "_1D-LTE_lines.csv"),
        note=("k 0.00564 against Kitt Peak's 0.00415, bootstrap sd 3.6%/3.7%. The same "
              "1.357 comes back on the laboratory-graded subset of each arm, where both "
              "k values are 1.90x larger — the line-selection factor cancels, the arm "
              "factor does not. Cross-checked by the err-vs-frac_rise knee ratio, 1.25. "
              "🔴 THIS RAISES THE HARPS BAR (1.0 -> 1.357). RYA-995 expected it to fall; "
              "at a given frac_rise a HARPS line is never better than a Kitt Peak one."),
    ),
}


def reference_arm() -> str:
    """The arm every scale is quoted against — DERIVED from the registry, never declared.

    The reference is the row that is its own denominator. Deriving it means the registry
    cannot disagree with a constant about which arm is the reference, and it means adding
    a second reference is a loud error rather than a silent re-anchoring of every cut.
    """
    refs = sorted(k for k, v in ARM_SCALE.items() if v.relative_to == k)
    if len(refs) != 1:
        raise ValueError(
            f"the stage-4 arm registry must have exactly ONE reference arm (a row whose "
            f"`relative_to` is itself); found {refs}. Every other scale is a ratio "
            f"against it, so two references — or none — makes every cut ambiguous.")
    return refs[0]


def synth_gof_cut(instrument: str) -> float:
    """The stage-4 frac_rise cut FOR THIS ARM.

    Loud on an unmeasured instrument. A grading run that cannot say which arm it is on
    cannot apply this gate, and guessing would reintroduce the transferred constant.
    """
    key = str(instrument)
    if key not in ARM_SCALE:
        raise KeyError(
            f"no measured stage-4 arm scale for instrument {key!r}; measured: "
            f"{sorted(ARM_SCALE)}. `frac_rise_weaker` carries a window-geometry factor "
            f"(k = frac_rise * sigma_A**2 = B**2/dof) that differs between instruments, "
            f"so the Kitt Peak cut of {SYNTH_GOF_REF_CUT} does not transfer. Measure the "
            f"scale with scripts/rya992_arm_scale.py against a pool of lines this arm "
            f"shares with {reference_arm()!r}, and add a row to ARM_SCALE.")
    return SYNTH_GOF_REF_CUT * ARM_SCALE[key].scale


def measure_arm_scale(frac_rise_weaker, sigma_A, ref_frac_rise_weaker, ref_sigma_A) -> float:
    """k(arm)/k(reference arm), the number that goes in an `ARM_SCALE` row.

    Pass the two arms' columns FOR THE SAME LINES. k moves with the line selection — the
    laboratory-graded subset runs 1.90x higher than the full pool on both arms — so a ratio
    taken across two different line sets measures the selection, not the instrument.

    🔒 RYA-161: takes fit metrics only. No abundance, published or otherwise, can reach it.
    """
    def _k(f, s):
        f = np.asarray(f, dtype=float)
        s = np.asarray(s, dtype=float)
        if f.shape != s.shape:
            raise ValueError(f"frac_rise and sigma_A differ in length: {f.shape} vs {s.shape}")
        m = np.isfinite(f) & np.isfinite(s) & (f > 0) & (s > 0)
        if m.sum() < 30:
            raise ValueError(
                f"only {int(m.sum())} rows carry a positive frac_rise AND a finite sigma_A; "
                f"a median k on that few is noise, and an arm scale is multiplied into "
                f"every line this arm will ever grade.")
        return float(np.median(f[m] * s[m] ** 2))

    return _k(frac_rise_weaker, sigma_A) / _k(ref_frac_rise_weaker, ref_sigma_A)


def synth_gof_ok(frac_rise_weaker: float | None, instrument: str) -> bool:
    """Did this synthesis fit pin the abundance well enough to grade on?

    `instrument` is REQUIRED and has no default. The cut is a property of the arm, and a
    default would be the fixed 1.0 coming back through the signature.

    `None` is NOT a pass. A row without the statistic is a row whose fit was never scored,
    and admitting it would let the un-scored population through the one stage that exists
    to score them (RYA-833: an absence is a hypothesis, never a conclusion).
    """
    cut = synth_gof_cut(instrument)      # first, so an unknown arm raises even on None
    if frac_rise_weaker is None:
        return False
    v = float(frac_rise_weaker)
    if not math.isfinite(v):
        return False
    return v >= cut


def synth_gof_reason(frac_rise_weaker: float | None, instrument: str) -> str:
    """The coded refusal, so a rejected line carries WHY (RYA-711 quarantine-not-cull)."""
    cut = synth_gof_cut(instrument)
    if frac_rise_weaker is None or not math.isfinite(float(frac_rise_weaker)):
        return ("SYNTH-GOF-UNSCORED: this row carries no frac_rise_weaker, so the fit was "
                "never scored. Not graded — an unscored fit is not a passing one.")
    v = float(frac_rise_weaker)
    return (f"SYNTH-GOF: chi2 rises only {v:.2f}x toward weaker abundances, below the "
            f"{cut:.3f}x this arm needs — {SYNTH_GOF_REF_CUT:.1f}x, the point at which the "
            f"objective merely doubles, scaled by {instrument}'s measured window-geometry "
            f"factor {ARM_SCALE[str(instrument)].scale:.3f}. The fit does not pin A(X) "
            f"tightly enough for its abundance to be graded (RYA-992; red_chi2 cannot make "
            f"this call on the synthesis route — RYA-981).")
