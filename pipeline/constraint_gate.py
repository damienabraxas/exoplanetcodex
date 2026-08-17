"""Should this fit enter the aggregate? — RYA-847, the DECIDER.

Deliberately a different module from `pipeline.fit_constraint`, which MEASURES. That
separation is enforced by a test: `fit_constraint` may not contain a threshold at all. The
reason is RYA-161 — when the measurer and the cut live together, the cut gets chosen from
whatever product is in front of whoever is editing, and the metric that justified it is
never re-examined.

    fit_constraint.py    measures sigma_A / frac_rise / edge_distance. No verdict.
    constraint_gate.py   applies the RATIFIED cut, or refuses to invent one.

WHAT THE GATE IS: A CORRECTNESS CHECK, AND NO THRESHOLD AT ALL
--------------------------------------------------------------
    frac_rise <= 0   the objective did not rise away from the point the optimizer
                     reported as its minimum. A minimum that is not a minimum is not a
                     measurement, whatever value it carries.

There is no number to choose. Zero is the boundary between "chi2 rose" and "it did not",
so the criterion has nothing to tune and no A_solar back-door.

🔴 AND THERE IS NO METRIC THRESHOLD, WHICH IS A MEASURED CONCLUSION, NOT AN OMISSION.
RYA-847 swept 9 cells and 581 synthesis lines looking for the universal gap item 3 asked
for. It does not exist:

    metric               per-band cut spread
    sigma_A                    x39.5
    frac_rise_weaker          x130.7
    edge_distance_dex          x24.4
    red_chi2                  x489.4

and the widest gap WITHIN almost every band is x1.02-2.88 — the distributions are
CONTINUOUS, not bimodal, so there is no separable population of unconstrained lines to
threshold against. RYA-342's clean-gap method (solar Fe II: clean <= 3.0, blends >= 128)
reproduces nowhere. RYA-843's 8%->47% gap in nine NIR lines was noise, exactly as that
ticket feared when it forbade setting the cut from that set.

"No defensible cut, measured" is a stronger statement than picking the least-bad number,
so `SYNTH_CONSTRAINT` stays None PERMANENTLY rather than pending. The machinery to apply
one remains, because a future band might genuinely be bimodal — but it must be earned the
same way, and the sweep is the standing evidence that today it cannot be.

WHAT THIS GATE DOES NOT CATCH, AND WHY THAT IS NOT FIXED HERE
--------------------------------------------------------------
A fit can find a genuine, sharp minimum of a model that cannot reproduce the data —
Fe I 11593.588 sits at A = 10.559 with red_chi2 = 1226 and a real minimum, so the sign
test passes it. Those are MODEL-INADEQUATE, equally non-measurements, and equally must not
be aggregated. But catching them needs a red_chi2 criterion, and red_chi2 is the metric
whose per-band cut spans x489 — the least transferable of the four. So they are routed to
`problem_children.csv` case by case under RYA-844 with a stated reason ("model inadequate,
red_chi2 >> band"), where a human names each one and a skeptic can reproduce the cut from
the reason. Auto-gating them would mean inventing exactly the threshold the sweep refuted.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

#: Which metric each name refers to on a `ConstraintMetrics`, and which direction fails.
#: Declared here so a ratified cut names its metric explicitly and cannot be applied to
#: the wrong one by position.
_DIRECTION = {
    "sigma_A": "above",              # dex; a large 1-sigma means A(X) was not pinned
    "frac_rise_weaker": "below",     # ratio; a small rise means chi2 is flat
    "edge_distance_dex": "below",    # dex; what edge_pinned used
    "red_chi2": "above",             # ratio; goodness of fit, NOT constraint
}


class ConstraintGateError(ValueError):
    """A cut was declared that this module cannot apply as stated."""


@dataclass(frozen=True)
class GateVerdict:
    """Whether the line may aggregate, and the sentence that says why.

    `reason` is empty ONLY when the line passes a ratified cut. When no cut is ratified
    it still carries text, because "not gated" is a fact about the product that belongs
    in the product (RYA-429: reported, never silently dropped).
    """
    ok: bool
    reason: str
    ratified: bool


def _cut():
    """The ratified cut as (metric, maximum) or None. Read lazily so a test can patch it."""
    from config import constants
    return getattr(constants, "SYNTH_CONSTRAINT", None)


def describe() -> str:
    """One sentence for the product provenance. Never silent about an unratified gate."""
    c = _cut()
    base = ("CONSTRAINT GATE: the NON-MINIMUM check is applied — a line whose chi2 does "
            "not rise away from its reported minimum (frac_rise <= 0) did not determine "
            "A(X) and is excluded from the aggregate, measured and retained (RYA-711). "
            "This is a correctness check with no threshold to choose.")
    if not c:
        return base + (
            " NO METRIC THRESHOLD IS APPLIED, and that is a measured conclusion rather "
            "than an omission: RYA-847 swept 9 cells and 581 synthesis lines and found "
            "no transferable cut for sigma_A, frac_rise, edge distance or red_chi2 — "
            "per-band cuts span x24 to x489 and every within-band distribution is "
            "continuous (widest gap x1.02-2.88). Model-inadequate fits that ARE real "
            "minima are routed to problem_children case by case under RYA-844 instead.")
    metric, maximum = c
    return (f"CONSTRAINT GATE: {metric} {_DIRECTION[metric]} {maximum} excludes the line "
            f"from the aggregate (measured and retained, never dropped — RYA-711).")


def verdict(metrics) -> GateVerdict:
    """Apply the ratified cut to one line's metrics.

    `metrics` is a `ConstraintMetrics`, or any object carrying the named attributes, or a
    mapping — the band-product route flattens them onto a `LineMeasurement` and the
    handler keeps the dataclass, and neither should have to convert for the other.
    """
    # ── 1. THE NON-MINIMUM CHECK — ALWAYS ON, AND NOT A THRESHOLD ────────────
    # This runs before, and independently of, any ratified cut, because it is a
    # CORRECTNESS check rather than a quality one: `frac_rise <= 0` means chi2 at a
    # bracket end is at or BELOW chi2 at the point the optimizer reported as the
    # minimum. A minimum that is not a minimum is not a measurement, whatever its value.
    #
    # There is no number to choose here. Zero is the boundary between "the objective rose
    # away from the answer" and "it did not", so the criterion cannot be tuned and has no
    # A_solar back-door — which is exactly why it survived the sweep that refuted every
    # threshold (RYA-847, 9 cells, 581 synthesis lines).
    #
    # `<= 0` and not `< 0`: a perfectly FLAT objective rises by exactly zero and is just
    # as unconstrained as one that dips. Measured on the sweep, no line ties at zero, so
    # the two agree on today's data — but the semantics are "chi2 did not rise", and that
    # is the form that stays correct when a tie eventually appears.
    fr = (metrics.get("frac_rise_weaker") if hasattr(metrics, "get")
          else getattr(metrics, "frac_rise_weaker", None))
    if fr is not None and np.isfinite(float(fr)) and float(fr) <= 0.0:
        return GateVerdict(
            False,
            f"NON-MINIMUM: chi2 at a bracket end is not above the reported minimum "
            f"(frac_rise = {float(fr):.3g} <= 0), so the optimizer returned a point that "
            f"is not a minimum and the fit did not determine A(X). A correctness "
            f"failure, not a quality threshold. Measured and retained; excluded from the "
            f"aggregate only (RYA-711).",
            ratified=True)

    # ── 2. the ratified cut, if one exists ────────────────────────────────────
    c = _cut()
    if not c:
        return GateVerdict(True, "", ratified=False)
    metric, maximum = c
    if metric not in _DIRECTION:
        raise ConstraintGateError(
            f"SYNTH_CONSTRAINT names metric {metric!r}, which this gate does not know how "
            f"to read. Known: {sorted(_DIRECTION)}. A cut on an unknown quantity must "
            f"fail loudly rather than pass everything.")
    v = metrics.get(metric) if hasattr(metrics, "get") else getattr(metrics, metric, None)
    v = None if v is None else float(v)

    # NOT MEASURABLE IS THE STRONGEST FAILURE, NOT A MISSING VALUE.
    # `sigma_A` is NaN exactly when the curvature could not be inverted: the fit is railed
    # against a bound, or the objective is flat. Treating that as "no data, let it
    # through" would admit precisely the fits this gate exists to catch — and RYA-848
    # showed what happens when an unmeasurable sigma is rendered as a number instead
    # (0.000 for a railed fit, the tightest possible bar for the least constrained fit).
    if v is None or not np.isfinite(v):
        return GateVerdict(
            False,
            f"UNCONSTRAINED FIT: {metric} is not measurable, so the data did not pin "
            f"A(X) at all — the fit is railed against a search bound or chi2 is flat in "
            f"A. Measured and retained; excluded from the aggregate only (RYA-711).",
            ratified=True)

    fails = v > maximum if _DIRECTION[metric] == "above" else v < maximum
    if fails:
        return GateVerdict(
            False,
            f"UNCONSTRAINED FIT (chi2 flat in A): {metric} = {v:.4g} "
            f"{_DIRECTION[metric]} the ratified {maximum}. The fit returned a number but "
            f"the objective barely moves with A(X), so it is not a measurement. Measured "
            f"and retained; excluded from the aggregate only (RYA-711).",
            ratified=True)
    return GateVerdict(True, "", ratified=True)
