"""Should this fit enter the aggregate? — RYA-847, the DECIDER.

Deliberately a different module from `pipeline.fit_constraint`, which MEASURES. That
separation is enforced by a test: `fit_constraint` may not contain a threshold at all. The
reason is RYA-161 — when the measurer and the cut live together, the cut gets chosen from
whatever product is in front of whoever is editing, and the metric that justified it is
never re-examined.

    fit_constraint.py    measures sigma_A / frac_rise / edge_distance. No verdict.
    constraint_gate.py   applies the RATIFIED cut, or refuses to invent one.

THE CUT IS NOT RATIFIED YET, AND THAT IS AN EXPLICIT STATE
----------------------------------------------------------
`SYNTH_CONSTRAINT` in `config/constants.py` is `None` until RYA-847's sweep across every
synthesis band chooses both the metric SHAPE and its value. While it is None this gate
carries the metrics onto every line and gates NOTHING, and `describe()` says so in the
product provenance.

That is the whole point of making it a state rather than a missing feature. The defect
RYA-843 found was not a wrong threshold, it was that nobody could see there was no
threshold: the band-product route reimplemented accept/reject as `status != 'ok'`, threw
`red_chi2` away, and quietly accepted fits whose chi2 moved 2.2% across eight dex of iron.
A gate that announces "no cut is ratified, N lines carry unexamined constraint metrics" is
honest. A gate that silently defaults to a plausible number is how the first one happened.

WHY THERE IS NO DEFAULT
-----------------------
The obvious default would be RYA-342's `SYNTH_CHI2_GATE = 10.0`, which is already applied
on the Engine-B path. Measured, it does not transfer: it was ratified against a bimodal
solar Fe II distribution (clean <= 3.0, blends >= 128, wide empty gap), while a
well-behaved NIR line sits at red_chi2 = 72 and all 40 near-UV lines sit between 27.7 and
999.5 on fits that are otherwise perfectly constrained. Inheriting it here would refuse
good lines in three bands to tidy one.
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
    if not c:
        return ("CONSTRAINT GATE: NOT RATIFIED. Every line carries sigma_A, "
                "frac_rise_weaker, edge_distance_dex and red_chi2, and NONE of them is "
                "gated on — RYA-847's sweep sets the metric and the cut across all "
                "synthesis bands, and RYA-161 forbids choosing either from one product.")
    metric, maximum = c
    return (f"CONSTRAINT GATE: {metric} {_DIRECTION[metric]} {maximum} excludes the line "
            f"from the aggregate (measured and retained, never dropped — RYA-711).")


def verdict(metrics) -> GateVerdict:
    """Apply the ratified cut to one line's metrics.

    `metrics` is a `ConstraintMetrics`, or any object carrying the named attributes, or a
    mapping — the band-product route flattens them onto a `LineMeasurement` and the
    handler keeps the dataclass, and neither should have to convert for the other.
    """
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
