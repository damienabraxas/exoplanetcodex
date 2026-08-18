"""Which harness residual does a product carry, and under whose name? — RYA-869.

The harness residual is the measurement handler's own systematic: what the handler
reproduces the known optical answer to, MEASURED rather than assumed zero
(`pipeline.measure.base.MeasurementHandler.systematic_dex`). It is a property of the
HANDLER THAT PRODUCED THE NUMBER and of nothing else -- not of the band, not of the
engine, and not of the treatment label the product happens to wear.

WHAT WENT WRONG, AND WHY A NAMED SET WOULD NOT HAVE BEEN ENOUGH
---------------------------------------------------------------
`scripts/derive_band_products.py` decided it like this::

    is_b = prod.treatment == "ENGINE-B"
    harness_residual_dex = 0.0 if is_b else PROFILE_FIT_RESIDUAL_DEX
    handler = "SynthesisHandler" if is_b else "ProfileFitHandler"

with a comment stating the rule correctly. The comment was right and the code was an
EQUALITY AGAINST A TREATMENT NAME, and that name gained a variant in RYA-798:
`ENGINE-B-NLTE` is the same `SynthesisHandler` flux fit on a different departure deck,
and `"ENGINE-B-NLTE" != "ENGINE-B"`. So every NLTE product was charged the profile
fitter's 0.0129 dex and LABELLED `ProfileFitHandler` in its own budget file, next to the
LTE product of the same handler charged 0.0000 and labelled `SynthesisHandler`. One
handler, two labels, two bars. Four published Fe matrix cells carried it (RYA-783/807):
VIS Fe I/II and red-optical Fe I/II, all ENGINE-B-NLTE, syst 0.1705 -> 0.1700 and
0.1731 -> 0.1726. The direction is that the bar was too LARGE, so nothing was ever going
to make it look wrong.

🔴 THE FIX IS NOT A BIGGER SET. Adding `"ENGINE-B-NLTE"` to the equality -- or writing
`treatment in {"ENGINE-B", "ENGINE-B-NLTE"}` -- rebuilds the same defect one variant
later, and `pipeline.band_products.TREATMENTS` has three more `X-VARIANT` names waiting
(`ENGINE-A-3DNLTE`, `1D-LTE-LABGF`). Worse, the mapping treatment -> handler DOES NOT
EXIST even today: the near-UV `1D-LTE` product is a flux fit (`derive_band_products`
`--synth-nearuv`, RYA-759) while the VIS `1D-LTE` product of the same treatment name is
an EW inversion of a profile fit. Any function from the label to the handler is wrong on
that pair before it meets a new variant.

So the product is asked what handler made it. `pipeline.band_products.Product` carries
`handler`, `build_product` refuses a product that does not declare one, and the residual
and the label are handed to `error_budget.build()` TOGETHER by `budget_kwargs()` -- the
RYA-855 shape, for the RYA-845 reason: a decision split across two arguments is a
decision that drifts between call sites, and the drift here is exactly what happened.

WHAT IS DECLARED HERE AND WHAT IS NOT
-------------------------------------
This module owns the number, ONCE. It was written out at `derive_band_products.py:77`
and again at `rya850_graded_products.py:73` -- two homes for one constant, the RYA-845
shape -- and both now import it.

⚠️ `SynthesisHandler` IS CHARGED 0.0 AND ITS BANKED CONTROL SAYS 0.0100 (RYA-870-adjacent
finding; see `UNCHARGED_CONTROL_RESIDUAL_DEX` below). That divergence is REPORTED, not
silently absorbed, and `tests/test_harness_residual_rya869.py` fails if it stops being
declared. It is deliberately NOT fixed here: this ticket moves four bars by 0.0005 dex
and folding a second bar-moving term into the same diff would make neither attributable
(RYA-848 -- prove a change with a same-inputs control).
"""
from __future__ import annotations

from dataclasses import dataclass

#: The MEASURED optical residual of each measurement handler, keyed by
#: `MeasurementHandler.name` (which is the class name).
#:
#: * `ProfileFitHandler` — 0.0129 dex. Reproduces the banked HARPS Fe I pool to a median
#:   EW ratio of 0.971 over 47 lines, MAD 0.060 (`pipeline/measure/profile_fit.py`
#:   module docstring, RYA-713/429). Charged as the harness systematic in every frontier
#:   band rather than assumed zero.
#: * `SynthesisHandler` — 0.0 dex AS PUBLISHED, which is not what its control measured.
#:   See `UNCHARGED_CONTROL_RESIDUAL_DEX`.
#:
#: 🔴 A HANDLER WITH NO ENTRY IS AN ERROR, NOT A ZERO. `for_handler` raises. A handler
#: that has not been controlled has no measured systematic and must not pretend to one —
#: that is `MeasurementHandler.systematic_dex`'s own rule, and defaulting to 0.0 here
#: would understate every frontier bar it touched while looking like a measurement.
HANDLER_RESIDUAL_DEX: dict[str, float] = {
    "ProfileFitHandler": 0.0129,
    "SynthesisHandler": 0.0,
}

#: Handlers whose CHARGED residual above disagrees with their own banked control, with
#: the control's value and the reason the charged number still stands. Declared so the
#: disagreement is a fact in the repo rather than an omission, and asserted in
#: `tests/test_harness_residual_rya869.py` so it cannot quietly stop being true.
#:
#: 🔴 `SynthesisHandler`'s control PASSED at dex_offset −0.0100
#: (`data/audit/synthesis_control/control_FeI.json`, n=18, tol 0.05, RYA-770/759) — and
#: every synthesis-route budget in the repo charges 0.0000 while `harness_term()` prints
#: the words "MEASURED against the known optical answer, not assumed zero". The prose and
#: the arithmetic disagree, which is the RYA-845 shape in a new place. Not fixed in
#: RYA-869: it moves the near-UV and every ENGINE-B cell, which is a different diff.
UNCHARGED_CONTROL_RESIDUAL_DEX: dict[str, dict] = {
    "SynthesisHandler": {
        "control_dex": 0.0100,
        "control_artifact": "data/audit/synthesis_control/control_FeI.json",
        "why_not_charged": (
            "published band products charge 0.0; charging the measured 0.0100 moves the "
            "near-UV 1D-LTE cell and every ENGINE-B / ENGINE-B-NLTE cell, so it is a "
            "separate diff from RYA-869's four-cell handler-attribution fix"),
    },
}


@dataclass(frozen=True)
class HarnessResidual:
    """One handler's harness term, and the pair of budget arguments that state it."""
    handler: str
    residual_dex: float

    def budget_kwargs(self) -> dict:
        """Exactly the harness arguments `error_budget.build()` must be handed.

        Returned as a mapping rather than as two fields for the reason RYA-869 exists:
        the value and the LABEL are one decision, and the published defect was a call
        site that got them from the same `if` and still managed to attribute the number
        to the wrong handler in both. A caller cannot pass one and forget the other.
        """
        return {"harness_residual_dex": float(self.residual_dex),
                "handler": self.handler}

    def describe(self) -> str:
        note = UNCHARGED_CONTROL_RESIDUAL_DEX.get(self.handler)
        extra = ("" if note is None else
                 f" [⚠️ its control measured {note['control_dex']:.4f} dex; "
                 f"{note['control_artifact']}]")
        return f"{self.handler} harness residual {self.residual_dex:.4f} dex{extra}"


def for_handler(handler: str) -> HarnessResidual:
    """The residual `handler` earns. Unknown handler -> loud, never a silent 0.0."""
    name = str(handler or "").strip()
    if name not in HANDLER_RESIDUAL_DEX:
        raise KeyError(
            f"no measured harness residual is declared for handler {handler!r}; known "
            f"handlers are {sorted(HANDLER_RESIDUAL_DEX)}. A handler with no control has "
            f"no measured systematic and must not be charged zero (RYA-869) — control it "
            f"and declare the number in pipeline.harness_residual, do not default it.")
    return HarnessResidual(handler=name, residual_dex=HANDLER_RESIDUAL_DEX[name])


def for_product(product) -> HarnessResidual:
    """The residual for one `pipeline.band_products.Product`, from ITS OWN handler.

    Reads `product.handler` and nothing else. It deliberately cannot see the treatment:
    a function of the label is the defect this module was written to remove, and the
    near-UV `1D-LTE` flux fit proves no such function exists even before a new variant
    is added.
    """
    handler = getattr(product, "handler", "")
    if not handler:
        raise ValueError(
            f"product {getattr(product, 'treatment', '?')!r} declares no handler, so the "
            f"harness residual it earns is unknown. Deriving it from the treatment name "
            f"is what RYA-869 removed — the producing route knows its handler and must "
            f"say so at build_product().")
    return for_handler(handler)


# ── archaeology: cells written before `Product.handler` existed ───────────────────────
#
# Every band product written from RYA-869 on carries its handler, in `Product.handler`
# and in the `handler` column of its `*_products.csv`. The artifacts BANKED BEFORE that
# do not, and three consumers still read them: `scripts/rya855_rung_audit.py` (per-line
# files), `scripts/rya850_graded_products.py` (the RYA-783 matrix) and this ticket's own
# before/after diff.
#
# 🔴 THIS IS THE ONLY PLACE IN THE REPO WHERE A HANDLER IS INFERRED, IT IS INFERRED FROM
# THE ROUTE AND NOT FROM THE ENGINE, AND IT IS FOR BANKED FILES ONLY. Nothing on a live
# path may call it: the producing route knows its handler and says so at `build_product`.
# `tests/test_harness_residual_rya869.py` holds it to that — it checks the table is
# EXHAUSTIVE over `band_products.TREATMENTS`, so a seventh treatment cannot be added
# without deciding which handler produces it, and it checks the table agrees with the
# handler the live deriver declares for every cell where both answers exist.
#
# The two keys are needed because neither alone decides it:
#
#   * the ROUTE is not enough — an `ENGINE-B` cell is written into a `*_PROFILEFIT_*`
#     stem, because that stem names where the LINE SET came from, not who measured it;
#   * the TREATMENT is not enough — `1D-LTE` is a profile fit in VIS and a flux fit in
#     the near-UV, which is the pair that makes a treatment->handler function impossible.

#: Routes that are a flux fit END TO END, whatever treatment name the cell wears:
#: `SYNTH` is RYA-759's near-UV synthesis route and `LABGF` is RYA-836's sub-pool of it.
_BANKED_FLUX_FIT_ROUTES = frozenset({"SYNTH", "LABGF"})

#: Within the EW-sourced `PROFILEFIT` route, which treatment re-measured the spectrum.
#: Exhaustive over `pipeline.band_products.TREATMENTS` by test, so a new treatment has to
#: state its handler here rather than inherit one by falling off the end of a set.
_BANKED_PROFILEFIT_HANDLER: dict[str, str] = {
    # inverted from a profile-fit EW ...
    "1D-LTE": "ProfileFitHandler",
    # ... and these two add a departure term to that same inversion (RYA-798/817)
    "ENGINE-A": "ProfileFitHandler",
    "ENGINE-A-3DNLTE": "ProfileFitHandler",
    # ... while these re-fit the FLUX and share only the line set (RYA-712/784/798)
    "ENGINE-B": "SynthesisHandler",
    "ENGINE-B-NLTE": "SynthesisHandler",
    # never written under this route — the lab-gf pool is the near-UV synthesis route —
    # but stated rather than omitted, so the table is total and the test can say so.
    "1D-LTE-LABGF": "SynthesisHandler",
}


def handler_of_banked_cell(*, route: str, treatment: str) -> str:
    """The handler that produced a band-product cell banked before RYA-869.

    ⚠️ ARCHAEOLOGY ONLY. For anything the deriver writes now, read `Product.handler` or
    the artifact's `handler` column; inferring it is what RYA-869 removed.
    """
    r, t = str(route).strip().upper(), str(treatment).strip()
    if r in _BANKED_FLUX_FIT_ROUTES:
        return "SynthesisHandler"
    if r != "PROFILEFIT":
        raise KeyError(
            f"unknown band-product route {route!r}; known routes are "
            f"{sorted(_BANKED_FLUX_FIT_ROUTES | {'PROFILEFIT'})}")
    if t not in _BANKED_PROFILEFIT_HANDLER:
        raise KeyError(
            f"no handler recorded for treatment {t!r} on the {r} route. A treatment "
            f"added to pipeline.band_products.TREATMENTS must state which handler "
            f"produces it (RYA-869) — inheriting one by default is the defect.")
    return _BANKED_PROFILEFIT_HANDLER[t]


def charged_in_banked_cell(*, route: str, treatment: str) -> HarnessResidual:
    """What a band-product cell banked BEFORE RYA-869 was actually CHARGED and LABELLED.

    ⚠️ THIS IS THE DEFECT, QUOTED. It exists so that a consumer reading a pre-RYA-869
    artifact can tell the difference between *"this cell does not reproduce because its
    inputs drifted"* and *"this cell does not reproduce because RYA-869 corrected its
    harness term"*. Those need opposite handling: the first invalidates a diff, the
    second is the diff.

    🔴 DO NOT CALL IT FROM ANYTHING THAT PRODUCES A NUMBER, and do not "fix" it to agree
    with `handler_of_banked_cell` — the day the two agree, every banked artifact has been
    re-derived and both this function and its caller can go.

    The rule, verbatim from `scripts/derive_band_products.py` as of `10016d5`::

        is_b = prod.treatment == "ENGINE-B"
        harness_residual_dex = 0.0 if is_b else PROFILE_FIT_RESIDUAL_DEX
        handler = "SynthesisHandler" if is_b else "ProfileFitHandler"

    with the synthesis route (`SYNTH`/`LABGF`) hardcoding its own 0.0 / SynthesisHandler
    pair at a separate call site — which is why it was right there and wrong here.
    """
    if str(route).strip().upper() in _BANKED_FLUX_FIT_ROUTES:
        return HarnessResidual("SynthesisHandler", 0.0)
    is_b = str(treatment).strip() == "ENGINE-B"          # ← the defect
    return HarnessResidual(
        "SynthesisHandler" if is_b else "ProfileFitHandler",
        0.0 if is_b else HANDLER_RESIDUAL_DEX["ProfileFitHandler"])
