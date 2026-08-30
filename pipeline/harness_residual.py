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

✅ `SynthesisHandler`'S 0.0 IS NOW MEASURED, NOT ASSUMED (RYA-875). It stood here as a
declared divergence — charged 0.0 while its banked control reported −0.0100 — until that
control was found to be comparing `median(fitted A)` over 18 lines against a HARDCODED
SCALAR that is the median of a DIFFERENT 23-line set. Paired against the banked engine's
own per-line answers the offset is **0.0000**, with 17 of 18 lines inside ±0.009 dex and a
MAD of 0.0060. The charged value never changed; what changed is that somebody measured it.
The divergence declaration was deleted BECAUSE THE NUMBERS AGREE, and
`tests/test_harness_residual_rya869.py` now enforces that undeclared handlers MATCH their
banked control — so the deletion is checked rather than merely permitted.

"""
from __future__ import annotations

from dataclasses import dataclass

from pipeline import treatment_axes  # RYA-1040: the ⟨3D⟩ pair is keyed off the
                                     # axis registry, never retyped here

#: How a charged harness residual came to be the number it is. The budget PRINTS this,
#: because the alternative was printing a claim: `harness_term` asserted "MEASURED
#: against the known optical answer, not assumed zero" beside whatever value it was
#: handed, including a zero nobody measured (RYA-873).
PROV_MEASURED = "measured"      # a control established this number
PROV_UNCHARGED = "uncharged"    # no residual is charged, and the reason is named
PROV_UNSTATED = "unstated"      # the caller did not say — printed as such, never assumed

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

#: Where each handler's control result is banked, so a charged number can be checked
#: against the thing that established it — RYA-875.
#:
#: 🔴 THIS EXISTS SO THAT DELETING A DIVERGENCE DECLARATION IS NOT A WAY TO STOP
#: CHECKING. `UNCHARGED_CONTROL_RESIDUAL_DEX` records handlers whose charged value
#: disagrees with their control, and the test that holds it only ever looked at handlers
#: IN that mapping — so removing an entry removed the check with it. That is precisely
#: the "relax the test to close the ticket" move RYA-873 forbade, available by accident.
#: With this map the other branch is real: a handler NOT declared as a divergence must
#: MATCH its artifact.
#:
#: A handler with no entry here has no banked control file to check against (the profile
#: fitter's 0.0129 is recorded in `pipeline/measure/profile_fit.py`'s module docstring,
#: not as JSON). That is a gap, not a pass, and it is named rather than hidden.
HANDLER_CONTROL_ARTIFACT: dict[str, str] = {
    "SynthesisHandler": "data/audit/synthesis_control/control_FeI.json",
}

#: Handlers whose CHARGED residual disagrees with their own banked control, each with the
#: control's value, the artifact, the ticket that owes the resolution, and why the number
#: is not simply charged. EMPTY, and that is a result rather than a default.
#:
#: 🔴 IT IS NOT A PLACE TO PARK A NUMBER NOBODY LIKES. `SynthesisHandler` lived here from
#: RYA-873 until RYA-875, and it left the only way an entry may: **the numbers agreed.**
#: Its control had reported −0.0100 by comparing a median over 18 lines to a scalar that
#: was the median of a different 23-line set; paired per line, the offset is 0.0000.
#: An entry is deleted when the disagreement is RESOLVED, never when the test that
#: surfaces it is relaxed — and `tests/test_harness_residual_rya869.py` now checks the
#: other branch too, so removing an entry re-arms a real assertion instead of dropping one.
UNCHARGED_CONTROL_RESIDUAL_DEX: dict[str, dict] = {}


@dataclass(frozen=True)
class HarnessResidual:
    """One handler's harness term, and the budget arguments that state it."""
    handler: str
    residual_dex: float
    #: RYA-873 — HOW this number was arrived at, carried so the budget prints the truth
    #: rather than a fixed claim. It travels with the value and the label for the same
    #: reason they travel with each other (RYA-869): they are one decision.
    provenance: str = PROV_MEASURED

    def budget_kwargs(self) -> dict:
        """Exactly the harness arguments `error_budget.build()` must be handed.

        Returned as a mapping rather than as two fields for the reason RYA-869 exists:
        the value and the LABEL are one decision, and the published defect was a call
        site that got them from the same `if` and still managed to attribute the number
        to the wrong handler in both. A caller cannot pass one and forget the other.
        """
        return {"harness_residual_dex": float(self.residual_dex),
                "handler": self.handler,
                "harness_provenance": self.provenance}

    def describe(self) -> str:
        note = UNCHARGED_CONTROL_RESIDUAL_DEX.get(self.handler)
        extra = ("" if note is None else
                 f" [⚠️ NOT CHARGED — its control is unresolved; abundance angle "
                 f"{note['control_dex']:.4f} dex, {note['control_artifact']}, "
                 f"{note.get('ticket', '')}]")
        return (f"{self.handler} harness residual {self.residual_dex:.4f} dex "
                f"({self.provenance}){extra}")


def for_handler(handler: str) -> HarnessResidual:
    """The residual `handler` earns. Unknown handler -> loud, never a silent 0.0."""
    name = str(handler or "").strip()
    if name not in HANDLER_RESIDUAL_DEX:
        raise KeyError(
            f"no measured harness residual is declared for handler {handler!r}; known "
            f"handlers are {sorted(HANDLER_RESIDUAL_DEX)}. A handler with no control has "
            f"no measured systematic and must not be charged zero (RYA-869) — control it "
            f"and declare the number in pipeline.harness_residual, do not default it.")
    return HarnessResidual(
        handler=name, residual_dex=HANDLER_RESIDUAL_DEX[name],
        # A handler listed as a known divergence has NOT had a residual established, and
        # the budget must not describe its number as measured (RYA-873).
        provenance=(PROV_UNCHARGED if name in UNCHARGED_CONTROL_RESIDUAL_DEX
                    else PROV_MEASURED))


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


def uncharged_note(handler: str) -> str:
    """One sentence for the budget: why this handler carries no residual. '' if it does."""
    note = UNCHARGED_CONTROL_RESIDUAL_DEX.get(str(handler).strip())
    if not note:
        return ""
    return (f"{note['why_not_charged']}. See "
            f"{note.get('ticket', 'the divergence declaration')} and "
            f"{note['control_artifact']}")


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
# EXHAUSTIVE over `band_products.TREATMENTS`, so a NEW treatment cannot be added
# without deciding which handler produces it,
# (RYA-1040: this said "a seventh" while there were six. A comment that names a
# count goes stale the moment the tuple grows -- the same defect shape as the
# spelling list it exists to replace, so it now states the invariant instead.) and it checks the table agrees with the
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
    # ... and this one adds a departure term to that same inversion (RYA-798)
    "ENGINE-A": "ProfileFitHandler",
    # ... while these re-fit the FLUX and share only the line set (RYA-712/784/798)
    "ENGINE-B": "SynthesisHandler",
    "ENGINE-B-NLTE": "SynthesisHandler",
    # 🔴 RYA-1106 — ENGINE-A-3DNLTE MOVED HERE FROM THE PROFILE-FIT BLOCK ABOVE. It was
    # listed as a profile-fit inversion on the strength of RYA-817's comment; RYA-1104
    # tested that against the data and refuted it -- the 1D-LTE base under this treatment
    # is the SYNTHESIS pool, 0e+00 dex on 50/50 lines, all four holdings. It re-fits the
    # flux like the ENGINE-B pair, so it takes their handler and their 0.0 residual. This
    # is what stops the un-earned +0.0129 dex being charged to the four Amarsi VIS cells.
    "ENGINE-A-3DNLTE": "SynthesisHandler",
    # never written under this route — the lab-gf pool is the near-UV synthesis route —
    # but stated rather than omitted, so the table is total and the test can say so.
    "1D-LTE-LABGF": "SynthesisHandler",
    # RYA-1040 — the ⟨3D⟩ pair. Both re-fit the FLUX on the ⟨3D⟩ STAGGER atmosphere, so
    # both are SynthesisHandler; the LTE member differs only in having the departures off,
    # which changes the physics and not the route. Keyed off the axis registry rather than
    # retyped, so this table cannot name a treatment the vocabulary does not have.
    treatment_axes.MEAN3D_NLTE_STAGGER.token: "SynthesisHandler",
    treatment_axes.MEAN3D_LTE_STAGGER.token: "SynthesisHandler",
    # RYA-1045 — the 1D rung's paired comparand. Same flux re-fit as `ENGINE-B-NLTE`
    # above and on the same MARCS.GES atmosphere; departures withheld is a physics
    # difference, not a route one, so it shares that handler.
    treatment_axes.GERBER1D_LTE_MARCS.token: "SynthesisHandler",
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
