"""The star's parameter systematics, as a budget term — RYA-1120 / RYA-282 §2.

🔴 THE DEFECT THIS CLOSES. There were TWO uncertainty budgets in this repo and nothing
joined them. `pipeline/uncertainty_stack.py` perturbs Teff, log g, xi and [Fe/H] by the
star's own 1-sigma uncertainties, re-derives A(X), and writes sigma_B per element to
`data/audit/uncertainty/`. `pipeline/error_budget.build()` assembles what a product
publishes. A `git grep` for that artifact returned only its own two stamping scripts:
the term was derived, sourced, and never travelled. RYA-282 §2 made it MANDATORY in
June 2026; every live VIS Fe product still published a `sigma_syst` with no Teff, no
log g and no microturbulence in it (RYA-1112 F1). **This module is that join.**

WHAT IT REFUSES TO DO
---------------------
🔴 It will not spread one element-level number across every product. RYA-1093 measured
xi-sensitivity as a STRONG-LINE phenomenon, so dA/dxi is a property of the LINE SET:
a DEEPGRADED pool (feature depth > 0.60, saturated by construction) and a GRADED pool
(<= 0.60) do not share one, and RYA-1089's -0.24 dex/(km/s) was measured on ONE 62-line
pool. That number is a reference MAGNITUDE for the size of the hole, never a per-product
value. A pool with no measured response of its own gets an UNMEASURED term — which makes
the published bar an honest floor — not a borrowed one.

🔴 And it will not decide applicability from an engine name. `"3D" in treatment` swept
the <3D> MEAN in with full 3D and exempted both (RYA-1092's name-vs-physics pattern).
They differ by measurement: RYA-1099 ran the <3D> mean at xi=0 and got +0.137 dex WORSE,
because a mean atmosphere averages the velocity structure OUT and the route still runs on
an inherited xi. So the caller states applicability with a reason, and `full_3d` is
matched against declared route tokens rather than a substring.

THE ARITHMETIC (RYA-282 §3)
---------------------------
    sigma_params = sqrt( sum_p (|dA/dp| * delta_p)^2 )   for p in Teff, logg, vmic, [Fe/H]

with `delta_p` the STAR's own parameter uncertainty, read from `uncertainty_stack` rather
than re-typed. For the Sun three of the four are ~0 by construction (Teff known to ~1 K,
log g fixed by the IAU nominal mass and radius, [Fe/H] = 0 by definition), so the solar
budget is xi-dominated — and RYA-1093 set that allowance to the unresolved method+selection
spread |1.0 - 0.709| rather than the (much smaller) formal error, deliberately, because
the formal one would have passed the gate (RYA-161). The value itself is deliberately not
spelled anywhere in this file: it has exactly one home, in `uncertainty_stack`, and a
number written down twice drifts (RYA-845).

It does NOT average down: perturbing Teff moves every line the same way, so more lines
cannot shrink it. That is the whole reason it must be published beside sigma/sqrt(N).
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

#: Route tokens whose model resolves the velocity field that microturbulence stands in
#: for. FULL 3D only — the <3D> MEAN is NOT here (RYA-1099 measured it as xi-sensitive),
#: and membership is by exact token, never by substring.
FULL_3D_TREATMENTS = frozenset({"ENGINE-A-3DNLTE"})

#: Recorded WITH its contrary measurement rather than asserted bare. `vmic` IS an input
#: axis of the Amarsi MLP (training box 0-3 km/s) and its CORRECTION responds at
#: d(aberr)/dxi = +0.0985 dex/(km/s) measured on its own optical Fe I set (112 of 117
#: lines, none exactly zero). The exemption rests on the published value being
#: A(1D-LTE)(xi) + aberr(xi), whose two halves move oppositely — the correction tracking
#: xi is the baseline's xi-dependence being undone, not a second one added.
FULL_3D_XI_EXEMPTION = (
    "NOT APPLICABLE — full 3D resolves the velocity field xi stands in for (Ryan, "
    "RYA-1120 2026-08-29). Recorded, not silent: the Amarsi MLP does take vmic as an "
    "input axis and its correction measures d(aberr)/dxi = +0.0985 dex/(km/s); the "
    "exemption rests on the published value being A(1D-LTE)(xi) + aberr(xi), whose two "
    "halves move oppositely.")


class UnknownParameter(KeyError):
    """A response was supplied for a parameter the star has no declared delta_p for."""


@dataclass(frozen=True)
class StellarParamSystematic:
    """One product's stellar-parameter term, and the budget arguments that state it.

    `responses` maps a parameter name (as `uncertainty_stack` spells it: 'teff_K',
    'logg', 'vturb_kms', 'feh') to |dA/dp| MEASURED ON THIS PRODUCT'S OWN POOL. A
    parameter absent from the mapping is UNMEASURED for this product; a parameter whose
    delta_p is 0 contributes exactly 0 whatever its response, and needs no measurement.
    """
    star_id: str
    deltas: dict           # delta_p, from uncertainty_stack -- never re-typed here
    responses: dict = field(default_factory=dict)   # |dA/dp| measured on THIS pool
    applicable: bool = True
    applicability_note: str = ""
    pool_note: str = ""

    def __post_init__(self):
        unknown = set(self.responses) - set(self.deltas)
        if unknown:
            raise UnknownParameter(
                f"{sorted(unknown)} has no declared delta_p for {self.star_id!r}. A "
                f"response cannot be charged against an uncertainty nobody declared "
                f"(RYA-282 §2: an undeclared delta_p is a STOP, not a zero).")
        if not self.applicable and self.responses:
            raise ValueError(
                "a NOT APPLICABLE stellar-parameter term cannot also carry measured "
                "responses -- if the parameters move this product, they apply to it.")

    # ── which parameters actually need a measurement ────────────────────────
    def required(self) -> list[str]:
        """Parameters with a non-zero delta_p. Only these can contribute."""
        return sorted(p for p, d in self.deltas.items() if d and float(d) > 0)

    def missing(self) -> list[str]:
        """Required parameters this product has NOT measured. Non-empty => floor."""
        return [p for p in self.required() if p not in self.responses]

    def sigma_dex(self) -> float | None:
        """sqrt(sum (|dA/dp| * delta_p)^2), or None if any required response is missing.

        None rather than a partial sum: a budget that quietly added the two terms it
        happened to have would publish a total that is smaller than the truth and look
        exactly like a complete one.
        """
        if not self.applicable or self.missing():
            return None
        return math.sqrt(sum((abs(float(self.responses[p])) * float(self.deltas[p])) ** 2
                             for p in self.required()))

    def contributions(self) -> dict:
        """Per-parameter dex, for the breakdown artifact (RYA-282 §3)."""
        out = {}
        for p, d in sorted(self.deltas.items()):
            dp = float(d or 0.0)
            if dp <= 0:
                out[p] = {"delta_p": dp, "response": None, "dex": 0.0,
                          "basis": "delta_p = 0 -> contributes exactly 0, no run needed"}
            elif p in self.responses:
                r = abs(float(self.responses[p]))
                out[p] = {"delta_p": dp, "response": r, "dex": r * dp,
                          "basis": "measured on this product's own pool"}
            else:
                out[p] = {"delta_p": dp, "response": None, "dex": None,
                          "basis": "UNMEASURED on this pool -- no perturb-and-re-derive"}
        return out

    def source(self) -> str:
        if not self.applicable:
            return self.applicability_note or FULL_3D_XI_EXEMPTION
        miss = self.missing()
        if miss:
            return (f"NOT MEASURED for this product: no perturb-and-re-derive on this "
                    f"pool for {', '.join(miss)} (RYA-282 §2). dA/dp is a property of "
                    f"the LINE SET (RYA-1093), so no element-level number may stand in "
                    f"for it. {self.pool_note}").strip()
        parts = ", ".join(
            f"|dA/d{p}|={abs(float(self.responses[p])):.4g} x delta={float(self.deltas[p]):.4g}"
            for p in self.required())
        return (f"RYA-282 §2 perturb-and-re-derive on this product's own pool "
                f"[{parts}]; delta_p from uncertainty_stack for {self.star_id}. "
                f"{self.pool_note}").strip()

    def budget_kwargs(self) -> dict:
        """Exactly the stellar-parameter arguments `error_budget.build()` must be handed.

        A mapping rather than three fields, for the reason `harness_residual` returns one:
        the value, its source and its applicability are ONE decision, and a caller must
        not be able to pass the number while forgetting to say where it came from.
        """
        return {"stellar_param_sigma_dex": self.sigma_dex(),
                "stellar_param_source": self.source(),
                "stellar_param_applicable": bool(self.applicable)}


def for_product(treatment: str, *, star_id: str = "solar",
                responses: dict | None = None, pool_note: str = "") \
        -> StellarParamSystematic:
    """The stellar-parameter term a product earns.

    `treatment` decides applicability against `FULL_3D_TREATMENTS` — by exact token.
    Everything else applies, including the <3D> mean.
    """
    from pipeline.uncertainty_stack import params_and_deltas
    _, deltas = params_and_deltas(star_id)
    if str(treatment) in FULL_3D_TREATMENTS:
        return StellarParamSystematic(star_id=star_id, deltas=deltas, responses={},
                                      applicable=False,
                                      applicability_note=FULL_3D_XI_EXEMPTION,
                                      pool_note=pool_note)
    return StellarParamSystematic(star_id=star_id, deltas=deltas,
                                  responses=dict(responses or {}), pool_note=pool_note)
