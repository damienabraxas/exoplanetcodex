"""RYA-906 — THE ONE PLACE THAT KNOWS WHAT A TREATMENT LABEL MEANS.

Engine naming caused three real defects in one cycle, one of them merged and published:

  * RYA-869   `prod.treatment == "ENGINE-B"` never matched `ENGINE-B-NLTE`, so four
              published Fe matrix bars carried the wrong handler's systematic. A string
              compare against a label whose variant set had silently grown.
  * RYA-896   the element tracker says the Fe II arbiter is "synthesis ... B x3" while the
              artifact's arbiter is ENGINE-A / ProfileFitHandler — an EW product. The
              label and the artifact disagree about which engine it even is.
  * RYA-906   keying on the deck directory (redundant with the label) inflated an audit
              to 53 cells against a real 31.

ROOT CAUSE: one token carries five orthogonal axes — route x scale x model x atmosphere x
gf-pool — collapsed inconsistently. `1D-LTE` names physics. `ENGINE-A` names a letter. The
same route differing only in scale sits on two different naming schemes. Because the axes
are collapsed into one string, code keys on the wrong axis and a new variant slips past an
`==`. "Everything NLTE" today means maintaining a list of spellings — which IS the RYA-869
bug class.

So: STORE THE AXES, DERIVE THE NAME. Code keys on fields. Nothing parses a display string,
and nothing types one.

    route   ew | synth                      HOW the abundance was extracted
    scale   1D-LTE | 1D-NLTE | 3D-NLTE      the physics of the model atmosphere/departures
    model   none | bergemann | amarsi | gerber    the NLTE model FAMILY
    atmos   atlas9 | marcs-ges              STORED, never derived (see below)
    gf      kurucz | lab                    the gf pool

🔴 ROUTE IS NEVER READ FROM THE LABEL. `1D-LTE` is used by BOTH routes: the VIS cells are
EW inversions, and the near-UV Fe I cell (n=40) is a RYA-759 synthesis flux fit wearing the
identical label. Measured on the committed tree: `1D-LTE` pairs with ProfileFitHandler on
some rows and SynthesisHandler on others, and so does `1D-LTE-LABGF`. **On those labels the
legacy string is not merely lossy, it is FALSE**, and any mapping that reads route off it
gets the near-UV cell backwards. That cell is the acceptance canary for this whole module.

🔴 `atmos` IS STORED, NOT DERIVED, even though it is currently predictable from `model`
(Gerber rides MARCS.GES; the EW/Bergemann/Amarsi routes ride ATLAS9). Deriving it would be
correct today and silently wrong the first day a model runs on a different atmosphere —
which is the entire lesson of the defects above. Ryan ratified storing it independently.

WHAT THIS MODULE DOES NOT DO: it does not rewrite the legacy `treatment` column. That
column is permanent (Ryan: dual-label forever). Rewriting it would make every historical
product incomparable to its own past and produce a diff that looks like values moved when
nothing did (the RYA-874 lesson). Axis columns are added ALONGSIDE it; every diff is a
pure addition.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict

# ── the axis vocabularies ────────────────────────────────────────────────────
ROUTES = ("ew", "synth")
SCALES = ("1D-LTE", "1D-NLTE", "3D-NLTE")
MODELS = ("none", "bergemann", "amarsi", "gerber")
ATMOSPHERES = ("atlas9", "marcs-ges")
GF_POOLS = ("kurucz", "lab")

#: The route each HANDLER implements. This is the authoritative witness (RYA-869 put the
#: handler on the product precisely so the question could be answered from data).
ROUTE_BY_HANDLER = {
    "ProfileFitHandler": "ew",
    "SynthesisHandler": "synth",
}

#: What each legacy label pins DOWN — scale, model and gf pool only. Route is deliberately
#: ABSENT from this table for the `1D-LTE` family: see the module docstring. For the
#: ENGINE-* family the route IS determined by the label, because those labels were coined
#: to name a route and pair with exactly one handler across the whole committed tree —
#: `_route_from_label` holds that, and `test_treatment_axes_rya906` re-measures it rather
#: than trusting this comment.
LEGACY = {
    "1D-LTE":          dict(scale="1D-LTE",  model="none",      gf="kurucz"),
    "1D-LTE-LABGF":    dict(scale="1D-LTE",  model="none",      gf="lab"),
    "ENGINE-A":        dict(scale="1D-NLTE", model="bergemann", gf="kurucz"),
    "ENGINE-A-3DNLTE": dict(scale="3D-NLTE", model="amarsi",    gf="kurucz"),
    "ENGINE-B":        dict(scale="1D-LTE",  model="none",      gf="kurucz"),
    "ENGINE-B-NLTE":   dict(scale="1D-NLTE", model="gerber",    gf="kurucz"),
}

#: Route for the labels that name one. `None` means "this label does not say" — which is
#: the honest answer for the `1D-LTE` family and the whole reason this module exists.
_ROUTE_BY_LABEL = {
    # 🔴 RYA-1002 — ENGINE-A NO LONGER PINS A ROUTE. It was EW-only until the synthesis
    # route gained its own ENGINE-A block, and it now appears on both. `None` is the
    # honest answer, the same one `1D-LTE` has always carried: the label does not say
    # which route ran, so the route must come from the product's HANDLER (RYA-869/906).
    # Leaving "ew" here would have made every synthesis-route NLTE product render as an
    # EW measurement — the exact mislabel this module exists to prevent.
    "ENGINE-A": None, "ENGINE-A-3DNLTE": "ew",
    "ENGINE-B": "synth", "ENGINE-B-NLTE": "synth",
    "1D-LTE": None, "1D-LTE-LABGF": None,
}

#: The atmosphere each model family actually ran on, used ONLY to fill legacy rows that
#: predate the stored column. Live emitters pass `atmos` explicitly. Kept separate from
#: LEGACY so that "we inferred this" and "the product said this" never look alike.
_LEGACY_ATMOS = {"gerber": "marcs-ges"}
_LEGACY_ATMOS_DEFAULT = "atlas9"

DISPLAY_MODEL = {"bergemann": "Bergemann", "amarsi": "Amarsi", "gerber": "Gerber"}
DISPLAY_ROUTE = {"ew": "EW", "synth": "Synth"}
#: 🔴 There is deliberately NO entry for an unknown route. A display name is the NAME OF A
#: PRODUCT, not a status report about it: a placeholder like "route?" is not the name of
#: anything, it reads as a rendering bug to anyone who sees it, and it duplicates a fact
#: that already has a proper home in `route_basis`. When the route is unknown the segment
#: is OMITTED and the name states exactly what is known and no more.


class UnknownTreatment(ValueError):
    """A treatment label this module has never heard of.

    Raised rather than defaulted. A new engine variant must be ADDED here, which is the
    single edit RYA-869 needed and did not get: a silent default is how `ENGINE-B-NLTE`
    slipped past an `== "ENGINE-B"` and reached publication.
    """


@dataclass(frozen=True)
class Axes:
    """The stored key. Code keys on THESE, never on a label or a display name."""
    route: str | None
    scale: str
    model: str
    atmos: str
    gf: str
    #: WHICH witness settled `route`. Carried so the derivation is auditable instead of
    #: silent: an inference and a reading must never look alike in an artifact.
    route_basis: str = "unknown"

    @property
    def is_nlte(self) -> bool:
        """'Everything NLTE' as a FIELD TEST, replacing a list of label spellings.

        This one property is the point of the refactor: `scale != "1D-LTE"` cannot be
        outgrown by a new variant, whereas `treatment in ("ENGINE-A", "ENGINE-B-NLTE")`
        silently stops being true the day someone adds one.
        """
        return self.scale != "1D-LTE"

    @property
    def display(self) -> str:
        """`route · scale · model` (+ `· lab-gf`). DERIVED — never typed, never parsed.

        An UNKNOWN route contributes NOTHING: the name becomes `1D-LTE` rather than
        `route? · 1D-LTE`. That is deliberately indistinguishable from the legacy label,
        because for those rows the legacy label is exactly how much is known — and
        `route_basis` says `unknown` beside it, which is where the caveat belongs. A name
        carrying a question mark is a rendering bug wearing a product's clothes.

        Atmosphere is stored but summarised OUT: it is currently redundant with `model`,
        so it clutters the label without adding a decision. It lives in the product detail
        and the download, for provenance and for the day that correlation breaks.
        """
        parts = []
        if self.route in DISPLAY_ROUTE:
            parts.append(DISPLAY_ROUTE[self.route])
        parts.append(self.scale)
        if self.model != "none":
            parts.append(DISPLAY_MODEL.get(self.model, self.model))
        if self.gf == "lab":
            parts.append("lab-gf")
        return " · ".join(parts)

    def as_columns(self) -> dict:
        """The axis columns, for emitters adding them ALONGSIDE `treatment`."""
        return asdict(self)


def _route_from_label(treatment: str) -> str | None:
    return _ROUTE_BY_LABEL.get(treatment)


def resolve_route(treatment: str, *, handler=None, ew_inversion=None) -> tuple[str | None, str]:
    """(route, basis) from the strongest witness available. Never from the label alone
    where the label is known to lie.

    Precedence, strongest first:
      1. `handler`       — RYA-869 put it on the product to answer exactly this.
      2. `ew_inversion`  — a per-line bool: False means the abundance came from a flux
                           fit, not an EW inversion. It is what proves the near-UV Fe I
                           cell is synthesis, on rows written before `handler` existed.
      3. the label       — ONLY for the ENGINE-* family, which was coined to name a route.
                           Never for the `1D-LTE` family.
      4. unknown         — stated, not guessed.
    """
    if handler:
        route = ROUTE_BY_HANDLER.get(str(handler).strip())
        if route:
            return route, "handler"
    if ew_inversion is not None and str(ew_inversion).strip() != "":
        s = str(ew_inversion).strip().lower()
        if s in ("true", "1"):
            return "ew", "ew_inversion"
        if s in ("false", "0"):
            return "synth", "ew_inversion"
    by_label = _route_from_label(treatment)
    if by_label:
        return by_label, "label-family"
    return None, "unknown"


def axes_for(treatment: str, *, handler=None, ew_inversion=None,
             model=None, atmos=None, gf=None) -> Axes:
    """Legacy label (+ whatever route evidence the product carries) -> stored axes.

    `atmos` should be PASSED by any live emitter that knows it. It is only inferred for
    rows that predate the stored column, and `route_basis` is not the place that records
    the difference — the inference is deliberately narrow and documented at `_LEGACY_ATMOS`.
    """
    key = (treatment or "").strip()
    if key not in LEGACY:
        raise UnknownTreatment(
            f"unknown treatment {treatment!r}. Add it to LEGACY and _ROUTE_BY_LABEL in "
            f"pipeline/treatment_axes.py — deliberately a hard failure, because a silent "
            f"default is how RYA-869 published four wrong systematics.")
    spec = LEGACY[key]
    route, basis = resolve_route(key, handler=handler, ew_inversion=ew_inversion)
    model = model or spec["model"]
    if model not in MODELS:
        raise ValueError(f"unknown model axis {model!r}; expected one of {MODELS}")
    return Axes(
        route=route, scale=spec["scale"], model=model,
        atmos=(atmos or _LEGACY_ATMOS.get(model, _LEGACY_ATMOS_DEFAULT)),
        gf=(gf or spec["gf"]), route_basis=basis)


def display_for(treatment: str, **kw) -> str:
    """Convenience: the derived display name for a legacy label plus its route evidence."""
    return axes_for(treatment, **kw).display
