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

#: 🔴 RYA-1040 — `<3D>` AND `3D` ARE DIFFERENT PHYSICS AND MUST NEVER COLLAPSE.
#: `3D-NLTE` is the Amarsi+2022 FULL 3D result: the radiative transfer is solved through
#: the inhomogeneous cube. `<3D>-NLTE` is a mean-3D result: the cube is averaged on
#: surfaces of constant tau500 FIRST, so the horizontal inhomogeneity — the term that
#: distinguishes 3D from 1D at all — is gone before the transfer runs. A ⟨3D⟩ number is
#: closer to 1D than to 3D in what it models, and reporting one under `3D-NLTE` would
#: claim the expensive result while having run the cheap one. The angle-bracket notation
#: is the field's own (Nordlander & Lind 2017 call theirs ⟨3D⟩ for exactly this reason).
SCALES = ("1D-LTE", "1D-NLTE", "<3D>-LTE", "<3D>-NLTE", "3D-NLTE")

#: Which scales are LTE. Stated as DATA rather than inferred from a substring, so that a
#: scale whose spelling does not contain a recognisable term fails loudly instead of being
#: silently classified. `is_nlte` reads this.
LTE_SCALES = ("1D-LTE", "<3D>-LTE")
NLTE_SCALES = ("1D-NLTE", "<3D>-NLTE", "3D-NLTE")

MODELS = ("none", "bergemann", "amarsi", "gerber")

#: `stagger-mean3d` is the ⟨3D⟩ STAGGER model (Magic et al. 2013) averaged on constant
#: tau500 surfaces -- `data/atmospheres/stagger_avg3d_rya442/sun_avg3d_stagger.mod`,
#: 101 depths, log tau500 -5..+5. It is a DIFFERENT atmosphere from `marcs-ges`, which is
#: why the axis is stored rather than derived from `model`: Gerber rides MARCS.GES on its
#: 1D decks and this one on its ⟨3D⟩ deck, so the correlation the legacy inference relies
#: on breaks here -- exactly the day the module's docstring said it would.
ATMOSPHERES = ("atlas9", "marcs-ges", "stagger-mean3d")

GF_POOLS = ("kurucz", "lab")

#: 🔴 RYA-1040 — THE DECK-PROVENANCE AXIS: WHOSE DEPARTURES THESE ARE.
#: Ratified by Ryan 2026-08-25. `stagger` is the vendor's STAGGERmean3D deck (Gerber et
#: al. 2023, fetched from the MPG Keeper share); `codex` is a deck we solve ourselves,
#: later. A `…gerber.codex` product is DISTINCT from a `…gerber.stagger` one, stands
#: beside it, and is NEVER blended with it (RYA-282/RYA-712).
#:
#: It is a separate axis from `model` because they answer different questions: `model` is
#: the NLTE model FAMILY (whose atom and rates), `deck` is who ran the statistical-
#: equilibrium solve and on which grid. Our own deck would still be `gerber`-family and
#: would still not be the vendor's numbers.
#:
#: ⚠️ It is also the fix for the deck-key bug: `for_node` resolves the deck from the
#: ELEMENT, so it picks `Fe` and can never reach `Fe@mean3D`. An explicit provenance axis
#: lets a call site SELECT the deck instead of guessing it from the element.
DECKS = ("none", "stagger", "codex")

#: The route each HANDLER implements. This is the authoritative witness (RYA-869 put the
#: handler on the product precisely so the question could be answered from data).
ROUTE_BY_HANDLER = {
    "ProfileFitHandler": "ew",
    "SynthesisHandler": "synth",
}

#: 🔴 RYA-1100 — THE ROUTE TOKEN A PUBLISHED PRODUCT RECORD CARRIES.
#: `scripts/publish_product.py` reads a products CSV, not a Product object, so it has the
#: stored `route` field and no `handler`. Without a witness it could resolve route ONLY
#: from the label -- and that is how 22 live ENGINE-A products, every one of them on
#: route=SYNTH, were published displaying "EW · 1D-NLTE · Bergemann". RYA-1002 removed the
#: route from `_ROUTE_BY_LABEL["ENGINE-A"]` for exactly this reason and the publisher's
#: own hand-typed map kept asserting it anyway.
#:
#: Stated as DATA, and unknown tokens are NOT defaulted: an unrecognised route token means
#: the route is unknown, which `resolve_route` already knows how to say.
ROUTE_BY_PRODUCT_TOKEN = {
    "PROFILEFIT": "ew",
    "SYNTH":      "synth",
    "EW-3D":      "ew",
}

#: 🔴 RYA-1100 — LABELS THAT NAME NOTHING NEW. A deprecated alias is a legacy spelling
#: whose axes are IDENTICAL to another label's; it survives only so historical products
#: stay readable (Ryan: dual-label forever -- rewriting the stored `treatment` column
#: would make every past product incomparable to its own past, the RYA-874 lesson).
#:
#: `ENGINE-B` is `1D-LTE` measured by the synthesis handler. Compare LEGACY: both declare
#: scale=1D-LTE, model=none, gf=kurucz. The ONLY thing the token added was a pinned route,
#: and RYA-906 moved route onto the handler -- which left the label naming nothing at all.
#: Measured on the live feed: 13 of 15 ENGINE-B products were byte-identical to a 1D-LTE
#: product in the same cell. One measurement, published twice under two names.
#:
#: The alias is NOT given a distinct display string. There is nothing to distinguish, and
#: inventing a name for it would re-create the collision one level down.
DEPRECATED_ALIASES = {"ENGINE-B": ("1D-LTE", "synth")}

#: 🔴 RYA-1102 — THE ERROR-PLOT ROW AXIS. Ordered, fixed, and DERIVED, never typed in JS.
#:
#: The plot used to emit ONE ROW PER PRODUCT, so a section's row count was whatever the
#: data happened to contain (Kitt Peak rendered 19). Ryan's rule is the opposite: every
#: section shows the SAME rows in the SAME order, and a model with no product for that
#: instrument/band says N/A. A missing measurement is information; a missing ROW is not.
#:
#: Ordered as a physics ladder -- 1D-LTE, then 1D-NLTE, then <3D>, then full 3D -- with
#: each scale's EW leg before its synthesis leg. EW and synthesis are SEPARATE PRODUCTS
#: (a profile fit and a flux fit are different measurements of different line pools), so
#: they are two rows, never one row to collapse.
#:
#: Each entry is (treatment, route_token, pending). The NAME is not stored: it is derived
#: through `display_for`, so this axis cannot drift from the labels the products carry.
#: Publishing it into the feed is what makes a NEW MODEL appear on the site with no site
#: change -- the "products should update live" requirement.
#:
#: ⚠️ `pending` rows are DECLARED BUT NOT RENDERED. The codex <3D>-NLTE deck is the model
#: we solve ourselves (`DECKS`); it is on the axis so the ladder is complete and so
#: turning it on is one flag, not a schema change. A rendered row for a model that does
#: not exist yet would read as a measurement we failed to make.
PLOT_ROW_AXIS = (
    ("1D-LTE",                           "PROFILEFIT", False),
    ("1D-LTE",                           "SYNTH",      False),
    ("synth-1D-LTE-gerber",              "SYNTH",      False),
    ("ENGINE-A",                         "PROFILEFIT", False),
    ("ENGINE-A",                         "SYNTH",      False),
    ("ENGINE-B-NLTE",                    "SYNTH",      False),
    ("synth-mean3D-LTE-gerber-stagger",  "SYNTH",      False),
    ("synth-mean3D-NLTE-gerber-stagger", "SYNTH",      False),
    ("ENGINE-A-3DNLTE",                  "EW-3D",      False),
    ("synth-mean3D-NLTE-gerber-codex",   "SYNTH",      True),
)


def plot_row_axis(*, include_pending: bool = False) -> list[dict]:
    """The ordered row axis as {name, treatment, route, pending} -- names DERIVED.

    A row whose treatment is registered but unnameable RAISES rather than being skipped:
    an axis that silently drops a row is how a section stops showing a model without
    anyone noticing, which is the failure this axis exists to prevent.
    """
    out = []
    for treatment, route, pending in PLOT_ROW_AXIS:
        if pending and not include_pending:
            continue
        if pending and treatment not in AXIS_NATIVE and treatment not in LEGACY:
            # A declared-but-unbuilt model has no axes yet, so it carries its intent
            # rather than a derived name, and never renders.
            out.append({"name": None, "treatment": treatment, "route": route,
                        "pending": True})
            continue
        out.append({"name": display_for(treatment, route_token=route),
                    "treatment": treatment, "route": route, "pending": pending})
    return out


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


def _nlte_from_scale(scale: str) -> bool:
    """Is this scale non-LTE? From the declared vocabularies, never from a substring."""
    if scale in LTE_SCALES:
        return False
    if scale in NLTE_SCALES:
        return True
    raise ValueError(
        f"unknown scale {scale!r}: it is in neither LTE_SCALES nor NLTE_SCALES, so "
        f"whether it is non-LTE is UNDECIDED. Add it to one of them in "
        f"pipeline/treatment_axes.py. Deliberately a hard failure -- the previous "
        f"implementation defaulted anything that was not '1D-LTE' to NLTE, which is how "
        f"an LTE comparand would have been counted as an NLTE product (RYA-1040).")


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
    #: RYA-1040 deck provenance. Defaults to `none` so every pre-existing construction and
    #: every committed row keeps its exact meaning: a product with no departure deck has
    #: no deck provenance, which is a fact and not a gap. New column, pure addition.
    deck: str = "none"

    @property
    def is_nlte(self) -> bool:
        """'Everything NLTE' as a FIELD TEST, replacing a list of label spellings.

        This one property is the point of the refactor: `treatment in ("ENGINE-A",
        "ENGINE-B-NLTE")` silently stops being true the day someone adds a variant.

        🔴 RYA-1040 — AND SO DID THE FIRST FIELD TEST. It was `scale != "1D-LTE"`, which
        is a spelling list of length one wearing a field test's clothes: it assumed the
        only LTE scale would ever be the 1D one. Adding `<3D>-LTE` — the comparand
        RYA-1040 REQUIRES, so that the NLTE effect is measured on ONE atmosphere — made
        an LTE product answer True. The module's own docstring claimed this "cannot be
        outgrown by a new variant"; it was outgrown by the second one, in the direction
        nobody was watching.

        Now it reads the scale's own LTE/NLTE term, which is what the axis actually
        encodes, and an unrecognised scale RAISES rather than defaulting to NLTE —
        because defaulting is how the original defect published four wrong systematics.
        """
        return _nlte_from_scale(self.scale)

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
        # RYA-1040: deck provenance is shown only when there IS a deck. `none` is not a
        # provenance, and a trailing "· none" would read as a missing value rather than
        # as the absence of the thing.
        if self.deck != "none":
            parts.append(self.deck)
        if self.gf == "lab":
            parts.append("lab-gf")
        return " · ".join(parts)

    @property
    def token(self) -> str:
        """The FILESYSTEM-SAFE product key. Derived from the same axes as `display`.

        🔴 RYA-1040 — WHY THIS IS NOT JUST `display`. The treatment string is not only a
        label: `derive_band_products` writes `f"{stem}_{prod.treatment}_lines.csv"`, so it
        becomes part of a FILENAME. Ryan's ratified science-facing label is
        `synth . <3D>-NLTE . gerber . stagger`, and `<` / `>` are shell redirection
        characters — a file called `…_<3D>-NLTE_lines.csv` needs quoting in every command
        that ever touches it, and silently truncates in any that forgets. The dots read as
        extension boundaries besides.

        So the LABEL keeps Ryan's form exactly and the KEY is a slug of the same axes.
        Both are DERIVED; neither is typed. `<3D>` becomes `mean3D`, which is the same
        distinction in characters a filename can hold.

        ⚠️ The two must never be independently edited. `test_treatment_axes_rya1040`
        asserts they carry the same axes, so a change to one that is not a change to the
        other is a test failure rather than a pair of names that drift apart.
        """
        parts = [self.route or "route-unknown", self.scale.replace("<3D>", "mean3D"),
                 self.model]
        if self.deck != "none":
            parts.append(self.deck)
        if self.gf == "lab":
            parts.append("lab-gf")
        return "-".join(parts)

    def as_columns(self) -> dict:
        """The axis columns, for emitters adding them ALONGSIDE `treatment`."""
        return asdict(self)


#: 🔴 RYA-1040 — AXIS-NATIVE TREATMENTS. NEW PRODUCTS ARE DEFINED AS AXES, NOT AS LABELS.
#:
#: `LEGACY` above maps a hand-coined label BACK to axes, which is the migration path for
#: products that already existed. Anything new goes here instead: the axes are the
#: definition and both names are derived from them, which is what RYA-906 asked for and
#: what `LEGACY` cannot deliver because its keys are the very letters RYA-906 retired.
#:
#: Ryan, 2026-08-25, ratifying the name: *"DROP the label ENGINE-B-MEAN3D -- it revives
#: the ENGINE-A/B vocabulary RYA-906 retired. The science-facing label follows
#: route.scale.model.deck-provenance."*
#:
#: 🔴 THE TWO MEMBERS ARE A MANDATORY PAIR, AND THAT IS PHYSICS, NOT TIDINESS. The NLTE
#: effect must be (⟨3D⟩-NLTE minus ⟨3D⟩-LTE) on ONE atmosphere. Differencing ⟨3D⟩-NLTE
#: against 1D-LTE instead folds the 1D→mean-3D ATMOSPHERE shift into the reported NLTE
#: effect -- the confound RYA-542 had to disentangle for Ti, where MARCS at +0.203
#: reproduced the deck's +0.221 and the difference turned out to be the atmosphere. Both
#: run the same deck; the LTE member simply has the departures OFF.
#:
#: ⚠️ They are two SEPARATE PRODUCTS, never a correction applied to one (RYA-712).
AXIS_NATIVE: dict[str, "Axes"] = {}


def _register_axis_native(**kw) -> "Axes":
    ax = Axes(**kw)
    if ax.scale not in SCALES:
        raise ValueError(f"unknown scale axis {ax.scale!r}; expected one of {SCALES}")
    if ax.model not in MODELS:
        raise ValueError(f"unknown model axis {ax.model!r}; expected one of {MODELS}")
    if ax.atmos not in ATMOSPHERES:
        raise ValueError(f"unknown atmos axis {ax.atmos!r}; expected one of {ATMOSPHERES}")
    if ax.deck not in DECKS:
        raise ValueError(f"unknown deck axis {ax.deck!r}; expected one of {DECKS}")
    AXIS_NATIVE[ax.token] = ax
    return ax


#: `synth · <3D>-NLTE · Gerber · stagger` — the Gerber solver's departures, solved on the
#: vendor STAGGERmean3D deck, applied through the ⟨3D⟩ STAGGER atmosphere they belong to.
MEAN3D_NLTE_STAGGER = _register_axis_native(
    route="synth", scale="<3D>-NLTE", model="gerber", atmos="stagger-mean3d",
    gf="kurucz", route_basis="axis-native", deck="stagger")

#: `synth · <3D>-LTE · Gerber · stagger` — THE SAME DECK WITH THE DEPARTURES OFF, and the
#: mandatory comparand for the member above. `model` stays `gerber` and `deck` stays
#: `stagger` because this product's identity is "the LTE limit OF THAT DECK'S SETUP", run
#: on the same atmosphere and line list; calling it `model="none"` would make it
#: indistinguishable from an unrelated LTE run and destroy the pairing.
MEAN3D_LTE_STAGGER = _register_axis_native(
    route="synth", scale="<3D>-LTE", model="gerber", atmos="stagger-mean3d",
    gf="kurucz", route_basis="axis-native", deck="stagger")

#: 🔴 THE 1D RUNG'S MISSING COMPARAND — RYA-1045.
#: The ⟨3D⟩ pair above exists so the NLTE effect is measured on ONE atmosphere. The 1D
#: rung never got the same treatment. `ENGINE-B-NLTE` runs the Gerber 1D deck on MARCS.GES
#: (the deck is computed there, RYA-798); the only LTE product on that route runs on
#: ATLAS9.Castelli. Differencing them reports the ATLAS9→MARCS.GES ATMOSPHERE change as
#: non-LTE physics — the RYA-542 confound, on the rung nobody checked. Measured on solar
#: Fe I: that mixed delta is +0.050, and it is NOT the 1D NLTE effect.
#:
#: ⚠️ ONLY THE COMPARAND IS NEW. The NLTE member already exists as the legacy label
#: `ENGINE-B-NLTE`, whose axes are exactly (synth, 1D-NLTE, gerber, marcs-ges, deck=none).
#: Registering a second axis-native token for the same axes would give one product two
#: names, so this pairs against the legacy token instead.
#:
#: `deck="none"` because the deck axis records ⟨3D⟩ SOLVE PROVENANCE (stagger vs codex);
#: the 1D MARCS deck is not a ⟨3D⟩ deck, and calling it `stagger` would claim it came from
#: the STAGGERmean3D solve. `model` stays `gerber` for the same reason the ⟨3D⟩ comparand
#: does: its identity is "the LTE limit OF THAT DECK'S SETUP".
GERBER1D_LTE_MARCS = _register_axis_native(
    route="synth", scale="1D-LTE", model="gerber", atmos="marcs-ges",
    gf="kurucz", route_basis="axis-native", deck="none")

NLTE_LTE_PAIRS = {MEAN3D_NLTE_STAGGER.token: MEAN3D_LTE_STAGGER.token,
                  # the legacy NLTE label, paired to its new axis-native comparand
                  "ENGINE-B-NLTE": GERBER1D_LTE_MARCS.token}


def comparand_for(token: str) -> str | None:
    """The LTE product a given NLTE product must be differenced against, or None.

    🔴 The NLTE effect is a DIFFERENCE OF TWO PRODUCTS ON ONE ATMOSPHERE, and this is the
    function that says which two. Reporting (⟨3D⟩-NLTE − 1D-LTE) as an NLTE effect
    attributes the atmosphere change to non-LTE physics (RYA-542).
    """
    return NLTE_LTE_PAIRS.get(token)


def _route_from_label(treatment: str) -> str | None:
    return _ROUTE_BY_LABEL.get(treatment)


def resolve_route(treatment: str, *, handler=None, ew_inversion=None,
                  route_token=None) -> tuple[str | None, str]:
    """(route, basis) from the strongest witness available. Never from the label alone
    where the label is known to lie.

    Precedence, strongest first:
      1. `handler`       — RYA-869 put it on the product to answer exactly this.
      2. `route_token`   — the stored `route` field of a published product record
                           (RYA-1100). Same provenance class as `handler`: written by the
                           route that actually ran. Ranked below it only because `handler`
                           names the class that did the work while this names its output.
      3. `ew_inversion`  — a per-line bool: False means the abundance came from a flux
                           fit, not an EW inversion. It is what proves the near-UV Fe I
                           cell is synthesis, on rows written before `handler` existed.
      4. the label       — ONLY for the ENGINE-* family, which was coined to name a route.
                           Never for the `1D-LTE` family.
      5. unknown         — stated, not guessed.
    """
    if handler:
        route = ROUTE_BY_HANDLER.get(str(handler).strip())
        if route:
            return route, "handler"
    if route_token:
        route = ROUTE_BY_PRODUCT_TOKEN.get(str(route_token).strip().upper())
        if route:
            return route, "route_token"
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
             model=None, atmos=None, gf=None, route_token=None) -> Axes:
    """Legacy label (+ whatever route evidence the product carries) -> stored axes.

    `atmos` should be PASSED by any live emitter that knows it. It is only inferred for
    rows that predate the stored column, and `route_basis` is not the place that records
    the difference — the inference is deliberately narrow and documented at `_LEGACY_ATMOS`.
    """
    key = (treatment or "").strip()
    # RYA-1040: an axis-native product already IS its axes -- there is nothing to infer,
    # and no legacy label to map. Checked first so a new product never depends on LEGACY.
    if key in AXIS_NATIVE:
        return AXIS_NATIVE[key]
    if key not in LEGACY:
        raise UnknownTreatment(
            f"unknown treatment {treatment!r}. Add it to LEGACY and _ROUTE_BY_LABEL in "
            f"pipeline/treatment_axes.py — or, if it is NEW, define it in AXIS_NATIVE as "
            f"axes rather than coining another label (RYA-906/RYA-1040). Deliberately a "
            f"hard failure, because a silent default is how RYA-869 published four wrong "
            f"systematics.")
    spec = LEGACY[key]
    route, basis = resolve_route(key, handler=handler, ew_inversion=ew_inversion,
                                 route_token=route_token)
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
