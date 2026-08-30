"""
pipeline/product_eligibility.py — RYA-1092
==========================================
THE GATE that decides whether a product may sit in the LIVE `products[]` list of an
element feed. One place, called by the sweep and by CI.

🔴 WHY THIS EXISTS: THE POOLS WERE BUILT, THE GATE WAS NOT.

`data/products/<star>/<El>.json` has had `quarantine[]`, `superseded[]` and `archive[]`
since RYA-711, and the discipline that governs them was ratified long ago. But nothing
ever *checked* `products[]` against that discipline, so ineligible products stayed live
and the public feed rendered them. A pool without a gate is a filing cabinet, not a
policy — and the failure mode is the RYA-1026 one: **wrong information that does not
announce itself**. An incomplete product renders a number and a bar and looks exactly
like a complete one.

WHAT THE GATE IS NOT
--------------------
⚠️ It is NOT a judgement about physics quality, and it is NOT allowed to become one.
Being 1D, or LTE, or having a large error bar, is not a disqualification: a 1D-LTE
product on a corrected holding with a real systematic budget is one of the SOLID
products and stays live. Every criterion below is a check on whether the record is
COMPLETE and COMPARABLE, never on whether its answer is liked.

⚠️ IT NEVER DELETES. Failing the gate moves a product to `quarantine[]` with its reason
code, its `quarantined_at`, and every published field byte-identical. Quarantine is loud
and reversible; deletion destroys the evidence for its own rejection (RYA-711).

THE CRITERIA
------------
1. `UNCORRECTED_HOLDING` — the holding is not telluric/continuum corrected. Read through
   `telluric_display_policy` / `telluric_policy.applied_state`, which is the ONE source
   for that fact; a hand-written clean-set is the RYA-845 defect shape and a first draft
   of that module got two of five entries wrong by trying.
2. `PRE_CONTINUUM_FIX` — the ARTIFACT was produced before the last of the RYA-933/1030
   continuum fixes landed. Same holding, different code: a product measured before the
   fix is not comparable with the grid measured after it, which is precisely why RYA-933
   withdrew two of its own.
3. `SYST_INCOMPLETE` — no real `sigma_syst`. A product without a systematic budget has no
   error bar, and rendering `A ± sigma_stat` for it silently states a precision nobody
   computed (RYA-968: admission must not manufacture precision).
4. `NOT_YET_DEFENSIBLE` — the species is outside the ratified real set {Fe I, Fe II, Al}.
   Ryan's ruling: everything else is a pilot, and a pilot presented as a product is a
   claim we cannot defend.
5. `ANOMALOUS_SCATTER` — the line-to-line scatter is grossly inconsistent with the other
   pools measuring the SAME species, band and holding. See the section below; the
   statistic is chosen so that it is not merely "this product has few lines".
6. `STAT_BASIS_MISMATCH` — `sigma_stat` is a different STATISTIC from the one its peers
   publish, and the feed carries no field saying so. The bar renders beside theirs and is
   read as comparable when it is not.
7. `SUPERSEDED` — the record itself declares a replacement. Derived from the record, never
   inferred from an engine's name (see the Amarsi note).

🔴 THE AMARSI CORRECTION (Ryan, 2026-08-28) — AND WHY IT IS A DESIGN CONSTRAINT.
RYA-1092 was written believing the EW · 3D-NLTE · Amarsi products were SUPERSEDED by the
redone mean-⟨3D⟩ STAGGER engine. They are not. Amarsi 3D-NLTE is the authoritative
reference the STAGGER deck is validated AGAINST (the ~44% amplitude-deficit check); the
redo was a PREREQUISITE for building STAGGER, not a replacement for Amarsi. The two
coexist and Amarsi is a keeper.

What is actually wrong with the Amarsi products in the feed is that they are INCOMPLETE:
`sigma_syst` is null and the line-to-line scatter is ~1.3 dex against ~0.18 for the 1D
pools on the same holdings. So they fail criteria 3 and 5 and are held out on those
reasons — never on supersession.

The design consequence is that **`SUPERSEDED` must never be reachable from an engine
name.** It is read from an explicit `superseded_by` / `superseded_reason` on the record.
Had the gate been allowed to encode "Amarsi is superseded", the correction above would
have required editing the gate rather than editing a record — and the wrong verdict would
have been baked into the enforcement mechanism, which is the worst place for it.

🔴 `sigma_stat` DOES NOT MEAN THE SAME THING IN EVERY PRODUCT, AND THAT IS ITS OWN DEFECT
------------------------------------------------------------------------------------------
The feed's `sigma_stat` is written by two different routes that compute two different
statistics, and nothing in the record says which:

    route              what `sigma_stat` IS                       written at
    SYNTH, PROFILEFIT  a STANDARD ERROR (scatter / sqrt(n), then  error_budget.py:609
                       RMS'd with the other averaging-down terms)
    EW-3D              the RAW per-line STANDARD DEVIATION        band_products.py:506
                                                                  (np.std(vals, ddof=1))

MEASURED, not inferred. Ten band products reproduce their published `sigma_stat` exactly
as `std/sqrt(n)` from their own per-line CSVs (e.g. Fe I VIS GRADED molecfit: std 0.1780,
se 0.0218, published 0.0218). The committed RYA-817 EW-3D artifact publishes sigma 0.3418,
which is exactly the per-line std of its 114 lines -- its standard error is 0.0320.

⚠️ THIS IS WHY THE 4 AMARSI PRODUCTS *LOOKED* LIKE 7x OUTLIERS AND ARE NOT. Treating their
`sigma_stat` as a standard error and multiplying by sqrt(n) inflated a raw scatter of 0.18
into a fictitious 1.3 dex. Their real line-to-line scatter is 0.1818-0.1861, which sits
right on top of the 1D pools on the same holdings (0.16-0.21). The ticket's framing, and
this module's first draft, both made that inference; the data does not support it.

The consequence for eligibility is not that they become fine. It is a DIFFERENT and better
reason: the feed renders `A +/- sigma_stat` uniformly, so a product whose `sigma_stat` is a
raw scatter sits beside products whose `sigma_stat` is a standard error and is read as
having an eight-times-worse bar when it does not. That is `STAT_BASIS_MISMATCH` -- an
uninterpretable published uncertainty, and exactly the RYA-1026 shape of wrong information
that does not announce itself. The band artifacts already carry a `stat_basis` column
(`error_budget.stat_basis`) saying where the number came from; the feed drops it. Carrying
it is the real fix and it is not this ticket's.

WHY `ANOMALOUS_SCATTER` IS STILL MEASURED ON RAW SCATTER
-------------------------------------------------------
`sigma_stat` carries a 1/sqrt(n) that has nothing to do with data quality: a clean 3-line
pool has a large standard error for a reason that is not a defect. `raw_scatter` recovers
the line-to-line scatter -- the quantity the ticket describes and the only one comparable
across pools of different size -- by dispatching on the route's declared basis rather than
assuming one. An UNKNOWN route abstains: a product whose statistic cannot be identified is
reported, never compared (RYA-907).

⚠️ AFTER THE BASIS CORRECTION, NO LIVE PRODUCT IS ANOMALOUS. The criterion is implemented
because the ticket specifies it and because it will bite the next time a genuinely broken
pool appears -- but it currently has NO measured instance, and saying so is the honest
report. `ANOMALY_RATIO = 3.0` is therefore justified against a null ALONE (the measured
maximum ratio across the 27 comparable products is 1.84, so the constant sits 1.6x above
the worst observed good pool) and NOT against a separation from bad ones, because there are no
bad ones to separate from. `test_product_eligibility_rya1092` asserts that null and fails
if any product ever approaches the threshold, which forces a re-derivation rather than
letting a stale constant decide (RYA-161, [[a tolerance needs a measured null]]).
"""
from __future__ import annotations

import math
import statistics
from dataclasses import dataclass

#: The canonical product identity. Imported concept, not a second definition -- the same
#: nine fields `scripts/publish_product.KEY_FIELDS` uses. Duplicated as a literal here
#: only because importing a script from a pipeline module would invert the dependency;
#: `test_product_eligibility_rya1092` asserts the two lists are equal.
#:
#: 🔴 RYA-1127 ADDS `line_set`, AND IT CLOSES A LATENT COLLISION RATHER THAN PREVENTING A
#: FUTURE ONE. RYA-1106 measured the Amarsi 3D-NLTE method on Asplund's own AGSS21 line
#: set across the four VIS holdings. Those products pass the eligibility gate, but their
#: key was IDENTICAL to the our-graded Amarsi products': same element, ion, band,
#: instrument, holding, tier, selector, route and treatment -- everything except the pool
#: of lines they were measured on, which the key did not carry.
#:
#: ⚠️ AND THE COLLISION WAS BEING MASKED BY A MISLABEL. Before RYA-1106, the our-graded
#: Amarsi leg stored `route=EW-3D` -- the stranded `ProfileFitHandler`'s route, refuted
#: line by line by RYA-1104 -- while the replication is `SYNTH`. So `route` was
#: accidentally supplying the distinguishing axis, and the two products did not collide
#: only because one of them was wrong. Correcting the label removed the accidental
#: distinguisher and exposed the real defect: line_set was never in the identity.
#:
#: That is the same shape as the lying departure column RYA-1106 also fixed -- a wrong
#: value quietly doing a job nobody assigned it, so that fixing it looks like it caused
#: the breakage it revealed.
KEY_FIELDS = ("element", "ion", "band", "instrument", "holding",
              "tier", "selector", "route", "treatment", "line_set")

#: The species a product may claim. Ryan's ruling: {Fe I, Fe II, Al}. `None` for the ion
#: means any ionisation stage of that element.
REAL_SPECIES: frozenset = frozenset({("Fe", "I"), ("Fe", "II"), ("Al", None)})

#: 🔴 THE CONTINUUM FIXES, BY COMMIT. The cutoff below is DERIVED from these, not typed
#: from a memory of when they landed -- `test_product_eligibility_rya1092` re-resolves
#: each sha with `git show` and re-checks the cutoff whenever the objects are reachable.
#: Listing the shas rather than only the instant is what makes the number auditable: a
#: reader can see which fixes it is the cutoff FOR.
CONTINUUM_FIX_COMMITS: dict[str, str] = {
    "fafe9246bcdc528425207166e49902dc19a24e49":
        "RYA-933: the corrected Kitt Peak holdings could not measure a band",
    "79ec5a3484763aaef20ad8b8ba0227e7cf1c48ea":
        "RYA-933/1026 (merge #371): Kurucz 2005 ships its own continuum -- we were "
        "fitting our own",
    "a61e7b8009ca2e58a295e780883c77c8f10408fe":
        "RYA-1030 (merge #373): normalisation is DETERMINED from the flux at intake, "
        "never asserted from a label",
}

#: The LAST of the above, in UTC. A product measured before this instant was produced by
#: code missing at least one continuum fix.
#: ⚠️ IT IS THE LATEST, NOT THE EARLIEST, AND THAT CHOICE MOVES PRODUCTS. Measured on the
#: committed feed: with the earliest fix (fafe924, 00:57:30Z) as the cutoff NO live
#: product is pre-fix; with this one, EIGHT are. The ticket says "post RYA-933/1030"
#: -- both -- so the cutoff is where both are in, and the conservative direction is also
#: the reversible one (quarantine is loud and undoable; shipping is not).
CONTINUUM_FIX_CUTOFF_UTC = "2026-08-24T03:51:12Z"

#: Flag a product whose line-to-line scatter exceeds this multiple of the leave-one-out
#: median of its peer group. NOT fitted -- see the module docstring and the null test.
ANOMALY_RATIO = 3.0

#: Below this many peers a "median of the others" is not a population, so the check
#: abstains rather than letting one product declare the other anomalous.
ANOMALY_MIN_PEERS = 3

#: The basis the feed as a whole publishes, and therefore the one a product must share to
#: be readable beside the others. 59 of 63 live products are on the band route.
_MAJORITY_STAT_BASIS = "standard_error"


class EligibilityError(RuntimeError):
    """A product could not be evaluated -- refuse, never pass by default."""


@dataclass(frozen=True)
class Ineligible:
    """One reason a product may not be live. Products carry a LIST: a record can be
    incomplete AND anomalous, and collapsing that to a single code would let the fix for
    one hide the other."""
    code: str
    detail: str


def line_set_of(product: dict) -> str:
    """This product's `line_set`, RESOLVED -- never read straight off the record.

    🔴 `product.get("line_set")` IS THE WRONG QUESTION AND WOULD HAVE SILENTLY BROKEN THE
    KEY. Our own products do not store the field: they state their pool in `tier`, and
    RYA-1111 derives `our-graded` / `our-deep-graded` from it precisely so the value is
    not written twice and free to disagree. A plain `.get` therefore returns "" for all 66
    live products, which would collapse every one of them onto an empty axis -- the key
    would gain a column that is blank everywhere, look like it was working, and still let
    the RYA-1106 replication collide with the our-graded leg.

    Resolving instead means an unrecognised pool RAISES here (RYA-869: never default). The
    import is deferred because this module is deliberately dependency-light and
    `reference_lineset` pulls in numpy/pandas; the module cache makes it free after the
    first call.
    """
    from pipeline.reference_lineset import line_set_for_product
    return line_set_for_product(product)


def key_of(product: dict) -> str:
    return "|".join(line_set_of(product) if k == "line_set"
                    else str(product.get(k) or "") for k in KEY_FIELDS)


def species_of(product: dict) -> tuple:
    return (str(product.get("element") or "").strip(),
            str(product.get("ion") or "").strip() or None)


def is_real_species(product: dict) -> bool:
    el, ion = species_of(product)
    return (el, ion) in REAL_SPECIES or (el, None) in REAL_SPECIES


#: What `sigma_stat` IS, per route. 🔴 NOT A CONVENTION -- A MEASURED FACT about two code
#: paths that disagree, each cited to the line that writes it. Anything not listed here
#: ABSTAINS rather than being assumed into one basis or the other.
STAT_BASIS_BY_ROUTE: dict[str, str] = {
    "SYNTH":      "standard_error",   # error_budget.py:609, scatter / sqrt(n)
    "PROFILEFIT": "standard_error",   # same budget assembly
    "EW-3D":      "line_scatter",     # band_products.py:506, np.std(vals, ddof=1)
}


def stat_basis_of(product: dict) -> str | None:
    """`standard_error` / `line_scatter` / None if it cannot be established.

    🔴 THE RECORD'S OWN `stat_basis` WINS, AND THE ROUTE MAP IS THE LEGACY FALLBACK.
    `publish_product` now carries the producing budget's `stat_basis` onto every product
    (RYA-1095), so what `sigma_stat` means is a property of the RECORD. Inferring it from
    the route was always a lookup that goes stale the moment a route changes what it
    writes -- and that happened: the EW-3D route now emits `stat_dex` (a standard error)
    where it used to emit a raw `sigma`. Keying on the route would have kept calling the
    corrected products mismatched.

    The map stays for records published BEFORE the field existed. Those are historical and
    their basis really is a property of the route that wrote them.
    """
    declared = str(product.get("stat_basis") or "")
    if declared:
        # `error_budget.stat_basis()` is prose, not a token -- it says "measured — RMS of
        # the random terms, ... dex at n_lines=N" or names the quantiser floor. Both are
        # STANDARD-ERROR constructions: every random term in that RMS has already been
        # divided by sqrt(n) (`error_budget.py:609`).
        if "RMS of the random terms" in declared or "quantiser-floor" in declared:
            return "standard_error"
        return None            # an unrecognised declaration is not silently classified
    return STAT_BASIS_BY_ROUTE.get(str(product.get("route") or ""))


def raw_scatter(product: dict) -> float | None:
    """Line-to-line scatter, in dex -- the size-independent quantity.

    🔴 THE CONVERSION DEPENDS ON THE ROUTE AND ASSUMING OTHERWISE IS THE BUG THIS FUNCTION
    WAS BORN WITH. A first version multiplied every product's `sigma_stat` by sqrt(n),
    which is right for the band route and wrong for EW-3D, where `sigma_stat` is already
    the scatter -- inflating 0.18 into a fictitious 1.3 dex and making four sound pools
    look like 7x outliers. Returns None when the basis is unknown: a statistic that cannot
    be identified must not be compared (RYA-907).
    """
    s, n = product.get("sigma_stat"), product.get("n_lines")
    if s is None or not n or float(n) <= 0:
        return None
    basis = stat_basis_of(product)
    if basis == "line_scatter":
        return float(s)
    if basis == "standard_error":
        return float(s) * math.sqrt(float(n))
    return None


def peer_group_of(product: dict) -> tuple:
    """Products that measure the SAME lines by different methods: same species, band,
    holding and TIER. Treatment and route deliberately vary -- comparing a pool only
    against itself would make the check vacuous.

    🔴 TIER IS IN THE KEY BECAUSE TIERS ARE DIFFERENT LINE POPULATIONS, NOT DIFFERENT
    QUALITIES. `GRADED` is the lab-gf lines at or below the depth gate; `DEEPGRADED` is
    the saturated ones above it (RYA-984). Their line-to-line scatters differ for reasons
    of physics, so comparing across them measures the tier split rather than a defect.
    Measured on the committed feed: leaving tier OUT put a perfectly sound HARPS
    DEEPGRADED pool at 2.79x its group median -- within touching distance of
    ANOMALY_RATIO -- purely because the surviving GRADED pools in that band are small and
    unsaturated. With tier in, the whole null tops out at 1.84x on 27 comparable products.
    """
    return (product.get("element"), product.get("ion"),
            product.get("band"), product.get("holding"), product.get("tier"))


def _artifact_mtime(product: dict) -> str | None:
    return (product.get("provenance") or {}).get("artifact_mtime")


def evaluate(product: dict, *, peers: list | None = None) -> tuple:
    """Every reason this product may not be live. Empty tuple == eligible.

    `peers` are the other products in its peer group (excluding itself); pass them to
    enable the scatter check, omit to skip it.
    """
    from pipeline import telluric_display_policy as tdp

    out: list = []

    holding = str(product.get("holding") or "")
    if not holding:
        raise EligibilityError(
            f"product {key_of(product)!r} has no holding, so its telluric state cannot "
            f"be read. A product whose input is unknown must not be judged eligible by "
            f"default (RYA-907).")
    #: The two states RYA-1026 lets render AS SCIENCE. CONTROL_ONLY renders too, but as a
    #: labelled control -- and a control in `products[]` is a science claim by placement,
    #: which is exactly the confusion RYA-1026 exists to stop. BLOCKED/UNREGISTERED are
    #: refusals there and here.
    state = tdp.display_state(holding)
    if state not in ("CLEAN", "CLEAN_WITH_ANOMALY"):
        out.append(Ineligible(
            "UNCORRECTED_HOLDING",
            f"holding {holding!r} is display_state={state!r} -- not a telluric/continuum-"
            f"corrected science basis (RYA-1026, read through "
            f"telluric_policy.applied_state, never inferred from the holding name)."))

    mt = _artifact_mtime(product)
    if mt is None:
        out.append(Ineligible(
            "PRE_CONTINUUM_FIX",
            "the record carries no provenance.artifact_mtime, so it cannot be shown to "
            "postdate the continuum fixes. Unknown is not 'after' (RYA-907)."))
    elif mt < CONTINUUM_FIX_CUTOFF_UTC:
        out.append(Ineligible(
            "PRE_CONTINUUM_FIX",
            f"artifact produced {mt}, before the last of the RYA-933/1030 continuum "
            f"fixes landed ({CONTINUUM_FIX_CUTOFF_UTC}). Not comparable with the grid "
            f"measured after them."))

    if product.get("sigma_syst") is None:
        out.append(Ineligible(
            "SYST_INCOMPLETE",
            "sigma_syst is null -- the product has no systematic budget, so it has no "
            "error bar to publish."))

    basis = stat_basis_of(product)
    if basis is None:
        out.append(Ineligible(
            "STAT_BASIS_MISMATCH",
            f"route {product.get('route')!r} is not in STAT_BASIS_BY_ROUTE, so what its "
            f"sigma_stat MEANS is unrecorded. An uncertainty whose definition is unknown "
            f"must not render beside ones whose is."))
    elif basis != _MAJORITY_STAT_BASIS:
        out.append(Ineligible(
            "STAT_BASIS_MISMATCH",
            f"sigma_stat here is a {basis!r} ({product.get('route')} route, "
            f"band_products.py:506) while the rest of the feed publishes a "
            f"{_MAJORITY_STAT_BASIS!r} (error_budget.py:609). The feed carries no "
            f"stat_basis field, so `A +/- sigma_stat` renders as comparable when it is "
            f"not -- this product's bar reads ~sqrt(n) times worse than it is."))

    if not is_real_species(product):
        el, ion = species_of(product)
        out.append(Ineligible(
            "NOT_YET_DEFENSIBLE",
            f"{el} {ion} is outside the ratified real set {{Fe I, Fe II, Al}}. It is a "
            f"pilot measurement, not a product."))

    if product.get("superseded_by") or product.get("superseded_reason"):
        out.append(Ineligible(
            "SUPERSEDED",
            f"the record declares a replacement: "
            f"{product.get('superseded_by') or product.get('superseded_reason')}"))

    if peers:
        r = scatter_ratio(product, peers)
        if r is not None and r > ANOMALY_RATIO:
            out.append(Ineligible(
                "ANOMALOUS_SCATTER",
                f"line-to-line scatter {raw_scatter(product):.4f} dex is {r:.2f}x the "
                f"median of the {len(peers)} other pools on the same species/band/"
                f"holding -- grossly inconsistent with them (threshold {ANOMALY_RATIO})."))

    return tuple(out)


def scatter_ratio(product: dict, peers: list) -> float | None:
    """This product's line-to-line scatter over the LEAVE-ONE-OUT median of its peers.

    Leave-one-out because a product must not be part of the baseline it is judged against;
    median because one bad pool in the group must not drag the baseline up to meet it.
    """
    mine = raw_scatter(product)
    if mine is None:
        return None
    # Only peers whose statistic means the same thing. Mixing bases here would compare a
    # scatter with a standard error and call the difference an anomaly.
    mine_basis = stat_basis_of(product)
    others = [raw_scatter(p) for p in peers if stat_basis_of(p) == mine_basis]
    others = [x for x in others if x is not None and x > 0]
    if len(others) < ANOMALY_MIN_PEERS:
        return None
    med = statistics.median(others)
    if med <= 0:
        return None
    return mine / med


def evaluate_feed(doc: dict) -> dict:
    """{key_of(product): (Ineligible, ...)} for every product in the live list."""
    live = doc.get("products") or []
    groups: dict = {}
    for p in live:
        groups.setdefault(peer_group_of(p), []).append(p)
    out = {}
    for p in live:
        peers = [q for q in groups[peer_group_of(p)] if q is not p]
        out[key_of(p)] = evaluate(p, peers=peers)
    return out


def duplicate_live_cells(doc: dict) -> dict:
    """{key: count} for identities that appear more than once in `products[]`.

    A cell with two live products cannot say which one IS the product for that cell, so
    the feed is ambiguous and a consumer picks arbitrarily. This is reported separately
    from `evaluate` because the remedy is a judgement -- which of the two is current --
    and the gate does not get to make it.
    """
    counts: dict = {}
    for p in doc.get("products") or []:
        counts[key_of(p)] = counts.get(key_of(p), 0) + 1
    return {k: v for k, v in counts.items() if v > 1}
