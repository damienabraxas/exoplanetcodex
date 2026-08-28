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
6. `SUPERSEDED` — the record itself declares a replacement. Derived from the record, never
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

WHY `ANOMALOUS_SCATTER` IS MEASURED ON RAW SCATTER, NOT ON `sigma_stat`
----------------------------------------------------------------------
`sigma_stat` is a STANDARD ERROR, so it carries a 1/sqrt(n) that has nothing to do with
data quality: a clean 3-line pool has a large `sigma_stat` for a reason that is not a
defect. Multiplying it back gives the line-to-line scatter, `sigma_stat * sqrt(n)`, which
is the quantity the ticket actually describes and the one that is comparable across pools
of different size.

MEASURED on the committed feed, under this module's own rule (leave-one-out median,
>= 3 peers, 37 comparable products, 4 of them the incomplete Amarsi pools):

    statistic      null max     flagged min     separation
    sigma_stat       3.137          4.641          1.48x
    raw scatter      1.889          5.085          2.69x

⚠️ SO THE HONEST STATEMENT IS NOT "sigma_stat CANNOT SEPARATE THEM" -- it can, on this
feed. It is that its window is 1.8x narrower, and that `ANOMALY_RATIO = 3.0` sits INSIDE
the raw-scatter window while falling BELOW the sigma_stat null maximum of 3.137 -- i.e.
the very same constant would have flagged a legitimate product had the statistic been the
standard error. (An earlier draft of this docstring claimed the two populations overlap on
`sigma_stat`. They do not, and the test that says so is what caught it.)

The threshold is not fitted to the gap. `test_product_eligibility_rya1092` asserts the
measured null, asserts that null and flagged do not overlap, asserts `ANOMALY_RATIO` lies
strictly between them, and asserts the window stays wider than 2x -- so no value in the
admissible range changes the verdict, and if new data narrows the window the test fails
and forces a re-derivation instead of letting a stale constant decide (RYA-161, and
[[a tolerance needs a measured null]]).
"""
from __future__ import annotations

import math
import statistics
from dataclasses import dataclass

#: The canonical product identity. Imported concept, not a second definition -- the same
#: nine fields `scripts/publish_product.KEY_FIELDS` uses. Duplicated as a literal here
#: only because importing a script from a pipeline module would invert the dependency;
#: `test_product_eligibility_rya1092` asserts the two lists are equal.
KEY_FIELDS = ("element", "ion", "band", "instrument", "holding",
              "tier", "selector", "route", "treatment")

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


class EligibilityError(RuntimeError):
    """A product could not be evaluated -- refuse, never pass by default."""


@dataclass(frozen=True)
class Ineligible:
    """One reason a product may not be live. Products carry a LIST: a record can be
    incomplete AND anomalous, and collapsing that to a single code would let the fix for
    one hide the other."""
    code: str
    detail: str


def key_of(product: dict) -> str:
    return "|".join(str(product.get(k) or "") for k in KEY_FIELDS)


def species_of(product: dict) -> tuple:
    return (str(product.get("element") or "").strip(),
            str(product.get("ion") or "").strip() or None)


def is_real_species(product: dict) -> bool:
    el, ion = species_of(product)
    return (el, ion) in REAL_SPECIES or (el, None) in REAL_SPECIES


def raw_scatter(product: dict) -> float | None:
    """Line-to-line scatter, recovered from the standard error: sigma_stat * sqrt(n).

    This is the quantity that is comparable across pools of different size. `sigma_stat`
    is not -- see the module docstring.
    """
    s, n = product.get("sigma_stat"), product.get("n_lines")
    if s is None or not n or float(n) <= 0:
        return None
    return float(s) * math.sqrt(float(n))


def peer_group_of(product: dict) -> tuple:
    """Products that measure the SAME thing by different methods: same species, same band,
    same holding. Treatment and tier deliberately vary -- comparing a pool only against
    itself would make the check vacuous."""
    return (product.get("element"), product.get("ion"),
            product.get("band"), product.get("holding"))


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
    others = [raw_scatter(p) for p in peers]
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
