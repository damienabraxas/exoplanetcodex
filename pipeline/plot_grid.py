"""RYA-1102 — the error-plot grid: fixed rows per section, N/A where absent.

🔴 THE DEFECT THIS REPLACES. `solar-report.js::plotRows` emitted ONE ROW PER PRODUCT,
grouped by band, from a build-time artifact. So a section's row count was whatever the
data happened to contain (Kitt Peak rendered 19), a model with no product for that
instrument simply vanished rather than reading N/A, and re-running a product could not
change the plot at all because the artifact was static.

Ryan's rule, implemented here rather than in JS so the site cannot drift from the physics:

  * GRADED only. DEEPGRADED is a secondary product -- documented in its own section, not
    showcased. EXCEPTION: a band where ONLY DEEPGRADED exists (near-UV) shows it, because
    "no graded product in this band" and "no product at all" are different facts.
  * 🔴 RYA-1213 ADDS REFERENCE TO THE SHOWCASE, AND GIVES IT ITS OWN SECTION. Reference
    Grade is the laboratory pool with the depth gate not applied; RYA-1213 requires the
    forest to draw it in every band that has one. It CANNOT share a section with the
    Codex product, because outside VIS the two measure the same holding, the same band
    and the same treatment -- the axis row name (`display`) is identical -- so one
    bucket would hold both and `_pick` would silently drop the one with fewer lines.
    In NIR and H the two pools are the SAME LINES, so the dropped row and the kept row
    would carry the same number under one grade label: a reader would see one product
    where there are two, and no trace of which. The section key therefore carries the
    resolved `line_set`.
  * EW and synthesis are SEPARATE PRODUCTS, both graded. A profile fit and a flux fit
    measure different line pools; they are two rows, never one row to collapse.
  * Section = ion x instrument x holding x band x line_set. The two telluric holdings stay apart:
    KPNO's kurucz2005 and molecfit differ by up to 0.26 dex, so collapsing them would
    silently choose a correction on the reader's behalf.

The ROWS come from `treatment_axes.plot_row_axis()` -- derived, ordered, never typed.
"""
from __future__ import annotations

import collections

from pipeline import product_eligibility as _pe
from pipeline import treatment_axes


def _sections_with_only_deepgraded(products: list[dict]) -> set:
    """Sections with NO Codex (GRADED) product, so the showcase rule would empty them.

    ⚠️ KEYED ON `section_of`, NOT ON THE DISPLAY SECTION, AND THE DIFFERENCE IS THE WHOLE
    EXCEPTION. The question this answers is "does this band on this holding have a Codex
    product" -- a question about the band. Asking it per line_set would make every deep
    section trivially deep-only, and the near-UV exception would quietly become a rule
    that showcases Deep Grade everywhere.

    🔴 RYA-1213 — THE TEST WAS `== {"DEEPGRADED"}` AND A REFERENCE PRODUCT EVICTED THE
    DEEP ONE. Publishing the near-UV Fe II Reference cells made that section's tier set
    {DEEPGRADED, REFERENCE}, which is no longer exactly {DEEPGRADED}, so the exception
    stopped firing and both Deep rows VANISHED from the plot -- a published product
    removed from the site by the arrival of a different one, silently. The intent was
    always "show Deep where there is no CODEX product", and a Reference product is not a
    Codex product, so the condition is now the absence of GRADED. near-UV Fe I is
    unaffected (still no GRADED, for the reason RYA-1213 Step 3 measured: the band has one
    lab line at or below the depth gate and a one-line pool has no scatter).
    """
    #: ⚠️ STILL ON `tier` HERE, DELIBERATELY. This asks "does this band on this holding
    #: have a CODEX product", and RYA-1213's fix was to test for the ABSENCE OF GRADED
    #: rather than for an exact tier set. Rewriting it in terms of grade would be a second
    #: change wearing the first one's clothes, and `is_displayable` consuming the result is
    #: what has to agree with the site, not this.
    tiers = collections.defaultdict(set)
    for p in products:
        tiers[section_of(p)].add(p.get("tier"))
    return {k for k, v in tiers.items() if "DEEPGRADED" in v and "GRADED" not in v}


def section_of(p: dict) -> tuple:
    return (p.get("ion"), p.get("instrument"), p.get("holding"), p.get("band"))


def display_section_of(p: dict) -> tuple:
    """The section a product is DRAWN in — `section_of` plus the pool it measured.

    RYA-1213: a Reference and a Codex product can agree on every other axis, including
    the row name, so the pool has to be part of what separates two sections or one of
    them is dropped without a word. Resolved through `reference_lineset`, never stored
    twice, so this cannot drift from the identity key the renderer joins on.
    """
    from pipeline.reference_lineset import line_set_for_product
    return section_of(p) + (line_set_for_product(p),)


#: 🔴 THE SHOWCASE RULE IS ON `grade`, NOT ON `tier`, AND THAT DISTINCTION IS THE BUG IT
#: FIXES. It read `tier in ("GRADED", "REFERENCE")`, and four products carry grade
#: "Reference Grade" while sitting in tier "ALL" -- the AGSS21 line-set 3D-NLTE cells, ONE
#: OF WHICH IS THE PUBLISHED Fe I HEADLINE. They were refused, and refused silently, so the
#: number the site reports at the top of the page was absent from the forest below it.
#:
#: `grade` is the published PROVENANCE CATEGORY -- what pool a product measured and how it
#: was selected -- which is exactly the question "should the showcase draw this" asks.
#: `tier` answers a narrower one and needed a new literal every time a tier appeared, which
#: is how ALL was missed. Verified on the live feed before the change: grade maps 1:1 onto
#: tier for every product that was already displayable (Codex Grade 59 <-> GRADED, Deep
#: Grade 29 <-> DEEPGRADED), so moving the rule to grade admits exactly those four and
#: changes nothing else. Ryan ruled them IN, 2026-09-25.
SHOWCASE_GRADES = ("Codex Grade", "Reference Grade")

#: Secondary by policy: "DEEPGRADED is a secondary product -- documented in its own
#: section, not showcased", with the deep-only exception below.
SECONDARY_GRADES = ("Deep Grade",)


def _why_not_displayed(p: dict) -> str:
    """Said in the terms the RULE uses, so a reader need not re-derive it."""
    grade = p.get("grade")
    if grade in SECONDARY_GRADES:
        return (f"grade {grade!r} is secondary and this section has a showcase product, so "
                f"the graded-only rule holds it back (tier {p.get('tier')!r})")
    if grade not in SHOWCASE_GRADES:
        return (f"grade {grade!r} is UNRECOGNISED -- it is in neither {SHOWCASE_GRADES} nor "
                f"{SECONDARY_GRADES}, so nothing decided whether the forest should draw it. "
                f"Held back and reported rather than dropped; give the grade a rule "
                f"(tier {p.get('tier')!r})")
    return (f"grade {grade!r} is a showcase grade, so this entry is a BUG in the reporter "
            f"itself -- it should have been drawn (tier {p.get('tier')!r})")


def is_displayable(p: dict, only_deep: set) -> bool:
    """A showcase grade, or a secondary one in a section that has nothing else.

    ⚠️ AN UNRECOGNISED GRADE IS NOT DRAWN, AND IT IS NOT AN EXCEPTION EITHER. My first cut
    RAISED here, and that was the wrong instrument twice over: `write_feed` calls `build()`
    on every publish, so a partial product dict could fail a PUBLISH over a field only the
    PLOT needs -- it broke `test_uncertainty_contract_rya587` immediately. The reason to be
    loud was that the old code dropped products silently; `not_displayed` now solves that
    directly, so an unknown grade is held back AND NAMED rather than crashing the run or
    vanishing. Report, do not crash, do not hide.
    """
    grade = p.get("grade")
    return grade in SHOWCASE_GRADES or section_of(p) in only_deep


def _pick(candidates: list[dict]) -> dict:
    """One product for a slot.

    🔴 A DEPRECATED ALIAS NEVER WINS A SLOT FROM THE LABEL IT ALIASES. Under RYA-1100's
    derived names the two `ENGINE-B` rows render "EW · 1D-LTE" and were BEATING the
    genuine product: harps showed 7.431/n=5 instead of 7.498/n=6, kpno 7.479/n=9 instead
    of 7.445/n=13. `ENGINE-B` is a retired spelling of `1D-LTE`, so it cannot represent
    1D-LTE on a plot where the real 1D-LTE product exists.

    Ties beyond that go to the larger line pool, and a tie there is reported by the
    caller rather than broken silently.
    """
    real = [c for c in candidates
            if c.get("treatment") not in treatment_axes.DEPRECATED_ALIASES]
    return sorted(real or candidates, key=lambda p: -(p.get("n_lines") or 0))[0]


def build(products: list[dict], *, include_pending: bool = False) -> dict:
    """{"axis": [...], "sections": [...]} — every section carries exactly len(axis) cells.

    A cell with no product is present with `product: None`. That is the whole point: the
    grid draws the empty cells too, so "we did not measure this" is visible instead of
    being indistinguishable from "this model does not exist".
    """
    axis = treatment_axes.plot_row_axis(include_pending=include_pending)
    only_deep = _sections_with_only_deepgraded(products)

    buckets: dict = collections.defaultdict(lambda: collections.defaultdict(list))
    # 🔴 `off_axis` DOES NOT CATCH A PRODUCT THE TIER RULE REFUSED, and for a long time
    # nothing did. It reports a product that entered a section and found no row; a product
    # `is_displayable` rejects never enters one, so it left the grid with NO TRACE ANYWHERE
    # in the feed -- 14 of 160 Fe products, and `off_axis` was empty on all 38 sections, so
    # the number looked like zero. Two of those 14 are the reason this matters: tier `ALL`
    # carries grade "Reference Grade" and one of them is the PUBLISHED Fe I headline, absent
    # from the forest with nothing saying so. Refusing to draw a product is a decision; not
    # recording the refusal is a defect (the same rule as RYA-711 and RYA-844).
    not_displayed = []
    for p in products:
        if is_displayable(p, only_deep):
            buckets[display_section_of(p)][p.get("display")].append(p)
        else:
            not_displayed.append({
                "product_key": _pe.key_of(p),
                "tier": p.get("tier"), "grade": p.get("grade"),
                "selector": p.get("selector"), "line_set": display_section_of(p)[-1],
                "reason": _why_not_displayed(p),
            })

    sections = []
    for key in sorted(buckets, key=lambda k: tuple(str(x) for x in k)):
        cells, unplaced = [], dict(buckets[key])
        for row in axis:
            cands = unplaced.pop(row["name"], []) if row["name"] else []
            hit = _pick(cands) if cands else None
            # 🔴 A KEY, NOT A COPY. Duplicating the product's numbers into the grid would
            # create a second place they can be wrong, and the two would drift the first
            # time one was re-published (RYA-353). The renderer joins on `products[]`.
            cells.append({"row": row["name"], "pending": row["pending"],
                          "product_key": _pe.key_of(hit) if hit else None,
                          "alternates": len(cands) - 1 if cands else 0})
        sections.append({
            "ion": key[0], "instrument": key[1], "holding": key[2], "band": key[3],
            # RYA-1213 — the pool this section drew, so the renderer can LABEL it. Two
            # sections now differ only here, and an unlabelled pair would read as a
            # duplicate rather than as two grades of one cell.
            "line_set": key[4],
            # ⚠️ THE FLAG IS ABOUT THIS SECTION, NOT ABOUT THE BAND. The renderer prints
            # "DEEPGRADED (no graded product in this band)" from it, and a Reference
            # section in a band that happens to have no Codex product would otherwise
            # wear that caption while drawing the lab pool. True only where this section
            # IS the deep exception.
            "only_deepgraded": key[:4] in only_deep and key[4] == "our-deep-graded",
            "cells": cells,
            # 🔴 REPORTED, NEVER DROPPED. A product whose display name is not on the axis
            # would otherwise disappear from the site with no trace. RYA-711's rule
            # applied to a renderer: if it is not shown, it must still be counted.
            "off_axis": sorted(unplaced),
        })
    return {
        "axis": [r["name"] for r in axis],
        # 🔴 EVERY PRODUCT THE TIER RULE HELD BACK, BY NAME. A reader can now subtract:
        # len(products) == (cells that resolved) + (alternates) + len(not_displayed), so a
        # product cannot leave the plot without appearing in one of the three. Nothing here
        # changes WHICH products are drawn -- it changes whether the omission is countable.
        "not_displayed": sorted(not_displayed, key=lambda r: r["product_key"]),
        # ⚠️ THE KEY FORMAT IS PUBLISHED, NOT RE-IMPLEMENTED BY THE READER. The renderer
        # joins cells to products on `product_key`; if it hard-coded the field order it
        # would silently stop matching the day KEY_FIELDS changed, and every cell would
        # render as N/A -- a join bug wearing the costume of a coverage gap.
        "key_fields": list(_pe.KEY_FIELDS),
        "sections": sections,
    }
