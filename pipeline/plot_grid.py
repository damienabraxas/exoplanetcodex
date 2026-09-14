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


def is_displayable(p: dict, only_deep: set) -> bool:
    """GRADED or REFERENCE, or DEEPGRADED in a section that has nothing else."""
    return (p.get("tier") in ("GRADED", "REFERENCE")
            or section_of(p) in only_deep)


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
    for p in products:
        if is_displayable(p, only_deep):
            buckets[display_section_of(p)][p.get("display")].append(p)

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
        # ⚠️ THE KEY FORMAT IS PUBLISHED, NOT RE-IMPLEMENTED BY THE READER. The renderer
        # joins cells to products on `product_key`; if it hard-coded the field order it
        # would silently stop matching the day KEY_FIELDS changed, and every cell would
        # render as N/A -- a join bug wearing the costume of a coverage gap.
        "key_fields": list(_pe.KEY_FIELDS),
        "sections": sections,
    }
