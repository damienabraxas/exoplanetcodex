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
  * EW and synthesis are SEPARATE PRODUCTS, both graded. A profile fit and a flux fit
    measure different line pools; they are two rows, never one row to collapse.
  * Section = ion x instrument x holding x band. The two telluric holdings stay apart:
    KPNO's kurucz2005 and molecfit differ by up to 0.26 dex, so collapsing them would
    silently choose a correction on the reader's behalf.

The ROWS come from `treatment_axes.plot_row_axis()` -- derived, ordered, never typed.
"""
from __future__ import annotations

import collections

from pipeline import product_eligibility as _pe
from pipeline import treatment_axes


def _sections_with_only_deepgraded(products: list[dict]) -> set:
    """Sections whose every product is DEEPGRADED, so the graded-only rule would empty them."""
    tiers = collections.defaultdict(set)
    for p in products:
        tiers[section_of(p)].add(p.get("tier"))
    return {k for k, v in tiers.items() if v == {"DEEPGRADED"}}


def section_of(p: dict) -> tuple:
    return (p.get("ion"), p.get("instrument"), p.get("holding"), p.get("band"))


def is_displayable(p: dict, only_deep: set) -> bool:
    """GRADED, or DEEPGRADED in a section that has nothing else."""
    return p.get("tier") == "GRADED" or section_of(p) in only_deep


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
            buckets[section_of(p)][p.get("display")].append(p)

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
            "only_deepgraded": key in only_deep,
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
