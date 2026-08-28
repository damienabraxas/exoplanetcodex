#!/usr/bin/env python3
"""RYA-1100 — withdraw the ENGINE-B products that duplicate a 1D-LTE product.

🔴 WHAT THIS IS NOT. It is not a de-duplication by value. Two products with the same
number are not necessarily the same product -- Fe I VIS DEEPGRADED and its
`-LOCALRENORM` sibling are different measurements that happen to share a cell, and a
value-only match would collapse real physics.

WHAT MAKES THESE ONE PRODUCT is the axis registry, not the numbers:

    pipeline/treatment_axes.py::LEGACY
        "1D-LTE":   dict(scale="1D-LTE", model="none", gf="kurucz")
        "ENGINE-B": dict(scale="1D-LTE", model="none", gf="kurucz")   # identical

Every declared axis agrees. The only thing the `ENGINE-B` token ever added was a pinned
route (`_ROUTE_BY_LABEL`), and RYA-906 moved route off the label onto the handler -- which
left the label naming nothing. An `ENGINE-B` product IS a `1D-LTE` product measured by the
synthesis handler. The identical VALUES are the CONSEQUENCE of that, and this script
requires BOTH: same axes AND byte-identical published numbers. Either alone is not enough.

WITHDRAWAL, NOT DELETION (RYA-711). Records move to `superseded[]` in full, carrying
`superseded_by` (the surviving twin's identity) and a reason. Reversible by hand, and the
pool is served, so nothing silently vanishes.

🔴 THE 2 ROWS THIS DELIBERATELY WILL NOT TOUCH -- AND WHY THEY ARE NOT A SCIENCE CALL.
Two `ENGINE-B` products sit on route=PROFILEFIT, contradicting `_ROUTE_BY_LABEL
["ENGINE-B"] = "synth"`. That first read as a line-selection disagreement. It is not.
The source CSV settles it:

    treatment  handler             route   A      n_lines  n_excluded
    1D-LTE     ProfileFitHandler   ew      7.498  6        0
    ENGINE-A   ProfileFitHandler   ew      7.433  5        1
    ENGINE-B   SynthesisHandler    synth   7.431  5        1   <-- declares SYNTH

The row declares `handler=SynthesisHandler` and `route=synth` IN ITS OWN FILE. The feed
publishes it as PROFILEFIT because `publish_product.py` takes the route from a `--route`
COMMAND-LINE FLAG and applies it file-wide, never reading the per-row `handler` / `route`
/ `route_basis` columns -- and this CSV is a PROFILEFIT-named file that happens to also
contain one synthesis row. RYA-869 put `handler` on the product to be the authoritative
witness for exactly this question, and the publisher overrides it with a filename.

So these two are a real, distinct SYNTH product wearing the wrong route -- not a duplicate
and not a disagreement. They stay LIVE and REPORTED here because correcting them changes
a published identity field (`route`), which moves the product to a different cell and is
Ryan's call to make, not a tidy-up. The publisher fix stops it recurring; see
`publish_product.normalise`.
"""
from __future__ import annotations

import argparse
import datetime as _dt
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline import treatment_axes  # noqa: E402

#: The identity of a product CELL -- what makes two rows the same measurement slot.
#: `treatment` is deliberately absent: it is the axis being collapsed.
CELL = ("element", "ion", "band", "instrument", "holding", "tier", "selector", "route")

#: The published numbers that must agree byte-for-byte. `n_excluded` IS included: two
#: products over the same lines that excluded different ones are different selections,
#: and that is exactly how the 2 PROFILEFIT rows are caught.
VALUES = ("A", "sigma_stat", "sigma_syst", "n_lines", "n_excluded")

ALIAS = "ENGINE-B"
REASON_CODE = "DUPLICATE_RETIRED_LABEL"


def _now() -> str:
    return _dt.datetime.now(_dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def cell_of(p: dict) -> tuple:
    return tuple(p.get(k) for k in CELL)


def values_of(p: dict) -> tuple:
    return tuple(p.get(k) for k in VALUES)


def assert_axes_are_identical(alias: str, canonical: str) -> None:
    """The PREMISE, re-measured at run time rather than trusted from the docstring.

    If someone gives ENGINE-B its own physics tomorrow, this script must stop working.
    Comparing the axis dicts is what makes that automatic -- a comment cannot.
    """
    a, c = treatment_axes.LEGACY[alias], treatment_axes.LEGACY[canonical]
    if a != c:
        raise SystemExit(
            f"REFUSING: {alias} and {canonical} no longer declare identical axes.\n"
            f"  {alias:10s} {a}\n  {canonical:10s} {c}\n"
            f"This script exists only because they were the same product. They are not "
            f"any more, so withdrawing one would destroy a real measurement.")


def classify(feed: dict) -> tuple[list, list, list]:
    """(duplicates, contradictions, orphans) — every alias product, accounted for."""
    live = feed["products"]
    by_cell: dict[tuple, list] = {}
    for p in live:
        by_cell.setdefault(cell_of(p), []).append(p)

    dupes, contra, orphans = [], [], []
    for p in live:
        if p.get("treatment") != ALIAS:
            continue
        canonical = treatment_axes.DEPRECATED_ALIASES[ALIAS][0]
        twins = [q for q in by_cell[cell_of(p)] if q.get("treatment") == canonical]
        if not twins:
            orphans.append((p, None))
        elif any(values_of(q) == values_of(p) for q in twins):
            dupes.append((p, next(q for q in twins if values_of(q) == values_of(p))))
        else:
            contra.append((p, twins[0]))
    return dupes, contra, orphans


def ident(p: dict) -> str:
    return (f"{p['element']} {p['ion']} {p['band']} {p['instrument']} {p['holding']} "
            f"{p['tier']}/{p['selector']} {p['route']} {p['treatment']}")


def withdraw(feed: dict, dupes: list) -> int:
    stamp = _now()
    canonical = treatment_axes.DEPRECATED_ALIASES[ALIAS][0]
    moved = 0
    for p, twin in dupes:
        rec = dict(p)
        rec["superseded_at"] = stamp
        rec["superseded_reason_code"] = REASON_CODE
        rec["superseded_by"] = {k: twin.get(k) for k in (*CELL, "treatment")}
        rec["superseded_reason"] = (
            f"`{ALIAS}` is a retired spelling of `{canonical}` measured by the synthesis "
            f"handler, NOT a separate engine: pipeline/treatment_axes.py::LEGACY declares "
            f"identical axes for the two (scale=1D-LTE, model=none, gf=kurucz), and the "
            f"only thing the token added was a pinned route that RYA-906 moved onto the "
            f"handler. This record and the surviving `{canonical}` product in the same "
            f"cell carry byte-identical {', '.join(VALUES)} -- one measurement published "
            f"twice under two names, which put two rows with the same display name on the "
            f"live element page. Withdrawn, not deleted (RYA-711): the record is kept in "
            f"full and `superseded_by` names the surviving twin, so this is reversible.")
        feed.setdefault("superseded", []).append(rec)
        moved += 1
    keep = {id(p) for p, _ in dupes}
    feed["products"] = [p for p in feed["products"] if id(p) not in keep]
    return moved


def rederive_display(feed: dict) -> list[tuple]:
    """Re-derive `display` on every surviving product. A REGENERATION, not a hand edit.

    🔴 WHY THIS IS ALLOWED TO CHANGE A PUBLISHED FIELD. RYA-1080 pins every field that
    existed at its baseline, and `display` is one -- but `display` is not a measurement.
    It is a LABEL COMPUTED FROM other published fields (`treatment`, `route`, `gf`) by
    code under version control, and those inputs do not move here. Correcting the code
    that derives it is a code fix; the label following is arithmetic, not a reconciliation.

    RYA-1100 takes OWNERSHIP of the field in exchange, with a check that is strictly
    stronger than the pin it replaces: every published `display` must EQUAL what the axis
    registry derives for that row. A pin can be satisfied by never touching the field; a
    derivability assertion cannot be satisfied by a hand edit at all.

    No measured field is touched here, and none is permitted to be -- see the caller.
    """
    from scripts.publish_product import display_name
    changed = []
    # ⚠️ EVERY POOL, not just `products`. The withdrawal pools are SERVED -- a quarantined
    # product is shown with its reason, which is the whole point of quarantining rather
    # than deleting (RYA-711) -- so a stale label there is just as visible as a live one.
    # Re-deriving only `products` left 7 rows in `quarantine`/`archive` still rendering
    # ENGINE-A synthesis products as "EW", which is the defect this fixes.
    for pool in ("products", "superseded", "quarantine", "archive"):
        for p in feed.get(pool) or []:
            was = p.get("display")
            now = display_name(p["treatment"], gf="kurucz", route=p.get("route"))
            if was != now:
                p["display"] = now
                changed.append((p["treatment"], p.get("route"), was, now))
    return changed


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--feed", default=str(ROOT / "data/products/solar/Fe.json"))
    ap.add_argument("--apply", action="store_true", help="write; default is a dry run")
    a = ap.parse_args()

    canonical = treatment_axes.DEPRECATED_ALIASES[ALIAS][0]
    assert_axes_are_identical(ALIAS, canonical)

    path = Path(a.feed)
    feed = json.loads(path.read_text())
    before = len(feed["products"])
    dupes, contra, orphans = classify(feed)

    print(f"feed {path}  v{feed.get('version')}  live={before}")
    print(f"\n{ALIAS} products: {len(dupes) + len(contra) + len(orphans)}")
    print(f"  exact duplicates of a {canonical} twin : {len(dupes)}  -> superseded[]")
    print(f"  route/selection contradictions        : {len(contra)}  -> LEFT LIVE, reported")
    print(f"  no {canonical} twin in the cell          : {len(orphans)}  -> LEFT LIVE, reported")

    if dupes:
        print("\n=== WITHDRAWN (duplicate of a live twin) ===")
        for p, t in dupes:
            print(f"  {ident(p)}\n      A={p['A']} n={p['n_lines']}  == twin {t['treatment']}"
                  f" A={t['A']} n={t['n_lines']}")
    if contra:
        print("\n🔴 === LEFT LIVE — RYA-1100 OPEN QUESTION, RYAN'S CALL ===")
        for p, t in contra:
            print(f"  {ident(p)}")
            print(f"      {ALIAS:9s} A={p['A']} n={p['n_lines']} n_excl={p['n_excluded']}")
            print(f"      {t['treatment']:9s} A={t['A']} n={t['n_lines']} n_excl={t['n_excluded']}")
            print(f"      published route={p['route']} but the SOURCE CSV row declares "
                  f"handler=SynthesisHandler / route=synth. The publisher stamped a "
                  f"file-wide --route flag over the row's own witness (RYA-869). This is "
                  f"a real SYNTH product, not a duplicate. Correcting `route` changes a "
                  f"published identity field -> Ryan's call.")
    if orphans:
        print("\n⚠️ === LEFT LIVE — no canonical twin, so nothing to supersede TO ===")
        for p, _ in orphans:
            print(f"  {ident(p)}  A={p['A']} n={p['n_lines']}")

    if not a.apply:
        print("\n(dry run — pass --apply to write)")
        return 0
    if not dupes:
        print("\nnothing to withdraw; feed untouched")
        return 0

    moved = withdraw(feed, dupes)

    # every measured field, snapshotted so the re-derive cannot touch one by accident
    MEASURED = ("A", "sigma_stat", "sigma_syst", "n_lines", "n_excluded", "treatment",
                "route", "tier", "selector", "holding", "instrument", "band", "ion")
    before_vals = [tuple(p.get(k) for k in MEASURED) for p in feed["products"]]

    relabelled = rederive_display(feed)
    after_vals = [tuple(p.get(k) for k in MEASURED) for p in feed["products"]]
    if before_vals != after_vals:
        raise SystemExit("REFUSING: re-deriving the display name moved a measured field")
    if relabelled:
        print(f"\n=== RE-DERIVED display on {len(relabelled)} surviving products ===")
        agg = {}
        for t, r, w, n in relabelled:
            agg.setdefault((t, r, w, n), 0)
            agg[(t, r, w, n)] += 1
        for (t, r, w, n), c in sorted(agg.items(), key=lambda x: -x[1]):
            print(f"  n={c:3d}  {t:32s} route={r:11s}\n           {w!r}\n        -> {n!r}")

    v = str(feed.get("version", "1.0"))
    feed["version"] = f"{v.split('.')[0]}.{int(v.split('.')[1]) + 1}"
    feed["updated_at"] = _now()
    # ⚠️ MATCH `publish_product`'s ENCODING, WHICH IS THE DEFAULT ensure_ascii=True.
    # Writing with ensure_ascii=False renders "·" literally instead of "\u00b7" and
    # reformats EVERY line carrying a display name -- a whole-file diff that hides the
    # real change and that the next publish would immediately revert.
    path.write_text(json.dumps(feed, indent=2) + "\n")
    print(f"\nwithdrew {moved}: live {before} -> {len(feed['products'])}, "
          f"superseded[] {len(feed['superseded'])}, v{feed['version']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
