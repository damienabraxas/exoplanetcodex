#!/usr/bin/env python3
"""RYA-1212 Part C — which published products does the 0.17 gate refuse, and what fixes each?

The gate (Part A) refuses at PUBLISH time by reading the budget. This asks the same
question of the feed as it stands: which LIVE products rest on the blanket, and what is
the resolution path for each.

🔴 THE ANSWER IS NOT "GO AND MEASURE N LABORATORY gf". `gf_graded` is a two-branch switch
that cannot say "56 of 57 lines are laboratory", so it returns the blanket on the strength
of the ones that are not — which is why RYA-855 moved 0 of 36 bars. So for every refused
product this reports the SPLIT, not just the verdict: how many lines are already graded and
how many actually block it. The near-UV DEEPGRADED pool is the worked example — 56 of 57
GF-LAB, blocked by ONE line, and RYA-1209 resolved it by dropping that line (57 -> 55),
which took its bar off the blanket without measuring anything new.
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline import error_budget as eb   # noqa: E402

FEED = ROOT / "data" / "products" / "solar" / "Fe.json"
OUT = ROOT / "data" / "results" / "rya1212"
BLANKET = eb.UNGRADED_GF_SYSTEMATIC_DEX


def budget_for(rec) -> Path | None:
    """The budget beside a product's artifact — the same resolution publish_product uses."""
    path = (rec.get("provenance") or {}).get("copied_to") \
        or (rec.get("provenance") or {}).get("path") or ""
    if not path.endswith("_products.csv"):
        return None
    stem = path[: -len("_products.csv")]
    for cand in (ROOT / (stem + "_budgets.txt"),
                 ROOT / (re.sub(r"_ENGINE-[A-Z-]+$", "", stem) + "_budgets.txt")):
        if cand.exists():
            return cand
    # products published from an off-repo scratch path still have the artifact in-repo
    tail = Path(stem).name
    for cand in (ROOT / "data/results/band_products" / (tail + "_budgets.txt"),
                 ROOT / "data/results/band_products"
                 / (re.sub(r"_ENGINE-[A-Z-]+$", "", tail) + "_budgets.txt")):
        if cand.exists():
            return cand
    return None


def split_from(text: str) -> dict:
    """The decider's own line/grade split, so the report says what BLOCKS each product."""
    m = re.search(r"MIXED POOL: (\d+) of (\d+) .*? lines are GF-LAB", text or "")
    if not m:
        return {}
    graded, total = int(m.group(1)), int(m.group(2))
    rest = re.search(r"the rest are ([^.]*)\.", text or "")
    return {"n_graded": graded, "n_lines": total, "n_blocking": total - graded,
            "blocking_grades": rest.group(1).strip() if rest else ""}


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    doc = json.loads(FEED.read_text())
    refused, checked, nobudget = [], 0, []

    for r in doc["products"]:
        checked += 1
        syst = float(r.get("sigma_syst") or 0.0)
        bud = budget_for(r)
        text = bud.read_text() if bud else ""
        by_value = abs(syst - BLANKET) < 5e-4 or syst >= BLANKET
        by_budget = eb.carries_ungraded_gf(text) if text else None
        if by_budget is None and not by_value:
            nobudget.append(r["holding"])
        if not (by_value or by_budget):
            continue
        row = {
            "holding": r["holding"], "ion": r["ion"], "band": r["band"],
            "tier": r["tier"], "treatment": r["treatment"],
            "line_set": r.get("line_set"), "grade": r.get("grade"),
            "sigma_syst": syst, "n_lines": r.get("n_lines"),
            "refused_by_published_value": by_value,
            "refused_by_budget_term": bool(by_budget),
            "budget": str(bud.relative_to(ROOT)) if bud else None,
        }
        row.update(split_from(text))
        refused.append(row)

    for row in refused:
        # The resolution path is a property of WHAT BLOCKS the pool, so it is derived from
        # the split rather than written per product.
        nb, ng = row.get("n_blocking"), row.get("n_graded")
        if row.get("line_set") == "asplund":
            row["resolution"] = (
                "RYA-1211 measured this pool: 5 of 21 GF-LAB and 12 GF-NIST = 17 of 21 "
                "carrying a citable per-line sigma (RMS 0.0475 dex), but CITED_COVERAGE_MIN "
                "is 90% so the cited route stays refused at 81%. 🔴 RYA-1211 does NOT "
                "resolve these, contrary to this ticket's Part C assumption. The route is "
                "RYA-968 per-line: 17 lines resolve as `cited`, and the remaining 4 need "
                "`fallback_sigma_dex`, which gf_empirical.Thresholds leaves UNDECLARED by "
                "design (a borrowed constant is not a control). Deriving it is the blocker.")
        elif nb is not None and nb <= 3:
            row["resolution"] = (
                f"{ng} of {row['n_lines']} lines are already GF-LAB and {nb} block the "
                f"pool ({row.get('blocking_grades')}). Cheapest path is drop/re-source "
                f"those {nb}, exactly as RYA-1209 did for the near-UV pool (57 -> 55).")
        else:
            row["resolution"] = (
                "wire RYA-968 per-line grading (empirical_gf_sigma_dex, now accepted by "
                "error_budget.build) once its thresholds are derived, or re-source the "
                "unresolved lines.")

    summary = {
        "ticket": "RYA-1212",
        "feed_version": doc.get("version"),
        "live_products_checked": checked,
        "refused": len(refused),
        "blanket_dex": BLANKET,
        "products": refused,
        "products_without_a_reachable_budget": len(nobudget),
        "note_on_budget_artifacts": (
            "20 of 92 committed *_budgets.txt still carry the blanket, but most are "
            "SUPERSEDED artifacts, not live products — e.g. the near-UV n=57 budget was "
            "replaced by RYA-1209's n=55/54 re-run, whose products publish 0.1122/0.1124. "
            "The gate fires on what is BEING published, so the live count is what matters."),
        "no_value_changed": True,
    }
    (OUT / "rya1212_refused_products.json").write_text(json.dumps(summary, indent=2) + "\n")

    print(f"feed v{doc.get('version')}: {checked} live products, {len(refused)} REFUSED by the gate")
    for row in refused:
        print(f"  {row['holding'][:33]:35}{row['treatment']:18}tier={row['tier']:5}"
              f"n={row['n_lines']:3} sigma_syst={row['sigma_syst']}  {row.get('grade')}")
    print(f"  wrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
