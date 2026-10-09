"""The definition of DONE for one element: every cell, every treatment, an outcome (RYA-1233).

Ryan, 2026-10-08: "did we hit every band? every instrument? every engine? with uncertainty?"
-- a question the orchestrator must answer itself. This table is that answer, built from the
orchestrator's own expansion (`run_matrix.expand`), its process holds, the last run report and
the published feed. Every (band x holding x ion x grade-pool x treatment) cell gets exactly one
outcome:

  PUBLISHED   a product in data/products/<star>/<El>.json, with its A, sigma and n
  HELD        the cell, or that treatment, cannot be produced -- with the stated reason
  NOT_RUN     applicable and not yet measured -- the only outcome that means work is owed

ENGINE-B on the ts-lte deck is the 1D-LTE fit relabelled and is not a separate cell.
An element is DONE when no cell is NOT_RUN. ELEMENT_PROTOCOL (docs/) "three honest outcomes"
map here: RESOLVED -> PUBLISHED, BOUNDED/UNRESOLVED -> HELD with its reason.
"""
from __future__ import annotations

import json
from pathlib import Path

from pipeline import band_policy
from pipeline import run_matrix as rm

#: Band order is the band policy's own (blue to red), never typed here.
BAND_ORDER = [p.name for p in band_policy.POLICIES]

#: The treatment families a synthesis cell emits on the ts-lte deck.
TREATMENTS = ("1D-LTE", "ENGINE-A")
#: Orchestrator pool -> the feed's `selector` token.
POOL_SELECTOR = {"codex": "GRADED", "deep": "DEEPGRADED", "reference": "REFERENCE"}
GRADE_OF_POOL = {"codex": "Codex Grade", "deep": "Deep Grade", "reference": "Reference (all lab gf)"}


def _selector(pool: str | None) -> str:
    if pool and pool.startswith(rm.SET_POOL_PREFIX):
        return "SET-" + pool[len(rm.SET_POOL_PREFIX):]
    return POOL_SELECTOR.get(pool or "", "")


def _grade(pool: str | None) -> str:
    return "Reference Grade" if pool and pool.startswith(rm.SET_POOL_PREFIX) else GRADE_OF_POOL.get(pool or "", "")


def coverage(star: str, element: str, *, report_dir: Path | None = None) -> dict:
    symbol, _ = rm.split_symbol(element)
    feed = rm.load_feed(star, symbol) or {}
    products = feed.get("products", [])
    index = {}
    for p in products:
        index[(p.get("band"), p.get("holding"), str(p.get("ion")), p.get("selector") or "",
               p.get("treatment"))] = p
    last = {}
    rep = (report_dir or rm.REPORT_DIR) / f"{star}_{element.replace(' ', '')}_latest.json"
    if rep.exists():
        for c in json.loads(rep.read_text()).get("cells", []):
            last[(c["band"], c["holding"], c["ion"], c.get("pool") or "")] = c
    rows = []
    for ion in rm.graded_ions(symbol) or ["I"]:
        # The THREE grades only (Ryan): published sets (Reference), Codex, Deep. The
        # orchestrator's `reference` pool (every lab line, no depth gate) is Codex + Deep
        # combined and is not a grade of its own -- not an expected cell.
        for d in rm.expand(star, symbol, ions=[ion], engines=["ts-lte"], methods=["synthesis"],
                           pools=["set", "codex", "deep"]):
            sel = _selector(d.pool)
            hold = rm.cell_process_hold(d)
            run = last.get((d.band, d.holding, ion, d.pool or ""))
            for t in TREATMENTS:
                p = index.get((d.band, d.holding, ion, sel, t))
                if p is not None:
                    out, why = "PUBLISHED", ""
                elif hold:
                    out, why = "HELD", hold
                elif run is not None and run.get("status") in (rm.FAILED, rm.HELD, rm.BLOCKED):
                    out, why = "HELD", f"{run['status']}: {run.get('reason', '')}"
                elif t == "ENGINE-A" and index.get((d.band, d.holding, ion, sel, "1D-LTE")):
                    out, why = "HELD", ("no NLTE correction for this pool's lines in the "
                                        "element's ENGINE-A table")
                else:
                    out, why = "NOT_RUN", ""
                rows.append({"band": d.band, "holding": d.holding, "ion": ion,
                             "grade": _grade(d.pool), "pool": d.pool, "treatment": t,
                             "outcome": out, "reason": why,
                             "A": p.get("A") if p else None,
                             "sigma": p.get("sigma_reported") if p else None,
                             "n_lines": p.get("n_lines") if p else None})
    counts = {k: sum(r["outcome"] == k for r in rows) for k in ("PUBLISHED", "HELD", "NOT_RUN")}
    return {"star": star, "element": element, "feed_version": feed.get("version"),
            "counts": counts, "done": counts["NOT_RUN"] == 0, "rows": rows}


def render(cov: dict) -> str:
    lines = [f"COVERAGE  {cov['star']} {cov['element']}  feed v{cov['feed_version']}  "
             + "  ".join(f"{k} {v}" for k, v in cov["counts"].items())
             + ("  -> DONE" if cov["done"] else "  -> NOT DONE (NOT_RUN cells owe work)")]
    for r in sorted(cov["rows"], key=lambda r: (BAND_ORDER.index(r["band"]) if r["band"] in BAND_ORDER else 99,
                                                r["holding"], r["ion"], r["grade"], r["treatment"])):
        val = (f"{r['A']:.3f} +/- {r['sigma']:.3f} (n={r['n_lines']})" if r["A"] is not None
               and r["sigma"] is not None else "")
        lines.append(f"  {r['outcome']:9} {r['band']:11} {r['holding'][:30]:30} {r['ion']:2} "
                     f"{r['grade'][:22]:22} {r['treatment']:8} {val or r['reason'][:70]}")
    return "\n".join(lines)
