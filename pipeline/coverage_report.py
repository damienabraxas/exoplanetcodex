"""The definition of DONE for one element: every cell, every treatment, an outcome (RYA-1233).

Ryan, 2026-10-08: "did we hit every band? every instrument? every engine? with uncertainty?"
-- a question the orchestrator must answer itself. This table is that answer, built from the
orchestrator's own expansion (`run_matrix.expand`), its process holds, the last run report and
the published feed. Every (band x holding x ion x grade-pool x treatment) cell gets exactly one
outcome:

  PUBLISHED   a product in data/products/<star>/<El>.json, with its A, sigma and n
  HELD        the cell, or that treatment, cannot be produced -- with the stated reason
  NOT_RUN     applicable and not yet measured -- work is owed
  STALE       published, but its inputs (lines, gf, spectrum, culls, NLTE table, code) changed
              since the build -- work is owed

A cell's treatments are its deck's emissions (run_matrix.DECK_EMITS); ENGINE-B on the
ts-lte deck is the 1D-LTE fit relabelled and is not a separate cell.
An element is DONE when no cell is NOT_RUN or STALE. ELEMENT_PROTOCOL (docs/) "three honest outcomes"
map here: RESOLVED -> PUBLISHED, BOUNDED/UNRESOLVED -> HELD with its reason.
"""
from __future__ import annotations

import json
from pathlib import Path

from pipeline import band_policy
from pipeline import run_matrix as rm

#: Band order is the band policy's own (blue to red), never typed here.
BAND_ORDER = [p.name for p in band_policy.POLICIES]

#: A cell's treatments are what its deck EMITS (run_matrix.DECK_EMITS): ts-lte -> 1D-LTE +
#: ENGINE-A; each Gerber deck -> its own token. ENGINE-B on ts-lte is not a cell.
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
            last[(c["band"], c["holding"], c["ion"], c.get("engine") or "ts-lte", c.get("pool") or "")] = c
    # The budget stage's own verdicts (element_publish): a product built but held by its
    # RYA-587 budget, or skipped as empty, is HELD with that reason -- not NOT_RUN.
    budget_held = {}
    bud = rm.REPORT_DIR / "budget" / f"{star}_{element.replace(' ', '')}" / "report.json"
    if bud.exists():
        for o in json.loads(bud.read_text()).get("products", []):
            why = o.get("skip") or o.get("verdict")
            if why:
                budget_held[(o.get("band"), o.get("holding"), o.get("selector") or "",
                             o.get("key_treatment"))] = f"budget: {why}"
    # A PUBLISHED cell is only done while its inputs are the ones it was built on: the ledger
    # records the inputs hash each build ran against; recompute it now (run()'s own hash).
    ledger = json.loads(rm.LEDGER.read_text()) if rm.LEDGER.exists() else {}
    manifests = rm._manifest_paths(star)

    def stale(d) -> str:
        rec = ledger.get(rm.ledger_key(star, d, rm._band_of(d, star).name_declared))
        if rec is None:
            # Published without an orchestrator build record (a hand run, or before the
            # ledger): nothing proves it was built on today's inputs.
            return "no orchestrator build record for this cell"
        try:
            from pipeline.run_descriptor import resolve
            now = rm.inputs_hash(d, resolve(d, interpreter=None, ispec_dir=None),
                                 manifest_path=manifests.get(d.holding))
        except Exception:
            return ""
        return "" if now == rec.get("inputs_hash") else "inputs changed since the build"

    rows = []
    for ion in rm.graded_ions(symbol) or ["I"]:
        # The THREE grades only (Ryan): published sets (Reference), Codex, Deep. The
        # orchestrator's `reference` pool (every lab line, no depth gate) is Codex + Deep
        # combined and is not a grade of its own -- not an expected cell.
        for d in rm.expand(star, symbol, ions=[ion], methods=["synthesis"],
                           pools=["set", "codex", "deep"]):
            sel = _selector(d.pool)
            hold = rm.cell_process_hold(d)
            # The run report stores the pool as the feed token (pool_tier: GRADED / DEEPGRADED /
            # SET-<NAME>), not the orchestrator's pool name.
            run = last.get((d.band, d.holding, ion, d.engine_deck, rm.pool_tier(d) or ""))
            for t in rm.DECK_EMITS[d.engine_deck]:
                p = index.get((d.band, d.holding, ion, sel, t))
                if p is not None:
                    why = stale(d)
                    out = "STALE" if why else "PUBLISHED"
                elif hold:
                    out, why = "HELD", hold
                elif run is not None and run.get("status") in (rm.FAILED, rm.HELD, rm.BLOCKED):
                    out, why = "HELD", f"{run['status']}: {run.get('reason', '')}"
                elif (d.band, d.holding, sel, t) in budget_held:
                    out, why = "HELD", budget_held[(d.band, d.holding, sel, t)]
                elif t == "ENGINE-A" and (index.get((d.band, d.holding, ion, sel, "1D-LTE"))
                                          or (d.band, d.holding, sel, "1D-LTE") in budget_held):
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
    counts = {k: sum(r["outcome"] == k for r in rows)
              for k in ("PUBLISHED", "HELD", "NOT_RUN", "STALE")}
    return {"star": star, "element": element, "feed_version": feed.get("version"),
            "counts": counts, "done": counts["NOT_RUN"] == 0 and counts["STALE"] == 0,
            "rows": rows}


def render(cov: dict) -> str:
    lines = [f"COVERAGE  {cov['star']} {cov['element']}  feed v{cov['feed_version']}  "
             + "  ".join(f"{k} {v}" for k, v in cov["counts"].items())
             + ("  -> DONE" if cov["done"] else "  -> NOT DONE (NOT_RUN / STALE cells owe work)")]
    for r in sorted(cov["rows"], key=lambda r: (BAND_ORDER.index(r["band"]) if r["band"] in BAND_ORDER else 99,
                                                r["holding"], r["ion"], r["grade"], r["treatment"])):
        val = (f"{r['A']:.3f} +/- {r['sigma']:.3f} (n={r['n_lines']})" if r["A"] is not None
               and r["sigma"] is not None else "")
        lines.append(f"  {r['outcome']:9} {r['band']:11} {r['holding'][:30]:30} {r['ion']:2} "
                     f"{r['grade'][:22]:22} {r['treatment']:8} {val or r['reason'][:70]}")
    return "\n".join(lines)
