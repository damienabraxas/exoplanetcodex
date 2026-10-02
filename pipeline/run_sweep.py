"""Runs every canonical element of one star through `run_matrix`, Fe first; decides no science.

RYA-1233 (the front door over RYA-1222's executor). Owns: the ORDER elements are
attempted in, the Fe gate, and one sweep report naming every element. Does NOT
own: what a cell is, whether it can run, whether its work is current, or how it is
reported -- every one of those is `run_matrix.run()`'s, called once per element
with its contract unchanged.

WHY FE FIRST AND WHY A GATE
---------------------------
Every other element's [X/Fe] product is relative to Fe, so a sweep that builds 25
elements on top of an Fe that produced nothing has built 25 numbers on nothing and
reports them beside one another as if they were comparable. The gate stops the
sweep there and SAYS so: the remaining elements are written into the report as
NOT_RUN_FE_GATE, a documented status rather than an absence (RYA-1187).

THE ELEMENT AXIS IS DERIVED, NEVER TYPED
----------------------------------------
The vocabulary and its order are `elements_master.json`'s. The ONLY rule this file
adds is `FE_FIRST`. `elements_master.json` lists `Fe` and `Fe II` as separate
targets; `run_matrix.run("Fe")` expands EVERY graded ion of Fe, so a sweep that ran
both master entries as written would build every Fe II cell twice under two report
names. The bare entry therefore runs the ions that no other master entry claims --
which for Fe is exactly Fe I -- read off the master list and the graded pool, not
off a literal.
"""
from __future__ import annotations

import json
import sys
from dataclasses import dataclass
from pathlib import Path

from pipeline import run_matrix

#: RYA-1233: these symbols run before everything else, and gate it. Fe because every
#: other element's [X/Fe] product is referenced to it. The order WITHIN this group,
#: and of everything after it, is elements_master.json's own.
FE_FIRST = ("Fe",)

#: The cell statuses that mean "this cell's work exists" -- the gate counts these.
#: run_matrix's own definition, not a second copy of it.
GATE_OK = run_matrix.TERMINAL_OK

#: Element-level statuses in the sweep report. RAN and NOT_RUN_FE_GATE are the
#: ticket's. The other two exist because an element whose matrix cannot even be
#: BUILT (no graded line of it, an unknown deck in --engine) would otherwise be a
#: hole in the report, which is the one outcome this layer is not allowed to have.
RAN = "RAN"
NOT_RUN_FE_GATE = "NOT_RUN_FE_GATE"
#: The graded pool holds no line of this element at the ions it would run. An HONEST
#: EMPTY (a pool finding -- RYA-945 owes the grading), not a failure: on 2026-10-01
#: 13 of the 27 solar targets were in this state, and folding them into MATRIX_ERROR
#: would have made every sweep exit non-zero for a fact no run can change.
NO_GRADED_POOL = "NO_GRADED_POOL"
MATRIX_ERROR = "MATRIX_ERROR"      # run_matrix refused to build the matrix; reason given
CRASHED = "CRASHED"                # an unexpected exception; loud, and the sweep continues
ELEMENT_STATUSES = (RAN, NOT_RUN_FE_GATE, NO_GRADED_POOL, MATRIX_ERROR, CRASHED)

FE_GATE_REASON = ("Fe produced no usable cell; every other element's [Fe/H]-dependent "
                  "product would be built on nothing.")


@dataclass(frozen=True)
class SweepEntry:
    """One element of the sweep: the master-list spelling, and the ions it runs."""
    target: str                 # as elements_master.json spells it ("Fe", "Fe II", "Si")
    ions: tuple[str, ...]       # () = run_matrix's default (every graded ion)
    label: str                  # what the report and stdout call it ("Fe I")
    plan_error: str = ""        # non-empty: the entry cannot be run; said why

    @property
    def symbol(self) -> str:
        return run_matrix.split_symbol(self.target)[0]


def sweep_plan() -> list[SweepEntry]:
    """Every canonical target, FE_FIRST symbols first, otherwise in master order."""
    targets = run_matrix.canonical_elements()
    # Ions that have their OWN master entry ("Fe II" -> Fe: {II}).
    claimed: dict[str, set[str]] = {}
    for t in targets:
        sym, ion = run_matrix.split_symbol(t)
        if ion:
            claimed.setdefault(sym, set()).add(ion)

    entries: list[SweepEntry] = []
    for t in targets:
        sym, ion = run_matrix.split_symbol(t)
        if ion or sym not in claimed:
            entries.append(SweepEntry(t, (), t))
            continue
        rest = tuple(i for i in run_matrix.graded_ions(sym) if i not in claimed[sym])
        if not rest:
            entries.append(SweepEntry(t, (), t, plan_error=(
                f"every graded ion of {sym} is claimed by its own master entry "
                f"({', '.join(f'{sym} {i}' for i in sorted(claimed[sym]))}), so the bare "
                f"{t!r} entry has no cell of its own. Running it unfiltered would rebuild "
                f"those ions' cells a second time.")))
            continue
        entries.append(SweepEntry(t, rest, f"{sym} {'+'.join(rest)}"))

    first = [e for e in entries if e.symbol in FE_FIRST]
    return first + [e for e in entries if e.symbol not in FE_FIRST]


def executed(doc: dict) -> int:
    """Stage dispatches this run actually made, from the element report's own cells.

    A dry run lists a WOULD_RUN cell's steps by bare name and a reused step is tagged
    `:reused`; only `:ok` / `:FAILED` were dispatched. This is the number an
    idempotent re-run must show as 0.
    """
    return sum(1 for c in doc.get("cells", []) for st in c.get("steps", [])
               if st.endswith((":ok", ":FAILED")))


def would_run_published(doc: dict) -> int:
    """Dry run only: WOULD_RUN cells whose product is ALREADY in the published feed.

    These are the cells a real run would most likely turn DONE -- run_matrix calls a
    cell DONE when its stages succeed and the feed holds its product. The gate does
    NOT count them (a dry run cannot know the stages will succeed); the number is
    reported beside the verdict so a projected HALT is not read as "Fe has nothing".
    """
    return sum(1 for c in doc.get("cells", [])
               if c.get("status") == run_matrix.WOULD_RUN and c.get("A") is not None)


def _pool_has(e: SweepEntry) -> bool:
    """Does the graded pool hold any line this entry would run? run_matrix's own answer."""
    graded = run_matrix.graded_ions(e.symbol)
    want = e.ions or ((run_matrix.split_symbol(e.target)[1],)
                      if run_matrix.split_symbol(e.target)[1] else tuple(graded))
    return any(i in graded for i in want)


def _zero_counts() -> dict[str, int]:
    return {s: 0 for s in run_matrix.STATUSES}


def _line(entry: dict) -> str:
    """One stdout line per element. Every status, zeroes included -- see render()."""
    k = entry["counts"]
    head = f"{entry['label']:<8}"
    if entry["status"] != RAN:
        return f"{head} {entry['status']}  {entry.get('reason', '')}".rstrip()
    return (head + " " + " ".join(f"{s} {k[s]}" for s in run_matrix.STATUSES)
            + f"  (executed {entry['executed']})")


def run_sweep(star: str, *, bands: list[str] | None = None,
              instruments: list[str] | None = None, engines: list[str] | None = None,
              methods: list[str] | None = None, dry_run: bool = False, interpreter: str | None = None,
              ispec_dir: str | None = None, step_timeout: int = 7200,
              report_dir: Path | None = None, echo: bool = True) -> dict:
    """Every canonical element of `star`, Fe first, through `run_matrix.run()`.

    Loud-fail-continue per element, exactly as per cell: an element with FAILED
    cells, or one whose matrix cannot be built, is recorded and the next is run.
    The ONE thing that stops the sweep is the Fe gate, and it stops it in the
    report, not by raising.
    """
    plan = sweep_plan()
    dest = report_dir or run_matrix.REPORT_DIR
    commit = run_matrix._code_commit()
    out: list[dict] = []
    fe_gate = None

    def say(msg: str) -> None:
        if echo:
            print(msg, flush=True)

    say(f"SWEEP  star={star}  {'DRY-RUN' if dry_run else 'EXECUTE'}  "
        f"{len(plan)} canonical elements  commit={commit}")
    for i, e in enumerate(plan):
        # ── the Fe gate: evaluated once, after the last FE_FIRST entry ─────────
        if fe_gate is None and e.symbol not in FE_FIRST:
            fe_ran = [r for r in out if r["symbol"] in FE_FIRST]
            usable = sum(r["counts"][s] for r in fe_ran for s in GATE_OK)
            fe_gate = "PASS" if usable else "HALT"
            if fe_gate == "HALT":
                say(f"\n!!! FE GATE HALT: {FE_GATE_REASON}")
                if dry_run:
                    pub = sum(r.get("would_run_with_published_product", 0) for r in fe_ran)
                    say(f"!!! DRY-RUN: 0 Fe cells are DONE/SKIP now; {pub} Fe WOULD_RUN "
                        f"cell(s) already have a published product and would turn DONE "
                        f"if their stages succeed.")
                    say("!!! DRY-RUN: nothing executes, so the remaining elements are "
                        "still expanded to show their matrices. A real sweep stops HERE.")

        entry = {"element": e.target, "label": e.label, "symbol": e.symbol,
                 "ions": list(e.ions) or None, "status": RAN, "report": None,
                 "counts": _zero_counts(), "cells_total": 0, "executed": 0}
        if fe_gate == "HALT" and not dry_run:
            entry.update(status=NOT_RUN_FE_GATE, reason=FE_GATE_REASON)
        elif not _pool_has(e):
            ion_txt = ", ".join(e.ions) if e.ions else (
                run_matrix.split_symbol(e.target)[1] or "any ion")
            entry.update(status=NO_GRADED_POOL, reason=(
                f"the graded pool in {run_matrix.CANONICAL_GF.name} holds no line of "
                f"{e.symbol} at {ion_txt}: no cell can be built. A POOL finding (grade "
                f"lines first, RYA-945), not a run failure."))
        elif e.plan_error:
            entry.update(status=MATRIX_ERROR, reason=e.plan_error)
            print(f"!!! {e.label}: {e.plan_error}", file=sys.stderr, flush=True)
        else:
            try:
                doc = run_matrix.run(
                    star, e.target, ions=list(e.ions) or None, bands=bands,
                    instruments=instruments, engines=engines, methods=methods,
                    dry_run=dry_run,
                    interpreter=interpreter, ispec_dir=ispec_dir,
                    step_timeout=step_timeout, report_dir=report_dir, echo=False)
                entry.update(report=doc.get("_report_path"), counts=doc["counts"],
                             cells_total=doc["cells_total"], executed=executed(doc))
                if dry_run:
                    entry["would_run_with_published_product"] = would_run_published(doc)
            except run_matrix.MatrixError as exc:
                entry.update(status=MATRIX_ERROR, reason=str(exc))
                print(f"!!! {e.label}: matrix not built -- {exc}", file=sys.stderr, flush=True)
            except Exception as exc:                                # noqa: BLE001
                # Loud-fail-continue: one element's crash must not erase the other
                # 26 from the report. It is recorded, printed, and fails the exit code.
                entry.update(status=CRASHED, reason=f"{type(exc).__name__}: {exc}")
                print(f"!!! {e.label}: CRASHED -- {type(exc).__name__}: {exc}",
                      file=sys.stderr, flush=True)
        if fe_gate == "HALT" and dry_run:
            entry["fe_gate_would_block"] = True
        out.append(entry)
        say(_line(entry))

    if fe_gate is None:                       # every target was in FE_FIRST
        fe_ran = [r for r in out if r["symbol"] in FE_FIRST]
        fe_gate = "PASS" if sum(r["counts"][s] for r in fe_ran for s in GATE_OK) else "HALT"

    doc = _write(star, out, fe_gate=fe_gate, dry_run=dry_run, commit=commit,
                 interpreter=interpreter, ispec_dir=ispec_dir, dest=dest)
    say(verdict(doc))
    return doc


def _write(star: str, entries: list[dict], *, fe_gate: str, dry_run: bool, commit: str,
           interpreter: str | None, ispec_dir: str | None, dest: Path) -> dict:
    totals = _zero_counts()
    for r in entries:
        for s, n in r["counts"].items():
            totals[s] += n
    canonical = run_matrix.canonical_elements()
    doc = {
        "schema": "codex.orchestrator_sweep/1",
        "star": star,
        "generated_at": run_matrix._utc(),
        "code_commit": commit,
        "dry_run": bool(dry_run),
        "environment": {"interpreter": interpreter or None, "ispec_dir": ispec_dir or None},
        "fe_first": list(FE_FIRST),
        "fe_gate": fe_gate,
        "fe_gate_reason": FE_GATE_REASON if fe_gate == "HALT" else None,
        # In a dry run nothing executes, so the gate cannot have stopped anything; it
        # reports what a real sweep WOULD do against the current state.
        "fe_gate_enforced": fe_gate == "HALT" and not dry_run,
        "canonical_count": len(canonical),
        "element_statuses": {s: sum(r["status"] == s for r in entries)
                             for s in ELEMENT_STATUSES},
        "elements": entries,
        "totals": totals,
        "cells_total": sum(r["cells_total"] for r in entries),
        "executed": sum(r["executed"] for r in entries),
    }
    # No silent gap, checked rather than trusted: every canonical target is in the
    # report exactly once, and the totals account for every cell.
    assert sorted(r["element"] for r in entries) == sorted(canonical), \
        "a canonical element is missing from (or duplicated in) the sweep report"
    assert sum(totals.values()) == doc["cells_total"], "a cell carries no terminal status"

    dest.mkdir(parents=True, exist_ok=True)
    stamp = doc["generated_at"].replace(":", "").replace("-", "")
    text = json.dumps(doc, indent=2) + "\n"
    (dest / f"{star}_SWEEP_{stamp}.json").write_text(text, encoding="utf-8")
    latest = dest / f"{star}_SWEEP_latest.json"
    latest.write_text(text, encoding="utf-8")
    try:
        doc["_report_path"] = str(latest.relative_to(run_matrix.ROOT))
    except ValueError:
        doc["_report_path"] = str(latest)
    return doc


def verdict(doc: dict) -> str:
    """The final line. ASCII only (Cloudflare WAF)."""
    t, es = doc["totals"], doc["element_statuses"]
    gate = doc["fe_gate"] + ("" if doc["fe_gate"] == "PASS" or doc["fe_gate_enforced"]
                             else " (dry-run: not enforced)")
    return (f"VERDICT  {doc['star']}  fe_gate {gate}  elements "
            + " ".join(f"{s} {es[s]}" for s in ELEMENT_STATUSES)
            + f" of {doc['canonical_count']}  cells "
            + " ".join(f"{s} {t[s]}" for s in run_matrix.STATUSES)
            + f" of {doc['cells_total']}  executed {doc['executed']}\nsweep report: {doc.get('_report_path', '')}")


def failed(doc: dict) -> bool:
    """Exit-code truth: any FAILED cell, any element not RAN, or an enforced gate halt."""
    return bool(doc["totals"][run_matrix.FAILED]
                or any(r["status"] in (MATRIX_ERROR, CRASHED) for r in doc["elements"])
                or doc["fe_gate_enforced"])
