"""
Refinement-debt join (RYA-676) — the registry, rendered.
========================================================
ONE QUESTION, ANSWERED PER ELEMENT: *this row is owed — WHICH TICKET WOULD FIX IT?*

WHY THIS EXISTS
---------------
RYA-524's 27-element audit spawned children. The architectural ones got executed;
RYA-581 (Ba deblend), RYA-585 (Zr rescue) and RYA-565 (Eu adjudication) sat in
Backlog through eight architecture tickets and were noticed only when RYA-672
went looking, post-hoc. Nothing was broken — every artifact was individually
correct. The gap was that no surface anywhere carried the sentence "this owed row
has a ticket, and it has not fired." So the tracker said `owed` and the ticket
said `Backlog` and the two never met.

This module is where they meet. It reads
``data/audit/element_refinement_registry.csv`` (hand-maintained SSOT, RYA-676 §2A)
and renders it into one tracker column and one informational CI report.

IT DECIDES NOTHING
------------------
Whether a debt exists is a registry question, answered by a human against the
admission rule in that file's header. This module only joins and formats. It has
no fallback: an unreadable registry raises rather than emitting blanks, because a
blank in this column reads as "nothing owed" on exactly the rows the ticket exists
to make visible.

THE JOIN IS BY ELEMENT, NOT BY SPECIES — DELIBERATELY
-----------------------------------------------------
The tracker carries the ion the Codex REPORTS; the registry carries the ion the
DEBT is about, and those genuinely differ (RYA-475 is about Y I lines while Y is
reported as Y II; gold v3 holds `Eu,I` while every Eu measurement is Eu II —
RYA-683). Joining on the pair would silently drop those rows, which is the exact
failure mode this whole ticket is about. So the join is by element and the ion
travels inside the rendered label, where a reader can see the mismatch.

VALUE FORMS (RYA-676 §2B)
-------------------------
    "Backlog: RYA-680 (Ba engine-B harness)"     a ticket exists and has not fired
    "In Progress: RYA-585 (Zr rescue)"           refinement work executing
    "Done: RYA-592 (Mg 5528)"                    retained for provenance until the
                                                 result folds into the phase_c verdict
    "TBD - no resolving ticket (second Ba II line)"
                                                 debt established, NO ticket filed —
                                                 the signal to file one
    ""                                           no registry row: no known refinement
                                                 path, or nothing owed

An element with several debts renders all of them, `; `-joined, most-actionable
first. Ba II is the worked example: an unwired Engine B *and* an unconfirmed single
line are two different fixes and neither substitutes for the other.
"""

from __future__ import annotations

import argparse
import csv
import sys
from dataclasses import dataclass
from pathlib import Path

# Runnable as `python -m pipeline.refinement_debt_join` AND as a plain file path -- see
# the same note in pipeline/ledger_consistency_guard.py.
if __package__ in (None, ""):                                   # pragma: no cover
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from pipeline import state_surfaces  # noqa: E402

REPO_ROOT = state_surfaces.REPO_ROOT

#: Registry path resolved from the RYA-659 surface registry, never re-spelled here.
REGISTRY_REL = state_surfaces.REFINEMENT_REGISTRY
REGISTRY_PATH = REPO_ROOT / REGISTRY_REL

REQUIRED_COLUMNS = ("element", "ion", "situation", "resolving_ticket", "ticket_state",
                    "provenance_ticket", "phase", "short_label", "notes")

#: The phase whose gate open debt is measured against. There is no machine source for
#: "which phase is active" — it is a Ryan decision carried on Linear labels — so it is
#: declared HERE, once, and every consumer imports it. Bump it when the phase closes;
#: do not let a second copy appear in a guard or a report.
ACTIVE_PHASE = "Solar Beta"

#: Rendered prefix per Linear state. Ordering is the render order: the most actionable
#: class first. `TBD` outranks `Done` because an unticketed debt is louder than a
#: resolved one, and `Done` rows are retained only for provenance.
_STATE_PREFIX = {
    "In Progress": "In Progress",
    "Todo": "Backlog",           # Linear's Todo and Backlog are one class to a reader
    "Backlog": "Backlog",
    "Done": "Done",
}
_RENDER_ORDER = ("In Progress", "Backlog", "TBD", "Done")

#: RYA-705 — the render classes that are OPEN DEBT.
#:
#: `In Progress` used to fall through. `open_debt` selected on
#: ``rendered.startswith("Backlog:")``, so a row whose ticket someone was actively
#: working on landed in neither bucket, vanished from the report, and was invisible to
#: the `--phase-close` gate. That inverts the intent: `_RENDER_ORDER` puts `In Progress`
#: FIRST because it is the most actionable class, and the selector then dropped exactly
#: that class. A phase could have closed over debt that was mid-flight.
#:
#: `Done` is deliberately NOT here — a discharged debt is not open — but it is no longer
#: silent either: `report()` prints the discharged rows, because "this debt went away"
#: is a claim that should be auditable rather than a row quietly leaving the file.
_OPEN_CLASSES = frozenset({"In Progress", "Backlog"})

#: The literal a consumer greps for. RYA-676 §2B fixes this string; the short label is
#: appended so an element with two unticketed debts does not render two identical cells.
TBD_TEXT = "TBD - no resolving ticket"


class RegistryError(RuntimeError):
    """The registry is missing or malformed. Never degraded to an empty registry."""


@dataclass(frozen=True)
class DebtRow:
    element: str
    ion: str
    situation: str
    resolving_ticket: str
    ticket_state: str
    provenance_ticket: str
    phase: str
    short_label: str
    notes: str

    @property
    def is_unticketed(self) -> bool:
        return self.resolving_ticket.strip().upper() == "TBD"

    @property
    def render_class(self) -> str:
        """One of _RENDER_ORDER."""
        if self.is_unticketed:
            return "TBD"
        prefix = _STATE_PREFIX.get(self.ticket_state.strip())
        if prefix is None:
            raise RegistryError(
                f"{self.element} {self.ion} / {self.situation}: ticket_state "
                f"{self.ticket_state!r} is not one of {sorted(_STATE_PREFIX)}. A state this "
                f"module cannot render must not be guessed at — fix the registry row.")
        return "In Progress" if prefix == "In Progress" else prefix

    def render(self) -> str:
        if self.is_unticketed:
            return f"{TBD_TEXT} ({self.short_label})"
        return f"{self.render_class}: {self.resolving_ticket} ({self.short_label})"


def load_registry(path: Path = REGISTRY_PATH) -> list[DebtRow]:
    """Read the registry. Loud on absence, on a missing column, and on a row that
    omits its `provenance_ticket` — the admission rule is that a debt without an
    establishing ticket is not a debt, so an unprovenanced row is a registry bug."""
    if not path.exists():
        raise RegistryError(
            f"refinement registry not found at {path} — it is the SSOT for "
            f"'this row is owed -> this ticket resolves it' (RYA-676 §2A). Restore it; "
            f"do not fall back to an empty registry, which would render every owed "
            f"element as having nothing to fix.")
    with path.open(encoding="utf-8") as fh:
        rows = list(csv.DictReader(_uncommented(fh)))
    if not rows:
        raise RegistryError(f"{path.name} has a header but no rows")
    missing = set(REQUIRED_COLUMNS) - set(rows[0])
    if missing:
        raise RegistryError(f"{path.name} is missing column(s) {sorted(missing)}")

    out: list[DebtRow] = []
    seen: set[tuple[str, str, str]] = set()
    for raw in rows:
        row = DebtRow(**{k: (raw[k] or "").strip() for k in REQUIRED_COLUMNS})
        if not row.provenance_ticket:
            raise RegistryError(
                f"{row.element} {row.ion} / {row.situation}: no provenance_ticket. The "
                f"admission rule is that a row exists only where a landed artifact or a "
                f"filed ticket ESTABLISHED the debt — an unprovenanced row is speculation.")
        if not row.short_label:
            raise RegistryError(
                f"{row.element} {row.ion} / {row.situation}: no short_label — it is what "
                f"renders in the tracker cell")
        key = (row.element, row.ion, row.situation)
        if key in seen:
            raise RegistryError(f"duplicate registry row {key}")
        seen.add(key)
        row.render_class          # validates ticket_state at load, not at render time
        out.append(row)
    return out


def _uncommented(lines):
    for line in lines:
        if not line.startswith("#"):
            yield line


def by_element(rows: list[DebtRow] | None = None) -> dict[str, list[DebtRow]]:
    """element -> its debt rows, in render order."""
    rows = load_registry() if rows is None else rows
    grouped: dict[str, list[DebtRow]] = {}
    for row in rows:
        grouped.setdefault(row.element, []).append(row)
    for element in grouped:
        grouped[element].sort(key=lambda r: (_RENDER_ORDER.index(r.render_class),
                                             r.situation))
    return grouped


def debt_cell(element: str, grouped: dict[str, list[DebtRow]]) -> str:
    """The `refinement_debt` tracker cell for one element. Empty when nothing is owed."""
    return "; ".join(r.render() for r in grouped.get(element, []))


def assert_elements_known(tracker_elements: set[str],
                          rows: list[DebtRow] | None = None) -> None:
    """Every registry element must exist in the tracker.

    A registry row naming an element the tracker does not carry would render nowhere
    and be invisible — the same silent-drop this module exists to prevent, one level up.
    """
    rows = load_registry() if rows is None else rows
    orphans = sorted({r.element for r in rows} - tracker_elements)
    if orphans:
        raise RegistryError(
            f"refinement registry names element(s) {orphans} that the element status "
            f"tracker does not carry. Reconcile them — a row that joins to nothing is a "
            f"debt nobody will ever see.")


# ─────────────────────────────────────────────────────────────────────────────
#  the informational report (RYA-676 §2C)
# ─────────────────────────────────────────────────────────────────────────────
def open_debt(tracker_tier_by_element: dict[str, str],
              rows: list[DebtRow] | None = None,
              phase: str = ACTIVE_PHASE) -> dict[str, list[dict]]:
    """The two classes of open debt in the active phase, keyed for the guard's --json.

    * ``refinement_debt_open`` — tier is owed, a ticket EXISTS, and it has not fired.
      "We know how to fix this and haven't." This is the RYA-672 class, and it is what
      RYA-676 §2C's informational check counts.
    * ``refinement_debt_unticketed`` — tier is owed and NO ticket exists. Strictly worse,
      and deliberately counted separately: it cannot be cleared by firing something, so
      folding it into the first count would make the first count un-actionable.

    `owed` here means the RYA-522 freeze tier, taken from the tracker, never re-derived.
    """
    rows = load_registry() if rows is None else rows
    owed_tiers = {"owed", "nlte_owed", "curation_owed"}
    out: dict[str, list[dict]] = {"refinement_debt_open": [],
                                  "refinement_debt_unticketed": [],
                                  "refinement_debt_discharged": []}
    for row in rows:
        if row.phase != phase:
            continue
        tier = (tracker_tier_by_element.get(row.element) or "").strip().lower()
        if tier not in owed_tiers:
            continue
        rendered = row.render()
        record = {"element": row.element, "ion": row.ion, "tier": tier,
                  "situation": row.situation, "refinement_debt": rendered,
                  "provenance_ticket": row.provenance_ticket}
        if row.is_unticketed:
            out["refinement_debt_unticketed"].append(record)
        elif row.render_class in _OPEN_CLASSES:
            out["refinement_debt_open"].append(record)
        else:                                  # "Done" — discharged, kept for provenance
            out["refinement_debt_discharged"].append(record)
    return out


def tracker_tiers() -> dict[str, str]:
    """element -> tier, read from the generated tracker (RYA-654)."""
    import pandas as pd
    path = REPO_ROOT / state_surfaces.TRACKER
    if not path.exists():
        raise RegistryError(f"element status tracker not found at {path}")
    df = pd.read_csv(path, comment="#")
    tiers: dict[str, str] = {}
    for _, r in df.iterrows():
        tier = str(r.get("tier") or "").strip()
        if tier == "diagnostic":          # the Fe II arbiter row asserts no verdict
            continue
        tiers[str(r["element"]).strip()] = tier
    return tiers


def report(phase: str = ACTIVE_PHASE) -> int:
    """Print the informational report. ALWAYS returns 0 — see RYA-676 §2C.

    This never blocks a merge. Seven of the registry's rows are `NO_MODEL_ATOM` Engine-B
    gaps that need an atom acquired; nobody can clear them this week, and a guard that is
    permanently red is a guard nobody reads. Visibility is the deliverable, not gating.
    Phase-close tickets escalate it explicitly (`--phase-close`).
    """
    rows = load_registry()
    tiers = tracker_tiers()
    assert_elements_known(set(tiers), rows)
    buckets = open_debt(tiers, rows, phase)

    backlog = buckets["refinement_debt_open"]
    unticketed = buckets["refinement_debt_unticketed"]
    print(f"REFINEMENT DEBT (RYA-676) — phase {phase!r}, INFORMATIONAL, never blocking")
    print(f"  registry: {len(rows)} row(s) over "
          f"{len({r.element for r in rows})} element(s) — {REGISTRY_REL}")
    print(f"  owed with a ticket that has NOT fired : {len(backlog)}")
    for r in backlog:
        print(f"    - {r['element']} {r['ion']}: {r['refinement_debt']}")
    print(f"  owed with NO ticket filed             : {len(unticketed)}"
          "   <- these need a ticket, not a run")
    for r in unticketed:
        print(f"    - {r['element']} {r['ion']}: {r['refinement_debt']}")
    # RYA-705: discharged rows are NOT debt and are not counted — but they are printed,
    # because a row silently leaving the report is how a debt gets lost, and "the ticket
    # closed" is a claim worth being able to audit against the element's actual state.
    # `.get` rather than `[]`: `open_debt` always populates all three keys, but this
    # report must survive a substituted bucket source — the one thing it may never do is
    # raise, since RYA-676 §2C makes exit 0 a property of the report, not of its content.
    discharged = buckets.get("refinement_debt_discharged", [])
    if discharged:
        print(f"  discharged (ticket Done, kept for provenance): {len(discharged)}"
              "   <- NOT counted as debt")
        for r in discharged:
            print(f"    - {r['element']} {r['ion']}: {r['refinement_debt']}")
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--report", action="store_true",
                    help="print the informational open-debt report (always exit 0)")
    ap.add_argument("--phase", default=ACTIVE_PHASE,
                    help=f"phase to scope the report to (default {ACTIVE_PHASE!r})")
    ap.add_argument("--phase-close", action="store_true",
                    help="ESCALATE: exit 1 if any owed element in the phase still carries "
                         "un-fired refinement debt. For phase-close / freeze tickets only "
                         "(RYA-677); never the default.")
    args = ap.parse_args(argv)

    try:
        rc = report(args.phase)
        if args.phase_close:
            tiers = tracker_tiers()
            buckets = open_debt(tiers, None, args.phase)
            n = len(buckets["refinement_debt_open"]) + len(buckets["refinement_debt_unticketed"])
            if n:
                print(f"\nPHASE-CLOSE ESCALATION: {n} open refinement debt row(s) in "
                      f"{args.phase!r}. A phase does not close over un-fired refinement "
                      f"work — fire it, or ratify each row as deferred and move its "
                      f"`phase` in {REGISTRY_REL}.", file=sys.stderr)
                return 1
        return rc
    except RegistryError as exc:
        print(f"REFINEMENT DEBT REGISTRY FAILED: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
