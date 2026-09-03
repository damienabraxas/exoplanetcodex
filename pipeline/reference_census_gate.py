"""
pipeline/reference_census_gate.py — RYA-1173
============================================
RYA-946's freeze gate, made enforceable:

    "No element is FROZEN_READY_FOR_MEASUREMENT until this cross-reference is complete or a
     documented, approved source-publication exception exists."

🔴 WHY THIS EXISTS. RYA-1132 stamped `intake_status = FROZEN` on 164 Al manifest rows while the
mandatory AGSS21 line-set census did not exist for Al and no exception was recorded. Nothing
noticed for four days; RYA-1141 found it by reading, as check D3. A sentence in a ticket is not a
gate — this module is.

🔴 IT ASKS THE REGISTRY, NOT THE FILESYSTEM. D3's own version of this test globbed
`data/reference/asplund*` for a directory whose name mentioned the element. An empty `mkdir` passes
that, a half-staged download passes it, and a set that exists under a name the glob does not
anticipate fails it while being perfectly real. `reference_lineset.sets_for_element` answers from
the registry instead, and `census_state` then LOADS the set — a set that cannot be read is not a
census, however many bytes are on disk.

⚠️ THE EXCEPTION ROUTE IS NARROW ON PURPOSE. RYA-946 allows a "documented, approved
source-publication exception", so refusing to implement one would be refusing the standard as
written. But an exception is a bypass of a science gate, so it must name a ticket, an approver and
a date, and say WHICH publication problem it covers. `data/reference/census_exceptions.csv` is the
only place one may live, and `census_state` reports every element that leans on one rather than
folding it silently into a PASS.

⚠️ WHERE THE GATE BITES. Not the element ledger — `intake_status_ledger.csv` currently carries no
FROZEN element at all, so a guard pointed only there would pass forever without ever evaluating
anything. It bites on the per-line intake manifests, which is where FROZEN is actually written.
`tests/test_reference_census_gate_rya1173.py` asserts the sweep is non-empty for exactly that
reason, and ships a synthetic element with no reference set to prove the gate can still fail.
"""
from __future__ import annotations

import csv
from pathlib import Path

import pandas as pd

from pipeline import reference_lineset as rls

ROOT = Path(__file__).resolve().parents[1]

#: The one place an approved exception may live.
EXCEPTIONS = ROOT / "data" / "reference" / "census_exceptions.csv"
EXCEPTION_FIELDS = ("element", "ticket", "approved_by", "approved_on",
                    "source_publication_problem", "scope", "notes")

#: Every `intake_status` value that asserts the intake is closed. RYA-946's gate names
#: FROZEN_READY_FOR_MEASUREMENT; RYA-1132 wrote plain FROZEN and FROZEN_SOURCE_CONTROL for the same
#: act, and reading the gate narrowly enough to miss those is how it came to be missed the first
#: time. Any status that begins FROZEN is a freeze.
def is_frozen(status: object) -> bool:
    return str(status or "").strip().upper().startswith("FROZEN")


#: Where per-line intake manifests live. A manifest is any CSV under an audit directory that
#: carries both an `intake_status` column and a `species` column.
MANIFEST_GLOB = "data/audit/*/*.csv"


def _rel(path: Path, root: Path) -> str:
    """`path` relative to `root` when it is under it, else absolute.

    ⚠️ PATHS HERE ARE NOT ALL UNDER ONE ROOT, AND A BARE `relative_to` RAISES WHEN THEY ARE NOT.
    `check_all(root=...)` takes an arbitrary root so the gate can be pointed at a synthetic tree,
    and a monkeypatched ReferenceSet can point outside the repo entirely. Both happened while
    writing the tests, and the first one raised FROM INSIDE THE ERROR PATH -- the guard crashed at
    exactly the moment it had correctly detected a violation, which reads as a broken test rather
    than a caught bug. Cosmetic formatting must never be able to do that.
    """
    try:
        return str(path.relative_to(root))
    except ValueError:
        return str(path)


class CensusGateError(RuntimeError):
    """An element is frozen without the RYA-946 census, or without a recorded exception."""


def read_exceptions() -> dict[str, dict]:
    """Approved source-publication exceptions, by element. Missing file = no exceptions."""
    if not EXCEPTIONS.exists():
        return {}
    out = {}
    with EXCEPTIONS.open(newline="") as fh:
        for row in csv.DictReader(fh):
            if not (row.get("element") or "").strip():
                continue
            missing = [f for f in EXCEPTION_FIELDS if not (row.get(f) or "").strip()]
            if missing:
                raise CensusGateError(
                    f"census exception for {row['element']!r} is incomplete: {missing} empty. An "
                    f"exception bypasses a science gate; an unattributed one is worse than none.")
            out[row["element"].strip()] = row
    return out


def census_state(element: str) -> dict:
    """Is RYA-946's census discharged for this element, and on what evidence?"""
    names = rls.sets_for_element(element)
    sets = []
    for n in names:
        spec = rls.SETS[n]
        entry = {"line_set": n, "path": _rel(spec.path, ROOT), "ticket": spec.ticket}
        try:
            ref = rls.load(n)
            entry.update(loads=True, n_rows=int(len(ref)),
                         n_used=int(len(rls.measurable(ref))),
                         n_excluded_by_source=int(len(rls.excluded_by_source(ref))))
        except Exception as exc:                      # a set that will not load is not a census
            entry.update(loads=False, n_rows=0, n_used=0, error=f"{type(exc).__name__}: {exc}")
        sets.append(entry)

    exc = read_exceptions().get(str(element).strip())
    # A set with zero measurable rows is not a census either -- it is an empty directory with a
    # registry entry, which is the failure mode the registry lookup was meant to rule out.
    usable = [s for s in sets if s.get("loads") and s.get("n_used", 0) > 0]
    return {
        "element": element,
        "reference_sets": sets,
        "n_usable_sets": len(usable),
        "census_complete": bool(usable),
        "exception": exc,
        "satisfied": bool(usable) or exc is not None,
        "basis": ("registered reference set(s) " + ", ".join(s["line_set"] for s in usable)
                  if usable else
                  (f"APPROVED EXCEPTION {exc['ticket']} ({exc['approved_by']}, {exc['approved_on']})"
                   if exc else "NOTHING -- no reference set and no exception")),
    }


def frozen_elements_in_manifests(root: Path = ROOT) -> dict[str, list[dict]]:
    """Every element carrying a FROZEN* row in any per-line intake manifest, and where."""
    found: dict[str, list[dict]] = {}
    for path in sorted(root.glob(MANIFEST_GLOB)):
        try:
            head = pd.read_csv(path, nrows=0)
        except Exception:
            continue
        if "intake_status" not in head.columns or "species" not in head.columns:
            continue
        d = pd.read_csv(path, low_memory=False)
        frozen = d[d["intake_status"].map(is_frozen)]
        if frozen.empty:
            continue
        # "Al I" / "Al II" -> Al. Multi-species cells ("C I|C II") split on the same rule.
        for sp in frozen["species"].astype(str):
            for one in sp.split("|"):
                el = one.strip().split(" ")[0].strip()
                if not el or el.lower() == "nan":
                    continue
                rows = found.setdefault(el, [])
                hit = next((r for r in rows if r["manifest"] == str(path.relative_to(root))), None)
                if hit is None:
                    rows.append({"manifest": str(path.relative_to(root)),
                                 "n_frozen_rows": int(len(frozen))})
    return found


def check_all(root: Path = ROOT) -> list[str]:
    """Complaints, one per offending element. Empty list = the gate is satisfied everywhere.

    Returns rather than raises so a caller can report every element in one pass.
    """
    bad = []
    for element, where in sorted(frozen_elements_in_manifests(root).items()):
        st = census_state(element)
        if st["satisfied"]:
            continue
        n = sum(w["n_frozen_rows"] for w in where)
        bad.append(
            f"{element}: {n} manifest row(s) stamped FROZEN in "
            f"{', '.join(w['manifest'] for w in where)}, but RYA-946's Solar reference-line-set "
            f"census is not discharged -- {st['basis']}. Either build the reference set (the "
            f"RYA-1109/RYA-1173 pattern) or record an approved exception in "
            f"{_rel(EXCEPTIONS, root)}.")
    return bad


def require_census(element: str) -> dict:
    """The state, or raise. For a builder to call BEFORE it writes FROZEN."""
    st = census_state(element)
    if not st["satisfied"]:
        raise CensusGateError(
            f"refusing to freeze {element}: RYA-946's Solar reference-line-set census is not "
            f"discharged -- {st['basis']}. 'No element is FROZEN_READY_FOR_MEASUREMENT until this "
            f"cross-reference is complete or a documented, approved source-publication exception "
            f"exists.'")
    return st
