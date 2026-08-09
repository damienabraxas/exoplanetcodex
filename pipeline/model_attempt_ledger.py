"""RYA-695 — every model/grid ATTEMPTED per element, and what became of it.

WHY THIS EXISTS
===============
Ryan, on the element tracker:

    "it should show which element was chosen, where we got the damn thing and why
     we chose it. But it should also reference any model we tried."
    "Hence the register, hence the element tracker, because we keep fucking this up.
     We keep repeating work."

The tracker answered the first three well enough and the fourth not at all. Which
models were TRIED for an element — and which were superseded, which were rejected,
which were never staged — lived only in ticket prose, so every session that wanted
to know re-derived it from scratch. That is the repetition, and it is not
hypothetical:

  * **Ti** — Bergemann-2011 was production until RYA-546's vintage audit found its
    scaled-Drawin H-collisions inflated the correction; RYA-544/545 replaced it with
    the Mallinson-2024 ab-initio atom (+0.0506, not +0.108/+0.20). Two tickets and a
    whole atmosphere-vs-deck investigation (RYA-535/542) to establish it.
  * **Ba** — RYA-559's EW→COG 2.410 was superseded by RYA-581's in-window deblend
    2.237 once the pool EW was shown to be blend-inflated.
  * **Zr** — RYA-560's synthesis was superseded by RYA-585's deblend refit, which
    fixed the red_chi2 and STILL could not clear the sensitivity floor.

None of that was readable from a generated artifact. All of it was already written
down somewhere.

SOURCES — JOINED, NEVER RE-DERIVED
==================================
Nothing here re-runs an experiment or re-decides a disposition. It joins records
that already exist:

  1. ``data/nlte_grids/*.csv``            the Engine-A grids that were actually built
  2. ``data/nlte_grids/*.prov.json``      their provenance, incl. the ``supersedes``
                                          key that names what each one replaced
  3. ``config.constants.NLTE_CORRECTION_ELEMENTS``  which grid is PRODUCTION today
  4. ``data/curation/engine_b_deck_availability.csv``  (RYA-695) the Gerber Engine-B
                                          deck: upstream / staged / RYA-534-ratified
  5. ``data/curation/nlte_grid_availability.csv``      role + wired flags, incl. the
                                          3D and offline-derivation legs

A grid that is present but NOT the registered production grid, and that some other
grid's provenance names in ``supersedes``, is reported SUPERSEDED-BY — the chain is
read from the data, not asserted here.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
GRID_DIR = ROOT / 'data' / 'nlte_grids'
DECK_CSV = ROOT / 'data' / 'curation' / 'engine_b_deck_availability.csv'
AVAIL_CSV = ROOT / 'data' / 'curation' / 'nlte_grid_availability.csv'


class LedgerError(RuntimeError):
    """A source this ledger joins is missing or unreadable (RYA-518: never silent)."""


def _element_of(grid_csv_name: str) -> str:
    """'Ti_Bergemann2011_MPIA.csv' -> 'Ti'. The naming convention is `<El>_<source>`."""
    return grid_csv_name.split('_', 1)[0]


def _engine_a_attempts() -> dict[str, list[dict]]:
    """element -> the Engine-A departure grids BUILT for it, with their fate."""
    from config.constants import NLTE_CORRECTION_ELEMENTS
    production = {el: str(spec.get('grid') or '')
                  for el, spec in NLTE_CORRECTION_ELEMENTS.items()}
    # Fe is production but is NOT in NLTE_CORRECTION_ELEMENTS: its leg is
    # ionization-balance-gated and resolves through `nlte_corrections._MPIA_FE_GRID`
    # instead of the per-element registry. Reading only the registry labels the Fe
    # anchor grid BUILT-NOT-REGISTERED, which is exactly the kind of wrong-but-
    # plausible line that sends someone to re-acquire a grid already in production.
    from pipeline import nlte_corrections as _nc
    production.setdefault('Fe', Path(_nc._MPIA_FE_GRID).name)

    grids, superseded_by = {}, {}
    for p in sorted(GRID_DIR.glob('*.csv')):
        el = _element_of(p.name)
        prov = {}
        for cand in (p.with_suffix('.prov.json'), Path(str(p) + '.prov.json')):
            if cand.exists():
                try:
                    prov = json.loads(cand.read_text(encoding='utf-8'))
                except Exception as exc:                   # noqa: BLE001 — reported
                    raise LedgerError(
                        f'{cand.name} is unreadable ({type(exc).__name__}: {exc}); a '
                        f'grid whose provenance cannot be read must not be reported as '
                        f'though its history were known') from exc
                break
        grids.setdefault(el, []).append((p.name, prov))
        # `supersedes` names the grid this one REPLACED — the attempt history, stated
        # by the replacement rather than inferred from dates.
        sup = str(prov.get('supersedes') or '').strip()
        if sup:
            superseded_by[sup.split()[0]] = (p.name, str(prov.get('ticket') or ''))

    out: dict[str, list[dict]] = {}
    for el, items in grids.items():
        for name, prov in items:
            if name == production.get(el):
                status, detail = 'PRODUCTION', 'registered in NLTE_CORRECTION_ELEMENTS'
            elif name in superseded_by:
                by, tick = superseded_by[name]
                status = 'SUPERSEDED'
                detail = f'replaced by {by}' + (f' ({tick})' if tick else '')
            elif production.get(el):
                # A different grid IS registered for this element, so this one was
                # replaced — but nothing recorded WHEN or WHY (no `supersedes` key
                # anywhere names it). Mg_Bergemann_MPIA, Si_Bergemann_MPIA and
                # Na_Lind2011_INSPECT are all in this state. The displacement is a fact
                # (read from the registry); the reason is genuinely not on record, and
                # saying so is the point — that missing rationale is what gets
                # re-litigated. Never inferred from file dates.
                status = 'REPLACED-IN-PRACTICE'
                detail = (f"not the registered grid for {el} ({production[el]} is); no "
                          f"`supersedes` key anywhere names it, so the REASON for the "
                          f"replacement is not on record — provenance gap, RYA-686 class")
            else:
                # Present, and the element has no registered production grid at all.
                status = 'BUILT-NOT-REGISTERED'
                detail = (f'on disk but {el} has no registered production grid in '
                          f'NLTE_CORRECTION_ELEMENTS; nothing names it in `supersedes`')
            out.setdefault(el, []).append({
                'engine': 'A', 'model': name, 'status': status, 'detail': detail,
                'ticket': str(prov.get('ticket') or ''),
                'reference': str(prov.get('reference') or '')[:120]})
    return out


def _engine_b_attempts() -> dict[str, list[dict]]:
    """element -> the Gerber Engine-B model atom and its staging/ratification state."""
    if not DECK_CSV.exists():
        raise LedgerError(
            f'{DECK_CSV.relative_to(ROOT)} is missing — regenerate it with '
            f'scripts/rya695_engine_b_deck_audit.py (Sirius-only). Without it this '
            f'ledger cannot tell an unstaged Engine-B atom from a nonexistent one, '
            f'which is the exact confusion it was built to end.')
    out: dict[str, list[dict]] = {}
    with DECK_CSV.open(encoding='utf-8') as fh:
        for row in csv.DictReader(fh):
            el = str(row['element']).strip()
            status = {
                'staged-and-ratified': 'PRODUCTION',
                'staged-NOT-ratified': 'STAGED-NOT-RATIFIED',
                'upstream-only-NOT-staged': 'AVAILABLE-NEVER-PULLED',
            }.get(str(row['status']).strip(), str(row['status']).strip())
            detail = str(row.get('rya534_gate') or '').strip()
            if not detail:
                detail = ('published in the Gerber/TSFitPy catalog but never staged on '
                          'Sirius' if status == 'AVAILABLE-NEVER-PULLED' else '')
            blocker = str(row.get('blocker_if_not') or '').strip()
            if blocker:
                detail = (detail + ' — staging does NOT unblock this element: '
                          + blocker[:200]).strip(' —')
            out.setdefault(el, []).append({
                'engine': 'B', 'model': f"gerber:{row['model_atom']}",
                'status': status, 'detail': detail,
                'ticket': 'RYA-534' if status == 'PRODUCTION' else 'RYA-540/695',
                'reference': 'Gerber et al. 2023, A&A 669, A43'})
    return out


def _other_legs() -> dict[str, list[dict]]:
    """element -> non-departure-grid legs on record (3D increments, offline .grd, …).

    These are attempts too: the RYA-399 3D metals increment and the offline PySME
    `.grd` re-derivation inputs are both things that were acquired and wired (or
    deliberately not), and both have been re-investigated more than once.
    """
    if not AVAIL_CSV.exists():
        raise LedgerError(f'{AVAIL_CSV.relative_to(ROOT)} is missing')
    out: dict[str, list[dict]] = {}
    with AVAIL_CSV.open(encoding='utf-8') as fh:
        for row in csv.DictReader(fh):
            role = str(row.get('role') or '').strip()
            if role in ('production',) and str(row.get('subsystem')) == 'registry-nlte':
                continue                       # already covered by _engine_a_attempts
            el = str(row['element']).strip()
            grid = str(row.get('grid_file') or '').strip()
            present = str(row.get('present')).strip().lower() == 'true'
            wired = str(row.get('wired')).strip().lower() == 'true'
            status = ('WIRED' if wired else
                      'ON-DISK-NOT-WIRED' if present else 'NOT-STAGED')
            out.setdefault(el, []).append({
                'engine': str(row.get('subsystem') or ''), 'model': grid,
                'status': status, 'detail': f'role={role}; ' + str(row.get('note') or '')[:150],
                'ticket': '', 'reference': ''})
    return out


def attempts_by_element() -> dict[str, list[dict]]:
    """element -> every model attempt on record, Engine A then Engine B then other legs."""
    merged: dict[str, list[dict]] = {}
    for part in (_engine_a_attempts(), _engine_b_attempts(), _other_legs()):
        for el, items in part.items():
            merged.setdefault(el, []).extend(items)
    return merged


def render_cell(element: str, attempts: dict[str, list[dict]] | None = None) -> str:
    """The `models_tried` tracker cell for one element — compact, one attempt per clause.

    Deliberately terse: the cell answers "has this been tried, and what happened",
    and points at the artifacts for the rest. An element with nothing on record says
    so explicitly, because a blank cell reads as "not investigated" when it may mean
    "investigated and there is nothing".
    """
    attempts = attempts if attempts is not None else attempts_by_element()
    items = attempts.get(element) or []
    if not items:
        return 'none on record'
    parts = []
    for a in items:
        bit = f"{a['model']} [{a['engine']}] {a['status']}"
        if a['status'] not in ('PRODUCTION', 'WIRED') and a['detail']:
            bit += f" ({a['detail'][:110]})"
        elif a['status'] == 'PRODUCTION' and a['engine'] == 'B' and a['detail']:
            bit += f" ({a['detail'][:70]})"
        parts.append(bit)
    return ' ; '.join(parts)
