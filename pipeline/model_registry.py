"""
The canonical roster of abundance-software models — RYA-1101.
=============================================================
`data/catalog/model_registry.csv` is the stable list of the models the pipeline can run.
It exists because the roster previously lived only in chat memory: the nicknames
"Frankenstein" (models 5+6) and "The Bride" (model 8) were in no ledger, and the roster
got re-litigated every session.

WHAT THIS FILE IS NOT. It carries no live product counts — those are volatile and owned by
RYA-1015's element x model matrix. Duplicating them here would create a second source of
truth for a number that changes every run, which is the defect this registry exists to end.

🔴 THE BARE `ENGINE-B` COLLISION. `ENGINE-B` is a live token, and it does not mean what it
used to. Read from `treatment_axes.LEGACY`:

    "ENGINE-B":      scale=1D-LTE,  model=none        -> model 1
    "ENGINE-B-NLTE": scale=1D-NLTE, model=gerber      -> model 4

Historically the bare string meant the Gerber NLTE synthesis, i.e. model 4. After RYA-906's
physics-axis renaming it resolves to model 1. **Same string, two physics, and nothing in the
old artifacts says which one was meant.** `resolve()` therefore REFUSES the bare token rather
than binding it, and names both candidates. The canonical binding is Ryan's decision, not
this module's — see `docs/catalog/model_registry_notes.md`.
"""
from __future__ import annotations

import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REGISTRY = ROOT / "data" / "catalog" / "model_registry.csv"

#: Tokens that are live in the code but are NOT a model in their own right. Listing them
#: here is what makes `resolve()` able to refuse instead of guessing.
AMBIGUOUS_TOKENS = {
    "ENGINE-B": (
        "bare 'ENGINE-B' is ambiguous and will not be bound. It resolves to model 1 "
        "(1D-LTE) under treatment_axes.LEGACY, but HISTORICALLY meant model 4 "
        "(1D-NLTE Gerber, now 'ENGINE-B-NLTE'). Same string, two physics — an artifact "
        "written before RYA-906 does not say which. Use the canonical token for the model "
        "you mean: '1D-LTE' for model 1, 'ENGINE-B-NLTE' for model 4. The canonical "
        "binding for legacy data is Ryan's call (RYA-1101); see "
        "docs/catalog/model_registry_notes.md."
    ),
}

#: The roster's status vocabulary. `not-emitted` is deliberately NOT folded into `in-dev`:
#: the Bride is being BUILT, while Frankenstein's Dog is a leg the science requires that
#: nothing currently emits. Collapsing them would hide which of the two is waiting on a
#: solver and which is waiting on a wiring decision — the same one-token-two-meanings
#: defect the bare `ENGINE-B` alias below exists to flag.
STATUSES = ("live", "in-dev", "not-emitted")

#: 🔴 THE `line_set` HOOK — RYA-1111 WIRES THE VALUES, RYA-1101 ONLY OPENS THE COLUMN.
#: Which pool of lines a measurement was made on is a PROVENANCE AXIS, not a property of
#: the model: any model here can in principle be measured on any of these sets, so every
#: roster row carries `-` ("not model-scoped") until RYA-1111 has products to key.
#:
#: The vocabulary is declared here so the column cannot accept a value typed from memory
#: the day it is first populated — `check_line_sets` refuses anything outside it. That is
#: the whole point of opening the column early rather than letting 1111 invent one.
#:
#: ⚠️ `consistent` IS ABSENT ON PURPOSE. RYA-1105 removes the Consistent tier from the
#: active pipeline and the website; the going-forward set is Asplund Grade / Our Grade /
#: Deep Grade. `--lines-tier consistent` is STILL LIVE in `derive_band_products.py` — that
#: is RYA-1105's to retire, not this ticket's (RYA-1101 forbids touching model code), so
#: the code and this vocabulary genuinely disagree today and the mismatch is recorded
#: rather than papered over.
#: ⚠️ `asplund-graded` WAS WRONG AND IS CORRECTED TO `asplund` (RYA-1111). RYA-1101 opened
#: this column before the axis had an owner, and guessed the spelling. RYA-1111 owns the
#: axis and names the values {asplund, gbs, our-graded, our-deep-graded}; that wins. This
#: tuple is now THE ONE DEFINITION -- `pipeline.reference_lineset` imports it rather than
#: restating it, so the registry's guard and the ingest loader cannot drift into
#: disagreeing about what a valid value is.
#:
#: ⚠️ A reference FILE may carry its own native spelling (the RYA-1109 artifact's column
#: says `asplund_agss21`). That is recorded by the adapter as `native_line_set` and mapped
#: to the canonical name -- never silently rewritten in the file.
LINE_SETS = (
    "-",                 # not model-scoped; the value every roster row carries today
    "asplund",           # the imported AGSS21 reference set (RYA-1109)
    "gbs",               # the Gaia FGK Benchmark Stars reference set (RYA-1110)
    "our-graded",        # our lab-gf graded pool, at or below the depth gate
    "our-deep-graded",   # the saturated population above the gate (RYA-984/954)
    # 🔴 RYA-1127 — ADDED BECAUSE THE KEY NOW NEEDS THEM, NOT SPECULATIVELY. Putting
    # `line_set` into the product identity key means EVERY product must resolve one, and
    # nine records in `quarantine[]` did not: six carry tier UNGRADED and three carry
    # tier ALL. Both are real `--lines-tier` values (derive_band_products.py:2102) naming
    # a real pool of OUR lines -- the non-lab-gf split, and the whole pool undivided -- so
    # the vocabulary was simply incomplete rather than these being unidentifiable.
    # Widening it here is the deliberate act `line_set_for_product` demands; the
    # alternative was to let `key_of` default, which is the RYA-869 class exactly.
    "our-ungraded",      # our lines BELOW the gf-grade bar (`--lines-tier ungraded`)
    "our-all",           # our whole pool, graded and ungraded together (`--lines-tier all`)
    # ⚠️ `consistent` REMAINS ABSENT, deliberately -- RYA-1105 retires that tier, and a
    # product carrying it must still fail loudly rather than acquire a name here.
)


class ModelRegistryError(LookupError):
    """A token could not be bound to exactly one model, SAFELY."""


def check_line_sets(rows: list[dict]) -> list[str]:
    """Every `line_set` value is in `LINE_SETS`, or say which are not.

    Returns the complaints rather than raising, so the verifier can report every bad row
    in one pass instead of stopping at the first.
    """
    bad = []
    for r in rows:
        v = (r.get("line_set") or "").strip()
        if v not in LINE_SETS:
            bad.append(f"model {r['model_id']} line_set {v!r} is not in LINE_SETS "
                       f"{LINE_SETS} — RYA-1111 must add the value to the vocabulary "
                       f"before writing it to the roster")
    return bad


def load(path: Path | None = None) -> list[dict]:
    with (path or REGISTRY).open(newline="") as fh:
        return list(csv.DictReader(fh))


def resolve(token: str, path: Path | None = None) -> dict:
    """The one model this token names, or RAISE.

    Refuses rather than guessing — for the bare `ENGINE-B` alias, and for any token that is
    absent or maps to more than one row.
    """
    tok = str(token).strip()
    if tok in AMBIGUOUS_TOKENS:
        raise ModelRegistryError(AMBIGUOUS_TOKENS[tok])
    hits = [r for r in load(path) if r["stored_token"].strip() == tok and tok]
    if len(hits) == 1:
        return hits[0]
    if not hits:
        raise ModelRegistryError(
            f"{tok!r} names no model in {REGISTRY.name}. If it is a real treatment, add it "
            f"to the registry with its source ticket — do not bind it here from memory.")
    raise ModelRegistryError(
        f"{tok!r} maps to {len(hits)} models "
        f"({', '.join(h['model_id'] for h in hits)}) — an undocumented collision.")
