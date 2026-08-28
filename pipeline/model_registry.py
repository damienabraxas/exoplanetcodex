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


class ModelRegistryError(LookupError):
    """A token could not be bound to exactly one model, SAFELY."""


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
