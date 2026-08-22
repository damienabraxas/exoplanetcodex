"""
RYA-964 — turn a raw data label into a star id at intake. One alias lookup, nothing more.

WHY THIS EXISTS. Our catalogues are already clean — `system_catalog`, `stars.yaml` and the
holdings registry all key on `alpha_cen_a`. The mess is entirely at the raw-data layer, where
whatever the telescope or scheduler typed comes through unfiltered: `alf Cen A`, `alf_Cen_A`,
`alpha-cen-a`, `HD128620`, `ALPHA-CEN-A`. RYA-952 found one star arriving under four different
`OBJECT` strings in the same archive (`eps Eri` / `HD 22049` / `Epsilon-Eridani` / `EPSERI`),
and 236 HARPS files split across them. This is the single gate a raw string passes through.

🔴 IT REFUSES RATHER THAN GUESSES, AND THAT IS THE WHOLE VALUE.
`STD`, `CAL_*`, `Star S5` are not aliases and never will be. RYA-952 found four tau Ceti
science frames labelled `OBJECT='STD'` — the pipeline wrote the ROLE where the name goes — and
a junk `Star S5` frame filed under Alpha Cen B. A resolver that tried to be clever about those
would have assigned them somewhere. `UNRESOLVED` sends them to a human, who adds one alias and
makes the label known forever.

⚠️ WHAT THIS IS NOT. It is not an identity CHECK. It reads a label and trusts that the label,
once recognised, is honest. Deciding whether a frame really points at the star it claims is a
different problem with a different method — astrometry (`pipeline.audit_crires`) or, for the
Alpha Cen A/B pair whose components no label can separate, RV/template physics (RYA-423).
Those stay separate and are not wired in here.

Aliases live ONLY in `data/catalog/system_catalog.csv`. Nothing in this module hardcodes a
name map, so adding a spelling is a one-line CSV edit and never a code change.
"""
from __future__ import annotations

import csv
from functools import lru_cache
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CATALOG = ROOT / "data" / "catalog" / "system_catalog.csv"

#: What `resolve_star` returns when a label matches nothing. A sentinel STRING rather than
#: None so it survives a CSV round-trip and shows up in an audit table as itself.
UNRESOLVED = "UNRESOLVED"


class StarLabelUnresolved(KeyError):
    """Raised by `resolve_star(..., strict=True)` for a label with no alias."""


def normalize_label(raw) -> str:
    """Collapse a raw label to its comparison form: lowercase alphanumerics only.

    Aggressive on purpose. `alf Cen A`, `alf_Cen_A`, `alpha-cen-a`, `ALPHA-CEN-A` and
    `HD 128620` differ only in separators and case, so folding those away means the CSV
    needs a handful of real spellings instead of the cartesian product of punctuation.
    It also folds SIMBAD's catalogue decoration — `* tau Cet`, `** STT 270A` — onto the
    same form, which is what lets the alias lists be pasted straight from SIMBAD.
    """
    return "".join(ch for ch in str(raw).lower() if ch.isalnum())


@lru_cache(maxsize=1)
def _alias_index(path: str = None) -> dict[str, str]:
    """{normalized alias -> system_id}, built from the catalogue and nowhere else."""
    p = Path(path) if path else CATALOG
    index: dict[str, str] = {}
    collisions: list[str] = []
    with open(p, encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            sid = (row.get("star_params_key") or "").strip()
            if not sid:
                continue          # a system we hold no parameters for cannot be a target yet
            names = [sid, row.get("system_name", "")]
            names += [a for a in (row.get("aliases") or "").split("|") if a.strip()]
            for n in names:
                k = normalize_label(n)
                if not k:
                    continue
                if k in index and index[k] != sid:
                    # 🔴 TWO STARS CLAIMING ONE LABEL IS A CATALOGUE BUG, NOT A TIE TO BREAK.
                    collisions.append(f"{k!r} claimed by both {index[k]} and {sid}")
                index[k] = sid
    if collisions:
        raise ValueError(
            "system_catalog.csv has ambiguous aliases — one label cannot mean two stars:\n  "
            + "\n  ".join(collisions))
    return index


def resolve_star(raw_label, *, strict: bool = False) -> str:
    """Raw label -> system_id, or `UNRESOLVED`.

    `strict=True` raises `StarLabelUnresolved` instead, for callers that would rather stop
    than carry a sentinel forward.
    """
    key = normalize_label(raw_label)
    sid = _alias_index().get(key, UNRESOLVED)
    if sid is UNRESOLVED or sid == UNRESOLVED:
        if strict:
            raise StarLabelUnresolved(
                f"no star matches the label {raw_label!r} (normalised {key!r}). If this is a "
                f"real target, add the spelling to the `aliases` column of "
                f"data/catalog/system_catalog.csv — one edit and it is known forever. If it "
                f"is a role or a placeholder (STD, CAL_*, Star S5), it SHOULD land here: "
                f"quarantine it.")
        return UNRESOLVED
    return sid


def known_aliases() -> dict[str, list[str]]:
    """{system_id: [normalized aliases]} — for reporting and tests."""
    out: dict[str, list[str]] = {}
    for alias, sid in sorted(_alias_index().items()):
        out.setdefault(sid, []).append(alias)
    return out
