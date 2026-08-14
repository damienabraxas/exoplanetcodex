#!/usr/bin/env python3
"""
pipeline/element_freeze.py — RYA-814 (child of RYA-812, Lever A)
================================================================
PER-ELEMENT FROZEN VERDICTS. Each element freezes independently with its own
version and hash; the Solar Gold Standard becomes the periodic ASSEMBLY of those
records rather than a single write-once CSV that every element close must
re-freeze.

WHAT THIS FIXES
---------------
Gold is one file, 26 rows, one `CURRENT`, one `hash_manifest.json`, one
whole-artifact `promote_solar_reference --apply`. So closing ANY element
re-freezes and re-gates ALL of them — the monolith tax RYA-812 exists to remove.
After this, freezing Fe cannot touch Ba's record or Ba's hash.

RYA-469 IS DECOUPLED, NOT WEAKENED
----------------------------------
A frozen element record is still WRITE-ONCE and hashed. `freeze_element` refuses
to overwrite an existing (element, version) — the immutability guarantee moves
from "the file never changes" to "each element's record never changes", which is
strictly finer-grained and equally absolute.

WHY EACH RECORD CARRIES A VERBATIM LINE
---------------------------------------
The acceptance gate is sha256 equality between the original monolith and the
re-assembled one. Round-tripping structured fields through a CSV writer will NOT
reproduce the original bytes in general: the notes carry commas, embedded quotes,
em-dashes and unicode arrows, and any writer's quoting/escaping choices are its
own. So each record stores BOTH:

  * `fields`   — the parsed, queryable record (what code should read), and
  * `verbatim` — the exact source line as frozen (what assembly emits).

`verify_consistency()` then asserts the two agree, so `fields` can never silently
drift away from the bytes it claims to describe. Storing only `verbatim` would
make the record opaque; storing only `fields` would make byte-identity a matter of
luck. Both, cross-checked, is the honest construction.

GOLD IS STILL OUTPUT-ONLY
-------------------------
Nothing here reads gold into a validation path. Assembly WRITES a view; RYA-813's
`validate_element` remains the only validator and it reads literature, never this.
"""
from __future__ import annotations

import hashlib
import json
import os
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional

ROOT = Path(__file__).resolve().parents[1]
SOLAR_DIR = ROOT / "data" / "reference" / "solar"
ELEMENTS_DIR = SOLAR_DIR / "elements"
PREAMBLE_NAME = "_preamble.json"
CURRENT_NAME = "_current.json"

#: RYA-812 principle 6 — the gold-complete bar: all 28 species VIS-validated.
#:
#: THE 26/28 GAP IS RESOLVED AND IT IS NOT AN ELEMENT-COUNT QUIRK. RYA-109 ratifies
#: the canonical 27 with **Fe I and Fe II COUNTED SEPARATELY** — Fe II is entry #2 at
#: priority 1, because it constrains log g via ionisation equilibrium while Fe I
#: constrains Teff ("must be tracked separately in the pipeline"). RYA-757 then adds
#: Zn as the 28th. So:
#:
#:      27 canonical (incl. Fe II)  +  Zn  =  28
#:      gold v4 carries 26          =  27 - Fe II
#:      missing: Fe II  and  Zn
#:
#: Fe II is MEASURED and ADOPTED already (arbiter 7.486-7.500, RYA-406/714 §3.2); it
#: simply has no frozen gold ROW yet.
#:
#: ⚠️ THEREFORE THE RECORD KEY IS (element, ion), NOT element. Keying by element
#: alone cannot represent Fe I and Fe II as distinct frozen records — they would
#: collide on one file and one would silently overwrite the other, which is the exact
#: class of defect RYA-469 exists to prevent. This is not hypothetical: gold v4
#: already carries **Sc II**, so the schema is ion-bearing today and Fe II is queued.
GOLD_COMPLETE_COUNT = 28


class ElementFreezeError(RuntimeError):
    """Write-once violated, or a record is internally inconsistent."""


def _sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


@dataclass(frozen=True)
class FrozenElement:
    element: str
    version: str
    fields: dict[str, Any]
    verbatim: str
    sha256: str
    frozen_utc: str
    source_commit: str
    provenance: str

    @property
    def ion(self) -> str:
        return str(self.fields.get("ion", "") or "")

    @property
    def key(self) -> str:
        return record_key(self.element, self.ion)

    @property
    def row_index(self) -> int:
        """Position in the assembled table. Preserved so assembly is byte-stable."""
        return int(self.fields.get("_row_index", -1))


def record_key(element: str, ion: str = "") -> str:
    """
    The per-record key. (element, ion) — NEVER element alone.

    RYA-109 counts Fe I and Fe II as separate canonical entries, and gold already
    carries Sc II, so an element-only key would collide two distinct frozen species
    onto one file. Ion is normalised and optional so records written before this fix
    (element-only names) still resolve.
    """
    ion = (ion or "").strip()
    return f"{element}_{ion}" if ion else element


def split_key(key: str) -> tuple[str, str]:
    """'Fe_II' -> ('Fe', 'II'). An un-suffixed legacy key returns ('Fe', '')."""
    el, _, ion = key.partition("_")
    return el, ion


def element_path(element: str, ion: str = "") -> Path:
    return ELEMENTS_DIR / f"{record_key(element, ion)}.json"


def preamble_path() -> Path:
    return ELEMENTS_DIR / PREAMBLE_NAME


def current_path() -> Path:
    return ELEMENTS_DIR / CURRENT_NAME


# ── serialisation that must reproduce the original bytes ─────────────────────
def serialise_fields(fields: dict[str, Any], columns: list[str]) -> str:
    """
    Re-serialise parsed fields to a CSV line, matching the frozen file's dialect.

    Deliberately hand-rolled rather than csv.writer: the target dialect quotes a
    field ONLY when it contains a comma, a quote or a newline, and doubles inner
    quotes. That is what the existing gold files do, and the round-trip gate is
    byte-level — a writer that quotes more (or less) eagerly fails it for reasons
    that have nothing to do with the science.
    """
    out = []
    for col in columns:
        v = fields.get(col, "")
        s = "" if v is None else str(v)
        if any(c in s for c in (",", '"', "\n")):
            s = '"' + s.replace('"', '""') + '"'
        out.append(s)
    return ",".join(out)


def verify_consistency(rec: FrozenElement, columns: list[str]) -> None:
    """`fields` must reproduce `verbatim`, or the record lies about its own bytes."""
    rebuilt = serialise_fields(rec.fields, columns)
    if rebuilt != rec.verbatim:
        raise ElementFreezeError(
            f"{rec.element}: parsed fields do not re-serialise to the frozen line.\n"
            f"  frozen  : {rec.verbatim[:160]}\n"
            f"  rebuilt : {rebuilt[:160]}\n"
            f"A record whose structured view disagrees with its bytes cannot be "
            f"trusted for either purpose.")


# ── read ─────────────────────────────────────────────────────────────────────
def read_preamble() -> dict[str, Any]:
    p = preamble_path()
    if not p.exists():
        raise ElementFreezeError(
            f"no preamble at {p} — the assembled file's comment header and column "
            f"order are part of its bytes and must be frozen alongside the rows.")
    return json.loads(p.read_text(encoding="utf-8"))


def read_element(element: str, ion: str = "") -> Optional[FrozenElement]:
    p = element_path(element, ion)
    if not p.exists() and not ion:
        # tolerate pre-fix element-only records
        cands = sorted(ELEMENTS_DIR.glob(f"{element}_*.json")) if ELEMENTS_DIR.is_dir() else []
        if len(cands) == 1:
            p = cands[0]
        elif len(cands) > 1:
            raise ElementFreezeError(
                f"{element} has {len(cands)} ion records "
                f"({', '.join(c.stem for c in cands)}) — read_element must be given "
                f"an ion to disambiguate. Fe I and Fe II are different species "
                f"(RYA-109), not one element read twice.")
    if not p.exists():
        return None
    d = json.loads(p.read_text(encoding="utf-8"))
    return FrozenElement(
        element=d["element"], version=d["version"], fields=d["fields"],
        verbatim=d["verbatim"], sha256=d["sha256"], frozen_utc=d["frozen_utc"],
        source_commit=d.get("source_commit", ""), provenance=d.get("provenance", ""))


def read_record(key: str) -> Optional["FrozenElement"]:
    """Read by the KEY that `frozen_elements()` returns (e.g. 'Fe_II')."""
    return read_element(*split_key(key))


def frozen_elements() -> list[str]:
    """Record KEYS (element_ion), one per frozen species."""
    if not ELEMENTS_DIR.is_dir():
        return []
    return sorted(p.stem for p in ELEMENTS_DIR.glob("*.json")
                  if not p.name.startswith("_"))


def read_current_map() -> dict[str, str]:
    p = current_path()
    if not p.exists():
        return {}
    return json.loads(p.read_text(encoding="utf-8"))


# ── write (write-once, per element) ──────────────────────────────────────────
def freeze_element(element: str, fields: dict[str, Any], verbatim: str, *,
                   version: str, source_commit: str = "", provenance: str = "",
                   columns: Optional[list[str]] = None,
                   allow_new_version: bool = False,
                   ion: Optional[str] = None) -> FrozenElement:
    """
    Freeze ONE element. Refuses to overwrite an existing record for the same
    version — RYA-469 immutability, decoupled to per-element granularity.

    `allow_new_version=True` permits writing a record for a DIFFERENT version of
    the same element (the normal re-freeze path). Same element + same version is
    always a refusal: that is the case where a value could silently change.
    """
    ELEMENTS_DIR.mkdir(parents=True, exist_ok=True)
    # `frozen_elements()` returns KEYS ('Fe_II'); these take an ELEMENT ('Fe').
    # Passing a key here silently produced 'Fe_II_II' -- a spurious extra record
    # that inflated the species count. Refuse it instead of accepting both shapes.
    if "_" in element:
        raise ElementFreezeError(
            f"freeze_element takes an ELEMENT ({element.split('_')[0]!r}) plus an "
            f"ion, not a record key ({element!r}). Use split_key() first, or "
            f"read_record() if you meant to read.")
    ion = ion if ion is not None else str(fields.get("ion", "") or "")
    key = record_key(element, ion)
    existing = read_element(element, ion)
    if existing is not None:
        if existing.version == version:
            raise ElementFreezeError(
                f"{key} is already frozen at {version} "
                f"(sha256 {existing.sha256[:12]}). A frozen element record is "
                f"WRITE-ONCE (RYA-469, decoupled per-element by RYA-814). To "
                f"supersede it, freeze a NEW version.")
        if not allow_new_version:
            raise ElementFreezeError(
                f"{key} is frozen at {existing.version}; refusing to write "
                f"{version} without allow_new_version=True.")

    rec = FrozenElement(
        element=element, version=version, fields=fields, verbatim=verbatim,
        sha256=_sha256_text(verbatim),
        frozen_utc=datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        source_commit=source_commit, provenance=provenance)
    if columns:
        verify_consistency(rec, columns)

    element_path(element, ion).write_text(
        json.dumps({
            "element": rec.element, "version": rec.version,
            "frozen_utc": rec.frozen_utc, "source_commit": rec.source_commit,
            "provenance": rec.provenance, "sha256": rec.sha256,
            "verbatim": rec.verbatim, "fields": rec.fields,
        }, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    cur = read_current_map()
    cur[key] = version
    current_path().write_text(
        json.dumps(dict(sorted(cur.items())), indent=2) + "\n", encoding="utf-8")
    return rec


# ── assembly (the periodic snapshot) ─────────────────────────────────────────
def assemble(version_label: Optional[str] = None) -> str:
    """
    Assemble the gold table from the per-element frozen records.

    Emits the frozen preamble, the frozen column header, then each element's
    VERBATIM line in its frozen row order. Byte-identical to the monolith it was
    exploded from — which is exactly what the RYA-814 round-trip gate asserts.
    """
    pre = read_preamble()
    els = frozen_elements()
    if not els:
        raise ElementFreezeError("no frozen element records to assemble")

    recs = []
    for key in els:
        p = ELEMENTS_DIR / f"{key}.json"
        d = json.loads(p.read_text(encoding="utf-8"))
        r = FrozenElement(
            element=d["element"], version=d["version"], fields=d["fields"],
            verbatim=d["verbatim"], sha256=d["sha256"], frozen_utc=d["frozen_utc"],
            source_commit=d.get("source_commit", ""),
            provenance=d.get("provenance", ""))
        verify_consistency(r, pre["columns"])
        if _sha256_text(r.verbatim) != r.sha256:
            raise ElementFreezeError(
                f"{key}: verbatim line does not match its recorded sha256 — the "
                f"frozen record has been altered on disk.")
        recs.append(r)
    recs.sort(key=lambda r: r.row_index)

    lines = list(pre["comment_lines"])
    if version_label:
        lines = [f"# version: {version_label}" if l.startswith("# version:") else l
                 for l in lines]
    lines.append(",".join(pre["columns"]))
    lines.extend(r.verbatim for r in recs)
    return "\n".join(lines) + ("\n" if pre.get("trailing_newline", True) else "")


def assembly_status() -> dict[str, Any]:
    """Is the assembly gold-complete (28) or partial? Partial must SAY so."""
    els = frozen_elements()
    n = len(els)
    return {
        "elements": n,
        "note": ("count is of SPECIES (element x ion): RYA-109 counts Fe I and Fe II "
                 "separately, and gold already carries Sc II."),
        "required_for_gold_complete": GOLD_COMPLETE_COUNT,
        "gold_complete": n >= GOLD_COMPLETE_COUNT,
        "missing_count": max(0, GOLD_COMPLETE_COUNT - n),
        "label": ("gold (complete, 28 species)" if n >= GOLD_COMPLETE_COUNT
                  else f"PARTIAL assembly — {n}/{GOLD_COMPLETE_COUNT} species. "
                       f"Below the RYA-812 principle-6 bar, so this is NOT 'gold'; "
                       f"it is a labelled partial snapshot. Missing: Fe II (measured "
                       f"and adopted as the RYA-406 arbiter but never frozen as a gold "
                       f"row) and Zn (RYA-757 intake, not yet measured)."),
        "frozen": els,
    }
