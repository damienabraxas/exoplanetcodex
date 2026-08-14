#!/usr/bin/env python3
"""
pipeline/litscan.py — RYA-813
=============================
THE LITERATURE SCAN, MACHINE-READABLE. The pass band for VIS validation comes from
HERE, per element, and never from a flat tolerance (RYA-812 principle 3: "O's
spread != Fe's spread").

WHAT THIS IS NOT
----------------
It is not gold. It is not a Codex measurement. Every value in a litscan is an
OUTSIDE number — a published determination with a citation — which is exactly what
makes it safe to validate against: `validate_element` compares our measurement to
the world, so the gold -> verdict -> gold loop cannot form (principle 5).

MISSING IS A FIRST-CLASS ANSWER
-------------------------------
`literature_range(el)` returns None when an element has no usable range, and the
caller MUST treat that as `un-anchorable` — report-only across all bands, flagged
loudly (principle 2). It is an honest gap, not a failure, and never a silent skip
(the RYA-786 class). A missing litscan therefore raises nothing here; it is data.

SCALE IS PART OF THE COMPARISON
-------------------------------
A range carries the scale it is quoted on. Comparing a 1D-NLTE measurement against
a 3D-NLTE range is a scale error dressed as a disagreement, so `LiteratureRange`
exposes `scale` and the validator is required to check it. This mirrors
`pipeline.solar_scale_provenance`, which owns the classification itself.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional

import yaml

LITSCAN_DIR = Path(__file__).resolve().parents[1] / "data" / "reference" / "litscan"


class LitscanError(RuntimeError):
    """A litscan exists but is malformed. Absent is NOT an error (see module docstring)."""


@dataclass(frozen=True)
class LiteratureEntry:
    source: str
    citation: str
    value: Optional[float]
    sigma: Optional[float]
    scale: Optional[str]
    role: str
    note: str = ""

    @property
    def contributes_value(self) -> bool:
        return self.value is not None


@dataclass(frozen=True)
class LiteratureRange:
    element: str
    central: float
    min: float
    max: float
    scale: Optional[str]
    basis: str
    derivation: str
    is_statistical_spread: bool
    sigma_external: Optional[float] = None
    deviate_beyond: Optional[float] = None
    best_external: str = ""
    entries: tuple[LiteratureEntry, ...] = field(default_factory=tuple)
    verification_owed: tuple[str, ...] = field(default_factory=tuple)

    def contains(self, value: float) -> bool:
        return self.min <= value <= self.max

    def offset(self, value: float) -> float:
        """Signed distance from the central literature value."""
        return value - self.central

    def is_deviate(self, value: float) -> bool:
        """
        RYA-714 §4 ratifies `DEVIATE (>2σ)` as a distinct code from a near miss.
        A value beyond 2σ_ext is a different claim from one 1.2σ out, and the
        dossier says so; collapsing them would throw away the distinction the
        litscan exists to make.
        """
        if self.deviate_beyond is None:
            return False
        return abs(self.offset(value)) > self.deviate_beyond

    def excess(self, value: float) -> float:
        """How far OUTSIDE the band, 0.0 if inside. The number an exception must justify."""
        if value < self.min:
            return self.min - value
        if value > self.max:
            return value - self.max
        return 0.0

    @property
    def citations(self) -> tuple[str, ...]:
        return tuple(f"{e.source} ({e.citation})" for e in self.entries if e.citation)


def litscan_path(element: str) -> Path:
    return LITSCAN_DIR / f"{element}.yaml"


def has_litscan(element: str) -> bool:
    return litscan_path(element).exists()


def load_litscan(element: str) -> Optional[dict[str, Any]]:
    """Raw litscan mapping, or None if this element has no litscan yet."""
    p = litscan_path(element)
    if not p.exists():
        return None
    try:
        with p.open() as fh:
            doc = yaml.safe_load(fh)
    except yaml.YAMLError as exc:
        raise LitscanError(f"litscan for {element} is not valid YAML: {exc}") from exc
    if not isinstance(doc, dict):
        raise LitscanError(f"litscan for {element} is not a mapping")
    return doc


def literature_range(element: str) -> Optional[LiteratureRange]:
    """
    The element's literature range, or None if it has no usable one.

    None means UN-ANCHORABLE and the caller must say so out loud. A litscan that
    exists but declares no `range` is also None — a scan can be started (entries
    collected) before it can define a band, and pretending otherwise would
    manufacture a pass band out of an incomplete survey.
    """
    doc = load_litscan(element)
    if doc is None:
        return None
    rng = doc.get("range")
    if not rng:
        return None
    for key in ("central", "min", "max"):
        if rng.get(key) is None:
            raise LitscanError(
                f"litscan for {element} declares a `range` but no `{key}`. An "
                f"incomplete band cannot be used as a pass criterion — either "
                f"complete it or remove `range` so the element reads as "
                f"un-anchorable.")
    lo, hi, mid = float(rng["min"]), float(rng["max"]), float(rng["central"])
    if not (lo <= mid <= hi):
        raise LitscanError(
            f"litscan for {element} has central {mid} outside [{lo}, {hi}]")

    entries = tuple(
        LiteratureEntry(
            source=str(e.get("source", "?")),
            citation=str(e.get("citation", "")),
            value=(None if e.get("value") is None else float(e["value"])),
            sigma=(None if e.get("sigma") is None else float(e["sigma"])),
            scale=(None if e.get("scale") is None else str(e["scale"])),
            role=str(e.get("role", "")),
            note=str(e.get("note", "")),
        )
        for e in (doc.get("entries") or [])
    )
    return LiteratureRange(
        element=element,
        central=mid, min=lo, max=hi,
        scale=(None if doc.get("scale") is None else str(doc["scale"])),
        basis=str(rng.get("basis", "")),
        derivation=str(rng.get("derivation", "unspecified")),
        is_statistical_spread=bool(rng.get("is_statistical_spread", False)),
        sigma_external=(None if rng.get("sigma_external") is None
                        else float(rng["sigma_external"])),
        deviate_beyond=(None if rng.get("deviate_beyond") is None
                        else float(rng["deviate_beyond"])),
        best_external=str(rng.get("best_external", "")),
        entries=entries,
        verification_owed=tuple(str(x) for x in (doc.get("verification_owed") or [])),
    )


def available_elements() -> tuple[str, ...]:
    if not LITSCAN_DIR.is_dir():
        return ()
    return tuple(sorted(p.stem for p in LITSCAN_DIR.glob("*.yaml")))
