#!/usr/bin/env python3
"""
pipeline/corrections_registry.py — RYA-674 §2A/§2B
==================================================
THE DICTIONARY FOR `corrections_applied`.

A frozen gold row records WHICH tabulated corrections its number already carries, as a
JSON list of correction identifiers in the `corrections_applied` column. This module
is what those identifiers MEAN: it loads `config/corrections_registry.yaml`, resolves
every `*_source` pointer into `config/constants.py`, and derives the quantities the
guards need.

Why pointers and not values
---------------------------
RYA-674 was opened against one defect class: *metadata asserted in one place about a
quantity that lives in another, with nothing keeping them in sync*. The ticket's own
§2B specified `pre_range: [7.50, 7.53]` / `post_range: [7.45, 7.48]` as literals in the
YAML — three hand-maintained copies of the Magic-2013 offset (the magnitude, the pre
band, the post band) that must all be edited together and none of which fails if only
one is. That is the defect, one level up.

So the YAML carries pointers; this module resolves them; and the "does this value look
already corrected?" bands are DERIVED from the magnitude via the RYA-681 construction
(`pipeline.solar_scale_provenance`): the two scale centres are separated by exactly
|magnitude|, so the two-hypothesis decision boundary is their midpoint and the
half-width is |magnitude| / 2. Revise the tabulated offset in `config/constants.py` and
every band moves with it, with no edit here.

What §2B asked for, and what it gets
------------------------------------
§2B's three-way outcome is implemented in
`pipeline.solar_scale_provenance.apply_reported_scale_correction`:

  value looks corrected  + record says applied      -> skip cleanly (idempotent)
  value looks corrected  + record says NOT applied  -> ScaleProvenanceError (loud)
  value looks pre-correction                        -> apply, and record the identifier

`value_check_required()` implements §2B's trigger condition: a correction whose
magnitude is at least half the acceptance gate's half-width cannot be caught by the
gate alone, so the guard MUST classify the value as well. It is COMPUTED from the two
pointers rather than asserted in prose.
"""
from __future__ import annotations

import re
from functools import lru_cache
from pathlib import Path
from typing import Any

import yaml

__all__ = [
    'CorrectionsRegistryError',
    'REGISTRY_PATH',
    'load_registry',
    'correction_ids',
    'correction',
    'corrections_for_element',
    'magnitude',
    'gate_half_width',
    'value_check_required',
    'scale_bands',
    'assert_registry_consistent',
]

REGISTRY_PATH = Path(__file__).resolve().parents[1] / 'config' / 'corrections_registry.yaml'

#: Only `config.constants` may be pointed at. A registry that could name an arbitrary
#: module would be an import side-channel, and every quantity a correction needs is a
#: declared constant by construction (RYA-659 single-sourcing).
_ALLOWED_MODULES = ('config.constants',)
_POINTER_RE = re.compile(r"^(?P<module>[A-Za-z_][\w.]*)\.(?P<attr>[A-Za-z_]\w*)"
                         r"(?:\[(?P<key>'[^']*'|\"[^\"]*\")\])?$")

_REQUIRED_KEYS = ('provenance_ticket', 'citation', 'scope', 'magnitude_source',
                  'pre_scale', 'post_scale')


class CorrectionsRegistryError(RuntimeError):
    """The corrections registry is missing, malformed, or names a quantity that does
    not exist. Raised instead of falling back — an unresolvable correction identifier
    is exactly how a correction gets applied twice (RYA-669)."""


def _resolve_pointer(pointer: str, *, what: str) -> Any:
    """Resolve a `config.constants.NAME` / `config.constants.NAME['key']` pointer.

    Deliberately not `eval`: a registry entry may name a declared constant and nothing
    else. Failure is loud — a dangling pointer means the registry and the constants
    have drifted, which is the very desync this file exists to prevent.
    """
    m = _POINTER_RE.match(str(pointer).strip())
    if not m:
        raise CorrectionsRegistryError(
            f"{what}: {pointer!r} is not a constants pointer of the form "
            f"'config.constants.NAME' or \"config.constants.NAME['key']\"")
    mod_name, attr, key = m.group('module'), m.group('attr'), m.group('key')
    if mod_name not in _ALLOWED_MODULES:
        raise CorrectionsRegistryError(
            f"{what}: {pointer!r} points at {mod_name!r}; only {_ALLOWED_MODULES} may be "
            f"named by the corrections registry")
    import importlib
    module = importlib.import_module(mod_name)
    if not hasattr(module, attr):
        raise CorrectionsRegistryError(f"{what}: {mod_name} has no attribute {attr!r}")
    value = getattr(module, attr)
    if key is not None:
        k = key[1:-1]
        try:
            value = value[k]
        except (KeyError, TypeError) as exc:
            raise CorrectionsRegistryError(
                f"{what}: {mod_name}.{attr} has no entry {k!r} ({exc})") from exc
    return value


@lru_cache(maxsize=1)
def load_registry() -> dict:
    """The parsed registry, keyed by correction identifier. Loud on absence."""
    if not REGISTRY_PATH.exists():
        raise CorrectionsRegistryError(
            f"corrections registry missing at {REGISTRY_PATH} — a correction whose "
            f"identifier cannot be resolved may not be recorded as applied (RYA-674)")
    data = yaml.safe_load(REGISTRY_PATH.read_text(encoding='utf-8'))
    if not isinstance(data, dict) or not data:
        raise CorrectionsRegistryError(f"{REGISTRY_PATH} did not parse to a non-empty mapping")
    return data


def correction_ids() -> tuple[str, ...]:
    return tuple(sorted(load_registry()))


def correction(correction_id: str) -> dict:
    """One registry entry. Raises on an unknown identifier — never returns a default."""
    reg = load_registry()
    if correction_id not in reg:
        raise CorrectionsRegistryError(
            f"unknown correction identifier {correction_id!r}; the registry declares "
            f"{list(reg)}. An identifier recorded in `corrections_applied` that this "
            f"registry cannot explain is an undeclared correction (RYA-674).")
    return reg[correction_id]


def corrections_for_element(element: str) -> tuple[str, ...]:
    """Every registered correction whose scope covers `element`, sorted."""
    out = []
    for cid, entry in load_registry().items():
        els = ((entry.get('scope') or {}).get('elements') or [])
        if element in els:
            out.append(cid)
    return tuple(sorted(out))


def magnitude(correction_id: str) -> float:
    """The tabulated magnitude in dex, resolved from `config/constants.py`."""
    entry = correction(correction_id)
    return float(_resolve_pointer(entry['magnitude_source'],
                                  what=f"{correction_id}.magnitude_source"))


def gate_half_width(correction_id: str) -> float | None:
    """Half-width of the acceptance gate this correction sits inside, or None if the
    entry declares no gate. Derived from the two gate-edge pointers."""
    entry = correction(correction_id)
    lo_p, hi_p = entry.get('gate_lower_source'), entry.get('gate_upper_source')
    if not lo_p or not hi_p:
        return None
    lo = float(_resolve_pointer(lo_p, what=f"{correction_id}.gate_lower_source"))
    hi = float(_resolve_pointer(hi_p, what=f"{correction_id}.gate_upper_source"))
    return abs(hi - lo) / 2.0


def value_check_required(correction_id: str) -> bool:
    """RYA-674 §2B's trigger, COMPUTED: is the correction large enough relative to the
    acceptance gate that the gate cannot catch a doubled application?

    True when ``|magnitude| >= 0.5 * gate_half_width`` — and also true when the entry
    declares no gate at all, because "no gate" is strictly less protection than a wide
    one. The Fe case is the archetype: |−0.05| against FE_GATE's 0.05 half-width, so a
    doubled correction lands exactly on the window edge and stays green (RYA-669).
    """
    half = gate_half_width(correction_id)
    if half is None:
        return True
    return abs(magnitude(correction_id)) >= 0.5 * half


def scale_bands(correction_id: str) -> dict[str, tuple[float, float]]:
    """The DERIVED pre/post "does this look already-corrected?" bands, per element.

    Returns ``{element: ...}`` is deliberately NOT the shape — a correction's scope may
    name several elements, so the caller asks per element via
    `pipeline.solar_scale_provenance.scale_centres`. This helper returns the bands for
    the entry's single declared post-scale centre, keyed by scale name:

        {'1D-NLTE': (lo, hi), '3D-NLTE': (lo, hi)}

    These are RYA-674 §2B's `pre_range` / `post_range`, computed from the magnitude
    rather than tabulated beside it.
    """
    entry = correction(correction_id)
    centre = float(_resolve_pointer(entry['post_scale_centre_source'],
                                    what=f"{correction_id}.post_scale_centre_source"))
    dex = magnitude(correction_id)
    half = abs(dex) / 2.0
    post_centre = centre
    pre_centre = round(centre - dex, 6)
    return {str(entry['post_scale']): (round(post_centre - half, 6), round(post_centre + half, 6)),
            str(entry['pre_scale']): (round(pre_centre - half, 6), round(pre_centre + half, 6))}


def assert_registry_consistent() -> None:
    """Every entry is well-formed and every pointer resolves. Loud on any breakage.

    Called by the test suite and by the emission-path guards, so a registry that has
    drifted from `config/constants.py` fails at the first use rather than silently
    resolving to a stale number.
    """
    reg = load_registry()
    for cid, entry in reg.items():
        if not isinstance(entry, dict):
            raise CorrectionsRegistryError(f"{cid}: entry is not a mapping")
        missing = [k for k in _REQUIRED_KEYS if k not in entry]
        if missing:
            raise CorrectionsRegistryError(f"{cid}: missing required key(s) {missing}")
        if not str(entry['provenance_ticket']).startswith('RYA-'):
            raise CorrectionsRegistryError(
                f"{cid}: provenance_ticket {entry['provenance_ticket']!r} is not an RYA-#; a "
                f"correction with no ratifying ticket may not be registered (RYA-674)")
        if not ((entry.get('scope') or {}).get('elements')):
            raise CorrectionsRegistryError(f"{cid}: scope.elements is empty")
        magnitude(cid)                      # resolves or raises
        if entry.get('post_scale_centre_source'):
            scale_bands(cid)                # resolves or raises
        gate_half_width(cid)                # resolves or raises
