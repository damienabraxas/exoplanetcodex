"""Method-selection policy resolver (RYA-315 / RYA-306).

For a given (star, species), which RT method do we trust: 'ew' or 'synthesis'?
Atomic species are EW-first; flip to synthesis only on evidence. Molecules are
synthesis-only a priori. A species we cannot classify is a SETUP ERROR
(fail loud) — never a silent default.

Species-key convention (matches the codebase): "<Element> <RomanIon>" — the
linelist (element, ion) pair joined by a space, e.g. "Fe I", "Fe II", "Ti I".
Molecular bands are keyed by formula, e.g. "CN". See data/method_policy.yaml.

This module stays star-agnostic (RYA-299): it does NOT own the canonical element
set — the caller passes `canonical_elements` (sourced from elements_master.json /
constants.TARGET_ELEMENTS). The policy file path is resolved via the same
config-path mechanism as everything else (config.constants.PATHS), not re-hardcoded.
"""
from __future__ import annotations

import functools

import yaml

from config.constants import PATHS

# Sourced from PATHS (RYA-292/288 pattern), not re-derived from __file__.
POLICY_PATH = PATHS['method_policy']


@functools.lru_cache(maxsize=1)
def _load_policy() -> dict:
    with open(POLICY_PATH) as f:
        return yaml.safe_load(f)


def _element_of(species: str) -> str:
    """'Fe II' -> 'Fe'; 'CN' -> 'CN'. Splits the codebase "<Element> <Ion>" key."""
    return species.split(" ")[0]


def get_method(star_id: str, species: str, *, canonical_elements: set[str]) -> dict:
    """Resolve the policy cell for (star_id, species).

    Returns {method, reason, evidence, confidence}.
    Fail-loud (ValueError) ONLY on a species we cannot classify — not a molecule
    and not a known atomic element. An un-overridden *known* atomic species
    correctly returns the EW default (that is the policy, not a gap).
    """
    policy = _load_policy()

    # 1. explicit override
    override = policy.get("overrides", {}).get(star_id, {}).get(species)
    if override is not None:
        return override

    # 2. molecular -> synthesis (a priori)
    if species in policy["molecular_species"] or _element_of(species) in policy["molecular_species"]:
        return {"method": policy["class_defaults"]["molecular"],
                "reason": "molecular", "evidence": "a-priori", "confidence": "a-priori"}

    # 3. known atomic -> EW default
    if _element_of(species) in canonical_elements:
        return {"method": policy["class_defaults"]["atomic"],
                "reason": "default", "evidence": "a-priori", "confidence": "default-pending-data"}

    # 4. unclassifiable -> SETUP ERROR (never silent-default)
    raise ValueError(
        f"method_policy: cannot classify species {species!r} for star {star_id!r} "
        f"- not a molecular band and not in the canonical element set. Setup error: "
        f"add it to elements_master / the policy; do not fall through to a silent default."
    )
