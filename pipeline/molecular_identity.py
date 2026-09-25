"""
RYA-1220 -- the molecular identity of a diagnostic, declared per diagnostic.

WHY THIS EXISTS. `nitrogen_cn_responses.py` built every route's pool as

    molecule = 'NH' if diag.key == 'NH_AX' else '12C14N'
    pool = {'molecule': molecule, 'system': 'A-X', 'band': '(0-0)', ...}

which special-cases one key and calls everything else CN A-X (0-0). OH_H (1.53-1.69 um)
and CO_K (2.31-2.46 um) were therefore both stamped 12C14N A-X (0-0) -- a band that has
no transitions anywhere near either window; CN A-X (0-0) lies at 1.09-1.32 um.

The label never reached the synthesis: the line list comes from `region_atomic_linelist`,
so the fits used the right transitions and no abundance is affected. But `pool` is hashed
into `indicator_id`, the identity key that `uncertainty_contract` and `publish_product`
join on, so two diagnostics were carrying a provenance claim about a molecule they do not
observe.

A ternary with a fallback is what made a partial fix look like a general one. Here every
diagnostic must be DECLARED; an unknown key raises rather than defaulting to CN.

Band designations are only asserted where the reference analysis pins them. Where we have
not established the vibrational assignment, the band is UNSPECIFIED -- an honest gap,
not an invented '(0-0)'.
"""
from __future__ import annotations

from typing import NamedTuple

UNSPECIFIED = "UNSPECIFIED"


class MolecularIdentity(NamedTuple):
    molecule: str
    system: str
    band: str
    element: str
    note: str


#: Declared identity per molecular diagnostic key. No fallback: adding a diagnostic
#: without adding it here is a hard error, which is the point.
IDENTITIES: dict[str, MolecularIdentity] = {
    "CN_AX_J": MolecularIdentity(
        "12C14N", "A-X", "(0-0)", "N",
        "CRIRES+ J, 11640-12076 A; inside the Amarsi 2021 used CN A-X range 10875-13208 A"),
    "CN_AX_IR": MolecularIdentity(
        "12C14N", "A-X", "(0-0)", "N",
        "IAG/KP NIR, 10871-11801 A; inside the Amarsi 2021 used CN A-X range"),
    "CN_red": MolecularIdentity(
        "12C14N", "A-X", UNSPECIFIED, "N",
        "6125-6200 A. The CN red system is A-X, but this is NOT the (0-0) band -- (0-0) "
        "lies at 1.09-1.32 um. The vibrational assignment here is not established, and "
        "the window is outside the Amarsi 2021 used set entirely."),
    "NH_AX": MolecularIdentity(
        "NH", "A-X", "(0-0)", "N",
        "3358-3373 A band head. Amarsi 2021 explicitly does not retain NH UV electronic "
        "lines; its NH is X-X vibration-rotation at 2.9-15 um."),
    "OH_H": MolecularIdentity(
        "OH", "X-X", UNSPECIFIED, "O",
        "CRIRES+ H, 15277-16910 A: OH vibration-rotation, NOT an electronic system. "
        "Previously stamped 12C14N A-X (0-0)."),
    "CO_K": MolecularIdentity(
        "12C16O", "X-X", UNSPECIFIED, "C",
        "CRIRES+ K, 23060-24572 A: CO vibration-rotation, NOT an electronic system. "
        "Previously stamped 12C14N A-X (0-0)."),
}


class MolecularIdentityError(KeyError):
    """The diagnostic has no declared identity; it must not default to CN."""


def identity_for(diagnostic_key: str) -> MolecularIdentity:
    try:
        return IDENTITIES[diagnostic_key]
    except KeyError:
        raise MolecularIdentityError(
            f"{diagnostic_key!r} has no declared molecular identity. Declare it in "
            f"pipeline.molecular_identity.IDENTITIES -- do not let it fall back to "
            f"12C14N, which is how OH_H and CO_K came to be stamped as CN A-X (0-0). "
            f"Declared: {sorted(IDENTITIES)}") from None


def pool_for(diagnostic_key: str, windows_air_A) -> dict:
    """The `pool` dict that gets hashed into indicator_id."""
    ident = identity_for(diagnostic_key)
    return {"molecule": ident.molecule, "system": ident.system, "band": ident.band,
            "diagnostic": diagnostic_key, "windows_air_A": windows_air_A}
