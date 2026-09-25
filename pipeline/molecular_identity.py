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
        "OH", "X-X", "(2-0)/(3-1)/(4-2)", "O",
        "CRIRES+ H, 15277-16910 A: OH vibration-rotation, NOT an electronic system. "
        "Previously stamped 12C14N A-X (0-0). Band from the diagnostic registry's own "
        "reference string (Brooke+2016 gf), not inferred."),
    "CO_K": MolecularIdentity(
        "12C16O", "X-X", "(2-0)/(3-1)", "C",
        "CRIRES+ K, 23060-24572 A: CO first overtone, NOT an electronic system. "
        "Previously stamped 12C14N A-X (0-0). Band from the diagnostic registry "
        "(Li2015 gf)."),
    "CH_Gband": MolecularIdentity(
        "CH", "A-X", "G-band 4290-4315", "C",
        "HARPS VIS 4303.5-4313.0 A; Masseron+2014/2022. The primary solar and Procyon "
        "carbon indicator."),
    "C2_Swan": MolecularIdentity(
        "C2", "Swan", "(0-0)", "C",
        "HARPS VIS 5160-5166 A, Swan (0,0) bandhead at 5165."),
    "OH_AX": MolecularIdentity(
        "OH", "A-X", "(0-0)/(1-1)", "O",
        "near-UV 3063-3125 A: OH A-X (0,0) 3064 + (1,1) 3123. ELECTRONIC, unlike OH_H. "
        "The registry marks it an UPPER BOUND, and it is outside the Amarsi 2021 used "
        "OH range (1528-12280 nm) exactly as NH_AX is."),
}


#: Diagnostics that are NOT molecular bands. The old
#: `molecule = 'NH' if key == 'NH_AX' else '12C14N'` stamped these as CN A-X (0-0) too:
#: two C I atomic lines and a forbidden [O I] blend, all labelled as a CN molecular band.
#: An atomic diagnostic has a species and a transition, not a molecule and a vibrational
#: band, so it gets a pool of the right SHAPE rather than a molecular pool with the
#: molecule fields left wrong.
ATOMIC_IDENTITIES: dict[str, dict] = {
    "CI_5052": {"species": "C", "ion": "I", "kind": "atomic",
                "note": "C I 5052 A, HARPS VIS cross-check"},
    "CI_5380": {"species": "C", "ion": "I", "kind": "atomic",
                "note": "C I 5380 A, HARPS VIS cross-check"},
    "OI_6300": {"species": "O", "ion": "I", "kind": "forbidden_blend",
                "note": "[O I] 6300 A forbidden line, Ni I blended -- a FORBIDDEN "
                        "transition, not a permitted atomic line and not a band"},
}


class MolecularIdentityError(KeyError):
    """The diagnostic has no declared identity; it must not default to CN."""


def identity_for(diagnostic_key: str) -> MolecularIdentity:
    try:
        return IDENTITIES[diagnostic_key]
    except KeyError:
        raise MolecularIdentityError(
            f"{diagnostic_key!r} has no declared identity. If it is a molecular band, "
            f"declare it in IDENTITIES; if it is atomic or forbidden, in "
            f"ATOMIC_IDENTITIES. Do not let it fall back to 12C14N -- that is how OH_H, "
            f"CO_K and the [O I] 6300 forbidden blend came to be stamped CN A-X (0-0). "
            f"Declared molecular: {sorted(IDENTITIES)}; "
            f"atomic: {sorted(ATOMIC_IDENTITIES)}") from None


def pool_for(diagnostic_key: str, windows_air_A) -> dict:
    """The `pool` dict that gets hashed into indicator_id, shaped by the diagnostic.

    A molecular band is identified by molecule/system/band; an atomic or forbidden
    diagnostic by species/ion. Giving an atomic line a molecular pool is what put
    `12C14N A-X (0-0)` on the [O I] 6300 forbidden blend.
    """
    atomic = ATOMIC_IDENTITIES.get(diagnostic_key)
    if atomic is not None:
        return {"species": atomic["species"], "ion": atomic["ion"],
                "kind": atomic["kind"], "diagnostic": diagnostic_key,
                "windows_air_A": windows_air_A}
    ident = identity_for(diagnostic_key)
    return {"molecule": ident.molecule, "system": ident.system, "band": ident.band,
            "diagnostic": diagnostic_key, "windows_air_A": windows_air_A}
