"""RYA-1172 — "NIST grade" is not one authority across CNO. Vintage and method differ.

A C grade and an O grade are not the same pedigree, and nothing in the store said so. The
chain is stated by the compilers themselves in Wiese & Fuhr, NASA LAW 2006
(`Reference documents/20060052476.pdf`), and every string below is quoted from that paper
or its reference list -- none is written from memory.

🔴 THE 2006 UPDATE IS EXPLICITLY PARTIAL, AND OXYGEN IS NOT IN IT.

    title      "New Critical Compilations of Atomic Transition Probabilities for Neutral
                and Singly Ionized Carbon, Nitrogen, and Iron"        <- no oxygen
    scope      "We have thus carried out a partial update for the transition probabilities
                of C I, C II, N I and N II, especially utilizing the results of
                sophisticated MCHF calculations that were performed by Froese Fischer and
                coworkers in the last six years"
    the base   "Our 1996 data volume [Wiese (1996)] for C, N, and O was primarily based on
                the very extensive calculational results of the OPACITY Project"
    why        "we found ... that the OPACITY data are often not as accurate as we had
                estimated ... This statement applies especially to neutral and
                singly-ionized carbon and nitrogen."

So C I / C II / N I / N II carry MCHF; O I / O II still rest on the 1996 Opacity Project
values. Same "NIST grade" label, ten years and a different method apart.

⚠️ BOTH ARE CRITICALLY-EVALUATED THEORY, NOT LABORATORY MEASUREMENTS. The Opacity Project
is a calculation and MCHF is a calculation. This is the accepted standard for the light
elements and it is NOT a lab gf -- tiering these LAB would repeat RYA-1005's Al mistake.
`is_laboratory` is False on every entry and there is no code path that sets it True.

⚠️ AND THE HIGHER IONISATION STAGES ARE NOT COVERED BY EITHER SOURCE. The store holds
C III-C V, N III-N V and O III-O VI. The 2006 paper addresses four spectra; the 1996
Monograph's stage coverage is not stated in the document we hold. Those species get
`NOT_ESTABLISHED_BY_HELD_SOURCES` -- an unknown pedigree recorded as unknown, never
extended from a neighbouring species (RYA-1072: an allow-list for the recognised case must
not launder the unrecognised one).
"""
from __future__ import annotations

from dataclasses import dataclass, asdict

#: The document every string here is quoted from.
SOURCE_DOCUMENT = "Reference documents/20060052476.pdf"
SOURCE_CITATION = ("Wiese, W. L. & Fuhr, J. R., 'New Critical Compilations of Atomic "
                   "Transition Probabilities for Neutral and Singly Ionized Carbon, "
                   "Nitrogen, and Iron', NASA LAW, February 14-16, 2006, UNLV, Las Vegas")


@dataclass(frozen=True)
class Pedigree:
    """One authority chain. Mirrors RYA-1130's provenance-first shape for the ATOMIC axis:
    RYA-1130 types molecular provenance as (database, release, method-ish conventions);
    the atomic axis is (compilation, method, vintage), plus the evidence for the claim."""
    key: str
    compilation: str
    method: str
    vintage: int
    is_laboratory: bool
    evidence: str

    def as_dict(self) -> dict:
        return asdict(self)


MCHF_UPDATE = Pedigree(
    key="NIST_MCHF_PARTIAL_UPDATE_2006",
    compilation=SOURCE_CITATION,
    method="MCHF (multiconfiguration Hartree-Fock), Froese Fischer and coworkers "
           "[Froese Fischer (2004); Zatsarinny (2002)]; Breit-Pauli relativistic terms",
    vintage=2006,
    is_laboratory=False,
    evidence="'We have thus carried out a partial update for the transition probabilities "
             "of C I, C II, N I and N II' — 20060052476.pdf",
)

WFD1996_OPACITY = Pedigree(
    key="WFD1996_MONOGRAPH7_OPACITY_PROJECT",
    compilation="Wiese, W. L. Fuhr, J. R. and Deters, T. M. 1996, 'Atomic Transition "
                "Probabilities of Carbon, Nitrogen, and Oxygen, A Critical Data "
                "Compilation', J. Phys. Chem. Ref. Data, Monograph No. 7",
    method="OPACITY Project (theory) — The Opacity Project Team, The Opacity Project "
           "Vol. I (Institute of Physics, Bristol, England, 1995)",
    vintage=1996,
    is_laboratory=False,
    evidence="'Our 1996 data volume [Wiese (1996)] for C, N, and O was primarily based on "
             "the very extensive calculational results of the OPACITY Project'; the 2006 "
             "update covers 'Carbon, Nitrogen, and Iron' and names only C I, C II, N I, "
             "N II — oxygen is not in it — 20060052476.pdf",
)

UNESTABLISHED = Pedigree(
    key="NOT_ESTABLISHED_BY_HELD_SOURCES",
    compilation="",
    method="",
    vintage=0,
    is_laboratory=False,
    evidence="Neither 20060052476.pdf (four spectra: C I, C II, N I, N II) nor the "
             "Monograph No. 7 reference as quoted there states a pedigree for this "
             "ionisation stage. Recorded as unknown rather than extended from a "
             "neighbouring species.",
)

#: Species -> pedigree, ONLY where a held source states it.
PEDIGREE_BY_SPECIES: dict[str, Pedigree] = {
    "C I": MCHF_UPDATE, "C II": MCHF_UPDATE,
    "N I": MCHF_UPDATE, "N II": MCHF_UPDATE,
    "O I": WFD1996_OPACITY, "O II": WFD1996_OPACITY,
}


def pedigree_for(species: str) -> Pedigree:
    """The authority chain for a species. Unknown stages are UNESTABLISHED, never guessed."""
    return PEDIGREE_BY_SPECIES.get((species or "").strip(), UNESTABLISHED)


def is_cno(species: str) -> bool:
    head = (species or "").strip().split()
    return bool(head) and head[0] in ("C", "N", "O")
