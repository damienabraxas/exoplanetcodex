#!/usr/bin/env python3
"""
pipeline/model_atom.py — RYA-818
================================
Reading Gerber TS-native `atom.*` model atoms, and — the part that matters —
telling whether an atom's levels are FINE-STRUCTURE or TERM-RESOLVED.

WHY THIS EXISTS
---------------
`scripts/rya763_level_mapping.py` already read these files, and RYA-776 imports
that reader. Both then identify a level by `(J, energy)`, where J is taken as
`(g - 1) / 2`. That is correct only for a fine-structure level, where g = 2J+1.

`atom.cr374` is not such an atom. Its levels are LS TERMS — every fine-structure
component of a term collapsed into one super-level — so its g is the term's TOTAL
statistical weight, not any single level's. `a5D` carries g=25, and the resolver
duly computes J=12 for it. No real Cr line has J=12, so every Cr line resolves as
ABSENT and Cr's Engine-B reach reads as zero.

That is a MANUFACTURED ABSENCE, and it is the same failure class as RYA-776's
super-level trap and RYA-763's "level index is a false coordinate". Measured: Cr
crossmatches at 2.1% by wavelength/J and 88.3% by term label. The 2% is an
artefact of asking the wrong question.

THE TEST THAT DECIDES IT
------------------------
Summing (2J+1) over the J of an LS term telescopes:

    sum_J (2J+1)  =  (2S+1)(2L+1)          for J = |L-S| .. L+S

So a term-resolved level satisfies `g == (2S+1)(2L+1)` EXACTLY, computable from
the label alone. Checked against `atom.cr374`: `a5D` 5x5=25, `a3H` 3x11=33,
`z7F*` 7x7=49, `a6D` 6x5=30 — 264 of its 374 levels confirm the identity. This is
an exact arithmetic check, not a heuristic threshold, which is what makes the
claim falsifiable rather than a judgement call.

WHAT THIS MODULE DOES NOT DO
----------------------------
It does not decide how to MATCH lines — that is
`pipeline/nlte_line_identification.py`. It reports what the atom is, so callers
stop assuming.
"""
from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import pandas as pd

CM1_PER_EV = 8065.543937

#: Orbital angular momentum letters, in spectroscopic order. J and P are skipped by
#: convention (P is used, J is not a term letter); this is the standard sequence.
L_OF_LETTER = {"S": 0, "P": 1, "D": 2, "F": 3, "G": 4, "H": 5, "I": 6,
               "K": 7, "L": 8, "M": 9, "N": 10, "O": 11, "Q": 12, "R": 13, "T": 14}

#: A bare LS term: optional series letter, multiplicity, L letter, optional parity.
#: e.g. a5D  z7F*  b4P  3F*
_BARE_TERM = re.compile(r"^([a-z]?)(\d+)([SPDFGHIKLMNOQRT])(\*?)$")

#: A J-resolved level: the same, plus J (integer or half-integer written as `.5`).
#: e.g. z3D3*  a2D1  a8S3*
_J_LEVEL = re.compile(r"^([a-z]?)(\d+)([SPDFGHIKLMNOQRT])(\d+(?:\.5)?)(\*?)$")

#: A level line in a Gerber atom: energy(cm-1), g, 'Level N = LABEL', ion stage.
_LEVEL_LINE = re.compile(r"\s*([-\d.Ee+]+)\s+([-\d.Ee+]+)\s+'([^']*)'\s+(\d+)")

TERM_RESOLVED = "term-resolved"
FINE_STRUCTURE = "fine-structure"
MIXED = "mixed"
UNDECIDED = "undecided"


class ModelAtomError(RuntimeError):
    """The atom file is not readable as a Gerber TS-native model atom."""


def term_statistical_weight(label: str) -> Optional[int]:
    """(2S+1)(2L+1) for a bare LS term label, else None.

    This is the total weight of the WHOLE term — what a super-level carries.
    Returns None for labels that are not bare terms (J-resolved or config-style),
    because for those the question does not apply.
    """
    m = _BARE_TERM.match(label.strip())
    if not m:
        return None
    mult, letter = int(m.group(2)), m.group(3)
    if letter not in L_OF_LETTER:
        return None
    return mult * (2 * L_OF_LETTER[letter] + 1)


def level_j(label: str) -> Optional[float]:
    """J parsed from a J-resolved label (`z3D3*` -> 3.0), else None.

    Deliberately parses J from the LABEL rather than deriving it from g. Deriving
    J as (g-1)/2 is exactly the step that fabricates J=12 for a super-level.
    """
    m = _J_LEVEL.match(label.strip())
    return float(m.group(4)) if m else None


def classify_label(label: str, g: float) -> str:
    """One level -> 'term' | 'j-resolved' | 'unparsed', confirmed against g.

    A label only counts as `term` when the (2S+1)(2L+1) identity HOLDS. A bare term
    whose g disagrees is left `unparsed`: it is not evidence either way, and
    counting it would let a threshold do work the arithmetic should do.
    """
    j = level_j(label)
    if j is not None:
        return "j-resolved" if abs(g - (2 * j + 1)) < 1e-6 else "unparsed"
    w = term_statistical_weight(label)
    if w is not None and abs(g - w) < 1e-6:
        return "term"
    return "unparsed"


@dataclass(frozen=True)
class AtomResolution:
    """What kind of level coordinates an atom actually uses."""
    verdict: str
    n_levels: int
    n_term: int
    n_j_resolved: int
    n_unparsed: int
    element: str

    @property
    def is_term_resolved(self) -> bool:
        return self.verdict == TERM_RESOLVED

    def describe(self) -> str:
        return (
            f"{self.element}: {self.n_levels} levels — {self.n_term} confirmed LS "
            f"terms (g == (2S+1)(2L+1)), {self.n_j_resolved} confirmed J-resolved "
            f"(g == 2J+1), {self.n_unparsed} unparsed -> {self.verdict.upper()}")


def read_gerber_atom(path: str | Path) -> pd.DataFrame:
    """Levels from a Gerber TS-native `atom.*` file.

    Columns: index (1-based, the atom's own level id), element, ion, term, g,
    energy_cm, energy_eV, kind.

    `element` is read from the file's FIRST LINE rather than passed in — the atom
    states its own element, and a caller-supplied one can disagree with the file.
    (`rya763_level_mapping.read_gerber_atom` hardcoded `species=f"Fe {ion}"`,
    which mislabels every non-Fe atom it is handed.)

    NOTE `g` is preserved as read. It is NOT converted to J here, because that
    conversion is only valid for fine-structure levels — see the module docstring.
    """
    p = Path(path)
    text = p.read_text(errors="replace").splitlines()
    if not text:
        raise ModelAtomError(f"{p}: empty file")

    element = text[0].strip().split()[0] if text[0].strip() else ""
    if not element:
        raise ModelAtomError(f"{p}: first line carries no element name")

    rows: list[dict] = []
    seen_header = False
    n_expect: Optional[int] = None
    for raw in text:
        line = raw.rstrip()
        if not line.strip() or line.lstrip().startswith(("*", "#")):
            continue
        parts = line.split()
        if not seen_header:
            # the counts line: n_levels n_transitions n_continua n_...
            if n_expect is None and len(parts) >= 3 and all(
                    x.lstrip("-").isdigit() for x in parts[:3]):
                n_expect = int(parts[0])
                seen_header = True
            continue
        m = _LEVEL_LINE.match(line)
        if not m:
            if len(rows) >= (n_expect or 0):
                break
            continue
        energy_cm, g, label, ion = m.groups()
        try:
            e, gg = float(energy_cm), float(g)
        except ValueError:
            continue
        term = label.split("=")[-1].strip() if "=" in label else label.strip()
        rows.append(dict(index=len(rows) + 1, element=element, ion=int(ion),
                         term=term, g=gg, energy_cm=e, energy_eV=e / CM1_PER_EV,
                         kind=classify_label(term, gg)))
        if n_expect and len(rows) >= n_expect:
            break

    if not rows:
        raise ModelAtomError(f"{p}: no level lines parsed")
    if n_expect and len(rows) != n_expect:
        raise ModelAtomError(
            f"{p}: header declares {n_expect} levels, parsed {len(rows)}. A short "
            f"read here silently shrinks the atom, so this refuses rather than "
            f"returning a partial level table.")
    return pd.DataFrame(rows)


def atom_resolution(levels: pd.DataFrame, *, term_threshold: float = 0.5
                    ) -> AtomResolution:
    """Is this atom term-resolved, fine-structure, or mixed?

    `term_threshold` gates only the SUMMARY verdict; the per-level `kind` column is
    exact arithmetic and carries no threshold. Callers that need to be careful
    should branch on the counts, not on the word.
    """
    n = len(levels)
    if n == 0:
        return AtomResolution(UNDECIDED, 0, 0, 0, 0, "")
    kinds = levels["kind"].value_counts()
    n_term = int(kinds.get("term", 0))
    n_j = int(kinds.get("j-resolved", 0))
    n_un = int(kinds.get("unparsed", 0))
    element = str(levels["element"].iloc[0])

    if n_term / n > term_threshold:
        verdict = TERM_RESOLVED
    elif n_j / n > term_threshold:
        verdict = FINE_STRUCTURE
    elif n_term and n_j:
        verdict = MIXED
    elif n_term:
        verdict = MIXED
    else:
        # Neither coordinate confirmed — usually a config-style label set
        # (`3s2.3p2P*`). Saying UNDECIDED is the honest answer; calling it
        # fine-structure by default is what lets a super-level through unnoticed.
        verdict = UNDECIDED
    return AtomResolution(verdict, n, n_term, n_j, n_un, element)


def ion_stage_histogram(levels: pd.DataFrame) -> dict[int, int]:
    """{ion_stage: n_levels}. Guard this before any per-stage join."""
    return {int(k): int(v) for k, v in levels["ion"].value_counts().sort_index().items()}
