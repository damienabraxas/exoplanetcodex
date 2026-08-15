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
    if n == 0 or "kind" not in levels.columns:
        # No `kind` means an Amarsi/PySME `label_*.txt` (Engine A), which carries no
        # `g` and so cannot be classified by the (2S+1)(2L+1) identity at all.
        # UNDECIDED is the honest answer and, crucially, is NOT term-resolved — so
        # callers branching on `is_term_resolved` leave Engine A exactly as it was.
        return AtomResolution(UNDECIDED, n, 0, 0, 0, "")
    kinds = levels["kind"].value_counts()
    n_term = int(kinds.get("term", 0))
    n_j = int(kinds.get("j-resolved", 0))
    n_un = int(kinds.get("unparsed", 0))
    # `element` is absent from the RYA-763-shaped frame (it carries `species`), and
    # this is called from `resolvable_j` on exactly that frame. The verdict does not
    # depend on the name, so a missing one is cosmetic, not a reason to fail.
    if "element" in levels.columns:
        element = str(levels["element"].iloc[0])
    elif "species" in levels.columns:
        element = str(levels["species"].iloc[0]).split()[0]
    else:
        element = ""

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


def resolvable_j(levels: pd.DataFrame) -> "pd.Series":
    """J per level, NaN where J does not exist — RYA-823.

    A fine-structure level has J = (g-1)/2. A SUPER-LEVEL does not have a J at all:
    its g is the whole term's (2S+1)(2L+1), so (g-1)/2 returns a number
    (`a5D` -> 12.0) that no line can ever carry. Handing that number to a
    (J, energy) matcher does not merely fail to match — it occasionally matches
    something by accident, and those accidents are worse than a clean miss:

      * Cr I VIS  4 of 5353 lines "resolve" (0.1%)
      * Cr II VIS 8 of 3660 (0.2%)

    which is enough to defeat a `nothing resolved, so the key must be wrong` guard
    and let a 0.1% figure be published as a measured reach. NaN here makes the
    absence clean, so that guard can see it.

    ⚠️ SUPPRESSION REQUIRES POSITIVE EVIDENCE, NOT MERE DOUBT. `kind == "unparsed"`
    means the label is in a form this module cannot read (`3s2.3p2P*` config style,
    `10p7P*` Rydberg series) — it does NOT mean super-level. A first cut nulled J for
    everything that was not confirmed `j-resolved`, and the regenerated coverage
    table showed what that costs: **Ca I VIS 55 -> 0 reach, Co I VIS 993 -> 0 and
    reclassified UNCOVERED**, 106 rows changed. Most atoms label most levels in a
    form the parser does not read, and for those `(g-1)/2` is both the best available
    key and demonstrably a working one.

    So J is withheld only in an atom that classifies as TERM-RESOLVED OVERALL — i.e.
    where the fabricated J is the only thing on offer and cannot be right. Cr is that
    atom (274 bare terms, 100 unreadable Rydberg labels, zero confirmed J-resolved
    levels), and withholding lets it fall cleanly to nothing-resolves.

    ⚠️ A MIXED atom keeps every J, including its super-levels', and that is
    deliberate restraint rather than an oversight. Withholding there costs real
    reach — Mn II VIS 171 -> 93, Mn I VIS 450 -> 421, Ti I VIS 1534 -> 1506 — and I
    cannot show those lost matches were WRONG. A low-g super-level's (g-1)/2 can
    coincide with a genuine J (`a6S` g=6 gives 2.5), the energy still had to agree
    within 1 meV to match at all, and where it does agree the coefficient returned is
    the term's, which is exactly what the label route would have returned anyway. So
    the match may be arrived at by luck but still be serviceable.

    Deciding that needs the raw-GES label route (RYA-818), because the iSpec frame's
    `nlte_label_*` columns are populated only for NLTE-tagged species and disagree
    with these atoms' own labels besides (Mn II label reach 0.5%). Until that is
    measured, this change stays MONOTONE: the label union can only ADD matches, and
    no element's reported reach is cut on a judgement call.
    """
    import numpy as np
    j = (levels["g"].astype(float) - 1.0) / 2.0
    if atom_resolution(levels).verdict == TERM_RESOLVED:
        return j.where(levels["kind"] == "j-resolved", other=np.nan)
    return j


def resolve_level(levels: pd.DataFrame, *, energy_eV: float = float("nan"),
                  j: float = float("nan"), term: str = "",
                  tol_eV: float = 0.001) -> tuple[str, int, int]:
    """Identify ONE level by whichever key its kind actually supports — RYA-823.

    Returns (verdict, level_index, n_candidates) with verdict in
    UNIQUE / AMBIGUOUS / ABSENT / NO-KEY, matching `rya763_level_mapping`'s
    vocabulary so callers need not learn a second one.

    THE UNION IS THE POINT, AND IT IS NOT A SWAP
    --------------------------------------------
    (J, energy) addresses fine-structure levels; the term label addresses
    super-levels. Neither dominates, measured per species:

        Ti I red-optical   (J,energy) 34.8%   label 56.2%   -> label wins
        Mn II VIS          (J,energy) 10.1%   label  0.5%   -> J wins, 20x

    So replacing one key with the other would have cost Mn II 95% of its reach
    while claiming to fix it. Each level is matched by the key appropriate to its
    OWN kind, and the results are unioned.

    Ambiguity refuses. Two keys agreeing on one level is one answer; two keys
    naming DIFFERENT levels is a contradiction, and picking either would attach a
    real departure coefficient to a guess.
    """
    import numpy as np
    hits: set[int] = set()
    used_a_key = False

    # A level table that carries no `kind` is an Amarsi/PySME `label_*.txt` (Engine A),
    # which is fine-structure throughout and has no `g`. Its own J is authoritative
    # there, so this falls back to it rather than demanding columns that deck never
    # had -- the union must not change Engine A's answers at all.
    if "kind" in levels.columns and "g" in levels.columns:
        jj = resolvable_j(levels)
        term_mask = levels["kind"] == "term"
    else:
        jj = levels["J"].astype(float)
        term_mask = None

    # NO-KEY must mean "this TABLE cannot be addressed", not "the caller passed
    # nothing". A term-resolved atom has NaN for every J, so a (J, energy) query
    # against it returns ABSENT -- which reads as a measured absence -- when the
    # truth is that the key does not apply. Requiring the table to offer at least one
    # real J is what turns that into an honest NO-KEY.
    j_offered = bool(np.isfinite(jj).any())
    if j_offered and np.isfinite(energy_eV) and np.isfinite(j):
        used_a_key = True
        c = levels[(np.abs(levels["energy_eV"].astype(float) - energy_eV) <= tol_eV)
                   & (np.abs(jj - j) < 0.01)]
        hits.update(int(x) for x in c["index"])

    t = (term or "").strip()
    if t and t.lower() != "none" and term_mask is not None:
        used_a_key = True
        c = levels[(levels["term"].astype(str).str.strip() == t) & term_mask]
        hits.update(int(x) for x in c["index"])

    if not used_a_key:
        return ("NO-KEY", -1, 0)
    if len(hits) == 1:
        return ("UNIQUE", hits.pop(), 1)
    if len(hits) > 1:
        return ("AMBIGUOUS", -1, len(hits))
    return ("ABSENT", -1, 0)
