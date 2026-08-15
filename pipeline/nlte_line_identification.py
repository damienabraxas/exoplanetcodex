#!/usr/bin/env python3
"""
pipeline/nlte_line_identification.py — RYA-818
==============================================
Give a Turbospectrum linelist the per-line NLTE level identifications a model atom
needs, for atoms whose levels are LS TERMS rather than fine-structure levels.

THE ARTIFACT BEING PRODUCED
---------------------------
An NLTE-tagged row in `nlte_ges_linelist_jmg17feb2022_I_II` carries six trailing
fields an LTE row does not:

    ... 'Fe I LS:... z3D* LS:... f3F'   82  0   'z3D3*'  'none'   'c' 'x'
                                        |   |   |        |        |    |
                             lower level id   |   lower label     |    upper flag
                                 upper level id     upper label   lower flag

The flags are MATCH-QUALITY CODES, kept compatible with the upstream TSFitPy
converter's vocabulary:

    'c'  identified by LEVEL LABEL
    'a'  identified by ENERGY, after the label route failed
    'b'  identified from a bound-bound transition (upstream's wavelength route;
         this module does not use it — see WHY WAVELENGTH IS THE WRONG KEY)
    'x'  NOT IDENTIFIED

A line flagged `x` on either endpoint runs in **LTE** even though it sits inside an
NLTE block. That is not a cosmetic detail: RYA-764 found 2,644 in-band Fe lines in
exactly that state, silently. So `x` is counted and reported here, never rounded
into "coverage".

WHY WAVELENGTH IS THE WRONG KEY HERE
------------------------------------
The upstream converter's primary route matches a line to a bound-bound transition
within +/-0.02 A. For a TERM-RESOLVED atom that route is structurally wrong: the
atom's transition wavelength is a term-to-term average, while the linelist's is one
fine-structure component. Measured on Cr, the median disagreement is 0.632 A (Cr I)
and 0.783 A (Cr II) — an order of magnitude past the tolerance. Wavelength matching
returns 2.1% / 1.8%; the term-label join returns 88.3% / 38.5%. The 2% is a
manufactured absence, and must never be quoted as reach.

THE PHYSICS BEING ASSUMED
-------------------------
Every fine-structure component of a term is assigned THAT TERM'S departure
coefficient. For a term-resolved atom this is the correct reading — the atom has no
finer coordinate to offer — but it IS an approximation, and it is invisible in the
output. `identification_provenance()` states it so the product carries it.

🔴 DO NOT MATCH ON MULTIPLICITY
-------------------------------
The upstream converter has a `match_multiplicity` option requiring the line's
2J+1 to equal the level's g. Against a term-resolved atom that is guaranteed to
fail: the line's g_up is ONE component's 2J+1, the super-level's g is the whole
term's sum. Measured agreement on Cr is 2-5%. Enabling it would collapse reach
from ~88% to ~3% while looking like a tightened cross-check.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Optional

import pandas as pd

from pipeline.model_atom import CM1_PER_EV

#: Match-quality flags, upstream-compatible.
FLAG_LABEL = "c"
FLAG_ENERGY = "a"
FLAG_TRANSITION = "b"
FLAG_UNMATCHED = "x"

#: First ionisation potentials (eV), for putting a GES excitation potential onto the
#: atom's GLOBAL energy scale. A Gerber atom numbers every stage from the NEUTRAL
#: ground state (atom.cr374's first Cr II level sits at 54580.617 cm-1 = 6.767 eV,
#: which IS the Cr I IP), while GES quotes excitation within the line's own stage.
#: Values: NIST ASD ionisation energies.
FIRST_IP_EV = {
    "Cr": 6.76651,
    "Fe": 7.9024681,
    "Y": 6.21726,
    "Eu": 5.670385,
    "Al": 5.985769,
    "Ti": 6.828120,
    "Mn": 7.4340380,
    "Co": 7.881010,
}

#: A quoted species header inside the linelist, e.g. `'Cr I    LTE'`.
_SPECIES_HEADER = re.compile(r"^'\s*([A-Za-z]+)\s+([IVX]+)")

_ROMAN = {"I": 1, "II": 2, "III": 3}


class LineIdentificationError(RuntimeError):
    """The inputs cannot support an honest identification."""


# ── air -> vacuum, and the upper-level energy ────────────────────────────────
def air_to_vacuum(wave_air_A: float) -> float:
    """Birch & Downs / Edlen dispersion, matching the upstream converter exactly.

    Kept bit-compatible with TSFitPy's constants deliberately: the upper-level
    energy derived here is compared against the same atoms upstream uses, so a
    different dispersion formula would shift every computed EC by a small amount
    and quietly change which levels fall inside an energy tolerance.
    """
    sigma2 = (1.0e4 / wave_air_A) ** 2
    fact = (1.0 + 8.336624212083e-5
            + 2.408926869968e-2 / (1.301065924522e2 - sigma2)
            + 1.599740894897e-4 / (3.892568793293e1 - sigma2))
    return wave_air_A * fact


def level_energies_cm(wave_air_A: float, ep_eV: float, ionisation_ev: float
                      ) -> tuple[float, float]:
    """(lower, upper) level energy in cm-1 on the atom's GLOBAL scale."""
    lower = CM1_PER_EV * (ep_eV + ionisation_ev)
    upper = lower + 1.0e8 / air_to_vacuum(wave_air_A)
    return lower, upper


def ionisation_offset_eV(element: str, ion: int) -> float:
    """Energy from the neutral ground state to this ion stage's ground state.

    Only stages I and II are supported: the offset for stage III needs the SECOND
    ionisation potential, and inventing one would silently misplace every level by
    several eV. Refuses instead.
    """
    if ion == 1:
        return 0.0
    if ion == 2:
        try:
            return FIRST_IP_EV[element]
        except KeyError:
            raise LineIdentificationError(
                f"no first ionisation potential recorded for {element}; add it to "
                f"FIRST_IP_EV from NIST ASD rather than defaulting to zero, which "
                f"would put every {element} II level on the wrong energy scale.")
    raise LineIdentificationError(
        f"ion stage {ion} is not supported: the offset needs the {ion - 1}th "
        f"ionisation potential, which is not recorded here.")


# ── reading the linelist ─────────────────────────────────────────────────────
def _terms_from_label(label: str) -> tuple[str, str]:
    """('a3G', '3F*') from `Cr I LS:3d5.(4G).4s a3G LS:3d4.(...) 3F*`.

    The term is the last whitespace-token of each `LS:` segment. Returns ('','')
    when the label does not carry two segments — an unlabelled line, which is a
    real and common case, not an error.
    """
    parts = label.split("LS:")
    if len(parts) < 3:
        return "", ""
    lo = parts[1].strip().split()[-1] if parts[1].strip() else ""
    up = parts[2].strip().split()[-1] if parts[2].strip() else ""
    return lo, up


def read_species_lines(path: str | Path, element: str, ion: int) -> pd.DataFrame:
    """Every line of one species from a Turbospectrum linelist.

    Columns: row (position within the species block), wave_A, ep_eV, loggf, g_up,
    term_low, term_up, raw.
    """
    want = (element.upper(), ion)
    rows: list[dict] = []
    current: Optional[tuple[str, int]] = None
    for raw in Path(path).read_text(errors="replace").splitlines():
        if raw.startswith("'"):
            m = _SPECIES_HEADER.match(raw)
            if m:
                stage = _ROMAN.get(m.group(2).upper())
                current = (m.group(1).upper(), stage) if stage else None
            continue
        if current != want or not raw.startswith(" "):
            continue
        fields = raw.split()
        if len(fields) < 5:
            continue
        try:
            wave, ep, loggf, g_up = (float(fields[0]), float(fields[1]),
                                     float(fields[2]), float(fields[4]))
        except ValueError:
            continue
        quoted = re.findall(r"'([^']*)'", raw)
        lo, up = _terms_from_label(quoted[-1] if quoted else "")
        rows.append(dict(row=len(rows), wave_A=wave, ep_eV=ep, loggf=loggf,
                         g_up=g_up, term_low=lo, term_up=up, raw=raw))
    return pd.DataFrame(rows)


# ── the identification itself ────────────────────────────────────────────────
@dataclass
class LevelResolver:
    """Resolves a level of ONE ion stage, by label first and energy second.

    Constructed per (element, ion) so a Cr I line can never be handed a Cr II
    level. That is not defensive tidiness: `ModelAtom`-style level tables expose
    the stage as `ionisation_stage`, and reading a differently-named attribute
    yields None for every level and pools the stages into one namespace. Making
    the stage a CONSTRUCTOR argument removes the chance to forget it.
    """
    levels: pd.DataFrame
    ion: int
    energy_tol_cm: float = 50.0

    _by_term: dict[str, list] = field(default_factory=dict, init=False)
    _energies: list = field(default_factory=list, init=False)

    def __post_init__(self) -> None:
        sub = self.levels[self.levels["ion"] == self.ion]
        if sub.empty:
            raise LineIdentificationError(
                f"the atom carries no levels for ion stage {self.ion}; refusing to "
                f"identify lines against an empty stage.")
        for rec in sub.to_dict("records"):
            self._by_term.setdefault(str(rec["term"]).strip(), []).append(rec)
        self._energies = sorted(sub.to_dict("records"),
                                key=lambda r: float(r["energy_cm"]))

    def by_label(self, term: str) -> tuple[Optional[dict], str]:
        """(level, reason). Ambiguity REFUSES rather than picking the first."""
        if not term:
            return None, "no-label"
        hits = self._by_term.get(term.strip())
        if not hits:
            return None, "label-absent"
        if len(hits) > 1:
            return None, f"label-ambiguous({len(hits)})"
        return hits[0], "label"

    def by_energy(self, energy_cm: float) -> tuple[Optional[dict], str]:
        """Nearest level within tolerance; ambiguity REFUSES."""
        cands = [r for r in self._energies
                 if abs(float(r["energy_cm"]) - energy_cm) <= self.energy_tol_cm]
        if not cands:
            return None, "energy-absent"
        if len(cands) > 1:
            return None, f"energy-ambiguous({len(cands)})"
        return cands[0], "energy"

    def resolve(self, term: str, energy_cm: float, *, allow_energy: bool = False
                ) -> tuple[Optional[dict], str, str]:
        """(level, flag, reason). Label first; energy only if explicitly allowed.

        🔴 `allow_energy` DEFAULTS TO FALSE, and that default is a finding, not
        caution. Against a term-resolved atom the energy route is measurably
        WRONG: run where the label route also answers (so the truth is known), it
        picks a different level 18% of the time for Cr I and 35% for Cr II. The
        cause is structural — the atom's energies are term AVERAGES, and
        neighbouring terms sit closer together than the fine-structure spread they
        average over, so a component's computed energy routinely lands nearest the
        wrong term. Cr I's a3P (23796 cm-1) and a3H (24080 cm-1) are 284 cm-1
        apart while the components spread further than that.

        Tightening does not fix it: sweeping tolerance 5-50 cm-1 against a
        runner-up separation margin of 1-10x, the best agreement reached was 86.7%
        (Cr I) / 86.2% (Cr II), and only by discarding ~two-thirds of matches. A
        route that is confidently wrong one time in seven must not run by default.
        """
        lvl, why = self.by_label(term)
        if lvl is not None:
            return lvl, FLAG_LABEL, why
        if not allow_energy:
            return None, FLAG_UNMATCHED, f"{why}->energy-route-disabled"
        lvl, why2 = self.by_energy(energy_cm)
        if lvl is not None:
            return lvl, FLAG_ENERGY, f"{why}->{why2}"
        return None, FLAG_UNMATCHED, f"{why}->{why2}"


def identify_lines(lines: pd.DataFrame, levels: pd.DataFrame, element: str,
                   ion: int, *, energy_tol_cm: float = 50.0,
                   energy_fallback: bool = False) -> pd.DataFrame:
    """Identify both endpoints of every line. Adds the six emitted fields.

    Columns added: level_low, level_up (atom level ids, 0 = unidentified),
    label_low, label_up ('none' = unidentified), flag_low, flag_up, reason_low,
    reason_up, and `nlte` (True only when BOTH endpoints are identified).

    `energy_fallback` defaults to FALSE — see `LevelResolver.resolve`. Turning it
    on buys reach at a measured ~1-in-6 chance of naming the wrong level, so any
    caller that enables it owns that trade and must report it.
    """
    resolver = LevelResolver(levels=levels, ion=ion, energy_tol_cm=energy_tol_cm)
    offset = ionisation_offset_eV(element, ion)

    out = []
    for rec in lines.to_dict("records"):
        e_lo, e_up = level_energies_cm(rec["wave_A"], rec["ep_eV"], offset)
        lo, f_lo, r_lo = resolver.resolve(rec["term_low"], e_lo,
                                          allow_energy=energy_fallback)
        up, f_up, r_up = resolver.resolve(rec["term_up"], e_up,
                                          allow_energy=energy_fallback)
        out.append(dict(
            rec,
            energy_low_cm=e_lo, energy_up_cm=e_up,
            level_low=int(lo["index"]) if lo else 0,
            level_up=int(up["index"]) if up else 0,
            label_low=str(lo["term"]) if lo else "none",
            label_up=str(up["term"]) if up else "none",
            flag_low=f_lo, flag_up=f_up, reason_low=r_lo, reason_up=r_up,
            nlte=bool(lo is not None and up is not None)))
    return pd.DataFrame(out)


def render_identification_fields(row) -> str:
    """The six trailing fields, in the linelist's own spacing."""
    return ("  %d %d  '%s' '%s'  '%s' '%s'"
            % (int(row["level_low"]), int(row["level_up"]),
               row["label_low"], row["label_up"],
               row["flag_low"], row["flag_up"]))


def reach_report(identified: pd.DataFrame) -> dict:
    """Reach, broken out PER FLAG. An `x` line runs in LTE and is counted as such."""
    n = len(identified)
    if n == 0:
        return {"n_lines": 0}
    both = int(identified["nlte"].sum())
    return {
        "n_lines": n,
        "n_nlte": both,
        "reach_pct": round(100.0 * both / n, 2),
        "n_lte_despite_nlte_block": n - both,
        "lower_by_flag": identified["flag_low"].value_counts().to_dict(),
        "upper_by_flag": identified["flag_up"].value_counts().to_dict(),
        "n_label_both": int(((identified["flag_low"] == FLAG_LABEL)
                             & (identified["flag_up"] == FLAG_LABEL)).sum()),
        "n_energy_assisted": int((identified["nlte"]
                                  & ((identified["flag_low"] == FLAG_ENERGY)
                                     | (identified["flag_up"] == FLAG_ENERGY))).sum()),
    }


def energy_route_agreement(lines: pd.DataFrame, levels: pd.DataFrame, element: str,
                           ion: int, *, energy_tol_cm: float = 50.0) -> dict:
    """POSITIVE CONTROL for the energy fallback.

    The fallback only ever fires where the label route failed, so on its own it is
    unfalsifiable: nothing it produces can be checked against a known answer. This
    runs it where the label route SUCCEEDED — cases it never sees in production —
    and asks whether it returns the same level.

    A high agreement rate means the energy route recovers the right level when the
    label is missing. A low one means it is confidently assigning wrong levels, and
    every `a` flag in the product is suspect. Reported as `agree / (agree+disagree)`
    with the disagreements kept, so the control can FAIL rather than merely reassure.
    """
    resolver = LevelResolver(levels=levels, ion=ion, energy_tol_cm=energy_tol_cm)
    offset = ionisation_offset_eV(element, ion)
    agree = disagree = unreachable = 0
    examples: list[dict] = []
    for rec in lines.to_dict("records"):
        e_lo, e_up = level_energies_cm(rec["wave_A"], rec["ep_eV"], offset)
        for term, energy, end in ((rec["term_low"], e_lo, "low"),
                                  (rec["term_up"], e_up, "up")):
            truth, why = resolver.by_label(term)
            if truth is None:
                continue                      # label route did not answer: no control
            guess, _ = resolver.by_energy(energy)
            if guess is None:
                unreachable += 1
            elif int(guess["index"]) == int(truth["index"]):
                agree += 1
            else:
                disagree += 1
                if len(examples) < 8:
                    examples.append(dict(wave_A=rec["wave_A"], end=end, term=term,
                                         label_level=int(truth["index"]),
                                         energy_level=int(guess["index"]),
                                         label_E=float(truth["energy_cm"]),
                                         energy_E=float(guess["energy_cm"]),
                                         computed_E=energy))
    decided = agree + disagree
    return {
        "controlled_endpoints": decided + unreachable,
        "agree": agree, "disagree": disagree,
        "energy_route_silent": unreachable,
        "agreement_pct": round(100.0 * agree / decided, 2) if decided else None,
        "disagreement_examples": examples,
    }


def identification_provenance(element: str, atom_name: str, resolution,
                              energy_tol_cm: float, *,
                              energy_fallback: bool = False) -> dict:
    """What a downstream reader must know to judge these identifications."""
    prov = {
        "element": element,
        "model_atom": atom_name,
        "atom_resolution": resolution.verdict,
        "matching": ("term-label join per ionisation stage"
                     + (", energy fallback ENABLED" if energy_fallback
                        else " (energy fallback disabled)")),
        "energy_fallback_used": energy_fallback,
        "energy_tolerance_cm-1": energy_tol_cm if energy_fallback else None,
        "energy_fallback_note": (
            "DISABLED by default. Controlled against the label route where both "
            "answer, the energy route names a DIFFERENT level 18% of the time for "
            "Cr I and 35% for Cr II, because this atom's energies are term averages "
            "and neighbouring terms sit closer together than the fine structure they "
            "average over. No tolerance/margin combination tested exceeded ~87% "
            "agreement (RYA-818)."),
        "wavelength_matching_used": False,
        "multiplicity_matching_used": False,
        "multiplicity_note": (
            "match_multiplicity is DISABLED deliberately. For a term-resolved atom "
            "the line's 2J+1 is one fine-structure component's weight while the "
            "level's g is the whole term's sum; requiring equality would reject "
            "~95% of correct matches (RYA-818)."),
    }
    if getattr(resolution, "is_term_resolved", False):
        prov["approximation"] = (
            "TERM-SHARED DEPARTURE COEFFICIENTS. This atom resolves LS terms, not "
            "fine-structure levels, so every component of a term is assigned that "
            "term's departure coefficient. This is the correct reading of a "
            "term-resolved atom, but it IS an approximation and it is not visible "
            "in the emitted linelist.")
    return prov
