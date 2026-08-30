"""RYA-1119 (the Bride, M4) — read a MULTI-format model atom, restrict it to a line set,
and emit a Lightweaver `AtomicModel`.

WHY THIS EXISTS. Lightweaver ingests CRTAF (`lightweaver/crtaf.py::from_crtaf`) and its own
`AtomicModel`; it does NOT read MULTI atoms — `lightweaver/multi.py` reads MULTI
*atmospheres*. `crtaf-py`, the reference implementation, ships only the Lightweaver round
trip (`from_lightweaver` / `to_lightweaver`). So the front half of MULTI -> Lightweaver is
the one missing hop, and this module is that hop. It targets `AtomicModel` directly rather
than going via CRTAF, because CRTAF buys nothing here and adds a dependency.

RYA-1137 established the two facts that make this worth doing:

  * every rate type `atom.fe607a` actually uses has a Lightweaver counterpart --
    CE (14,001), CH (15,090), CI (548), CH0 (344) -> `CE`, `CH`, `CI`,
    `ChargeExchangeNeutralH`. No missing physics, so this is parsing, not research.
  * the FULL atom is not runnable (607 levels / 12,635 lines -> hours-to-days per column),
    but the ASSIGNED line set is: the 40 AGSS21 Fe I lines touch only **58 levels / 298
    transitions**, and a measured proxy of that size (`MgI_66_atom`, 66 levels / 315 lines)
    converges in 32.15 s/column -- ~57 core-hours for a 6400-column cube.

🔴 RESTRICTION IS NOT FREE, AND THIS MODULE DOES NOT PRETEND OTHERWISE. Dropping 490 levels
changes the statistical-equilibrium solution: the omitted levels carry population and
collisional/radiative coupling. Selecting from a vetted atom is not building one, but the
RESULT is not vetted either. `restrict()` records what it dropped in `RestrictionReport` so
the eventual validation -- reproduce fe607a's own departures on the 40 AGSS21 lines, and
land A(Fe) on 7.46 -- is checkable rather than assumed. A restricted atom that converges to
a different departure pattern is a wrong atom that runs (RYA-1118's striding lesson: every
value/range/finite check passed on a transposed cube; only STRUCTURE caught it).

FORMAT, read off the file rather than from documentation:

    line 1        element symbol
    line 2        abundance, atomic weight
    line 3        nlevel nline ncont nfix
    levels        E[cm^-1]  g  'label'  stage        <- stage is 1-BASED (1 = neutral)
    lines         j i f nq qmax iw iqo gRad gVdw gStark PROFILE
    continua      j i alpha0 nlamb wavelength_dep     <- then nlamb (lambda, alpha) pairs
    GENCOL / TEMP ntemp T1..Tn
    collisions    <KEYWORD> \n j i rate1..rate_ntemp

Comment lines start with '*' and are stripped before any indexing.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field

import numpy as np

#: MULTI writes level energies in cm^-1 and so does Lightweaver's `AtomicLevel.E`, so no
#: conversion happens anywhere in this module. Kept as a named constant only for the eV
#: views used in reports, so a reader never has to wonder which unit a number is in.
CM1_TO_EV = 1.0 / 8065.54429

#: Collision keywords this module understands. Anything else is carried through as
#: UNMAPPED rather than dropped silently -- an unrecognised rate is a physics change.
_KNOWN_COLLISIONS = {"CE", "CI", "CH", "CH0", "CP", "OHM", "CR"}


@dataclass
class MultiLevel:
    E_cm: float
    g: float
    label: str
    stage: int          # 1-based, as written in the file
    index: int          # 1-based, as referenced by lines/continua/collisions

    @property
    def E_eV(self) -> float:
        return self.E_cm * CM1_TO_EV


@dataclass
class MultiLine:
    j: int              # upper, 1-based
    i: int              # lower, 1-based
    f: float
    nq: int
    qmax: float
    gRad: float
    gVdw: float
    gStark: float
    profile: str


@dataclass
class MultiContinuum:
    j: int
    i: int
    alpha0: float
    nlamb: int
    wavelength: np.ndarray
    alpha: np.ndarray


@dataclass
class MultiCollision:
    kind: str
    j: int
    i: int
    rates: np.ndarray


@dataclass
class MultiAtom:
    symbol: str
    abundance: float
    weight: float
    levels: list[MultiLevel]
    lines: list[MultiLine]
    continua: list[MultiContinuum]
    temperature: np.ndarray
    collisions: list[MultiCollision]

    def level(self, idx1: int) -> MultiLevel:
        return self.levels[idx1 - 1]


@dataclass
class RestrictionReport:
    """What a restriction kept and what it threw away — so the cost is visible.

    `restrict()` is a physics change, not a filter. This is the record that makes the
    change auditable instead of implicit.
    """
    kept_levels: list[int]
    n_levels_before: int
    n_lines_before: int
    n_continua_before: int
    n_collisions_before: int
    n_levels_after: int
    n_lines_after: int
    n_continua_after: int
    n_collisions_after: int
    dropped_levels: list[int] = field(default_factory=list)
    unmapped_collision_kinds: dict = field(default_factory=dict)

    def describe(self) -> str:
        L = ["RESTRICTION — a physics change, recorded so it can be checked:",
             f"  levels      {self.n_levels_before:>6} -> {self.n_levels_after:<6}"
             f" ({self.n_levels_before - self.n_levels_after} dropped)",
             f"  lines       {self.n_lines_before:>6} -> {self.n_lines_after}",
             f"  continua    {self.n_continua_before:>6} -> {self.n_continua_after}",
             f"  collisions  {self.n_collisions_before:>6} -> {self.n_collisions_after}"]
        if self.unmapped_collision_kinds:
            L.append(f"  ⚠️ collision kinds with no Lightweaver mapping: "
                     f"{self.unmapped_collision_kinds}")
        L.append("  ⚠️ dropped levels carried population and coupling. The restricted atom "
                 "is NOT validated by construction — it must reproduce the source atom's "
                 "departures on the target lines before any number from it is trusted.")
        return "\n".join(L)


def _strip(path) -> list[str]:
    with open(path, errors="replace") as fh:
        return [l for l in fh.read().splitlines() if not l.lstrip().startswith("*")]


def read_multi_atom(path) -> MultiAtom:
    """Parse a MULTI model atom. Raises rather than guessing on any structural surprise."""
    rows = _strip(path)
    symbol = rows[0].split()[0]
    abundance, weight = (float(x) for x in rows[1].split()[:2])
    nlev, nlin, ncont, _nfix = (int(x) for x in rows[2].split()[:4])

    levels = []
    for k, line in enumerate(rows[3:3 + nlev], start=1):
        m = re.match(r"\s*([-\d.Ee+]+)\s+([-\d.Ee+]+)\s+'(.*?)'\s+(\d+)", line)
        if not m:
            raise ValueError(f"{path}: level {k} does not parse: {line!r}")
        levels.append(MultiLevel(float(m.group(1)), float(m.group(2)),
                                 m.group(3).strip(), int(m.group(4)), k))
    if len(levels) != nlev:
        raise ValueError(f"{path}: header declares {nlev} levels, parsed {len(levels)}")

    lines = []
    for line in rows[3 + nlev:3 + nlev + nlin]:
        t = line.split()
        lines.append(MultiLine(int(t[0]), int(t[1]), float(t[2]), int(t[3]), float(t[4]),
                               float(t[7]), float(t[8]), float(t[9]),
                               t[10] if len(t) > 10 else "VOIGT"))
    if len(lines) != nlin:
        raise ValueError(f"{path}: header declares {nlin} lines, parsed {len(lines)}")

    # 🔴 CONTINUA ARE VARIABLE LENGTH — header then NLAMB (lambda, alpha) pairs — so this
    # walks sequentially. Indexing them at a fixed offset silently reads wavelength pairs
    # as continuum headers and every subsequent block is misaligned.
    i = 3 + nlev + nlin
    continua = []
    for _ in range(ncont):
        t = rows[i].split()
        j_, i_, alpha0, nlamb = int(t[0]), int(t[1]), float(t[2]), int(t[3])
        pairs = np.array([[float(x) for x in rows[i + 1 + k].split()[:2]]
                          for k in range(nlamb)], dtype=float)
        continua.append(MultiContinuum(j_, i_, alpha0, nlamb, pairs[:, 0], pairs[:, 1]))
        i += 1 + nlamb

    while i < len(rows) and rows[i].split()[0].upper() != "GENCOL":
        i += 1
    if i >= len(rows):
        raise ValueError(f"{path}: no GENCOL block — collisional rates are missing")
    i += 1
    if rows[i].split()[0].upper() != "TEMP":
        raise ValueError(f"{path}: GENCOL is not followed by TEMP, got {rows[i]!r}")
    t = rows[i + 1].split()
    ntemp = int(t[0])
    temps = [float(x) for x in t[1:]]
    i += 2
    while len(temps) < ntemp:                     # the grid may wrap onto further rows
        temps += [float(x) for x in rows[i].split()]
        i += 1
    temperature = np.array(temps[:ntemp], dtype=float)

    collisions = []
    while i < len(rows):
        w = rows[i].split()
        if not w:
            i += 1
            continue
        kind = w[0].upper()
        if kind == "END":
            break
        if i + 1 >= len(rows):
            break
        vals = rows[i + 1].split()
        try:
            j_, i_ = int(vals[0]), int(vals[1])
        except (ValueError, IndexError):
            i += 1
            continue
        rates = np.array([float(x) for x in vals[2:2 + ntemp]], dtype=float)
        collisions.append(MultiCollision(kind, j_, i_, rates))
        i += 2

    return MultiAtom(symbol, abundance, weight, levels, lines, continua,
                     temperature, collisions)


def levels_touched(atom: MultiAtom, lines_elo_eup_eV, tol_eV: float = 0.002):
    """1-based level indices touched by (Elo, Eup) pairs, plus the lines that did not match.

    Matching is on energy because a MULTI atom's labels are free text. `tol_eV` is small on
    purpose: fe607a is J-resolved, so an exact fine-structure match is the right test here —
    unlike the bundled RH atoms, whose term-averaged levels make an energy match meaningless
    (that mistake cost a wrong 'Fe is dead' verdict on RYA-1137).
    """
    import bisect
    neutral = [(l.index, l.E_eV) for l in atom.levels if l.stage == 1]
    Es = sorted(E for _, E in neutral)
    by_E = {E: idx for idx, E in neutral}

    def nearest(E):
        k = bisect.bisect_left(Es, E)
        best = None
        for c in (k - 1, k):
            if 0 <= c < len(Es) and abs(Es[c] - E) <= tol_eV:
                if best is None or abs(Es[c] - E) < abs(best - E):
                    best = Es[c]
        return by_E[best] if best is not None else None

    touched, missed = set(), []
    for elo, eup in lines_elo_eup_eV:
        a, b = nearest(elo), nearest(eup)
        if a is None or b is None:
            missed.append((elo, eup))
        else:
            touched.update((a, b))
    return sorted(touched), missed


def restrict(atom: MultiAtom, keep: list[int], *, keep_ion_ground: bool = True):
    """Keep `keep` (1-based level indices) and everything that connects only those.

    🔴 The ion ground state is added unless explicitly refused. Without it there is no
    continuum and no ionisation balance, and the atom silently becomes a closed neutral
    system — which converges perfectly happily to the wrong answer.
    """
    keep = set(keep)
    if keep_ion_ground:
        higher = [l for l in atom.levels if l.stage > 1]
        if not higher:
            raise ValueError("no ionised stage in this atom — cannot close the continuum")
        top = min(l.stage for l in higher)
        ground = min((l for l in higher if l.stage == top), key=lambda l: l.E_cm)
        keep.add(ground.index)

    kept = sorted(keep)
    remap = {old: new for new, old in enumerate(kept)}       # -> 0-based, Lightweaver
    lines = [l for l in atom.lines if l.j in keep and l.i in keep]
    continua = [c for c in atom.continua if c.j in keep and c.i in keep]
    cols, unmapped = [], {}
    for c in atom.collisions:
        if c.j in keep and c.i in keep:
            cols.append(c)
            if c.kind not in _KNOWN_COLLISIONS:
                unmapped[c.kind] = unmapped.get(c.kind, 0) + 1

    report = RestrictionReport(
        kept_levels=kept,
        n_levels_before=len(atom.levels), n_lines_before=len(atom.lines),
        n_continua_before=len(atom.continua), n_collisions_before=len(atom.collisions),
        n_levels_after=len(kept), n_lines_after=len(lines),
        n_continua_after=len(continua), n_collisions_after=len(cols),
        dropped_levels=[l.index for l in atom.levels if l.index not in keep],
        unmapped_collision_kinds=unmapped)

    sub = MultiAtom(atom.symbol, atom.abundance, atom.weight,
                    [atom.level(k) for k in kept], lines, continua,
                    atom.temperature, cols)
    return sub, remap, report


def to_lightweaver(atom: MultiAtom, remap: dict | None = None):
    """Build a Lightweaver `AtomicModel`. Import is local so parsing needs no Lightweaver."""
    from lightweaver.atomic_model import (AtomicLevel, AtomicModel, ExplicitContinuum,
                                          LineType, VoigtLine)
    from lightweaver.atomic_table import Element, PeriodicTable
    from lightweaver.broadening import (LineBroadening, RadiativeBroadening,
                                        VdwUnsold, QuadraticStarkBroadening)
    from lightweaver.collisional_rates import CE, CH, CI, ChargeExchangeNeutralH

    if remap is None:
        remap = {l.index: n for n, l in enumerate(atom.levels)}

    # 🔴 STAGE IS 1-BASED IN MULTI AND 0-BASED IN LIGHTWEAVER. An off-by-one here does not
    # crash: it produces a perfectly convergent atom with the wrong ionisation balance.
    levels = [AtomicLevel(E=l.E_cm, g=l.g, label=l.label, stage=l.stage - 1, J=None)
              for l in atom.levels]

    lines = []
    for L in atom.lines:
        broad = LineBroadening(
            natural=[RadiativeBroadening(L.gRad)] if L.gRad > 0 else [],
            elastic=([VdwUnsold(vals=[1.0, 1.0])] if L.gVdw != 0 else [])
                    + ([QuadraticStarkBroadening(1.0)] if L.gStark != 0 else []))
        lines.append(VoigtLine(j=remap[L.j], i=remap[L.i], f=L.f, type=LineType.CRD,
                               quadrature=_quadrature(L), broadening=broad))

    # 🔴 TWO UNIT/ORDER TRAPS IN ONE LINE, both of which produce a well-formed atom.
    #
    #  1. MULTI writes continuum grids in DESCENDING wavelength; Lightweaver requires
    #     ascending and raises. Sort BOTH arrays by the SAME permutation — reversing only
    #     the wavelengths would pair every cross-section with the wrong wavelength and
    #     invert the photoionisation edge without any error.
    #  2. UNITS, both axes, and Lightweaver's own docstring states them:
    #     `wavelengthGrid [nm]`, `alphaGrid [m2]`. MULTI writes ANGSTROM and CM^2.
    #       - Angstrom -> nm (x0.1). Getting this wrong put every grid point ~10x above
    #         `lambdaEdge`, the grid filtered to EMPTY, and it surfaced far away as
    #         `IndexError: index -1 is out of bounds` inside `compute_wavelength_grid`.
    #       - cm^2 -> m^2 (x1e-4). This one does NOT raise: it makes every photoionisation
    #         cross-section 10,000x too strong and still converges. Caught only by
    #         comparing magnitudes against the bundled atoms (theirs 1e-24..1e-21 m^2,
    #         mine 1e-20..1e-16 cm^2 — four orders apart, which IS the cm^2/m^2 factor).
    CM2_TO_M2 = 1e-4
    continua = []
    for c in atom.continua:
        order = np.argsort(c.wavelength)
        continua.append(ExplicitContinuum(
            j=remap[c.j], i=remap[c.i],
            wavelengthGrid=list(c.wavelength[order] * 0.1),
            alphaGrid=list(c.alpha[order] * CM2_TO_M2)))

    # 🔴 CGS -> SI. Every Lightweaver rate class documents `m^3`; MULTI is CGS (`cm^3`).
    # Passing cm^3 makes every collisional rate 1e6 too large, which does not raise and
    # does not diverge — it THERMALISES the atom completely and returns b = 1.0000 at
    # every level and every depth. A converged LTE answer wearing an NLTE label, i.e.
    # exactly the failure this work exists to avoid. It was caught by checking the
    # departure STRUCTURE, never by convergence (RYA-1118's lesson, again).
    CM3_TO_M3 = 1e-6

    # ⚠️ CH0 is CHARGE EXCHANGE with neutral H (downward only), NOT the same thing as CH
    # (collisions with neutral H, both directions). Mapping CH0 -> CH adds an upward rate
    # the physics does not have.
    kind_map = {"CE": CE, "CI": CI, "CH": CH, "CH0": ChargeExchangeNeutralH}
    collisions = []
    for c in atom.collisions:
        cls = kind_map.get(c.kind)
        if cls is None:
            continue
        collisions.append(cls(j=remap[c.j], i=remap[c.i],
                              temperature=list(atom.temperature),
                              rates=list(c.rates * CM3_TO_M3)))

    Z = PeriodicTable[atom.symbol.capitalize()].Z
    return AtomicModel(element=Element(Z=Z), levels=levels, lines=lines,
                       continua=continua, collisions=collisions)


def _quadrature(L: MultiLine):
    from lightweaver.atomic_model import LinearCoreExpWings
    return LinearCoreExpWings(qCore=max(L.qmax / 4.0, 1.0), qWing=max(L.qmax, 10.0),
                              Nlambda=max(L.nq, 5))
