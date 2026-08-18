"""Which gf rung is a PRODUCT entitled to? — RYA-855.

`pipeline.error_budget` has carried a three-rung oscillator-strength ladder since
RYA-850:

    rung 1  `gf scale (UNGRADED)`      0.17    Kurucz semi-empirical (RYA-161)
    rung 2  `gf scale (NIST-graded)`   0.041   a BOUND -- the worst grade we accept
    rung 3  `gf scale (cited lab)`     measured  the pool's OWN published per-line sigmas

`error_budget.build()` can reach all three. `scripts/derive_band_products.py` could reach
only the first: it passed `gf_graded=False` unconditionally at BOTH call sites and never
passed a cited sigma at all, so every band product it emitted charged the ungraded 0.17
regardless of what its lines actually were. RYA-850's promotion lives in the REPORTING
layer, a separate script, so regenerating any band product silently reverted to rung 1.

WHY THIS IS A MODULE AND NOT TWO `if` STATEMENTS
------------------------------------------------
Because the defect being fixed is a SECOND DECLARATION, and fixing it by writing the
decision twice would reproduce it in a form that is harder to see. RYA-845 (the
pseudo-continuum counted twice) and RYA-847 (the constraint decider bypassed on one of
two routes) are the same shape: a rule that lives at its call sites drifts between them,
and each copy is internally consistent while the pair is wrong. So the rung is decided
HERE, once, and both routes hand the answer straight to `error_budget.build()` via
`GfRung.budget_kwargs()` -- neither route is able to state a rung of its own.

THE DECISION, AND THE THREE REFUSALS IN IT
------------------------------------------
The rung is a property of THE LINES THAT ENTERED THE AGGREGATE -- not of the band, not
of the element, and not of the best line in the pool:

1. **A MIXED POOL IS UNGRADED.** A pool is graded only if EVERY line in it is graded.
   Letting a pool claim rung 2/3 on the strength of a subset would attribute a
   laboratory pedigree to lines that do not have one, and the resulting bar would
   describe a pool nobody measured. This is the rule that decides every real case in
   the Fe matrix: the VIS 1D-LTE pool has 5 primary-lab lines out of 152 and the
   red-optical one 12 out of 101, so both land on rung 1 -- correctly, and for a stated
   reason instead of by hardcode. Measured, not assumed: see
   `data/results/rya855/rya855_rung_by_cell.csv`.

2. **A SPECIES WITH NO PRIMARY-LAB TABLE IS UNGRADED, AND IS NEVER ASKED.**
   `pipeline.gf_grades` is Fe I: its lab table is Fe I and `canonical_fe1()` filters to
   `species == 'Fe I'`. Neither checks the species of the line handed to it, so grading
   an Al or Fe II line through it would match Fe I laboratory rows on wavelength and
   excitation potential alone and MANUFACTURE a graded pool out of a coincidence. The
   species gate is therefore here, in front of the call, not inside it.

3. **A LINE WE CANNOT PRICE IS UNGRADEABLE, NOT UNGRADED-BY-DEFAULT.** A grade describes
   the log gf THE POOL ACTUALLY USED (RYA-799), so the value is read from the loaded line
   list -- the same object the inversion and the flux fit ran on, canonical-gf
   substitutions included. A line absent from that list has no stateable gf, so it cannot
   be graded and it forces the whole pool to rung 1. Counted and reported, never dropped.

WHY RUNG 3 CAN WIDEN THE BAR
----------------------------
For the Fe I pools the cited laboratory sigma is 0.052-0.060 dex, LARGER than the 0.041
bound (RYA-850). The bound was optimistic. Nothing here clamps to it: clamping would turn
a measurement back into the assumption it supersedes, and would make a genuinely grade-A
pool report a bar it does not deserve. Rung 3 is not "the smaller number" -- it is the
MEASURED one, and it moves either way.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from collections import Counter

import numpy as np
import pandas as pd

from pipeline import gf_grades
from pipeline.error_budget import GRADED_GF_SYSTEMATIC_DEX

#: A pool whose cited sigmas cover only PART of it is not described by their RMS: the
#: unmatched lines would silently inherit the matched ones' uncertainty. Below this the
#: cited term is refused and the graded BOUND stands, which is the honest fallback.
#: Declared here and imported by `scripts/rya850_graded_products.py` rather than written
#: in both -- a threshold with two homes is the RYA-845 defect shape.
CITED_COVERAGE_MIN = 0.90

#: Species for which a PRIMARY LABORATORY gf table exists in-repo, i.e. the species
#: `pipeline.gf_grades` is actually about. Everything else is rung 1 by construction and
#: is never handed to the grader -- see refusal 2 in the module docstring.
#:
#: 🔴 THIS IS NOT A TODO LIST. Adding a species here without adding its lab table would
#: make `grade_line` referee it against Fe I laboratory rows.
LAB_GRADED_SPECIES: frozenset[tuple[str, str]] = frozenset({("Fe", "I")})

#: Wavelength tolerance when resolving a measured line back to its row in the loaded
#: line list, and the reason it is TIGHT.
#:
#: On the synthesis route the pool is keyed at the list's own wavelength, so this is a
#: rounding tolerance and nothing more. On the EW route it is not: those wavelengths come
#: from the MEASUREMENT, and two catalogues quoting one transition differ by up to ~0.03 A
#: (RYA-704). Measured on the VIS 1D-LTE pool: 16 of 152 lines do not resolve here, and
#: for most of them the nearest Fe I row sits 0.006-0.02 A away.
#:
#: ⚠️ WIDENING IT DOES NOT FIX THAT — it converts "absent" into "ambiguous". At 0.02 A
#: several of those lines have TWO Fe I rows straddling them (5620.3945 has rows at
#: +0.0055 and -0.0175), so a wider window buys a choice, not an identification. The
#: second key that would settle it is the excitation potential, which `gf_grades` uses and
#: which the EW route's per-line artifact does not carry.
#:
#: So an unresolvable line stays UNGRADEABLE and forces rung 1. That is the conservative
#: direction and it flips no answer today — every EW pool is mixed several times over —
#: but it is a real ceiling: a pool that was otherwise entirely lab-gf would be held at
#: rung 1 by lines nobody can identify. Recorded in `rya855_summary.json` under `caveats`.
LINELIST_MATCH_TOL_A = 0.005

#: RYA-871 — the tolerance used WHEN THE MEASURED LINE CARRIES ITS EXCITATION POTENTIAL,
#: and the reason it is four times wider without being four times looser.
#:
#: The paragraph above is still true of a wavelength-only match and is why this constant
#: is a SEPARATE one rather than a widening of it. What changed is that the EW route now
#: carries `ep_eV` (RYA-871), so a wider window can be disambiguated instead of guessed:
#: two REAL transitions at one wavelength differ by whole eV (RYA-855's 3125.651/3125.683
#: pair sits at 0.990 and 2.404 eV), so the EP separates exactly the cases the wider
#: window creates.
#:
#: 🔴 MEASURED, NOT CHOSEN — `data/results/rya871/rya871_ep_resolution_probe.csv`, scored
#: on all 35 banked EW-route cells against the control that matters: does a variant
#: re-identify a line the current rule ALREADY identified? Compared by which ROW it lands
#: on, because swapping one identification for another scores as a tie on a count.
#:
#:     tol    EP   unique  absent  ambiguous   cells that RE-IDENTIFY something
#:     0.005  no     1601     163         40   0     <- today
#:     0.005  yes    1641     163          0   0
#:     0.010  no     1619     115         70   8     <- widening ALONE already breaks it
#:     0.010  yes    1679     125          0   0
#:     0.020  yes    1764      23         17   0     <- here
#:     0.030  yes    1764      23         17   0     <- IDENTICAL: a plateau, not a knob
#:     0.050  yes    1755      23         26   9     <- and it ends
#:     0.020  no     1556      12        236  25
#:
#: The answer is FLAT across 0.020-0.030 and regresses on both sides, so this is not a
#: tuned threshold: any value inside the plateau gives the same identifications, and
#: nothing downstream changes if it moves within it (RYA-161/847). On the ticket's own
#: cell — VIS Fe I 1D-LTE, 152 lines — it takes 136 resolved to **148**, leaving 3 absent
#: and 1 genuinely (wavelength AND EP) degenerate, with 136 of 136 prior identifications
#: unchanged.
LINELIST_MATCH_TOL_EPKEY_A = 0.020

#: How closely the carried EP must agree with the line list's. The accounting table
#: rounds to 4 dp and the list carries full precision, so this is a ROUNDING tolerance,
#: not a physical one. Nothing depends on where inside it the cut sits: distinct
#: transitions at one wavelength differ by ~1 eV, which is 200x this.
EP_MATCH_TOL_EV = 0.005

_ROMAN = {"I": 1, "II": 2, "III": 3, "IV": 4, "V": 5, "VI": 6}


@dataclass(frozen=True)
class GfRung:
    """The gf rung one product is entitled to, and the evidence for it."""
    rung: int                        # 1 ungraded | 2 graded bound | 3 cited lab sigma
    gf_graded: bool
    cited_sigma_dex: float | None
    cited_source: str
    n_lines: int
    n_graded: int
    n_unresolved: int                # lines absent from the loaded line list
    coverage: float                  # fraction of the pool carrying a cited sigma
    grade_counts: dict[str, int] = field(default_factory=dict)
    reason: str = ""

    @property
    def term_name(self) -> str:
        return {1: "gf scale (UNGRADED)", 2: "gf scale (NIST-graded)",
                3: "gf scale (cited lab)"}[self.rung]

    def budget_kwargs(self) -> dict:
        """Exactly the gf arguments `error_budget.build()` must be handed.

        Returned as a mapping rather than as separate fields so a caller cannot pass the
        graded flag and forget the cited sigma -- the pair is the decision, and splitting
        it is how one call site drifts from the other.
        """
        kw: dict = {"gf_graded": self.gf_graded}
        if self.cited_sigma_dex is not None:
            kw["cited_gf_sigma_dex"] = float(self.cited_sigma_dex)
            kw["cited_gf_source"] = self.cited_source
        return kw

    def describe(self) -> str:
        return (f"gf rung {self.rung} ({self.term_name}): {self.reason}")


def _species_label(element: str, ion: str) -> str:
    """The line list's own species spelling, e.g. ('Fe', 'I') -> 'Fe 1' (RYA-759)."""
    n = _ROMAN.get(str(ion).strip().upper())
    if n is None:
        raise ValueError(f"unrecognised ionisation stage {ion!r}")
    return f"{str(element).strip()} {n}"


def linelist_frame(linelist) -> pd.DataFrame:
    """The loaded iSpec line list as (species, wavelength, EP, log gf).

    `loggf` here is the value the run ACTUALLY used -- `_load_synth_resources` applies the
    canonical-gf substitutions before handing the list over, and RYA-799's whole point is
    that a grade describes the number that was used and not the best number available.
    Reading it from anywhere else would re-derive it and could drift.
    """
    names = getattr(linelist, "dtype", None)
    names = names.names if names is not None else None
    if names is None:                       # already a frame
        d = pd.DataFrame(linelist)
        names = list(d.columns)
    else:
        d = None
    def col(key):
        return np.asarray(linelist[key]) if d is None else d[key].to_numpy()
    w = (col("wave_A").astype(float) if "wave_A" in names
         else col("wave_nm").astype(float) * 10.0)
    return pd.DataFrame({
        "species": [str(x).strip() for x in col("element")],
        "wavelength_air_A": w,
        "ep_eV": col("lower_state_eV").astype(float),
        "log_gf": col("loggf").astype(float)})


def resolve_lines(element: str, ion: str, wavelengths, linelist,
                  measured_ep_eV=None) -> pd.DataFrame:
    """Attach the EP and the log gf THE POOL USED to each measured wavelength.

    A wavelength that resolves to more than one row of the SAME species inside the
    tolerance is returned unresolved rather than settled by `iloc[0]`: picking the first
    match is exactly how RYA-853 manufactured 12-dex "defects", and here it would let a
    pool inherit a pedigree from the wrong transition.

    RYA-871 — `measured_ep_eV`, when supplied, is the excitation potential the MEASURED
    line carries (`LineMeasurement.ep_eV`). A line that has one is matched on wavelength
    AND EP inside `LINELIST_MATCH_TOL_EPKEY_A`; a line that has none keeps the narrow
    wavelength-only rule at `LINELIST_MATCH_TOL_A`.

    🔴 THE TWO TOLERANCES ARE NOT A PREFERENCE, THEY ARE THE SAME RULE. The window may be
    wide only because the EP can settle what falls into it; widening it for a line with
    no EP would buy a choice rather than an identification, which is precisely what
    RYA-855 refused and what the RYA-871 probe measured breaking (7 of 136 already-
    identified lines change row at 0.020 A with no EP key). So the tolerance travels with
    the key, per line, and a route that does not carry an EP is not silently widened.
    """
    ll = linelist_frame(linelist)
    ll = ll[ll.species == _species_label(element, ion)]
    lw = ll.wavelength_air_A.to_numpy()
    lep = ll.ep_eV.to_numpy()
    # Materialised once: `wavelengths` may be a generator, and consuming it to size the
    # EP list would leave the loop below with nothing to iterate.
    ws = [float(x) for x in wavelengths]
    eps = ([None] * len(ws) if measured_ep_eV is None
           else [None if x is None or not np.isfinite(float(x)) else float(x)
                 for x in measured_ep_eV])
    if len(eps) != len(ws):
        raise ValueError(
            f"{len(ws)} wavelengths but {len(eps)} excitation potentials — these are "
            f"per-line parallel arrays and a length mismatch would silently key a line "
            f"on its neighbour's EP (RYA-871)")
    rows = []
    for i, w in enumerate(ws):
        ep = eps[i]
        tol = LINELIST_MATCH_TOL_A if ep is None else LINELIST_MATCH_TOL_EPKEY_A
        m = np.abs(lw - w) <= tol
        if ep is not None:
            m = m & (np.abs(lep - ep) <= EP_MATCH_TOL_EV)
        n = int(m.sum())
        key = ("wavelength alone" if ep is None
               else f"wavelength+EP ({ep:.4f} eV)")
        if n == 1:
            r = ll.iloc[int(np.flatnonzero(m)[0])]
            rows.append({"wavelength_air_A": w, "ep_eV": float(r.ep_eV),
                         "log_gf": float(r.log_gf), "resolved": True,
                         "unresolved_why": ""})
        else:
            rows.append({"wavelength_air_A": w, "ep_eV": np.nan, "log_gf": np.nan,
                         "resolved": False,
                         "unresolved_why": (
                             f"absent: no same-species row within {tol} A on {key}"
                             if n == 0 else
                             f"ambiguous: {n} same-species rows within {tol} A on {key}")})
    return pd.DataFrame(rows)


def decide(element: str, ion: str, lines: pd.DataFrame) -> GfRung:
    """The rung this pool is entitled to. `lines` needs wavelength_air_A / ep_eV / log_gf.

    Rows may carry `resolved=False` (see `resolve_lines`); those are counted as
    ungradeable and force rung 1, because a line whose gf we cannot state cannot be
    said to be graded.
    """
    n = int(len(lines))
    species = f"{element} {ion}"
    if n == 0:
        return GfRung(1, False, None, "", 0, 0, 0, 0.0, {},
                      "empty pool — no line to grade, so the ungraded systematic stands")

    if (str(element), str(ion)) not in LAB_GRADED_SPECIES:
        return GfRung(
            1, False, None, "", n, 0, 0, 0.0, {},
            f"no primary-laboratory gf table exists for {species} — the graded ladder is "
            f"Fe I only (RYA-799/824/836), and grading {species} through it would referee "
            f"it against Fe I lab rows on wavelength and EP alone")

    resolved = lines.get("resolved")
    unresolved = int((~resolved.astype(bool)).sum()) if resolved is not None else 0

    grades, sigmas, cites = [], [], []
    for r in lines.itertuples():
        if resolved is not None and not bool(getattr(r, "resolved")):
            grades.append("UNRESOLVED")
            continue
        v = gf_grades.grade_line(float(r.wavelength_air_A), float(r.ep_eV),
                                 float(r.log_gf))
        grades.append(v.gf_grade)
        if v.is_graded:
            if np.isfinite(v.gf_sigma_dex):
                sigmas.append(float(v.gf_sigma_dex))
            cites.append(v.gf_grade_source)
    counts = dict(sorted(Counter(grades).items()))
    n_graded = counts.get(gf_grades.GRADE_LAB, 0)

    if n_graded < n:
        other = ", ".join(f"{k} x{v}" for k, v in counts.items()
                          if k != gf_grades.GRADE_LAB)
        return GfRung(
            1, False, None, "", n, n_graded, unresolved, 0.0, counts,
            f"MIXED POOL: {n_graded} of {n} {species} lines are {gf_grades.GRADE_LAB} "
            f"(primary laboratory gf); the rest are {other}. A pool is graded only if "
            f"every line in it is — inheriting the best rung a subset carries would "
            f"attribute a laboratory pedigree to lines that do not have one")

    coverage = len(sigmas) / n
    if coverage < CITED_COVERAGE_MIN or not sigmas:
        return GfRung(
            2, True, None, "", n, n_graded, unresolved, coverage, counts,
            f"every one of the {n} {species} lines is {gf_grades.GRADE_LAB}, but only "
            f"{len(sigmas)} ({coverage:.0%} < {CITED_COVERAGE_MIN:.0%}) carry a citable "
            f"per-line sigma — an RMS over part of a pool does not describe the pool, so "
            f"the graded BOUND stands")

    sig = float(np.sqrt(np.mean(np.asarray(sigmas, dtype=float) ** 2)))
    # NAME THE PAPERS. `error_budget.cited_gf_term` refuses an unsourced sigma, and the
    # citation is built from the lines actually in this pool rather than from a constant
    # list, so it cannot name a paper the pool did not use.
    source = "; ".join(sorted(set(cites)))
    return GfRung(
        3, True, sig, source, n, n_graded, unresolved, coverage, counts,
        f"every one of the {n} {species} lines is {gf_grades.GRADE_LAB} and "
        f"{coverage:.0%} carry a published per-line sigma; RMS {sig:.4f} dex "
        f"({sig / GRADED_GF_SYSTEMATIC_DEX:.2f}x the generic bound "
        f"{GRADED_GF_SYSTEMATIC_DEX}) MEASURES the same quantity the bound only "
        f"estimates, so it supersedes it — larger or smaller (RYA-850)")


def for_lines(element: str, ion: str, measurements, *, linelist) -> GfRung:
    """The rung for a list of `LineMeasurement`, counting only the ones IN the aggregate.

    The membership rule is `build_product`'s own -- `in_aggregate` AND an abundance --
    written once here so a line the product excluded can never price the product's bar.
    A quarantined line is not part of the measurement, so its pedigree is not part of the
    measurement's uncertainty either.
    """
    used = [l for l in measurements if l.in_aggregate and l.abundance is not None]
    # RYA-871 — hand the resolver the EP each measurement carries. `getattr` rather than
    # attribute access so a caller passing some other line-like object still works, and
    # None (the default) keeps the narrow wavelength-only rule for it.
    lines = resolve_lines(element, ion, [l.wavelength_air_A for l in used], linelist,
                          measured_ep_eV=[getattr(l, "ep_eV", None) for l in used])
    return decide(element, ion, lines)


def for_product(product, *, linelist) -> GfRung:
    """The rung for one `pipeline.band_products.Product`. Delegates to `for_lines`."""
    return for_lines(product.element, product.ion, product.lines, linelist=linelist)
