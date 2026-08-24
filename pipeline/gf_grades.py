"""Per-line gf accuracy: grade it, or bound it honestly — RYA-799.

Ryan, 2026-08-12 (ratified): *"we still do our best job, and if there is no grade, no way
to tie it to a graded line, then we note the systematic."* A larger bar is a result; a
hidden bar is a defect (RYA-713).

WHAT A gf GRADE HAS TO SURVIVE
------------------------------
Two rules, and the second is the one that is easy to get wrong.

1. **The grade names its subject** (RYA-711). Our values are prefixed `GF-` because they
   describe the ATOMIC DATUM, not our measurement of the line — `MQ-A..D` already means
   measurement quality and is a different claim about a different thing. NIST's own
   letters stay UNPREFIXED wherever they are quoted: renaming someone else's published
   scale misquotes the source, and RYA-711 asserts we never do it.

2. **A grade describes the NUMBER THAT WAS USED, not the best number available.** This is
   the trap this module exists for. `data/linelists/canonical_gf.csv` carries a real
   per-line `loggf_reference` for the IR Fe I lines — but the measured pool was inverted
   on the VALD Kurucz-2014 value, and for a large minority of lines those two gf differ,
   by up to ~0.6 dex. Attaching the good reference's grade to an abundance derived from a
   different gf would be a fabricated pedigree: the number would look graded and would
   not be. So every tie is confirmed against the log gf the pool actually used, and a tie
   that fails that check becomes `SCALE-MISMATCH` — still bounded by the systematic, but
   flagged, because "a better gf exists in-repo and was not used" is actionable and must
   not be silently absorbed.

THE THREE TERMINAL STATES
-------------------------
  `GF-LAB`                      a PRIMARY laboratory measurement (Ruffoni 2014 /
                                Den Hartog 2014 / Belmonte 2017) covers this line AND its
                                value is the one the pool used. sigma = the paper's own
                                per-line uncertainty in dex. This is the only state with
                                a measured bar.
  `systematic:K07`              no graded tie. sigma = the RYA-161 semi-empirical
                                Kurucz systematic. Never blank, never averaged down.
  `systematic:K07/SCALE-MISMATCH`  a graded or primary-referenced gf for this line exists
                                in `canonical_gf.csv`, but it is NOT the value the pool
                                was measured with. Same bar as above — the bar cannot be
                                improved by a number we did not use — plus the flag.

WHY NOT NIST ASD LETTERS
------------------------
They were the intended grade axis and the endpoint is down: `lines1.pl` has returned
`Can't use an undefined value as an ARRAY reference` on every query form since at least
2026-08-11, re-verified 2026-08-15 (HTTP 500 on the exact RYA-760 recipe). Recorded as an
external blocker rather than quietly skipped. It is also the weaker axis: RYA-760 showed
`FMW` IS a NIST compilation and VALD copies it, so an ASD agreement proves only that
nobody fat-fingered a transcription. Only a primary lab measurement can referee a gf.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

_REPO_ROOT = Path(__file__).resolve().parents[1]
CANONICAL_GF = _REPO_ROOT / "data" / "linelists" / "canonical_gf.csv"

#: Per-species primary-laboratory gf tables — RYA-1002.
#:
#: This was a single module constant pointing at the Fe file, which meant `GF-LAB` was
#: STRUCTURALLY UNREACHABLE for every other element: a species with a perfectly good lab
#: measurement could not be graded, because there was nowhere to put its table. Al is the
#: case that surfaced it (Burheim 2023 grades 8 of our Al I lines at 1.5-11%), but the
#: defect was never Al-specific.
#:
#: A table joins this registry by matching the Fe schema — `source, wavelength_air_A,
#: elo_cm1, eup_cm1, elo_eV, eup_eV, loggf, e_loggf_dex`. Extra columns are ignored, so a
#: table may carry its own provenance beside the shared ones. Every source named in a
#: table MUST have a CITATIONS entry: `grade_line` cites the paper by name in the verdict,
#: and an uncited lab value would be an unattributable pedigree.
LAB_TABLES: dict[str, Path] = {
    "Fe I": _REPO_ROOT / "data" / "reference" / "fe_gf_lab" / "fe1_lab_loggf.csv",
    "Al I": _REPO_ROOT / "data" / "reference" / "al_gf_lab" / "al1_lab_loggf.csv",
    # RYA-953: Fe II. The table itself has existed since RYA-945 -- full Den Hartog 2019
    # Table 6, 131 lines, with its own provenance JSON. What was missing was an entry
    # HERE, so `gf_rung.LAB_GRADED_SPECIES` (derived from this dict) refused Fe II and
    # every Fe II pool sat at rung 1 carrying the 0.17 dex Kurucz placeholder, while the
    # measurements sat on disk. Same shape RYA-1002 fixed for Al.
    #
    # ⚠️ 2249.2-4583.8 A. VERIFIED as the SOURCE's own limit, not an ingest truncation:
    # DH19 Table 6 ends there. So Fe II laboratory gf reaches near-UV (12 lines) and the
    # blue half of VIS (9, none above 4583.8), and there is NONE in red-optical, NIR or
    # IR -- for those bands the best available is the NIST-C+ compilation (147 lines in
    # red-optical), which is rung 2, not rung 3. A "Fe II VIS graded" product is a BLUE
    # subset wearing a band name, and must say so.
    "Fe II": _REPO_ROOT / "data" / "reference" / "fe_gf_lab" / "fe2_lab_loggf_dh19.csv",
}
#: How to rebuild each, quoted in the not-found error so the message is actionable.
LAB_REGEN = {
    "Fe I": "python3 scripts/rya799_fetch_fe_gf_lab.py",
    "Al I": "python3 scripts/rya1002_fetch_al_gf_lab.py",
    "Fe II": "python3 scripts/rya945_fetch_fe2_gf_dh19.py",
}
#: The default species. Every pre-RYA-1002 caller passes no species and must keep getting
#: exactly the Fe I behaviour it got before — the generalisation is additive, never a
#: silent re-grade of the Fe pool.
DEFAULT_SPECIES = "Fe I"
#: Back-compat alias. Several scripts import the Fe path directly.
LAB_CSV = LAB_TABLES[DEFAULT_SPECIES]

#: RYA-161. The semi-empirical Kurucz gf systematic carried when nothing better ties.
K07_SYSTEMATIC_DEX = 0.20
#: The tag `canonical_gf.loggf_reference` uses for that semi-empirical source.
K07_TAG = "K07"

GRADE_LAB = "GF-LAB"
GRADE_SYSTEMATIC = "systematic:K07"
GRADE_MISMATCH = "systematic:K07/SCALE-MISMATCH"
#: RYA-822. A value carrying a NIST ASD accuracy class — a citable per-line sigma, but from
#: a COMPILATION, not a primary laboratory measurement. Kept distinct from GF-LAB for that
#: reason (see this module's RYA-822 note).
GRADE_NIST = "GF-NIST"

#: NIST's accuracy ladder, worst-case percentage on the transition probability.
#: ⚠️ ENUMERATED IN FULL, '+' TIERS INCLUDED. Omitting them once put a B+ line (<=7%)
#: BELOW a B line (<=10%) — an inverted ladder that silently demoted the better
#: measurement (RYA-592, found in `curate_nonfe_pools.NIST_GRADE_HIGH`).
NIST_ACC_PCT = {
    "AAA": 0.3, "AA": 1.0, "A+": 2.0, "A": 3.0, "B+": 7.0, "B": 10.0,
    "C+": 18.0, "C": 25.0, "D+": 40.0, "D": 50.0, "E": 100.0,
}


def nist_sigma_dex(grade: str) -> float:
    """A NIST accuracy class as a log-space sigma: dex = log10(1 + pct/100).

    The same percent->dex bridge `error_budget.gf_term` already documents ("NIST grade
    B = 10% on the transition probability = log10(1.10) dex") and the one Belmonte 2017
    Eq. 1 states independently. Returns NaN for a class not in the ladder rather than
    guessing — an unmapped grade must not silently become an accuracy.
    """
    pct = NIST_ACC_PCT.get(str(grade).strip())
    return float(np.log10(1.0 + pct / 100.0)) if pct is not None else float("nan")

#: Match windows. Wavelength is a rounding tolerance, not a search radius — the pool is
#: measured AT the line list's own wavelength. The EP window is the second key: RYA-780
#: ratified matching on wavelength AND excitation potential, because a same-species
#: neighbour does not cancel (RYA-785).
WAVE_TOL_A = 0.02
EP_TOL_EV = 0.02
#: How close the reference's log gf must be to the value the pool used for the reference's
#: accuracy to describe that value. 0.02 dex is below the 2-decimal quantisation of the
#: VALD delivery, so this admits "same number, different rounding" and nothing else.
LOGGF_MATCH_TOL = 0.02

CITATIONS = {
    "Ruffoni2014": ("Ruffoni et al. 2014, MNRAS 441, 3127", "10.1093/mnras/stu780"),
    "DenHartog2014": ("Den Hartog et al. 2014, ApJS 215, 23", "10.1088/0067-0049/215/2/23"),
    "Belmonte2017": ("Belmonte et al. 2017, ApJ 848, 125", "10.3847/1538-4357/aa8cd3"),
    # RYA-1002. Experimental Al I: FTS branching fractions x radiative lifetimes, 12
    # lines at 2-11%. The only primary lab source for Al in the repo.
    "Burheim2023": ("Burheim, Hartman & Nilsson 2023, A&A 672, A197",
                    "10.1051/0004-6361/202245394"),
    # RYA-953. Experimental Fe II: the only primary laboratory gf source for Fe II in
    # the repo. DOI carried on every row RYA-945 ingested, and asserted by the build
    # script rather than assumed -- `grade_line` looks the source up HERE, so a lab
    # table whose source is missing from this dict raises KeyError on the first graded
    # line instead of quietly grading it uncited.
    "DenHartog2019": ("Den Hartog et al. 2019, ApJS 243, 33",
                      "10.3847/1538-4365/ab322e"),
}
NO_TIE_SOURCE = "no-tie->K07 (RYA-161 semi-empirical Kurucz systematic)"


@dataclass
class GradeVerdict:
    gf_grade: str
    gf_grade_source: str
    gf_sigma_dex: float
    gf_reference_tag: str        # what canonical_gf says the source is, evidence only
    gf_ref_loggf: float          # the reference value, for the mismatch audit
    gf_delta_dex: float          # reference minus the value the pool used
    note: str

    @property
    def is_graded(self) -> bool:
        """Primary LABORATORY measurement. Unchanged by RYA-822 on purpose — widening it
        to include compilations would restate 799's headline without saying so."""
        return self.gf_grade == GRADE_LAB

    @property
    def has_cited_sigma(self) -> bool:
        """Carries a citable per-line sigma — a lab measurement OR a NIST accuracy class.
        This is what an error budget can consume instead of the blanket systematic."""
        return self.gf_grade in (GRADE_LAB, GRADE_NIST)


_cache: dict = {}


def lab_lines(species: str = DEFAULT_SPECIES) -> pd.DataFrame:
    """The primary-laboratory gf table for `species`.

    A species with NO registered table raises rather than returning empty. An empty
    frame would send every line of that species to the systematic and look exactly like
    "no lab measurement exists" — the absence-is-a-hypothesis failure (RYA-833). Not
    having a table and there being no measurement are different facts.
    """
    key = f"lab::{species}"
    if key not in _cache:
        path = LAB_TABLES.get(species)
        if path is None:
            raise KeyError(
                f"no primary-lab gf table registered for {species!r}. Registered: "
                f"{sorted(LAB_TABLES)}. Add one to LAB_TABLES (and a CITATIONS entry for "
                f"every source it names) rather than letting the species fall silently to "
                f"the systematic — 'we hold no table' and 'no lab measurement exists' are "
                f"different claims and must not look alike.")
        if not path.exists():
            raise FileNotFoundError(
                f"primary-lab {species} gf table missing at {path} — regenerate with "
                f"`{LAB_REGEN.get(species, '(no regen command registered)')}`. Without it "
                f"every line would fall to the systematic and the run would look like a "
                f"result.")
        _cache[key] = pd.read_csv(path)
    return _cache[key]


def canonical_species(species: str = DEFAULT_SPECIES) -> pd.DataFrame:
    """`canonical_gf.csv` filtered to one species."""
    key = f"cgf::{species}"
    if key not in _cache:
        df = pd.read_csv(CANONICAL_GF, comment="#", low_memory=False)
        _cache[key] = df[df["species"].astype(str) == species].reset_index(drop=True)
    return _cache[key]


def canonical_fe1() -> pd.DataFrame:
    """Back-compat alias — the Fe I slice. Pre-RYA-1002 name, kept because callers and
    docstrings across the repo refer to it."""
    return canonical_species("Fe I")


def _nearest(df: pd.DataFrame, wave: float, ep: float,
             wcol: str, ecol: str) -> pd.Series | None:
    m = df[(np.abs(df[wcol] - wave) <= WAVE_TOL_A) & (np.abs(df[ecol] - ep) <= EP_TOL_EV)]
    if m.empty:
        return None
    return m.iloc[int(np.abs(m[wcol] - wave).values.argmin())]


def grade_line(wavelength_air_A: float, ep_eV: float, log_gf_used: float,
               species: str = DEFAULT_SPECIES) -> GradeVerdict:
    """Grade one line of `species`, or bound it. Never returns a blank bar.

    `species` defaults to Fe I so every pre-RYA-1002 caller is unchanged.
    """
    lab = _nearest(lab_lines(species), wavelength_air_A, ep_eV,
                   "wavelength_air_A", "elo_eV")
    cgf = _nearest(canonical_species(species), wavelength_air_A, ep_eV,
                   "wavelength_air_A", "excitation_potential_eV")
    tag = str(cgf["loggf_reference"]) if cgf is not None else ""
    ref_loggf = float(cgf["log_gf"]) if cgf is not None else np.nan

    if lab is not None:
        d_lab = float(lab["loggf"]) - float(log_gf_used)
        if abs(d_lab) <= LOGGF_MATCH_TOL:
            cite, doi = CITATIONS[str(lab["source"])]
            return GradeVerdict(
                GRADE_LAB, f"{cite} (DOI {doi})", float(lab["e_loggf_dex"]),
                tag, float(lab["loggf"]), d_lab,
                "primary laboratory measurement; the value the pool used is this value")
        # A lab measurement EXISTS and the pool did not use it. That is the mismatch,
        # and it is the most actionable form of it.
        cite, doi = CITATIONS[str(lab["source"])]
        return GradeVerdict(
            GRADE_MISMATCH, NO_TIE_SOURCE, K07_SYSTEMATIC_DEX, tag,
            float(lab["loggf"]), d_lab,
            f"a PRIMARY lab gf exists ({cite}, DOI {doi}, sigma "
            f"{float(lab['e_loggf_dex']):.2f} dex) but differs from the value the pool "
            f"used by {d_lab:+.3f} dex — the lab accuracy does NOT describe this "
            f"measurement, so the systematic stands")

    # RYA-822: a NIST accuracy class on the canonical row IS a citable per-line sigma —
    # but only for the value it describes, so the value the pool used must match it.
    if cgf is not None:
        ngrade = str(cgf.get("nist_grade", "") or "").strip()
        nsig = nist_sigma_dex(ngrade) if ngrade else float("nan")
        if ngrade and np.isfinite(nsig) and np.isfinite(ref_loggf):
            d_nist = ref_loggf - float(log_gf_used)
            if abs(d_nist) <= LOGGF_MATCH_TOL:
                return GradeVerdict(
                    GRADE_NIST, f"NIST ASD accuracy class {ngrade} "
                                f"(<={NIST_ACC_PCT[ngrade]:.1f}% on A_ki)",
                    nsig, tag, ref_loggf, d_nist,
                    "graded against NIST ASD, which for Fe I in this regime is largely a "
                    "COMPILATION — agreement proves no transcription error, not "
                    "independence (RYA-760: FMW *is* NIST and VALD copies it)")
            return GradeVerdict(
                GRADE_MISMATCH, NO_TIE_SOURCE, K07_SYSTEMATIC_DEX, tag, ref_loggf,
                d_nist,
                f"canonical_gf carries a NIST-graded ({ngrade}) gf for this line but it "
                f"differs from the value used by {d_nist:+.3f} dex — the NIST accuracy "
                f"does NOT describe this measurement, so the systematic stands")

    if cgf is not None and tag and tag != K07_TAG:
        d_ref = ref_loggf - float(log_gf_used)
        if abs(d_ref) <= LOGGF_MATCH_TOL:
            return GradeVerdict(
                GRADE_SYSTEMATIC, NO_TIE_SOURCE, K07_SYSTEMATIC_DEX, tag,
                ref_loggf, d_ref,
                f"value confirmed against canonical_gf reference {tag!r}, but that "
                f"reference publishes no per-line uncertainty we can cite and the line "
                f"carries no NIST accuracy class — attributed, not graded. (NIST "
                f"ASD itself is reachable again via astroquery; the lines1.pl CGI is "
                f"broken only in its `unit` path — RYA-822.)")
        return GradeVerdict(
            GRADE_MISMATCH, NO_TIE_SOURCE, K07_SYSTEMATIC_DEX, tag, ref_loggf, d_ref,
            f"canonical_gf carries a {tag!r}-referenced gf for this line but it differs "
            f"from the value the pool used by {d_ref:+.3f} dex")

    return GradeVerdict(
        GRADE_SYSTEMATIC, NO_TIE_SOURCE, K07_SYSTEMATIC_DEX, tag, ref_loggf,
        (ref_loggf - float(log_gf_used)) if np.isfinite(ref_loggf) else np.nan,
        "no primary-lab and no non-K07 reference tie" if tag in ("", K07_TAG)
        else "no tie")


def grade_pool(pool: pd.DataFrame, *, wave_col: str = "wavelength_air_A",
               ep_col: str = "ep_eV", loggf_col: str = "log_gf") -> pd.DataFrame:
    """Attach gf_grade / gf_grade_source (+ evidence) to a measured pool.

    Row count in == row count out, asserted here rather than trusted: a merge that drops
    or duplicates a line would change the reported graded fraction, which is the whole
    output of this ticket.
    """
    n_in = len(pool)
    out = pool.copy()
    cols = {"gf_grade": [], "gf_grade_source": [], "gf_sigma_dex": [],
            "gf_reference_tag": [], "gf_ref_loggf": [], "gf_delta_dex": [], "gf_note": []}
    for r in pool.itertuples():
        v = grade_line(float(getattr(r, wave_col)), float(getattr(r, ep_col)),
                       float(getattr(r, loggf_col)))
        cols["gf_grade"].append(v.gf_grade)
        cols["gf_grade_source"].append(v.gf_grade_source)
        cols["gf_sigma_dex"].append(v.gf_sigma_dex)
        cols["gf_reference_tag"].append(v.gf_reference_tag)
        cols["gf_ref_loggf"].append(v.gf_ref_loggf)
        cols["gf_delta_dex"].append(v.gf_delta_dex)
        cols["gf_note"].append(v.note)
    for k, v in cols.items():
        out[k] = v
    if len(out) != n_in:
        raise ValueError(f"row count changed {n_in} -> {len(out)}; refusing to write")
    blank = out[out["gf_grade"].astype(str).str.strip() == ""]
    if len(blank):
        raise ValueError(f"{len(blank)} row(s) carry a blank gf_grade — never allowed "
                         f"(RYA-799: never leave the bar blank)")
    return out
