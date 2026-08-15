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
LAB_CSV = _REPO_ROOT / "data" / "reference" / "fe_gf_lab" / "fe1_lab_loggf.csv"
CANONICAL_GF = _REPO_ROOT / "data" / "linelists" / "canonical_gf.csv"

#: RYA-161. The semi-empirical Kurucz gf systematic carried when nothing better ties.
K07_SYSTEMATIC_DEX = 0.20
#: The tag `canonical_gf.loggf_reference` uses for that semi-empirical source.
K07_TAG = "K07"

GRADE_LAB = "GF-LAB"
GRADE_SYSTEMATIC = "systematic:K07"
GRADE_MISMATCH = "systematic:K07/SCALE-MISMATCH"

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
        return self.gf_grade == GRADE_LAB


_cache: dict = {}


def lab_lines() -> pd.DataFrame:
    if "lab" not in _cache:
        if not LAB_CSV.exists():
            raise FileNotFoundError(
                f"primary-lab Fe I gf table missing at {LAB_CSV} — regenerate with "
                f"`python3 scripts/rya799_fetch_fe_gf_lab.py`. Without it every line "
                f"would fall to the systematic and the run would look like a result.")
        _cache["lab"] = pd.read_csv(LAB_CSV)
    return _cache["lab"]


def canonical_fe1() -> pd.DataFrame:
    if "cgf" not in _cache:
        df = pd.read_csv(CANONICAL_GF, comment="#", low_memory=False)
        _cache["cgf"] = df[df["species"].astype(str) == "Fe I"].reset_index(drop=True)
    return _cache["cgf"]


def _nearest(df: pd.DataFrame, wave: float, ep: float,
             wcol: str, ecol: str) -> pd.Series | None:
    m = df[(np.abs(df[wcol] - wave) <= WAVE_TOL_A) & (np.abs(df[ecol] - ep) <= EP_TOL_EV)]
    if m.empty:
        return None
    return m.iloc[int(np.abs(m[wcol] - wave).values.argmin())]


def grade_line(wavelength_air_A: float, ep_eV: float, log_gf_used: float) -> GradeVerdict:
    """Grade one Fe I line, or bound it. Never returns a blank bar."""
    lab = _nearest(lab_lines(), wavelength_air_A, ep_eV, "wavelength_air_A", "elo_eV")
    cgf = _nearest(canonical_fe1(), wavelength_air_A, ep_eV,
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

    if cgf is not None and tag and tag != K07_TAG:
        d_ref = ref_loggf - float(log_gf_used)
        if abs(d_ref) <= LOGGF_MATCH_TOL:
            return GradeVerdict(
                GRADE_SYSTEMATIC, NO_TIE_SOURCE, K07_SYSTEMATIC_DEX, tag,
                ref_loggf, d_ref,
                f"value confirmed against canonical_gf reference {tag!r}, but that "
                f"reference publishes no per-line uncertainty we can cite (NIST ASD is "
                f"externally down) — attributed, not graded")
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
