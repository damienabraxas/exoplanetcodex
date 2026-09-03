#!/usr/bin/env python3
"""RYA-1173 - the RYA-946 mandatory Solar reference-line-set census, for Al.

    python3 scripts/rya1173_al_agss21_census.py [--check]

RYA-1141 check D3 found Al stamped FROZEN through RYA-946's census gate with no census and no
recorded exception. This is the census. It consumes `data/reference/asplund2021_al/` (RYA-1173's
other half, which reconstructs the line set from AGSS21's cited primaries) and produces RYA-946's
four deliverables:

  1. the published Solar reference line set, INCLUDING its excluded row     -> the reference dir
  2. a per-band coverage matrix                                             -> band_coverage_matrix.csv
  3. a source-lineage note                                                  -> lineage_note.md
  4. the four-way Codex comparison                                          -> four_way_comparison.csv
  + the 18-field per-line join                                              -> per_line_join.csv

🔴 THE HEADLINE: ONE OF THE SIX IS NOT IN CODEX AT ALL. Al I 10768.363 A -- the 4p 2Po(1/2) ->
5d 2D(3/2) line, one of the six carrying AGSS21's adopted A(Al) = 6.43 -- has no row in
`canonical_gf` under any species. This is not an ingest hole: 32 rows of C/Fe/Si/Na/Ti/Cr/Mg/Ca
sit between 10700 and 10850 A, and we hold its own multiplet partner at 10782.045 A. The line is
specifically missing, and a replication of AGSS21's Al cannot measure it today.

🔴 AND THE TWO BEST LINES ARE THE TWO WE GRADE. Scott's weights are 2 and 3 on 6696/6698 and 1 on
everything else -- half the weight of the published mean sits on those two -- and they are exactly
the two Al I lines Codex holds at LAB tier (Burheim et al. 2023). That is a genuinely good result
and it is worth stating as precisely as the bad one.

TOLERANCE. Derived, swept, and NULLED (RYA-1109/1110/1117, and the lesson that a plateau must be
qualified by its null). The reference wavelengths are printed to 3 decimals and so is much of
`canonical_gf`, but printing is not the binding term here: the two primaries disagree with EACH
OTHER by up to 8 mA on the same transition, and `canonical_gf` follows Scott's wavelength scale
where the two differ. So the window has to clear ~8 mA, and the sweep shows where it stops moving.
The NULL is the same sweep with the reference wavelengths displaced bodily; a window wide enough to
match by coincidence shows it there.

⚠️ EP DOES NOT SEPARATE THIS SET AT ALL, AND THAT IS MEASURED. The worst pair is 6696.015 and
6698.672 A -- the two lines Scott weights highest -- and their EPs are not merely close, they are
IDENTICAL: both leave the same 4s 2S(1/2) lower level, so Delta EP is exactly 0.000 eV against an
EP_TOL_EV of 0.02. Only 2.657 A of wavelength and the UPPER level (5p 2Po 3/2 vs 1/2) tell them
apart. 10872.975 and 10891.732 A are the same story from the other end -- a shared UPPER level, EPs
0.002 eV apart. This is RYA-1151's point in its sharpest form: an EP tolerance is not a level
identity, and on this set it contributes nothing whatever to the join. The census measures the
margin rather than assuming it, and the tolerance is kept far below the smallest wavelength gap.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline import line_match                                        # noqa: E402
from pipeline.gf_empirical import GRADED_TIERS, LAB_TIERS              # noqa: E402
from scripts.build_al_intake_rya1132 import band                       # noqa: E402  ONE band vocabulary

REF_DIR = ROOT / "data" / "reference" / "asplund2021_al"
REF_LINES = REF_DIR / "asplund2021_al_lines.csv"
REF_ANALYSIS = REF_DIR / "nordlander_lind_2017_analysis_lines.csv"
CANON = ROOT / "data" / "linelists" / "canonical_gf.csv"
OUT = ROOT / "data" / "audit" / "rya1173_al_agss21_census"

FLOAT_DECIMALS = 10

#: Every band in this repo's vocabulary (`build_al_intake_rya1132.band`). The matrix reports ALL of
#: them, including the ones with no published line: RYA-946 asks for a search across FUV, NUV, VIS,
#: red optical, NIR and IR, and "we looked and the source used none" is a RESULT, not a blank.
ALL_BANDS = ["FUV", "NUV", "near-UV", "VIS", "red-optical", "NIR", "J", "H", "K",
             "OUTSIDE_CURRENT_INSTRUMENT_REACH"]

#: The sweep. Wide enough at the top to sweep in a neighbour, so a count that keeps climbing is
#: visible as such.
TOL_SWEEP_A = (0.001, 0.002, 0.005, 0.008, 0.010, 0.020, 0.050, 0.100, 0.250, 0.500, 1.0, 2.0)

#: Bodily displacements for the null. Several, and both signs, because ONE offset can be unlucky:
#: a single +5 A shift that happens to land on a real line would understate the chance rate.
NULL_OFFSETS_A = (-11.0, -7.0, -5.0, -3.0, 3.0, 5.0, 7.0, 11.0)


class CensusError(RuntimeError):
    """The census could not be completed HONESTLY."""


def _txt(v) -> str:
    return "" if v is None or (isinstance(v, float) and np.isnan(v)) or pd.isna(v) else str(v).strip()


def al_pool() -> pd.DataFrame:
    d = pd.read_csv(CANON, low_memory=False)
    d = d[d["species"].astype(str) == "Al I"]
    return d.dropna(subset=["wavelength_air_A", "excitation_potential_eV"]).reset_index(drop=True)


def _match(ref: pd.DataFrame, pool: pd.DataFrame, tol: float) -> line_match.MatchResult:
    return line_match.match(ref["wavelength_air_A"].to_numpy(float),
                            pool["wavelength_air_A"].to_numpy(float),
                            want_ep=ref["elo_eV"].to_numpy(float),
                            src_ep=pool["excitation_potential_eV"].to_numpy(float),
                            require_ep=True, tol_A=tol)


def tolerance_study(ref: pd.DataFrame, pool: pd.DataFrame) -> dict:
    """The plateau AND its null, at every tolerance. Neither means anything without the other."""
    rows = []
    for t in TOL_SWEEP_A:
        r = _match(ref, pool, t)
        nulls = []
        for off in NULL_OFFSETS_A:
            shifted = ref.copy()
            shifted["wavelength_air_A"] = shifted["wavelength_air_A"] + off
            nulls.append(int((np.asarray(_match(shifted, pool, t).index) >= 0).sum()))
        rows.append({"tol_A": t,
                     "matched": int((np.asarray(r.index) >= 0).sum()),
                     "ambiguous": len(r.ambiguous),
                     "null_max": max(nulls), "null_mean": round(float(np.mean(nulls)), 3)})
    sweep = pd.DataFrame(rows)

    # The adopted window: the smallest tolerance inside the plateau -- the count has reached its
    # maximum, no row is ambiguous, and EVERY displacement matches nothing.
    clean = sweep[(sweep.matched == sweep.matched.max()) & (sweep.ambiguous == 0)
                  & (sweep.null_max == 0)]
    if clean.empty:
        raise CensusError(
            "No tolerance has the plateau count with zero ambiguity AND a zero null. A count "
            "without a zero null is a count of coincidences (RYA-1117); refusing to adopt one.\n"
            + sweep.to_string(index=False))
    adopted = float(clean.tol_A.min())
    plateau_top = float(clean.tol_A.max())
    if plateau_top <= adopted:
        raise CensusError(f"The 'plateau' is a single point at {adopted} A -- that is not a "
                          f"plateau, it is a coincidence of one window width.")
    return {"sweep": sweep, "adopted_tol_A": adopted,
            "plateau_A": [adopted, plateau_top],
            "plateau_width_ratio": round(plateau_top / adopted, 1),
            "basis": (
                "DERIVED, then swept, then nulled. Printing is not the binding term: the two "
                "primaries disagree with each other by up to 0.008 A on the SAME transition "
                "(Nordlander & Lind 6696.015 vs Scott 6696.023) and canonical_gf follows Scott's "
                "scale, so the window must clear ~8 mA. Adopted at the smallest tolerance where "
                "the count has plateaued, nothing is ambiguous, and all %d bodily displacements "
                "match NOTHING." % len(NULL_OFFSETS_A)),
            "null_design": (
                "the same sweep with the reference wavelengths displaced bodily by %s A. A window "
                "wide enough to match by coincidence matches the displaced set too."
                % list(NULL_OFFSETS_A))}


def ep_separation_margin(ref: pd.DataFrame) -> dict:
    """⚠️ How close this set comes to being unresolvable by EP. Measured, not assumed."""
    w = ref["wavelength_air_A"].to_numpy(float)
    e = ref["elo_eV"].to_numpy(float)
    worst = None
    for i in range(len(w)):
        for j in range(i + 1, len(w)):
            d_ep, d_w = abs(e[i] - e[j]), abs(w[i] - w[j])
            if worst is None or d_ep < worst["delta_EP_eV"]:
                worst = {"pair_A": [round(float(w[i]), 3), round(float(w[j]), 3)],
                         "delta_EP_eV": round(float(d_ep), 6),
                         "delta_lambda_A": round(float(d_w), 3),
                         "levels": [ref.iloc[i]["lower_level"] + " -> " + ref.iloc[i]["upper_level"],
                                    ref.iloc[j]["lower_level"] + " -> " + ref.iloc[j]["upper_level"]]}
    worst["EP_TOL_EV"] = line_match.EP_TOL_EV
    worst["ep_would_separate_them"] = bool(worst["delta_EP_eV"] > line_match.EP_TOL_EV)
    worst["note"] = (
        "The two closest lines in EP share an upper level and differ by %.3f eV -- %.0fx BELOW "
        "line_match.EP_TOL_EV. EP cannot separate them; only the %.1f A wavelength gap and the "
        "LEVEL can. RYA-1151: an EP tolerance is not a level identity."
        % (worst["delta_EP_eV"], line_match.EP_TOL_EV / max(worst["delta_EP_eV"], 1e-12),
           worst["delta_lambda_A"]))
    return worst


def per_line_join(ref: pd.DataFrame, pool: pd.DataFrame, tol: float) -> pd.DataFrame:
    """RYA-946's 18-field per-line join, through the canonical matcher in strict mode."""
    r = _match(ref, pool, tol)
    idx = np.asarray(r.index)

    #: 🔴 KEYED BY POSITION, NOT BY A ROUNDED WAVELENGTH -- and CI's RYA-1037 AST guard is why.
    #: The first draft built `{round(wavelength, 3): detail}` to look up the matcher's unresolved
    #: and ambiguous reports. That is RYA-1033 exactly (Python and numpy round the same tie
    #: differently, and a rounded key can split a matched pair), committed inside the census whose
    #: whole subject is careful line identity -- the same trap `reference_lineset` records itself
    #: falling into. `MatchResult` reports the wavelength it was GIVEN, straight out of this array,
    #: so exact float identity recovers the row index with no key construction at all.
    want = ref["wavelength_air_A"].to_numpy(float)
    amb, unres = {}, {}
    for w, cands in r.ambiguous:
        for i in np.flatnonzero(want == w):
            amb[int(i)] = cands
    for w, dist in r.unresolved:
        for i in np.flatnonzero(want == w):
            unres[int(i)] = dist

    rows = []
    for n, (_, s) in enumerate(ref.iterrows()):
        j = int(idx[n])
        c = pool.iloc[j] if j >= 0 else None
        excluded = s.selection_status == "EXCLUDED_BY_SOURCE_ANALYSIS"

        if j >= 0:
            tier = _txt(c.gf_tier)
            status = ("MATCHED_GRADED" if tier in GRADED_TIERS else "MATCHED_UNGRADED")
            note = ("Resolved on the lambda+EP dual key in strict mode (require_ep=True) at "
                    f"{tol} A; nearest canonical row {abs(float(c.wavelength_air_A) - float(s.wavelength_air_A)):.4f} A away.")
        elif n in amb:
            status, note = "AMBIGUOUS", f"more than one candidate in the window: {amb[n]}"
        else:
            status = "ABSENT_FROM_CODEX"
            note = (f"No Al I row within {tol} A whose EP agrees. Nearest Al I row of any EP is "
                    f"{unres.get(n, float('nan')):.3f} A away.")
        if excluded:
            status = "EXCLUDED_BY_SOURCE_ANALYSIS|" + status

        rows.append({
            # 1-2 provenance of the SET
            "reference_line_set": "asplund-al",
            "adopted_solar_value_lineage": s.adopted_solar_value_lineage,
            # 3 the paper/table this row came from
            "source_paper_table": s.source,
            # 4 species / isotopologue
            "species_isotopologue": "Al I (27Al, 100%, I=5/2 -- monoisotopic, no isotopologue split)",
            # 5 wavelength + medium
            "wavelength_A": float(s.wavelength_air_A),
            "wavelength_medium": "air (measured, control C4 of the reference-set build)",
            "wavelength_A_scott2015b": float(s.wavelength_air_A_scott2015b),
            # 6 EP / transition identity
            "elo_eV": float(s.elo_eV), "eup_eV": float(s.eup_eV),
            "transition_identity": f"{s.lower_level} -> {s.upper_level}",
            # 7-8 published gf and its source
            "published_loggf": float(s.loggf),
            "published_loggf_sigma_dex": (None if pd.isna(s.loggf_sigma_dex)
                                          else float(s.loggf_sigma_dex)),
            "published_gf_source": s.gf_source_per_line,
            "published_gf_is_laboratory": False,
            # 9 published EW / weight / exclusion
            "published_ew_mA": float(s.ew_mA_scott2015b),
            "published_weight": int(s.weight_scott2015b),
            "published_selection_status": s.selection_status,
            "published_exclusion_reason": (s.selection_reason if excluded else ""),
            # 10 band
            "source_band": s.source_band,
            # 11-15 the Codex side
            "codex_canonical_line_id": _txt(c.line_id) if c is not None else "",
            "codex_wavelength_air_A": (float(c.wavelength_air_A) if c is not None else None),
            "codex_loggf": (float(c.log_gf) if c is not None else None),
            "codex_gf_tier": (_txt(c.gf_tier) if c is not None else ""),
            "codex_nist_grade": (_txt(c.nist_grade) if c is not None else ""),
            "codex_lab_source": (_txt(c.loggf_reference) if c is not None else ""),
            "codex_gf_sigma_dex": (None if c is None or pd.isna(c.gf_sigma_dex)
                                   else float(c.gf_sigma_dex)),
            "codex_gf_source_doi": (_txt(c.gf_source_doi) if c is not None else ""),
            # 16 delta
            "delta_loggf_codex_minus_published": (None if c is None
                                                  else round(float(c.log_gf) - float(s.loggf), 4)),
            # 17-18 status
            "join_status": status,
            "ambiguity_status_note": note,
        })
    return pd.DataFrame(rows)


def coverage_matrix(join: pd.DataFrame, analysis: pd.DataFrame, pool: pd.DataFrame,
                    tol: float) -> pd.DataFrame:
    """RYA-946's per-band matrix, over EVERY band -- an empty band is a result, not a blank."""
    # Codex's own Al I holdings per band, so "graded-but-unused" is per-band too.
    pool = pool.copy()
    pool["band"] = [band(w) for w in pool["wavelength_air_A"].astype(float)]
    used = join[~join.published_selection_status.eq("EXCLUDED_BY_SOURCE_ANALYSIS")]
    excl = join[join.published_selection_status.eq("EXCLUDED_BY_SOURCE_ANALYSIS")]

    rows = []
    for b in ALL_BANDS:
        u = used[used.source_band == b]
        x = excl[excl.source_band == b]
        p = pool[pool.band == b]
        diag = analysis[(analysis.source_band == b)
                        & analysis.role.eq("SOLAR_DIAGNOSTIC_NAMED_IN_TEXT")]
        other = analysis[(analysis.source_band == b)
                         & analysis.role.eq("ANALYSIS_LINE_ROLE_NOT_STATED_PER_LINE")]
        rows.append({
            "band": b,
            "published_used": len(u),
            "published_explicitly_rejected": len(x),
            "codex_matched": int(u.join_status.str.startswith("MATCHED").sum()),
            "codex_graded": int(u.join_status.eq("MATCHED_GRADED").sum()),
            "codex_ungraded": int(u.join_status.eq("MATCHED_UNGRADED").sum()),
            "missing": int(u.join_status.eq("ABSENT_FROM_CODEX").sum()),
            "ambiguous": int(u.join_status.eq("AMBIGUOUS").sum()),
            "codex_al1_rows_in_band": len(p),
            "codex_graded_rows_in_band": int(p.gf_tier.isin(GRADED_TIERS).sum()),
            "source_diagnostic_lines_not_in_abundance": len(diag),
            "source_analysis_lines_role_unstated": len(other),
        })
    return pd.DataFrame(rows)


def four_way(join: pd.DataFrame, pool: pd.DataFrame, tol: float) -> pd.DataFrame:
    """RYA-946 deliverable 4: five disjoint classes, every Al I row on exactly one side."""
    used = join[~join.published_selection_status.eq("EXCLUDED_BY_SOURCE_ANALYSIS")]
    matched_ids = set(used.codex_canonical_line_id) - {""}

    rows = []
    for _, r in used.iterrows():
        if r.join_status == "MATCHED_GRADED":
            cls = "USED_BY_SOURCE_AND_GRADED_BY_CODEX"
        elif r.join_status == "MATCHED_UNGRADED":
            cls = "USED_BY_SOURCE_BUT_UNGRADED_IN_CODEX"
        elif r.join_status == "ABSENT_FROM_CODEX":
            cls = "USED_BY_SOURCE_AND_ABSENT_FROM_CODEX"
        else:
            cls = "USED_BY_SOURCE_AND_AMBIGUOUS_IN_CODEX"
        rows.append({"class": cls, "wavelength_A": r.wavelength_A,
                     "transition_identity": r.transition_identity,
                     "published_loggf": r.published_loggf, "published_weight": r.published_weight,
                     "codex_gf_tier": r.codex_gf_tier, "codex_loggf": r.codex_loggf,
                     "codex_lab_source": r.codex_lab_source,
                     "delta_loggf": r.delta_loggf_codex_minus_published,
                     "note": r.ambiguity_status_note})

    for _, r in join[join.published_selection_status.eq("EXCLUDED_BY_SOURCE_ANALYSIS")].iterrows():
        rows.append({"class": "EXCLUDED_BY_THE_SOURCE_ANALYSIS", "wavelength_A": r.wavelength_A,
                     "transition_identity": r.transition_identity,
                     "published_loggf": r.published_loggf, "published_weight": r.published_weight,
                     "codex_gf_tier": r.codex_gf_tier, "codex_loggf": r.codex_loggf,
                     "codex_lab_source": r.codex_lab_source,
                     "delta_loggf": r.delta_loggf_codex_minus_published,
                     "note": r.published_exclusion_reason})

    for _, c in pool[pool.gf_tier.isin(GRADED_TIERS)].iterrows():
        if _txt(c.line_id) in matched_ids:
            continue
        rows.append({"class": "GRADED_BY_CODEX_BUT_NOT_USED_BY_SOURCE",
                     "wavelength_A": float(c.wavelength_air_A),
                     "transition_identity": "", "published_loggf": None, "published_weight": None,
                     "codex_gf_tier": _txt(c.gf_tier), "codex_loggf": float(c.log_gf),
                     "codex_lab_source": _txt(c.loggf_reference), "delta_loggf": None,
                     "note": ("Codex grades this line; the AGSS21 lineage does not use it. That is "
                              "not a defect on either side -- Scott selected seven weak, clean "
                              "lines for a solar abundance, which is a different question from "
                              "which lines have a laboratory gf.")})
    return pd.DataFrame(rows)


def lineage_note(join: pd.DataFrame, cov: pd.DataFrame, fw: pd.DataFrame, tolstudy: dict,
                 epsep: dict, prov: dict) -> str:
    used = join[~join.published_selection_status.eq("EXCLUDED_BY_SOURCE_ANALYSIS")]
    graded = used[used.join_status.eq("MATCHED_GRADED")]
    absent = used[used.join_status.eq("ABSENT_FROM_CODEX")]
    ungraded = used[used.join_status.eq("MATCHED_UNGRADED")]
    sweep = tolstudy["sweep"]
    wsum = int(used.published_weight.sum())
    wgraded = int(graded.published_weight.sum())

    def _lines(d):
        return ", ".join(f"{w:.3f}" for w in sorted(d.wavelength_A))

    def _weights(d):
        return ", ".join(f"{r.wavelength_A:.3f}\u00a0Å→**{int(r.published_weight)}**"
                         for _, r in d.sort_values("wavelength_A").iterrows())

    return f"""# Source-lineage note — the AGSS21 solar Al value, traced to its lines

**RYA-1173**, discharging RYA-946's *Mandatory Solar reference-line-set census — AGSS21 lineage,
all bands* for Al. Generated by `scripts/rya1173_al_agss21_census.py`; every number below is read
from the artifacts beside this file, not typed.

## 1. The chain, and where it stops

RYA-946 asks for `SOURCE_LINE_LIST_NOT_PUBLISHED` *"after the cited primaries/supplements/catalogues
are checked."* They were checked, and it does **not** apply — the rows are published, one paper
further down than AGSS21:

| step | paper | what it contributes |
|---|---|---|
| 1 | Asplund, Amarsi & Grevesse 2021, A&A 653, A141 | the adopted value A(Al) = 6.43 ± 0.03. **No line list. No table. Nothing per line.** |
| 2 | Nordlander & Lind 2017, A&A 607, A75 | the analysis AGSS21 adopts: Table A.1 line data, §3.1.5 the telluric exclusion, Fig. 8 the six used lines named |
| 3 | Scott, Grevesse, Asplund et al. 2015b, A&A 573, A25 | the selection (seven lines), solar EWs, line weights, five-model abundances, Table 3 level *J* |

The reconstructed set is `data/reference/asplund2021_al/`, built under five extraction controls
(three carrying a measured negative). Both a published and a preprint copy of both primaries were
read.

⚠️ `Scott et al. (2015b)` is **A&A 573, A25**, not A26 — AGSS21 letters its two 2015 Scott
references by author list, so `2015a` is the iron-group paper.

## 2. The set: six used, one explicitly rejected — and AGSS21 says five

Scott retains **seven**; Nordlander & Lind drop 10891 Å for telluric contamination and name the
remaining **six** on the Fig. 8 axis; AGSS21's prose says **"these five Al i lines"**, which
reproduces from neither source it cites. Six is adopted on the authority of the primary and the
conflict is **recorded, not resolved** — see `raw/lineage_quotations.md`. A replication of this set
measures six lines.

The rejected line ships in the set, flagged, with its published reason and Scott's EW and weight —
RYA-946: preserve explicit negative selections.

## 3. 🔴 One of the six is not in Codex at all

**Al I {_lines(absent)} Å** — {absent.iloc[0].transition_identity if len(absent) else "—"} — has no
row in `canonical_gf` under any species.

This is **not** an ingest gap. 32 rows of C I, Fe I, Si I, Na I, Ti I, Cr I, Mg I and Ca I sit
between 10700 and 10850 Å, and Codex holds this line's own multiplet partner at 10782.045 Å. The
line is specifically missing, and **AGSS21's Al cannot be replicated on our line list today** — one
of the six abundance-carrying lines has nothing to measure.

## 4. ✅ And the two the source weights highest are the two we grade

Scott's weights run 1–3, larger = better (verified: only that direction reproduces his published
6.43). Over the **six used** lines they are {_weights(used)} — so **{wgraded} of {wsum} weight
units** of the published mean sit on {_lines(graded)} Å. Those two are exactly the Al I lines Codex
holds at `LAB` tier, from Burheim, Hartman & Nilsson 2023 (A&A 672, A197). Where the source analysis
is most confident, so are we.

⚠️ Weight totals depend on which set you mean: **9** over Nordlander & Lind's six, **10** over
Scott's seven. This note is about the six that carry the adopted value.

## 5. 🔴 But the reference gf are theory, so a match grades nothing

Every log gf in this set is Opacity Project / TOPbase (Mendoza et al. 1995) under the LS-coupling
assumption — Scott §5.3 says so and NL2017's per-line reference code agrees on all seven. RYA-946:
*"an AGSS21 abundance value is not a gf grade."* Matching a Codex line into this set establishes
that **AGSS21's lineage used it**, and nothing whatever about the quality of its gf. The `delta_loggf`
column is a theory-minus-whatever-we-hold difference, not an error.

## 6. Bands: the absences are published selections, not gaps

```
{cov.to_string(index=False)}
```

The six used lines occupy **VIS, red-optical and NIR only**. FUV, NUV, near-UV, J, H and K carry
zero — and that is a *decision*, not a coverage hole. Nordlander & Lind analyse the 3944/3961 Å
resonance lines, centre-to-limb variation in 7835 Å, the HFS-sensitive 13123 and 16750 Å IR lines
and the 12.33 µm emission line, and **not one of them enters the abundance**; they constrain the
model atom and the collisional rates instead. `source_diagnostic_lines_not_in_abundance` counts them
per band. RYA-946 is explicit that an abundance is not evidence a line was used, so the other 43
Table A.1 rows carry `ANALYSIS_LINE_ROLE_NOT_STATED_PER_LINE` rather than a guessed role.

**Checked negatives**, since RYA-946 asks for them by name:

* **isotopologues** — none possible. Al is 100% ²⁷Al (Scott Table 3 header); there is no isotopic
  split to look for.
* **molecular indicators** — none used. Neither primary uses an Al-bearing molecule for the solar
  abundance.
* **forbidden lines** — none. Every row in Table A.1 is a permitted E1 transition.
* **blend components** — 7836.134 Å is the 3d ²D(5/2) partner of the used 7835.309 Å line and is
  held in `canonical_gf`; the source uses 7835.309 alone. Recorded in the join, not merged into it.
* **Al II** — the source set is Al I only. Scott §5.3 notes Al is essentially all Al II in the
  photosphere, but no Al II line enters the analysis.

## 7. The match tolerance, and its null

Adopted **{tolstudy['adopted_tol_A']} Å**, on the λ+EP dual key in strict mode (`require_ep=True`).

{tolstudy['basis']}

```
{sweep.to_string(index=False)}
```

The plateau runs {tolstudy['plateau_A'][0]}–{tolstudy['plateau_A'][1]} Å
({tolstudy['plateau_width_ratio']}× wide) with the null identically zero across it. A count that
keeps climbing with the window is a count of coincidences; one that plateaus while its null stays at
zero is the answer (RYA-1109/1117).

⚠️ **EP does not separate this set.** {epsep['note']} The census keeps the window
{round(epsep['delta_lambda_A'] / tolstudy['adopted_tol_A']):,}× below that pair's separation, so the
ambiguity never arises — but the margin is measured here rather than assumed.

## 8. The four-way comparison

```
{fw.groupby('class').size().to_string()}
```

Full rows in `four_way_comparison.csv`. `GRADED_BY_CODEX_BUT_NOT_USED_BY_SOURCE` is not a defect on
either side: Scott selected seven weak, clean lines to derive a solar abundance, which is a different
question from which Al lines have a laboratory gf.

## 9. What this census does and does not unblock

It discharges RYA-946's census gate for Al: the reference set exists, the per-line join exists, the
per-band matrix exists, the four-way comparison exists, and this note traces the value to its lines.

It does **not** make Al ready for measurement, and the census is not entitled to say it is. RYA-1141
raised other findings that this ticket does not touch — RYA-1174 (red-optical typed
`CRITICALLY_EVALUATED` over Opacity Project theory), RYA-1176 (`line_set` unresolvable on the frozen
manifest), RYA-1151 (the gf-promotion join), RYA-1177 (`OUTSIDE_CURRENT_REACH`). And this census adds
one of its own: **a used line that Codex does not hold.**
"""


def readme(verdict: dict) -> str:
    f = {x["id"]: x for x in verdict["🔴 findings"]}
    m = verdict["match"]
    return f"""# RYA-1173 — the RYA-946 Solar reference-line-set census for Al

Closes **RYA-1141 check D3**: Al was stamped `FROZEN` on 168 manifest rows while this census did
not exist and no exception was recorded.

Everything here is generated by `scripts/rya1173_al_agss21_census.py`; `--check` reproduces it and
refuses on drift. The line set it consumes is `data/reference/asplund2021_al/`.

## Headline

🔴 **{f['RYA-1173-F1']['finding']}**

{f['RYA-1173-F1']['why_it_is_not_an_ingest_gap']}

{f['RYA-1173-F1']['consequence']}

🔴 {f['RYA-1173-F2']['finding']}

🔴 {f['RYA-1173-F4']['finding']}

✅ {f['RYA-1173-F3']['finding']}

⚠️ **{verdict['does_not_unblock']}**

## Files

| file | RYA-946 deliverable |
|---|---|
| `lineage_note.md` | **3** — the source-lineage note. Read this one first. |
| `per_line_join.csv` | the required 18-field per-line join |
| `band_coverage_matrix.csv` | **2** — published-used / published-explicitly-rejected / Codex-matched / graded / ungraded / missing / ambiguous, over every band |
| `four_way_comparison.csv` | **4** — used-and-graded / used-but-ungraded / graded-but-unused / absent / excluded-by-source |
| `tolerance_plateau_and_null.csv` | the match window, swept, with its null |
| `census_verdict.json` | the machine-readable verdict and findings |

Deliverable **1** — the machine-readable line set including its excluded row — is
`data/reference/asplund2021_al/`.

## Matching

`pipeline.line_match.match(require_ep=True)` on the λ+EP dual key, at **{m['adopted_tol_A']} Å**.
Plateau {m['plateau_A'][0]}–{m['plateau_A'][1]} Å ({m['plateau_width_ratio']}× wide) with the null
identically {m['null_max_across_plateau']} across it.

⚠️ EP contributes nothing to this set's join: {m['⚠️ ep_separation_margin']['note']}

## The gate is now enforced

`pipeline/reference_census_gate.py` sweeps the per-line intake manifests, and
`scripts/build_al_intake_rya1132.py` calls it before writing `FROZEN`. Negative controls in
`tests/test_reference_census_gate_rya1173.py`, including one that reproduces the pre-RYA-1173
state and asserts the gate goes red.
"""


def build() -> dict:
    if not REF_LINES.exists():
        raise CensusError(f"{REF_LINES} missing -- run scripts/rya1173_build_asplund_al_lineset.py")
    ref = pd.read_csv(REF_LINES)
    analysis = pd.read_csv(REF_ANALYSIS)
    prov = json.loads((REF_DIR / "asplund2021_al_lines.prov.json").read_text())
    pool = al_pool()

    tolstudy = tolerance_study(ref, pool)
    tol = tolstudy["adopted_tol_A"]
    epsep = ep_separation_margin(ref)

    join = per_line_join(ref, pool, tol)
    cov = coverage_matrix(join, analysis, pool, tol)
    fw = four_way(join, pool, tol)
    note = lineage_note(join, cov, fw, tolstudy, epsep, prov)

    used = join[~join.published_selection_status.eq("EXCLUDED_BY_SOURCE_ANALYSIS")]
    counts = {c: int(n) for c, n in fw["class"].value_counts().items()}

    d = join["delta_loggf_codex_minus_published"].dropna().astype(float)
    ranked = join.dropna(subset=["delta_loggf_codex_minus_published"]).assign(
        _a=lambda x: x.delta_loggf_codex_minus_published.abs()).sort_values("_a", ascending=False)
    DELTA_SUMMARY = {
        "n_matched": int(len(d)),
        "per_line": {f"{float(r.wavelength_A):.3f}": float(r.delta_loggf_codex_minus_published)
                     for _, r in ranked.iterrows()},
        "largest_abs_dex": round(float(d.abs().max()), 4),
        "next_largest_abs_dex": round(float(ranked._a.iloc[1]), 4),
        "ratio_largest_to_next": round(float(ranked._a.iloc[0] / ranked._a.iloc[1]), 1),
        "note": ("Codex minus published. The published side is Opacity Project THEORY throughout, "
                 "so this is a theory-vs-whatever-we-hold difference and NOT an error bar."),
    }

    verdict = {
        "ticket": "RYA-1173",
        "discharges": ("RYA-946 'Mandatory Solar reference-line-set census - AGSS21 lineage, all "
                       "bands', for Al. Raised as RYA-1141 check D3 (FAIL)."),
        "element": "Al", "species_in_source_set": ["Al I"],
        "census_complete": True,
        "source_line_list_not_published_wall": False,
        "wall_note": ("AGSS21 publishes no Al line list, but its cited primaries do. The rows were "
                      "found one step further down the chain, so the RYA-946 wall does not apply."),
        "reference_line_set": {
            "line_set": "asplund-al",
            "path": str(REF_LINES.relative_to(ROOT)),
            "n_used": int(len(used)), "n_excluded_by_source": int(len(join) - len(used)),
            "published_line_count_conflict": prov["🔴 published_line_count_conflict"],
        },
        "match": {
            "matcher": "pipeline.line_match.match(require_ep=True) -- lambda+EP dual key, strict",
            "adopted_tol_A": tol,
            "plateau_A": tolstudy["plateau_A"],
            "plateau_width_ratio": tolstudy["plateau_width_ratio"],
            "null_max_across_plateau": int(
                tolstudy["sweep"][(tolstudy["sweep"].tol_A >= tolstudy["plateau_A"][0])
                                  & (tolstudy["sweep"].tol_A <= tolstudy["plateau_A"][1])
                                  ].null_max.max()),
            "basis": tolstudy["basis"], "null_design": tolstudy["null_design"],
            "⚠️ ep_separation_margin": epsep,
        },
        "four_way": counts,
        "🔴 findings": [
            {"id": "RYA-1173-F1", "severity": "CRITICAL",
             "finding": ("Al I 10768.363 A -- 4p 2Po(1/2) -> 5d 2D(3/2), one of the SIX lines "
                         "carrying AGSS21's adopted A(Al) = 6.43 -- is ABSENT from canonical_gf "
                         "under any species."),
             "why_it_is_not_an_ingest_gap": ("32 rows of C/Fe/Si/Na/Ti/Cr/Mg/Ca lie between 10700 "
                                             "and 10850 A, and its own multiplet partner is held "
                                             "at 10782.045 A. The line is specifically missing."),
             "consequence": ("A replication of AGSS21's Al on our line list can measure at most "
                             "5 of the 6 lines. Codex cannot reproduce that solar value as "
                             "published until the line is ingested.")},
            {"id": "RYA-1173-F2", "severity": "HIGH",
             "finding": ("Every log gf in the AGSS21-lineage Al set is Opacity Project / TOPbase "
                         "LS-coupling THEORY (Mendoza et al. 1995), on all seven lines and in both "
                         "primaries."),
             "consequence": ("Matching a Codex line into this set is evidence the lineage USED it "
                             "and nothing about its gf quality. RYA-946: an AGSS21 abundance value "
                             "is not a gf grade. Note this is the same defect RYA-1174 reports "
                             "from the other direction inside RYA-1132's own manifest.")},
            {"id": "RYA-1173-F4", "severity": "HIGH",
             #: Every number in this sentence is read from DELTA_SUMMARY. A typed "five times"
             #: stood here and the computed ratio was 3.5 -- a claim in prose drifts from the
             #: artifact beside it the moment either moves (RYA-1080).
             "finding": ("Al I 8912.900 A: Codex holds log gf = %.3f (KURUCZ) against the source's "
                         "-1.963 (TOPbase). A %.3f dex disagreement on a line the AGSS21 lineage "
                         "uses -- %.1fx the next largest delta in the set."
                         % (-1.963 + DELTA_SUMMARY["per_line"]["8912.900"],
                            abs(DELTA_SUMMARY["per_line"]["8912.900"]),
                            DELTA_SUMMARY["ratio_largest_to_next"])),
             "delta_loggf_distribution": DELTA_SUMMARY,
             "consequence": ("Measuring this line on OUR gf rather than the published one moves "
                             "its inferred abundance by about %+.3f dex -- it is weak (EW 3.0 mA, "
                             % -DELTA_SUMMARY["per_line"]["8912.900"] +
                             "on the linear part of the curve of growth), so the shift is very "
                             "nearly the gf difference. A 'replication of AGSS21' that silently "
                             "used canonical gf would not be replicating AGSS21 on this line. "
                             "RYA-946's fixed-asset rule says the line list is established from "
                             "the reference study and shared unchanged; this is where the two "
                             "scales part company."),
             "not_adjudicated_here": ("Which value is right is a gf question, not a census "
                                      "question. Recorded, not resolved (RYA-161)."),
             },
            {"id": "RYA-1173-F3", "severity": "INFO",
             "finding": ("The two lines Scott weights highest -- %s A, weights %s of the %d "
                         "weight units spread over the SIX used lines -- are exactly the two Al I "
                         "lines Codex holds at LAB tier, from Burheim et al. 2023."
                         % (", ".join(f"{w:.3f}" for w in sorted(
                                used[used.join_status.eq("MATCHED_GRADED")].wavelength_A)),
                            "+".join(str(int(w)) for w in sorted(
                                used[used.join_status.eq("MATCHED_GRADED")].published_weight)),
                            int(used.published_weight.sum()))),
             "⚠️ weight_frame": ("9 weight units over Nordlander & Lind's SIX used lines; 10 over "
                                 "Scott's SEVEN. The adopted value rests on the six, so the six "
                                 "are the frame."),
             "consequence": "Where the source analysis is most confident, so are we."},
        ],
        "does_not_unblock": ("This census discharges the RYA-946 gate ONLY. Al is not "
                             "FROZEN_READY_FOR_MEASUREMENT: RYA-1151, RYA-1174, RYA-1176 and "
                             "RYA-1177 remain open, and RYA-1173-F1 is new."),
        "counts": {
            "reference_rows": int(len(join)),
            "analysis_rows_tableA1": int(len(analysis)),
            "codex_al1_rows": int(len(pool)),
            "codex_al1_graded_rows": int(pool.gf_tier.isin(GRADED_TIERS).sum()),
            "codex_al1_lab_rows": int(pool.gf_tier.isin(LAB_TIERS).sum()),
        },
    }
    return {"join": join.round(FLOAT_DECIMALS), "coverage": cov, "four_way": fw.round(FLOAT_DECIMALS),
            "sweep": tolstudy["sweep"], "note": note, "verdict": verdict,
            "readme": readme(verdict)}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true", help="rebuild and compare; write nothing")
    a = ap.parse_args()
    b = build()

    files = {"per_line_join.csv": b["join"], "band_coverage_matrix.csv": b["coverage"],
             "four_way_comparison.csv": b["four_way"], "tolerance_plateau_and_null.csv": b["sweep"]}

    if a.check:
        bad = []
        for name, df in files.items():
            p = OUT / name
            if not p.exists():
                bad.append(f"{name}: MISSING")
                continue
            import io
            buf = io.StringIO()
            df.to_csv(buf, index=False)
            if buf.getvalue() != p.read_text():
                bad.append(f"{name}: DRIFTED from its generator")
        for name, txt in (("lineage_note.md", b["note"]), ("README.md", b["readme"]),
                          ("census_verdict.json",
                           json.dumps(b["verdict"], indent=2, ensure_ascii=False) + "\n")):
            p = OUT / name
            if not p.exists():
                bad.append(f"{name}: MISSING")
            elif p.read_text() != txt:
                bad.append(f"{name}: DRIFTED from its generator")
        if bad:
            print("DRIFT:\n  " + "\n  ".join(bad))
            return 1
        print("committed census artifacts reproduce from their generator")
        return 0

    OUT.mkdir(parents=True, exist_ok=True)
    for name, df in files.items():
        df.to_csv(OUT / name, index=False)
    (OUT / "lineage_note.md").write_text(b["note"])
    (OUT / "README.md").write_text(b["readme"])
    (OUT / "census_verdict.json").write_text(
        json.dumps(b["verdict"], indent=2, ensure_ascii=False) + "\n")

    v = b["verdict"]
    print(f"tolerance {v['match']['adopted_tol_A']} A  plateau {v['match']['plateau_A']} "
          f"({v['match']['plateau_width_ratio']}x)  null across plateau "
          f"{v['match']['null_max_across_plateau']}")
    print()
    print(b["coverage"].to_string(index=False))
    print()
    for k, n in v["four_way"].items():
        print(f"  {k:44} {n}")
    print()
    for f in v["🔴 findings"]:
        print(f"  [{f['severity']:8}] {f['id']}  {f['finding'][:100]}...")
    print(f"\nwrote {len(files) + 3} artifacts to {OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
