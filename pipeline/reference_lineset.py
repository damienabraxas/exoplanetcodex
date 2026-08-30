"""RYA-1111 - ingest an EXTERNAL reference line list and measure it on our spectra.

🔴 WHY THIS EXISTS. "Replicate Asplund" and "replicate GBS" are not re-slices of our own
pool. RYA-1109 measured it: of AGSS21's 40 Fe I lines only 19 sit in our graded (lab-gf)
tier, and the 11 below our 2.85 eV floor are held at NIST-C+/VALD3/OTHER, not LAB. So a
replication has to MEASURE SOMEONE ELSE'S LINES, ON THEIR gf, IN OUR SPECTRA - and the
resulting product has to say so. `line_set` is that statement, and it is a first-class
product axis, never collapsed into the others (RYA-712).

🔴 ONE VOCABULARY, DEFINED ONCE. Three spellings of this axis were live when this module
was written: `model_registry.LINE_SETS` said `asplund-graded` (my own guess in RYA-1101),
the RYA-1109 artifact's own column says `asplund_agss21`, and RYA-1111's spec says
`asplund`. The spec wins - it is the ticket that owns the axis - and the vocabulary is
imported from `model_registry` rather than restated here, so the registry's guard and this
module's loader cannot drift into disagreeing about what a valid value is.

⚠️ An artifact's NATIVE value is recorded, never silently rewritten. `native_line_set` says
what the file itself carries; the canonical name is what products are tagged with. A reader
who finds `asplund_agss21` in the reference file and `asplund` on a product can see, from
the register, that those are the same axis value and not two different sets.

🔴 THE MATCH TOLERANCE IS A PROPERTY OF THE SOURCE, NOT A CONSTANT. This is the RYA-1109
trap, and it cost a wrong published number: `line_match.MATCH_TOL_A` is 0.005 A, right for
a table printed to 0.01 A. AGSS21 prints lambda in NANOMETRES to two decimals - 0.1 A
resolution - so that default is 20x too tight for it and silently discards real matches
(the Fe I overlap read 2/40 on the default and 19/40 at the derived 0.05 A). Every set
here therefore declares `match_tol_A` WITH ITS BASIS, and `coverage()` ships a plateau
sweep so a count can be seen not to be an artifact of the window.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from pipeline import line_match
from pipeline.model_registry import LINE_SETS  # the ONE vocabulary

ROOT = Path(__file__).resolve().parents[1]

#: EP agreement window. The same 0.02 eV `line_match` uses; the coincident-line clusters
#: it has to separate differ by 0.1-2.3 eV, so this is not the binding constraint.
EP_TOL_EV = line_match.EP_TOL_EV

#: 🔴 A REPRESENTATIONAL EPSILON, NOT A WIDER WINDOW. A wavelength printed to a step S
#: lies within +/-S/2 of truth INCLUSIVE -- the rounding interval is CLOSED at its edge.
#: Comparing `|d| <= tol` in floating point fails exactly there: AGSS21 prints 524.70 nm,
#: the GES list holds 5247.05 A, and the distance evaluates to 0.0500000000001819 against
#: a 0.05 tolerance. That is 1.8e-13 A of float representation error, and on it the gf
#: override declared a line ABSENT that is plainly present.
#:
#: ⚠️ This is NOT the RYA-1109 move of widening a tolerance until a count improves. 1e-9 A
#: is physically meaningless -- a millionth of a mA - and cannot admit a different line;
#: the nearest distinct Fe line anywhere in the GES list is ~1e-3 A away. It exists solely
#: to make a closed interval actually closed.
_CLOSED_EDGE_EPS_A = 1e-9


class ReferenceLineSetError(RuntimeError):
    """A reference set could not be loaded or measured HONESTLY."""


@dataclass(frozen=True)
class ReferenceSet:
    """One external reference line list, and everything needed to read it faithfully."""
    name: str                      # canonical `line_set` value (model_registry.LINE_SETS)
    path: Path
    source: str                    # the citation
    ticket: str
    native_line_set: str           # what the file's own line_set column says
    wl_col: str
    ep_col: str
    gf_col: str
    species_col: str | None        # None => derive from an `ion` column
    gf_source_cols: tuple[str, ...]
    match_tol_A: float
    tol_basis: str
    notes: str = ""
    extra: dict = field(default_factory=dict)


#: 🔴 EXPLICIT ADAPTERS, NOT COLUMN SNIFFING. The two sets genuinely disagree on names
#: (`elo_eV` vs `excitation_potential_eV`, `loggf` vs `log_gf_gbs`, `ion` vs `species`).
#: A heuristic that guessed would bind the wrong column the first day a third set arrives
#: with a plausible-but-different name, and a gf column bound wrongly moves every
#: abundance while looking completely normal.
SETS: dict[str, ReferenceSet] = {
    "asplund": ReferenceSet(
        name="asplund",
        path=ROOT / "data" / "reference" / "asplund2021_fe" / "asplund2021_fe_lines.csv",
        source="Asplund, Amarsi & Grevesse 2021, A&A 653, A141, Table A.2",
        ticket="RYA-1109",
        native_line_set="asplund_agss21",
        wl_col="wavelength_air_A", ep_col="elo_eV", gf_col="loggf",
        species_col=None,                       # has `ion` = I / II
        gf_source_cols=("gf_source_per_line", "gf_source_collective"),
        match_tol_A=0.05,
        tol_basis=("half AGSS21 Table A.2's printed resolution: lambda is printed in "
                   "NANOMETRES to 2 dp = 0.1 A, so half-width 0.05 A (RYA-1109)"),
        notes=("gf_source_per_line is EMPTY BY DESIGN - Table A.2 publishes no per-line "
               "gf source column, so the attribution is collective (RYA-161: flag, do "
               "not guess)."),
    ),
    "gbs": ReferenceSet(
        name="gbs",
        path=ROOT / "data" / "linelists" / "reference_sets" / "gbs_solar_fe_rya1110.csv",
        source="Jofre et al. 2014, A&A 564, A133, Tables 4/5 (the GBS 'golden' Fe set)",
        ticket="RYA-1110",
        native_line_set="gbs",
        wl_col="wavelength_air_A", ep_col="excitation_potential_eV", gf_col="log_gf_gbs",
        species_col="species",
        gf_source_cols=("loggf_ref_gbs_resolved", "gf_provenance_gbs"),
        match_tol_A=0.015,
        tol_basis=("derived in RYA-1110 and re-checked in RYA-1117: Jofre Tables 4/5 "
                   "print lambda in ANGSTROMS to 2 dp (half-width 0.005 A); 0.015 A is "
                   "3x that, and RYA-1117 confirmed the counts were NOT computed on the "
                   "0.005 A default"),
    ),
}


def _normalise_species(d: pd.DataFrame, spec: ReferenceSet) -> pd.Series:
    if spec.species_col:
        return d[spec.species_col].astype(str).str.strip()
    if "ion" not in d.columns:
        raise ReferenceLineSetError(
            f"{spec.name}: no species column and no `ion` column to derive one from")
    return "Fe " + d["ion"].astype(str).str.strip()


def load(name: str) -> pd.DataFrame:
    """A reference set in ONE normalised schema, with its provenance carried per row."""
    if name not in SETS:
        raise ReferenceLineSetError(
            f"unknown line set {name!r}; known: {sorted(SETS)}. Add a ReferenceSet with "
            f"its adapter and its DERIVED match tolerance - do not reuse another set's.")
    spec = SETS[name]
    if not spec.path.exists():
        raise ReferenceLineSetError(
            f"{name}: reference list missing at {spec.path}. It is built by "
            f"{spec.ticket}; there is nothing to substitute.")
    d = pd.read_csv(spec.path)

    native = d["line_set"].astype(str).str.strip().unique().tolist()
    if native != [spec.native_line_set]:
        raise ReferenceLineSetError(
            f"{name}: file's own line_set column is {native}, expected "
            f"[{spec.native_line_set!r}]. The adapter and the artifact disagree about "
            f"what this file IS - refusing rather than relabelling it.")

    out = pd.DataFrame({
        "line_set": spec.name,
        "native_line_set": spec.native_line_set,
        "species": _normalise_species(d, spec),
        "wavelength_air_A": d[spec.wl_col].astype(float),
        "elo_eV": d[spec.ep_col].astype(float),
        "loggf": d[spec.gf_col].astype(float),
        "source": spec.source,
        "source_ticket": spec.ticket,
    })
    # gf provenance: first non-empty of the declared columns, per line. Recorded as-is.
    prov = pd.Series([""] * len(d), index=d.index, dtype=object)
    for c in spec.gf_source_cols:
        if c in d.columns:
            v = d[c].fillna("").astype(str).str.strip()
            prov = prov.where(prov != "", v)
    out["gf_source"] = prov.to_numpy()
    out["gf_source_is_per_line"] = bool(
        spec.gf_source_cols and spec.gf_source_cols[0] in d.columns
        and d[spec.gf_source_cols[0]].fillna("").astype(str).str.strip().ne("").any())

    # 🔴 A LINE WITH NO PUBLISHED gf IS FLAGGED HERE AND REFUSED AT MEASUREMENT -- it is
    # NOT dropped, and it is NOT quietly given a gf from somewhere else.
    #
    # GBS has 21 such lines (20 Fe I, 1 Fe II): their `gf_provenance_gbs` reads "NOT
    # PUBLISHED IN TABLES 4/5". A `heiter2021_log_gf` IS staged for all 21 (Jofre's list
    # is GES-v3, so Heiter+2021 is the plausible source), and adopting it is an OPEN
    # DECISION RYA-1110 flagged for Ryan. Taking it here would silently turn "GBS's own
    # scale" into "GBS where published, GES elsewhere" and report one number for the
    # mixture -- the confound a replication exists to remove (RYA-161/429).
    # 🔴 ADOPTED gf ARE JOINED HERE, FROM A SIDECAR, AND NEVER FROM `log_gf_gbs`.
    # Ryan ratified adopting 12 of GBS's 21 unpublished lines (RYA-1110, 2026-08-30) --
    # the ones whose reference code reproduces Jofre's PUBLISHED gf exactly on >= 5
    # held-out lines. The value is Heiter+2021's, so it lives in its own file with its own
    # provenance columns: writing it into `log_gf_gbs` would make that column assert
    # something false about a publication and make the 12 indistinguishable from the 138.
    #
    # ⚠️ The remaining 9 (1 THIN + 8 RISKY) stay REFUSED. `gf_adopted` records which side
    # of that line each row is on, so a product can never quietly include one.
    out["gf_adopted"] = False
    out["gf_adopted_source"] = ""
    adopt_path = spec.path.parent / "gbs_gf_adoption_rya1110.csv"
    if spec.name == "gbs" and adopt_path.exists():
        ad = pd.read_csv(adopt_path)
        # 🔴 THE CANONICAL MATCHER, NOT A ROUNDED KEY. An earlier draft of this join keyed
        # on `round(wavelength, 2)` -- and CI's RYA-1037 AST guard caught it, in a module
        # whose entire subject is careful line identity. A rounded wavelength is not an
        # identity (RYA-1033): it splits a matched pair, and Python and numpy round the
        # same tie differently. Matched here on the lambda+EP dual key at THIS set's
        # derived tolerance, exactly as every other join in this module.
        idx = line_match.match(
            ad["wavelength_air_A"].to_numpy(float),
            out["wavelength_air_A"].to_numpy(float),
            want_ep=ad["elo_eV"].to_numpy(float),
            src_ep=out["elo_eV"].to_numpy(float),
            require_ep=True, tol_A=spec.match_tol_A).index
        for k, j in enumerate(np.asarray(idx)):
            if j < 0:
                raise ReferenceLineSetError(
                    f"adopted gf row {k} ({ad.iloc[k]['species']} "
                    f"{ad.iloc[k]['wavelength_air_A']} A) does not resolve into the {name} "
                    f"set on the lambda+EP key -- refusing to place a gf by position")
            if not pd.isna(out.at[j, "loggf"]):
                raise ReferenceLineSetError(
                    f"adopted gf would OVERRIDE a published value at "
                    f"{out.at[j, 'wavelength_air_A']} A -- refusing (RYA-161)")
            hit = ad.iloc[k]
            out.at[j, "loggf"] = float(hit["log_gf_adopted"])
            out.at[j, "gf_adopted"] = True
            out.at[j, "gf_adopted_source"] = str(hit["gf_adopted_source"])
            out.at[j, "gf_source"] = (
                f"ADOPTED (not published by this set): {hit['gf_adopted_source']} | "
                f"{hit['gf_adopted_basis']} | {hit['gf_adopted_ticket']}")

    out["gf_missing"] = out["loggf"].isna()
    if out["gf_missing"].any():
        out.loc[out["gf_missing"], "gf_source"] = (
            out.loc[out["gf_missing"], "gf_source"].astype(str)
            + "  [NO PUBLISHED gf -- not measurable on this set's own scale]")
    return out


def measurable(ref: pd.DataFrame) -> pd.DataFrame:
    """The rows a replication may actually measure: those with a published gf."""
    return ref[~ref["gf_missing"]].reset_index(drop=True)


def gf_gap(name: str) -> dict:
    """The lines this set publishes without a gf, named. A coverage fact, not a footnote."""
    ref = load(name)
    miss = ref[ref["gf_missing"]]
    return {
        "line_set": name,
        "n_ref": int(len(ref)),
        "n_without_published_gf": int(len(miss)),
        "n_measurable": int(len(ref) - len(miss)),
        "by_species": {k: int(v) for k, v in miss["species"].value_counts().items()},
        "lines": [{"wavelength_air_A": round(float(r.wavelength_air_A), 3),
                   "elo_eV": round(float(r.elo_eV), 4), "species": r.species}
                  for r in miss.itertuples()],
        "disposition": ("REFUSED at measurement, not dropped and not substituted. "
                        "Adopting another table's gf for these is an OPEN DECISION "
                        "(RYA-1110) and is Ryan's, not this module's."),
    }


def match_to(ref: pd.DataFrame, pool: pd.DataFrame, spec: ReferenceSet, *,
             pool_wl: str = "wavelength_air_A",
             pool_ep: str = "excitation_potential_eV") -> np.ndarray:
    """Index into `pool` per reference line, or -1. Strict lambda+EP (RYA-1037)."""
    if ref.empty or pool.empty:
        return np.full(len(ref), -1, int)
    r = line_match.match(ref["wavelength_air_A"].to_numpy(float),
                         pool[pool_wl].to_numpy(float),
                         want_ep=ref["elo_eV"].to_numpy(float),
                         src_ep=pool[pool_ep].to_numpy(float),
                         require_ep=True, tol_A=spec.match_tol_A)
    return np.asarray(r.index)


def plateau(ref: pd.DataFrame, pool: pd.DataFrame, spec: ReferenceSet,
            tols=(0.005, 0.01, 0.015, 0.02, 0.05, 0.10, 0.25, 0.50)) -> dict:
    """Match count vs window. A real overlap plateaus; coincidences keep climbing.

    Shipped WITH every coverage report, because the count is only meaningful once the
    reader can see it is not a function of the tolerance (RYA-1109).
    """
    out = {}
    for t in tols:
        s = ReferenceSet(**{**spec.__dict__, "match_tol_A": t})
        out[str(t)] = int((match_to(ref, pool, s) >= 0).sum())
    return out


def coverage(name: str, pool: pd.DataFrame, *, pool_label: str,
             pool_wl: str = "wavelength_air_A",
             pool_ep: str = "excitation_potential_eV") -> dict:
    """What fraction of a reference set a pool can serve, and WHICH lines it cannot.

    🔴 LOUD, NEVER SILENT. Every unmatched line is named with its wavelength and EP. A
    reference line our data cannot serve is a RESULT about coverage (RYA-429/711); it is
    never dropped quietly and never replaced by a neighbour.
    """
    spec = SETS[name]
    ref = load(name)
    idx = match_to(ref, pool, spec, pool_wl=pool_wl, pool_ep=pool_ep)
    hit = idx >= 0
    per_species = {}
    for sp in sorted(ref["species"].unique()):
        m = (ref["species"] == sp).to_numpy()
        per_species[sp] = {"n_ref": int(m.sum()), "matched": int((hit & m).sum()),
                           "unmatched": int(((~hit) & m).sum())}
    return {
        "line_set": name,
        "source": spec.source,
        "pool": pool_label,
        "match_tol_A": spec.match_tol_A,
        "tol_basis": spec.tol_basis,
        "ep_tol_eV": EP_TOL_EV,
        "n_ref": int(len(ref)),
        "matched": int(hit.sum()),
        "unmatched": int((~hit).sum()),
        "per_species": per_species,
        "unmatched_lines": [
            {"wavelength_air_A": round(float(r.wavelength_air_A), 3),
             "elo_eV": round(float(r.elo_eV), 4), "species": r.species}
            for r in ref[~hit].itertuples()],
        "plateau": plateau(ref, pool, spec),
    }


def apply_gf_override(linelist, targets: pd.DataFrame, spec: ReferenceSet, *,
                      on_missing: str = "raise") -> dict:
    """Put the REFERENCE's own log gf on the target lines of an in-memory synthesis list.

    Generalised from RYA-1106's Asplund-specific version; the reasoning is unchanged.

    🔴 IN MEMORY, NOT VIA A WRITTEN LIST. The fitter is handed the same array object, so
    there is no second file to drift and no write/read round-trip to normalise a value
    behind our back (RYA-1084's lesson about round-tripping data files applies to line
    lists too). `canonical_gf.csv` is NEVER written.

    ⚠️ ONLY THE TARGET LINES MOVE. Blends inside the window keep canonical gf exactly as
    production set them: the blend context is OUR modelling choice and is not part of what
    the reference contributed. Substituting their scale into lines they never measured
    would invent a third thing that is neither their analysis nor ours.

    🔴 REFUSES A PARTIAL OVERRIDE. Measuring some lines on their scale and the rest on
    ours, then reporting one number for the mixture, is exactly the confound a replication
    exists to remove (RYA-429).
    """
    names = linelist.dtype.names
    w_A = np.asarray(linelist["wave_A"] if "wave_A" in names
                     else linelist["wave_nm"] * 10.0, dtype=float)
    ep = np.asarray(linelist["lower_state_eV"], dtype=float)
    el = np.asarray([str(x).strip() for x in linelist["element"]])

    applied, missing, ambiguous = [], [], []
    for r in targets.itertuples():
        want_el = str(r.species).replace(" I", " 1").replace(" II", " 2")
        sel = ((el == want_el)
               & (np.abs(w_A - r.wavelength_air_A)
                  <= spec.match_tol_A + _CLOSED_EDGE_EPS_A)
               & (np.abs(ep - r.elo_eV) <= EP_TOL_EV))
        hit = np.flatnonzero(sel)
        if hit.size == 0:
            missing.append({"wavelength_air_A": round(float(r.wavelength_air_A), 4),
                            "elo_eV": round(float(r.elo_eV), 4),
                            "species": str(r.species)})
            continue
        if hit.size > 1:
            # A hyperfine cluster is legitimately several rows for ONE physical line, and
            # the gf then belongs on every row of it. But that is only true while the rows
            # agree on being the same transition -- never resolved by argmin.
            spread = float(np.ptp(w_A[hit]))
            if spread > spec.match_tol_A + _CLOSED_EDGE_EPS_A:
                ambiguous.append({"wavelength_air_A": round(float(r.wavelength_air_A), 4),
                                  "n_rows": int(hit.size),
                                  "wavelength_spread_A": round(spread, 5)})
                continue
        before = np.asarray(linelist["loggf"], dtype=float)[hit]
        linelist["loggf"][hit] = float(r.loggf)
        applied.append({"wavelength_air_A": round(float(r.wavelength_air_A), 4),
                        "elo_eV": round(float(r.elo_eV), 4),
                        "species": str(r.species), "n_rows": int(hit.size),
                        "loggf_canonical": round(float(np.mean(before)), 4),
                        "loggf_reference": round(float(r.loggf), 4),
                        "delta_dex": round(float(r.loggf - np.mean(before)), 4),
                        "gf_source": str(getattr(r, "gf_source", ""))[:200]})
    # 🔴 TWO DIFFERENT FAILURES, AND COLLAPSING THEM WAS WRONG.
    #   * AMBIGUOUS -- the line IS in the list but we cannot say which row is it. Placing
    #     a gf there could move the wrong transition, so this always refuses.
    #   * MISSING   -- the line is NOT in the synthesis list at all (AGSS21's 9786.6 A sits
    #     past the GES list's 9200 A red edge). Nothing can be measured there by any
    #     scale, so it is not a "partial override" at all: it is COVERAGE, and the ticket
    #     asks for coverage to be REPORTED, never silently dropped (RYA-429/711).
    # `on_missing="report"` is what a replication run passes, and it must then EXCLUDE and
    # NAME those lines. The default stays "raise" so no existing caller changes behaviour.
    if on_missing not in ("raise", "report"):
        raise ReferenceLineSetError(f"on_missing must be 'raise' or 'report', not {on_missing!r}")
    if ambiguous or (missing and on_missing == "raise"):
        raise ReferenceLineSetError(
            f"the synthesis line list cannot carry {spec.name}'s gf for the whole pool: "
            f"{len(missing)} absent on the lambda+EP dual key (+/-{spec.match_tol_A} A / "
            f"+/-{EP_TOL_EV} eV) {missing[:4]}, {len(ambiguous)} ambiguous "
            f"{ambiguous[:3]}. A PARTIAL override would measure some lines on their scale "
            f"and the rest on ours and report one number for the mixture - refusing "
            f"(RYA-429).")
    deltas = np.array([a["delta_dex"] for a in applied], dtype=float)
    if deltas.size == 0:
        raise ReferenceLineSetError(
            f"{spec.name}: not one target line is in the synthesis list -- nothing to "
            f"measure. A leg that fits ZERO lines must fail, not emit an empty product.")
    return {
        "line_set": spec.name, "n_targets": int(len(targets)), "n_applied": len(applied),
        "not_in_synthesis_list": missing,
        "per_line": applied,
        "delta_vs_canonical_dex": {
            "mean": round(float(deltas.mean()), 4),
            "median": round(float(np.median(deltas)), 4),
            "min": round(float(deltas.min()), 4),
            "max": round(float(deltas.max()), 4),
            "n_exact": int((np.abs(deltas) < 1e-9).sum()),
        },
        "note": (f"{spec.name}'s published log gf on the target lines; canonical gf "
                 f"everywhere else in the window. RYA-353 single-sourcing is DECLARED "
                 f"OFF for the targets and canonical_gf.csv is not written."),
    }


# ── wiring the axis to the feed ───────────────────────────────────────────────

#: Our own products already say which pool they were measured on, in `tier`. The axis is
#: therefore DERIVED for them rather than stored twice.
_TIER_TO_LINE_SET = {"GRADED": "our-graded", "DEEPGRADED": "our-deep-graded"}


def line_set_for_product(product: dict) -> str:
    """The `line_set` axis value for one feed product.

    🔴 STORED FOR A REPLICATION, DERIVED FOR OUR OWN -- and the asymmetry is deliberate.

    A replication product is measured on someone else's list, so nothing already in the
    record implies which; it carries an explicit `line_set` and this returns it.

    Our own products already state their pool in `tier` (GRADED / DEEPGRADED), and those
    map one-to-one onto `our-graded` / `our-deep-graded`. Storing it again would create a
    second source of truth for a value the record already carries -- and the two would be
    free to disagree, which is the defect the model registry exists to end.

    ⚠️ IT REFUSES AN UNKNOWN TIER RATHER THAN DEFAULTING. `treatment_axes` warns that a
    derivation which is correct today goes silently wrong the first day the correlation
    breaks; the guard against that is to fail loudly the day a third tier appears, not to
    pick the nearest value. A caller that needs a new tier supported must SAY so here.

    ⚠️ It does NOT write anything. Back-stamping an explicit `line_set` onto the 66 live
    products is a feed edit governed by RYA-1080's published-value guard and is a separate
    ticket; RYA-1111 adds the axis and the replication products, and changes no existing
    product value (RYA-161).
    """
    stored = str(product.get("line_set") or "").strip()
    if stored:
        if stored not in LINE_SETS:
            raise ReferenceLineSetError(
                f"product carries line_set {stored!r}, which is not in the vocabulary "
                f"{LINE_SETS}. Add it to model_registry.LINE_SETS deliberately, or fix "
                f"the tag -- do not widen the reader to accept it.")
        return stored
    tier = str(product.get("tier") or "").strip().upper()
    if tier not in _TIER_TO_LINE_SET:
        raise ReferenceLineSetError(
            f"cannot derive line_set: product has no explicit `line_set` and its tier is "
            f"{tier!r}, which is not one of {sorted(_TIER_TO_LINE_SET)}. A product "
            f"measured on an unrecognised pool must SAY which -- refusing to guess.")
    return _TIER_TO_LINE_SET[tier]


def tag_product(product: dict, name: str, gf_source_note: str = "") -> dict:
    """Stamp a replication product with its line set and whose gf it used.

    Returns a NEW dict; the caller decides where it lands. Nothing here writes the feed.
    """
    if name not in SETS:
        raise ReferenceLineSetError(f"unknown line set {name!r}; known: {sorted(SETS)}")
    spec = SETS[name]
    out = dict(product)
    out["line_set"] = spec.name
    out["line_set_source"] = spec.source
    out["line_set_ticket"] = spec.ticket
    out["gf_source"] = gf_source_note or spec.notes or spec.source
    return out
