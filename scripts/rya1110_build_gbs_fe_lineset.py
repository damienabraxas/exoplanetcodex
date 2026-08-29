#!/usr/bin/env python3
"""RYA-1110 — the Gaia FGK Benchmark Stars (Jofré+ 2014) SOLAR Fe reference line set.

    python3 scripts/rya1110_build_gbs_fe_lineset.py

Builds, from the RYA-1110 holding under `data/reference/jofre2014_gbs/`:

    data/linelists/reference_sets/gbs_solar_fe_rya1110.csv           the line set
    data/linelists/reference_sets/gbs_solar_fe_coverage_rya1110.csv  per-VIS-holding coverage

This is the reference line set for the GBS/GAIA replication — the sister of RYA-1109's
Asplund set. It is a TRANSCRIPTION plus a JOIN, never a re-derivation: no gf is
substituted, no value is edited, and everything the published record does not say is left
empty and reported as empty (RYA-161 firewall).

WHAT "THE GBS SOLAR LINE SET" MEANS HERE
----------------------------------------
Jofré+2014 Table 3 states the Sun's selected-line counts outright: N(Fe I) = 150,
N(Fe II) = 9. VizieR `table6.dat` carries a per-star usage flag and `ew.dat` a per-star
measurement block; both give 150 + 9 = 159 for the Sun, and the same identity holds for
all 34 benchmark stars (checked, not assumed — `_check_published_counts`). So the
published solar set is those 159 lines, and this file carries all 159.

🔴 THE PAPER'S OWN -4.8 CUT DOES NOT SELECT THE SET IT PUBLISHES
----------------------------------------------------------------
Sect. 6.1 states the first selection step exactly: "we selected those lines with
log(EW/λ) ≤ -4.8", to stay on the linear part of the curve of growth. Applied to the
Sun's own published equivalent widths, **14 of those 159 lines exceed -4.8 on EVERY ONE of
the six methods that measured them** — not marginally: Fe I 6393.6 sits at -4.66. It is
not a solar quirk; over all 4252 published star-line rows, 713 (16.8%) violate the stated
cut on every method, and the rate tracks stellar type (0.0% for HD 140283 and HD 84937,
32.0% for μ Leo) exactly as it would if the cut had not been applied to what is tabulated.

The ticket asks for THEIR selection, and their two statements disagree. Neither is
overridden here. The file carries both:

  * `gbs_selected_sun` = 1 on all 159 rows — the published selection, as published;
  * `rew_class` — this ticket's application of the published -4.8 rule to the published
    solar EWs.

WHY `rew_class` HAS THREE VALUES AND NOT TWO
--------------------------------------------
Each line has up to six independent EW measurements (EPINARBO, UCM, Porto, Bologna, ULB,
LUMBA) and the paper does not say which one the cut was applied to. For most lines it does
not matter — every method lands on the same side. For three it does. Widening the passing
set to absorb those three would launder an ambiguity into a decision (RYA-1072), so there
are two closed disjoint sets and an explicitly ambiguous remainder:

    pass       every method that measured the line gives log(EW/λ) ≤ -4.8   (142)
    excluded   every method gives log(EW/λ) > -4.8                          ( 14)
    ambiguous  the methods straddle -4.8                                    (  3)

The GBS replication set is `rew_class == "pass"`. The other 17 rows stay in the file
because a dropped row cannot be audited (RYA-931: quarantine, never cull).

THE JOIN TO OUR POOL
--------------------
λ+EP dual key through `pipeline.line_match` (RYA-1037), `require_ep=True`. The tolerance
is DERIVED, not chosen, and its null is MEASURED (RYA-1070's lesson — a tolerance without
a measured null is not a tolerance):

  * Jofré prints λ to 2 dp, so agreement can be no better than ±5 mÅ from rounding alone.
  * The measured worst genuine disagreement over the 159 lines is 5.0 mÅ (Fe I) and
    14.0 mÅ (Fe II, 2 of 9 rows).
  * `_MATCH_TOL_A` = 0.015 Å resolves all 159 with ZERO ambiguous candidates. At 0.030 Å
    a genuine fork appears, so the window is not near the ambiguity edge.
  * NULL CONTROL, asserted at build time: re-run the same match with the GBS wavelengths
    DISPLACED by ±0.2/0.3/0.5 Å and the resolved count must be 0. It is — the EP half of
    the key is what makes it 0, which is the whole argument for the dual key.

THE DECISION FLAG (spec item 4 — NOT resolved here)
---------------------------------------------------
Both gf columns ship, so either can be selected downstream:

    log_gf_gbs    Jofré Tables 4/5, as published   (138 of 159 rows have one)
    log_gf_ours   canonical_gf's adopted value     (159 of 159)
    gf_synth_ges  our GES seed, for the v3-vs-v5 question

`docs/design/rya1110_gbs_reference_lineset.md` records the measured consequences of the
choice. It is Ryan's, and nothing here presumes it.
"""
from __future__ import annotations

import csv
import math
import re
import statistics
import sys
import unicodedata
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from pipeline import line_match                                    # noqa: E402
from pipeline import band_policy                                   # noqa: E402
from pipeline.telluric_policy import exclusion as telluric_exclusion  # noqa: E402

HOLDING = ROOT / "data" / "reference" / "jofre2014_gbs"
VIZIER = HOLDING / "vizier"
GES = ROOT / "data" / "reference" / "heiter2021_ges"
OUT = ROOT / "data" / "linelists" / "reference_sets"
LINESET = OUT / "gbs_solar_fe_rya1110.csv"
COVERAGE = OUT / "gbs_solar_fe_coverage_rya1110.csv"

LINE_SET_TAG = "gbs"           # pipeline.model_registry.LINE_SETS
BAND = "VIS"
STAR = "Sun"

#: The published criterion, verbatim from Sect. 6.1: log(EW/λ) ≤ -4.8.
REW_CUT = -4.8

#: See the module docstring. DERIVED from the 2-dp print precision and the MEASURED
#: table-to-table spread; its chance-match rate is asserted to be 0 at build time.
_MATCH_TOL_A = 0.015

#: Displacements for the null control. Large enough to break every real identification,
#: small enough that the line DENSITY the match sees is unchanged.
_NULL_SHIFTS_A = (-0.5, -0.3, -0.2, 0.2, 0.3, 0.5)

#: `table6.dat` prints λ to 1 dp where `ew.dat` prints 2. Half a step is the most that
#: rounding can move a value, and a full step is the closest two lines may sit before the
#: ordinal pairing in `_bind_table6` stops being provably correct.
_PRINT_STEP_A = 0.1
_PRINT_HALF_STEP_A = 0.05
#: 4873.8 - 4873.75 is 0.050000000000454747 in binary floating point, so the half-step
#: bound has to be compared with representation slack or the exactly-on-the-bound pairs
#: (there are several, and 0.05 IS the measured maximum) read as violations.
_FP_SLACK_A = 1.0e-9

#: Jofré Table 3's own solar counts — the published number this build must reproduce.
PUBLISHED_SUN_COUNTS = {"Fe I": 150, "Fe II": 9}

#: Heiter+2021 keys species as (Element, Ion); we key it as a name. One mapping, here.
_GES_ION = {"Fe I": 1, "Fe II": 2}

#: A gf value and a gf SOURCE agree only if the values match to the precision the source
#: table prints (3 dp on both sides). Anything looser would let a revised value keep the
#: old value's pedigree, which is the `gf_grades` SCALE-MISMATCH defect.
_GF_SAME_DEX = 0.0005

_SPECIES = {"FeI": "Fe I", "FeII": "Fe II"}
_METHODS = ("EPINARBO", "UCM", "Porto", "Bologna", "ULB", "LUMBA")
_NULL = -999.00


class BuildError(RuntimeError):
    """A published fact this build depends on did not hold."""


# ── the published holding ────────────────────────────────────────────────────
def _read_measurement_table(path: Path) -> list[dict]:
    """ew.dat / abund.dat — identical byte layout, six per-method columns."""
    rows = []
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        vals = [float(line[37 + 10 * i:44 + 10 * i]) for i in range(6)]
        rows.append(dict(star=line[0:9].strip(),
                         species=_SPECIES[line[10:14].strip()],
                         wavelength_air_A=float(line[17:24]),
                         excitation_potential_eV=float(line[29:34]),
                         values=[v for v in vals if v != _NULL],
                         methods=[m for m, v in zip(_METHODS, vals) if v != _NULL]))
    return rows


def _read_table6(path: Path) -> dict:
    """table6 as a flat list: species, λ (1 dp AS PUBLISHED), golden flag, Sun usage flag.

    λ is kept exactly as printed and is never re-rounded to meet ew.dat's 2 dp. The two
    tables genuinely disagree on 7 lines through rounding alone (ew.dat 4985.55 is printed
    4985.5 here, not 4985.6), and a re-rounded key would silently drop those 7 — RYA-1033,
    a rounded number is not an identity. `_bind_table6` pairs the two ordinally instead.
    """
    out = []
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        out.append(dict(species=_SPECIES[line[0:4].strip()],
                        wavelength_air_A=float(line[5:11]),
                        golden=line[12],
                        sun=line[78]))
    return out


def _bind_table6(ew_rows: list[dict], t6: list[dict]) -> dict:
    """(species, ew λ) -> its table6 row.

    🔴 THIS IS NOT A NEAREST-WAVELENGTH SEARCH, AND MUST NOT BECOME ONE. `table6.dat` and
    `ew.dat` are the SAME 242 lines of the same paper printed at two precisions (1 dp and
    2 dp), and `table6` carries no excitation potential at all — so there is no second key
    to disambiguate with, and a tolerance match here would be the exact λ-only join
    RYA-1037 exists to end. Instead the pairing is ORDINAL, and the one condition under
    which an ordinal pairing could be wrong is ASSERTED rather than handled:

      * the two tables must list the same NUMBER of lines per species (242 = 224 + 18);
      * no two lines of a species may sit within one print step of each other, because
        that is the only way 1-dp rounding could reorder them. Measured: the closest pair
        is 0.19 Å apart (Fe I) and 0.16 Å (Fe II), against a 0.1 Å step — so no rounding
        can cross a neighbour and sorted order is preserved exactly;
      * and each pair must then differ by no more than HALF a print step, which is all the
        rounding is allowed to move it. Measured maximum: exactly 0.05 Å, the bound.

    Every one of those is checked below, so a changed holding breaks loudly instead of
    silently pairing the wrong rows.
    """
    bound: dict = {}
    for species in sorted({r["species"] for r in ew_rows}):
        want = sorted({r["wavelength_air_A"] for r in ew_rows
                       if r["species"] == species})
        have = sorted(t["wavelength_air_A"] for t in t6 if t["species"] == species)
        if len(want) != len(have):
            raise BuildError(
                f"{species}: ew.dat lists {len(want)} lines and table6 lists {len(have)}. "
                f"They are the same table at two precisions; unequal counts mean the "
                f"holding changed and the pairing is not the identity.")
        gaps = [b - a for a, b in zip(want, want[1:])]
        if gaps and min(gaps) <= _PRINT_STEP_A:
            i = gaps.index(min(gaps))
            raise BuildError(
                f"{species}: {want[i]} and {want[i + 1]} are {min(gaps):.3f} Å apart, "
                f"within the {_PRINT_STEP_A} Å print step. 1-dp rounding could reorder "
                f"them, so the ordinal pairing is no longer provably the identity.")
        if len(set(have)) != len(have):
            raise BuildError(
                f"{species}: table6 prints the same wavelength twice. The ordinal pairing "
                f"assumes one row per line; two rows on one λ would collapse silently.")
        by_lam = {t["wavelength_air_A"]: t for t in t6 if t["species"] == species}
        for a, b in zip(want, have):
            dev = b - a
            if not (-_PRINT_HALF_STEP_A - _FP_SLACK_A <= dev
                    <= _PRINT_HALF_STEP_A + _FP_SLACK_A):
                raise BuildError(
                    f"{species}: ew.dat {a} pairs ordinally with table6 {b}, which is "
                    f"{dev:+.3f} Å away — more than rounding to 1 dp can move it.")
            bound[(species, a)] = by_lam[b]
    return bound


def _read_paper_tables() -> dict:
    """(species, λ) -> the paper's atomic data. λ is 2 dp on BOTH sides, so exact."""
    out = {}
    for species, name in (("Fe I", "paper_table4_fe1_golden.tsv"),
                          ("Fe II", "paper_table5_fe2_golden.tsv")):
        with (HOLDING / name).open() as fh:
            for row in csv.DictReader((l for l in fh if not l.startswith("#")),
                                      delimiter="\t"):
                key = (species, float(row["wavelength_air_A"]))
                if key in out:
                    raise BuildError(f"duplicate paper row for {key}")
                out[key] = dict(elow_eV=float(row["elow_eV"]),
                                log_gf=float(row["log_gf"]),
                                vdw_abo=float(row["vdw_abo"]),
                                ref_code=int(row["loggf_ref_code"]))
    return out


def _read_tsv(path: Path) -> list[dict]:
    """A staged decoder TSV, comment lines dropped."""
    with path.open() as fh:
        return list(csv.DictReader((l for l in fh if not l.startswith("#")),
                                   delimiter="\t"))


def _read_ges_lines() -> pd.DataFrame:
    df = pd.DataFrame(_read_tsv(GES / "geslines_Fe_4700_6900.tsv"))
    if df.empty:
        raise BuildError(f"{GES / 'geslines_Fe_4700_6900.tsv'} is empty — run "
                         f"scripts/rya1110_stage_heiter2021.py")
    for c in ("lambda", "loggf", "e_loggf", "Elow"):
        df[c] = df[c].astype(float)
    df["Ion"] = df["Ion"].astype(int)
    return df


def _read_ges_refs() -> dict:
    """GES reference code -> (author, bibcode). The decoder proper."""
    return {r["Ref"]: (r["Aut"], r["BibCode"], r["Com"])
            for r in _read_tsv(GES / "refs.tsv")}


def _read_jofre_codes() -> dict:
    """Jofré integer code -> (table, sources as published, first-author surnames)."""
    return {int(r["code"]): (r["table"], r["sources_as_published"],
                             tuple(r["first_author_surnames"].split(",")))
            for r in _read_tsv(HOLDING / "paper_table45_refcodes.tsv")}


def _fold(text: str) -> str:
    """ASCII-fold, for surname comparison ONLY — never for anything that gets written."""
    return "".join(c for c in unicodedata.normalize("NFD", text)
                   if unicodedata.category(c) != "Mn").lower().replace("'", "")


def _decode_ges_code(code: str, refs: dict) -> tuple[str, tuple]:
    """A GES `r_loggf` -> a human string and the first-author surnames it names.

    Heiter+2021 Note (3): a code may combine labels with `+` (the adopted value is an
    AVERAGE over those sources) or `|` (relative values from the first source put on an
    absolute scale using lifetimes from the second). Both are preserved in the rendered
    string — collapsing a composite to its first label would silently drop a source that
    contributed to the number, which is the opposite of provenance.
    """
    parts, names = [], []
    joiner = "+" if "+" in code else ("|" if "|" in code else "")
    for tok in re.split(r"[+|]", code):
        tok = tok.strip()
        if tok in refs:
            aut, bib, com = refs[tok]
            # Not every refs.dat row has a BibCode — BWL, K07 and K13 carry the citation in
            # the free-text Com field instead. Fall through to it rather than emitting a
            # bare year, which would look like a citation and identify nothing.
            cite = bib or com or tok
            parts.append(f"{aut} [{cite}]")
            names.append(_fold(re.split(r"\s+(?:and|et)\b|,", aut)[0].strip()))
        else:
            parts.append(f"{tok} (NOT IN refs.dat)")
            names.append(f"?{tok}")
    sep = {"+": " + ", "|": " | ", "": ""}[joiner]
    text = sep.join(parts) if len(parts) > 1 else parts[0]
    if joiner == "+":
        text += "  (adopted value is the AVERAGE of these sources)"
    elif joiner == "|":
        text += "  (relative gf from the first, absolute scale from the second)"
    return text, tuple(names)


# ── controls on the published record ─────────────────────────────────────────
def _check_published_counts(ew_rows: list[dict]) -> None:
    """Every star's ew.dat block must equal Table 3's N(Fe I)+N(Fe II) for that star.

    Only the Sun's counts are pinned here (that is the set this ticket ships), but the
    identity is what licenses reading `ew.dat` as "the selected lines" at all, so the
    solar half is asserted rather than assumed.
    """
    got = {sp: sum(1 for r in ew_rows if r["star"] == STAR and r["species"] == sp)
           for sp in PUBLISHED_SUN_COUNTS}
    if got != PUBLISHED_SUN_COUNTS:
        raise BuildError(
            f"solar line counts {got} do not reproduce Jofré Table 3 "
            f"{PUBLISHED_SUN_COUNTS} — the holding is not the published set")


def _check_ep_agreement(sun: list[dict], paper: dict) -> None:
    """VizieR's 2-dp EP against the paper's 4-dp Elow, for every joined line.

    This is the CONTROL on the λ join: two tables that agree on wavelength but describe
    different transitions would disagree here. Nothing else in this build cross-checks
    the paper join, so without it the gf column would rest on a wavelength alone —
    exactly the RYA-853 crosscheck_nist failure.
    """
    bad = [(r["species"], r["wavelength_air_A"], r["excitation_potential_eV"],
            paper[(r["species"], r["wavelength_air_A"])]["elow_eV"])
           for r in sun if (r["species"], r["wavelength_air_A"]) in paper
           and abs(paper[(r["species"], r["wavelength_air_A"])]["elow_eV"]
                   - r["excitation_potential_eV"]) > 0.006]
    if bad:
        raise BuildError(f"VizieR EP and paper Elow disagree on {len(bad)} lines "
                         f"(first: {bad[0]}) — the paper join is not this transition")


def _null_control(sun: list[dict], cg: pd.DataFrame) -> list[int]:
    """Resolved counts for wavelength-DISPLACED GBS lines. Must be all zero."""
    counts = []
    for shift in _NULL_SHIFTS_A:
        n = 0
        for species in PUBLISHED_SUN_COUNTS:
            want = [r for r in sun if r["species"] == species]
            src = cg[cg["species"] == species]
            n += line_match.match(
                [r["wavelength_air_A"] + shift for r in want],
                src["wavelength_air_A"].values,
                want_ep=[r["excitation_potential_eV"] for r in want],
                src_ep=src["excitation_potential_eV"].values,
                tol_A=_MATCH_TOL_A, require_ep=True).n_resolved
        counts.append(int(n))
    return counts


# ── the build ────────────────────────────────────────────────────────────────
def _rew(ew_mA: float, lam_A: float) -> float:
    return math.log10(ew_mA * 1.0e-3 / lam_A)


#: 🔴 SOLAR-FITTED gf SOURCES — a gf calibrated ON the Sun cannot referee a solar
#: abundance (RYA-161). Keyed by the first-author surname the decoders print. Meléndez &
#: Barbuy 2009 (A&A 497, 611) is `melendez2009` in data/refs/bibliography.csv: *"multiplets
#: are globally normalised on laboratory data but individual values are partly solar-fitted,
#: so it must never referee a solar abundance"*. Decoding the provenance is what made this
#: checkable, so the check ships with the decode rather than waiting for a reader to notice.
FIREWALLED_SOURCES = {
    "melendez": "Meléndez & Barbuy 2009 (A&A 497, 611) — partly solar-fitted Fe II gf, "
                "FIREWALLED by RYA-161/852 (bibliography key `melendez2009`)",
}

#: 🔴 RATIFIED DISPOSITION — Ryan, RYA-1110, 2026-08-29. FLAG AND KEEP.
#:
#: The three GBS Fe II lines whose gf decodes to Meléndez & Barbuy 2009 stay in the
#: 142-line replication set. Dropping them would break replication fidelity: we replicate
#: Jofré's PUBLISHED set and flag its properties — the same principle already applied to
#: the −4.8 quirk, where their stated rule and their published selection disagree and BOTH
#: are carried.
#:
#: WHAT THE FLAG MEANS, and it must travel with the line: on these lines the gf was itself
#: calibrated on the Sun, so a GBS solar Fe II value derived from them is a
#: **METHOD-REPRODUCTION CHECK, NOT AN INDEPENDENT VALIDATION**. That sentence is the
#: point of the flag; a flag that only names the paper leaves the reader to rediscover the
#: consequence.
#:
#: ⚠️ NOTHING CONSUMES THIS YET. RYA-1111 is the measurement path that wires `line_set` to
#: products, and it does not exist. `SOLAR_CIRCULAR` below is the importable form so 1111
#: binds to it rather than retyping three wavelengths — the RYA-845 shape (two declarations
#: of one fact) is how a flag like this goes missing.
DISPOSITION = "FLAG-AND-KEEP (Ryan, RYA-1110, 2026-08-29)"

#: The controlled vocabulary for `gbs_solar_validity`. `not-flagged` is deliberately NOT
#: called "independent": this check tests one thing — whether the GBS gf decodes to a
#: solar-fitted source — and passing it is not a certificate of anything else.
SOLAR_VALIDITY = ("not-flagged", "method-reproduction-only")

_CIRCULAR_MEANING = (
    "METHOD-REPRODUCTION CHECK, NOT AN INDEPENDENT VALIDATION: this line's gf was itself "
    "calibrated on the Sun, so a solar Fe II value derived from it reproduces the method "
    "rather than testing it")


def solar_circular_lines(df: pd.DataFrame) -> pd.DataFrame:
    """The rows RYA-1110 flags circular for a SOLAR number, for RYA-1111 to bind to.

    DERIVED from the built line set, never a typed list of wavelengths — the flag and the
    lines it names cannot drift apart if there is only one of them.
    """
    return df[df["gbs_solar_validity"] == "method-reproduction-only"]


def _provenance(gbs_gf, our_gf, paper_row, h, refs: dict, jofre_codes: dict) -> dict:
    """Decode WHERE THE GBS gf VALUE CAME FROM, and say which decoder answered.

    🔴 THE ORDERING IS THE WHOLE POINT. Heiter+2021 gives a per-line source and is the
    finer instrument — but it describes GES **v6**, and Jofré used **v3**. Where the two
    versions carry different log gf, Heiter's `r_loggf` is the provenance of a DIFFERENT
    NUMBER. Attaching it to Jofré's value would manufacture a pedigree: the line would
    look sourced to a 2014 laboratory paper whose value it does not carry. So:

      heiter2021-exact    Heiter's log gf EQUALS the published GBS value -> Heiter's
                          per-line code, decoded through refs.dat. The strong case.
      jofre2014-footnote  the values differ (GES revised it after v3), so the only
                          provenance valid for the GBS value is the paper's own footnote —
                          a source LIST, coarser, but attached to the right number.
      no-gbs-value        Jofré published no gf for this line (Tables 4/5 carry golden
                          lines only). There is nothing to attribute; Heiter's v6 value and
                          source are still recorded, in their own columns.
      unresolved          neither route answers. NEVER guessed.
    """
    h_text, h_names = _decode_ges_code(h["r_loggf"], refs)
    code = None if paper_row is None else paper_row["ref_code"]
    j = jofre_codes.get(code) if code is not None else None
    j_text = "" if j is None else j[1]
    same = gbs_gf is not None and abs(gbs_gf - h["loggf"]) <= _GF_SAME_DEX

    if gbs_gf is None:
        basis = "no-gbs-value"
        text = ("NO GBS gf PUBLISHED (line is not in Jofré Tables 4/5). GES v6 value "
                f"{h['loggf']:+.3f} is from: {h_text}")
        names = ()
    elif same:
        basis = "heiter2021-exact"
        text = (f"{h_text}  [Heiter et al. 2021, A&A 645 A106, geslines.dat r_loggf="
                f"{h['r_loggf']}; GES v6 log gf equals the published GBS value]")
        names = h_names
    elif j is not None:
        basis = "jofre2014-footnote"
        text = (f"{j_text}  [Jofré et al. 2014, A&A 564 A133, {j[0]} footnote, code "
                f"{code}. Heiter+2021 carries {h['loggf']:+.3f} for this line against the "
                f"GBS {gbs_gf:+.3f}, so the GES v6 per-line source describes a different "
                f"number and is recorded separately]")
        names = j[2]
    else:
        basis = "unresolved"
        text = (f"UNRESOLVED — Jofré code {code} is not defined in the published Table 4/5 "
                f"footnote, and Heiter+2021's {h['loggf']:+.3f} differs from the GBS "
                f"{gbs_gf:+.3f} so its per-line source is not this value's source")
        names = ()

    fired = sorted({FIREWALLED_SOURCES[n] for n in names if n in FIREWALLED_SOURCES})
    if fired:
        # Does "use our gf instead" actually escape the circularity on THIS line? Ryan's
        # decision notes it does not on two of the three, and that is a per-line fact, so
        # it is derived per line rather than stated once in prose that a filter cannot see.
        escapes = our_gf is not None and abs(our_gf - gbs_gf) > _GF_SAME_DEX
        fired.append(
            _CIRCULAR_MEANING + ". DISPOSITION " + DISPOSITION + ": the line STAYS in the "
            "replication set — dropping it would break replication fidelity — and carries "
            "this flag wherever it feeds a solar number. " +
            ("Our adopted gf differs, so the our-gf arm of the do-both comparison is not "
             "circular here." if escapes else
             "🔴 OUR ADOPTED gf IS THE SAME NUMBER, so 'use our gf' does NOT escape the "
             "circularity on this line — neither arm of the do-both comparison is "
             "independent here."))
    return {
        "gf_source_per_line": text,
        "gf_source_basis": basis,
        "gf_source_firewalled": "; ".join(fired),
        "gbs_solar_validity": "method-reproduction-only" if fired else "not-flagged",
        "jofre_refcode_sources": j_text,
        "heiter2021_source": h_text,
    }


def build() -> pd.DataFrame:
    ew_rows = _read_measurement_table(VIZIER / "ew.dat")
    ab_rows = _read_measurement_table(VIZIER / "abund.dat")
    _check_published_counts(ew_rows)

    sun = [r for r in ew_rows if r["star"] == STAR]
    sun.sort(key=lambda r: (r["species"], r["wavelength_air_A"]))
    ab = {(r["species"], r["wavelength_air_A"]): r
          for r in ab_rows if r["star"] == STAR}

    # The pairing is built from ALL 242 published lines, not just the Sun's 159 — it is
    # a property of the two TABLES, and restricting it first would make the count check
    # compare a subset against the whole and fail for the wrong reason.
    t6 = _bind_table6(ew_rows, _read_table6(VIZIER / "table6.dat"))
    paper = _read_paper_tables()
    _check_ep_agreement(sun, paper)

    cg = pd.read_csv(ROOT / "data" / "linelists" / "canonical_gf.csv", low_memory=False)
    cg = cg[cg["species"].isin(PUBLISHED_SUN_COUNTS)].reset_index(drop=True)

    nulls = _null_control(sun, cg)
    if any(nulls):
        raise BuildError(
            f"the displaced-null control resolved {nulls} lines at tol "
            f"{_MATCH_TOL_A} Å — a tolerance with a non-zero chance rate is not a "
            f"tolerance (RYA-1070); do not widen it, narrow it")

    # λ+EP dual key into canonical_gf, per species.
    ours: dict = {}
    for species in PUBLISHED_SUN_COUNTS:
        want = [r for r in sun if r["species"] == species]
        src = cg[cg["species"] == species].reset_index(drop=True)
        res = line_match.match([r["wavelength_air_A"] for r in want],
                               src["wavelength_air_A"].values,
                               want_ep=[r["excitation_potential_eV"] for r in want],
                               src_ep=src["excitation_potential_eV"].values,
                               tol_A=_MATCH_TOL_A, require_ep=True)
        # RAISE on an unresolved or ambiguous line rather than letting it become an empty
        # gf column. RYA-833/1033: a line with no atomic-data row that travels on as NaN is
        # indistinguishable in the output from a genuinely low-tier line.
        idx = line_match.require_resolved(res, what="the GBS solar set (RYA-1110)",
                                          species=species)
        for k, r in enumerate(want):
            ours[(r["species"], r["wavelength_air_A"])] = (src.iloc[int(idx[k])],
                                                           float(res.distance_A[k]))

    # ── the gf-PROVENANCE decode (RYA-1110 second pass) ──────────────────────
    ges = _read_ges_lines()
    refs = _read_ges_refs()
    jofre_codes = _read_jofre_codes()
    heiter: dict = {}
    for species, ion in _GES_ION.items():
        want = [r for r in sun if r["species"] == species]
        src = ges[ges["Ion"] == ion].reset_index(drop=True)
        res = line_match.match([r["wavelength_air_A"] for r in want],
                               src["lambda"].values,
                               want_ep=[r["excitation_potential_eV"] for r in want],
                               src_ep=src["Elow"].values,
                               tol_A=_MATCH_TOL_A, require_ep=True)
        idx = line_match.require_resolved(
            res, what="the GBS solar set vs Heiter+2021 geslines (RYA-1110)",
            species=species, source="Heiter+2021 geslines")
        for k, r in enumerate(want):
            heiter[(r["species"], r["wavelength_air_A"])] = (src.iloc[int(idx[k])],
                                                            float(res.distance_A[k]))
    if _null_control(sun, ges.rename(columns={"lambda": "wavelength_air_A",
                                              "Elow": "excitation_potential_eV"})
                     .assign(species=np.where(ges["Ion"] == 1, "Fe I", "Fe II"))) != [0] * len(_NULL_SHIFTS_A):
        raise BuildError("the displaced-null control against Heiter+2021 resolved lines — "
                         "the provenance join has a non-zero chance rate and is not evidence")

    out = []
    for r in sun:
        key = (r["species"], r["wavelength_air_A"])
        lam, ep = r["wavelength_air_A"], r["excitation_potential_eV"]
        rews = [_rew(e, lam) for e in r["values"]]
        cls = ("pass" if max(rews) <= REW_CUT
               else "excluded" if min(rews) > REW_CUT else "ambiguous")
        p = paper.get(key)
        o = ours[key]
        row = o[0]
        gbs_gf = None if p is None else p["log_gf"]
        our_gf = row["log_gf"]
        ges_gf = None if pd.isna(row["gf_synth_ges"]) else float(row["gf_synth_ges"])
        h, hd = heiter[key]
        prov = _provenance(gbs_gf, our_gf, p, h, refs, jofre_codes)
        band = band_policy.resolve(lam).name
        if band != BAND:
            raise BuildError(
                f"{r['species']} {lam} resolves to band {band!r}, not {BAND!r}. The "
                f"ticket scopes this set to VIS; a line outside it is not a rounding "
                f"question, it is a scope question and belongs to Ryan, not to a filter.")
        out.append({
            "line_set": LINE_SET_TAG,
            "band": band,
            "species": r["species"],
            "wavelength_air_A": f"{lam:.2f}",
            "excitation_potential_eV": f"{ep:.2f}",
            # ---- as published by Jofré+2014 ----
            "gbs_golden": t6[key]["golden"],
            "gbs_selected_sun": t6[key]["sun"],
            "elow_eV_paper": "" if p is None else f"{p['elow_eV']:.4f}",
            "log_gf_gbs": "" if gbs_gf is None else f"{gbs_gf:.3f}",
            "vdw_abo_gbs": "" if p is None else f"{p['vdw_abo']:.3f}",
            "loggf_ref_code_gbs": "" if p is None else p["ref_code"],
            # DERIVED, not asserted. It was hardcoded False in the first pass because the
            # arXiv copy did not typeset the footnote; the PUBLISHED PDF does, so this now
            # says whether THIS row's code actually decodes. A status column that keeps
            # saying "unresolved" after the decode lands is the RYA-1110-shaped lie.
            "loggf_ref_gbs_resolved": str(p is not None
                                          and p["ref_code"] in jofre_codes),
            "gf_provenance_gbs": (
                "NOT PUBLISHED IN TABLES 4/5" if p is None else
                f"Jofre+2014 ({'Table 4' if r['species'] == 'Fe I' else 'Table 5'}) "
                f"code {p['ref_code']}; line list GES-v3 (Heiter+ 2014, in prep); "
                f"code decoder did not typeset in arXiv:1309.1099v2"),
            # ---- the published solar measurement the cut is applied to ----
            "n_methods_ew": len(r["values"]),
            "n_methods_abund": len(ab[key]["values"]) if key in ab else "",
            "ew_mA_min": f"{min(r['values']):.2f}",
            "ew_mA_max": f"{max(r['values']):.2f}",
            "ew_mA_mean": f"{statistics.fmean(r['values']):.2f}",
            "rew_min": f"{min(rews):.4f}",
            "rew_max": f"{max(rews):.4f}",
            "rew_mean": f"{_rew(statistics.fmean(r['values']), lam):.4f}",
            "rew_class": cls,
            # ---- our pool, joined on λ+EP ----
            "our_wavelength_air_A": f"{row['wavelength_air_A']:.4f}",
            "our_excitation_potential_eV": f"{row['excitation_potential_eV']:.4f}",
            "match_distance_mA": f"{o[1] * 1000:.2f}",
            "log_gf_ours": f"{our_gf:.4f}",
            "loggf_reference_ours": row["loggf_reference"],
            "gf_tier_ours": row["gf_tier"],
            "nist_grade_ours": "" if pd.isna(row["nist_grade"]) else row["nist_grade"],
            "gf_synth_ges": "" if ges_gf is None else f"{ges_gf:.4f}",
            "delta_gbs_minus_ours": ("" if gbs_gf is None
                                     else f"{gbs_gf - our_gf:+.4f}"),
            "delta_gbs_minus_ges": ("" if gbs_gf is None or ges_gf is None
                                    else f"{gbs_gf - ges_gf:+.4f}"),
            # ---- gf PROVENANCE, decoded (RYA-1110 second pass) ----
            "gf_source_per_line": prov["gf_source_per_line"],
            "gf_source_basis": prov["gf_source_basis"],
            "gf_source_firewalled": prov["gf_source_firewalled"],
            "gbs_solar_validity": prov["gbs_solar_validity"],
            "jofre_refcode_sources": prov["jofre_refcode_sources"],
            "heiter2021_r_loggf": h["r_loggf"],
            "heiter2021_source": prov["heiter2021_source"],
            "heiter2021_log_gf": f"{h['loggf']:.3f}",
            "heiter2021_gfflag": h["gfflag"],
            "heiter2021_synflag": h["synflag"],
            "heiter2021_match_distance_mA": f"{hd * 1000:.2f}",
            "delta_gbs_minus_heiter2021": ("" if gbs_gf is None
                                           else f"{gbs_gf - h['loggf']:+.4f}"),
            # ---- reachability ----
            "telluric_exclusion": telluric_exclusion(lam) or "",
        })
    return pd.DataFrame(out)


# ── coverage ─────────────────────────────────────────────────────────────────
#: VERIFIED spans for the two Kitt Peak 1984 holdings, whose HoldingSpec declares
#: `span_A=None` because the reader inventories its own segments. `covers()` therefore
#: answers True for ANY wavelength there, which is not a coverage claim — it is the
#: absence of one (RYA-767). The number used instead is the VERIFIED span in
#: data/catalog/solar_reference_holdings_rya708.csv, and the report says which it used.
_REGISTRY_SPAN = {"kpno_solar_atlas": (2960.0, 13000.0)}


def coverage(df: pd.DataFrame) -> pd.DataFrame:
    from measure_band_ew import _INSTRUMENT_HOLDINGS

    vis = next(b for b in band_policy.POLICIES if b.name == BAND)
    sel = df[df["rew_class"] == "pass"]
    lam = sel["wavelength_air_A"].astype(float).to_numpy()
    # RYA-1110 disposition: the flag has to travel to anything that could feed a solar
    # number, and the coverage report is the only such surface that exists today. A
    # holding that reaches these lines reaches a method-reproduction check, not an
    # independent one, and the row says so rather than leaving it two files away.
    circ = (sel["gbs_solar_validity"] == "method-reproduction-only").to_numpy()

    rows = []
    for instrument, specs in _INSTRUMENT_HOLDINGS.items():
        for h in specs:
            if h.span_A is not None:
                lo, hi = h.span_A
                source = "HoldingSpec.span_A"
            elif instrument in _REGISTRY_SPAN:
                lo, hi = _REGISTRY_SPAN[instrument]
                source = "solar_reference_holdings_rya708.csv (VERIFIED)"
            else:
                lo = hi = float("nan")
                source = "UNDECLARED — reader inventories its own coverage"
            if math.isnan(lo):
                inside = np.zeros(lam.shape, dtype=bool)
            else:
                inside = (lam >= lo) & (lam <= hi)
            tell = np.array([bool(telluric_exclusion(float(w), instrument))
                             for w in lam])
            rows.append({
                "line_set": LINE_SET_TAG,
                "band": BAND,
                "band_lo_A": f"{vis.lo_A:.1f}",
                "band_hi_A": f"{vis.hi_A:.1f}",
                "instrument": instrument,
                "holding_id": h.holding_id,
                "span_lo_A": "" if math.isnan(lo) else f"{lo:.2f}",
                "span_hi_A": "" if math.isnan(hi) else f"{hi:.2f}",
                "span_source": source,
                "n_lines_selected": int(len(lam)),
                "n_in_span": int(inside.sum()),
                "n_telluric_excluded": int((inside & tell).sum()),
                "n_reachable": int((inside & ~tell).sum()),
                "n_reachable_solar_circular": int((inside & ~tell & circ).sum()),
                "pre_normalised": str(h.pre_normalised),
            })
    return pd.DataFrame(rows)


def main() -> int:
    df = build()
    OUT.mkdir(parents=True, exist_ok=True)
    df.to_csv(LINESET, index=False, lineterminator="\n")
    cov = coverage(df)
    cov.to_csv(COVERAGE, index=False, lineterminator="\n")

    n = len(df)
    counts = df["rew_class"].value_counts().to_dict()
    print(f"{LINESET.relative_to(ROOT)}  {n} rows "
          f"({int((df.species == 'Fe I').sum())} Fe I + "
          f"{int((df.species == 'Fe II').sum())} Fe II) — reproduces Jofré Table 3's "
          f"published solar N = {PUBLISHED_SUN_COUNTS}")
    print(f"  REW <= {REW_CUT}: pass {counts.get('pass', 0)}  "
          f"ambiguous {counts.get('ambiguous', 0)}  excluded {counts.get('excluded', 0)}")
    lo = df["rew_min"].astype(float).min()
    hi = df["rew_max"].astype(float).max()
    print(f"  REW range {lo:.4f} .. {hi:.4f}   "
          f"EP range {df.excitation_potential_eV.astype(float).min():.2f} .. "
          f"{df.excitation_potential_eV.astype(float).max():.2f} eV   "
          f"lambda {df.wavelength_air_A.astype(float).min():.2f} .. "
          f"{df.wavelength_air_A.astype(float).max():.2f} A")
    print(f"  gf published by Jofré: {int((df.log_gf_gbs != '').sum())}/{n};  "
          f"joined to canonical_gf: {int((df.log_gf_ours != '').sum())}/{n}")
    print(f"{COVERAGE.relative_to(ROOT)}  {len(cov)} holdings")
    for _, r in cov[cov.n_reachable.astype(int) > 0].iterrows():
        print(f"  {r.holding_id:38s} {r.n_reachable:3d}/{r.n_lines_selected}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
