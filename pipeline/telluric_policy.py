"""Telluric policy — the single source for a question that keeps getting re-asked (RYA-786).

THE RECURRENCE THIS ENDS
------------------------
The telluric decision for a reference atlas was already made and built three times over:

  * **RYA-424** — telluric correction as a standing data-input stage, with
    instrument-aware routing (molecfit / cr2res / APERO) and an analysis-ready flag.
  * **`data/catalog/instrument_catalog.csv`** — `kpno_solar_atlas` is registered
    `telluric_required = no`; **RYA-380** is the molecfit/GDAS recipe for the instruments
    that do need it.
  * **RYA-460 / `config/physics_regime_rya400.yaml`** — per-line KPNO handling, e.g.
    K I "7665 stays in the O2 A-band, 7699 is the clean line".

Despite that, `SynthesisHandler.prepare` refused a KPNO run on a BAND-level
`telluric_required=True`, so every new band run re-collided with a settled question. That
is a single-source-of-truth defect, not a science gap, and this module is the source.

THE DISTINCTION THAT KEEPS GETTING LOST
---------------------------------------
`telluric_required = no` for a reference atlas does **NOT** mean "this atlas has been
telluric-divided". The Kurucz 1984 KPNO atlas HAS telluric absorption in it. It means
**the tellurics are handled by per-line CLEAN-LINE SELECTION rather than by a correction
stage** — the standard method for a reference atlas that has a clean alternative line and
a second arm (IAG) to cross-check against.

So the honest basis for running KPNO is *"instrument flag + per-line selection"*, cited.
It is NOT a `telluric_corrected` declaration, and fabricating one to satisfy a gate is
forbidden (RYA-786): it asserts a correction that was never applied.

WHAT THIS MODULE DECIDES
------------------------
  1. `requires_correction(instrument)` — from the catalog, the single registry of what an
     instrument needs. `yes` routes to the RYA-380/424 molecfit path; the correction
     machinery stays real and nothing is "avoided" as architecture.
  2. `exclusion(wave_A)` — the O2/H2O band set, enumerated ONCE below. A line inside a
     band is QUARANTINED-TELLURIC: a valid physics exclusion (RYA-777), not a cull.
  3. `gate(instrument, analysis_ready)` — what a handler should ask instead of carrying
     its own band flag.
  4. `applied_state(holding)` / `gate_holding(holding)` — the SECOND axis (RYA-806).
  5. `correction_requirement(instrument)` / `reconcile_axes(instrument)` — the THREE-way
     read of the requirement, and the check that the catalog's two telluric columns do
     not state opposite rules (RYA-1072).

THE TWO AXES (RYA-806) — ORTHOGONAL, NEVER COLLAPSED
-----------------------------------------------------
    telluric_basis    per-INSTRUMENT   does this band NEED correction?
                      `instrument_catalog.csv`, RYA-786
    telluric_applied  per-HOLDING      has correction been APPLIED to THIS product?
                      `holdings_manifest_registry.csv`, RYA-806

The second is NOT derivable from the first: it is a fact about the downloaded PRODUCT
LEVEL, and the same exposures can exist at two levels. Both are consumed together in
`gate_holding` so no caller has to combine them itself — combining them by hand in three
places is how the RYA-786 defect happened in the first place.

    state         instrument needs it?    outcome
    applied       either                  serve as-is, SKIP the RYA-424 stage
    not-applied   yes                     REFUSE -> route through RYA-424
    not-applied   no  (line_selection)    run, on the per-line selection basis
    not-applied   UNRESOLVED              REFUSE, always (`TelluricStateUnknown`)
    unknown       either                  REFUSE, always (`TelluricStateUnknown`)

⚠️ THE INSTRUMENT COLUMN ABOVE HAS THREE VALUES, NOT TWO (RYA-1072)
-------------------------------------------------------------------
"does this band need correction?" answers REQUIRED, NOT_REQUIRED **or AMBIGUOUS**, and
that third row is not decoration — `instrument_catalog.csv` carries
`telluric_required=mode_dependent` for UVES today. `requires_correction` used to test
`v in ('yes','true','1','required')` and return False for everything else, so the
ambiguous case took the identical branch to a declared `no` and `gate_holding` served
the holding on the clean-line basis with no mode ever consulted. Reading a three-valued
fact through a two-valued test, with the unresolved case landing on the permissive side,
is the same shape as the RYA-786 defect this module was written to end.

    requires_correction(instrument)     -> bool, RAISES on AMBIGUOUS
    correction_requirement(instrument)  -> REQUIRED | NOT_REQUIRED | AMBIGUOUS
    reconcile_axes(instrument)          -> raises if the two columns contradict

Use the three-outcome accessor where a caller must survive an unresolved instrument (a
sweep over every registered holding, a status surface); use the boolean where the answer
is about to decide whether data may be measured, because there it must refuse.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
CATALOG = ROOT / "data" / "catalog" / "instrument_catalog.csv"
HOLDINGS = ROOT / "data" / "catalog" / "holdings_manifest_registry.csv"

# The per-holding axis (RYA-806). Anything else in the column reads as `unknown`, so a
# typo refuses rather than silently passing.
_APPLIED_VALUES = ("applied", "not-applied", "unknown")

# ── THE authoritative telluric exclusion set (RYA-786) ───────────────────────
#
# Enumerated once, here, and consumed everywhere. It was previously a three-entry list
# inside a measurement script — O2 A-band 7600-7640 plus two H2O bands — which is both
# incomplete and too narrow: the A-band runs to ~7685, the O2 B-band was absent entirely,
# and two H2O complexes in the red-optical/NIR were missing. A line sitting in an
# unlisted band is measured as if it were clean, which is the silent version of this bug.
#
# Ranges are air wavelengths in Angstrom, inclusive.
TELLURIC_BANDS: tuple[tuple[float, float, str], ...] = (
    (6867.0, 6884.0, "O2 B-band"),
    (7160.0, 7340.0, "H2O"),
    (7594.0, 7685.0, "O2 A-band"),
    (8100.0, 8400.0, "H2O"),
    (9280.0, 9600.0, "H2O"),
    (11120.0, 11560.0, "H2O"),
)

QUARANTINE_TAG = "QUARANTINED-TELLURIC"


class TelluricStateUnknown(RuntimeError):
    """We hold this product in an UNVERIFIED telluric state, OR the catalog does not
    resolve whether this instrument needs a correction stage at all.

    Deliberately distinct from "this band needs correction" and from "we lack this
    band" — the RYA-796 `RestFrameNotConditioned` discipline: a refusal names the state
    it is refusing, so the fix is unambiguous. `unknown` is never defaulted either way:
    defaulting to `applied` fabricates a correction (forbidden, RYA-786) and defaulting
    to `not-applied` sends a corrected product through a second correction.

    RYA-1072 widened it to the INSTRUMENT axis. It was already the right refusal for an
    undetermined per-holding state; an undetermined per-instrument REQUIREMENT is the
    same kind of not-knowing and earns the same loud stop.
    """


class TelluricCatalogContradiction(RuntimeError):
    """`telluric_required` and `telluric_basis` name OPPOSITE operative rules.

    A separate type from `TelluricStateUnknown` on purpose: not-knowing is fixed by
    DETERMINING a value, and this is fixed by deciding WHICH OF TWO RECORDED VALUES is
    stale. Serving the holding on either reading would settle that question silently
    (RYA-1069 ledger warning 2 — `delbouille_liege_intensity`).
    """


# ── THE INSTRUMENT AXIS HAS THREE OUTCOMES, NOT TWO (RYA-1072) ───────────────
#
# 🔴 THE DEFECT THIS ENDS. `requires_correction` was `return v in ('yes','true','1',
# 'required')`. That is an allow-list for the TRUE case with EVERYTHING ELSE falling into
# FALSE — so `mode_dependent`, a value the catalog actually uses and which means "the
# answer depends on a mode nobody has resolved", took the identical branch to a declared
# `no`. `gate_holding` then served UVES on the clean-line-selection basis without any
# mode ever being consulted. UVES RED860 reaches 10427 A, well inside the registered H2O
# complexes, so the first UVES product would have been measured under-corrected with no
# refusal anywhere in the chain.
#
# The failure shape is the recurring one: a two-valued reading imposed on a three-valued
# fact, with the ambiguous case laundered into the permissive answer. `pipeline/telluric/
# routing.py` (RYA-927) already got this right on the BASIS axis — `mode_dependent` and
# `unspecified` fall through to `audit_required` there — so the repo held the correct
# reading in one module and the collapsed one in this module's own gate.
#
# BOTH SETS ARE CLOSED AND EXPLICIT, and anything outside them is AMBIGUOUS — that is
# what stops a NEW catalog value from silently joining the permissive side.
# ⚠️ DO NOT "fix" a future refusal by adding its value to `_NOT_REQUIRED_VALUES`.
# Widening the recognised set to swallow an ambiguity is the original bug wearing a patch.
REQUIRED = "required"
NOT_REQUIRED = "not-required"
AMBIGUOUS = "ambiguous"

# ── ONE VOCABULARY FOR BOTH TELLURIC COLUMNS (RYA-1078) ─────────────────────
#
# 🔴 RYA-1072 fixed how THIS module reads an ambiguous value and left
# `pipeline/telluric/routing.py` reading it separately. Two modules implementing "how to
# read an ambiguous telluric value" is the defect one level up from the defect 1072 fixed:
# `unspecified` resolved to audit-required in routing and to run-on-clean-line-selection
# here, and `route_for` settled the delbouille contradiction by preferring the requirement
# column while `reconcile_axes` refused it. The reading now lives HERE, once, and routing
# CALLS it -- so a third caller cannot reintroduce the split by writing a fourth `in`-test.
#
# The two columns share one vocabulary because they answer the SAME question in different
# words -- is a correction STAGE owed:
#
#     telluric_required   yes | no                     | mode_dependent
#     telluric_basis      correction_required          | line_selection | corrected
#                                                      | not_applicable | unspecified
#
# so `correction_required` classifies REQUIRED exactly as `yes` does, and `line_selection`
# / `corrected` / `not_applicable` classify NOT_REQUIRED exactly as `no` does. That shared
# reading is what makes `reconcile_axes` expressible as "the two columns must not classify
# differently", and it is why `unspecified` and `mode_dependent` cannot diverge again: they
# are in neither set, in one place, for both columns.
#
# ⚠️ BOTH SETS STAY CLOSED. Anything outside them is AMBIGUOUS. Do NOT add
# `unspecified` or `mode_dependent` to either set to quiet a refusal -- that is the
# RYA-1072 bug wearing a patch, and now it would be the bug in two modules at once.
_REQUIRED_VALUES = frozenset({
    "yes", "true", "1", "required",     # telluric_required
    "correction_required",              # telluric_basis
})
_NOT_REQUIRED_VALUES = frozenset({
    "no", "false", "0", "not-required", "not_required",   # telluric_required
    "line_selection", "corrected", "not_applicable",      # telluric_basis
})


def classify(value) -> str:
    """THE reading of a telluric catalog value -- either column. RYA-1078.

    REQUIRED / NOT_REQUIRED / AMBIGUOUS. Pure, so the three-way logic is testable without
    a catalog and so a caller that must BRANCH on ambiguity (rather than stop on it) can
    do so without catching an exception.

    This is the single entry point the spec requires: `pipeline/telluric/routing.py`
    calls it rather than re-deriving the reading, so `mode_dependent` and `unspecified`
    resolve identically in both modules BY CONSTRUCTION rather than by two lists agreeing.
    """
    v = str(value if value is not None else "").strip().lower()
    if v in _REQUIRED_VALUES:
        return REQUIRED
    if v in _NOT_REQUIRED_VALUES:
        return NOT_REQUIRED
    return AMBIGUOUS


def correction_requirement_from_value(value) -> str:
    """RYA-1072's name for `classify`, kept so its call sites and tests still read.

    One implementation, two names: the 1072 name says WHICH COLUMN it was written for,
    the 1078 name says it reads either. Aliased rather than duplicated -- a second body
    here would be the exact split this ticket closed.
    """
    return classify(value)


def requires_correction_from_value(value) -> bool:
    """True / False for a RESOLVED value; LOUD on an ambiguous one.

    `mode_dependent`, `unspecified`, a blank cell and any unrecognised string all raise.
    There is no third boolean, and manufacturing one by returning False is the defect.
    """
    if correction_requirement_from_value(value) == AMBIGUOUS:
        raise TelluricStateUnknown(
            f"telluric_required={str(value)!r} is a DECLARED AMBIGUITY, not a declared "
            f"answer: it resolves to neither {sorted(_REQUIRED_VALUES)} nor "
            f"{sorted(_NOT_REQUIRED_VALUES)}. It must NOT be read as 'no correction "
            f"needed' — reading it that way is what served UVES (telluric_required="
            f"mode_dependent, RED860 reaching 10427 A) on the clean-line basis with no "
            f"mode ever consulted (RYA-1072). Resolve it on the instrument/mode axis and "
            f"record the resolved value; do not widen the recognised sets to absorb it.")
    return correction_requirement_from_value(value) == REQUIRED


def basis(instrument: str) -> str:
    """`telluric_basis` for this instrument, from the catalog. Loud on unknown."""
    if "df" not in _catalog_cache:
        _catalog_cache["df"] = pd.read_csv(CATALOG)
    df = _catalog_cache["df"]
    hit = df[df.instrument_id.astype(str) == str(instrument)]
    if not len(hit):
        raise KeyError(
            f"instrument {instrument!r} is not in {CATALOG.name}; its telluric basis is "
            f"unknown and must not be assumed (RYA-786).")
    return str(hit.iloc[0].get("telluric_basis", "unspecified")).strip().lower()


def exclusion(wave_A: float, instrument: str | None = None) -> str:
    """Reason string if this line must be excluded for tellurics, else ''.

    ⚠️ SUPERSEDED AS A DECISION INPUT (RYA-1079), and kept as the enumeration it always
    was. The band-membership rule below answers "is this line inside a registered telluric
    complex", which was the right question while correction was expensive and rare. It is
    not the right question now: it excludes a 0.03-deep line and a saturated core alike,
    and on our own Kitt Peak red-optical pool that throws away 12 of 32 in-band graded
    Fe I lines that measure out RECOVERABLE.

    The successor is `pipeline.telluric_observability`, which decides per line from the
    MEASURED transmission depth over that line's own window -- CLEAN / RECOVERABLE /
    SATURATED -- and routes the recoverable class to correction instead of to the bin.
    `telluric_basis = line_selection` is descriptive from here on; it does not decide.

    This function is NOT removed and its callers are NOT rewired here: the two remaining
    consumers sit inside the measurement path, and moving them is a measurement change,
    not a policy change. What RYA-1079 settles is which of the two answers is
    authoritative when they disagree -- the measured one.

    ⚠️ THE BASIS DECIDES, NOT THE WAVELENGTH. `telluric_required = no` was carrying two
    states that behave OPPOSITELY at the line level, and treating them alike threw away
    real data:

        corrected        the tellurics are REMOVED in the data product, so a line inside
                         an O2/H2O band is ordinary solar spectrum and IS measurable
        line_selection   the tellurics are PRESENT and we exclude the affected lines

    Measured on our own two solar atlases, same windows:

        region            KP min  KP mean  KP<0.5    IAG min  IAG mean  IAG<0.5
        O2 A-band core    -0.003    0.453   51.3%      0.458     0.988     0.1%
        H2O 9280-9600     -0.001    0.663   23.1%      0.547     0.971     0.0%
        clean continuum    0.396    0.954    0.2%      0.462     0.971     0.1%

    Kitt Peak flux is driven to ZERO inside the bands — saturated telluric absorption, not
    a stellar spectrum — so excluding there is correct. IAG sits at continuum in the same
    bands while showing the SAME solar line depths in clean regions, so its correction
    removed tellurics without erasing solar structure and those lines are usable.

    An instrument whose basis we do not positively know is treated as `line_selection`:
    excluding a good line costs coverage, measuring a telluric one costs the number.
    """
    for lo, hi, name in TELLURIC_BANDS:
        if lo <= float(wave_A) <= hi:
            b = basis(instrument) if instrument else "unspecified"
            if b in ("corrected", "not_applicable"):
                return ""
            return (f"{QUARANTINE_TAG}: inside the {name} ({lo:.0f}-{hi:.0f} A) and "
                    f"{instrument or 'this instrument'} is telluric_basis={b}, so the "
                    f"observed flux there is not stellar. Excluded by per-line selection "
                    f"(RYA-460/786), not corrected.")
    return ""


def in_telluric_band(wave_A: float, instrument: str | None = None) -> bool:
    return bool(exclusion(wave_A, instrument))


_catalog_cache: dict = {}


def _catalog_row(instrument: str):
    """This instrument's catalog row. Loud on an unknown id.

    The telluric state of an instrument the registry does not know cannot be asserted,
    and guessing it is how the fabricated declaration got in.
    """
    if "df" not in _catalog_cache:
        _catalog_cache["df"] = pd.read_csv(CATALOG)
    df = _catalog_cache["df"]
    hit = df[df.instrument_id.astype(str) == str(instrument)]
    if not len(hit):
        raise KeyError(
            f"instrument {instrument!r} is not in {CATALOG.name}; its telluric requirement "
            f"is unknown and must not be assumed. Register it first (RYA-786).")
    return hit.iloc[0]


def correction_requirement(instrument: str) -> str:
    """REQUIRED / NOT_REQUIRED / AMBIGUOUS for this instrument, from the catalog.

    The three-outcome accessor. Use it wherever a caller must KEEP GOING past an
    unresolved instrument — a sweep over every registered holding, a status surface —
    and report the ambiguity rather than stop on it. Where the answer is about to decide
    whether data may be measured, use `requires_correction`, which refuses.
    """
    return correction_requirement_from_value(
        _catalog_row(instrument).get("telluric_required", ""))


def requires_correction(instrument: str) -> bool:
    """Does this instrument need a telluric CORRECTION STAGE? From the catalog only.

    🔴 RAISES `TelluricStateUnknown` on an AMBIGUOUS requirement (RYA-1072). It used to
    return False there, which is the CRITICAL defect this signature now makes impossible:
    a caller can no longer receive "no correction needed" for an instrument whose
    requirement nobody has determined. A caller that must survive an unresolved
    instrument asks `correction_requirement()` for the three-way answer instead.

    Loud on an unknown instrument id, as before.
    """
    return requires_correction_from_value(
        _catalog_row(instrument).get("telluric_required", ""))


def reconcile_axes(instrument: str) -> None:
    """Raise if `telluric_required` and `telluric_basis` contradict each other.

    The two columns answer different questions — DOES this band need a correction stage,
    and BY WHAT METHOD are its tellurics handled — but they are not independent: a basis
    of `line_selection`, `corrected` or `not_applicable` asserts no correction stage is
    owed, and `correction_required` asserts one is. When that assertion disagrees with
    `telluric_required`, the row states two opposite operative rules and a consumer
    silently picks one by reading whichever column it happens to consult first.

    Measured, not supposed: across `instrument_catalog.csv` exactly one row disagrees —
    `delbouille_liege_intensity` is `telluric_required=yes` with
    `telluric_basis=line_selection`, while every other `yes` row is
    `correction_required` and every other `line_selection` row is `no`. This function
    does NOT decide which of its two values is stale. That is a data decision, and
    guessing it here would be the RYA-786 fabrication in a new place (RYA-1069 ledger
    warning 2). It refuses, names both values, and hands the choice to a human.

    A basis outside both sets (`unspecified`, `mode_dependent`) asserts nothing, so it
    cannot contradict anything and is left to the requirement axis to refuse or permit.
    """
    row = _catalog_row(instrument)
    raw_req = str(row.get("telluric_required", "")).strip()
    raw_basis = str(row.get("telluric_basis", "")).strip()
    req, b = classify(raw_req), classify(raw_basis)

    # RYA-1078: expressed as ONE comparison over the shared reading, rather than two
    # hand-written membership tests per direction. A column that classifies AMBIGUOUS
    # asserts nothing, so it cannot contradict anything -- that is what keeps UVES
    # (mode_dependent on BOTH columns) a consistent ambiguity rather than a contradiction.
    if AMBIGUOUS in (req, b) or req == b:
        return
    raise TelluricCatalogContradiction(
        f"{instrument}: telluric_required={raw_req!r} reads {req} and "
        f"telluric_basis={raw_basis!r} reads {b}. The two columns state OPPOSITE "
        f"operative rules about whether a correction stage is owed, and this row states "
        f"both — so serving it would settle the question by whichever column the consumer "
        f"read first, which is precisely what `route_for` used to do (RYA-1078). Decide "
        f"which value is stale and record it in {CATALOG.name}; it is not derivable from "
        f"either column (RYA-1069/RYA-1072).")


def gate(instrument: str, analysis_ready: bool = False) -> tuple[bool, str]:
    """(may_run, basis). What a handler asks instead of carrying a band flag.

    A reference atlas registered `telluric_required=no` runs on the per-line selection
    basis. An instrument that DOES require correction runs only once the RYA-424
    analysis-ready flag says the correction was applied and verified.

    🔴 REFUSES (RYA-1072) rather than returning a verdict when the catalog does not
    resolve the requirement, or when its two telluric columns contradict each other.
    Both were previously served as "no correction needed".
    """
    reconcile_axes(instrument)
    if not requires_correction(instrument):
        return True, (f"{instrument} is registered telluric_required=no; tellurics are "
                      f"handled by per-line clean-line selection over "
                      f"{len(TELLURIC_BANDS)} enumerated bands (RYA-460/786), not by a "
                      f"correction stage. No telluric_corrected declaration is made.")
    if analysis_ready:
        return True, (f"{instrument} requires correction and the RYA-424 analysis-ready "
                      f"flag is set: molecfit/GDAS applied and verified (RYA-380).")
    return False, (f"{instrument} is registered telluric_required=yes and the RYA-424 "
                   f"analysis-ready flag is not set. Route it through the RYA-380 "
                   f"molecfit path; before correction the observed flux is not stellar.")


# ── THE SECOND AXIS: per-HOLDING telluric_applied (RYA-806) ──────────────────
#
# ⚠️ THESE TWO AXES ARE ORTHOGONAL AND MUST NEVER BE COLLAPSED.
#
#   telluric_basis   (RYA-786, instrument_catalog.csv)  per-INSTRUMENT:
#                    does this band NEED correction?
#   telluric_applied (RYA-806, holdings_manifest_registry.csv)  per-HOLDING:
#                    has correction been APPLIED to THIS product?
#
# The second is not derivable from the first. The same exposures can exist as a raw
# cr2res IDP and as a corrected ADP, and the instrument axis cannot tell them apart.
# Measured, not argued: alpha Cen CRIRES+ and alpha Cen NIRPS are BOTH registered
# `telluric_required=yes`, and their products differ — CRIRES+ `OBS_NODDING_EXTRACTC_IDP`
# is not corrected, NIRPS `S1D_FINAL_A` is (its `FLUX_TELL_*` / `FLUX_EL` ratio tracks
# `1/ATM_TRANSM` at r=1.0000). Collapsing the axes would either refuse good NIRPS data or
# run a second correction over it.

def applied_state(holding_id: str) -> str:
    """`telluric_applied` for this holding, from the registry only. Loud on unknown id.

    Determined AT INTAKE from the product's own headers — see `telluric_intake` — never
    inferred from the instrument, which answers a different question.
    """
    if "holdings" not in _catalog_cache:
        _catalog_cache["holdings"] = pd.read_csv(HOLDINGS)
    df = _catalog_cache["holdings"]
    hit = df[df.holding_id.astype(str) == str(holding_id)]
    if not len(hit):
        raise KeyError(
            f"holding {holding_id!r} is not in {HOLDINGS.name}; its telluric state is "
            f"unknown and must not be assumed. Register it at intake first (RYA-806).")
    v = str(hit.iloc[0].get("telluric_applied", "")).strip().lower()
    return v if v in _APPLIED_VALUES else "unknown"


def gate_holding(holding_id: str, instrument: str | None = None) -> tuple[bool, str]:
    """(may_run, reason) for a SPECIFIC held product. Consumes BOTH axes, in one place.

        applied      -> serve as-is; skip the RYA-424 telluric stage
        not-applied  -> route through the RYA-424 engine first, IF the band needs it
        unknown      -> refuse, always

    Note the asymmetry that keeps the axes honest: `not-applied` is only a refusal when
    the INSTRUMENT axis says the band needs correction. Kitt Peak is `not-applied` and
    runs anyway, on the per-line clean-line selection basis it was always on (RYA-460/786)
    — it is not "uncorrected data sneaking through", it is a different, stated method.
    """
    if "holdings" not in _catalog_cache:
        _catalog_cache["holdings"] = pd.read_csv(HOLDINGS)
    df = _catalog_cache["holdings"]
    hit = df[df.holding_id.astype(str) == str(holding_id)]
    if not len(hit):
        raise KeyError(
            f"holding {holding_id!r} is not in {HOLDINGS.name}; register it at intake "
            f"(RYA-806) rather than assuming its telluric state.")
    inst = instrument or str(hit.iloc[0].instrument_id)
    state = applied_state(holding_id)
    b = basis(inst)

    # ── the INSTRUMENT axis is consulted FIRST, and can short-circuit ────────
    # Order matters, and this is not `unknown` leaking through a default. Where the band
    # needs no correction stage, the applied-state is not load-bearing: the ratified
    # method is per-line clean-line selection (RYA-460/786), which is defined ON
    # uncorrected data and so is valid whatever the product turns out to be. Refusing
    # there would ground HARPS/ESPRESSO — instruments with no telluric question — on a
    # switch built for the IR. The ticket's own smoke test scopes the raise to "an IR
    # dataset", which is exactly what `requires_correction` names.
    # ── the CATALOG must agree with itself before it can gate anything (RYA-1072) ──
    # Ordered first deliberately: a contradictory row would otherwise be resolved by
    # whichever of its two columns the branches below happen to consult, which is
    # exactly the silent settling this refusal exists to prevent.
    reconcile_axes(inst)

    if b == "not_applicable":
        return True, (f"{holding_id}: {inst} is telluric_basis=not_applicable (no "
                      f"terrestrial atmosphere in this data), so telluric_applied="
                      f"{state} does not gate it.")

    # ── an UNRESOLVED requirement is a refusal, never a pass (RYA-1072) ──────
    # 🔴 THE CRITICAL FIX. `requires_correction` used to return False here for
    # `mode_dependent`, so control fell into the clean-line-selection branch below and
    # this function returned (True, "registered telluric_required=no ...") for an
    # instrument whose requirement nobody had determined — a fabricated permission
    # quoting a value the catalog does not contain. It is named, not merely denied, so
    # the fix is unambiguous: RESOLVE THE MODE.
    # RYA-1078: EITHER column being unresolved is an unresolved instrument. Before this,
    # `unspecified` basis fell straight through to the clean-line branch here while
    # `routing.route_for` sent it to `audit_required` -- the same value, two answers, in
    # the two modules that both claim to own this decision.
    if AMBIGUOUS in (correction_requirement(inst), classify(b)):
        raise TelluricStateUnknown(
            f"{holding_id}: {inst} is telluric_basis={b!r} with telluric_required="
            f"{str(_catalog_row(inst).get('telluric_required', '')).strip()!r} — reading "
            f"{classify(b)} and {correction_requirement(inst)}. At least one column "
            f"resolves to neither required nor not-required. This is NOT 'we hold it "
            f"uncorrected' (route it through RYA-424) and NOT 'this band needs no "
            f"correction' (run it on clean-line selection) — it is that nobody has "
            f"determined which of those two this instrument/mode is. Serving it on the "
            f"clean-line basis is what RYA-1072 fixed: UVES RED860 reaches 10427 A, "
            f"inside the registered H2O complexes. Resolve the mode "
            f"(instrument_modes.csv) and record the value.")

    if not requires_correction(inst):
        extra = ("" if state != "unknown" else
                 " Its telluric_applied is UNKNOWN and stays unknown — that is recorded, "
                 "not resolved, and it costs only coverage: per-line selection may "
                 "exclude a line this product had already corrected (RYA-786).")
        return True, (f"{holding_id}: {inst} is registered telluric_required=no on basis "
                      f"'{b}', so tellurics are handled by per-line clean-line selection "
                      f"over {len(TELLURIC_BANDS)} enumerated bands (RYA-460/786), not by "
                      f"a correction stage. telluric_applied={state}. No "
                      f"telluric_corrected declaration is made.{extra}")

    # ── the band DOES need correction: the holding axis now decides ──────────
    if state == "unknown":
        raise TelluricStateUnknown(
            f"holding {holding_id!r} ({inst}) has telluric_applied=unknown and {inst} is "
            f"registered telluric_required=yes: we hold an IR product in an UNVERIFIED "
            f"telluric state. This is not 'we lack this band' (a coverage answer) and not "
            f"'this band needs correction' (an instrument answer) — it is that nobody has "
            f"determined, from this product's own headers, whether a correction was "
            f"already applied. Determine it with "
            f"pipeline.telluric_intake.from_headers() and record the value on "
            f"{HOLDINGS.name}. It is never defaulted: assuming 'applied' fabricates a "
            f"correction (forbidden, RYA-786), assuming 'not-applied' risks correcting an "
            f"already-corrected product twice.")

    if state == "applied":
        return True, (f"{holding_id}: telluric_applied=applied — the correction is "
                      f"already in this product, so the RYA-424 stage is SKIPPED. "
                      f"(Instrument axis: {inst} telluric_basis={b}.)")

    return False, (f"{holding_id}: telluric_applied=not-applied and {inst} is registered "
                   f"telluric_required=yes (basis '{b}'). Route it through the RYA-424 "
                   f"telluric stage (molecfit / cr2res / APERO+Wapiti) BEFORE any IR "
                   f"abundance; before correction the observed flux is not stellar.")
