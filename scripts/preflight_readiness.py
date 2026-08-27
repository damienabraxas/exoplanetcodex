#!/usr/bin/env python3
"""RYA-1069: preflight data-readiness orchestrator.

Per-system conductor. Answers, for every holding of a system, whether it is
measurement-ready, and if not, WHICH gate is missing. Runs no science, selects
no product, normalises nothing, corrects no tellurics -- it CONDUCTS the existing
audit/policy stages and derives one measurement_ready verdict per holding x band.

Sibling, NOT duplicate, of `scripts/preflight_check.py` (RYA-905, per-ELEMENT-RUN
advisory, WARN-not-block, runs immediately before a measurement). This one is
per-SYSTEM and GATING, runs once after download, and its verdict is what the
measurement orchestrator (RYA-767 `pipeline/run_descriptor.py`) keys on.
Origin: the RYA-1064 STOP -- `docs/science/rya1064_orchestration_stop.md`.

WHAT IT CONDUCTS (nothing here is new science; every gate binds to a standing stage)
-----------------------------------------------------------------------------------
    evidence       `holdings_manifest_registry.csv::evidence_state`
    product        DERIVED from the registry's own manifest -- see `product_selection`
    normalization  ADAPTER 1: `rya1030_backfill_normalisation.determine()`, which probes
                   THROUGH THE HARNESS READER and classifies the FLUX. Never a label.
    telluric       `pipeline.telluric_policy` -- BOTH axes (RYA-786 instrument basis,
                   RYA-806 per-holding applied), combined here by `telluric_satisfied`.
    reader         ADAPTER 2: `measure_band_ew.holdings_for()` + `HoldingSpec.covers()`
    line_pool      `canonical_gf.csv` filtered by `gf_empirical.GRADED_TIERS`, then by
                   what the reader can actually serve.

THE RULE THAT DECIDES EVERY AMBIGUOUS CASE (RYA-161 / loud-fail / no-silent-fallback)
-------------------------------------------------------------------------------------
An UNKNOWN on any gate makes the holding NOT-READY **with that gate named**. `unknown`
is never defaulted to ready. A scientifically correct NO-GO carrying the exact blocking
gate is the SUCCESS output of this script, exactly as RYA-1064 stopped -- it never
invents a coadd, a continuum, or a telluric state to manufacture a GO.

No band edge, path, or threshold is written down here. Band edges come from
`config.synth_bands`, instrument reach from `instrument_catalog.csv`, the graded tier set
from `pipeline.gf_empirical`, the telluric enumeration from `pipeline.telluric_policy`,
and the normalisation probe geometry from the RYA-1030 backfill that owns it.

WHAT IT DOES NOT DO. It does not write back into `holdings_manifest_registry.csv` (that
file is generated / hand-curated; the rollup wiring is a follow-on once this verdict is
trusted), and it does not select a product, normalise, or run a telluric correction. It
DETECTS and REPORTS which conditioning stage is owed.
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from dataclasses import asdict, dataclass, fields
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from config import synth_bands                                        # noqa: E402
from pipeline import telluric_policy                                  # noqa: E402
from pipeline.gf_empirical import GRADED_TIERS                        # noqa: E402
from pipeline.normalization_intake import NORMALISED, UNKNOWN         # noqa: E402

REGISTRY = REPO / "data" / "catalog" / "holdings_manifest_registry.csv"
INSTRUMENTS = REPO / "data" / "catalog" / "instrument_catalog.csv"
CANONICAL_GF = REPO / "data" / "linelists" / "canonical_gf.csv"
OUT_DIR = REPO / "data" / "audit" / "readiness"

#: The `evidence_state` the registry writes when the identity/provenance evidence has
#: been checked. Read as a constant rather than retyped at each comparison.
EVIDENCE_VERIFIED = "verified"

#: The red end of the project's ENUMERATED telluric contamination -- derived from
#: `telluric_policy.TELLURIC_BANDS`, never written down again. Past it the per-line
#: clean-line selection basis has nothing to select against, so an absence of enumerated
#: bands there is not evidence of a clean band (see `telluric_satisfied`).
TELLURIC_DOMAIN_MAX_A = max(hi for _lo, hi, _name in telluric_policy.TELLURIC_BANDS)

#: Instrument telluric bases that positively declare "no correction STAGE is owed".
#: `unspecified` is included on RYA-786's ratified reading -- an instrument registered
#: `telluric_required=no` whose basis is not stated is treated as `line_selection`,
#: because excluding a good line costs coverage while measuring a telluric one costs the
#: number. `mode_dependent` is deliberately ABSENT: it is a declared ambiguity, not a
#: declared basis, and resolving it needs the mode. See the ledger note in this ticket.
NON_CORRECTION_BASES = ("line_selection", "corrected", "unspecified")


@dataclass
class ReadinessRow:
    system_id: str
    holding_id: str
    instrument_id: str
    band: str
    coverage_A: str
    evidence_state: str
    product_selected: bool
    normalization_state: str
    telluric_basis: str
    telluric_applied: str
    telluric_satisfied: bool
    reader_wired: bool
    line_pool_reachable: bool
    measurement_ready: str
    blocking_gate: str
    source_issue_ids: str


# ── ADAPTER 2: the measurement harness's HoldingSpec table ───────────────────
#
# `scripts/measure_band_ew.py` resolves the Kitt Peak atlas AT IMPORT and `SystemExit`s
# when it is absent -- RYA-1064 friction #3, reproduced here: a pure Alpha Cen preflight
# stops on unrelated SOLAR bytes before it can report the Alpha Cen refusal. This
# conductor must NEVER inherit that: an unimportable harness is reported as a NAMED
# reader refusal, distinct from "this holding has no HoldingSpec", because collapsing
# those two would let a wiring failure wear a data verdict's clothes (RYA-833).

_HARNESS: dict = {}


def harness():
    """The measurement harness module, or None. Loaded ONCE, never fatally."""
    if "mod" not in _HARNESS:
        import importlib.util
        spec = importlib.util.spec_from_file_location(
            "_mbe_rya1069", REPO / "scripts" / "measure_band_ew.py")
        mod = importlib.util.module_from_spec(spec)
        sys.modules["_mbe_rya1069"] = mod
        try:
            spec.loader.exec_module(mod)
            _HARNESS["mod"], _HARNESS["why"] = mod, ""
        except BaseException as exc:            # SystemExit included, deliberately
            sys.modules.pop("_mbe_rya1069", None)
            _HARNESS["mod"] = None
            _HARNESS["why"] = f"{type(exc).__name__}: {str(exc).splitlines()[0]}"
            print(f"⚠️  ADAPTER 2 UNAVAILABLE -- scripts/measure_band_ew.py could not be "
                  f"imported: {_HARNESS['why']}\n"
                  f"    Every `reader` gate below reports harness-unavailable, which is "
                  f"NOT the same finding as 'this holding has no HoldingSpec'. Resolve "
                  f"the harness (CODEX_KP_ATLAS, or stage the atlas) and re-run before "
                  f"reading any reader verdict as a fact about the DATA.",
                  file=sys.stderr)
    return _HARNESS["mod"]


def holding_spec(holding_id: str):
    """The `HoldingSpec` wired for this holding, or None. ADAPTER 2's lookup."""
    h = harness()
    if h is None:
        return None
    for specs in h._INSTRUMENT_HOLDINGS.values():
        for spec in specs:
            if spec.holding_id == holding_id:
                return spec
    return None


# ── ADAPTER 1: the codex-data-audit stages ───────────────────────────────────

def _backfill():
    """`scripts/rya1030_backfill_normalisation.py`, loaded by PATH.

    Imported by file rather than by name so that `scripts/` never joins `sys.path` --
    that directory holds ~90 modules whose names would shadow real packages on the next
    import anyone adds.
    """
    if "backfill" not in _HARNESS:
        import importlib.util
        spec = importlib.util.spec_from_file_location(
            "_rya1030_rya1069", REPO / "scripts" / "rya1030_backfill_normalisation.py")
        mod = importlib.util.module_from_spec(spec)
        sys.modules["_rya1030_rya1069"] = mod
        spec.loader.exec_module(mod)
        _HARNESS["backfill"] = mod
    return _HARNESS["backfill"]


@dataclass
class AuditResult:
    """What the intake stages said about the HOLDING (not about one band)."""
    product_selected: bool
    product_tag: str
    product_evidence: str
    normalization_state: str
    normalization_evidence: str


def _tokens(value) -> set[str]:
    return set(re.split(r"[^A-Za-z0-9_]+", str(value).lower()))


def _squash(value) -> str:
    return re.sub(r"[^a-z0-9]", "", str(value).lower())


#: A cell that NAMES A FILE: it carries a path separator, or ends in a data-file
#: extension. Content test, not a column-name guess -- the four manifests in the registry
#: call the same column `file`, `filename`, `filepath` and `product`.
_FILEISH = re.compile(r"[/\\]|\.(fits|fit|dat|csv|tsv|txt|gz|ecsv|json)(\.[a-z0-9]+)?$",
                      re.IGNORECASE)


def _read_manifest_rows(path: Path) -> list[dict]:
    """Rows of a tabular manifest, with leading `#` comment lines skipped.

    `solar_reference_holdings_rya708.csv` opens with four commented lines. Feeding those
    to `DictReader` makes the comment the header and the real header a data row, which
    silently turns a one-row-per-instrument catalogue into an unreadable blob.
    """
    with path.open(encoding="utf-8", errors="replace", newline="") as fh:
        lines = [ln for ln in fh if not ln.lstrip().startswith("#")]
    return list(csv.DictReader(lines))


_ALIAS: dict = {}


def _instrument_alias(cell: str) -> str | None:
    """The catalog `instrument_id` a manifest cell names, or None.

    Manifests write the FACILITY's name for an instrument -- the alpha Cen inventory says
    `CRIRES` for products this repo registers as `crires_plus`. An exact squashed match is
    tried first; a prefix match is accepted only when EXACTLY ONE registered instrument
    extends it, so `HARPS` can never be silently attributed to `harps_n`.
    """
    if not _ALIAS:
        _ALIAS["ids"] = [r["instrument_id"] for r in
                         csv.DictReader(INSTRUMENTS.open(encoding="utf-8"))]
    key = _squash(cell)
    if not key:
        return None
    exact = [i for i in _ALIAS["ids"] if _squash(i) == key]
    if exact:
        return exact[0]
    starts = [i for i in _ALIAS["ids"] if _squash(i).startswith(key)]
    return starts[0] if len(starts) == 1 else None


def _attribute(rows: list[dict], holding_id: str, system_id: str,
               instrument_id: str) -> tuple[list[dict], str]:
    """Which manifest rows belong to this holding, and how we decided.

    Most narrowly first. A manifest shared by six holdings (the RYA-479 alpha Cen
    inventory is one file for HARPS, ESPRESSO, NIRPS and CRIRES+ on both components)
    must not be counted whole against each of them -- that would report the same
    291 products six times and lose the fact that the HARPS arm alone holds 88.
    """
    by_hold = [r for r in rows if any(holding_id in _tokens(v) for v in r.values() if v)]
    by_sys = [r for r in rows if any(system_id in _tokens(v) for v in r.values() if v)]
    by_inst = [r for r in rows
               if any(_instrument_alias(v) == instrument_id for v in r.values() if v)]
    both = [r for r in by_sys if r in by_inst]
    for rows_here, how in (
            (by_hold, f"row(s) naming holding {holding_id}"),
            (both, f"row(s) naming system {system_id} AND instrument {instrument_id}"),
            (by_inst, f"row(s) naming instrument {instrument_id}"),
            (by_sys, f"row(s) naming system {system_id}")):
        if rows_here:
            return rows_here, how
    return rows, "NO row attributes itself to this holding, so ALL rows"


def product_selection(holding_id: str, system_id: str, instrument_id: str,
                      manifest_path: str) -> tuple[bool, str, str]:
    """Does the registry declare ONE science product for this holding, or an inventory?

    🔴 THIS IS THE RYA-1064 GATE, and it is deliberately NOT the reader gate. 1064 stopped
    because `alpha_cen_a_harps` resolves to 205 loose source products with "no declared
    run product/coadd" -- so no checksum, coverage, resolving power or continuum contract
    could be pinned to anything. That is a fact about the REGISTRY, answerable before any
    reader exists, and keeping the two gates independent is what lets the verdict say
    which of the two follow-ups (condition/register vs. wire) is owed.

    Derived from the manifest's own CONTENT, in this order:

      * the manifest does not resolve            -> NOT selected (nothing to check)
      * JSON declaring `holding_id` == this one  -> selected (a per-holding declaration)
      * a tabular manifest -> attribute rows to this holding (by holding_id token, else
        by instrument_id), then count the DISTINCT files those rows name:
            1 file      -> selected
            >1 files    -> NOT selected: an N-file source inventory
            no file column at all -> the manifest_path IS the product (the RYA-794 and
                RYA-1054 conditioned CSVs are spectra, not inventories) -> selected
      * anything else (prose, a checksum list)   -> NOT selected: it declares no product
    """
    if not manifest_path:
        return False, "no-manifest", "the registry names no manifest_path for this holding"
    path = REPO / manifest_path
    if not path.exists():
        return False, "manifest-missing", (
            f"manifest_path {manifest_path} does not resolve, so nothing declares which "
            f"product this holding is")

    if path.suffix.lower() == ".json":
        try:
            doc = json.loads(path.read_text(encoding="utf-8", errors="replace"))
        except ValueError as exc:
            return False, "manifest-unreadable", f"{manifest_path} is not readable JSON ({exc})"
        if isinstance(doc, dict) and str(doc.get("holding_id", "")) == holding_id:
            return True, "declared", (
                f"{manifest_path} is a per-holding intake declaration naming "
                f"holding_id={holding_id}")
        return False, "not-declared", (
            f"{manifest_path} is JSON but declares no holding_id={holding_id}, so it does "
            f"not select a product for this holding")

    if path.suffix.lower() not in (".csv", ".tsv"):
        return False, "not-a-manifest", (
            f"{manifest_path} is not a tabular manifest and declares no single science "
            f"product -- it cannot pin a checksum, coverage or continuum contract to one "
            f"file (RYA-1064)")

    rows = _read_manifest_rows(path)
    if not rows:
        return False, "manifest-empty", f"{manifest_path} holds no rows"

    cols = [c for c in rows[0] if c]
    file_cols = [c for c in cols
                 if sum(bool(_FILEISH.search(str(r.get(c) or ""))) for r in rows)
                 >= max(1, len(rows) // 2)]

    rows_here, how = _attribute(rows, holding_id, system_id, instrument_id)

    if not file_cols:
        return True, "is-the-product", (
            f"{manifest_path} names no product files -- the manifest_path IS the single "
            f"conditioned product ({len(rows)} data rows, columns {cols})")

    products = {tuple(str(r.get(c) or "") for c in file_cols) for r in rows_here}
    if len(products) == 1:
        return True, "selected", (
            f"{manifest_path} declares exactly ONE product for this holding, from "
            f"{len(rows_here)} {how}: {'|'.join(next(iter(products)))}")
    return False, f"n-file-inventory({len(products)})", (
        f"{manifest_path} is an N-FILE SOURCE INVENTORY: {len(products)} distinct "
        f"products across {len(rows_here)} {how}. No single approved run product/coadd is "
        f"selected, so its checksum, coverage, resolving power and continuum contract are "
        f"all unpinned (RYA-1064). Choosing or coadding one here would be inventing "
        f"science -- that is follow-up #1's job, not this conductor's.")


def run_data_audit(h: dict) -> AuditResult:
    """ADAPTER 1. The codex-data-audit stages, CALLED, never re-implemented.

    `normalization_state` comes from `scripts/rya1030_backfill_normalisation.determine`,
    which is the ratified flux-measured entry point: it probes at four fractions of the
    holding's declared span THROUGH `measure_band_ew.load_window_ex` -- the same call the
    measurement harness makes, on the same holding -- and classifies what comes back with
    `pipeline.normalization_intake.detect`. Reading the registry's own
    `normalization_state` column instead would be reading a LABEL, which the skill
    forbids and which is exactly how KP2005 ran mis-routed for months (RYA-929/933/1026).
    """
    product_ok, product_tag, product_why = product_selection(
        h["holding_id"], h["system_id"], h["instrument_id"], h.get("manifest_path", ""))

    mod = harness()
    if mod is None:
        return AuditResult(
            product_ok, product_tag, product_why, UNKNOWN,
            f"NOT MEASURED: the harness reader is unavailable ({_HARNESS['why']}), so "
            f"the flux was never read. `unknown` here is a fact about THIS RUN, not "
            f"about the product (RYA-833).")

    value, evidence = _backfill().determine(mod, h["holding_id"], h["instrument_id"])
    return AuditResult(product_ok, product_tag, product_why, value, evidence)


# ── bands ────────────────────────────────────────────────────────────────────

def _instrument_span_A(instrument_id: str) -> tuple[float, float]:
    """The instrument's declared reach, in Angstrom, from the catalog. Loud on unknown."""
    for row in csv.DictReader(INSTRUMENTS.open(encoding="utf-8")):
        if row["instrument_id"] == instrument_id:
            return (float(row["wavelength_min_nm"]) * 10.0,
                    float(row["wavelength_max_nm"]) * 10.0)
    raise SystemExit(
        f"instrument {instrument_id!r} is not in {INSTRUMENTS.name}; its wavelength reach "
        f"is unknown and must not be assumed. Register it first (RYA-786).")


def bands_for(instrument_id: str) -> list[tuple[str, float, float]]:
    """(band_name, lo_A, hi_A) -- the science bands this instrument actually reaches.

    The band EDGES are `config/synth_bands.yaml` (RYA-967's single source); the REACH is
    `instrument_catalog.csv`. Each row is the OVERLAP, so `coverage_A` is what the
    instrument can serve of that band rather than the band's nominal extent. Nothing here
    is written down: an edge typed into this file would become the second home RYA-967
    just removed.

    The overlap is computed on WAVELENGTH, not on `bands_supported`, whose vocabulary
    (`red_optical`, `NUV`, `Y`/`J`/`H`) is the instrument axis and does not map onto the
    four science bands one-for-one.
    """
    lo_i, hi_i = _instrument_span_A(instrument_id)
    out = []
    for name, band in sorted(synth_bands.SYNTH_BANDS.items(), key=lambda kv: kv[1].lo_A):
        lo, hi = max(band.lo_A, lo_i), min(band.hi_A, hi_i)
        if hi > lo:
            out.append((name, lo, hi))
    return out


# ── the derived telluric verdict ─────────────────────────────────────────────

def _reader_reach(spec, lo: float, hi: float, pad: float) -> tuple[float, bool]:
    """(probe centre, does the reader reach INTO this band?) -- ADAPTER 2's question.

    🔴 ASKED AT THE READER'S OWN OVERLAP, NOT AT THE BAND CENTRE. `HoldingSpec.covers`
    demands TOTAL coverage of a window, which is right (half a window inside a product is
    a truncated window). Asking it at the band's midpoint therefore answers a different
    question -- "does this holding span the whole band" -- and the RYA-794 Y arm
    (10280-10680 A) would report `reader: out-of-span` for the NIR band it demonstrably
    serves, while `line_pool_reachable` simultaneously found lines it CAN serve. Two
    gates disagreeing about the same reader is a defect, not a nuance.

    A `span_A` of None means the reader inventories its own coverage (the Kitt Peak
    segment list, the CRIRES+ IDP comb), so `covers()` answers True and the band centre
    is the honest probe point.
    """
    if spec is None:
        return 0.5 * (lo + hi), False
    if spec.span_A is None:
        centre = 0.5 * (lo + hi)
        return centre, spec.covers(centre, pad)
    o_lo, o_hi = max(spec.span_A[0], lo), min(spec.span_A[1], hi)
    if o_hi <= o_lo:
        return 0.5 * (lo + hi), False
    centre = 0.5 * (o_lo + o_hi)
    return centre, spec.covers(centre, pad)


def telluric_satisfied(instrument_id: str, band: str, lo: float, hi: float,
                       basis: str, applied: str,
                       evidence_state: str) -> tuple[bool, str, str]:
    """Does this band's telluric state permit measurement? BOTH axes, in order.

    🔴 THE IR RULE, no exceptions: a band that reaches past the reddest ENUMERATED
    telluric complex is NOT satisfied unless the correction was actually applied to this
    product. Per-line clean-line selection is the ratified basis for the optical
    (RYA-460/786) precisely because the O2/H2O complexes there are enumerated and a clean
    alternative line exists. Past `TELLURIC_DOMAIN_MAX_A` the enumeration stops, so an
    absence of listed bands is an absence of LOOKING, not evidence of a clean band, and
    treating it as satisfied would measure terrestrial absorption as stellar flux.

    `unknown` never reaches a True by any route. Neither does `mode_dependent`: this
    conductor refused it from the start, and RYA-1072 made `telluric_policy` refuse it
    too -- `requires_correction` now RAISES `TelluricStateUnknown` on an unresolved
    requirement rather than reading it as "no", so the two now agree at the source
    instead of only here.
    """
    if basis == "not_applicable":
        return True, "", (f"{instrument_id} is telluric_basis=not_applicable -- no "
                          f"terrestrial atmosphere is in this light path (RYA-806)")
    if applied == "applied" and evidence_state == EVIDENCE_VERIFIED:
        return True, "", (f"telluric_applied=applied on {EVIDENCE_VERIFIED} evidence; the "
                          f"correction is in this product, RYA-424 is SKIPPED")
    if hi > TELLURIC_DOMAIN_MAX_A:
        return False, "ir-band-uncorrected", (
            f"THE IR RULE: this band reaches {hi:.0f} A, past the reddest ENUMERATED "
            f"telluric complex ({TELLURIC_DOMAIN_MAX_A:.0f} A), and telluric_applied="
            f"{applied}. Beyond the enumeration, per-line clean-line selection has "
            f"nothing to select against -- an absence of listed bands there is an absence "
            f"of LOOKING, not a clean band. No exceptions (RYA-1069 spec).")
    if applied == "applied":
        return False, "applied-unverified", (
            f"telluric_applied=applied but evidence_state={evidence_state!r}, not "
            f"{EVIDENCE_VERIFIED!r} -- an unverified correction is not a correction")
    requirement = telluric_policy.correction_requirement(instrument_id)
    if requirement == telluric_policy.AMBIGUOUS:
        return False, f"basis-{basis}", (
            f"{instrument_id} is telluric_required=mode_dependent/unresolved: the "
            f"catalog does not say whether a correction stage is owed, and RYA-1072 "
            f"makes telluric_policy refuse rather than read that as 'no'. Resolve the "
            f"mode (instrument_modes.csv) and record the value.")
    if requirement == telluric_policy.REQUIRED:
        return False, f"needs-correction({applied})", (
            f"{instrument_id} is registered telluric_required=yes and telluric_applied="
            f"{applied}; route it through the RYA-424 telluric stage first")
    if basis in NON_CORRECTION_BASES:
        return True, "", (
            f"{instrument_id} is telluric_required=no on basis {basis!r}; tellurics are "
            f"handled by per-line clean-line selection over the "
            f"{len(telluric_policy.TELLURIC_BANDS)} enumerated bands (RYA-460/786), which "
            f"is defined for this band because it ends at {hi:.0f} A, inside the "
            f"enumerated domain. No telluric_corrected declaration is made.")
    return False, f"basis-{basis}", (
        f"telluric_basis={basis!r} is a DECLARED AMBIGUITY, not a declared basis, so "
        f"per-line clean-line selection cannot be claimed for this band. Resolve it on "
        f"the instrument/mode axis; `unknown` is never defaulted to ready (RYA-161).")


# ── the graded line pool ─────────────────────────────────────────────────────

_POOL: dict = {}


def graded_lines(element: str | None, ion: str | None) -> list[float]:
    """Air wavelengths of every GRADED line in `canonical_gf.csv`.

    `GRADED_TIERS` is imported from `pipeline.gf_empirical` (RYA-945's tier vocabulary),
    never re-enumerated here -- a second copy of that set is how a pool definition drifts.
    """
    key = (element, ion)
    if key not in _POOL:
        import pandas as pd
        df = pd.read_csv(CANONICAL_GF, low_memory=False)
        df = df[df.gf_tier.astype(str).isin(GRADED_TIERS)]
        if element:
            df = df[df.species.astype(str).str.split().str[0] == element]
        if ion:
            df = df[df.species.astype(str).str.split().str[-1] == ion]
        _POOL[key] = sorted(
            float(w) for w in pd.to_numeric(df.wavelength_air_A, errors="coerce").dropna())
    return _POOL[key]


def pool_reachable(spec, lo: float, hi: float, pad: float,
                   element: str | None,
                   ion: str | None) -> tuple[bool, int, int]:
    """Is at least one graded line in this band SERVABLE by the wired reader?

    Two conditions, and they are different questions: the line has to be in the band, and
    the reader has to cover a window around it. With no HoldingSpec the answer is zero
    reachable lines -- not "the pool is empty", which would blame the linelist for a
    wiring gap (RYA-833).
    """
    in_band = [w for w in graded_lines(element, ion) if lo <= w <= hi]
    if spec is None:
        return False, len(in_band), 0
    servable = [w for w in in_band if spec.covers(w, pad)]
    return bool(servable), len(in_band), len(servable)


# ── the conductor ────────────────────────────────────────────────────────────

def holdings_for_system(system_id: str, instrument: str | None) -> list[dict]:
    rows = list(csv.DictReader(REGISTRY.open(encoding="utf-8")))
    hits = [r for r in rows if r["system_id"] == system_id
            and (instrument is None or r["instrument_id"] == instrument)]
    if not hits:
        known = sorted({r["system_id"] for r in rows})
        raise SystemExit(
            f"NO HOLDINGS for system={system_id!r} instrument={instrument!r} in "
            f"{REGISTRY.name} -- refuse, do not guess a holding. "
            f"Registered systems: {', '.join(known)}")
    return hits


def assess(h: dict, element: str | None, ion: str | None,
           evidence: dict) -> list[ReadinessRow]:
    inst = h["instrument_id"]
    audit = run_data_audit(h)
    spec = holding_spec(h["holding_id"])
    evidence_state = (h.get("evidence_state") or "unknown").strip()

    basis = telluric_policy.basis(inst)
    try:
        applied = telluric_policy.applied_state(h["holding_id"])
    except telluric_policy.TelluricStateUnknown:
        applied = "unknown"

    bands = bands_for(inst)
    reddest = max(b.hi_A for b in synth_bands.SYNTH_BANDS.values())
    _lo_i, hi_i = _instrument_span_A(inst)
    if hi_i > reddest:
        pool = graded_lines(element, ion)
        beyond = sum(1 for w in pool if reddest < w <= hi_i)
        evidence[(h["holding_id"], "DECLARED GAP")] = [
            f"reach:{reddest:.0f}-{hi_i:.0f} A -- {inst} reaches {hi_i:.0f} A but "
            f"config/synth_bands.yaml defines no band past {reddest:.0f} A, so that "
            f"stretch is NOT ASSESSED here. Named rather than dropped: silence would "
            f"read as 'no coverage owed' (RYA-833). {beyond} graded line(s) sit in it "
            f"(canonical_gf reaches {max(pool):.0f} A), so the gap is real work, not a "
            f"formality."]

    out: list[ReadinessRow] = []
    for band, lo, hi in bands:
        pad = synth_bands.SYNTH_BANDS[band].half_width_A
        centre, reader_wired = _reader_reach(spec, lo, hi, pad)
        tsat, tag, twhy = telluric_satisfied(
            inst, band, lo, hi, basis, applied, evidence_state)
        reachable, n_in_band, n_servable = pool_reachable(
            spec, lo, hi, pad, element, ion)

        if harness() is None:
            reader_why, reader_tag = (
                f"the harness could not be imported ({_HARNESS['why']}), so nothing was "
                f"asked about a reader for this holding", "harness-unavailable")
        elif spec is None:
            reader_why, reader_tag = (
                f"no HoldingSpec for {h['holding_id']} in "
                f"measure_band_ew._INSTRUMENT_HOLDINGS -- the holding is registered but "
                f"the measurement harness cannot address it (RYA-1030/1064)", "unwired")
        else:
            reader_why, reader_tag = (
                f"{spec.holding_id} declares span_A={spec.span_A}, which does not reach "
                f"into {lo:.0f}-{hi:.0f} A far enough to hold one {pad} A window",
                "out-of-band")

        if spec is None:
            pool_why, pool_tag = (
                f"{n_in_band} graded line(s) lie in {lo:.0f}-{hi:.0f} A but no reader "
                f"serves them -- a wiring gap, not an empty pool (RYA-833)",
                f"unreachable({n_in_band}-in-band)")
        else:
            pool_why, pool_tag = (
                f"{n_in_band} graded line(s) in {lo:.0f}-{hi:.0f} A, {n_servable} of them "
                f"servable by {spec.holding_id}", f"0-servable-of-{n_in_band}")

        # (passed, SHORT tag for `blocking_gate`, FULL evidence for the operator)
        gates = {
            "evidence": (evidence_state == EVIDENCE_VERIFIED, evidence_state,
                         f"registry evidence_state={evidence_state!r}, not "
                         f"{EVIDENCE_VERIFIED!r}"),
            "product": (audit.product_selected, audit.product_tag,
                        audit.product_evidence),
            "normalization": (audit.normalization_state == NORMALISED,
                              audit.normalization_state, audit.normalization_evidence),
            "telluric": (tsat, tag, twhy),
            "reader": (reader_wired, reader_tag, reader_why),
            "line_pool": (reachable, pool_tag, pool_why),
        }
        ready = all(ok for ok, _t, _w in gates.values())
        blocking = ";".join(f"{k}:{t}" for k, (ok, t, _w) in gates.items() if not ok)

        # ── RYA-1079 §5: an uncorrected but CORRECTABLE band is WORK, not a wall ──
        # The three-state per-line policy (pipeline.telluric_observability) says a
        # RECOVERABLE telluric is corrected and then measured; only a SATURATED core is
        # genuinely lost. The conductor has to carry that distinction or an operator
        # reads "NO-GO" and abandons a band that needs one conditioning run. The gate
        # tag is what separates them: `needs-correction(...)` names an available stage,
        # while `ir-band-uncorrected`, `basis-*` and `applied-unverified` do not.
        # ACTIONABLE means a correction stage can be NAMED for this cell. That is true
        # whenever the requirement is RESOLVED and the correction simply has not been run
        # -- `needs-correction(...)` (the instrument is registered telluric_required=yes)
        # and `ir-band-uncorrected` (the band reaches past the enumerated complexes and
        # carries no applied correction) are both that case. Neither is terminal: RYA-963
        # molecfit-corrected exactly the first, RYA-940 exactly the second.
        #
        # It is NOT actionable when the requirement or basis is UNRESOLVED (`basis-*`),
        # because you cannot route what you cannot classify -- resolve the mode first
        # (RYA-1072/1078).
        #
        # ⚠️ SATURATION IS NOT DECIDABLE AT BAND LEVEL. Whether a given line is genuinely
        # lost is a PER-LINE measurement (pipeline.telluric_observability), and a band
        # holds clean, recoverable and saturated lines at once. So the band verdict here
        # is NEEDS-CORRECTION and the per-line census decides what survives it; calling a
        # whole band terminal on a telluric ground would exclude the recoverable majority
        # to spare the saturated few, which is the data loss RYA-1079 exists to stop.
        _ACTIONABLE_TELLURIC = ("needs-correction", "ir-band-uncorrected")
        only_telluric = {k for k, (ok, _t, _w) in gates.items() if not ok} == {"telluric"}
        actionable = gates["telluric"][1].startswith(_ACTIONABLE_TELLURIC)
        if only_telluric and actionable:
            verdict = "NEEDS-CORRECTION"
        elif ready:
            verdict = "GO"
        else:
            verdict = "NO-GO"
        evidence[(h["holding_id"], band)] = [
            f"{k}:{t} -- {w}" for k, (ok, t, w) in gates.items() if not ok]
        if ready and (audit.normalization_state == UNKNOWN
                      or (applied == "unknown" and basis != "not_applicable")):
            raise SystemExit(
                f"CRITICAL: {h['holding_id']} {band} was about to be reported GO while "
                f"normalization_state={audit.normalization_state} "
                f"telluric_applied={applied}. An `unknown` reached a READY verdict, which "
                f"is the silent-ready fallback this conductor exists to prevent "
                f"(RYA-161). Refusing rather than emitting the row.")

        out.append(ReadinessRow(
            system_id=h["system_id"], holding_id=h["holding_id"], instrument_id=inst,
            band=band, coverage_A=f"{lo:.1f}-{hi:.1f}",
            evidence_state=evidence_state,
            product_selected=audit.product_selected,
            normalization_state=audit.normalization_state,
            telluric_basis=basis, telluric_applied=applied,
            telluric_satisfied=tsat, reader_wired=reader_wired,
            line_pool_reachable=reachable,
            measurement_ready=verdict,
            blocking_gate=blocking,
            source_issue_ids=h.get("source_issue_ids", "")))

    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description="Per-system data-readiness conductor (RYA-1069, RYA-1064 follow-up).")
    ap.add_argument("--system", required=True, help="system_id as the registry writes it")
    ap.add_argument("--instrument", default=None, help="limit to one instrument_id")
    ap.add_argument("--element", default=None,
                    help="limit the graded line pool to one element (e.g. Fe)")
    ap.add_argument("--ion", default=None,
                    help="limit the graded line pool to one ion (e.g. I)")
    a = ap.parse_args(argv)

    evidence: dict = {}
    rows = [r for h in holdings_for_system(a.system, a.instrument)
            for r in assess(h, a.element, a.ion, evidence)]

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    dest = OUT_DIR / f"{a.system}_readiness.csv"
    with dest.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=[f.name for f in fields(ReadinessRow)])
        w.writeheader()
        for r in rows:
            w.writerow(asdict(r))

    ready = sum(r.measurement_ready == "GO" for r in rows)
    owed = sum(r.measurement_ready == "NEEDS-CORRECTION" for r in rows)
    print(f"{a.system}: {ready}/{len(rows)} holding x band cells READY"
          + (f", {owed} NEEDS-CORRECTION (actionable, not terminal)" if owed else "")
          + f" -> "
          f"{dest.relative_to(REPO)}")
    for r in rows:
        mark = ("READY" if r.measurement_ready == "GO"
                else f"{r.measurement_ready} [{r.blocking_gate}]")
        print(f"  {r.holding_id:<34} {r.band:<12} {r.coverage_A:<16} {mark}")

    if evidence:
        print("\nWHY -- the evidence behind every blocked gate. `blocking_gate` in the "
              "CSV is the\nshort form; this is what each verdict was read off.")
        for (hid, band), lines in evidence.items():
            print(f"\n  {hid} | {band}")
            for line in lines:
                print(f"    - {line}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
