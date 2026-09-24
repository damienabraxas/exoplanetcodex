"""Expands (star, element) into RunDescriptors and executes the runnable ones; decides no science.

RYA-1222 (RYA-767 Seam 2's executor). Owns: the CELL AXES of a run matrix, the
order cells are attempted in, whether a cell's work is already current, and one
report that names every cell. Does NOT own: whether a run is possible (that is
`run_descriptor.resolve`), whether the data is ready (that is RYA-1069's
conductor), or any measurement (that is the stage scripts). It runs nothing it
was not handed by the resolver and it never decides a number.

WHY THIS EXISTS
---------------
`run_descriptor.resolve()` answers ONE descriptor. Producing an element meant a
person writing the descriptors out, remembering which holdings exist, which
bands each instrument reaches, which decks are live, and which of those cells
were already built last week. That is the loop this closes: the judgment is in
the ledgers, and the executor makes no decisions at all.

THE FOUR RULES THIS FILE ENFORCES
---------------------------------
1. NO AXIS IS TYPED HERE. Holdings and instruments come from
   `holdings_manifest_registry.csv`, band edges from `config/synth_bands.yaml`
   via RYA-1069's `bands_for` (so the matrix and the preflight agree cell for
   cell), instrument reach from `instrument_catalog.csv`, the engine axis from
   `model_registry.csv`, the ion axis from the graded pool in
   `canonical_gf.csv`, and the element vocabulary from `elements_master.json`.
   A literal list of any of them would be a second declaration of a fact that
   already has a home, which is how RYA-845's double-count survived.

2. EVERY CELL GETS A TERMINAL STATUS. A cell that cannot run is BLOCKED or
   NOT_READY *with the reason*, never dropped. RYA-1187: a silent gap and an
   honest empty read identically to every downstream consumer, and only one of
   them is a finding. The sole cells never produced are the ones the wired
   `HoldingSpec` positively declares out of span -- and that is a declaration,
   not an absence.

3. ONE FAILURE IS ONE CELL. A stage error is captured, the cell goes FAILED,
   and the matrix continues. An aborted matrix reports nothing about the cells
   it never reached, which is the same silent gap as rule 2 wearing a crash.

4. WORK ALREADY CURRENT IS NOT REDONE. A cell is SKIP when the published feed
   still holds this cell's products AND the inputs-hash recorded for the cell
   equals the one this run computes. The hash is over the descriptor, the
   ordered stage scripts' own bytes, and the named input artifacts -- so a code
   edit or a re-staged spectrum invalidates it and a clock tick does not. The
   hash lives in a sidecar ledger, NOT on the product; see `LEDGER` for the
   measurement that forced that and why it changes nothing about the guarantee.

WHAT A "CELL" IS
----------------
(band x instrument/holding x ion x engine-deck) for one (star, element). The
holding is on the axis and the instrument is not enough: RYA-933/934 collided
`solar_harps` with `solar_harps_molecfit_corrected`, two products differing by
exactly whether tellurics were removed.
"""
from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import os
import subprocess
import sys
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

from pipeline.run_descriptor import RunDescriptor, resolve  # noqa: E402
from pipeline import band_policy  # noqa: E402

MODEL_REGISTRY = ROOT / "data" / "catalog" / "model_registry.csv"
HOLDINGS_REGISTRY = ROOT / "data" / "catalog" / "holdings_manifest_registry.csv"
CANONICAL_GF = ROOT / "data" / "linelists" / "canonical_gf.csv"
ELEMENTS_MASTER = ROOT / "data" / "config" / "elements_master.json"
PRODUCTS = ROOT / "data" / "products"
REPORT_DIR = ROOT / "data" / "results" / "orchestrator"
#: What lands in REPORT_DIR and why the hash is not on the product:
#: docs/orchestrator_run_reports_rya1222.md

# ── the terminal statuses ────────────────────────────────────────────────────
#
# SKIP and DONE are the only ones that mean "this cell's work exists". The rest
# are all honest empties and each names a DIFFERENT owner: BLOCKED is the
# resolver's (this run cannot be expressed), NOT_READY is the data's (RYA-1069
# names the conditioning stage owed), FAILED is ours (a stage ran and broke).
# Collapsing any two of them would put a wiring fault in a data verdict's
# clothes -- the RYA-833 shape.
DONE, SKIP, FAILED, BLOCKED, NOT_READY = "DONE", "SKIP", "FAILED", "BLOCKED", "NOT_READY"
#: --dry-run only. A cell that passed every gate and WOULD have been executed.
#: It is its own status because the four above are all claims about the world and
#: this one is a claim about an intention: calling it DONE would report work that
#: did not happen, and calling it NOT_READY would blame the data for a run we
#: chose not to start.
WOULD_RUN = "WOULD_RUN"
#: 🔴 The cell RAN, every stage exited 0, and NOTHING reached the published feed.
#:
#: `derive_band_products.py` writes an artifact under `data/results/band_products/`
#: and then tells the operator to publish it -- publication is `publish_product.py`,
#: which is human-gated on purpose: RYA-1034 requires a stated `--reason` to move a
#: value and RYA-772 calls copying a foreign artifact without provenance laundering.
#: An orchestrator that auto-published would be doing exactly that.
#:
#: So this state is real and it must not be called DONE. Currency is defined against
#: the PUBLISHED product, so a cell in this state is not current and WILL be re-run
#: on the next pass -- which is the loop this ticket exists to kill, and the only
#: honest thing to do is say so in the report rather than let it spin silently.
UNPUBLISHED = "UNPUBLISHED"
STATUSES = (DONE, SKIP, FAILED, BLOCKED, NOT_READY, WOULD_RUN, UNPUBLISHED)
#: A cell whose work exists. Everything else lands in the report's `resume` list.
TERMINAL_OK = (DONE, SKIP)


class MatrixError(RuntimeError):
    """The matrix itself cannot be built -- distinct from any cell failing."""


# ═════════════════════════════════════════════════════════════════════════════
# The engine axis
# ═════════════════════════════════════════════════════════════════════════════
#
# The DISPATCHABLE engine axis is the deck -- `derive_band_products.py
# --engine-b-deck` -- because that is the only engine knob a run actually turns.
# `model_registry.csv` is the axis of TREATMENTS, which are what a deck EMITS: a
# single `ts-lte` run writes both a 1D-LTE and an ENGINE-A product. The two are
# not the same list and pretending they were is what would fabricate cells.
#
# So the binding below is the one fact this file declares, and `assert_registry_covered`
# makes it impossible to hold it wrong quietly: every live token in the registry
# must be claimed by a deck or explicitly named un-dispatchable WITH A REASON, and
# every token a deck claims must still be live. A new registry row therefore fails
# loud here instead of silently never being run.
DECK_EMITS: dict[str, tuple[str, ...]] = {
    "ts-lte": ("1D-LTE", "ENGINE-A"),
    "gerber-1d-lte": ("synth-1D-LTE-gerber",),
    "gerber-nlte": ("ENGINE-B-NLTE",),
    "gerber-mean3d-lte": ("synth-mean3D-LTE-gerber-stagger",),
    "gerber-mean3d": ("synth-mean3D-NLTE-gerber-stagger",),
}

#: RYA-1040 -- `gerber-mean3d` and `gerber-mean3d-lte` are a MANDATORY PAIR: the
#: NLTE effect is (<3D>-NLTE minus <3D>-LTE) on ONE atmosphere, because
#: differencing against 1D-LTE would report the 1D -> mean-3D ATMOSPHERE shift as
#: non-LTE physics (RYA-542). Selecting one of them alone is refused.
DECK_PAIRS: dict[str, str] = {"gerber-mean3d": "gerber-mean3d-lte",
                              "gerber-mean3d-lte": "gerber-mean3d"}

NOT_DISPATCHABLE: dict[str, str] = {
    "ENGINE-A-3DNLTE": (
        "the Amarsi 3D-NLTE leg. It is a per-line grid interpolation applied to an "
        "EW measurement, not a `--engine-b-deck` synthesis, so no descriptor this "
        "layer emits can dispatch it. It is named here rather than omitted so that "
        "its absence from the matrix is a DECLARATION and not a gap (RYA-1187)."),
}


def _registry_rows() -> list[dict]:
    with MODEL_REGISTRY.open(encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def live_treatments() -> set[str]:
    """Every treatment token the model registry currently declares live.

    `status` is the registry's own gate: `in-dev` and `not-emitted` rows carry no
    stored token and must not become matrix cells.
    """
    return {r["stored_token"].strip() for r in _registry_rows()
            if r["status"].strip() == "live"
            and r["stored_token"].strip() not in ("", "-")}


def assert_registry_covered() -> None:
    """The deck binding and the model registry agree, or we stop.

    Loud both ways on purpose. An unclaimed live token is a treatment that would
    never be produced and never be missed; a claimed dead token is a deck that
    would be dispatched for a product nothing accepts.
    """
    live = live_treatments()
    claimed = {t for toks in DECK_EMITS.values() for t in toks}
    unclaimed = live - claimed - set(NOT_DISPATCHABLE)
    if unclaimed:
        raise MatrixError(
            f"model_registry.csv declares {sorted(unclaimed)} live, and no deck in "
            f"DECK_EMITS emits them and NOT_DISPATCHABLE does not name them. A live "
            f"treatment that no cell produces is a silent gap: bind it to a deck, or "
            f"name it un-dispatchable WITH THE REASON.")
    stale = claimed - live
    if stale:
        raise MatrixError(
            f"DECK_EMITS claims {sorted(stale)}, which model_registry.csv no longer "
            f"lists as live. Dispatching a deck for a retired treatment writes a "
            f"product with nowhere to land.")


def decks() -> tuple[str, ...]:
    """The engine axis, in a stable order. Checked against the registry first."""
    assert_registry_covered()
    return tuple(DECK_EMITS)


# ═════════════════════════════════════════════════════════════════════════════
# Adapters -- the two modules that must never be imported by name
# ═════════════════════════════════════════════════════════════════════════════

_ADAPTERS: dict = {}


def preflight():
    """RYA-1069's conductor, loaded BY PATH and never fatally.

    By path because `scripts/` holds ~90 modules whose names would shadow real
    packages on the next import anyone adds -- the reason RYA-1069 loads its own
    adapters the same way. Never fatally because `measure_band_ew` resolves the
    Kitt Peak atlas AT IMPORT and `SystemExit`s when it is absent: that is
    RYA-1064 friction #3, and inheriting it would stop an Alpha Cen matrix on
    unrelated SOLAR bytes.
    """
    if "preflight" not in _ADAPTERS:
        spec = importlib.util.spec_from_file_location(
            "_preflight_rya1222", ROOT / "scripts" / "preflight_readiness.py")
        mod = importlib.util.module_from_spec(spec)
        sys.modules["_preflight_rya1222"] = mod
        try:
            spec.loader.exec_module(mod)
            _ADAPTERS["preflight"], _ADAPTERS["preflight_why"] = mod, ""
        except BaseException as exc:            # SystemExit included, deliberately
            sys.modules.pop("_preflight_rya1222", None)
            _ADAPTERS["preflight"] = None
            _ADAPTERS["preflight_why"] = f"{type(exc).__name__}: {str(exc).splitlines()[0]}"
    return _ADAPTERS["preflight"]


def _publisher():
    """`scripts/publish_product.py` -- the SSOT for what a product's IDENTITY is.

    Imported rather than re-derived so this layer and the publisher cannot drift
    about which fields key a product (RYA-1127 added `line_set` to that key after
    two legs collided on it).
    """
    if "publisher" not in _ADAPTERS:
        spec = importlib.util.spec_from_file_location(
            "_publish_rya1222", ROOT / "scripts" / "publish_product.py")
        mod = importlib.util.module_from_spec(spec)
        sys.modules["_publish_rya1222"] = mod
        spec.loader.exec_module(mod)
        _ADAPTERS["publisher"] = mod
    return _ADAPTERS["publisher"]


# ═════════════════════════════════════════════════════════════════════════════
# The element and ion axes
# ═════════════════════════════════════════════════════════════════════════════

def canonical_elements() -> list[str]:
    """The 27 targets, from `elements_master.json`. Never a literal here."""
    doc = json.loads(ELEMENTS_MASTER.read_text(encoding="utf-8"))
    return [e["symbol"] for e in doc["elements"]]


def validate_element(element: str) -> str:
    """The element, or a refusal naming the vocabulary.

    `elements_master.json` lists `Fe` and `Fe II` as separate targets because the
    linelist writes both as element='Fe' and tells them apart by the ion column.
    This layer takes the SYMBOL and expands the ion axis itself, so `Fe II` is
    accepted and read as (element=Fe, ion=II) rather than being refused for a
    spelling the master list itself uses.
    """
    raw = element.strip()
    known = canonical_elements()
    if raw in known:
        return raw
    raise MatrixError(
        f"{element!r} is not one of the {len(known)} canonical targets in "
        f"{ELEMENTS_MASTER.relative_to(ROOT)}. Add it there before running it; do not "
        f"fall through to a silent default. Known: {', '.join(known)}")


def split_symbol(element: str) -> tuple[str, str | None]:
    """`Fe II` -> ('Fe', 'II'); `Si` -> ('Si', None). The master list's own spelling."""
    parts = element.split()
    return (parts[0], parts[1]) if len(parts) == 2 else (parts[0], None)


_ION_ORDER = ("I", "II", "III", "IV", "V", "VI")


def graded_ions(element: str) -> list[str]:
    """Ions of this element the GRADED pool actually holds a line for.

    The ion axis is the pool's, not a guess: running Si III because the periodic
    table permits it would fabricate a cell that no line can ever fill. The tier
    vocabulary is imported from `gf_empirical` rather than re-enumerated -- a
    second copy of that set is how a pool definition drifts (RYA-945).
    """
    key = ("ions", element)
    if key not in _ADAPTERS:
        import pandas as pd
        from pipeline.gf_empirical import GRADED_TIERS
        df = pd.read_csv(CANONICAL_GF, low_memory=False,
                         usecols=["species", "gf_tier"])
        df = df[df.gf_tier.astype(str).isin(GRADED_TIERS)]
        sp = df.species.astype(str).str.split()
        hit = {i for e, i in zip(sp.str[0], sp.str[-1]) if e == element}
        _ADAPTERS[key] = [i for i in _ION_ORDER if i in hit]
    return _ADAPTERS[key]


# ═════════════════════════════════════════════════════════════════════════════
# The expander
# ═════════════════════════════════════════════════════════════════════════════

def holdings_for(star: str) -> list[dict]:
    """Every registered holding of this system. Refuses rather than guessing one."""
    with HOLDINGS_REGISTRY.open(encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    hits = [r for r in rows if r["system_id"] == star]
    if not hits:
        known = sorted({r["system_id"] for r in rows})
        raise MatrixError(
            f"no holdings for system {star!r} in {HOLDINGS_REGISTRY.name} -- refuse, do "
            f"not guess one. Registered systems: {', '.join(known)}")
    return hits


def _bands_for(instrument: str) -> list[tuple[str, float, float]]:
    """(band, lo_A, hi_A) overlaps, from RYA-1069's own function.

    Deliberately NOT re-derived: the matrix and the preflight verdict must be
    keyed on the same (holding, band) pairs, or a NOT_READY row would be looked
    up for a cell that does not exist and every cell would read ready by default.
    """
    p = preflight()
    if p is None:
        raise MatrixError(
            f"the RYA-1069 conductor could not be imported "
            f"({_ADAPTERS.get('preflight_why')}), so the band axis has no source. "
            f"Refusing to substitute one.")
    return p.bands_for(instrument)


def _clip_to_holding(spec, lo: float, hi: float) -> tuple[float, float]:
    """Narrow the band window to what this HOLDING can serve. Empty means no cell.

    🔴 THIS IS WHY THE MATRIX IS CLIPPED AND NOT FILTERED. `HoldingSpec.covers()`
    demands TOTAL coverage -- a window half inside the product is a truncated
    window and serving it quietly would be worse than refusing. Asking it about
    the whole band therefore DELETED the single most productive holding in the
    repo: `bands_for('harps')` returns VIS as 3780-6910 A from the instrument
    catalogue, `solar_harps` declares its span as 3782.6-6910.0, and 2.4 A of
    disagreement at the blue edge made `covers()` false. Every published HARPS
    VIS product -- the Fe anchor among them -- sat in a cell the expander had
    silently dropped, which is precisely the silent gap this ticket exists to
    close (RYA-1187).

    The instrument axis already works this way: RYA-1069's `bands_for` returns
    the OVERLAP of the band and the instrument's reach rather than testing for
    containment. This applies the same rule one level down, so a cell's window is
    band edges INTERSECT instrument reach INTERSECT holding span. A holding that
    inventories its own coverage (`span_A is None` -- the Kitt Peak segment list,
    the CRIRES+ IDP comb) is left alone: it has no declared span to clip to, and
    the resolver still asks it `covers()` before anything runs.
    """
    if spec is None or getattr(spec, "span_A", None) is None:
        return lo, hi
    return max(lo, spec.span_A[0]), min(hi, spec.span_A[1])


def expand(star: str, element: str, *, ions: list[str] | None = None,
           bands: list[str] | None = None, instruments: list[str] | None = None,
           engines: list[str] | None = None) -> list[RunDescriptor]:
    """Every applicable (band x holding x ion x deck) cell, in a stable sorted order.

    A cell is omitted ONLY where the wired `HoldingSpec` positively declares a
    span that does not contain the band. That is a declaration, not an absence,
    and the resolver would reach the same verdict -- omitting it here keeps the
    report about cells that could conceivably run. Everything else is kept and
    documented, including a holding with no HoldingSpec at all: "registered and
    unreachable" is the RYA-904 shape and it must not read as "no data".
    """
    element = validate_element(element)
    symbol, pinned_ion = split_symbol(element)

    want_ions = [pinned_ion] if pinned_ion else (ions or graded_ions(symbol))
    if not want_ions:
        raise MatrixError(
            f"the graded pool in {CANONICAL_GF.name} holds no line of {symbol!r} at any "
            f"ion, so there is no cell to build. This is a POOL finding, not a run "
            f"failure -- grade some lines first (RYA-945).")

    want_decks = list(engines) if engines else list(decks())
    unknown = [d for d in want_decks if d not in DECK_EMITS]
    if unknown:
        raise MatrixError(
            f"unknown engine deck(s) {unknown}; the axis is {sorted(DECK_EMITS)} "
            f"(bound to model_registry.csv by DECK_EMITS).")
    missing_pair = [d for d in want_decks
                    if d in DECK_PAIRS and DECK_PAIRS[d] not in want_decks]
    if missing_pair:
        raise MatrixError(
            f"{missing_pair} selected without its comparand "
            f"{[DECK_PAIRS[d] for d in missing_pair]}. RYA-1040/542: the <3D> decks are a "
            f"MANDATORY PAIR -- the NLTE effect is (<3D>-NLTE minus <3D>-LTE) on ONE "
            f"atmosphere, and differencing the survivor against 1D-LTE would report the "
            f"1D -> mean-3D ATMOSPHERE shift as non-LTE physics.")

    p = preflight()
    out: list[RunDescriptor] = []
    for h in holdings_for(star):
        inst, hid = h["instrument_id"], h["holding_id"]
        if instruments and inst not in instruments and hid not in instruments:
            continue
        spec = p.holding_spec(hid) if p is not None else None
        for band, lo, hi in _bands_for(inst):
            if bands and band not in bands:
                continue
            lo, hi = _clip_to_holding(spec, lo, hi)
            if hi <= lo:
                continue                      # declared out of span -- rule 2's one exit
            for ion in want_ions:
                for deck in want_decks:
                    out.append(RunDescriptor(
                        element=symbol, ion=ion, instrument=inst, holding=hid,
                        lo_A=lo, hi_A=hi, engine_deck=deck))
    out.sort(key=lambda d: (d.lo_A, d.instrument, d.holding, d.ion, d.engine_deck))
    return out


# ═════════════════════════════════════════════════════════════════════════════
# The idempotency key
# ═════════════════════════════════════════════════════════════════════════════

_HASH_CACHE: dict = {}


def _file_fingerprint(path: Path) -> str:
    """sha256 of the file's bytes, or a NAMED absence. Never silently nothing.

    Content, not mtime: a re-checkout or a `touch` moves the mtime without moving
    the input, and re-running the whole matrix for that is precisely the waste
    this ticket exists to end. Memoised on (path, mtime, size) so a 24 MB pool
    file is read once per run rather than once per cell.
    """
    try:
        st = path.stat()
    except OSError as exc:
        return f"ABSENT:{type(exc).__name__}"
    key = (str(path), st.st_mtime_ns, st.st_size)
    if key not in _HASH_CACHE:
        h = hashlib.sha256()
        with path.open("rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                h.update(chunk)
        _HASH_CACHE[key] = h.hexdigest()
    return _HASH_CACHE[key]


def input_fingerprints(descriptor: RunDescriptor, resolved, *,
                       manifest_path: str | None) -> list[dict]:
    """Every named input this cell depends on, each with the bytes it was hashed on.

    Returned as a LIST OF NAMED ROWS rather than one opaque digest so that a
    changed hash can be explained -- "it re-ran and I cannot tell you why" is how
    an idempotency layer stops being trusted. What is in here:

      * the descriptor itself (star is not a field of it; the caller's holding is)
      * each ordered stage SCRIPT's own bytes -- this is the `code_commit of the
        stage modules` the spec asks for, done exactly. Repo HEAD would move for
        edits to files this cell never loads, and would not move at all for an
        uncommitted edit that changes the answer.
      * the holding's manifest, which is the registry's own name for this
        holding's input artifact
      * the graded line pool and the three ledgers the axes were read off
    """
    rows: list[dict] = [
        {"kind": "descriptor", "name": "run_descriptor",
         "digest": hashlib.sha256(
             json.dumps(descriptor.as_dict(), sort_keys=True).encode()).hexdigest()},
        {"kind": "method", "name": "resolved_method", "digest": resolved.method},
    ]
    for step in resolved.steps:
        script = ROOT / step["script"]
        rows.append({"kind": "stage_script", "name": step["script"],
                     "digest": _file_fingerprint(script)})
    if manifest_path:
        rows.append({"kind": "holding_manifest", "name": manifest_path,
                     "digest": _file_fingerprint(ROOT / manifest_path)})
    for ledger in (CANONICAL_GF, MODEL_REGISTRY, HOLDINGS_REGISTRY,
                   ROOT / "data" / "catalog" / "instrument_catalog.csv"):
        rows.append({"kind": "ledger", "name": str(ledger.relative_to(ROOT)),
                     "digest": _file_fingerprint(ledger)})
    return rows


def inputs_hash(descriptor: RunDescriptor, resolved, *,
                manifest_path: str | None) -> str:
    rows = input_fingerprints(descriptor, resolved, manifest_path=manifest_path)
    blob = json.dumps(rows, sort_keys=True).encode()
    return hashlib.sha256(blob).hexdigest()


# ═════════════════════════════════════════════════════════════════════════════
# The published feed -- what "already current" means
# ═════════════════════════════════════════════════════════════════════════════

def load_feed(star: str, element: str) -> dict:
    """The published element feed, via the publisher's own loader."""
    return _publisher().load(element, star)


def cell_prefix_key(descriptor: RunDescriptor) -> str:
    """The part of a product's identity a CELL fixes.

    A cell is not one product: a `ts-lte` run emits both a 1D-LTE and an ENGINE-A
    product, and the tier/selector/route the run lands on are OUTPUTS, not inputs.
    So the cell is keyed on the identity fields the descriptor actually pins, and
    the treatment axis is matched separately against the deck's emitted set.
    """
    return "|".join((descriptor.element, descriptor.ion, descriptor.band,
                     descriptor.instrument, descriptor.holding))


def cell_products(feed: dict, descriptor: RunDescriptor, band: str) -> list[dict]:
    """The published products this cell owns, matched on identity + the deck's tokens."""
    emits = set(DECK_EMITS[descriptor.engine_deck])
    return [p for p in feed.get("products", [])
            if str(p.get("element")) == descriptor.element
            and str(p.get("ion")) == descriptor.ion
            and str(p.get("band")) == band
            and str(p.get("instrument")) == descriptor.instrument
            and str(p.get("holding")) == descriptor.holding
            and str(p.get("treatment")) in emits]


#: Where the inputs-hash lives.
#:
#: 🔴 DEVIATION FROM THE RYA-1222 SPEC, WITH THE MEASUREMENT THAT FORCED IT. The
#: spec says to store the hash ON the product in `data/products/<star>/<El>.json`.
#: That is refused by the RYA-587 publication contract. `publish_product.write_feed`
#: calls `uncertainty_contract.assert_publication_feed(doc, previous=...)`, which
#: re-validates every product that is not BYTE-IDENTICAL to the row already on
#: disk. All 160 committed solar Fe products are legacy rows with no `uncertainty`
#: block; they pass today only via the identity retention clause. Measured on the
#: real feed: the unchanged document passes, and the same document with ONE
#: `inputs_hash` key added on ONE product is refused with
#:   "Fe|II|VIS|harps|solar_harps_molecfit_corrected|...: star: provenance is required".
#: So writing the hash where the spec asks would either fail every DONE cell or
#: require weakening a science gate to let an orchestrator bookkeeping field in.
#: Neither is acceptable, so the hash lives in this sidecar instead: same key, same
#: comparison, tracked in git, and NOTHING of this layer's writes touches the
#: science feed. The spec's actual requirement -- that a second run does no work --
#: is met exactly.
LEDGER = REPORT_DIR / "inputs_hashes.json"


def _product_key(product: dict) -> str:
    """The publisher's own identity key, or a NAMED fallback -- never an exception.

    `publish_product.key_of` resolves `line_set` through `reference_lineset`, which
    REFUSES to guess when a product carries neither an explicit `line_set` nor a
    recognised `tier`. That refusal is right for a publisher and wrong here: this
    list is bookkeeping recorded AFTER a cell's stages have already succeeded, and
    letting it raise would turn a completed run into a crash and lose the cell.
    The unresolvable case is recorded as what it is.
    """
    try:
        return _publisher().key_of(product)
    except Exception as exc:                                   # noqa: BLE001
        return (f"UNRESOLVED-KEY({type(exc).__name__}): "
                + "|".join(str(product.get(k) or "") for k in
                           ("element", "ion", "band", "instrument", "holding", "treatment")))


def ledger_key(star: str, descriptor: RunDescriptor, band: str) -> str:
    """The cell's identity in the sidecar. Carries the deck, which the product key
    does not: two decks over one window are two runs and two hashes."""
    return "|".join((star, descriptor.element, descriptor.ion, band,
                     descriptor.instrument, descriptor.holding, descriptor.engine_deck))


def load_ledger() -> dict:
    if not LEDGER.exists():
        return {}
    try:
        return json.loads(LEDGER.read_text(encoding="utf-8"))
    except (ValueError, OSError) as exc:
        # Loud, and treated as EMPTY rather than as "everything is current": a
        # corrupt ledger must cost a re-run, never a skipped one.
        print(f"WARNING: {LEDGER.name} is unreadable ({type(exc).__name__}: {exc}); "
              f"every cell will be treated as NOT current and re-run.", file=sys.stderr)
        return {}


def record_inputs_hash(key: str, want_hash: str, *, product_keys: list[str],
                       code_commit: str) -> None:
    """Record what this cell was built from, after it was built."""
    led = load_ledger()
    led[key] = {"inputs_hash": want_hash, "recorded_at": _utc(),
                "code_commit": code_commit, "products": sorted(product_keys)}
    LEDGER.parent.mkdir(parents=True, exist_ok=True)
    LEDGER.write_text(json.dumps(led, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def is_current(products: list[dict], want_hash: str,
               recorded: dict | None) -> tuple[bool, str]:
    """(is this cell's work current, why-not).

    BOTH halves are required. The published products must still be there -- a
    ledger entry for a product someone quarantined would skip a cell that no
    longer has an answer -- and the recorded hash must match. A cell with no
    ledger entry is NOT current: it was built before this layer existed, so
    nothing recorded what it was built from, and treating "we cannot tell" as
    "yes" is exactly the silent fallback the permanent rules forbid. It is re-run
    once, and that run records the hash.
    """
    if not products:
        return False, "no published product for this cell"
    if not recorded:
        return False, (f"{len(products)} published product(s) exist but no inputs_hash was "
                       f"ever recorded for this cell -- built before this layer, so what "
                       f"they were built from is unknown")
    if recorded.get("inputs_hash") != want_hash:
        return False, (f"recorded inputs_hash {str(recorded.get('inputs_hash'))[:12]} != "
                       f"{want_hash[:12]} -- an input or a stage script moved since "
                       f"{recorded.get('recorded_at')}")
    return True, (f"{len(products)} published product(s) current at inputs_hash "
                  f"{want_hash[:12]} recorded {recorded.get('recorded_at')}")


# ═════════════════════════════════════════════════════════════════════════════
# The executor
# ═════════════════════════════════════════════════════════════════════════════

@dataclass
class CellResult:
    """One cell's terminal verdict. Every field the report prints comes from here."""
    band: str
    instrument: str
    holding: str
    ion: str
    engine: str
    status: str
    reason: str
    A: float | None = None
    n_lines: int | None = None
    product_key: str = ""
    inputs_hash: str = ""
    lo_A: float = 0.0
    hi_A: float = 0.0
    steps_run: list[str] = field(default_factory=list)

    @property
    def cell_key(self) -> str:
        return "|".join((self.band, self.instrument, self.holding, self.ion, self.engine))


def _utc() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _code_commit() -> str:
    try:
        out = subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"],
                             capture_output=True, text=True, timeout=30)
        sha = out.stdout.strip()
        dirty = subprocess.run(["git", "-C", str(ROOT), "status", "--porcelain"],
                               capture_output=True, text=True, timeout=60).stdout.strip()
        return f"{sha}{'-dirty' if dirty else ''}" if sha else "unknown"
    except Exception as exc:                                   # noqa: BLE001
        return f"unknown ({type(exc).__name__})"


def verify_numpy_ceiling(interpreter: str) -> tuple[bool, str]:
    """The RYA-682 ceiling, CHECKED rather than asserted.

    `run_descriptor.resolve` records that the caller pinned an interpreter and says
    in as many words that the EXECUTOR must verify the version. This is that check.
    Above the ceiling iSpec writes a zero-row artifact and exits 0 -- a control run
    that reports success on no data -- so it is checked before dispatch, never
    discovered when a product comes back empty.
    """
    from pipeline.run_descriptor import NUMPY_CEILING
    ceiling = ".".join(map(str, NUMPY_CEILING))
    try:
        out = subprocess.run([interpreter, "-c", "import numpy; print(numpy.__version__)"],
                             capture_output=True, text=True, timeout=120)
    except Exception as exc:                                   # noqa: BLE001
        return False, f"{interpreter} could not be run: {type(exc).__name__}: {exc}"
    if out.returncode != 0:
        return False, (f"{interpreter} cannot import numpy: "
                       f"{(out.stderr or out.stdout).strip().splitlines()[-1:] or ['(silent)']}")
    got = out.stdout.strip()
    try:
        parts = tuple(int(x) for x in got.split(".")[:2])
    except ValueError:
        return False, f"{interpreter} reported an unparseable numpy version {got!r}"
    if parts >= NUMPY_CEILING:
        return False, (f"numpy {got} on {interpreter} is AT OR ABOVE the {ceiling} ceiling "
                       f"(RYA-682): iSpec would write a ZERO-ROW artifact and exit 0, so "
                       f"this refusal is the only thing between us and a silent empty "
                       f"product reported as a success")
    return True, f"numpy {got} < {ceiling} on {interpreter}"


def _readiness_index(star: str, element: str, ion: str | None):
    """{(holding, band): ReadinessRow} from RYA-1069, or (None, why).

    One call per (star, element, ion): the conductor re-reads every manifest and the
    graded pool, so calling it per cell would multiply a ~7 s scan by the matrix.
    """
    p = preflight()
    if p is None:
        return None, (f"the RYA-1069 conductor could not be imported: "
                      f"{_ADAPTERS.get('preflight_why')}")
    key = ("readiness", star, element, ion)
    if key not in _ADAPTERS:
        evidence: dict = {}
        rows = [r for h in p.holdings_for_system(star, None)
                for r in p.assess(h, element, ion, evidence)]
        _ADAPTERS[key] = {(r.holding_id, r.band): r for r in rows}
    return _ADAPTERS[key], ""


def _manifest_paths(star: str) -> dict[str, str]:
    return {h["holding_id"]: h.get("manifest_path", "") for h in holdings_for(star)}


def _accepts_star(script: Path) -> bool:
    """Does this stage script declare a `--star` flag?

    Asked of the script's own source rather than of `--help`, because several
    stage scripts do work at import and `--help` is not a safe probe on them.
    It matters: `derive_band_products.py` takes `--star` and DEFAULTS IT TO
    'solar', while `measure_band_profilefit.py` has no star axis at all -- the
    instrument and holding already pin whose spectrum it is. Passing the flag
    blindly would make every profile-fit step die on an unrecognised argument,
    and omitting it from the derive step would silently measure the Sun for
    whichever star was asked for. Neither is allowed to happen quietly.
    """
    key = ("accepts_star", str(script))
    if key not in _ADAPTERS:
        try:
            _ADAPTERS[key] = '"--star"' in script.read_text(encoding="utf-8")
        except OSError:
            _ADAPTERS[key] = False
    return _ADAPTERS[key]


def _dispatch_signature(step: dict, star: str) -> str:
    """Exactly what this step will execute: interpreter, script, args and env overrides.

    Two steps with this signature equal do the same work on the same inputs, so the
    second one is free. Nothing weaker is safe -- see the reuse check in `run`.
    """
    return json.dumps([step.get("interpreter") or sys.executable, step["script"], star,
                       list(step["args"]), sorted((step.get("env") or {}).items())],
                      sort_keys=True)


def _run_step(step: dict, star: str, *, timeout: int) -> tuple[bool, str]:
    """Dispatch one ordered step. Returns (ok, detail). Raises nothing."""
    interp = step.get("interpreter") or sys.executable
    script = ROOT / step["script"]
    cmd = [interp, str(script)]
    if _accepts_star(script):
        cmd += ["--star", star]
    cmd += list(step["args"])
    env = dict(os.environ)
    env.update({k: v for k, v in (step.get("env") or {}).items() if v})
    try:
        out = subprocess.run(cmd, capture_output=True, text=True, cwd=str(ROOT),
                             env=env, timeout=timeout)
    except subprocess.TimeoutExpired:
        return False, f"{step['name']}: TIMEOUT after {timeout}s"
    except Exception as exc:                                   # noqa: BLE001
        return False, f"{step['name']}: {type(exc).__name__}: {exc}"
    if out.returncode != 0:
        tail = (out.stderr or out.stdout or "").strip().splitlines()[-4:]
        return False, f"{step['name']}: exit {out.returncode}: {' | '.join(tail)}"
    return True, f"{step['name']}: exit 0"


def run(star: str, element: str, *, ions: list[str] | None = None,
        bands: list[str] | None = None, instruments: list[str] | None = None,
        engines: list[str] | None = None, dry_run: bool = False,
        interpreter: str | None = None, ispec_dir: str | None = None,
        step_timeout: int = 7200, report_dir: Path | None = None,
        echo: bool = True) -> dict:
    """Drive the full matrix for one (star, element). Returns the run report.

    Loud-fail-CONTINUE is the whole contract: a cell that blocks, is not ready, or
    breaks is recorded with its reason and the next cell is attempted. The only
    thing that raises out of here is a failure to BUILD the matrix at all, which
    is not a cell verdict and must not be reported as one.
    """
    element = validate_element(element)
    symbol, pinned_ion = split_symbol(element)
    interpreter = interpreter or os.environ.get("CODEX_INTERPRETER") or ""
    ispec_dir = ispec_dir or os.environ.get("ISPEC_DIR") or ""

    descriptors = expand(star, element, ions=ions, bands=bands,
                         instruments=instruments, engines=engines)

    # The RYA-682 ceiling once per run, not once per cell -- but its verdict is
    # carried INTO every cell, because a failed check must block dispatch rather
    # than warn beside it.
    numpy_ok, numpy_why = (True, "no interpreter pinned")
    if interpreter:
        numpy_ok, numpy_why = verify_numpy_ceiling(interpreter)

    ready_ix, ready_why = _readiness_index(star, symbol, pinned_ion or (ions[0] if ions else None))
    manifests = _manifest_paths(star)
    ledger = load_ledger()
    commit = _code_commit()
    #: Artifacts this invocation has already built successfully, by the resolver's own
    #: `produces` path. Scoped to the process on purpose -- see the reuse check below.
    produced_this_run: set[str] = set()

    cells: list[CellResult] = []
    for d in descriptors:
        band = _band_of(d, star)
        cell = CellResult(band=band.name_declared, instrument=d.instrument,
                          holding=d.holding, ion=d.ion, engine=d.engine_deck,
                          status=BLOCKED, reason="", lo_A=d.lo_A, hi_A=d.hi_A,
                          product_key=cell_prefix_key(d))
        cells.append(cell)

        # ── 0. the two band tables must agree about where this cell IS ───────
        if band.disagreement:
            cell.status, cell.reason = BLOCKED, band.disagreement
            continue

        # ── 1. can this run be expressed at all (RYA-767) ────────────────────
        try:
            resolved = resolve(d, interpreter=interpreter or None,
                               ispec_dir=ispec_dir or None)
        except Exception as exc:                               # noqa: BLE001
            cell.status, cell.reason = BLOCKED, f"resolver raised {type(exc).__name__}: {exc}"
            continue
        if resolved.blocked_reason:
            cell.status, cell.reason = BLOCKED, resolved.blocked_reason
            continue
        unmet = [p for p in resolved.preconditions if not p.satisfied]
        if unmet:
            cell.status = BLOCKED
            cell.reason = "; ".join(f"{p.name}: {p.detail}" for p in unmet)
            continue

        # ── 2. is the DATA ready (RYA-1069) ──────────────────────────────────
        if ready_ix is None:
            cell.status, cell.reason = NOT_READY, ready_why
            continue
        row = ready_ix.get((d.holding, cell.band))
        if row is None:
            cell.status, cell.reason = NOT_READY, (
                f"the RYA-1069 conductor returned no verdict for "
                f"({d.holding}, {cell.band}). An unassessed cell is NOT a ready one.")
            continue
        if row.measurement_ready != "GO":
            cell.status, cell.reason = NOT_READY, (
                f"{row.measurement_ready} [{row.blocking_gate}]")
            continue

        # ── 3. is the work already current (the loop-killer) ─────────────────
        want = inputs_hash(d, resolved, manifest_path=manifests.get(d.holding))
        cell.inputs_hash = want
        feed = load_feed(star, symbol)
        published = cell_products(feed, d, cell.band)
        lkey = ledger_key(star, d, cell.band)
        current, why = is_current(published, want, ledger.get(lkey))
        if published:
            cell.A = published[0].get("A")
            cell.n_lines = published[0].get("n_lines")
        if current:
            cell.status, cell.reason = SKIP, why
            continue

        # ── 4. execute, loud per cell, never aborting the matrix ─────────────
        if not numpy_ok:
            cell.status, cell.reason = BLOCKED, f"numpy ceiling: {numpy_why}"
            continue
        if dry_run:
            cell.status, cell.reason = WOULD_RUN, f"would execute -- {why}"
            cell.steps_run = [s["name"] for s in resolved.steps]
            continue
        ok, detail = True, ""
        for step in resolved.steps:
            # 🔴 WITHIN ONE RUN, DO NOT RE-DISPATCH AN IDENTICAL COMMAND.
            # The EW step is DECK-INDEPENDENT: `descriptor.key` carries no deck, so
            # all five decks over one (holding, band) invoke `measure_band_profilefit`
            # with byte-identical arguments. Measured on the first full Si run:
            # `measure_ew` executed 15 times for 3 distinct commands, and the 12 extra
            # runs rewrote bytes that were already correct.
            #
            # Keyed on the COMMAND, not on the step's declared `produces`. Keying on
            # `produces` is what an earlier version of this did and it was WRONG, for
            # a reason worth writing down: `run_descriptor.resolve` builds the derive
            # step's `produces` from `descriptor.key` too, so ALL FIVE DECKS DECLARE
            # THE SAME PRODUCTS ARTIFACT even though each writes a different product.
            # Reusing on that string silently skipped `derive_products` for four decks
            # out of five while the report said they were handled -- a silent gap
            # produced by the thing built to abolish silent gaps. An identical command
            # is identical work by construction; a shared output path is not.
            # (The `produces` collision itself is a defect in the resolver, not here.)
            sig = _dispatch_signature(step, star)
            if sig in produced_this_run:
                cell.steps_run.append(f"{step['name']}:reused")
                continue
            ok, detail = _run_step(step, star, timeout=step_timeout)
            cell.steps_run.append(f"{step['name']}:{'ok' if ok else 'FAILED'}")
            if not ok:
                break
            produced_this_run.add(sig)
        if not ok:
            cell.status, cell.reason = FAILED, detail
            continue
        feed = load_feed(star, symbol)
        published = cell_products(feed, d, cell.band)
        if not published:
            cell.status = UNPUBLISHED
            cell.reason = (
                f"all {len(resolved.steps)} stage(s) exited 0 and wrote their artifact, but "
                f"no product for this cell is in data/products/{star}/{symbol}.json. "
                f"Publication is `scripts/publish_product.py` and is human-gated by design "
                f"(RYA-1034 needs a stated reason; RYA-772 forbids copying an artifact "
                f"without provenance). Until it is published this cell is NOT current and "
                f"WILL re-run.")
            continue
        cell.A = published[0].get("A")
        cell.n_lines = published[0].get("n_lines")
        record_inputs_hash(lkey, want, code_commit=commit,
                           product_keys=[_product_key(p) for p in published])
        ledger[lkey] = {"inputs_hash": want}
        cell.status = DONE
        cell.reason = (f"ran {len(resolved.steps)} step(s); {len(published)} published "
                       f"product(s) recorded at inputs_hash {want[:12]}")

    report = _report(star, element, cells, dry_run=dry_run, numpy_why=numpy_why,
                     interpreter=interpreter, ispec_dir=ispec_dir,
                     code_commit=commit, report_dir=report_dir)
    if echo:
        print(render(report))
    return report


@dataclass
class _Band:
    name_declared: str
    disagreement: str = ""


def _band_of(d: RunDescriptor, star: str) -> _Band:
    """The cell's band, and a refusal when the repo's two band tables disagree.

    🔴 `config/synth_bands.yaml` and `pipeline/band_policy.py` do NOT share edges.
    synth_bands ends red-optical at 9199 A and band_policy at 10000 A. A holding
    whose reach clips the NIR band to 9199-10650 therefore has a MIDPOINT of
    9924.5 A, which `band_policy.resolve` -- and with it `method_for` -- calls
    red-optical, where profile-fit is permitted. The cell would have been
    dispatched under the wrong regime's method with nothing in the output saying
    so. This layer does not get to pick a winner between two SSOTs, so it names
    the disagreement and blocks the cell (RYA-1187: a documented cell, not a
    silent gap). The fix is a band-edge reconciliation ticket, not a tiebreak here.
    """
    # By CONTAINMENT, not equality: the window was clipped to the holding's span,
    # so it is a sub-interval of the band it belongs to rather than the band itself.
    from config import synth_bands
    declared = ""
    for name, b in sorted(synth_bands.SYNTH_BANDS.items(), key=lambda kv: kv[1].lo_A):
        if b.lo_A <= d.lo_A and d.hi_A <= b.hi_A:
            declared = name
            break
    policy = band_policy.resolve(0.5 * (d.lo_A + d.hi_A)).name
    if declared and declared != policy:
        return _Band(declared, (
            f"BAND-EDGE DISAGREEMENT: config/synth_bands.yaml places "
            f"{d.lo_A:g}-{d.hi_A:g} A in {declared!r}, while band_policy.resolve() on its "
            f"midpoint {0.5 * (d.lo_A + d.hi_A):g} A returns {policy!r} -- and the METHOD "
            f"comes from band_policy. Dispatching would run this cell under the "
            f"{policy!r} regime's rules while every artifact called it {declared!r}. "
            f"Refusing to tiebreak between two source-of-truth tables."))
    return _Band(declared or policy)


# ═════════════════════════════════════════════════════════════════════════════
# The report
# ═════════════════════════════════════════════════════════════════════════════

def _report(star: str, element: str, cells: list[CellResult], *, dry_run: bool,
            numpy_why: str, interpreter: str, ispec_dir: str, code_commit: str,
            report_dir: Path | None) -> dict:
    counts = {s: sum(c.status == s for c in cells) for s in STATUSES}
    doc = {
        "schema": "codex.orchestrator_run/1",
        "star": star,
        "element": element,
        "generated_at": _utc(),
        "code_commit": code_commit,
        "dry_run": bool(dry_run),
        "environment": {"interpreter": interpreter or None,
                        "ispec_dir": ispec_dir or None,
                        "numpy_ceiling_check": numpy_why},
        "cells_total": len(cells),
        "counts": counts,
        "cells": [{"band": c.band, "instrument": c.instrument, "holding": c.holding,
                   "ion": c.ion, "engine": c.engine, "status": c.status,
                   "reason": c.reason, "A": c.A, "n_lines": c.n_lines,
                   "product_key": c.product_key, "inputs_hash": c.inputs_hash,
                   "lo_A": c.lo_A, "hi_A": c.hi_A, "steps": c.steps_run}
                  for c in cells],
        "resume": [c.cell_key for c in cells if c.status not in TERMINAL_OK],
    }
    # Rule 2, checked rather than trusted: the counts must account for every cell
    # and every cell must carry one of the five statuses. A cell that fell through
    # a branch without a verdict is a BUILD DEFECT and says so here, loudly.
    assert sum(counts.values()) == len(cells), "a cell carries no terminal status"
    dest = report_dir or REPORT_DIR
    dest.mkdir(parents=True, exist_ok=True)
    stamp = doc["generated_at"].replace(":", "").replace("-", "")
    safe = element.replace(" ", "")
    (dest / f"{star}_{safe}_{stamp}.json").write_text(
        json.dumps(doc, indent=2) + "\n", encoding="utf-8")
    (dest / f"{star}_{safe}_latest.json").write_text(
        json.dumps(doc, indent=2) + "\n", encoding="utf-8")
    latest = dest / f"{star}_{safe}_latest.json"
    # Repo-relative where it is inside the repo, absolute where a caller redirected it
    # (a test's tmpdir). RYA-1096: an absolute path is a MACHINE fact, so the committed
    # shape stays relative -- but refusing to report a path at all would be worse.
    try:
        doc["_report_path"] = str(latest.relative_to(ROOT))
    except ValueError:
        doc["_report_path"] = str(latest)
    return doc


def render(doc: dict) -> str:
    """The compact ASCII table plus the one-line verdict. ASCII only (Cloudflare WAF)."""
    head = (f"{'band':<12} {'instrument':<24} {'holding':<34} {'ion':<4} "
            f"{'engine':<18} {'status':<10} {'A':>8} {'n':>5}  reason")
    lines = [
        f"ORCHESTRATOR  star={doc['star']}  element={doc['element']}  "
        f"{'DRY-RUN' if doc['dry_run'] else 'EXECUTE'}  commit={doc['code_commit']}",
        f"numpy ceiling: {doc['environment']['numpy_ceiling_check']}",
        "",
        head, "-" * len(head),
    ]
    for c in doc["cells"]:
        a = f"{c['A']:.4f}" if isinstance(c["A"], (int, float)) else "-"
        n = str(c["n_lines"]) if c["n_lines"] is not None else "-"
        reason = (c["reason"] or "").replace("\n", " ")
        if len(reason) > 96:
            reason = reason[:93] + "..."
        lines.append(f"{c['band']:<12} {c['instrument']:<24} {c['holding']:<34} "
                     f"{c['ion']:<4} {c['engine']:<18} {c['status']:<10} {a:>8} {n:>5}  {reason}")
    # Every status in STATUSES, including the zeroes: a verdict line that prints only
    # the non-zero counts cannot be read as "and nothing else happened", which is the
    # one thing an operator needs it to say.
    k = doc["counts"]
    lines += [
        "-" * len(head),
        " / ".join(f"{s} {k[s]}" for s in STATUSES) + f" across {doc['cells_total']} cells",
        f"report: {doc.get('_report_path', '')}",
    ]
    return "\n".join(lines)
