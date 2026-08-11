#!/usr/bin/env python3
"""
pipeline/coverage.py — RYA-708
==============================
"IS THIS LINE OBSERVED?" IS A QUESTION ABOUT THE CODEX, NOT ABOUT ONE FILE.

The defect this closes
----------------------
RYA-707's unresolved-element appendix reported Al I 7835/7836 and 8772/8773 as
**NO DATA, outside our wavelength coverage, zero pixels**. That verdict reached the
state register (v48) and a published artifact.

It was false. It was true of HARPS (3782.6-6910.0 A) and false of the project: the IAG
Reiners+2016 FTS atlas covers 4047.4-10649.9 A, and the deblend harnesses have been
loading it as a SECOND ARM since RYA-551. Al I 8772/8773 is 36.7% deep there and carries
log gf -0.161 — the strongest line in the element, and the appendix said we could not
see it.

The cause was not a bad measurement. `plot_unmeasured_lines_rya707.py` defined "our
coverage" as `data/processed/solar_normalized.csv` — one hardcoded file — and printed a
claim about that file as a claim about the Codex. The IAG path was independently
hardcoded in four more scripts (RYA-551/564/581/592), so nothing in the repo knew what
we own.

What this module is
-------------------
The single source for that question. It answers, for a star and a wavelength, WHICH
instruments cover it — never whether one particular file does.

IT OWNS NO INSTRUMENT DATA. Ryan, 2026-08-09: *"we also have an instruments.csv that is
also our source for the instruments and the wavelengths they cover."* He was right, and
the first version of this module shipped a three-row registry of its own — duplicating a
25-instrument, 32-column catalog that already had a loader, a validator and a test.
`kpno_solar_atlas` was already in it at 296-1300 nm, which is exactly the 2960-13000 A
this module then went and measured. Building a second source while fixing a
single-source bug is the joke writing itself, so the duplicate is deleted.

The division of labour, which is the catalog's own and not invented here:

  * `data/catalog/instrument_catalog.csv` — WHICH INSTRUMENTS EXIST and what they cover.
  * `data/catalog/instrument_modes.csv`   — per-mode coverage where it differs.
  * `data/catalog/holdings_manifest_registry.csv` — WHAT WE HOLD, per system.

This module joins them. Where it needs something the catalog genuinely does not carry —
the path to a solar spectrum on this machine and how to read it — that goes in the
HOLDINGS registry, because that is what holdings means.

Three distinctions it refuses to collapse, because collapsing them is how the wrong
answer got out:

  * **not covered** — no instrument we own reaches this wavelength. A real data gap,
    and the only state that justifies "NO DATA" in an appendix.
  * **covered, not reachable here** — an instrument covers it, but the file lives on
    the other machine. IAG is Sirius-only; a line is not a data gap because you asked
    from the Mac. Callers get `covered=True` and `reachable=False`.
  * **covered and loadable** — the file is present and can be read now.

Declared vs verified spans
--------------------------
A registry row carries a span and a `span_status`. `--verify` measures the file where it
is reachable and REPORTS disagreement; it never rewrites a declared span with a measured
one. A registry that silently self-corrects cannot be audited, and the point of this
module is that its answers are checkable.

Two coverage questions, one discipline (RYA-776)
------------------------------------------------
"What covers this wavelength" has two independent halves, and this module now answers
both with the same refusal to collapse states:

  * WHICH INSTRUMENT sees it   — the first half, above (RYA-708).
  * WHICH ENGINE/GRID models it — the second half, at the bottom of this file, reading
    the generated `data/catalog/engine_coverage.csv` (RYA-776).

They are complementary and neither substitutes for the other: an instrument can see a
line no departure grid reaches, and a grid can reach a line no instrument we hold
observes. See the ENGINE LAYER banner below for its states.
"""
from __future__ import annotations

import argparse
import os
import sys
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
CATALOG  = ROOT / "data" / "catalog" / "instrument_catalog.csv"
MODES    = ROOT / "data" / "catalog" / "instrument_modes.csv"
HOLDINGS = ROOT / "data" / "catalog" / "holdings_manifest_registry.csv"

VALID_HOSTS = frozenset({"mac", "sirius", "both"})
VALID_STATUS = frozenset({"VERIFIED", "DECLARED"})


class CoverageError(RuntimeError):
    """A catalog or holdings file is missing or malformed. Never degraded to 'no
    coverage' — an unreadable registry must not be indistinguishable from a genuine
    data gap, which is exactly the confusion this module exists to end."""


@dataclass(frozen=True)
class Instrument:
    """One instrument WE HOLD DATA FOR, for one star: the catalog row joined to the
    holdings row. Spans come from the holdings side because that is the file we would
    actually open; the catalog span is carried alongside so the two can disagree
    visibly rather than silently."""
    star: str
    instrument: str
    instrument_id: str
    host: str
    path: str
    loader: str
    wave_min_A: float
    wave_max_A: float
    span_status: str
    catalog_min_A: float
    catalog_max_A: float
    resolving_power: float
    provenance_ticket: str
    notes: str

    @property
    def abs_path(self) -> Path:
        q = Path(self.path)
        return q if q.is_absolute() else ROOT / q

    @property
    def reachable(self) -> bool:
        return self.abs_path.exists()

    def covers(self, wave_A: float) -> bool:
        return bool(self.wave_min_A <= float(wave_A) <= self.wave_max_A)


@dataclass(frozen=True)
class CoverageAnswer:
    """Deliberately not a bool. A bool is what produced the false NO-DATA verdict."""
    wave_A: float
    star: str
    covering: tuple[Instrument, ...]
    reachable_here: tuple[Instrument, ...]

    @property
    def covered(self) -> bool:
        return bool(self.covering)

    @property
    def is_data_gap(self) -> bool:
        """The ONLY state that justifies 'NO DATA' in an appendix or a ledger."""
        return not self.covering

    def why(self) -> str:
        if not self.covering:
            return (f"{self.wave_A:.3f} A is outside every instrument we hold data for "
                    f"on {self.star} — a real data gap")
        names = ", ".join(i.instrument for i in self.covering)
        if self.reachable_here:
            here = ", ".join(i.instrument for i in self.reachable_here)
            return f"{self.wave_A:.3f} A is covered by {names}; loadable here from {here}"
        return (f"{self.wave_A:.3f} A IS covered by {names}, but none of those files is "
                f"present on this machine — covered, not reachable. NOT a data gap.")


def load_registry(star: str | None = None) -> list[Instrument]:
    """Join the instrument catalog to the holdings registry. Neither file is owned here."""
    for q in (CATALOG, HOLDINGS):
        if not q.exists():
            raise CoverageError(
                f"{q} not found — it is the project's single source and this module "
                f"deliberately keeps no copy of it (RYA-708).")
    cat = {r["instrument_id"]: r for r in
           pd.read_csv(CATALOG, comment="#").to_dict("records")}
    hold = pd.read_csv(HOLDINGS, comment="#").to_dict("records")
    out: list[Instrument] = []
    for h in hold:
        if star is not None and str(h["system_id"]).strip() != star:
            continue
        man = ROOT / str(h["manifest_path"]).strip()
        if not man.exists():
            continue                      # holdings whose manifest is a per-target audit
        rows = pd.read_csv(man, comment="#")
        if "loader" not in rows.columns:
            continue                      # not a spectrum-location manifest
        iid = str(h["instrument_id"]).strip()
        for _, r in rows[rows.instrument_id == iid].iterrows():
            c = cat.get(iid)
            if c is None:
                raise CoverageError(
                    f"holdings names instrument_id {iid!r} which is not in the catalog. "
                    f"Register it there first — the catalog is the source for what an "
                    f"instrument IS.")
            host = str(r.host).strip().lower()
            if host not in VALID_HOSTS:
                raise CoverageError(f"{iid}: host {host!r} not in {sorted(VALID_HOSTS)}")
            out.append(Instrument(
                star=str(h["system_id"]).strip(), instrument=str(c["instrument_name"]).strip(),
                instrument_id=iid, host=host, path=str(r.path).strip(),
                loader=str(r.loader).strip(), wave_min_A=float(r.wave_min_A),
                wave_max_A=float(r.wave_max_A), span_status=str(r.span_status).strip().upper(),
                catalog_min_A=float(c["wavelength_min_nm"]) * 10.0,
                catalog_max_A=float(c["wavelength_max_nm"]) * 10.0,
                resolving_power=float(c["resolving_power_max"]),
                provenance_ticket=str(h["source_issue_ids"]).strip(),
                notes=str(r.get("notes") or "").strip()))
    if not out:
        raise CoverageError("no instrument holdings resolved — check the manifests")
    return out


def instruments_for(star: str, registry: list[Instrument] | None = None) -> list[Instrument]:
    reg = registry if registry is not None else load_registry(star)
    hits = [i for i in reg if i.star == star]
    if not hits:
        raise CoverageError(
            f"no instrument registered for star {star!r}. Register it before asking a "
            f"coverage question about it — an empty answer here would read as a data gap.")
    return hits


def coverage_at(wave_A: float, star: str = "solar",
                registry: list[Instrument] | None = None) -> CoverageAnswer:
    """THE question. Which instruments cover this wavelength for this star."""
    inst = instruments_for(star, registry)
    covering = tuple(i for i in inst if i.covers(wave_A))
    return CoverageAnswer(wave_A=float(wave_A), star=star, covering=covering,
                          reachable_here=tuple(i for i in covering if i.reachable))


def verify(star: str | None = None) -> int:
    """Measure the span of every reachable file and REPORT disagreement with the
    registry. Never rewrites — see the module docstring."""
    rows = load_registry(star)
    bad = 0
    print(f"{'star':7s}{'instrument':12s}{'host':8s}{'status':10s}"
          f"{'declared span':>24s}   measured")
    for i in rows:
        if not i.reachable:
            # Two different things, and conflating them hides a real defect: a file
            # absent from a machine it was never on is fine; a file absent from the
            # machine the registry SAYS it lives on is a broken row.
            import platform
            here = "sirius" if platform.node().lower().startswith("sirius") else "mac"
            elsewhere = i.host not in (here, "both")
            verdict = (f"on {i.host} only — not checked here, NOT a failure" if elsewhere
                       else f"MISSING at {i.path} though host={i.host} claims this machine")
            if not elsewhere:
                bad += 1
            print(f"{i.star:7s}{i.instrument[:26]:28s}{i.host:8s}{i.span_status:10s}"
                  f"{i.wave_min_A:10.1f}-{i.wave_max_A:.1f}   {verdict}")
            continue
        try:
            lo, hi, n = measure_span(i)
        except Exception as exc:                       # noqa: BLE001 — reported, not raised
            bad += 1
            print(f"{i.star:7s}{i.instrument[:26]:28s}{i.host:8s}{i.span_status:10s}"
                  f"{i.wave_min_A:10.1f}-{i.wave_max_A:.1f}   UNREADABLE: {exc}")
            continue
        agree = abs(lo - i.wave_min_A) < 1.0 and abs(hi - i.wave_max_A) < 1.0
        # The catalog and the holdings manifest are two independent statements about the
        # same instrument. They agreed here (kpno 296-1300 nm vs a measured 2960-13000 A),
        # but nothing was checking, and an unchecked agreement is a coincidence.
        cat_ok = (i.catalog_min_A - 1.0) <= lo and hi <= (i.catalog_max_A + 1.0)
        if not cat_ok:
            print(f"{'':7s}{'':28s}{'':8s}{'':10s}{'':24s}   note: measured span falls "
                  f"outside the CATALOG span {i.catalog_min_A:.0f}-{i.catalog_max_A:.0f} A")
        if not agree:
            bad += 1
        print(f"{i.star:7s}{i.instrument[:26]:28s}{i.host:8s}{i.span_status:10s}"
              f"{i.wave_min_A:10.1f}-{i.wave_max_A:.1f}   {lo:.1f}-{hi:.1f} "
              f"({n} pts) {'OK' if agree else '<-- DISAGREES with the registry'}")
    if bad:
        print(f"\n{bad} row(s) disagree with the registry or could not be read. The "
              f"registry is NOT auto-corrected: decide which is right and edit it "
              f"deliberately (RYA-708).")
    return 1 if bad else 0


def measure_span(inst: Instrument) -> tuple[float, float, int]:
    """Read a spectrum's actual air-wavelength span. Loaders are named in the registry
    so a new instrument is a row plus a loader, never an edit to every caller."""
    if inst.loader == "csv_normalized":
        d = pd.read_csv(inst.abs_path)
        col = next((c for c in d.columns if "wave" in c.lower()), None)
        if col is None:
            raise ValueError(f"no wavelength column in {inst.abs_path.name}")
        w = d[col].astype(float).values
    elif inst.loader == "iag_wavenumber_gz":
        # col0 is VACUUM WAVENUMBER in cm^-1, not a wavelength. Reading it as one gives
        # 9387-24700 and a confident wrong answer — I did exactly that before checking
        # the harness's own loader.
        import gzip
        import numpy as np
        from pipeline.wavelength_util import vac_to_air
        wn = []
        with gzip.open(inst.abs_path, "rt") as fh:
            for ln in fh:
                p = ln.split()
                if len(p) < 2:
                    continue
                try:
                    wn.append(float(p[0]))
                except ValueError:
                    continue
        w = vac_to_air(1e8 / np.asarray(wn))
    elif inst.loader == "kittpeak_lm_nm_dir":
        # 251 files named lmNNNN, NNNN = START WAVELENGTH IN NM. Column 0 is nm, not A.
        # The span is read from the file names and the first/last file's contents rather
        # than by loading 250 MB, because this runs on every --verify.
        import re
        import numpy as np
        d = inst.abs_path
        names = sorted(f for f in os.listdir(d) if re.match(r"^lm\d+$", f))
        if not names:
            raise ValueError(f"no lmNNNN files under {d}")
        def _edge(fname, last=False):
            vals = []
            for ln in open(d / fname):
                q = ln.split()
                if q:
                    try:
                        vals.append(float(q[0]) * 10.0)
                    except ValueError:
                        pass
            return (max(vals) if last else min(vals))
        w = np.array([_edge(names[0]), _edge(names[-1], last=True)])
        return float(w.min()), float(w.max()), len(names)
    else:
        raise ValueError(f"unknown loader {inst.loader!r} for {inst.instrument}")
    return float(w.min()), float(w.max()), int(len(w))


# ─────────────────────────────────────────────────────────────────────────────
# THE ENGINE LAYER — RYA-776
# ─────────────────────────────────────────────────────────────────────────────
#
# Everything above answers "WHICH INSTRUMENT sees this wavelength". That is half the
# question. The other half — "which ENGINE/GRID can model this wavelength" — was being
# re-derived from scratch every time it came up (RYA-763 existed to answer it for Fe I
# once; the answer then had no durable home). `element_status_tracker.csv` names WHICH
# grid and WHAT STATE per engine, but carries no WAVELENGTH REACH, so
# "do we have Engine A on Fe in the IR?" is not answerable from it.
#
# This half reads `data/catalog/engine_coverage.csv`, which is GENERATED from the grids
# and model atoms by scripts/generate_engine_coverage_rya776.py. It is never hand-edited:
# a hand-maintained coverage map that drifts is worse than no map at all, because a
# reader trusts it and gets more lost than they started.
#
# THE STATES IT REFUSES TO COLLAPSE
# ---------------------------------
# The instrument half above refuses to collapse "no instrument reaches this" into
# "this file does not have it". The engine half refuses the analogous collapse, which
# has bitten repeatedly as one blurry "no coverage":
#
#   * SERVED                    the per-line extract resolves the wavelength. Usable now.
#   * REACHABLE-NOT-EXTRACTED   the .grd / model atom LEVELS cover it, but the per-line
#                               extract does not expose it. The cheap-to-unlock class —
#                               a derivation away, not a data gap. This is what RYA-763
#                               found for Fe in the IR.
#   * UNCOVERED                 no atom/grid carries the levels at this wavelength.
#                               Genuinely absent. The ONLY state that justifies saying
#                               we cannot model here.
#
#   * REACH-UNKNOWN             (RYA-776, added on measurement — see below)
#
# WHY THERE IS A FOURTH STATE
# ---------------------------
# The three above presume reach is LOCALLY DECIDABLE — that some file on disk states
# which levels the engine carries. For a large part of the deck it is not:
#
#   * Fe's Engine A is the Bergemann/MPIA WEB SERVICE. There is no local Fe departure
#     grid to interrogate, so no local file can say whether MPIA reaches 8000 A. The
#     committed extract stops at 6843.7 A, but RYA-763 measured the live service still
#     answering at 46.7% in 6910-9199 A — so recording UNCOVERED there would be FALSE,
#     and it is exactly the false "no coverage" this table exists to end.
#   * Where the GES linelist itself carries no line in a band, there is nothing to
#     resolve, so a zero reach measures OUR CATALOGUE's span and not the grid's. The
#     9199.9 A wall is a linelist limit, not physics (atom.fe607a reaches 20000 A).
#
# Collapsing either case into UNCOVERED would manufacture the precise wrong answer this
# module was built to prevent, so it gets its own state and says which of the two it is.
# REACH-UNKNOWN is not a hedge: it is the difference between "we know there is nothing"
# and "we have not got a local file that could tell us".

ENGINE_COVERAGE = ROOT / "data" / "catalog" / "engine_coverage.csv"

SERVED = "SERVED"
REACHABLE_NOT_EXTRACTED = "REACHABLE-NOT-EXTRACTED"
UNCOVERED = "UNCOVERED"
REACH_UNKNOWN = "REACH-UNKNOWN"

VALID_ENGINE_STATES = (SERVED, REACHABLE_NOT_EXTRACTED, REACH_UNKNOWN, UNCOVERED)

# Strength order for reducing several grids to one answer for a species. SERVED beats
# everything; UNCOVERED is LAST on purpose — asserting "genuinely absent" requires every
# grid to have been decidable, so a single REACH-UNKNOWN outranks it and the species
# answer degrades to "not established" rather than to a false absence.
_STATE_RANK = {s: i for i, s in enumerate(VALID_ENGINE_STATES)}

VALID_ENGINES = frozenset({"A", "B"})


@dataclass(frozen=True)
class EngineCoverage:
    """One (element, ion, engine, grid, band) cell of the generated reference."""
    element: str
    ion: str
    engine: str
    grid_id: str
    band: str
    band_lo_A: float
    band_hi_A: float
    state: str
    n_lines_served: int
    n_lines_reachable: int
    n_lines_catalogued: int
    level_asset: str
    grid_asset: str
    note: str

    @property
    def species(self) -> str:
        return f"{self.element} {self.ion}"

    def covers(self, wave_A: float) -> bool:
        return bool(self.band_lo_A <= float(wave_A) < self.band_hi_A)


@dataclass(frozen=True)
class EngineReachAnswer:
    """Deliberately not a bare string, for the same reason CoverageAnswer is not a bool.

    An element can hold SEVERAL grids for one engine (Mg and Si each have both a
    Bergemann/MPIA and an Amarsi/PySME extract), and they do not have the same reach.
    `.state` reduces them for the common lookup; `.rows` keeps the grid-by-grid detail
    so the reduction is never the only thing on offer.
    """
    element: str
    ion: str
    engine: str
    wave_A: float
    rows: tuple[EngineCoverage, ...]

    @property
    def state(self) -> str:
        if not self.rows:
            return REACH_UNKNOWN
        return min((r.state for r in self.rows), key=lambda s: _STATE_RANK[s])

    @property
    def band(self) -> str:
        return self.rows[0].band if self.rows else ""

    @property
    def is_data_gap(self) -> bool:
        """True only when EVERY grid for this species/engine is decidably UNCOVERED."""
        return bool(self.rows) and all(r.state == UNCOVERED for r in self.rows)

    def why(self) -> str:
        head = (f"{self.element} {self.ion}, Engine {self.engine}, "
                f"{self.wave_A:.3f} A ({self.band or 'no declared band'})")
        if not self.rows:
            return (f"{head}: NO ROW in {ENGINE_COVERAGE.name} — this species/engine pair "
                    f"was never generated. That is an unbuilt reference, NOT a coverage "
                    f"verdict; regenerate before reading anything into it.")
        best = self.state
        det = "; ".join(f"{r.grid_id} [{r.state}]"
                        f"{f' {r.n_lines_served} served' if r.n_lines_served else ''}"
                        for r in self.rows)
        if best == SERVED:
            n = sum(r.n_lines_served for r in self.rows)
            return f"{head}: SERVED — {n} line(s) in band across {det}"
        if best == REACHABLE_NOT_EXTRACTED:
            n = sum(r.n_lines_reachable for r in self.rows)
            return (f"{head}: REACHABLE-NOT-EXTRACTED — the levels cover it "
                    f"({n} line(s) resolvable) but no per-line extract exposes it. "
                    f"A derivation away, NOT a data gap. {det}")
        if best == REACH_UNKNOWN:
            return (f"{head}: REACH-UNKNOWN — no local asset can decide this engine's "
                    f"reach here (service-only supplier, or no catalogued line in band). "
                    f"NOT a claim of absence. {det}")
        return (f"{head}: UNCOVERED — no atom or grid carries the levels here. "
                f"A real modelling gap. {det}")


def load_engine_coverage(path: Path | None = None) -> list[EngineCoverage]:
    """Read the generated engine x wavelength reference. Owns no data of its own."""
    p = path or ENGINE_COVERAGE
    if not p.exists():
        raise CoverageError(
            f"{p} not found. It is GENERATED — run "
            f"scripts/generate_engine_coverage_rya776.py on Sirius (the grids and model "
            f"atoms live there only). Never hand-write it: a stale coverage map is worse "
            f"than none, because it is trusted (RYA-776).")
    rows: list[EngineCoverage] = []
    for r in pd.read_csv(p, comment="#").to_dict("records"):
        state = str(r["state"]).strip().upper()
        if state not in _STATE_RANK:
            raise CoverageError(
                f"{p}: state {state!r} is not one of {list(VALID_ENGINE_STATES)}. An "
                f"unrecognised state must not be silently read as a coverage verdict.")
        engine = str(r["engine"]).strip().upper()
        if engine not in VALID_ENGINES:
            raise CoverageError(f"{p}: engine {engine!r} not in {sorted(VALID_ENGINES)}")
        rows.append(EngineCoverage(
            element=str(r["element"]).strip(), ion=str(r["ion"]).strip(),
            engine=engine, grid_id=str(r["grid_id"]).strip(),
            band=str(r["band"]).strip(), band_lo_A=float(r["band_lo_A"]),
            band_hi_A=float(r["band_hi_A"]), state=state,
            n_lines_served=int(r["n_lines_served"]),
            n_lines_reachable=int(r["n_lines_reachable"]),
            n_lines_catalogued=int(r["n_lines_catalogued"]),
            level_asset=str(r["level_asset"]).strip(),
            grid_asset=str(r["grid_asset"]).strip(),
            note=str(r.get("note") or "").strip()))
    if not rows:
        raise CoverageError(f"{p} has no rows — an empty reference is not 'no coverage'")
    return rows


def engine_reach(element: str, ion: str, engine: str, wave_A: float,
                 table: list[EngineCoverage] | None = None) -> EngineReachAnswer:
    """THE engine question: can `engine` model this species at this wavelength?

    A LOOKUP, not a re-derivation — that is the whole point of RYA-776. Returns the
    answer object; `.state` is the one-word reduction and `.why()` explains it.
    """
    eng = str(engine).strip().upper()
    if eng not in VALID_ENGINES:
        raise CoverageError(f"engine {engine!r} not in {sorted(VALID_ENGINES)}")
    tab = table if table is not None else load_engine_coverage()
    el, io = str(element).strip(), str(ion).strip()
    hits = tuple(r for r in tab
                 if r.element == el and r.ion == io and r.engine == eng
                 and r.covers(wave_A))
    return EngineReachAnswer(element=el, ion=io, engine=eng,
                             wave_A=float(wave_A), rows=hits)


def engine_summary(element: str, ion: str,
                   table: list[EngineCoverage] | None = None) -> str:
    """Compact per-species reach line for surfacing in the element status tracker.

    e.g. `A:VIS · B:VIS,red-optical?` — bands the engine SERVES, then bands it only
    REACHES marked `?`. Deliberately short: the tracker references this table, it does
    not absorb it.
    """
    tab = table if table is not None else load_engine_coverage()
    el, io = str(element).strip(), str(ion).strip()
    out = []
    for eng in sorted(VALID_ENGINES):
        rows = [r for r in tab if r.element == el and r.ion == io and r.engine == eng]
        if not rows:
            continue
        by_band: dict[str, str] = {}
        for r in sorted(rows, key=lambda r: r.band_lo_A):
            cur = by_band.get(r.band)
            if cur is None or _STATE_RANK[r.state] < _STATE_RANK[cur]:
                by_band[r.band] = r.state
        served = [b for b, s in by_band.items() if s == SERVED]
        reach = [b for b, s in by_band.items() if s == REACHABLE_NOT_EXTRACTED]
        bits = served + [f"{b}?" for b in reach]
        out.append(f"{eng}:{','.join(bits) if bits else 'none'}")
    return " · ".join(out) if out else "(no engine rows)"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--star", default="solar")
    ap.add_argument("--at", type=float, action="append",
                    help="ask about a wavelength in A (repeatable)")
    ap.add_argument("--verify", action="store_true",
                    help="measure every reachable file and report registry disagreement")
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--engine", choices=sorted(VALID_ENGINES),
                    help="ask the ENGINE question instead: does this engine reach --at? "
                         "Requires --element (and --ion).")
    ap.add_argument("--element")
    ap.add_argument("--ion", default="I")
    args = ap.parse_args(argv)

    if args.engine:
        if not (args.element and args.at):
            ap.error("--engine needs --element and at least one --at")
        tab = load_engine_coverage()
        for w in args.at:
            print(engine_reach(args.element, args.ion, args.engine, w, tab).why())
        return 0
    if args.verify:
        return verify(args.star)
    if args.list or not args.at:
        for i in instruments_for(args.star):
            mark = "loadable" if i.reachable else f"on {i.host} only"
            print(f"  {i.instrument[:26]:28s} {i.wave_min_A:8.1f}-{i.wave_max_A:.1f} A  "
                  f"{i.span_status:9s} R~{i.resolving_power:,.0f}  [{mark}]  {i.provenance_ticket}")
        return 0
    for w in args.at:
        print(coverage_at(w, args.star).why())
    return 0


if __name__ == "__main__":
    sys.exit(main())
