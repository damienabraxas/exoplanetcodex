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


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--star", default="solar")
    ap.add_argument("--at", type=float, action="append",
                    help="ask about a wavelength in A (repeatable)")
    ap.add_argument("--verify", action="store_true",
                    help="measure every reachable file and report registry disagreement")
    ap.add_argument("--list", action="store_true")
    args = ap.parse_args(argv)

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
