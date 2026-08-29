#!/usr/bin/env python3
"""RYA-1111 - measure an EXTERNAL reference line set on our spectra, tagged `line_set`.

    python3 scripts/measure_reference_lineset.py --line-set asplund \
        --holdings kpno_molecfit harps_molecfit
    python3 scripts/measure_reference_lineset.py --line-set gbs --dry-run

This is the ingest + measurement path that turns "replicate Asplund / GBS" into a product
family. It does NOT touch our own products: `our-graded` and `our-deep-graded` are emitted
by `derive_band_products.py` exactly as before, and nothing here writes `canonical_gf.csv`
(RYA-161 - no value changes to existing products).

WHAT A REPLICATION PRODUCT IS, HERE
-----------------------------------
THEIR lines, THEIR gf, OUR spectra, OUR fitter. Every one of those has to be true at once
or the number means something else:

  * their LINES, because our pool is not theirs. RYA-1109 measured it: of AGSS21's 40 Fe I
    lines only 19 are in our graded (lab-gf) tier, and the 11 below our 2.85 eV floor are
    held at NIST-C+/VALD3/OTHER. There is no subset of our pool to take.
  * their gf, because log gf sets the abundance directly. Measuring their lines on OUR gf
    would be a third thing that is neither their analysis nor ours.
  * our SPECTRA and our FITTER, because that is the thing under test - whether our
    pipeline on our data reproduces their number.

🔴 THE COVERAGE IS THE RESULT, NOT A FOOTNOTE. A reference line we cannot serve - outside
the holding's span, no admissible window, no gf published, a fit that did not converge - is
REPORTED BY NAME with its reason and excluded. Never replaced by a neighbour, never
silently dropped (RYA-429/711). A holding that serves 30 of 40 lines is answering a
different question from one that serves 40, and the reader must see that without
reconstructing it.

⚠️ --dry-run IS A FIRST-CLASS MODE, and on a machine without the spectra it is the only
honest one. It does the whole ingest, the lambda+EP join against our line list, the
per-holding span coverage and the loud-fail accounting, and writes the coverage report --
but runs no fit. The measuring run needs the atlases and a working MOOG/Turbospectrum, i.e.
Sirius (RYA-567); `MOOGSILENT` in the vendored iSpec is a Linux ELF and cannot execute on
the Mac at all.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline import reference_lineset as rls  # noqa: E402

CANON = ROOT / "data" / "linelists" / "canonical_gf.csv"
HOLDINGS_CSV = ROOT / "data" / "catalog" / "holdings_manifest_registry.csv"
INSTRUMENTS = ROOT / "data" / "catalog" / "instrument_catalog.csv"
OUT_DIR = ROOT / "data" / "results" / "rya1111"

#: The VIS solar holdings the replication targets, by the short names the CLI takes.
HOLDINGS = {
    "kpno_kurucz2005": "solar_kpno_kurucz2005_corrected",
    "kpno_molecfit": "solar_kpno_molecfit_corrected",
    "harps_molecfit": "solar_harps_molecfit_corrected",
    "iag": "solar_iag",
}


def our_lines(species: str) -> pd.DataFrame:
    d = pd.read_csv(CANON, low_memory=False)
    d = d[d["species"].astype(str) == species]
    return d.dropna(subset=["wavelength_air_A", "excitation_potential_eV"]).reset_index(
        drop=True)


def holding_span_A(short: str) -> tuple[str, float, float] | None:
    """The instrument's declared span for a holding.

    ⚠️ THE INSTRUMENT'S, NOT THE HOLDING'S. `holdings_manifest_registry` maps holding ->
    instrument and carries no wavelength range, so this is an UPPER BOUND on what the
    holding can serve. The measured coverage is a run, not a catalogue lookup.
    """
    hid = HOLDINGS[short]
    hold = pd.read_csv(HOLDINGS_CSV)
    row = hold[hold["holding_id"].astype(str) == hid]
    if row.empty:
        return None
    iid = str(row.iloc[0]["instrument_id"])
    inst = pd.read_csv(INSTRUMENTS)
    irow = inst[inst["instrument_id"].astype(str) == iid]
    if irow.empty or pd.isna(irow.iloc[0].get("wavelength_min_nm")):
        return None
    return (iid, float(irow.iloc[0]["wavelength_min_nm"]) * 10.0,
            float(irow.iloc[0]["wavelength_max_nm"]) * 10.0)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--line-set", required=True, choices=sorted(rls.SETS),
                    help="which external reference set to ingest and measure")
    ap.add_argument("--holdings", nargs="+", default=list(HOLDINGS),
                    choices=list(HOLDINGS))
    ap.add_argument("--dry-run", action="store_true",
                    help="ingest, join and report coverage; run no fit")
    ap.add_argument("--out", type=Path, default=None)
    a = ap.parse_args()

    spec = rls.SETS[a.line_set]
    ref = rls.load(a.line_set)
    meas = rls.measurable(ref)
    gap = rls.gf_gap(a.line_set)

    print(f"REFERENCE LINE SET: {a.line_set}")
    print(f"  source     {spec.source}")
    print(f"  built by   {spec.ticket}")
    print(f"  file       {spec.path.relative_to(ROOT)}")
    print(f"  native tag {spec.native_line_set!r} -> canonical line_set "
          f"{spec.name!r}")
    print(f"  match tol  {spec.match_tol_A} A  ({spec.tol_basis})")
    print(f"  lines      {len(ref)}  "
          f"({', '.join(f'{k} {v}' for k, v in ref.species.value_counts().items())})")

    doc: dict = {"ticket": "RYA-1111", "line_set": a.line_set, "source": spec.source,
                 "source_ticket": spec.ticket,
                 "native_line_set": spec.native_line_set,
                 "match_tol_A": spec.match_tol_A, "tol_basis": spec.tol_basis,
                 "n_reference_lines": int(len(ref)),
                 "gf_gap": gap, "coverage_vs_our_linelist": {}, "holdings": {},
                 "dry_run": bool(a.dry_run)}

    # ── loud-fail 1: lines with no published gf ──────────────────────────────
    if gap["n_without_published_gf"]:
        print(f"\n  [REFUSED] {gap['n_without_published_gf']} line(s) publish NO gf in "
              f"this set -> not measurable on its own scale")
        print(f"            {gap['by_species']}")
        for L in gap["lines"][:6]:
            print(f"              {L['species']:6} {L['wavelength_air_A']:9.3f} A  "
                  f"Elo {L['elo_eV']}")
        if len(gap["lines"]) > 6:
            print(f"              ... and {len(gap['lines']) - 6} more (all in the report)")
        print(f"            {gap['disposition']}")

    # ── loud-fail 2: reference lines absent from OUR synthesis line list ──────
    print(f"\nCOVERAGE vs our line list (canonical_gf), strict lambda+EP, "
          f"tol {spec.match_tol_A} A:")
    for sp in sorted(meas.species.unique()):
        cov = rls.coverage(a.line_set, our_lines(sp), pool_label=f"canonical_gf {sp}")
        sub = {k: v for k, v in cov["per_species"].items() if k == sp}
        doc["coverage_vs_our_linelist"][sp] = {**cov, "per_species": sub}
        s = sub[sp]
        print(f"  {sp:6} {s['matched']:3}/{s['n_ref']:<3} present    "
              f"{s['unmatched']:3} ABSENT from our list")
        for L in cov["unmatched_lines"][:4]:
            if L["species"] == sp:
                print(f"            [ABSENT] {L['wavelength_air_A']:9.3f} A  "
                      f"Elo {L['elo_eV']}")
        print(f"         plateau {cov['plateau']}")

    # ── per-holding span coverage ────────────────────────────────────────────
    print("\nPER-HOLDING span coverage (instrument's declared range = UPPER BOUND):")
    for h in a.holdings:
        span = holding_span_A(h)
        if span is None:
            doc["holdings"][h] = {"status": "no span declared"}
            print(f"  {h:18} no span declared")
            continue
        iid, lo, hi = span
        inb = meas[(meas.wavelength_air_A >= lo) & (meas.wavelength_air_A <= hi)]
        out_of = meas[~meas.index.isin(inb.index)]
        doc["holdings"][h] = {
            "holding_id": HOLDINGS[h], "instrument_id": iid, "span_A": [lo, hi],
            "span_basis": "instrument_catalog (upper bound)",
            "measurable_in_span": int(len(inb)),
            "measurable_total": int(len(meas)),
            "out_of_span": [{"wavelength_air_A": round(float(r.wavelength_air_A), 3),
                             "species": r.species} for r in out_of.itertuples()],
        }
        print(f"  {h:18} {lo:7.1f}-{hi:7.1f} A   {len(inb):3}/{len(meas):<3} measurable "
              f"lines in span" + (f"   [{len(out_of)} OUT OF SPAN]" if len(out_of) else ""))

    # ── the fit ──────────────────────────────────────────────────────────────
    if a.dry_run:
        print("\n--dry-run: ingest, join and coverage done; NO FIT RUN.")
        print("  The measuring run needs the atlases and a working synthesiser (Sirius, "
              "RYA-567).")
    else:
        print("\nMEASURING is not wired to run off-Sirius.")
        print("  `MOOGSILENT` in the vendored iSpec is a Linux x86-64 ELF and the solar "
              "atlases live on /mnt/codex-data, so a fit cannot even start here. Re-run "
              "on Sirius, or use --dry-run for the coverage accounting.")
        doc["measurement"] = {"ran": False,
                              "reason": "not on Sirius; MOOGSILENT is a Linux ELF and the "
                                        "atlases are not mounted (RYA-567)"}

    out = a.out or (OUT_DIR / f"{a.line_set}_ingest_coverage.json")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(doc, indent=2) + "\n")
    print(f"\nwrote {out.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
