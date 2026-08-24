#!/usr/bin/env python3
"""RYA-1030: determine `normalization_state` per holding FROM THE FLUX, and record it.

Sibling of `rya806_backfill_telluric_applied.py`, which does the same job for the
telluric axis. Run with no flags to report; `--write` to update the registry.

🔴 IT MEASURES WHAT THE READER RETURNS, NOT WHAT A FILE GLOB FINDS. This is the whole
point and it is the one design decision worth defending. RYA-806's backfill resolves a
holding to its files and reads their headers, which is right for a question the headers
answer. Normalisation is not that question: the failure this axis exists to catch is a
MIS-ROUTE -- a holding declared normalised whose reader opens the raw file (KP2005 before
RYA-933; `iag_fts_solar_atlas` on the telluric axis, still open). Globbing the staged
directory would find the normalised file sitting there unused and cheerfully report
`normalised`, reproducing the exact blindness that let KP2005 run mis-routed for months.
So every probe goes through `measure_band_ew.load_window_ex` -- the same call the
measurement harness makes, on the same holding -- and classifies what actually comes back.

MULTIPLE PROBE WINDOWS, AND DISAGREEMENT IS AN ANSWER. One window in a line-crowded
region can hold no true continuum and drag the envelope low (KP1984 col1 reads 0.98, not
1.00). Several windows spread across the holding's declared span are classified
independently: unanimity is the verdict, and a SPLIT is reported `unknown` with the split
shown, because one flag cannot honestly describe a product that answers two ways. That
mirrors RYA-806's MIXED handling rather than inventing a second convention.
"""
from __future__ import annotations

import argparse
import csv
import io
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from pipeline.normalization_intake import (  # noqa: E402
    NORMALISED, UNKNOWN, UN_NORMALISED, detect)

HOLDINGS = ROOT / "data" / "catalog" / "holdings_manifest_registry.csv"
COLUMN = "normalization_state"

#: Half-width of each probe, in Angstrom. Wide enough that the envelope statistic has
#: many windows to work with (the detector needs >= 5 x 2 A), narrow enough to load fast.
PROBE_PAD_A = 25.0

#: Fractions across a holding's declared span to probe at. Spread deliberately: a blaze
#: or SED shows itself BETWEEN windows as much as within one, and three points that agree
#: on the level while disagreeing on nothing else is a much stronger statement than one.
PROBE_FRACTIONS = (0.25, 0.50, 0.75)


def _instrument_span_A(instrument: str) -> tuple[float, float] | None:
    """Declared reach of the arm, from `instrument_catalog.csv`, in Angstrom.

    The fallback when a HoldingSpec declares no `span_A` of its own -- which most do not
    (only HARPS and KP2005 carry one, RYA-767). Taken from the catalogue rather than
    guessed: a hard-coded "somewhere in the visible" would silently probe outside a
    holding's reach and report `unknown` for a product that is perfectly readable, which
    is manufacturing an absence.
    """
    import pandas as pd
    cat = pd.read_csv(ROOT / "data" / "catalog" / "instrument_catalog.csv")
    row = cat[cat.instrument_id.astype(str) == str(instrument)]
    if not len(row):
        return None
    lo = pd.to_numeric(row.iloc[0].get("wavelength_min_nm"), errors="coerce")
    hi = pd.to_numeric(row.iloc[0].get("wavelength_max_nm"), errors="coerce")
    if pd.isna(lo) or pd.isna(hi):
        return None
    return float(lo) * 10.0, float(hi) * 10.0


def _probe_centres(spec, instrument: str) -> list[float]:
    span = spec.span_A or _instrument_span_A(instrument)
    if not span:
        return []
    lo, hi = float(span[0]) + PROBE_PAD_A, float(span[1]) - PROBE_PAD_A
    if hi <= lo:
        return []
    return [lo + frac * (hi - lo) for frac in PROBE_FRACTIONS]


def determine(harness, holding_id: str, instrument: str) -> tuple[str, str]:
    """(value, evidence). Probes through the HARNESS's reader; see the module docstring."""
    specs = {h.holding_id: h
             for hs in harness._INSTRUMENT_HOLDINGS.values() for h in hs}
    spec = specs.get(holding_id)
    if spec is None:
        return UNKNOWN, (
            f"RYA-1030: no HoldingSpec -- `{holding_id}` is registered but the measurement "
            f"harness has no reader for it, so there is nothing to scan. This is not "
            f"'un-normalised'; it is unmeasured, and the two must not collapse (RYA-833).")

    centres = _probe_centres(spec, instrument)
    if not centres:
        return UNKNOWN, (
            f"RYA-1030: neither the HoldingSpec nor `instrument_catalog.csv` declares a "
            f"usable wavelength span for `{instrument}`, so no probe window could be "
            f"placed. Declare `span_A` (RYA-767) or the catalogue range, then re-run -- "
            f"do not guess a window, because probing outside a holding's reach reports "
            f"`unknown` for a product that reads perfectly (RYA-833).")

    verdicts, notes, failures = [], [], []
    for centre in centres:
        try:
            win = harness.load_window_ex(instrument, centre, PROBE_PAD_A,
                                         holding=holding_id, allow_uncorrected=True)
        except Exception as exc:                       # unreachable data, refused gate
            failures.append(f"{centre:.0f} A: {type(exc).__name__}: {exc}")
            continue
        ev = detect(win.wave, win.flux)
        verdicts.append(ev.value)
        notes.append(f"{centre:.0f} A -> {ev.value} "
                     f"(P95={ev.median_p95:.4f}, slope={ev.slope_across_band:+.4f})"
                     if ev.median_p95 is not None else f"{centre:.0f} A -> {ev.value}")

    if not verdicts:
        return UNKNOWN, (
            f"RYA-1030: no probe window could be READ, so the flux never spoke. "
            f"{'; '.join(failures)}. Product not reachable from this machine, or the "
            f"holding is gated -- either way `unknown` is the honest record, and a "
            f"consumer refuses it rather than assuming (RYA-833).")

    served = ""
    distinct = set(verdicts)
    if len(distinct) > 1:
        return UNKNOWN, (
            f"RYA-1030 INTAKE (through the harness reader{served}): ⚠️ SPLIT across probe "
            f"windows {dict(Counter(verdicts))} => the holding is `unknown`: one flag "
            f"cannot honestly describe a product that answers two ways. "
            f"{'; '.join(notes)}")

    value = verdicts[0]
    note = (f"RYA-1030 INTAKE (through the harness reader, {len(verdicts)} probe window"
            f"{'s' if len(verdicts) != 1 else ''}): {'; '.join(notes)}")
    if failures:
        note += f" | unread: {'; '.join(failures)}"
    return value, note


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--write", action="store_true", help="write the registry")
    ap.add_argument("--only", help="limit to one holding_id (for a quick check)")
    a = ap.parse_args(argv)

    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "_mbe_rya1030", ROOT / "scripts" / "measure_band_ew.py")
    harness = importlib.util.module_from_spec(spec)
    sys.modules["_mbe_rya1030"] = harness
    spec.loader.exec_module(harness)

    raw = HOLDINGS.read_bytes()
    crlf = b"\r\n" in raw
    rows = list(csv.reader(io.StringIO(raw.decode("utf-8"), newline="")))
    hdr = rows[0]
    if COLUMN not in hdr:
        hdr.append(COLUMN)
        for r in rows[1:]:
            if r:
                r.append("")
    i_id, i_inst = hdr.index("holding_id"), hdr.index("instrument_id")
    i_val, i_notes = hdr.index(COLUMN), hdr.index("notes")

    print(f"{'holding_id':<34} {'value':<14} evidence")
    print("-" * 118)
    changed = 0
    for r in rows[1:]:
        if not r:
            continue
        hid, inst = r[i_id], r[i_inst]
        if a.only and hid != a.only:
            continue
        value, note = determine(harness, hid, inst)
        print(f"{hid:<34} {value:<14} {note[:150]}")
        if r[i_val] != value:
            changed += 1
        r[i_val] = value
        if note not in r[i_notes]:
            r[i_notes] = (r[i_notes] + " " + note).strip()

    print(f"\n{changed} row(s) would change.")
    if not a.write:
        print("Dry run. Pass --write to record.")
        return 0

    out = io.StringIO(newline="")
    csv.writer(out, lineterminator="\r\n" if crlf else "\n").writerows(rows)
    HOLDINGS.write_text(out.getvalue(), newline="")
    print(f"WROTE {HOLDINGS}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
