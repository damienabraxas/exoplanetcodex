"""RYA-806 — determine `telluric_applied` for every holding and write it to the registry.

This is the INTAKE step made repeatable. It does not accept an argued value: where the
product files are reachable it reads their headers via `pipeline.telluric_intake`, and
where they are not it takes a value only if a prior ticket MEASURED it, carrying that
citation into the registry. Anything else stays `unknown`, which consumers refuse.

    python3 scripts/rya806_backfill_telluric_applied.py            # report only
    python3 scripts/rya806_backfill_telluric_applied.py --write    # write the registry

⚠️ The registry is edited with `csv`, never round-tripped through pandas: pandas rewrites
quoting and line endings across every untouched row and buries the real change (RYA-786).
"""
from __future__ import annotations

import argparse
import csv
import glob
import io
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from pipeline.telluric_intake import APPLIED, NOT_APPLIED, UNKNOWN, from_many  # noqa: E402
from config.constants import codex_path  # RYA-810 path register

HOLDINGS = ROOT / "data" / "catalog" / "holdings_manifest_registry.csv"
DATA = codex_path('data.spectra_local')
ACEN = DATA / "Alpha Centauri (vetted)"

# holding_id -> glob(s) of the product files. Determined from headers when they resolve.
FILE_SOURCES: dict[str, tuple[str, ...]] = {
    "alpha_cen_a_crires_plus": (str(ACEN / "Alpha Cen A" / "CRIRES" / "*.fits"),),
    "alpha_cen_b_crires_plus": (str(ACEN / "Alpha Cen B" / "CRIRES" / "*.fits"),),
    "alpha_cen_b_nirps": (str(ACEN / "Alpha Cen B" / "NIRPS" / "*.fits"),
                          str(ACEN / "Alpha Cen A" / "NIRPS" / "*.fits")),
    "alpha_cen_a_harps": (str(ACEN / "Alpha Cen A" / "HARPS" / "*.fits"),),
    "alpha_cen_b_harps": (str(ACEN / "Alpha Cen B" / "HARPS" / "*.fits"),),
    "alpha_cen_b_espresso": (str(ACEN / "Alpha Cen B" / "ESPRESSO" / "*.fits"),),
    "procyon_uves": (str(DATA / "Procyon" / "Procyon UVES" / "**" / "*.fits"),),
    "solar_harps": (str(DATA / "Solar Calibration" / "**" / "*HARPS*" / "**" / "*.fits"),),
}

# Holdings whose telluric state was MEASURED by a prior ticket. The citation is the
# evidence; these are not opinions and each names the number that settled it.
CITED: dict[str, tuple[str, str]] = {
    "solar_vesta_crires_plus_idp": (NOT_APPLIED, (
        "RYA-805 (headers, all 18): PRO CATG=OBS_NODDING_EXTRACTC_IDP, the PRO REC chain "
        "is exactly ONE recipe cr2res_obs_nodding with no REC2, no transmission HDU, no "
        "FLUX_TELL column, zero genuine telluric keywords. Confirmed in the flux: 2.02-"
        "2.81% of px below 0.7 continuum in the O2 1.27um band vs 0.54% for the "
        "telluric-CORRECTED Elgueta control, and H-band H2O depth tracks header IWV x "
        "airmass at Pearson r=+0.996")),
    "solar_crires_plus_y_rya794": (APPLIED, (
        "RYA-794/805: this is a DERIVED product built from Elgueta+2026 sp/Sun_Y_rv.dat, "
        "not from the raw IDPs. Elgueta's reduced spectra are telluric-corrected - "
        "measured at 0.10% of window points below 0.5 in the O2 A-band vs 51.3% for Kitt "
        "Peak. Inherits `applied` from its source, NOT from the crires_plus instrument "
        "axis")),
    "elgueta2026_vizier": (APPLIED, (
        "RYA-789/794: Elgueta+2026 published REDUCED spectra (sp/Sun_*_rv.dat, and the "
        "GBS star spectra); telluric-corrected by their pipeline - 0.10% of window points "
        "below 0.5 in the O2 A-band vs 51.3% for Kitt Peak in the same window")),
    "solar_kpno": (NOT_APPLIED, (
        "RYA-783/786 MEASURED: the Kurucz 1984 Kitt Peak atlas RETAINS its tellurics - "
        "51.3% of O2 A-band core points below 0.5 and flux driven to zero inside the "
        "bands. It runs anyway on telluric_basis=line_selection (per-line clean-line "
        "selection, RYA-460/786), which is a stated method, not uncorrected data "
        "sneaking through")),
    "solar_iag": (APPLIED, (
        "RYA-783/786 MEASURED: the IAG FTS atlas IS telluric-corrected - 0.1% below 0.5 "
        "in the O2 A-band core and 0.0% in H2O 9280-9600, while showing the SAME solar "
        "line depths as Kitt Peak in clean regions, so the correction removed tellurics "
        "without erasing solar structure. Matches telluric_basis=corrected")),
    "alpha_cen_a_stis": (NOT_APPLIED, "space-based UV"),
    "alpha_cen_b_stis": (NOT_APPLIED, "space-based UV"),
    "procyon_stis": (NOT_APPLIED, "space-based UV"),
    "procyon_cos": (NOT_APPLIED, "space-based UV"),
}
_SPACE_UV = (
    "HST is above the atmosphere, so no telluric correction was applied and none is "
    "possible or needed. Recorded as not-applied because that is literally true; it does "
    "NOT gate, because the instrument axis carries telluric_basis=not_applicable. This "
    "is exactly why the two axes are kept separate (RYA-806)")


def resolve(holding_id: str):
    files: list[str] = []
    for pat in FILE_SOURCES.get(holding_id, ()):
        files += glob.glob(pat, recursive=True)
    return sorted(set(files))


def determine(holding_id: str) -> tuple[str, str, int]:
    """(value, evidence, n_files_read)."""
    files = resolve(holding_id)
    if files:
        value, evs = from_many(files)
        e = evs[0]
        note = (f"RYA-806 INTAKE (read {len(files)} file"
                f"{'s' if len(files) != 1 else ''} from headers): {e.citation()}")
        if len({x.value for x in evs}) > 1:
            from collections import Counter
            note += (f" ⚠️ MIXED across files {dict(Counter(x.value for x in evs))} => the "
                     f"holding is unknown: one flag cannot honestly describe it")
        return value, note, len(files)
    if holding_id in CITED:
        value, cite = CITED[holding_id]
        if cite == "space-based UV":
            cite = _SPACE_UV
        return value, f"RYA-806: {cite}", 0
    return UNKNOWN, (
        "RYA-806: product files were not reachable at backfill time and no prior ticket "
        "measured this holding's telluric state, so it stays `unknown` and consumers "
        "refuse it. Determine with pipeline.telluric_intake.from_headers() and record."), 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--write", action="store_true", help="write the registry")
    a = ap.parse_args(argv)

    raw = HOLDINGS.read_bytes()
    crlf = b"\r\n" in raw
    rows = list(csv.reader(io.StringIO(raw.decode("utf-8"), newline="")))
    hdr = rows[0]
    if "telluric_applied" not in hdr:
        hdr.append("telluric_applied")
        for r in rows[1:]:
            if r:
                r.append("")
    i_id, i_notes = hdr.index("holding_id"), hdr.index("notes")
    i_app, i_src = hdr.index("telluric_applied"), hdr.index("source_issue_ids")

    print(f"{'holding_id':<30} {'value':<12} {'files':>6}  evidence")
    print("-" * 110)
    changed = 0
    for r in rows[1:]:
        if not r:
            continue
        hid = r[i_id]
        value, note, n = determine(hid)
        print(f"{hid:<30} {value:<12} {n:>6}  {note[:70]}")
        if r[i_app] != value:
            r[i_app] = value
            changed += 1
        if "RYA-806" not in r[i_notes]:
            r[i_notes] = (r[i_notes].rstrip() + " " + note).strip()
        if "RYA-806" not in r[i_src]:
            r[i_src] = (r[i_src] + ";RYA-806") if r[i_src] else "RYA-806"

    from collections import Counter
    tally = Counter(r[i_app] for r in rows[1:] if r)
    print("-" * 110)
    print(f"tally: {dict(tally)}   rows changed: {changed}")

    if a.write:
        buf = io.StringIO(newline="")
        csv.writer(buf, lineterminator="\r\n" if crlf else "\n",
                   quoting=csv.QUOTE_MINIMAL).writerows(rows)
        HOLDINGS.write_bytes(buf.getvalue().encode("utf-8"))
        print(f"wrote {HOLDINGS} ({'CRLF' if crlf else 'LF'} preserved)")
    else:
        print("(report only; pass --write to update the registry)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
