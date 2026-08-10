#!/usr/bin/env python3
"""Ingest per-line gf SOURCES from a VALD long-format extract — RYA-713.

    python3 scripts/ingest_vald_references.py \
        --raw ../exoplanetcodex-rya558/data/linelists/vald_solar_nearuv_2000_3780_hfson_raw.txt \
        --element Fe --lo 3000 --hi 3800

Ryan: *"what about the VALD references in the IR? Not everything needs to be NIST graded"*
→ *"re-ingest the near-UV VALD with references."*

WHY THIS EXISTS
---------------
Every near-UV Fe line in `linelist_solar.csv` carries `loggf_source = VALD3` — a bare tag
with no source behind it. All 5135 of them. That is not a statement that the oscillator
strengths are poor; it is a statement that **we did not record where they came from**, so
no uncertainty can be assigned to any of them.

The information was never missing. VALD's long-format extract carries the source inline,
as a fourth line per record:

    'Ti 2',  2000.00849, -2.199, 4.0088, 2.5, ...
    '  LS                                3d3 d2D1'
    '  LS                                3d2.(3F).5d 4D*'
    '_   Kurucz TiII 2016   1 wl:K16   1 gf:K16   1 K16 ...   Ti+'
                                          ^^^^^^ the gf source

This reads that tag and writes it per line, which is exactly what made the IR product
possible: in the red-optical the same information was already carried through
`canonical_gf` (which spans 3780.3–9199.9 Å and therefore **starts above the near-UV**),
and splitting on it moved the IR abundance 0.106 dex without touching a measurement.

WHAT IT DOES NOT DO
-------------------
It does not assign an accuracy. It records the SOURCE verbatim. Mapping a source to an
uncertainty is a separate, citable judgement (`GF_SOURCE_ACCURACY` below) and is kept
separate so a wrong accuracy estimate can be corrected without re-parsing 68 MB.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

OUT = ROOT / "data" / "linelists"

# Published accuracies for the sources VALD cites, in dex. Laboratory measurements first.
# These are ESTIMATES OF THE SOURCE, not of any individual line, and each is citable.
GF_SOURCE_ACCURACY = {
    # --- laboratory ---
    "BWL":   (0.02, "Blackwell/Wells/Lynas-Gray — Oxford Fe I furnace absorption"),
    "BKK":   (0.04, "Bard, Kock & Kock 1991 — laboratory emission"),
    "BK":    (0.04, "Bard & Kock 1994 — laboratory emission"),
    "RHL":   (0.03, "Ruffoni/Heiter/Lawler — FTS laboratory"),
    "RU":    (0.03, "Ruffoni et al. — FTS laboratory"),
    "MRW":   (0.08, "May, Richter & Wichelmann 1974 — laboratory"),
    "FMW":   (0.08, "Fuhr, Martin & Wiese — NIST critical compilation"),
    "BGHL":  (0.05, "Bard/Gurell/Hartman/Lawler-family FTS"),
    "OBRIAN": (0.03, "O'Brian et al. 1991 — Fe I laboratory"),
    # --- semi-empirical / theoretical ---
    "K":     (0.20, "Kurucz semi-empirical (K07/K14/K16/K88…) — 0.1–0.3 dex, RYA-161"),
}

_REC = re.compile(r"^'([A-Z][a-z]?)\s+(\d+)',\s*([\d.]+),\s*(-?[\d.]+),")
_GF = re.compile(r"gf:(\S+)")


def source_accuracy(tag: str) -> tuple[float, str]:
    """Map a VALD gf tag to (dex, citation). Kurucz codes are matched last, as a family."""
    t = str(tag).strip().upper()
    for k in sorted((k for k in GF_SOURCE_ACCURACY if k != "K"), key=len, reverse=True):
        if t.startswith(k):
            return GF_SOURCE_ACCURACY[k]
    if re.match(r"^K\d", t):            # K07, K14, K16, K88 ...
        return GF_SOURCE_ACCURACY["K"]
    return (float("nan"), f"unrecognised VALD source tag {tag!r} — accuracy NOT assumed")


def parse(raw: Path, element: str, lo: float, hi: float) -> pd.DataFrame:
    """Walk the extract. A record is a data line plus up to three continuation lines; the
    gf tag lives on the one beginning with `'_`."""
    rows, pending = [], None
    with raw.open(errors="replace") as f:
        for line in f:
            m = _REC.match(line)
            if m:
                el, ion, wl, lgf = m.group(1), int(m.group(2)), float(m.group(3)), float(m.group(4))
                pending = dict(element=el, ion=ion, wavelength_air_A=wl, log_gf=lgf) \
                    if (el == element and lo <= wl < hi) else None
                continue
            if pending is not None and line.lstrip().startswith("'_"):
                g = _GF.search(line)
                pending["vald_gf_source"] = g.group(1) if g else ""
                rows.append(pending)
                pending = None
    return pd.DataFrame(rows)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--raw", required=True)
    ap.add_argument("--element", required=True)
    ap.add_argument("--lo", type=float, required=True)
    ap.add_argument("--hi", type=float, required=True)
    ap.add_argument("--out", default=str(OUT))
    a = ap.parse_args()

    raw = Path(a.raw)
    if not raw.exists():
        raise SystemExit(f"VALD extract not found: {raw}")
    print(f"parsing {raw.name} ({raw.stat().st_size/1e6:.1f} MB) for {a.element} "
          f"in {a.lo:.0f}-{a.hi:.0f} A ...")
    df = parse(raw, a.element, a.lo, a.hi)
    if not len(df):
        raise SystemExit("no matching records — check the element symbol and range")

    df[["gf_dex", "gf_citation"]] = df.vald_gf_source.apply(
        lambda t: pd.Series(source_accuracy(t)))

    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    dest = out / f"vald_gf_sources_{a.element}_{int(a.lo)}_{int(a.hi)}.csv"
    df.to_csv(dest, index=False)

    print(f"\n  {len(df)} {a.element} lines with a recorded gf source\n")
    print(f"{'source':>14s}{'n':>7s}{'dex':>8s}  citation")
    for tag, n in df.vald_gf_source.value_counts().head(15).items():
        d, cite = source_accuracy(tag)
        print(f"{tag or '(none)':>14s}{n:7d}{d:8.3f}  {cite[:58]}")
    lab = df[df.gf_dex <= 0.08]
    print(f"\n  laboratory / compilation (<=0.08 dex): {len(lab)} of {len(df)} "
          f"({100*len(lab)/len(df):.1f}%)")
    print(f"  semi-empirical (>0.10 dex)           : {int((df.gf_dex > 0.10).sum())}")
    print(f"  unrecognised, accuracy NOT assumed   : {int(df.gf_dex.isna().sum())}")
    print(f"\n  wrote {dest}")


if __name__ == "__main__":
    main()
