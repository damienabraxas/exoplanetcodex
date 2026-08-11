#!/usr/bin/env python3
"""RYA-760 — adjudicate the gf sources against VALD, from data already in the repo.

    python3 scripts/rya760_ges_vs_vald.py

WHY NOT NIST
------------
The ticket says the FMW offset is "checkable against NIST ASD directly". Two problems with
that, and Ryan flagged the second:

  1. ASD is currently returning a server-side error for every query (lines1.pl:701),
     including the project's own previously-recorded working recipe. Externally blocked.
  2. **It would have been close to circular anyway.** FMW is Fuhr, Martin & Wiese — itself
     a NIST compilation. Asking NIST whether NIST's numbers are right is not an independent
     test; it can only catch a transcription error.

And it is unnecessary. `canonical_gf.csv` already carries, per line, the gf from MORE THAN
ONE catalogue: `gf_synth_ges` (29300 Fe rows) and `gf_linelist_vald` (8034 Fe rows). Where
both exist, the difference is an independent cross-catalogue check on the SAME transition,
with no network and a sample two orders of magnitude larger than the 12 measured FMW lines.

THE TEST
--------
Group by `loggf_reference` and compare `gf_synth_ges - gf_linelist_vald`.

A source with a genuine scale offset should show a systematic, non-zero difference against
VALD, while the laboratory sources (BWL / Ruffoni / Bard-Kock) should sit near zero. That
is the same claim the ticket wanted to put to NIST, asked of a catalogue that is actually
independent of it.

Caveat kept in view: VALD is itself a compilation and inherits its own sources, so a
non-zero difference localises a DISAGREEMENT between catalogues rather than proving which
one is wrong. What it can do — and what matters here — is say whether FMW is anomalous
*relative to the sources we already trust*, which is exactly the adjudication asked for.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

CANON = ROOT / "data" / "linelists" / "canonical_gf.csv"
SOURCE_DEX = {
    "BWL": 0.02, "2014MNRAS.": 0.03, "GESHRL14": 0.03, "BK": 0.04, "BKK": 0.04,
    "FMW": 0.08, "GESB82c": 0.08, "MRW": 0.08, "K07": 0.20, "K13": 0.20, "K75": 0.20,
}


def base(ref) -> str:
    return str(ref).strip().split("+")[0]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--element", default="Fe")
    ap.add_argument("--ion", type=int, default=1)
    ap.add_argument("--min-n", type=int, default=5)
    a = ap.parse_args()

    c = pd.read_csv(CANON, low_memory=False)
    c["el"] = c.species.astype(str).str.split().str[0]
    d = c[(c.el == a.element) & (c.ion == a.ion)].copy()
    d["base"] = d.loggf_reference.map(base)
    d["acc"] = d.base.map(SOURCE_DEX)

    both = d[d.gf_synth_ges.notna() & d.gf_linelist_vald.notna()].copy()
    both["delta"] = both.gf_synth_ges - both.gf_linelist_vald
    print(f"\n{a.element} {'I'*a.ion}: {len(d)} canonical rows, "
          f"{len(both)} with BOTH a GES and a VALD log gf")

    print("\n" + "=" * 78)
    print("GES - VALD log gf, by source   (a scale offset shows up here, no network)")
    print("=" * 78)
    print(f"{'source':<12}{'acc':>6}{'n':>6}{'median d':>11}{'MAD':>9}"
          f"{'|d|>0.05':>10}{'|d|>0.10':>10}")
    rows = []
    for b, g in both.groupby("base"):
        if len(g) < a.min_n:
            continue
        med = g.delta.median()
        mad = (g.delta - med).abs().median()
        rows.append(dict(source=b, acc=g.acc.iloc[0], n=len(g), median=med, mad=mad,
                         f05=float((g.delta.abs() > 0.05).mean()),
                         f10=float((g.delta.abs() > 0.10).mean())))
    r = pd.DataFrame(rows).sort_values(["acc", "source"], na_position="last")
    for _, x in r.iterrows():
        acc = f"{x.acc:.2f}" if np.isfinite(x.acc) else "  ?"
        print(f"{x.source:<12}{acc:>6}{int(x.n):>6}{x['median']:>11.4f}{x.mad:>9.4f}"
              f"{x.f05:>10.3f}{x.f10:>10.3f}")

    lab = r[r.acc <= 0.04]
    loose = r[(r.acc > 0.04) & (r.acc <= 0.08)]
    print("\n" + "-" * 78)
    if len(lab) and len(loose):
        lm = np.average(lab["median"], weights=lab.n)
        gm = np.average(loose["median"], weights=loose.n)
        print(f"  laboratory sources (acc <= 0.04), n-weighted median d = {lm:+.4f}")
        print(f"  loose sources      (0.04-0.08)  , n-weighted median d = {gm:+.4f}")
        print(f"  separation                                            = {gm - lm:+.4f} dex")
        print("\n  A gf too LOW makes the derived abundance too HIGH, so a NEGATIVE")
        print("  delta on the loose sources is what the +0.29 dex A offset predicts.")
    print("\n  NOTE: VALD is itself a compilation, so a non-zero delta localises a")
    print("  DISAGREEMENT between catalogues, not which one is wrong. What it settles")
    print("  is whether a source is anomalous relative to the ones we already trust.")


if __name__ == "__main__":
    main()
