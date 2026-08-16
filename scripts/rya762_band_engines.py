#!/usr/bin/env python3
"""RYA-762 — which ENGINES can serve Fe I in 9199.9-13000 A, per line.

    python3 scripts/rya762_band_engines.py            # the table
    python3 scripts/rya762_band_engines.py --scan-tol # measure the tolerance, don't pick it

THE TICKET'S PREMISE WAS BACKWARDS, AND THIS IS THE MEASUREMENT THAT SHOWS IT
-----------------------------------------------------------------------------
RYA-762 was written with this product table:

    1D-LTE    ✓        Engine A  ✓ to 11973 A        Engine B  ✗ no level IDs past 9199.9

The Engine-B cell is wrong, and it was wrong for a reason worth keeping: **9199.9 A is
where the GES linelist ENDS**, not where the physics stops. The Gerber model atom
`atom.fe607a` is cut at lambda < 20000 A and carries 548 Fe I levels, so the deck already
reaches this band. What was missing was never the atom — it was a way to identify a
line's levels without GES's level-index column. RYA-763 built exactly that: identify by
**(J, energy)**, which VALD carries natively for both endpoints.

So Engine B is available here, and by a wide margin:

    Engine A (MPIA, RYA-763's live probe)   0.183
    Engine B (atom.fe607a, this script)     0.813

⚠️ THE RATES ARE NOT COMPARABLE POPULATIONS. RYA-763 measured 0.276 for Engine B in
6910-9199 A against the GES list, which carries a long tail of weak lines whose upper
levels the 607-level atom does not have. This band's candidates come from a VALD
extraction cut at depth 0.05, so they are the strong, low-lying-level lines the atom is
most likely to carry. A higher rate here is a property of the LINE SELECTION, not
evidence that the atom serves the IR better than the optical. Quote it as "0.813 of the
239 catalogued at depth>=0.05", never as a bare "81%".

⚠️ AND IT DOES NOT MAKE ENGINE B "PRIMARY". Ryan, 2026-08-11: *"there is no Primary. All
Engines, LTE and NLTE, are products that get presented."* A wider reach makes an engine
BROADER, not better. Each engine is its own product keyed (instrument x band x engine)
per RYA-712; the cross-engine delta is a FINDING to report, and it stops being readable
the moment one engine is treated as the answer the others correct.

WHAT THIS DOES NOT DO
---------------------
It does not wavelength-match into the grid. Every line here is level-identified on both
endpoints or it is REFUSED — which is the distinction RYA-762's "do NOT extend Engine B"
instruction was actually protecting (a wavelength match without level identification is
the silent-LTE-fallback failure RYA-713 exists to prevent). The instruction's letter said
9199.9 A; its intent was "never guess a level". This obeys the intent.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from data.linelists.vald_parse import parse_vald_long                 # noqa: E402
from pipeline import band_policy                                      # noqa: E402
from scripts.rya763_level_mapping import (                            # noqa: E402
    GERBER_DIR, read_gerber_atom, resolve_by_label,
)

ACCOUNTING = ROOT / "data" / "audit" / "line_accounting" / "per_line.csv"
VALD_RAW = ROOT / "data" / "linelists" / "vald_solar_ir_6910_14000_raw.txt"
OUT = ROOT / "data" / "audit" / "rya762" / "band_engines.csv"

LO_A, HI_A = 9199.9, 13000.0
TOL_EV = 0.001          # RYA-763's MEASURED optimum; --scan-tol re-measures it here
MATCH_TOL_A = 0.01      # accounting <-> VALD wavelength join


def candidates() -> pd.DataFrame:
    """The band's Fe I candidates, from the committed accounting pool.

    The pool is the SSOT for "what lines exist and which instruments see them" — it
    already reaches 24985 A, so this band needed no new inventory, only the engine
    question answering. Depth and instrument filters match the band harness's defaults
    so the two agree on what a candidate is.
    """
    acc = pd.read_csv(ACCOUNTING)
    sel = acc[(acc.element == "Fe") & (acc.ion == "I") &
              (acc.wave_air_A > LO_A) & (acc.wave_air_A <= HI_A) &
              acc.predicted_depth.between(0.05, 0.60) &
              acc.instruments.notna()].copy()
    return sel.sort_values("wave_air_A").reset_index(drop=True)


def vald_levels() -> pd.DataFrame:
    """Fe I level endpoints from the VALD extract — the (J, energy) key, natively.

    This is why the band is reachable at all: GES stops at 9200 A and would have had to
    supply a level INDEX anyway, which is deck-native and hazardous across decks
    (RYA-763). VALD supplies the physical coordinate instead.
    """
    recs, report = parse_vald_long(str(VALD_RAW))
    if report["n_failures"]:
        raise SystemExit(f"VALD parse dropped {report['n_failures']} line(s) — "
                         f"no-silent-drop (RYA-429); fix before trusting this table")
    rows = [r for r in recs
            if r["element"] == "Fe" and r["ion"] == "I" and LO_A < r["wavelength"] <= HI_A]
    return pd.DataFrame(rows)


def resolve(cands: pd.DataFrame, vald: pd.DataFrame, lab: pd.DataFrame,
            tol_eV: float) -> pd.DataFrame:
    out = []
    for _, c in cands.iterrows():
        w = float(c.wave_air_A)
        near = vald[(vald.wavelength - w).abs() <= MATCH_TOL_A]
        if near.empty:
            # Stated, not dropped: the accounting pool and this VALD extract are two
            # different holdings and need not agree line-for-line.
            out.append(dict(wave_air_A=w, engine_b="NO-VALD-LEVELS",
                            lower="", upper="", band=band_policy.resolve(w).name))
            continue
        v = near.iloc[(near.wavelength - w).abs().argmin()]
        vlo, _, _ = resolve_by_label(lab, float(v.e_low_eV), float(v.j_low), tol_eV)
        vup, _, _ = resolve_by_label(lab, float(v.e_up_eV), float(v.j_up), tol_eV)
        out.append(dict(
            wave_air_A=w, engine_b=("AVAILABLE" if vlo == "UNIQUE" and vup == "UNIQUE"
                                    else "REFUSED"),
            lower=vlo, upper=vup, band=band_policy.resolve(w).name))
    df = pd.DataFrame(out)
    return df.merge(cands[["wave_air_A", "log_gf", "ep_eV", "predicted_depth",
                           "instruments"]], on="wave_air_A", how="left")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--scan-tol", action="store_true",
                    help="measure resolution rate vs tolerance instead of assuming it")
    a = ap.parse_args(argv)

    cands = candidates()
    vald = vald_levels()
    atom = sorted(GERBER_DIR.glob("atom.fe*"))
    if not atom:
        raise SystemExit(f"no Gerber Fe atom under {GERBER_DIR} — this must run on Sirius")
    lab_all = read_gerber_atom(atom[0])
    lab = lab_all[lab_all["ion"].astype(int) == 1]

    print(f"Fe I  {LO_A:.1f}-{HI_A:.0f} A")
    print(f"  candidates (accounting pool, depth 0.05-0.60, instrument-covered): {len(cands)}")
    print(f"  VALD Fe I records in band: {len(vald)}")
    print(f"  {atom[0].name}: {len(lab)} Fe I levels")

    # The band STRADDLES a policy boundary, and that decides the method per line --
    # this is not cosmetic: profile-fit is FORBIDDEN above 10000 A.
    print("\n  band policy split (RYA-713 intake rules):")
    for p in sorted(band_policy.POLICIES, key=lambda p: p.lo_A):
        n = int(((cands.wave_air_A >= p.lo_A) & (cands.wave_air_A < p.hi_A)).sum())
        if n:
            print(f"    {p.name:12s} {n:4d} lines   permitted={p.permitted_methods}"
                  f"{'  telluric REQUIRED' if p.telluric_required else ''}")

    if a.scan_tol:
        print("\n  Engine-B resolution vs tolerance (MEASURED, not picked):")
        print(f"    {'tol (eV)':>10} {'available':>10} {'rate':>8}")
        for t in (0.0005, 0.001, 0.002, 0.005, 0.01, 0.05):
            r = resolve(cands, vald, lab, t)
            k = int((r.engine_b == "AVAILABLE").sum())
            print(f"    {t:10.4f} {k:10d} {k/max(len(r),1):8.3f}")

    df = resolve(cands, vald, lab, TOL_EV)
    k = int((df.engine_b == "AVAILABLE").sum())
    print(f"\n  ENGINE B at tol={TOL_EV} eV: {k} of {len(df)}  ({k/max(len(df),1):.3f})")
    for st, n in df.engine_b.value_counts().items():
        print(f"    {st:<16} {n}")
    print("\n  per policy band:")
    for b, g in df.groupby("band"):
        kk = int((g.engine_b == "AVAILABLE").sum())
        print(f"    {b:12s} {kk:4d} of {len(g):4d} Engine-B available")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    df.sort_values("wave_air_A").to_csv(OUT, index=False)
    print(f"\n  wrote {OUT.relative_to(ROOT)}")
    print("\n  NOTE: Engine A is NOT re-derived here — RYA-763 measured it live at 0.183\n"
          "  in this band, and the committed Fe_Bergemann_MPIA.csv extract has ZERO nodes\n"
          "  past 6910 A, so any Engine-A count must name which source it used.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
