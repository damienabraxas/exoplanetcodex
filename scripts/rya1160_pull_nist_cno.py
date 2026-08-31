#!/usr/bin/env python3
"""
RYA-1160 — pull NIST ASD C/N/O gf: the first graded gf source CNO has ever had.

WHY THIS EXISTS.  Before this pull, `canonical_gf` carried 2,606 CNO atomic rows and
essentially no grade on any of them: `gf_tier` LAB = 0 (the repo holds 453 LAB rows and
CNO owns none), `lab_source_tag` blank on all 2,606, `gf_sigma_dex` blank on all 2,606,
`gf_source_doi` blank on all 2,606, and exactly five rows with a NIST grade.  RYA-1129
recorded the same finding per element and left the gate closed on it.

WHAT THE SOURCE ACTUALLY IS, from the compilers' own prose rather than from memory.
`Reference documents/20060052476.pdf` (Wiese & Fuhr, NASA LAW 2006) states the chain:

  * the base compilation for all three elements is Wiese, Fuhr & Deters 1996, "Atomic
    Transition Probabilities of Carbon, Nitrogen, and Oxygen: A Critical Data
    Compilation", JPCRD Monograph 7, which was "primarily based on the very extensive
    calculational results of the OPACITY Project";
  * Wiese & Fuhr then published a PARTIAL update, and its scope is C I, C II, N I and
    N II ONLY, driven by new MCHF calculations (Froese Fischer & Tachiev 2004) which
    showed "the OPACITY data are often not as accurate as we had estimated ... This
    statement applies especially to neutral and singly-ionized carbon and nitrogen."

🔴 OXYGEN WAS NOT UPDATED.  O I and O II still rest on the 1996 Opacity Project values
while C and N carry the newer MCHF ones.  The two are not the same vintage or the same
method, and nothing downstream should treat "NIST grade" as one homogeneous authority
across CNO.  The per-line `ref_transition_probability` (TP) code is carried for exactly
this reason -- it names the underlying source per line.

🔴 THIS IS NOT A LABORATORY MEASUREMENT AND MUST NOT BE TIERED `LAB`.  For C, N and O
these are critically evaluated, largely THEORETICAL values (Opacity Project close
coupling, then MCHF).  That is the accepted standard for the light elements, where
theory is most accurate and lab measurement hardest -- but the honest tier is the
existing `NIST-C+` / graded family, never `LAB`.  RYA-1005's lesson on Al applies
verbatim: a truncated bibcode was once read as "lab" and built an
experimental-vs-semi-empirical comparison out of two theory sources.

This script only ACQUIRES and records.  It writes a primary_gf holding plus provenance
and touches no gf, no tier and no `canonical_gf` row.  Adjudication is a separate step.

Usage:
    python3 scripts/rya1160_pull_nist_cno.py             # all six species
    python3 scripts/rya1160_pull_nist_cno.py --species "O I"
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import sys
import time
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "linelists" / "primary_gf"

#: The pull mechanics are RYA-822's and are deliberately REUSED, not reimplemented:
#: `wavelength_type='vac+air'` (astroquery defaults to vacuum, and our lists are air
#: above 2000 A -- a ~1 A silent offset), log gf = log10(gi * fik) with NaN kept honest,
#: Ei already in eV, and the TP/Line reference codes preserved.
RYA822 = ROOT / "scripts" / "rya822_pull_nist_nearuv.py"

#: 25000 A IS A HARD CEILING: the Codex holds no stellar spectra beyond it, so a line
#: above 25000 A can never be measured here no matter how good its atomic data. Do not
#: pull past it. 3000 A is the matching floor -- `canonical_gf` holds zero rows below it
#: for any element (RYA-1158), though NIST does have graded CNO down to 900 A.
LO_A, HI_A = 3000.0, 25000.0
CEILING_A = 25000.0

SPECIES = ("C I", "C II", "N I", "N II", "O I", "O II")

#: Scope of the C/N update, per the compilers. Recorded PER SPECIES so the vintage
#: asymmetry survives into the artifact instead of living only in this docstring.
UPDATE_SCOPE = {
    "C I": "Wiese & Fuhr partial update (MCHF); base Wiese, Fuhr & Deters 1996 Monograph 7",
    "C II": "Wiese & Fuhr partial update (MCHF); base Wiese, Fuhr & Deters 1996 Monograph 7",
    "N I": "Wiese & Fuhr partial update (MCHF); base Wiese, Fuhr & Deters 1996 Monograph 7",
    "N II": "Wiese & Fuhr partial update (MCHF); base Wiese, Fuhr & Deters 1996 Monograph 7",
    "O I": "NOT UPDATED -- Wiese, Fuhr & Deters 1996 Monograph 7 (Opacity Project) only",
    "O II": "NOT UPDATED -- Wiese, Fuhr & Deters 1996 Monograph 7 (Opacity Project) only",
}


def safe_pull(rya822, lo, hi, step, species, pause):
    """RYA-822's pull, with one fix that only shows up over a wide band.

    astroquery returns masked columns whose dtype depends on what a chunk happens to
    contain, so a column that is Int64 in one chunk and float64 in the next makes
    pd.concat raise `cannot safely cast non-equivalent float64 to int64`. It surfaced
    only above 25000 A. Every column is cast to object BEFORE concat; tidy() re-parses
    with pd.to_numeric anyway, so nothing is lost and no value is coerced.
    """
    import pandas as pd, numpy as np, time
    from astroquery.nist import Nist
    import astropy.units as u

    frames = []
    edges = np.arange(lo, hi, step)
    for i, a in enumerate(edges):
        b = min(a + step, hi)
        t = None
        for attempt in range(3):
            try:
                t = Nist.query(a * u.AA, b * u.AA, linename=species,
                               wavelength_type="vac+air")
                break
            except Exception as e:
                if attempt == 2:
                    print(f"  {a:.0f}-{b:.0f} FAILED: {type(e).__name__}")
                    t = None
                else:
                    time.sleep(2 + 3 * attempt)
        if t is None or len(t) == 0:
            print(f"  {a:9.1f}-{b:9.1f} A   0 rows")
            continue
        df = t.to_pandas().astype(object)
        df["chunk_lo_A"], df["chunk_hi_A"] = a, b
        frames.append(df)
        print(f"  {a:9.1f}-{b:9.1f} A {len(df):5d} rows ({i+1}/{len(edges)})", flush=True)
        time.sleep(pause)
    if not frames:
        raise SystemExit("NIST returned nothing across the whole band — refusing to "
                         "write an empty pull that would read as 'no data exists'")
    return pd.concat(frames, ignore_index=True)


def load_rya822():
    spec = importlib.util.spec_from_file_location("rya822", RYA822)
    module = importlib.util.module_from_spec(spec)
    sys.modules["rya822"] = module
    spec.loader.exec_module(module)
    return module


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--species", action="append", default=None)
    ap.add_argument("--lo-A", type=float, default=LO_A)
    ap.add_argument("--hi-A", type=float, default=HI_A)
    ap.add_argument("--step-A", type=float, default=2000.0)
    ap.add_argument("--pause-s", type=float, default=0.4)
    a = ap.parse_args()
    if a.hi_A > CEILING_A:
        raise SystemExit(f"refusing: --hi-A {a.hi_A:.0f} exceeds the {CEILING_A:.0f} A "
                         "instrument ceiling; the Codex has no spectra beyond it")
    species = a.species or list(SPECIES)

    rya822 = load_rya822()
    import astroquery

    OUT.mkdir(parents=True, exist_ok=True)
    summary = {}
    for sp in species:
        print(f"\n=== {sp}  {a.lo_A:.0f}-{a.hi_A:.0f} A ===")
        raw = safe_pull(rya822, a.lo_A, a.hi_A, a.step_A, sp, a.pause_s)
        tid = rya822.tidy(raw)
        slug = sp.replace(" ", "")
        dest = OUT / f"nist_asd_{slug}_{int(a.lo_A)}_{int(a.hi_A)}.tsv"
        tid.to_csv(dest, sep="\t", index=False)

        graded = int(tid.nist_grade.notna().sum())
        with_gf = int(tid.log_gf.notna().sum())
        prov = {
            "ticket": "RYA-1160",
            "source": "NIST Atomic Spectra Database (ASD), lines",
            "access": f"astroquery.nist {astroquery.__version__} (Nist.query)",
            "reused_from": "scripts/rya822_pull_nist_nearuv.py (pull + tidy)",
            "species": sp,
            "band_A": [a.lo_A, a.hi_A],
            "chunk_A": a.step_A,
            "underlying_compilation": UPDATE_SCOPE[sp],
            "compilation_chain": (
                "Wiese, Fuhr & Deters 1996, JPCRD Monograph 7 (base, Opacity Project); "
                "Wiese & Fuhr partial update for C I/C II/N I/N II only (MCHF, Froese "
                "Fischer & Tachiev 2004). Chain read from the compilers' own text in "
                "Reference documents/20060052476.pdf, not from memory."),
            "not_laboratory": (
                "These are CRITICALLY EVALUATED, LARGELY THEORETICAL values (Opacity "
                "Project close coupling, then MCHF) -- the accepted standard for C/N/O, "
                "but NOT laboratory measurements. Must never be tiered LAB. The "
                "per-line TP code names the underlying source."),
            "wavelength_frame": "vac+air (air above 2000 A, matching canonical_gf)",
            "n_rows_raw": int(len(raw)),
            "n_rows_tidy": int(len(tid)),
            "n_with_log_gf": with_gf,
            "n_graded": graded,
            "grade_counts": {k: int(v) for k, v in tid.nist_grade.value_counts().items()},
            "adjudicated": False,
            "canonical_gf_modified": False,
            "pulled_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        }
        (OUT / f"nist_asd_{slug}_{int(a.lo_A)}_{int(a.hi_A)}.prov.json").write_text(
            json.dumps(prov, indent=2) + "\n")
        summary[sp] = {"rows": len(tid), "with_log_gf": with_gf, "graded": graded}
        print(f"  rows {len(tid)}   log gf {with_gf}   graded {graded}")
        print(f"  [out] {dest.relative_to(ROOT)}")

    print("\n=== RYA-1160 summary ===")
    tot = 0
    for sp, s in summary.items():
        print(f"  {sp:6s} rows {s['rows']:5d}   log gf {s['with_log_gf']:5d}   graded {s['graded']:5d}")
        tot += s["graded"]
    print(f"  TOTAL graded CNO transitions acquired: {tot}")
    print("  canonical_gf NOT modified; adjudication is a separate step.")


if __name__ == "__main__":
    main()
