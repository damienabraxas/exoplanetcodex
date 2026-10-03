#!/usr/bin/env python3
"""RYA-1233 -- pull NIST ASD Si I / Si II gf (governing process step 6, "secure correct
VALD/NIST lines"), for the solar Si orchestrator test.

WHY. canonical_gf carries 2,473 Si I/II rows and exactly ONE NIST grade (Si I 5793.073,
B); 2,050 Si I rows are Kurucz and 165 VALD3. So process step 7 (graded lines) holds
every Si band but VIS, which has 3 graded lines.

WHAT THE SOURCE IS. NIST ASD's Si I/II transition probabilities are the critical
compilation of Kelleher & Podobedova 2008, "Atomic Transition Probabilities of Silicon.
A Critical Compilation", J. Phys. Chem. Ref. Data 37, 1285 (doi:10.1063/1.2734566). For
Si I it rests on Garz 1973 renormalised to the O'Brian & Lawler 1991 lifetimes -- the
same gf base Amarsi & Asplund 2017 / Asplund 2021 use (RYA-725). These are critically
EVALUATED values: the honest tier is the NIST graded family (NIST-C+ at C+ or better),
never LAB.

This script only ACQUIRES and records (a primary_gf holding + provenance). The
adjudication into canonical_gf is scripts/rya1233_adjudicate_si_nist.py.
The pull mechanics are RYA-1160's `safe_pull` and RYA-822's `tidy`, reused, not copied.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "linelists" / "primary_gf"
SPECIES = ("Si I", "Si II")
SOURCE = ("Kelleher & Podobedova 2008, J. Phys. Chem. Ref. Data 37, 1285 "
          "(doi:10.1063/1.2734566) -- the NIST ASD Si critical compilation")


def _load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--species", action="append", default=None)
    ap.add_argument("--lo-A", type=float, default=3000.0)
    ap.add_argument("--hi-A", type=float, default=25000.0)
    ap.add_argument("--step-A", type=float, default=2000.0)
    ap.add_argument("--pause-s", type=float, default=0.4)
    a = ap.parse_args()
    r1160 = _load(ROOT / "scripts" / "rya1160_pull_nist_cno.py", "rya1160_pull")
    if a.hi_A > r1160.CEILING_A:
        raise SystemExit(f"refusing: --hi-A beyond the {r1160.CEILING_A:.0f} A ceiling")
    rya822 = r1160.load_rya822()
    import astroquery
    OUT.mkdir(parents=True, exist_ok=True)
    for sp in a.species or list(SPECIES):
        print(f"\n=== {sp}  {a.lo_A:.0f}-{a.hi_A:.0f} A ===")
        raw = r1160.safe_pull(rya822, a.lo_A, a.hi_A, a.step_A, sp, a.pause_s)
        tid = rya822.tidy(raw)
        dest = OUT / f"nist_asd_{sp.replace(' ', '')}_{int(a.lo_A)}_{int(a.hi_A)}.tsv"
        tid.to_csv(dest, sep="\t", index=False)
        prov = {
            "ticket": "RYA-1233",
            "source": "NIST Atomic Spectra Database (ASD), lines",
            "compilation": SOURCE,
            "access": f"astroquery.nist {astroquery.__version__} (Nist.query)",
            "reused_from": "scripts/rya1160_pull_nist_cno.py safe_pull + "
                           "scripts/rya822_pull_nist_nearuv.py tidy",
            "species": sp, "band_A": [a.lo_A, a.hi_A], "chunk_A": a.step_A,
            "wavelength_type": "vac+air (air above 2000 A)",
            "pulled_at_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "n_rows": int(len(tid)), "n_with_log_gf": int(tid.log_gf.notna().sum()),
            "n_graded": int(tid.nist_grade.notna().sum()),
            "grade_counts": {str(k): int(v) for k, v in tid.nist_grade.value_counts().items()},
        }
        dest.with_suffix(".prov.json").write_text(json.dumps(prov, indent=2) + "\n")
        print(f"  -> {dest.relative_to(ROOT)}: {prov['n_rows']} rows, "
              f"{prov['n_graded']} graded {prov['grade_counts']}")


if __name__ == "__main__":
    main()
