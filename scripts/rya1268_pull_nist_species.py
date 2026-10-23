#!/usr/bin/env python3
"""RYA-1268 -- pull NIST ASD lines for every remaining canonical species (the lab-gf half of
Reference Grade: "Asplund set PLUS any lab-referenced line").

WHY. canonical_gf carries NIST grades / lab tags for Fe, Si, Al and CNO only; every other
TARGET_ELEMENTS species has ZERO graded rows (RYA-1268 coverage matrix). The Asplund-cited
tables (Scott+2015a/b, Grevesse+2015) print no per-line gf sigma, so a published set cannot
be graded without a per-line error -- NIST's accuracy class supplies it where the transition
is in the critical compilation.

This script only ACQUIRES and records (a primary_gf holding + provenance), exactly as
scripts/rya1233_pull_nist_si.py does for Si. Adjudication into canonical_gf is a separate,
reviewed step. The pull mechanics are RYA-1160's `safe_pull` and RYA-822's `tidy`, reused.
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
SPECIES = ("Li I", "Na I", "P I", "S I", "K I", "Ca I", "Ca II", "Sc I", "Sc II",
           "Ti I", "Ti II", "V I", "V II", "Cr I", "Cr II", "Mn I", "Mn II", "Co I", "Co II",
           "Ni I", "Ni II", "Cu I", "Zn I", "Sr I", "Sr II", "Y I", "Y II", "Zr I", "Zr II",
           "Ba I", "Ba II", "Eu II")


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
    ap.add_argument("--skip-existing", action="store_true")
    a = ap.parse_args()
    r1160 = _load(ROOT / "scripts" / "rya1160_pull_nist_cno.py", "rya1160_pull")
    if a.hi_A > r1160.CEILING_A:
        raise SystemExit(f"refusing: --hi-A beyond the {r1160.CEILING_A:.0f} A ceiling")
    rya822 = r1160.load_rya822()
    import astroquery
    OUT.mkdir(parents=True, exist_ok=True)
    for sp in a.species or list(SPECIES):
        # El_Ion, never ElIon: on a case-insensitive filesystem (macOS) "NiI" IS RYA-1160's
        # "NII" (N II) -- --skip-existing skipped Ni I, and without it the pull overwrote N II.
        dest = OUT / f"nist_asd_{sp.replace(' ', '_')}_{int(a.lo_A)}_{int(a.hi_A)}.tsv"
        if a.skip_existing and dest.exists():
            continue
        print(f"\n=== {sp}  {a.lo_A:.0f}-{a.hi_A:.0f} A ===", flush=True)
        try:
            raw = r1160.safe_pull(rya822, a.lo_A, a.hi_A, a.step_A, sp, a.pause_s)
        except Exception as e:  # one species failing must not lose the rest
            print(f"  FAILED {sp}: {e!r}", flush=True)
            continue
        tid = rya822.tidy(raw)
        tid.to_csv(dest, sep="\t", index=False)
        prov = {
            "ticket": "RYA-1268",
            "source": "NIST Atomic Spectra Database (ASD), lines",
            "compilation": "NIST ASD; per-line TP/Line reference codes retained (tp_ref)",
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
              f"{prov['n_graded']} graded {prov['grade_counts']}", flush=True)


if __name__ == "__main__":
    main()
