#!/usr/bin/env python3
"""RYA-1079: per-line telluric observability census over a real band.

Classifies every graded line in a band CLEAN / RECOVERABLE / SATURATED from the MEASURED
telluric depth (see `pipeline.telluric_observability`), reports the census, and emits the
`problem_children` rows for the saturated class.

    python3 scripts/rya1079_observability_census.py --instrument kpno_solar_atlas \\
        --holding solar_kpno --band red-optical --element Fe --ion I

Emits nothing to the registry unless `--write-problem-children` is given: the verdict is
reported first and adopted second (RYA-807).
"""
from __future__ import annotations

import argparse
import csv
import importlib.util
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from config import synth_bands                                       # noqa: E402
from pipeline import telluric_observability as obs                   # noqa: E402
from pipeline.gf_empirical import GRADED_TIERS                       # noqa: E402
from pipeline.problem_children import SCHEMA_COLUMNS                 # noqa: E402

OUT_DIR = ROOT / "data" / "audit" / "rya1079_observability"


def harness():
    spec = importlib.util.spec_from_file_location(
        "_mbe_rya1079", ROOT / "scripts" / "measure_band_ew.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules["_mbe_rya1079"] = mod
    spec.loader.exec_module(mod)
    return mod


def graded_lines(band: str, element: str | None, ion: str | None) -> list[float]:
    import pandas as pd
    b = synth_bands.SYNTH_BANDS[band]
    df = pd.read_csv(ROOT / "data/linelists/canonical_gf.csv", low_memory=False)
    df = df[df.gf_tier.astype(str).isin(GRADED_TIERS)]
    if element:
        df = df[df.species.astype(str).str.split().str[0] == element]
    if ion:
        df = df[df.species.astype(str).str.split().str[-1] == ion]
    w = pd.to_numeric(df.wavelength_air_A, errors="coerce").dropna()
    return sorted(float(x) for x in w if b.lo_A <= x <= b.hi_A)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Per-line telluric observability (RYA-1079).")
    ap.add_argument("--instrument", required=True)
    ap.add_argument("--holding", required=True)
    ap.add_argument("--band", required=True, choices=sorted(synth_bands.SYNTH_BANDS))
    ap.add_argument("--element", default=None)
    ap.add_argument("--ion", default=None)
    ap.add_argument("--write-problem-children", action="store_true",
                    help="append the SATURATED rows to data/registry/problem_children.csv")
    a = ap.parse_args(argv)

    band = synth_bands.SYNTH_BANDS[a.band]
    lines = graded_lines(a.band, a.element, a.ion)
    print(f"{a.holding} / {a.band} ({band.lo_A:.0f}-{band.hi_A:.0f} A): "
          f"{len(lines)} graded lines, half-width {band.half_width_A} A")

    p = obs.partition_pool(harness(), a.instrument, a.holding, lines,
                           band.half_width_A, (band.lo_A, band.hi_A))
    print(f"\nOBSERVABILITY CENSUS (measured, before any tier is assigned)\n  {p.summary()}")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    stem = f"{a.holding}_{a.band}"
    with (OUT_DIR / f"{stem}_observability.csv").open("w", newline="",
                                                      encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["wavelength_air_A", "verdict", "disposition", "depth",
                    "transmission_min", "clean_max_depth", "saturated_min_depth",
                    "snr_continuum", "corrected", "reason"])
        for o in p.measurable + p.needs_correction + p.unmeasured:
            w.writerow([f"{o.wavelength_A:.4f}", o.verdict, o.disposition,
                        f"{o.depth:.6f}", f"{o.transmission_min:.6f}",
                        f"{o.thresholds.clean_max_depth:.6f}",
                        f"{o.thresholds.saturated_min_depth:.6f}",
                        f"{o.thresholds.snr_continuum:.1f}", o.corrected, o.reason])
        for wl, why in p.uncharacterised:
            w.writerow([f"{wl:.4f}", "UNCHARACTERISED", "REFUSE", "", "", "", "", "",
                        "", why[:200]])
    print(f"  -> {(OUT_DIR / f'{stem}_observability.csv').relative_to(ROOT)}")

    species = f"{a.element} {a.ion}" if a.element and a.ion else "unknown"
    rows = [obs.problem_child_row(o, species) for o in p.unmeasured]
    if rows:
        dest = OUT_DIR / f"{stem}_problem_children.csv"
        with dest.open("w", newline="", encoding="utf-8") as fh:
            wr = csv.DictWriter(fh, fieldnames=SCHEMA_COLUMNS)
            wr.writeheader()
            wr.writerows(rows)
        print(f"  -> {dest.relative_to(ROOT)}  ({len(rows)} saturated-core exclusions)")

    if p.needs_correction:
        print(f"\nRESCUED CLASS -- {len(p.needs_correction)} line(s) that the RYA-460/786 "
              f"band rule would have SKIPPED.\nThey are RECOVERABLE: route this holding "
              f"through the RYA-424 correction stage and measure them.")
        for o in p.needs_correction[:5]:
            print(f"    {o.wavelength_A:10.3f}  depth={o.depth:.4f}  T_min="
                  f"{o.transmission_min:.4f}  -> {o.disposition}")
    if p.unmeasured:
        print(f"\nGENUINELY LOST -- {len(p.unmeasured)} line(s) under a saturated core:")
        for o in p.unmeasured[:5]:
            print(f"    {o.wavelength_A:10.3f}  depth={o.depth:.4f}  T_min="
                  f"{o.transmission_min:.4f}  reason={o.reason}")

    if a.write_problem_children and rows:
        reg = ROOT / "data" / "registry" / "problem_children.csv"
        raw = reg.read_bytes()
        crlf = b"\r\n" in raw
        with reg.open("a", newline="", encoding="utf-8") as fh:
            wr = csv.DictWriter(fh, fieldnames=SCHEMA_COLUMNS,
                                lineterminator="\r\n" if crlf else "\n")
            wr.writerows(rows)
        print(f"\nWROTE {len(rows)} row(s) to {reg.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
