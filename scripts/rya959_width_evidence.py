#!/usr/bin/env python3
"""RYA-959 — the evidence behind the width ceiling, derived rather than typed.

`data/measured/band_ew/` is gitignored (RYA-469): the pools are regenerable from the
committed reference table plus the committed harness, so what gets committed is the
EVIDENCE — the verdict census, the fitted-sigma histogram that shows the ceiling lands in
a valley, the implied-width distribution that shows no impossible width survives, and the
fate of the laboratory-graded lines.

Every number quoted in `docs/science/rya959_fe_vis_remeasure.md` and in the Linear report
comes out of here. Re-run it and the doc can be checked against the artifact:

    python3 scripts/rya959_width_evidence.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from pipeline.line_width import (  # noqa: E402
    physical_ceiling_sigma_A, voigt_fwhm, max_stellar_sigma_kms, irreducible_sigma_kms)

EW_DIR = ROOT / "data" / "measured" / "band_ew"
OUT = ROOT / "data" / "audit" / "rya959_width_ceiling"

#: The two arms of the RYA-959 VIS run. HARPS is the science holding (RYA-931
#: telluric-corrected); Kitt Peak is the control, and it is a CONTROL precisely because it
#: is a different instrument on the same Sun — an agreement between them is evidence, an
#: agreement of one arm with itself is not.
ARMS = {
    "harps_molecfit_corrected":
        "FeI_3780_6910_harps_solar_harps_molecfit_corrected_PROFILEFIT_ew.csv",
    "kpno_solar_atlas":
        "FeI_3780_6910_kpno_solar_atlas_PROFILEFIT_ew.csv",
}

#: Resolving power per arm, read from the catalogue rather than typed (RYA-713).
CATALOG = ROOT / "data" / "catalog" / "instrument_catalog.csv"
FWHM_TO_SIGMA = 1.0 / 2.35482

#: Histogram edges for the fitted sigma. Chosen to RESOLVE the valley the ceiling sits in,
#: not to make one appear: the ceiling range (0.057-0.104 A) spans three of these bins.
SIGMA_EDGES = [0.0, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10,
               0.15, 0.20, 0.30, 0.39, 0.3999, 0.4001]


def _verdict(reason: object) -> str:
    """Collapse an exclusion sentence to its coded verdict.

    The REW branch spells its own number into the reason (`REW -4.578 above ...`), so a
    naive `value_counts` returns one row per line rather than one row per cause — which is
    how a 303-line verdict reads as 200 distinct verdicts.
    """
    if not isinstance(reason, str) or not reason:
        return "IN-AGGREGATE"
    if reason.startswith("REW "):
        return "REW-SATURATED"
    return reason.split(":")[0]


def _sigma_inst(R: float, wave_A: np.ndarray) -> np.ndarray:
    return wave_A / R * FWHM_TO_SIGMA


def _lab_wavelengths() -> np.ndarray:
    cg = pd.read_csv(ROOT / "data" / "linelists" / "canonical_gf.csv", low_memory=False)
    lab = cg[(cg.species == "Fe I")
             & cg.gf_tier.astype(str).str.contains("LAB", na=False)]
    return np.sort(lab.wavelength_air_A.values)


def _is_lab(waves: np.ndarray, lab: np.ndarray, tol: float = 0.005) -> np.ndarray:
    """Nearest-within-tolerance, never a rounded key (the RYA-703/704 mistake)."""
    i = np.searchsorted(lab, waves)
    best = np.full(waves.shape, np.inf)
    for off in (-1, 0, 1):
        j = np.clip(i + off, 0, len(lab) - 1)
        best = np.minimum(best, np.abs(lab[j] - waves))
    return best <= tol


def audit() -> dict:
    cat = pd.read_csv(CATALOG)
    idcol = cat.columns[0]
    lab = _lab_wavelengths()
    out: dict = {
        "broadening_kms": {
            "irreducible_sigma_floor": round(irreducible_sigma_kms(), 4),
            "max_stellar_sigma_ceiling": round(max_stellar_sigma_kms(), 4),
        },
        "arms": {},
    }

    for arm, fname in ARMS.items():
        path = EW_DIR / fname
        if not path.exists():
            raise SystemExit(
                f"{path} missing. Re-run:\n"
                f"  python3 scripts/measure_band_profilefit.py --element Fe --ion I "
                f"--lo 3780 --hi 6910 --instrument <arm>\n"
                f"Refusing to emit an evidence file describing a run that did not happen.")
        d = pd.read_csv(path)
        instrument = str(d.instrument.iloc[0])
        row = cat[cat[idcol].astype(str) == instrument]
        if row.empty:
            raise SystemExit(f"{instrument!r} absent from {CATALOG.name}")
        R = float(row.iloc[0]["resolving_power_max"])

        verdicts = d.excluded_reason.map(_verdict)
        fits = d.dropna(subset=["profile_sigma_A"]).copy()
        si = _sigma_inst(R, fits.wavelength_air_A.values)
        ceil_sigma = np.array([physical_ceiling_sigma_A(w, s)
                               for w, s in zip(fits.wavelength_air_A.values, si)])
        sig = fits.profile_sigma_A.values
        hist, _ = np.histogram(sig, bins=SIGMA_EDGES)

        ok = d[d.in_aggregate]
        lab_mask = _is_lab(d.wavelength_air_A.values, lab)

        out["arms"][arm] = {
            "instrument": instrument,
            "resolving_power": R,
            "candidates": int(len(d)),
            "verdicts": {k: int(v) for k, v in verdicts.value_counts().items()},
            "sigma_histogram": {f"{a}-{b}": int(n) for a, b, n
                                in zip(SIGMA_EDGES[:-1], SIGMA_EDGES[1:], hist)},
            "sigma_pinned_at_optimiser_bound": int(np.sum(np.abs(sig - 0.40) < 1e-6)),
            "ceiling_sigma_A": {"min": round(float(ceil_sigma.min()), 4),
                                "max": round(float(ceil_sigma.max()), 4)},
            # The property that makes the ceiling a MEASURED split rather than a chosen
            # cut: almost nothing lives in the decade around it.
            "fraction_between_1x_and_2x_ceiling":
                round(float(np.mean((sig > ceil_sigma) & (sig < 2 * ceil_sigma))), 4),
            "in_aggregate": int(len(ok)),
            "in_aggregate_implied_width_A": {
                "median": round(float(ok.implied_width_A.median()), 4),
                "p95": round(float(ok.implied_width_A.quantile(0.95)), 4),
                "max": round(float(ok.implied_width_A.max()), 4),
            },
            "in_aggregate_ew_mA": {
                "median": round(float(ok.ew_mA.median()), 2),
                "max": round(float(ok.ew_mA.max()), 2),
            },
            "lab_graded_lines_reaching_the_fitter": int(lab_mask.sum()),
            "lab_graded_fate": {k: int(v) for k, v
                                in verdicts[lab_mask].value_counts().items()},
        }
    return out


def main() -> int:
    res = audit()
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "width_ceiling_evidence.json").write_text(json.dumps(res, indent=2) + "\n")

    rows = []
    for arm, a in res["arms"].items():
        for verdict, n in sorted(a["verdicts"].items(), key=lambda kv: -kv[1]):
            rows.append(dict(arm=arm, instrument=a["instrument"], verdict=verdict,
                             n=n, of_candidates=a["candidates"]))
    pd.DataFrame(rows).to_csv(OUT / "verdict_census.csv", index=False)

    print(f"solar broadening sigma: floor {res['broadening_kms']['irreducible_sigma_floor']}"
          f" km/s, ceiling {res['broadening_kms']['max_stellar_sigma_ceiling']} km/s")
    for arm, a in res["arms"].items():
        iw = a["in_aggregate_implied_width_A"]
        print(f"\n{arm} (R={a['resolving_power']:.0f}) — {a['candidates']} candidates")
        print(f"  in aggregate {a['in_aggregate']}; pinned at the 0.40 A bound "
              f"{a['sigma_pinned_at_optimiser_bound']}")
        print(f"  ceiling {a['ceiling_sigma_A']['min']}-{a['ceiling_sigma_A']['max']} A; "
              f"{a['fraction_between_1x_and_2x_ceiling']:.1%} of fits sit within 1-2x it")
        print(f"  in-aggregate implied width median {iw['median']} A, max {iw['max']} A")
        print(f"  lab-graded reaching the fitter {a['lab_graded_lines_reaching_the_fitter']}"
              f": {a['lab_graded_fate']}")
    print(f"\nwrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
