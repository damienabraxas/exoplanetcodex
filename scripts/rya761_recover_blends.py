#!/usr/bin/env python3
"""RYA-761 — recover quarantined lines with a two-component fit. Measured, not assumed.

    python3 scripts/rya761_recover_blends.py --element Fe --ion I --lo 6910 --hi 9199
    python3 scripts/rya761_recover_blends.py --element Al --ion I --lo 6910 --hi 9199

WHAT IS RECOVERED, AND WHAT IS NOT
----------------------------------
Only **FIT-PINNED** and **BLEND-DOMINATED**. Both are one failure wearing two faces: the
window holds more than one line and a one-component model was asked to describe it.
Saturated lines, GF-GHOST and GF-GHOST-ABSENT are NOT touched -- a saturated line carries
no EW information and a ghost is an atomic-data fault (RYA-354), neither of which a
richer profile model can fix.

THE PIN IS NOT RELAXED. `permits_profile_fit_for_line` and the sigma floor are untouched.
What changes is the number of components, not the threshold (RYA-761: *"the pin is the
detector working. Fix the model, not the threshold."*).

THE SEPARABILITY TIERS ARE MEASURED
-----------------------------------
`scripts/rya761_separability.py` injected 7488 synthetic blends and recovered them. With
a PERFECT interloper wavelength the floor is 0.06 A (1.33 sigma), and 0.08 A for the hard
case where the neighbour is deeper than the target -- the BLEND-DOMINATED class.

But the fit holds the interloper centre FIXED (freeing it lets the optimiser put both
components on one line), so a wrong catalogue wavelength pins the second component in the
wrong place, and that -- not the optimiser -- is what actually limits real data:

    catalogue error   median EW error at separation 0.10 / 0.20 / 0.30 / 0.60 A
        0.00 A            0.001   0.001   0.001   0.001
        0.01 A            0.055   0.029   0.013   0.006
        0.03 A            0.535   0.110   0.059   0.032
        0.05 A            0.565   0.211   0.117   0.065

Two catalogues quoting the same transition disagree by up to ~0.03 A (RYA-704). So a
recovery at 0.10 A separation is worthless unless the wavelength is known to better than
0.01 A, while one at 0.30 A survives a 0.03 A error. Every recovered line is therefore
emitted with its separation and a tier, and the tier is the honest read:

    ROBUST       sep >= 0.30 A -- tolerates the RYA-704 scale of catalogue disagreement
    PROVISIONAL  0.08 - 0.30 A -- needs the wavelength good to ~0.01 A; verify per line
    UNRECOVERABLE sep < 0.08 A -- below the measured floor even with a perfect catalogue

RYA-711: a recovered line keeps BOTH facts -- its original quarantine reason and its
recovery. Recovery is not promotion; the line still faces the RYA-760 gf tier split like
any other before it can enter a product.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline.lines_fit import _fit_blend, _fit_profile, _voigt_abs  # noqa: E402
from scripts.measure_band_ew import kp_segments, load_kp_window      # noqa: E402
from scripts.measure_band_profilefit import (CONT_HALF_A, FIT_HALF_A,  # noqa: E402
                                             instrument_sigma)

ACCOUNTING = ROOT / "data" / "audit" / "line_accounting" / "per_line.csv"
EW_DIR = ROOT / "data" / "measured" / "band_ew"
OUT_DIR = ROOT / "data" / "results" / "rya761"

# Measured by scripts/rya761_separability.py — see the module docstring. These are
# separations in Angstroms, not guesses.
FLOOR_PERFECT_CATALOGUE_A = 0.06
FLOOR_HARD_CASE_A = 0.08
FLOOR_ROBUST_A = 0.30

RECOVERABLE = ("FIT-PINNED", "BLEND-DOMINATED")


def tier(sep_A: float) -> str:
    if sep_A >= FLOOR_ROBUST_A:
        return "ROBUST"
    if sep_A >= FLOOR_HARD_CASE_A:
        return "PROVISIONAL"
    return "UNRECOVERABLE"


def quarantine_class(reason: str) -> str | None:
    """Which recoverable class a line is in, or None if it is not ours to fix."""
    if not isinstance(reason, str):
        return None
    for k in RECOVERABLE:
        if k in reason:
            return k
    return None


def component_ew_mA(depth: float, sigma: float, gamma: float) -> float:
    """EW of ONE fitted component, integrated from its own parameters, in mA."""
    x = np.linspace(-3.0, 3.0, 4001)
    return float(np.trapezoid(1.0 - _voigt_abs(x, 0.0, depth, sigma, gamma), x)) * 1000.0


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--element", required=True)
    ap.add_argument("--ion", default="I")
    ap.add_argument("--lo", type=float, required=True)
    ap.add_argument("--hi", type=float, required=True)
    ap.add_argument("--instrument", default="kpno_solar_atlas")
    ap.add_argument("--max-neighbours", type=int, default=2,
                    help="components beyond the target (default 2 -> up to 3 total)")
    ap.add_argument("--ew-csv", help="read the measured EWs from this path instead of "
                                     "data/measured/band_ew/. Lets the recovery run "
                                     "against a measurement that has not landed on main "
                                     "yet, WITHOUT staging a foreign artifact here.")
    a = ap.parse_args()

    stem = f"{a.element}{a.ion}_{int(a.lo)}_{int(a.hi)}_{a.instrument}_PROFILEFIT"
    src = Path(a.ew_csv) if a.ew_csv else EW_DIR / f"{stem}_ew.csv"
    if not src.exists():
        raise SystemExit(f"no measured EWs at {src}")
    ew = pd.read_csv(src)
    acc = pd.read_csv(ACCOUNTING)
    allw = np.sort(acc.wave_air_A.values.astype(float))

    ew["quarantine_class"] = ew.excluded_reason.map(quarantine_class)
    pool = ew[ew.quarantine_class.notna()].copy()
    print(f"{a.element} {a.ion} {a.lo:.0f}-{a.hi:.0f} — {len(ew)} measured, "
          f"{len(pool)} in the recoverable quarantine "
          f"({', '.join(f'{k}:{int((pool.quarantine_class == k).sum())}' for k in RECOVERABLE)})")
    if pool.empty:
        raise SystemExit("nothing to recover")

    segs = kp_segments()
    rows = []
    for _, r in pool.iterrows():
        c = float(r.wavelength_air_A)
        d = np.abs(allw - c)
        near = allw[(d > 1e-4) & (d <= FIT_HALF_A)]
        near = near[np.argsort(np.abs(near - c))][: a.max_neighbours]
        sep = float(np.min(np.abs(near - c))) if len(near) else np.nan

        rec = dict(wavelength_air_A=c, quarantine_class=r.quarantine_class,
                   original_reason=str(r.excluded_reason)[:120],
                   ew_single_mA=float(r.ew_mA), n_neighbours=len(near),
                   nearest_sep_A=sep, tier=(tier(sep) if len(near) else "NO-NEIGHBOUR"))

        if not len(near):
            # No catalogued neighbour: this is the "missing real solar line" atomic-data
            # fault (98 such windows in RYA-713), not a fitting problem. Different ticket.
            rec.update(recovered=False,
                       outcome="NO-CATALOGUED-NEIGHBOUR: interloper absent from our list "
                               "-> atomic-data fault (RYA-713/354), not recoverable here")
            rows.append(rec); continue

        if rec["tier"] == "UNRECOVERABLE":
            rec.update(recovered=False,
                       outcome=f"BELOW-MEASURED-FLOOR: separation {sep:.3f} A < "
                               f"{FLOOR_HARD_CASE_A} A, not identifiable even with a "
                               f"perfect catalogue wavelength")
            rows.append(rec); continue

        try:
            w, f, _ = load_kp_window(segs, c, pad=CONT_HALF_A * 1.3)
            m = np.abs(w - c) <= FIT_HALF_A
            wf, ff = w[m], f[m]
            if wf.size < 40:
                raise ValueError(f"only {wf.size} points in the fit window")
            s_init, s_min, s_max = instrument_sigma(a.instrument, c)
            kw = dict(sigma_init=s_init, sigma_min=s_min, sigma_max=s_max)

            centres = [c] + [float(x) for x in near]
            p2, _, kind2, chi2_2 = _fit_blend(wf, ff, centres, **kw)
            p1, _, kind1, chi2_1 = _fit_profile(wf, ff, c, **kw)

            if p2 is None:
                rec.update(recovered=False, outcome="FIT-FAILED: blend fit did not converge")
                rows.append(rec); continue

            sig2 = float(p2[2]); gam2 = float(p2[3]) if kind2 == 'voigt' else 0.0
            pinned = sig2 <= s_min * 1.001
            ew2 = component_ew_mA(float(p2[1]), sig2, gam2)

            rec.update(profile_type=kind2, chi2_two=chi2_2, chi2_one=chi2_1,
                       sigma_fit_A=sig2, still_pinned=bool(pinned),
                       ew_two_mA=ew2,
                       interloper_depths=";".join(f"{x:.4f}" for x in p2[4 if kind2 == 'voigt' else 3:]))

            if pinned:
                rec.update(recovered=False,
                           outcome="STILL-PINNED: sigma remains on the instrumental floor "
                                   "with the richer model — the feature is not two clean "
                                   "components either")
            elif not np.isfinite(chi2_2) or chi2_2 > chi2_1:
                rec.update(recovered=False,
                           outcome=f"NO-IMPROVEMENT: chi2 {chi2_2:.2e} not better than the "
                                   f"single-component {chi2_1:.2e}")
            else:
                rec.update(recovered=True,
                           outcome=f"RECOVERED ({rec['tier']}): EW {ew2:.2f} mA from the "
                                   f"target's own component; chi2 {chi2_1:.2e} -> {chi2_2:.2e}")
        except Exception as e:
            rec.update(recovered=False, outcome=f"ERROR: {type(e).__name__}: {str(e)[:70]}")
        rows.append(rec)

    df = pd.DataFrame(rows)
    n_rec = int(df.recovered.fillna(False).sum())
    print(f"\n{'='*78}\nRECOVERY — measured, not projected\n{'='*78}")
    print(f"  attempted on {len(df)} quarantined lines")
    print(f"  RECOVERED   {n_rec}  ({100.0*n_rec/max(len(df),1):.1f}%)")
    if n_rec:
        rr = df[df.recovered.fillna(False)]
        for t in ("ROBUST", "PROVISIONAL"):
            k = int((rr.tier == t).sum())
            print(f"    {t:<12} {k}")
        print(f"  by quarantine class:")
        for k in RECOVERABLE:
            sub = df[df.quarantine_class == k]
            print(f"    {k:<16} {int(sub.recovered.fillna(False).sum())} of {len(sub)}")
    print("\n  not recovered, by reason:")
    for reason, k in (df[~df.recovered.fillna(False)].outcome
                      .str.split(":").str[0].value_counts().items()):
        print(f"    {reason:<26} {k}")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    f = OUT_DIR / f"{stem}_blend_recovery.csv"
    df.to_csv(f, index=False)
    print(f"\n  wrote {f.relative_to(ROOT)}")
    print("  NOTE: recovery is not promotion — every recovered line still faces the "
          "RYA-760 gf tier split before it can enter a product.")


if __name__ == "__main__":
    main()
