#!/usr/bin/env python3
"""Measure EWs for one element in one band, on every instrument that covers it — RYA-713.

    python3 scripts/measure_band_ew.py --element Fe --ion I --lo 6910 --hi 9199

This is the MEASUREMENT half. It produces equivalent widths and nothing else: no
abundances, no engines, no corrections. The three products (1D-LTE, Engine A, Engine B)
are all built from this one set of EWs downstream, which is what makes them comparable —
they differ only in treatment, never in what was measured.

Element, ion, band and instrument are all arguments. There is no element symbol in the
logic, because the Ba->Al copy proved that a hand-adapted harness keeps its source's
identity in places nobody looks (RYA-701).
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline.band_products import (  # noqa: E402
    LineMeasurement, equivalent_width, assert_single_element)

ACCOUNTING = ROOT / "data" / "audit" / "line_accounting" / "per_line.csv"
OUT = ROOT / "data" / "measured" / "band_ew"

# Kitt Peak FTS: lm#### files, NNNN = segment start wavelength in NM. Three whitespace
# columns (nm air, residual flux, irradiance) -- same reader as RYA-459's intake.
KP_DIR = (Path("/Users/ryanschmitt/Documents/Exoplanet Codex/data/spectra/"
               "exoplanetcodex-data/Solar Calibration/Kitt Peak Flux Atlas"))

# Regions where the terrestrial atmosphere, not the Sun, sets the flux. A line here is
# not measurable from the ground without a telluric correction we have not applied.
# Which instruments arrive ALREADY continuum-normalised. This is a property of the data
# product, not a preference: Kitt Peak column 1 is residual flux (Kurucz divided by his
# continuum); HARPS arrives un-normalised and our pipeline sets the continuum. Two
# normalisation histories is a reason two instruments can disagree methodologically.
PRE_NORMALISED = {"kpno_solar_atlas": True, "harps": False, "iag_fts_solar_atlas": True}

TELLURIC = [(7600, 7640, "O2 A-band"), (9280, 9600, "H2O"), (11120, 11560, "H2O")]


def telluric_reason(wave: float) -> str:
    for lo, hi, name in TELLURIC:
        if lo <= wave <= hi:
            return f"inside the {name} telluric band ({lo}-{hi} A); not measurable here"
    return ""


def kp_segments() -> list[tuple[float, float, Path]]:
    """Inventory the atlas as (lo_A, hi_A, path). Reads each file's ACTUAL span rather
    than trusting the filename -- the lm#### stem is a start hint, not a guarantee."""
    segs = []
    for p in sorted(KP_DIR.glob("lm[0-9]*")):
        if not p.is_file():
            continue
        try:
            head = np.loadtxt(p, max_rows=1)
            tail = np.loadtxt(p, skiprows=max(0, sum(1 for _ in open(p)) - 2))
        except Exception:
            continue
        lo = float(np.atleast_2d(head)[0, 0]) * 10.0
        hi = float(np.atleast_2d(tail)[-1, 0]) * 10.0
        segs.append((lo, hi, p))
    return segs


def load_kp_window(segs, centre: float, pad: float) -> tuple[np.ndarray, np.ndarray, str]:
    """Load the atlas around one line. Spans segment boundaries when a window straddles
    two files -- a line near a seam is a real line, not a missing one."""
    lo, hi = centre - pad, centre + pad
    hits = [p for (a, b, p) in segs if not (b < lo or a > hi)]
    if not hits:
        raise LookupError(f"no Kitt Peak segment covers {centre:.3f} A")
    W, F = [], []
    for p in hits:
        arr = np.loadtxt(p)
        w = arr[:, 0] * 10.0
        m = (w >= lo) & (w <= hi)
        if m.any():
            W.append(w[m]); F.append(arr[m, 1])
    if not W:
        raise LookupError(f"segments cover {centre:.3f} but hold no points in the window")
    w = np.concatenate(W); f = np.concatenate(F)
    o = np.argsort(w)
    return w[o], f[o], ",".join(p.name for p in hits)


def window_half_width(waves: np.ndarray, centre: float,
                      floor: float = 0.12, cap: float = 0.45) -> float:
    """Half-width from the distance to the NEAREST NEIGHBOURING LINE, not a constant.

    A fixed window is what makes a crowded line swallow its neighbour and an isolated
    line clip its own wings. Half the gap to the nearest catalogued line, bounded.
    """
    other = waves[np.abs(waves - centre) > 1e-4]
    if not len(other):
        return cap
    gap = float(np.min(np.abs(other - centre)))
    return float(np.clip(gap / 2.0, floor, cap))


def verify_feature(w: np.ndarray, f: np.ndarray, centre: float, half_width: float,
                   predicted_depth: float) -> tuple[bool, str]:
    """Is the thing we just integrated actually the line we asked for?

    Window-integrated EW answers "how much absorption is in this interval", which is only
    the line's EW if the line is (a) present, at (b) the catalogued position, and (c) the
    dominant feature there. In the crowded IR none of those can be assumed, and our line
    inventory is too sparse to flag the blends -- an empty neighbour list means our
    CATALOGUE is empty there, not the spectrum.

    Three checks, each returning a named reason. A line that fails is still measured and
    still reported; it is marked so its number is never mistaken for a clean EW.
    """
    cont = float(np.percentile(f, 95))
    m = np.abs(w - centre) <= half_width
    if m.sum() < 3:
        return False, "too few points in the window to verify"

    # (a) is the deepest thing in this window AT the catalogued position?
    i = int(np.argmin(f[m]))
    peak_at = float(w[m][i])
    offset = peak_at - centre
    depth = 1.0 - f[m][i] / cont
    if abs(offset) > max(0.05, half_width * 0.25):
        return False, (f"deepest feature sits {offset:+.3f} A from the catalogued "
                       f"position (depth {depth:.3f}) — the window is dominated by a "
                       f"DIFFERENT feature, so this EW is an upper bound on a blend, "
                       f"not this line's EW")

    # (b) is there anything there at all?
    if depth < 0.02:
        return False, (f"no absorption at the catalogued position (depth {depth:.3f}) — "
                       f"the line is absent from the spectrum or misplaced in the list")

    # (c) does the observed depth resemble what the line parameters predict?
    if predicted_depth and predicted_depth > 0:
        ratio = depth / predicted_depth
        if not (0.25 <= ratio <= 4.0):
            return False, (f"observed depth {depth:.3f} vs predicted {predicted_depth:.3f} "
                           f"(x{ratio:.1f}) — the line parameters and the spectrum "
                           f"disagree; gf or identification is suspect")
    return True, ""


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--element", required=True)
    ap.add_argument("--ion", default="I")
    ap.add_argument("--lo", type=float, required=True, help="band start, Angstrom")
    ap.add_argument("--hi", type=float, required=True, help="band end, Angstrom")
    ap.add_argument("--instrument", default="kpno_solar_atlas")
    ap.add_argument("--max-lines", type=int, default=None)
    ap.add_argument("--depth-min", type=float, default=0.05)
    ap.add_argument("--depth-max", type=float, default=0.60)
    ap.add_argument("--out", default=str(OUT))
    a = ap.parse_args()

    if a.instrument != "kpno_solar_atlas":
        raise SystemExit(f"instrument {a.instrument!r} has no reader wired here yet. "
                         f"Add it to this driver rather than copying this file.")

    acc = pd.read_csv(ACCOUNTING)
    sel = acc[(acc.element == a.element) & (acc.ion == a.ion) &
              (acc.wave_air_A >= a.lo) & (acc.wave_air_A <= a.hi) &
              acc.predicted_depth.between(a.depth_min, a.depth_max) &
              acc.instruments.notna()].copy()
    sel = sel.sort_values("wave_air_A").reset_index(drop=True)
    if a.max_lines:
        # Even sampling across the band, not the first N -- the first N is one corner.
        idx = np.unique(np.linspace(0, len(sel) - 1, a.max_lines).astype(int))
        sel = sel.iloc[idx].reset_index(drop=True)

    print(f"{a.element} {a.ion}  band {a.lo:.0f}-{a.hi:.0f} A  "
          f"instrument {a.instrument}  candidates {len(sel)}")

    segs = kp_segments()
    print(f"  atlas segments inventoried: {len(segs)}")
    allw = acc[(acc.element.notna())].wave_air_A.values

    rows, skipped = [], []
    for _, r in sel.iterrows():
        why = telluric_reason(r.wave_air_A)
        if why:
            skipped.append(dict(wave=r.wave_air_A, reason=why)); continue
        hw = window_half_width(allw, float(r.wave_air_A))
        try:
            w, f, src = load_kp_window(segs, float(r.wave_air_A), pad=hw * 3.0)
            # Kitt Peak ships residual flux -- already normalised. Say so explicitly
            # rather than re-normalising on top of it (see band_products docstring).
            ew, method, concern = equivalent_width(
                w, f, float(r.wave_air_A), hw, pre_normalised=PRE_NORMALISED[a.instrument])
            ok, why = verify_feature(w, f, float(r.wave_air_A), hw,
                                     float(r.predicted_depth))
        except Exception as e:
            skipped.append(dict(wave=r.wave_air_A, reason=f"{type(e).__name__}: {e}"))
            continue
        lm = LineMeasurement(
            element=a.element, ion=a.ion, wavelength_air_A=float(r.wave_air_A),
            instrument=a.instrument, ew_mA=ew,
            ew_method=f"{method}; segment(s) {src}; half-width {hw:.3f} A from line separation")
        if concern:
            lm.ew_method += f" | CONCERN: {concern}"
        if not ok:
            # Measured, kept, reported -- and barred from the aggregate with its reason.
            # Quarantine, not a cull (RYA-711).
            lm.in_aggregate = False
            lm.excluded_reason = f"FEATURE-VERIFICATION: {why}"
        rows.append(lm)

    assert_single_element(rows, a.element)

    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    stem = f"{a.element}{a.ion}_{int(a.lo)}_{int(a.hi)}_{a.instrument}"
    df = pd.DataFrame([{k: v for k, v in vars(l).items()} for l in rows])
    df.to_csv(out / f"{stem}_ew.csv", index=False)
    (out / f"{stem}_skipped.json").write_text(json.dumps(skipped, indent=2))

    print(f"\n  measured {len(rows)}, skipped {len(skipped)}")
    if len(df):
        sat = df[df.in_aggregate == False]  # noqa: E712
        print(f"  EW range {df.ew_mA.min():.1f} - {df.ew_mA.max():.1f} mA")
        print(f"  quarantined by saturation ceiling: {len(sat)} "
              f"(measured and kept, excluded from aggregate only)")
    for s in skipped[:6]:
        print(f"    skip {s['wave']:.3f}: {s['reason'][:88]}")
    print(f"\n  wrote {out / (stem + '_ew.csv')}")


if __name__ == "__main__":
    main()
