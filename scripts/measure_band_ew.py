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
from pipeline.band_policy import check_intake, resolve, BandPolicyError  # noqa: E402

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


# ── Root-cause attribution ───────────────────────────────────────────────────
# Ryan, 2026-08-09: "In QA, we want to find root causes. Why did it fail? What is the
# mechanism? Is it our model? the data? Something wonky?"
#
# A symptom is not a cause. "GF ghost" says what we SEE; it does not say whether the
# fault lives in the atomic data, the observation, our physics, or our method -- and
# those four have different owners and different fixes. Every failure therefore carries
# a FAULT DOMAIN, the MECHANISM that produces the symptom, the DISCRIMINATOR that
# distinguished it from the alternatives, and the FIX that would resolve it.
#
# When the evidence does not separate two candidates we say UNKNOWN and name both.
# A confidently wrong root cause is worse than an honest undetermined one.
FAULT_DOMAINS = (
    "ATOMIC-DATA",   # our line list: wrong gf, wrong wavelength, wrong species
    "OBSERVATION",   # the spectrum: telluric, coverage gap, S/N, upstream normalisation
    "MODEL",         # our physics: predicted depth from the wrong atmosphere/abundance
    "METHOD",        # our measurement: window, continuum policy, EW-inversion regime
    "UNKNOWN",
)


def attribute_root_cause(w, f, centre, half_width, predicted_depth, symptom,
                         catalogue_waves) -> dict:
    """Given a failed line, work out WHERE the fault lives and HOW it produces the symptom."""
    cont = 1.0
    j = int(np.argmin(np.abs(w - centre)))
    depth_at = 1.0 - float(f[j]) / cont
    # The caller passes the stored reason, which carries a "FEATURE-VERIFICATION: "
    # prefix. Match on the tag itself rather than the start of the string.
    symptom = symptom.replace("FEATURE-VERIFICATION: ", "", 1)

    if symptom.startswith("GF-GHOST-ABSENT"):
        # Discriminator: is there a feature of ABOUT THE RIGHT DEPTH nearby? If yes, the
        # line is real and our WAVELENGTH is wrong. If nothing of that depth exists
        # anywhere near, the gf is wrong or the species is misassigned. These are both
        # ATOMIC-DATA faults but they have completely different fixes.
        near = np.abs(w - centre) <= 0.6
        if near.sum() and predicted_depth:
            depths = 1.0 - f[near] / cont
            k = int(np.argmax(depths))
            best_d, best_w = float(depths[k]), float(w[near][k])
            if 0.5 <= best_d / predicted_depth <= 2.0:
                return dict(
                    fault_domain="ATOMIC-DATA", mechanism="wavelength error in our line list",
                    discriminator=(f"a feature of depth {best_d:.3f} (predicted "
                                   f"{predicted_depth:.3f}) sits at {best_w:.3f}, "
                                   f"{best_w - centre:+.3f} A away — the line is REAL, "
                                   f"our position is wrong"),
                    fix="re-source this line's wavelength from a graded reference (NIST)")
        return dict(
            fault_domain="ATOMIC-DATA",
            mechanism="log gf far too strong, or the transition is assigned to the wrong species",
            discriminator=(f"nothing within 0.6 A has a depth resembling the predicted "
                           f"{predicted_depth:.3f}; observed at the position is {depth_at:.3f}"),
            fix="re-adjudicate log gf against a NIST-graded source; if it survives, the "
                "species assignment is suspect")

    if symptom.startswith("GF-GHOST"):
        # Present, correctly positioned, wrong STRENGTH. The direction discriminates:
        # too deep means our gf is too weak OR an unrecognised blend is adding depth --
        # both possible, so both are named rather than guessing. Too shallow has no
        # blend explanation (a blend cannot REMOVE absorption), so it is unambiguously
        # our atomic data.
        ratio_txt = f"{depth_at:.3f} observed vs {predicted_depth:.3f} predicted"
        if predicted_depth and depth_at > predicted_depth:
            return dict(
                fault_domain="ATOMIC-DATA",
                mechanism="log gf too weak, or an uncatalogued blend adds depth at this position",
                discriminator=(f"{ratio_txt} — the line is present and correctly placed, "
                               f"so this is not a positional error. Deeper than predicted "
                               f"has TWO explanations and this test does not separate "
                               f"them; a synthesis fit would"),
                fix="re-adjudicate log gf against NIST; if it holds, look for an "
                    "uncatalogued blend by fitting the profile")
        return dict(
            fault_domain="ATOMIC-DATA",
            mechanism="log gf too strong",
            discriminator=(f"{ratio_txt} — correctly placed and present, and SHALLOWER "
                           f"than predicted. A blend can only add absorption, never "
                           f"remove it, so a blend cannot explain this; the gf is wrong"),
            fix="re-adjudicate log gf against a NIST-graded source")

    if symptom.startswith("BLEND-DOMINATED"):
        # Discriminator: is the interloper in OUR catalogue? If yes, we knew about it and
        # our window was simply too wide -- a METHOD fault we own. If no, our line list is
        # missing a real solar line -- an ATOMIC-DATA/coverage fault.
        m = np.abs(w - centre) <= half_width
        i = int(np.argmin(f[m]))
        peak = float(w[m][i])
        known = catalogue_waves[np.abs(catalogue_waves - peak) < 0.05]
        if len(known):
            return dict(
                fault_domain="METHOD",
                mechanism="integration window wide enough to swallow a KNOWN neighbour",
                discriminator=(f"the dominant feature at {peak:.3f} IS in our catalogue "
                               f"({len(known)} entry/entries within 0.05 A)"),
                fix="narrow the window, or measure by profile fitting/synthesis which "
                    "models the neighbour instead of integrating over it")
        return dict(
            fault_domain="ATOMIC-DATA",
            mechanism="a real solar line missing from our list dominates the window",
            discriminator=(f"the dominant feature at {peak:.3f} (depth "
                           f"{1.0 - float(f[m][i]):.3f}) has NO catalogue entry within "
                           f"0.05 A — absence of a neighbour in our list is not absence "
                           f"in the spectrum"),
            fix="extend the IR line list from a graded source before measuring this region")

    if "saturation ceiling" in symptom:
        return dict(
            fault_domain="METHOD",
            mechanism="EW->abundance inversion runs on the flat part of the curve of growth",
            discriminator=f"REW above {-4.9}; the line itself is real and well measured",
            fix="measure by synthesis, which uses the profile shape rather than inverting EW")

    return dict(fault_domain="UNKNOWN", mechanism="", discriminator="", fix="")


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

    i = int(np.argmin(f[m]))
    peak_at = float(w[m][i])
    offset = peak_at - centre
    depth = 1.0 - f[m][i] / cont
    misplaced = abs(offset) > max(0.05, half_width * 0.25)

    # Depth AT the catalogued position, which is a different question from the depth of
    # whatever happens to be deepest in the window. A ghost is diagnosed here.
    j = int(np.argmin(np.abs(w - centre)))
    depth_at = 1.0 - f[j] / cont
    ratio = (depth_at / predicted_depth) if (predicted_depth and predicted_depth > 0) else None

    # GF-GHOST-ABSENT: the catalogue promises a line and the Sun shows nothing there.
    # Checked BEFORE the position test -- otherwise an absent line gets blamed on
    # whatever neighbour happened to be deepest, which mislabels the real fault.
    if depth_at < 0.02:
        return False, (f"GF-GHOST-ABSENT: no absorption at the catalogued position "
                       f"(depth {depth_at:.3f}, predicted {predicted_depth:.3f}) — the "
                       f"line is absent from the spectrum or misplaced in the list")

    # GF-GHOST: the line IS there, at the right place, but nothing like the strength the
    # line parameters claim. That is an atomic-data fault, not a measurement fault.
    if ratio is not None and not (0.25 <= ratio <= 4.0) and not misplaced:
        return False, (f"GF-GHOST: observed depth {depth_at:.3f} vs predicted "
                       f"{predicted_depth:.3f} (x{ratio:.1f}) at the correct position — "
                       f"the spectrum and the line parameters disagree; gf or "
                       f"identification is suspect, not the measurement")

    # BLEND-DOMINATED: something else owns this window.
    if misplaced:
        return False, (f"BLEND-DOMINATED: deepest feature sits {offset:+.3f} A from the "
                       f"catalogued position (depth {depth:.3f} vs {depth_at:.3f} at the "
                       f"line) — the EW is an upper bound on a blend, not this line")

    # Position right, present, but strength still off -- report it as a ghost too.
    if ratio is not None and not (0.25 <= ratio <= 4.0):
        return False, (f"GF-GHOST: observed depth {depth_at:.3f} vs predicted "
                       f"{predicted_depth:.3f} (x{ratio:.1f})")
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
    ap.add_argument("--diagnostic-only", action="store_true",
                    help="acknowledge that the EWs are interval-integrated absorption, "
                         "not equivalent widths, and are for verification/root-cause only")
    a = ap.parse_args()

    if a.instrument != "kpno_solar_atlas":
        raise SystemExit(f"instrument {a.instrument!r} has no reader wired here yet. "
                         f"Add it to this driver rather than copying this file.")

    # INTAKE CHECK (RYA-713). This driver measures by interval integration, which the
    # optical control proved is not an equivalent width -- median EW ratio 0.773 against
    # the HARPS pool, 5x spread. The band policy now refuses it everywhere, which is the
    # correct outcome: this script is a DIAGNOSTIC harness (is the line there? is it a
    # ghost? is it blended?) and must not be mistaken for a measurement harness again.
    for edge in (a.lo, a.hi - 1e-6):
        try:
            check_intake(edge, "interval-integration", instrument=a.instrument)
        except BandPolicyError as e:
            if not a.diagnostic_only:
                raise SystemExit(
                    f"{e}\n\n  This driver only produces interval-integrated absorption. "
                    f"Re-run with --diagnostic-only to use it for feature verification and "
                    f"root-cause work, where the EW value is not the product.")
    if a.diagnostic_only:
        pol = resolve((a.lo + a.hi) / 2.0)
        print(f"  DIAGNOSTIC ONLY — band {pol.name}: interval integration is forbidden here")
        print(f"  for measurement. EW values are absorption-in-window, NOT equivalent widths.")
        print(f"  permitted for measurement: {pol.permitted_methods}"
              + (f" · telluric correction REQUIRED" if pol.telluric_required else ""))

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

    rows, skipped, causes = [], [], []
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
        # Root-cause every quarantine, including the saturation ones set in __post_init__.
        if not lm.in_aggregate:
            rc = attribute_root_cause(w, f, float(r.wave_air_A), hw,
                                      float(r.predicted_depth or 0.0),
                                      lm.excluded_reason, allw)
            causes.append(dict(wave=float(r.wave_air_A), symptom=lm.excluded_reason[:90], **rc))
        rows.append(lm)

    assert_single_element(rows, a.element)

    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    stem = f"{a.element}{a.ion}_{int(a.lo)}_{int(a.hi)}_{a.instrument}"
    df = pd.DataFrame([{k: v for k, v in vars(l).items()} for l in rows])
    df.to_csv(out / f"{stem}_ew.csv", index=False)
    (out / f"{stem}_skipped.json").write_text(json.dumps(skipped, indent=2))
    pd.DataFrame(causes).to_csv(out / f"{stem}_root_causes.csv", index=False)

    print(f"\n  measured {len(rows)}, skipped {len(skipped)}")
    if len(df):
        sat = df[df.in_aggregate == False]  # noqa: E712
        print(f"  EW range {df.ew_mA.min():.1f} - {df.ew_mA.max():.1f} mA")
        print(f"  quarantined (ALL causes, not only saturation): {len(sat)} "
              f"— measured and kept, excluded from the aggregate only")
    for s in skipped[:6]:
        print(f"    skip {s['wave']:.3f}: {s['reason'][:88]}")
    if causes:
        cf = pd.DataFrame(causes)
        print("\n  ROOT CAUSE — where the fault actually lives:")
        for dom, g in cf.groupby("fault_domain"):
            print(f"    {dom:12s} {len(g):2d}")
            for mech, gg in g.groupby("mechanism"):
                print(f"        {len(gg):2d} x {mech}")
    print(f"\n  wrote {out / (stem + '_ew.csv')}")


if __name__ == "__main__":
    main()
