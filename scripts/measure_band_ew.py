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
import os
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import NamedTuple

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline.band_products import (  # noqa: E402
    carried_ep,
    LineMeasurement, equivalent_width, assert_single_element)
from pipeline.band_policy import check_intake, resolve, BandPolicyError  # noqa: E402
from pipeline import kp_atlas_integrity as kp_integrity  # noqa: E402  RYA-938
from pipeline.prenormalised_guard import (  # noqa: E402  RYA-1026
    assert_data_matches_declaration)
from config.constants import codex_path, codex_root, PATHS  # RYA-810 path register

ACCOUNTING = ROOT / "data" / "audit" / "line_accounting" / "per_line.csv"

OUT = ROOT / "data" / "measured" / "band_ew"

# Kitt Peak FTS: lm#### files, NNNN = segment start wavelength in NM. Three whitespace
# columns (nm air, residual flux, irradiance) -- same reader as RYA-459's intake.
#
# The atlas lives on BOTH machines (Ryan, 2026-08-09: "the spectra should be on Sirius as
# well"), so the path is resolved rather than hardcoded to one host. Compute runs on
# Sirius; hardcoding the Mac path is what made the first synthesis control fail with
# "no Kitt Peak segment covers 4065.381 A" -- a coverage-shaped error message for what
# was really a missing directory.
_KP_CANDIDATES = (
    os.environ.get("CODEX_KP_ATLAS", ""),
    str(codex_path('data.spectra_local') / 'Solar Calibration' / 'Kitt Peak Flux Atlas'),
    str(codex_path('data.spectra_kitt_peak')),
)


def _resolve_kp_dir() -> Path:
    for c in _KP_CANDIDATES:
        if c and Path(c).is_dir() and any(Path(c).glob("lm[0-9]*")):
            return Path(c)
    raise SystemExit(
        "Kitt Peak atlas not found. Looked in:\n  " + "\n  ".join(x for x in _KP_CANDIDATES if x)
        + "\nSet CODEX_KP_ATLAS, or stage the atlas. Failing here rather than reporting "
          "every line as 'no segment covers' -- that message described coverage when the "
          "real fault was a missing directory.")


KP_DIR = _resolve_kp_dir()

# Regions where the terrestrial atmosphere, not the Sun, sets the flux. A line here is
# not measurable from the ground without a telluric correction we have not applied.
# Which instruments arrive ALREADY continuum-normalised. This is a property of the data
# product, not a preference: Kitt Peak column 1 is residual flux (Kurucz divided by his
# continuum); HARPS arrives un-normalised and our pipeline sets the continuum. Two
# normalisation histories is a reason two instruments can disagree methodologically.
# crires_plus: FALSE and it matters. These are reduced but UN-normalised IDPs in adu
# (TUNIT2 = "adu", flux runs to ~1.9e5), unlike Kitt Peak residual flux or Elgueta's
# sp/. Treating them as pre-normalised is the RYA-713 continuum defect, which measured
# EWs low by a median 11.7 % and up to 71.4 %.
# crires_plus: BOTH, and that is the whole of RYA-904. The raw Vesta IDPs are reduced
# but UN-normalised, in adu (TUNIT2 = "adu", flux runs to ~1.9e5); the RYA-794 corrected
# Y arm arrives already continuum-normalised (median 0.9985). One instrument, two
# holdings, opposite answers -- so this cannot be an instrument-keyed dict and stay true.
# It is derived from `_INSTRUMENT_HOLDINGS` below, keyed by HOLDING, so the flag and the
# loader that serves the data cannot drift apart. Treating a normalised product as
# un-normalised (or the reverse) is the RYA-713 continuum defect, which measured EWs low
# by a median 11.7 % and up to 71.4 %.
#
# Defined after `_INSTRUMENT_HOLDINGS`; see `PRE_NORMALISED` further down.

# RYA-786: the band set lives in pipeline/telluric_policy.py and is imported, never
# re-typed. The list that used to sit here was both incomplete and too narrow — it had
# the O2 A-band at 7600-7640 (it runs to ~7685), no O2 B-band at all, and was missing the
# ~7160-7340 and ~8100-8400 H2O complexes. A line in an unlisted band gets measured as if
# it were clean, which is the silent version of this bug.
from pipeline.telluric_policy import (TELLURIC_BANDS as TELLURIC,  # noqa: E402
                                      exclusion as _telluric_exclusion)


def telluric_reason(wave: float, instrument: str | None = None) -> str:
    """RYA-786: the INSTRUMENT decides, not just the wavelength — a telluric-corrected
    atlas (IAG) keeps lines that an uncorrected one (KPNO) must drop."""
    return _telluric_exclusion(wave, instrument)


# ── the IAG arm (RYA-783) ────────────────────────────────────────────────────
#
# The second instrument. Every Fe band product so far is Kitt Peak ONLY, so the product
# key (instrument x band x engine, RYA-712) has had a single instrument in it and there
# has been no cross-instrument check anywhere — the very test that settled Al 6696 vs 6698.
#
# IAG is a genuinely independent arm, not a second copy of the same thing:
#
#   * different telescope and spectrograph — Goettingen VVT FTS at R ~ 700 000 against
#     Kitt Peak's ~300-500 k;
#   * different TELLURIC BASIS — Baker, Blake & Reiners 2020 (ApJS 247, 24;
#     Zenodo 10.5281/zenodo.3598136) is telluric-CORRECTED, where the Kitt Peak atlas
#     carries its tellurics and we exclude lines per-line (RYA-786). So IAG can reach
#     lines KPNO must quarantine, and disagreement between the arms inside a telluric
#     band is diagnostic rather than noise.
#
# THE TRAP, RECORDED IN THE CATALOG BEFORE I GOT HERE: column 0 is VACUUM WAVENUMBER in
# cm^-1, not a wavelength. `instrument_catalog.reduction_requirements` says it outright —
# "convert 1e8/wn then vac_to_air. Reading col0 as a wavelength returns 9387-24700 and a
# confident wrong answer." Our line lists are AIR wavelengths, so both steps are required
# and `vac_to_air` is imported from pipeline.uv_conditioning rather than re-typed here.

IAG_FITS = Path(str(codex_path('data.solar_iag_baker2020')))
_iag_cache: dict = {}


def iag_atlas() -> tuple[np.ndarray, np.ndarray]:
    """The IAG telluric-corrected atlas as (wave_air_A, flux), ascending. Cached."""
    if "wf" in _iag_cache:
        return _iag_cache["wf"]
    from astropy.io import fits
    from pipeline.uv_conditioning import vac_to_air
    if not IAG_FITS.exists():
        raise LookupError(f"IAG atlas not staged at {IAG_FITS} (Sirius-only, RYA-485)")
    with fits.open(IAG_FITS) as h:
        d = h[1].data
        wn = np.asarray(d["v"], dtype=float)      # VACUUM WAVENUMBER, cm^-1
        fl = np.asarray(d["s"], dtype=float)
    w_air = vac_to_air(1.0e8 / wn)                # -> vacuum A -> air A
    o = np.argsort(w_air)
    _iag_cache["wf"] = (w_air[o], fl[o])
    return _iag_cache["wf"]



#: The OTHER IAG atlas. Two files with DIFFERENT telluric states and DIFFERENT reach have
#: lived under one instrument_id since RYA-944 flagged it, which is the RYA-904 defect
#: shape and is why "is IAG corrected?" had no single answer. They are separate holdings
#: now, each declaring its own span, so selection can never silently pick the wrong one.
#:
#: WHY THE RAW ATLAS IS KEPT AND NOT RETIRED: Baker+2020 starts at 5002.5 A and Reiners
#: reaches 4048.6 A, so ~954 A of the blue is REINERS-ONLY. Baker cannot serve it at all.
#:
#: 🔴 AND IN THAT BLUE STRETCH THERE IS NOTHING TO CORRECT. The bluest registered telluric
#: complex is the O2 B-band at 6867 A, so no band this reader can reach carries one.
#: MEASURED, not argued -- 4100-5000 A against KP2005's residual atlas, which is
#: telluric-free AT SOURCE and does reach the blue:
#:     Reiners        median 0.9256, 10.48% of pixels below 0.5
#:     KP2005 residual median 0.9049, 10.35% below 0.5
#: and Reiners is not systematically DEEPER (mean difference +0.109, i.e. shallower).
#: Telluric contamination shows as one-sided deepening concentrated in bands; this is
#: symmetric scatter from a resolution difference. So `telluric_applied` for the blue
#: holding is `applied` in the only sense that matters -- there is no absorption to remove.
IAG_REINERS = Path(str(codex_path('data.solar_iag_atlas')))
_iag_reiners_cache: dict = {}

#: Where Baker+2020 takes over. Below this the blue holding is the ONLY IAG option; above
#: it, the corrected atlas covers the same wavelengths and must be preferred.
#: 🔴 AIR, not vacuum. Measured 5002.5 A VACUUM on the file and converted -- both atlases
#: store VACUUM WAVENUMBER, while `covers()` and every line list here are AIR. Declaring a
#: vacuum number as an air span would misplace the boundary by ~1.4 A and silently hand
#: ~1.4 A of blue to the atlas that cannot serve it.
IAG_BAKER_BLUE_EDGE_A = 5001.10


def iag_reiners_atlas() -> tuple[np.ndarray, np.ndarray]:
    """Reiners+2016 IAG as (wave_air_A, flux), ascending. Cached.

    Same vacuum-wavenumber convention as the Baker atlas -- col0 is cm^-1, NOT a
    wavelength, and reading it as one returns a confident wrong answer.
    """
    if "wf" in _iag_reiners_cache:
        return _iag_reiners_cache["wf"]
    import gzip
    from pipeline.uv_conditioning import vac_to_air
    if not IAG_REINERS.exists():
        raise LookupError(
            f"IAG Reiners+2016 atlas not staged at {IAG_REINERS} (Sirius-only). This is "
            f"the ONLY atlas reaching 4047.5-5001.1 A (air); Baker+2020 starts at "
            f"{IAG_BAKER_BLUE_EDGE_A} A and cannot substitute for it.")
    wn, fl = [], []
    with gzip.open(IAG_REINERS, "rt") as fh:
        for line in fh:
            q = line.split()
            if len(q) < 2:
                continue
            try:
                v, y = float(q[0]), float(q[1])
            except ValueError:
                continue                          # header / comment
            wn.append(v); fl.append(y)
    w_air = vac_to_air(1.0e8 / np.asarray(wn, float))
    f = np.asarray(fl, float)
    o = np.argsort(w_air)
    _iag_reiners_cache["wf"] = (w_air[o], f[o])
    return _iag_reiners_cache["wf"]


MIN_WINDOW_POINTS = 12


def _slice_window(w: np.ndarray, f: np.ndarray, centre: float, pad: float,
                  what: str, min_points: int = MIN_WINDOW_POINTS):
    """Mask an ascending (wave, flux) pair to a window, or say why it cannot.

    Shared rather than re-typed per arm (RYA-701: one Ba->Al copy produced 13 defects).
    The failure message carries the arm's own span, because "0 points near 21000 A" and
    "this arm stops at 6910 A" are different diagnoses and only the second is actionable.
    """
    m = (w >= centre - pad) & (w <= centre + pad)
    n = int(m.sum())
    if n < min_points:
        raise LookupError(f"{what} holds only {n} points near {centre:.3f} A "
                          f"(spans {w.min():.1f}-{w.max():.1f} A)")
    return w[m], f[m]


def load_iag_window(centre: float, pad: float) -> tuple[np.ndarray, np.ndarray, str]:
    w, f = iag_atlas()
    w, f = _slice_window(w, f, centre, pad, "IAG (Baker+2020, telluric-corrected)")
    return w, f, IAG_FITS.name


def load_iag_reiners_window(centre: float, pad: float):
    """The blue IAG arm, 4048.6-5002.5 A -- the stretch Baker+2020 cannot reach."""
    w, f = iag_reiners_atlas()
    w, f = _slice_window(w, f, centre, pad, "IAG (Reiners+2016, blue arm)")
    return w, f, IAG_REINERS.name


# ── the direct-solar HARPS arm (RYA-897, ported to the RYA-904 holding table) ─
#
# THE ARM THE WHOLE SCALE IS DIFFERENTIAL TO, and it was half-wired from day one.
# `PRE_NORMALISED` carried a "harps" entry from RYA-713 while `load_window` never got the
# dispatch branch, so the config declared an instrument the loader could not reach.
# Nothing scientific excluded HARPS; the path was simply never written (RYA-896/897).
#
# Why it matters more than a fourth data point: every Codex number is
# [X/H] = A(X)star - A(X)sun, and program stars are observed with HARPS. A solar cell
# measured on THIS harness, on THIS instrument, is the only anchor whose instrumental and
# methodological systematics cancel against the program stars.
#
# ⚠️ THE DOUBLE-NORMALISATION TRAP, avoided by construction. This file carries FOUR
# columns -- wavelength_air_A, flux_raw, continuum, flux_normalized -- so it is possible
# to hand the harness an already-normalised flux and have it set a continuum on top: the
# RYA-713 defect. We read `flux_raw` and the holding declares pre_normalised=False, so the
# harness sets the continuum ONCE, per window.
#
# 🔴 AND THAT CONTINUUM IS THE SUBJECT OF RYA-911. This product carries the pipeline's
# OWN fitted `continuum` column, so for the first time the harness's per-window continuum
# has a same-file reference to be read against -- see `reference_continuum()` below and
# the instrumentation in `measure_band_profilefit`. It is a REFERENCE, not an authority:
# which of the two is right is a question this module records the evidence for and does
# not answer by itself.
HARPS_CSV = Path(str(PATHS['solar_normalized']))
#: RYA-931 — the telluric-corrected sibling product, same columns, same normaliser.
#: A SECOND PATH, not a swap: `solar_harps` keeps serving what it always served.
HARPS_TELLCORR_CSV = HARPS_CSV.with_name('solar_normalized_tellcorr.csv')
_harps_cache: dict = {}


def harps_spectrum(csv: Path | None = None) -> tuple[np.ndarray, np.ndarray]:
    """The direct-solar HARPS spectrum as (wave_air_A, flux_NORMALISED), ascending.

    Air wavelengths already (the column says so, and RYA-643 conditioned this product to
    rest frame); no vacuum conversion, unlike IAG.

    🔴 RYA-911 — THIS RETURNS `flux_normalized`, AND THAT IS THE FIX. RYA-897 wired it to
    `flux_raw` with `pre_normalised=False`, which routed every window through
    `pipeline.lines_fit._local_renorm` — a degree-1 fit through the top 20% of pixels in
    the outer 25% of each +/-1.2 A edge strip, i.e. through 0.3 A strips of the crowded
    solar optical. MEASURED, on the harness's OWN pixels with no extrapolation and no
    reconstruction: that continuum sits BELOW THE OBSERVED FLUX on 8 of the 16
    in-aggregate Fe II lines, worst case max(F/C) = 1.204. Flux above the continuum is
    not a small error; it is impossible. The product's own column never exceeds 1.007.
    Isolated in a same-inputs control (same fitter, same window, only the continuum
    swapped) the local re-fit costs a MEDIAN 23.8% OF THE EQUIVALENT WIDTH -- which is
    the -0.34 dex RYA-897 measured and could not explain.
    
    ⚠️ THE TEST IS ONE-SIDED, DELIBERATELY. `max(F/C) > 1` convicts a continuum outright.
    `F/C < 1` everywhere does NOT convict one: in the blanketed blue there may be no
    unabsorbed pixel in the window at all, so a correct continuum SHOULD sit above every
    observed point. So this evidence rules the local re-fit out; it does not certify the
    product's column, which the blue lines (max F/C 0.90-0.94 at 4138/4413/4620) suggest
    may itself sit a little high there. That is recorded, not resolved.

    This is NOT the RYA-713 double-normalisation trap in reverse: the holding declares
    `pre_normalised=True` to match, so the harness uses unity as the continuum and sets
    none of its own. Read one and declare the other and you get RYA-713 exactly.
    """
    # RYA-931 — keyed by PRODUCT, not by instrument. HARPS now serves two holdings
    # (uncorrected and telluric-corrected) from the same reader, and a single-slot
    # cache would hand the second caller the first caller's pixels. That is the
    # RYA-904 defect shape exactly: one instrument silently standing for one holding.
    csv = HARPS_CSV if csv is None else Path(csv)
    key = str(csv)
    if key in _harps_cache:
        return _harps_cache[key]["wf"]
    if not csv.exists():
        raise LookupError(
            f"HARPS solar spectrum not staged at {csv}. This is the direct-solar "
            f"arm (Dumusque ESO 1102.D-0954); stage it rather than falling back to an "
            f"atlas, which would label an atlas measurement 'harps'.")
    d = pd.read_csv(csv)
    missing = {"wavelength_air_A", "flux_normalized", "continuum"} - set(d.columns)
    if missing:
        raise LookupError(
            f"{csv.name} lacks {sorted(missing)}; got {list(d.columns)}. The band "
            f"harness consumes this product's OWN fitted continuum (RYA-911) and will "
            f"not fall back to re-fitting one, because that re-fit is the defect.")
    w = np.asarray(d["wavelength_air_A"], dtype=float)
    f = np.asarray(d["flux_normalized"], dtype=float)
    cont = np.asarray(d["continuum"], dtype=float)
    ok = np.isfinite(w) & np.isfinite(f)
    w, f, cont = w[ok], f[ok], cont[ok]
    o = np.argsort(w)
    _harps_cache[key] = {"wf": (w[o], f[o])}
    #: kept for provenance: the ABSOLUTE continuum this normalisation divided by, so the
    #: per-line record can name a number rather than only the word "unity".
    _harps_cache[key]["cont"] = (w[o], cont[o])
    return _harps_cache[key]["wf"]


def load_harps_window(centre: float, pad: float,
                      csv: Path | None = None) -> tuple[np.ndarray, np.ndarray, str]:
    csv = HARPS_CSV if csv is None else Path(csv)
    w, f = harps_spectrum(csv)
    w, f = _slice_window(w, f, centre, pad, "HARPS")
    return w, f, csv.name


def reference_continuum(spec: "HoldingSpec", centre: float) -> float | None:
    """The SOURCE PRODUCT's own continuum at this wavelength, or None — RYA-911.

    🔴 THE POINT OF THIS FUNCTION IS THAT IT IS A READING. RYA-897's RCA compared the
    harness continuum against the pipeline's by REFITTING the window edges to approximate
    what the harness does — an inference about the code, not a measurement of it, and it
    said so. This module now hands the number out, and the harness hands out the one it
    actually placed, so the comparison downstream is two readings.

    ⚠️ NOT AN AUTHORITY. A holding whose product carries a continuum column gives us a
    second opinion, nothing more. Disagreement is a finding to explain, and the direction
    of the fix is decided by measuring which one sits at the unabsorbed flux level — not
    by preferring whichever was written down first.
    """
    if not spec.reference_continuum:
        return None
    if spec.reader in ("harps", "harps_tellcorr"):
        csv = HARPS_TELLCORR_CSV if spec.reader == "harps_tellcorr" else HARPS_CSV
        harps_spectrum(csv)                    # populates the cache
        cw, cc = _harps_cache.get(str(csv), {}).get("cont", (None, None))
        if cw is None:
            return None
        j = int(np.argmin(np.abs(cw - centre)))
        if abs(cw[j] - centre) > 0.05:
            return None
        return float(cc[j])
    return None


# ── the CRIRES+ arm (RYA-796) ────────────────────────────────────────────────
#
# 18 Vesta solar IDPs, 13 grating settings, 9479-24855 A. Our OWN IR spectrum of the Sun,
# as opposed to a published atlas -- which is the point: it is the same instrument class
# the science targets are observed with.
#
# FOUR PROPERTIES OF THESE FILES THAT A NAIVE READER GETS WRONG, all measured, not assumed:
#
#  1. WAVE IS NANOMETRES (TUNIT1 = 'nm'). Read as Angstrom it returns 949-2485 and a
#     confident wrong answer -- the same shape of trap as IAG's vacuum wavenumber column.
#     Converted once, here, and asserted against the catalog span rather than trusted.
#
#  2. SPECSYS = 'TOPOCENT', and `ESO TEL TARG RADVEL` is 0.0 -- a placeholder, not a
#     measurement. Vesta is a MOVING REFLECTOR, so the shift to remove is the two-leg
#     Sun->Vesta->observer rate (RYA-372), which no stellar BERV keyword can supply. A
#     wavelength measured off these frames as delivered is wrong by tens of km/s. This
#     loader therefore REFUSES raw frames by default (see the rest-frame gate below) --
#     RYA-643's defect class was exactly a rest-frame slip reaching a line fit, in four
#     copies.
#
#  3. DETECTOR GAPS ARE REAL. A setting is a comb of echelle orders across three
#     detectors, not a filled band: selection is on the actual WAVE array with QUAL == 0
#     and finite non-zero flux, never on WAVELMIN/WAVELMAX (RYA-377's warning). Measured:
#     807 488 good pixels of 830 000, and two genuine inter-band gaps at 1349-1439 nm and
#     1796-1946 nm that no setting covers.
#
#  4. CO-ADDING EPOCHS IS A ROTATION QUESTION, NOT AN SNR ONE. Five settings have two
#     epochs, and RYA-372 wrote its conditioned frames "per-frame only -- never coadded
#     (Vesta 5.3 h rotation)". Measured against the 5.342 h period, the two kinds are not
#     alike: Y1029 / J1232 / H1575 pairs are 28 / 32 / 4 minutes apart, so the sub-observer
#     longitude moves 32 / 36 / 5 degrees -- the same face. H1559 and H1582 are ~23.8 h
#     apart: 166 and 163 degrees, the OPPOSITE face. Co-adding the second kind averages
#     two different hemispheres of a reflecting body and calls it signal-to-noise.
#     So co-adding is gated on rotation phase, and never spans settings (their continua
#     and blaze differ -- RYA-794 Step 2).

CRIRES_ROTATION_H = 5.342          # Vesta sidereal rotation period
CRIRES_COADD_MAX_LON_DEG = 45.0    # co-add only within this sub-observer longitude change

_CRIRES_CANDIDATES = (
    os.environ.get("CODEX_CRIRES_VESTA", ""),
    str(codex_path('data.spectra_vesta_crires')),
)
#: Rest-frame conditioned output of pipeline/reflected_solar_rv.py `write_set`.
_CRIRES_REST_CANDIDATES = (
    os.environ.get("CODEX_CRIRES_VESTA_REST", ""),
    str(codex_root('data') / 'spectra' / 'vesta' / 'CRIRESPlus_rest' / 'vesta_crires'),
)
_crires_cache: dict = {}


class RestFrameNotConditioned(LookupError):
    """Raised when a topocentric frame is asked for as if it were science-ready.

    A LookupError subclass so a driver's existing "this arm cannot serve this line"
    handling still catches it -- but a distinct type, because "we do not hold this
    wavelength" and "we hold it in the wrong velocity frame" are different problems with
    different fixes, and collapsing them is how the second one gets ignored.
    """


def _resolve_dir(candidates, pattern: str):
    for c in candidates:
        if c and Path(c).is_dir() and any(Path(c).glob(pattern)):
            return Path(c)
    return None


def crires_frames() -> list[dict]:
    """Inventory the staged IDPs: one record per FRAME, with its measured good-pixel span.

    Cached. Reads the actual arrays rather than the span keywords, because the keywords
    describe the grating setting and the arrays describe the data (property 3).
    """
    if "frames" in _crires_cache:
        return _crires_cache["frames"]
    from astropy.io import fits
    d = _resolve_dir(_CRIRES_CANDIDATES, "*.fits")
    if d is None:
        raise LookupError(
            "CRIRES+ Vesta IDPs not staged. Looked in:\n  "
            + "\n  ".join(x for x in _CRIRES_CANDIDATES if x)
            + "\nSet CODEX_CRIRES_VESTA. Sirius-only (RYA-567).")
    out = []
    for p in sorted(d.glob("*.fits")):
        with fits.open(p) as h:
            ph, t = h[0].header, h[1].data
            unit = str(h[1].header.get("TUNIT1", "")).strip().lower()
            wave = np.asarray(t["WAVE"][0], dtype=float)
            flux = np.asarray(t["FLUX"][0], dtype=float)
            qual = np.asarray(t["QUAL"][0])
            # Property 1: convert ONCE, from the unit the file declares. An unexpected
            # unit is a loud stop, not a guess -- guessing is what produced RYA-794's
            # inverted nm/Angstrom confusion in the first place.
            if unit == "nm":
                wave_A = wave * 10.0
            elif unit in ("angstrom", "a", "0.1nm"):
                wave_A = wave
            else:
                raise LookupError(f"{p.name}: unexpected TUNIT1 {unit!r}; refusing to "
                                  f"guess the wavelength unit.")
            good = (qual == 0) & np.isfinite(flux) & (flux != 0)
            if not good.any():
                continue
            out.append(dict(path=p, name=p.name, setting=str(ph.get("ESO INS WLEN ID")),
                            mjd=float(ph.get("MJD-OBS", np.nan)),
                            snr=float(ph.get("SNR", np.nan)),
                            specsys=str(ph.get("SPECSYS", "")).strip().upper(),
                            wave_A=wave_A[good], flux=flux[good],
                            lo=float(wave_A[good].min()), hi=float(wave_A[good].max())))
    if not out:
        raise LookupError(f"no usable CRIRES+ frames under {d}")
    # Property 1, asserted rather than trusted: the catalog says 950-5300 nm.
    lo = min(r["lo"] for r in out); hi = max(r["hi"] for r in out)
    if not (9000.0 < lo < 30000.0 and 9000.0 < hi < 60000.0):
        raise LookupError(
            f"CRIRES+ wavelengths land at {lo:.0f}-{hi:.0f} A after unit conversion, "
            f"outside the instrument's catalogued 9500-53000 A. The unit handling is "
            f"wrong -- refusing to measure off it.")
    _crires_cache["frames"] = out
    return out


def crires_rest_frames() -> dict:
    """Rest-frame conditioned frames keyed by source stem, if RYA-372/373 has produced
    them. Empty dict when the conditioning leg has not been run."""
    if "rest" in _crires_cache:
        return _crires_cache["rest"]
    d = _resolve_dir(_CRIRES_REST_CANDIDATES, "*_rest.csv")
    rest = {}
    if d is not None:
        for p in sorted(d.glob("*_rest.csv")):
            rest[p.name.replace("_rest.csv", "")] = p
    _crires_cache["rest"] = rest
    return rest


def _coaddable(group: list[dict]) -> list[list[dict]]:
    """Split one setting's frames into rotation-safe co-add clusters (property 4)."""
    clusters: list[list[dict]] = []
    for fr in sorted(group, key=lambda r: r["mjd"]):
        placed = False
        for cl in clusters:
            dt_h = abs(fr["mjd"] - cl[0]["mjd"]) * 24.0
            if (360.0 * dt_h / CRIRES_ROTATION_H) % 360.0 <= CRIRES_COADD_MAX_LON_DEG:
                cl.append(fr); placed = True; break
        if not placed:
            clusters.append([fr])
    return clusters


def load_crires_window(centre: float, pad: float, *, allow_topocentric: bool = False
                       ) -> tuple[np.ndarray, np.ndarray, str]:
    """One CRIRES+ Vesta window as (wave_air_A, flux, provenance).

    `allow_topocentric` exists for the CONDITIONING leg itself (RYA-372/373), which must
    read raw frames to measure the shift it is going to remove. It is not a convenience
    flag for measurement, and the provenance string says so when it is used.
    """
    frames = crires_frames()
    hits = [r for r in frames if not (r["hi"] < centre - pad or r["lo"] > centre + pad)]
    if not hits:
        raise LookupError(
            f"no CRIRES+ setting covers {centre:.3f} A (arm spans "
            f"{min(f['lo'] for f in frames):.1f}-{max(f['hi'] for f in frames):.1f} A, "
            f"with real inter-band gaps -- coverage is a comb of settings, not a span)")

    # ── the rest-frame gate ──────────────────────────────────────────────────
    # These are TOPOCENT and Vesta is a moving reflector. Refuse unless the frame has
    # been conditioned, or the caller is the conditioning leg saying so explicitly.
    rest = crires_rest_frames()
    conditioned = [r for r in hits if Path(r["name"]).stem in rest]
    if not conditioned and not allow_topocentric:
        raise RestFrameNotConditioned(
            f"CRIRES+ frames covering {centre:.3f} A are SPECSYS="
            f"{sorted({r['specsys'] for r in hits})} and have not been rest-frame "
            f"conditioned. Vesta is a moving reflector: the shift to remove is the "
            f"two-leg Sun->Vesta->observer rate, which no BERV keyword carries "
            f"(ESO TEL TARG RADVEL is 0.0 in every one of these files). Condition them "
            f"with pipeline/reflected_solar_rv.py (write_set) and stage the output at "
            f"{_CRIRES_REST_CANDIDATES[-1]}, or pass allow_topocentric=True if you ARE "
            f"the conditioning leg. Measuring a wavelength off an unconditioned frame is "
            f"the RYA-643 defect.")

    # ── choose ONE setting, then co-add only within a rotation-safe cluster ──
    # Never across settings: their continua and blaze differ and none of this is
    # normalised (RYA-794 Step 2), so a cross-setting co-add manufactures a continuum.
    by_setting: dict[str, list[dict]] = {}
    for r in hits:
        by_setting.setdefault(r["setting"], []).append(r)

    def in_window(r):
        m = (r["wave_A"] >= centre - pad) & (r["wave_A"] <= centre + pad)
        return int(m.sum())

    best_setting = max(by_setting,
                       key=lambda s: (max(in_window(r) for r in by_setting[s]),
                                      max(r["snr"] for r in by_setting[s])))
    group = by_setting[best_setting]
    cluster = max(_coaddable(group),
                  key=lambda cl: (sum(in_window(r) for r in cl),
                                  sum(r["snr"] for r in cl)))

    ws, fs = [], []
    for r in cluster:
        w, f = _slice_window(r["wave_A"], r["flux"], centre, pad, f"CRIRES+ {r['name']}")
        ws.append(w); fs.append(f)
    if len(cluster) == 1:
        w, f = ws[0], fs[0]
    else:
        # SNR-weighted co-add onto the highest-SNR frame's grid. Same setting, same
        # rotation phase, so this is the same spectrum measured twice.
        ref = int(np.argmax([r["snr"] for r in cluster]))
        w = ws[ref]
        stack, wts = [], []
        for i, (wi, fi) in enumerate(zip(ws, fs)):
            stack.append(np.interp(w, wi, fi)); wts.append(max(cluster[i]["snr"], 1e-6))
        f = np.average(np.vstack(stack), axis=0, weights=np.asarray(wts))

    o = np.argsort(w)
    dropped = sorted(set(by_setting) - {best_setting})
    prov = (f"crires_plus {best_setting} "
            + "+".join(r["name"] for r in cluster)
            + (f" [co-added {len(cluster)} epochs within "
               f"{CRIRES_COADD_MAX_LON_DEG:.0f} deg rotation]" if len(cluster) > 1 else "")
            + (f" [also covered by {','.join(dropped)}, not merged]" if dropped else "")
            + ("" if conditioned else " [RAW TOPOCENTRIC -- conditioning leg only]"))
    return w[o], f[o], prov


# ── the CRIRES+ TELLURIC-CORRECTED Y arm (RYA-794 product, wired RYA-904) ────
#
# A SECOND CRIRES+ holding, and a different KIND of thing from the IDPs above: not raw
# frames we condition here, but the finished science-ready product RYA-794 built from
# Elgueta+2026's telluric-corrected `sp/Sun_Y_rv.dat`. Two columns, air wavelength and
# normalised flux, 13707 rows over 10280-10680 A.
#
# THREE PROPERTIES, EACH ALREADY MEASURED BY RYA-794 RATHER THAN ASSUMED:
#
#  1. TELLURIC-CORRECTED, measured: 0.10 % of window points below 0.5, against 51.3 %
#     for the Kitt Peak O2 A-band. Registered `telluric_applied=applied` on holding
#     `solar_crires_plus_y_rya794`, and it is the REGISTRY the gate reads -- this loader
#     never asserts a telluric state of its own (RYA-786).
#
#  2. ALREADY CONTINUUM-NORMALISED (median 0.9985). `pre_normalised` is TRUE here and
#     FALSE for the raw IDPs of the SAME INSTRUMENT. That disagreement is the reason
#     RYA-904 had to move normalisation to the holding key at the same time as the
#     loader key: wiring one without the other sets a continuum on an already-normalised
#     spectrum, which is the RYA-713 defect.
#
#  3. AIR WAVELENGTHS, IN THE SOLAR REST FRAME -- and this is checked, not inherited
#     from a column name. The instrument catalog says CRIRES+ delivers VACUUM and
#     TOPOCENT, so `wavelength_air_A` in a derived product is a claim. RYA-794's own
#     cross-match is the positive control: its 5 Fe I candidates sit within 0.001 A of
#     the VALD AIR wavelengths and the observed features within 0.04 A, where the
#     air-vacuum difference at 10500 A is ~2.9 A and Elgueta's `_rv` files are
#     RV-corrected. `_assert_air_rest_frame` below re-runs that check on load rather
#     than trusting the header, because a silently vacuum product would shift every line
#     by ~100 pixels and still fit something.
CRIRES_Y_CSV = Path(str(codex_path('repo.crires_plus_solar_y_rya794')))
#: Air wavelengths of Fe I lines RYA-794 identified in this product, with the measured
#: offset it found. Used as the load-time frame control (property 3). These are the
#: ticket's own published numbers, not new ones.
_CRIRES_Y_FRAME_CONTROL = (10535.709, 10577.139, 10611.686, 10616.721, 10674.070)
_CRIRES_Y_FRAME_TOL_A = 0.10        # RYA-794 measured max |offset| 0.039 A


def crires_y_spectrum() -> tuple[np.ndarray, np.ndarray]:
    """The RYA-794 corrected Y product as (wave_air_A, flux_normalised), ascending."""
    if "y" in _crires_cache:
        return _crires_cache["y"]
    if not CRIRES_Y_CSV.exists():
        raise LookupError(
            f"the RYA-794 corrected CRIRES+ Y product is not at {CRIRES_Y_CSV}. It is "
            f"git-tracked under data/results/rya794/ and registered as path key "
            f"'repo.crires_plus_solar_y_rya794'; rebuild it with "
            f"scripts/normalize_vesta_ir.py rather than falling back to the raw IDPs, "
            f"which are telluric-uncorrected and refused by the gate.")
    d = pd.read_csv(CRIRES_Y_CSV)
    missing = {"wavelength_air_A", "flux_normalized"} - set(d.columns)
    if missing:
        raise LookupError(f"{CRIRES_Y_CSV.name} is missing column(s) {sorted(missing)}; "
                          f"refusing to guess which column is the wavelength.")
    w = np.asarray(d["wavelength_air_A"], dtype=float)
    f = np.asarray(d["flux_normalized"], dtype=float)
    o = np.argsort(w)
    w, f = w[o], f[o]
    _assert_air_rest_frame(w, f)
    _crires_cache["y"] = (w, f)
    return _crires_cache["y"]


def _assert_air_rest_frame(w: np.ndarray, f: np.ndarray) -> None:
    """Refuse the product if its wavelength solution is not the AIR, rest-frame one.

    An ABSENCE-SHAPED failure otherwise (RYA-833): a vacuum or topocentric version of
    this file has the same column name, the same row count and the same flux, and every
    line simply fits at the wrong place. So the check is POSITIVE -- find the deepest
    point near each of RYA-794's published Fe I positions and require it to be close.
    ~2.9 A separates air from vacuum at 10500 A, and the tolerance here is 0.10 A, so the
    test discriminates by a factor of ~29.
    """
    off = []
    for c in _CRIRES_Y_FRAME_CONTROL:
        m = np.abs(w - c) <= 0.25
        if m.sum() < 5:
            raise LookupError(
                f"{CRIRES_Y_CSV.name}: the frame control line {c:.3f} A has only "
                f"{int(m.sum())} points near it. The product does not span what RYA-794 "
                f"says it spans ({w.min():.1f}-{w.max():.1f} A) -- refusing to measure.")
        off.append(float(w[m][int(np.argmin(f[m]))] - c))
    worst = max(abs(o) for o in off)
    if worst > _CRIRES_Y_FRAME_TOL_A:
        raise LookupError(
            f"{CRIRES_Y_CSV.name}: RYA-794's Fe I features sit {worst:.3f} A from their "
            f"AIR wavelengths (offsets {['%+.3f' % o for o in off]}), against a "
            f"{_CRIRES_Y_FRAME_TOL_A:.2f} A tolerance. RYA-794 measured at most 0.039 A. "
            f"CRIRES+ delivers VACUUM/TOPOCENT natively, so this is what a product that "
            f"lost its air conversion or its RV correction looks like -- and it would "
            f"still fit lines, at the wrong abundance. Refusing.")


def load_crires_y_window(centre: float, pad: float) -> tuple[np.ndarray, np.ndarray, str]:
    w, f = crires_y_spectrum()
    w, f = _slice_window(w, f, centre, pad, "CRIRES+ Y (RYA-794 telluric-corrected)")
    return w, f, CRIRES_Y_CSV.name


# ── HOLDINGS, not instruments — RYA-904 ──────────────────────────────────────
#
# 🔴 THE DEFECT THIS REPLACES. Two dicts were keyed by INSTRUMENT: `_LOADER_HOLDING`
# (which holding an arm serves) and `PRE_NORMALISED` (whether it arrives normalised).
# One instrument therefore got one answer to each. `crires_plus` has THREE solar
# holdings and they disagree on BOTH axes:
#
#   solar_crires_plus_y_rya794    telluric applied      normalised      MAY RUN
#   elgueta2026_vizier            telluric applied      (upstream)      MAY RUN
#   solar_vesta_crires_plus_idp   telluric not-applied  un-normalised   REFUSED
#
# The map named the third. So `load_window("crires_plus", ...)` asked for the one holding
# that should be refused, got refused -- correctly -- and the two corrected holdings were
# unreachable, not by any telluric decision but because NOTHING COULD NAME THEM. The gate
# was right the whole time; the dispatch one level up was wrong.
#
# ⚠️ BOTH KEYS MOVE TOGETHER OR NEITHER DOES. Wiring the corrected Y arm while leaving
# normalisation instrument-keyed would set a continuum on an already-normalised spectrum
# -- the RYA-713 double-normalisation defect, EWs low by a median 11.7 %. `PRE_NORMALISED`
# is now DERIVED from this table rather than written a second time, so the pair cannot
# come apart again.
@dataclass(frozen=True)
class HoldingSpec:
    """One HOLDING an instrument can serve a window from.

    `span_A` is the holding's own extent where it is a fixed product (the RYA-794 Y arm
    is exactly 10280-10680 A) and None where the reader inventories its own coverage
    (the Kitt Peak segment list, the CRIRES+ IDP comb). Coverage is required to be TOTAL,
    not overlapping: a window half inside the corrected product is a truncated window,
    and quietly serving it would be worse than falling through to the next candidate and
    saying so.
    """
    holding_id: str
    reader: str
    pre_normalised: bool
    span_A: tuple[float, float] | None = None
    caveat: str = ""
    note: str = ""
    #: RYA-911 — does the SOURCE PRODUCT carry a continuum of its own that the harness's
    #: placement can be read against? Only true where the product ships that column. It
    #: is a REFERENCE, never an authority: two continua disagreeing is evidence, and
    #: which one is right is decided by measurement, not by which is written down first.
    reference_continuum: bool = False

    def covers(self, centre: float, pad: float) -> bool:
        if self.span_A is None:
            return True
        return self.span_A[0] <= centre - pad and centre + pad <= self.span_A[1]


#: RYA-794's finding, carried to every product measured off that holding (RYA-904 spec 6).
#: It is a property of the DATA, so it travels with the holding rather than being retyped
#: by each consumer -- a caveat that lives only in a ticket does not reach the page.
GDSAT_CAVEAT = (
    "ROBUSTNESS CAVEAT (RYA-794): Elgueta+2026's own G-dwarf saturation flag GDSat=Y is "
    "set on NO Fe I line anywhere in atomicy.dat -- 42 lines are certified robust for a "
    "solar-type star and every one is another species. Fe I in the Y band is high-EP and "
    "only strengthens in cooler stars, which is the physical reason. RYA-794 therefore "
    "quoted NO ABUNDANCE from this arm. Reaching and measuring the pool does not "
    "pre-decide that it is publishable.")

#: Instrument -> its holdings IN PREFERENCE ORDER. Selection takes the first candidate
#: that (a) covers the window and (b) PASSES `gate_holding()`. It never falls back to a
#: refused holding: a refusal is reported with every candidate's reason, so "we hold this
#: window in no state that may be measured" is distinguishable from "we do not hold it".
_INSTRUMENT_HOLDINGS: dict[str, tuple[HoldingSpec, ...]] = {
    "harps": (
        HoldingSpec("solar_harps", reader="harps", pre_normalised=True,
                    reference_continuum=True,
                    # RYA-767 -- DECLARE the extent, do not leave it to be discovered.
                    # The span is already VERIFIED in solar_reference_holdings_rya708.csv
                    # (3782.6-6910.0) and was simply not carried here, so `covers()`
                    # answered True for every window ever asked -- including the NIR.
                    # A planner enumerating holdings x bands therefore proposed eight
                    # HARPS runs above the 6912 A cutoff, each of which would have failed
                    # one line at a time deep inside a run instead of being ruled out
                    # before it started.
                    span_A=(3782.6, 6910.0),
                    note="Direct-solar HARPS (Dumusque ESO 1102.D-0954), RYA-897, "
                         "continuum contract FIXED by RYA-911. PRE-NORMALISED: the "
                         "product ships its own fitted continuum and we consume it. "
                         "🔴 It was pre_normalised=False and re-fitting locally, which "
                         "put the continuum BELOW the observed flux on 8 of 16 lines "
                         "(worst max(F/C) 1.204 -- impossible) and cost a median 23.8% "
                         "of the EW in a same-inputs control. That was the -0.34 dex."),
        HoldingSpec("solar_harps_molecfit_corrected", reader="harps_tellcorr",
                    pre_normalised=True, reference_continuum=True,
                    span_A=(3782.6, 6910.0),
                    note="RYA-931 telluric-corrected sibling: the same ten exposures "
                         "through the same normaliser, with the O2 B band divided out "
                         "per exposure by its own molecfit/GDAS transmission. Listed "
                         "SECOND ON PURPOSE -- it is reachable by name, but selection "
                         "order is unchanged, so no existing measurement silently "
                         "switches product. Choosing it for the affected red-edge "
                         "windows is RYA-936's decision to make and record, not a "
                         "side effect of this wiring. 256 saturated-core pixels are "
                         "NaN (quarantined, not divided), so a window overlapping them "
                         "loses those pixels rather than being served a fabricated "
                         "flux."),
    ),
    "kpno_solar_atlas": (
        HoldingSpec("solar_kpno", reader="kpno", pre_normalised=True,
                    note="Kurucz/Brault FTS residual flux -- unity IS the continuum. "
                         "🔴 STANDING RULE, RYA-1026: DO NOT NORMALISE ANY KITT PEAK "
                         "ATLAS. Both KP products ship their own continuum and a second "
                         "one adds a spurious TILT that follows the saturated bands "
                         "down. It has bitten twice -- RYA-940 here, RYA-929 on the 2005 "
                         "sibling. Enforced by pipeline.prenormalised_guard. The ONLY "
                         "thing done to a KP atlas on the way in is TELLURIC "
                         "CORRECTION; never a continuum refit."),
        HoldingSpec("solar_kpno_molecfit_corrected", reader="kpno_1984_composite",
                    pre_normalised=True,
                    note="RYA-933: serves the WHOLE band -- corrected flux inside the "
                         "six RYA-940 bands, the untouched 1984 atlas outside them, and "
                         "a REFUSAL inside a registered telluric band RYA-940 could not "
                         "fit (H2O 7160-7340). Until then this reader returned only the "
                         "six corrected windows and nothing else, so every graded VIS "
                         "line fell outside it and the band measured zero lines -- the "
                         "holding had never produced a single product. "
                         "RYA-940 telluric-corrected 1984 atlas. Same conventions as "
                         "solar_kpno -- air, residual flux, unity IS the continuum -- "
                         "differing only by the six corrected telluric bands. Listed "
                         "AFTER solar_kpno on purpose: reachable by name, selection "
                         "order unchanged, so no existing measurement silently switches "
                         "product. NaN marks a quarantined saturated core."),
        HoldingSpec("solar_kpno_kurucz2005_corrected", reader="kurucz2005",
                    pre_normalised=True, span_A=(2990.0, 10010.0),
                    note="Kurucz 2005, telluric-corrected at source. Served from "
                         "irradrelwl.dat, the RESIDUAL atlas Kurucz ships alongside the "
                         "irradiance file -- so this holding DOES carry its own continuum "
                         "and pre_normalised is True. It was False until RYA-933 because "
                         "`0irrad.readme` does not list the residual file and the RYA-929 "
                         "intake never took it; placing our own continuum instead tilted "
                         "the band 4% blue-to-red and cost 0.022 dex. VACUUM grid "
                         "(gravitational redshift included), converted to air on read. "
                         "Spans 2990-10010 A; nothing telluric-free reaches the IR. "
                         "🔴 RYA-1026 ratified this for the WHOLE KITT PEAK CLASS and "
                         "made it ENFORCED, not remembered: pipeline.prenormalised_guard "
                         "refuses to fit/apply a continuum here, and cross-checks the "
                         "flag against the FLUX -- a declared flag and a mis-routed file "
                         "agree with each other perfectly and are both wrong, which is "
                         "exactly what happened here for months."),
    ),
    "iag_fts_solar_atlas": (
        HoldingSpec("solar_iag", reader="iag", pre_normalised=True,
                    span_A=(5001.10, 11083.46),   # AIR (vac 5002.5-11086.5)
                    note="Baker+2020 TELLURIC-CORRECTED atlas, normalised. 🔴 THE SPAN IS "
                         "DECLARED NOW (RYA-767): it was absent, so covers() answered True "
                         "for every window ever asked -- including the 954 A of blue this "
                         "atlas does not reach at all. MEASURED telluric state on the "
                         "served flux: O2 A-band 0.18% of pixels below 0.5, H2O 9280-9600 "
                         "0.00% -- against 46.25%/51.63% for the raw Reiners sibling, a "
                         "250x separation. This holding IS corrected; the long-standing "
                         "doubt was a stale CATALOGUE ROW pointing at the other file, "
                         "never the reader."),
        HoldingSpec("solar_iag_reiners2016", reader="iag_reiners", pre_normalised=True,
                    span_A=(4047.46, IAG_BAKER_BLUE_EDGE_A),   # AIR (vac 4048.6-)
                    note="Reiners+2016 IAG, THE BLUE ARM. Listed second and SPAN-CAPPED at "
                         "5001.1 A on purpose: where Baker+2020 reaches, the corrected "
                         "atlas wins, so this can never be selected in preference to it. "
                         "Below that it is the ONLY IAG option -- Baker starts at 5001.1 A "
                         "and cannot substitute. NOT A COVERAGE GAP: together the two span "
                         "4047.5-11083.5 A (air). 🔴 AND THERE IS NOTHING TO CORRECT HERE -- the "
                         "bluest registered telluric complex is O2 B at 6867 A, and "
                         "measured against KP2005's telluric-free-at-source residual atlas "
                         "over 4100-5000 A the two agree (10.48% vs 10.35% of pixels below "
                         "0.5) with Reiners not systematically deeper. 5.6x denser "
                         "sampling than Baker (4.06M vs 728K points)."),
    ),
    "crires_plus": (
        HoldingSpec("solar_crires_plus_y_rya794", reader="crires_y", pre_normalised=True,
                    span_A=(10280.0, 10680.0), caveat=GDSAT_CAVEAT,
                    note="RYA-794 science-ready Y arm: telluric-corrected (measured), "
                         "continuum-normalised, air, solar rest frame. PREFERRED over "
                         "the raw IDPs wherever it covers the window."),
        HoldingSpec("solar_vesta_crires_plus_idp", reader="crires_idp",
                    pre_normalised=False,
                    note="Raw Vesta IDPs: adu, un-normalised, TOPOCENT, telluric "
                         "not-applied. Correctly refused by the gate for measurement; "
                         "reachable only by the correction/conditioning legs."),
    ),
}

#: DERIVED, never written twice (RYA-845's defect shape: two declarations of one fact).
PRE_NORMALISED: dict[str, bool] = {
    h.holding_id: h.pre_normalised
    for specs in _INSTRUMENT_HOLDINGS.values() for h in specs}


class TelluricNotCorrected(LookupError):
    """We hold this window, but not in a telluric state that may be measured.

    A LookupError subclass so a driver's existing "this arm cannot serve this line"
    handling still catches it -- but a distinct type, because "we do not hold this
    wavelength", "we hold it in the wrong velocity frame" (`RestFrameNotConditioned`) and
    "we hold it uncorrected" are three different problems with three different fixes.
    """


class Window(NamedTuple):
    """A loaded window AND the holding it came from — RYA-904.

    The holding travels WITH the data. Before this, a caller asked `load_window` for the
    flux and a separate dict for the normalisation, which is how those two could describe
    different products. Anything a caller needs to know about how this spectrum was made
    is on `.holding`, which is the object selection actually used.
    """
    wave: np.ndarray
    flux: np.ndarray
    provenance: str
    holding: HoldingSpec

    @property
    def pre_normalised(self) -> bool:
        return self.holding.pre_normalised


def holdings_for(instrument: str) -> tuple[HoldingSpec, ...]:
    """Every holding wired for an instrument, in preference order. Loud on an unknown."""
    try:
        return _INSTRUMENT_HOLDINGS[instrument]
    except KeyError:
        raise LookupError(
            f"no window loader for instrument {instrument!r}. Add its holdings to "
            f"_INSTRUMENT_HOLDINGS here rather than letting a driver fall back to "
            f"another arm's data.") from None


def select_holding(instrument: str, centre: float, pad: float, *,
                   allow_uncorrected: bool = False,
                   holding: str | None = None) -> HoldingSpec:
    """Which HOLDING serves this window — the dispatch RYA-904 exists to fix.

    Order: cover the window, then pass `gate_holding()`, then first-wins. A candidate that
    covers but is refused is REPORTED, never silently skipped past and never fallen back
    to; if nothing passes, the raise carries every candidate's own reason so the caller
    can tell a coverage answer from a telluric one.

    `holding=` names one explicitly and is how a correction/conditioning leg asks for the
    product it is about to correct -- and how the RYA-904 control points the instrument
    back at only the raw IDPs. `allow_uncorrected=True` suspends the GATE, not the
    preference order; it exists for the leg that must read uncorrected flux in order to
    correct it (mirrors `allow_topocentric` on the rest-frame gate) and is not a general
    escape hatch.
    """
    specs = holdings_for(instrument)
    if holding is not None:
        named = [h for h in specs if h.holding_id == holding]
        if not named:
            raise LookupError(
                f"holding {holding!r} is not wired for {instrument!r}; wired: "
                f"{[h.holding_id for h in specs]}")
        specs = tuple(named)

    covering = [h for h in specs if h.covers(centre, pad)]
    if not covering:
        spans = "; ".join(
            f"{h.holding_id} {h.span_A[0]:.1f}-{h.span_A[1]:.1f} A" for h in specs
            if h.span_A is not None)
        raise LookupError(
            f"no {instrument} holding covers {centre:.3f} +/- {pad:.3f} A ({spans}). "
            f"Coverage must be TOTAL -- a window half inside a fixed product is a "
            f"truncated window, not a measurement.")

    from pipeline.telluric_policy import gate_holding
    refusals, unknown = [], None
    for h in covering:
        if allow_uncorrected:
            return h
        try:
            ok, why = gate_holding(h.holding_id, instrument)
        except Exception as e:           # TelluricStateUnknown and friends
            unknown = unknown or e
            refusals.append(f"{h.holding_id}: {type(e).__name__}: {e}")
            continue
        if ok:
            return h
        refusals.append(why)
    if unknown is not None:
        raise unknown
    raise TelluricNotCorrected(
        f"{instrument} covers {centre:.3f} A but no holding may be measured there. "
        + " | ".join(refusals)
        + " If you ARE the telluric correction leg and need the uncorrected flux in "
          "order to correct it, pass allow_uncorrected=True.")


#: holding.reader -> the function that actually reads it. Split from the specs so the
#: preference table stays readable and a reader can be shared by two holdings.
def _reader(spec: HoldingSpec, centre: float, pad: float, segs):
    if spec.reader == "kpno":
        return load_kp_window(segs if segs is not None else kp_segments(), centre, pad)
    if spec.reader == "kpno_1984_corrected":
        return load_kp1984_corrected_window(centre, pad)
    if spec.reader == "kpno_1984_composite":
        return load_kp1984_composite_window(centre, pad, segs)
    if spec.reader == "kurucz2005":
        return load_kurucz2005_window(centre, pad)
    if spec.reader == "iag":
        return load_iag_window(centre, pad)
    if spec.reader == "iag_reiners":
        return load_iag_reiners_window(centre, pad)
    if spec.reader == "harps":
        return load_harps_window(centre, pad)
    if spec.reader == "harps_tellcorr":
        return load_harps_window(centre, pad, HARPS_TELLCORR_CSV)
    if spec.reader == "crires_y":
        return load_crires_y_window(centre, pad)
    if spec.reader == "crires_idp":
        return load_crires_window(centre, pad)
    raise LookupError(f"holding {spec.holding_id} names reader {spec.reader!r}, which "
                      f"is not implemented here.")


def load_window_ex(instrument: str, centre: float, pad: float, segs=None,
                   allow_uncorrected: bool = False,
                   holding: str | None = None) -> Window:
    """One window, WITH the holding that served it — RYA-904.

    Loud on an unknown instrument: silently defaulting to Kitt Peak is how a product gets
    labelled with an instrument it was not measured on. Loud, too, on a telluric state
    that forbids measurement (RYA-806) -- checked BEFORE any data is read, so a refusal
    costs nothing and cannot be half-completed.

    ⚠️ THE PROVENANCE NAMES THE HOLDING, not just the instrument. A CRIRES+ product could
    have come from either the corrected Y arm or the raw IDPs, and "crires_plus" does not
    say which. That was true before this ticket too -- the string simply could not have
    been wrong, because only one holding was reachable.
    """
    spec = select_holding(instrument, centre, pad,
                          allow_uncorrected=allow_uncorrected, holding=holding)
    w, f, prov = _reader(spec, centre, pad, segs)
    return Window(w, f, f"holding={spec.holding_id} · {prov}", spec)


def load_window(instrument: str, centre: float, pad: float, segs=None,
                allow_uncorrected: bool = False, holding: str | None = None):
    """`load_window_ex` as the historical (wave, flux, provenance) triple."""
    win = load_window_ex(instrument, centre, pad, segs,
                         allow_uncorrected=allow_uncorrected, holding=holding)
    return win.wave, win.flux, win.provenance



def kp_segments(allow_corrupt: bool = False) -> list[tuple[float, float, Path]]:
    """Inventory the atlas as (lo_A, hi_A, path). Reads each file's ACTUAL span rather
    than trusting the filename -- the lm#### stem is a start hint, not a guarantee.

    RYA-938 -- THIS USED TO SWALLOW A PARSE FAILURE. The body was `except Exception:
    continue`, so a segment that would not parse was dropped from the inventory and the
    next question answered "no Kitt Peak segment covers 8420.000 A". A CORRUPT FILE
    PRESENTED AS MISSING COVERAGE, which is the RYA-833 shape and exactly the failure
    `_resolve_kp_dir` was written to prevent one level up. It was not hypothetical:
    `lm0840` in both staged copies was a saved HTTP 500 page, hiding 8400-8441 A --
    where the solar line list holds 33 Fe lines alone.

    `allow_corrupt=True` is the deliberate escape for a caller that wants degraded
    operation; it still returns only good segments, but the caller has SAID SO.
    """
    reports = [kp_integrity.inspect_segment(p)
               for p in sorted(KP_DIR.glob("lm[0-9]*")) if p.is_file()]
    if not allow_corrupt:
        kp_integrity.require_parseable(reports)
    return [(r.lo_A, r.hi_A, Path(r.path)) for r in reports if r.ok]


#: RYA-940's corrected 1984 products, and the Kurucz 2005 reference. Both are
#: HOLDINGS of the kpno_solar_atlas instrument, and until now neither could be
#: NAMED by the harness -- the RYA-904 shape, where a holding nobody can reach
#: reads to every caller exactly like having no data.
KP1984_CORRECTED_DIR = codex_root('repo') / 'data' / 'processed' / 'kp1984_telluric_corrected'


def load_kp1984_corrected_window(centre: float, pad: float):
    """RYA-940 telluric-corrected 1984 segments. Air, residual flux, NaN = quarantined.

    Only the six corrected bands exist; everywhere else this holding simply does not
    reach, and says so rather than falling back to the uncorrected atlas -- which would
    label an uncorrected measurement 'corrected'.
    """
    lo, hi = centre - pad, centre + pad
    W, F, used = [], [], []
    for path in sorted(KP1984_CORRECTED_DIR.glob("kp1984_corrected_*.txt")):
        a, b = (float(x) for x in path.stem.split("_")[-2:])
        if b < lo or a > hi:
            continue
        arr = np.loadtxt(path)
        m = (arr[:, 0] >= lo) & (arr[:, 0] <= hi) & np.isfinite(arr[:, 1])
        if m.any():
            W.append(arr[m, 0]); F.append(arr[m, 1]); used.append(path.name)
    if not W:
        raise LookupError(
            f"no RYA-940 corrected 1984 band covers {centre:.3f} A. Only six telluric "
            f"bands were corrected; this window is not one of them. Use solar_kpno for "
            f"the uncorrected atlas -- do NOT relabel it corrected.")
    w = np.concatenate(W); f = np.concatenate(F)
    o = np.argsort(w)
    return w[o], f[o], ",".join(used)


def load_kp1984_composite_window(centre: float, pad: float, segs=None):
    """RYA-933: the 1984 atlas with RYA-940's corrected bands SUBSTITUTED IN.

    `solar_kpno_molecfit_corrected` holds ONLY the six fitted telluric windows, so every
    line outside them is unserved and a VIS band measures nothing. This composes the two
    without relabelling anything:

      * inside a band RYA-940 corrected -> the corrected flux;
      * outside every REGISTERED telluric band -> the original 1984 atlas, which is
        telluric-clean there by the project's own enumeration (`TELLURIC_BANDS`);
      * inside a registered band with NO corrected file -> REFUSE.

    🔴 THE REFUSAL IS THE POINT, and it is why this is not the silent fallback
    `load_kp1984_corrected_window` rightly forbids. H2O 7160-7340 is registered and
    RYA-940 got NO ADMISSIBLE FIT for it, so a window there is genuinely uncorrected and
    must not be served under a corrected name. Falling through everywhere would do
    exactly that; falling through only where nothing needs correcting does not.
    """
    from pipeline.telluric_policy import TELLURIC_BANDS
    lo, hi = centre - pad, centre + pad
    corrected = []
    for path in sorted(KP1984_CORRECTED_DIR.glob("kp1984_corrected_*.txt")):
        a, b = (float(x) for x in path.stem.split("_")[-2:])
        corrected.append((a, b, path))

    # a registered band this window touches, for which no corrected file exists -> refuse
    for blo, bhi, bname in TELLURIC_BANDS:
        if bhi < lo or blo > hi:
            continue
        if not any(a <= bhi and b >= blo for a, b, _ in corrected):
            raise LookupError(
                f"window {lo:.2f}-{hi:.2f} A overlaps REGISTERED telluric band "
                f"{blo:.0f}-{bhi:.0f} A ({bname}), and RYA-940 produced no admissible "
                f"correction for it. Refusing to serve uncorrected flux under a "
                f"corrected holding -- measure it on solar_kpno and label it uncorrected.")

    W, F, used = [], [], []
    for a, b, path in corrected:                       # corrected flux where it exists
        if b < lo or a > hi:
            continue
        arr = np.loadtxt(path)
        m = (arr[:, 0] >= lo) & (arr[:, 0] <= hi) & np.isfinite(arr[:, 1])
        if m.any():
            W.append(arr[m, 0]); F.append(arr[m, 1]); used.append(path.name)

    aw, af, aprov = load_kp_window(segs if segs is not None else kp_segments(), centre, pad)
    keep = np.ones(np.asarray(aw).size, bool)          # atlas flux everywhere else
    for a, b, _ in corrected:
        keep &= ~((aw >= a) & (aw <= b))
    if keep.any():
        W.append(np.asarray(aw)[keep]); F.append(np.asarray(af)[keep])
        used.append(f"{aprov}[uncorrected: no telluric band registered]")
    if not W:
        raise LookupError(f"no 1984 atlas or corrected coverage at {centre:.3f} A")
    w = np.concatenate(W); f = np.concatenate(F)
    o = np.argsort(w)
    return w[o], f[o], ",".join(used)


def _read_kurucz2005_residual(path: Path, lo: float, hi: float):
    """(air_A, residual_flux) from irradrelwl.dat. Vacuum nm in, air Angstrom out."""
    from pipeline.uv_conditioning import vac_to_air
    w, f = [], []
    with Path(path).open(errors="replace") as fh:
        for line in fh:
            q = line.split()
            if len(q) < 2:
                continue
            try:
                x, y = float(q[0]) * 10.0, float(q[1])
            except ValueError:
                continue                      # header lines
            if lo - 5.0 <= x <= hi + 5.0:
                w.append(x); f.append(y)
    w, f = np.asarray(w), np.asarray(f)
    if not w.size:
        return w, f
    return np.asarray(vac_to_air(w)), f


def load_kurucz2005_window(centre: float, pad: float):
    """Kurucz 2005 RESIDUAL irradiance atlas -- the continuum Kurucz shipped.

    🔴 IT DOES SHIP A CONTINUUM, and reading the wrong file cost 0.022 dex. The
    distribution carries `irradthu.dat` (absolute irradiance, W/m2/nm) AND
    `irradrelwl.dat`, the "KURUCZ RESIDUAL IRRADIANCE ATLAS 2005" -- the same spectrum
    with Kurucz's own continuum divided out. **`0irrad.readme` does not list the residual
    file**, so the RYA-929 intake took only the flux file and this holding was recorded as
    shipping no continuum. Placing our own instead put a 4% blue-to-red tilt on the band
    (1.0238 at 4400 A falling to 0.9848 at 6800 A, measured against the 1984 atlas) and
    biased A(Fe I) low by 0.0218 +/- 0.0040 dex, correlated with wavelength (r=+0.373).
    The shipped continuum tilts 12x less: -0.0033 across the same band. This is exactly
    the RYA-911/938 rule -- do not re-fit a continuum where the product ships one.

    Still VACUUM ("vacuum wavelength including gravitational red shift"), so the same
    vac->air conversion applies; verified at 0.000 A shift against the air 1984 atlas
    (r=0.9997). Residual flux, so the holding is `pre_normalised=True` and needs no
    `--place-continuum`.
    """
    path = codex_path('data.kurucz2005_residual')
    w, f = _read_kurucz2005_residual(Path(path), centre - pad, centre + pad)
    if w.size < 5:
        raise LookupError(f"Kurucz 2005 does not cover {centre:.3f} A (it spans "
                          f"2990-10010 A); nothing beyond 10010 A is telluric-free here.")
    assert_data_matches_declaration(
        "solar_kpno_kurucz2005_corrected", w, f, declared=True,
        where=f"load_kurucz2005_window({centre:.3f} A) -> {Path(path).name}")
    return w, f, Path(path).name


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
            # RYA-904 — through the holding dispatch, so the normalisation flag comes
            # from the SAME object that served the flux. Same reader, same data as the
            # direct `load_kp_window` call this replaces; what changes is that `src` now
            # names the holding and `pre_normalised` can no longer describe a different
            # product from the one measured.
            win = load_window_ex(a.instrument, float(r.wave_air_A), hw * 3.0, segs)
            w, f, src = win.wave, win.flux, win.provenance
            # Kitt Peak ships residual flux -- already normalised. Say so explicitly
            # rather than re-normalising on top of it (see band_products docstring).
            ew, method, concern = equivalent_width(
                w, f, float(r.wave_air_A), hw, pre_normalised=win.pre_normalised)
            ok, why = verify_feature(w, f, float(r.wave_air_A), hw,
                                     float(r.predicted_depth))
        except Exception as e:
            skipped.append(dict(wave=r.wave_air_A, reason=f"{type(e).__name__}: {e}"))
            continue
        lm = LineMeasurement(
            element=a.element, ion=a.ion, wavelength_air_A=float(r.wave_air_A),
            instrument=a.instrument, ew_mA=ew,
            # RYA-871 — the same carry as the profile-fit driver. This harness is
            # diagnostic-only for the EW VALUE, but the line IDENTITY it writes is read
            # by the same consumers, so it must be identifiable too.
            ep_eV=carried_ep(r, wavelength_A=float(r.wave_air_A), element=a.element,
                             ion=a.ion),
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
    if a.max_lines:
        # A SUBSET run never writes over a full one. A --max-lines smoke test
        # silently clobbered a documented 445-line result; the fixed output path
        # made a partial run indistinguishable from a complete one on disk.
        stem += f"_SUBSET{a.max_lines}"
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
