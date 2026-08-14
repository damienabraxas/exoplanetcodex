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
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline.band_products import (  # noqa: E402
    LineMeasurement, equivalent_width, assert_single_element)
from pipeline.band_policy import check_intake, resolve, BandPolicyError  # noqa: E402
from config.constants import codex_path  # RYA-810 path register

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
    "/mnt/codex-data/spectra/Solar Calibration/Kitt Peak Flux Atlas",
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
PRE_NORMALISED = {"kpno_solar_atlas": True, "harps": False, "iag_fts_solar_atlas": True,
                  "crires_plus": False}

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

IAG_FITS = Path("/srv/codex/solar_reference/iag_baker2020/iag_telfree_solaratlas.fits")
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
    w, f = _slice_window(w, f, centre, pad, "IAG")
    return w, f, IAG_FITS.name


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
    "/mnt/codex-data/spectra/vesta/CRIRESPlus",
)
#: Rest-frame conditioned output of pipeline/reflected_solar_rv.py `write_set`.
_CRIRES_REST_CANDIDATES = (
    os.environ.get("CODEX_CRIRES_VESTA_REST", ""),
    "/mnt/codex-data/spectra/vesta/CRIRESPlus_rest/vesta_crires",
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


#: instrument -> the holding each arm actually serves (RYA-806). The telluric gate is a
#: per-HOLDING question and this loader dispatches per INSTRUMENT, so the two are joined
#: here, once, rather than by every caller. An instrument absent from this map is gated
#: on its instrument axis alone.
_LOADER_HOLDING = {
    "kpno_solar_atlas": "solar_kpno",
    "iag_fts_solar_atlas": "solar_iag",
    # load_crires_window serves the staged Vesta IDPs specifically -- not the Elgueta
    # reduced spectra, which are a different holding at a different product level.
    "crires_plus": "solar_vesta_crires_plus_idp",
}


def _assert_telluric_state(instrument: str, allow_uncorrected: bool = False) -> None:
    """Refuse a window whose telluric state forbids it (RYA-806).

    THE GAP THIS CLOSES. RYA-805 found this module consuming only `TELLURIC_BANDS` and
    `exclusion()` while never calling the policy gate — and `exclusion()` returns '' for
    every CRIRES+ wavelength because the enumerated band set stops at 11560 A, so the
    whole J/H/K arm fell off the end of the list and every IR line read as clean. The
    only thing refusing the arm was the REST-FRAME gate, which would stop refusing the
    moment RYA-372/373 conditioning ran, leaving telluric-uncorrected H-band flux to be
    measured with no telluric objection at all.

    ⚠️ THIS FIRES BEFORE THE REST-FRAME GATE, and the order is the physics, not an
    accident: tellurics are stationary in the TOPOCENTRIC frame, so the correction must
    happen there and the RV shift comes after (RYA-373). Telluric is the earlier blocker,
    so it is the one a caller is told about first. This does change what RYA-796's
    `load_window` raises for the staged IDPs — `TelluricNotCorrected` now, where it was
    `RestFrameNotConditioned` — because a second, earlier defect was found in the same
    data. Both refusals remain reachable and both are still tested.

    `allow_uncorrected=True` is for the CORRECTION LEG ITSELF. RYA-373's molecfit driver
    has to read uncorrected flux in order to correct it, so a gate with no door would
    lock out the only thing that can clear it. Exactly mirrors `allow_topocentric` on the
    rest-frame gate, and is equally not a general escape hatch.
    """
    if allow_uncorrected:
        return
    holding = _LOADER_HOLDING.get(instrument)
    if holding is None:
        return
    from pipeline.telluric_policy import gate_holding
    ok, why = gate_holding(holding, instrument)
    if not ok:
        raise TelluricNotCorrected(
            f"{why} If you ARE the telluric correction leg and need the uncorrected "
            f"flux in order to correct it, pass allow_uncorrected=True.")


class TelluricNotCorrected(LookupError):
    """We hold this window, but not in a telluric state that may be measured.

    A LookupError subclass so a driver's existing "this arm cannot serve this line"
    handling still catches it -- but a distinct type, because "we do not hold this
    wavelength", "we hold it in the wrong velocity frame" (`RestFrameNotConditioned`) and
    "we hold it uncorrected" are three different problems with three different fixes.
    """


def load_window(instrument: str, centre: float, pad: float, segs=None,
                allow_uncorrected: bool = False):
    """One entry point per instrument, so a driver does not hardcode an arm.

    Loud on an unknown instrument: silently defaulting to Kitt Peak is how a product gets
    labelled with an instrument it was not measured on. Loud, too, on a telluric state
    that forbids measurement (RYA-806) -- checked BEFORE any data is read, so a refusal
    costs nothing and cannot be half-completed.
    """
    _assert_telluric_state(instrument, allow_uncorrected)
    if instrument == "kpno_solar_atlas":
        return load_kp_window(segs if segs is not None else kp_segments(), centre, pad)
    if instrument == "iag_fts_solar_atlas":
        return load_iag_window(centre, pad)
    if instrument == "crires_plus":
        return load_crires_window(centre, pad)
    raise LookupError(
        f"no window loader for instrument {instrument!r}. Add one here rather than "
        f"letting a driver fall back to another arm's data.")


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
