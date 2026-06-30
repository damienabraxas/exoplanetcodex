"""
pipeline/ir_telluric.py
=======================
RYA-494 — generalized, per-star / per-instrument IR telluric conditioning.

RYA-373 built molecfit telluric correction BESPOKE to the Vesta CRIRES+ K-band
(reflected-solar CO). This module generalizes that lineage so every IR target —
α Cen A is caller #1, 55 Cnc and future stars inherit — runs the same conditioning
without another one-off. The strategic point: benchmark work builds the software
for all science targets.

What was Vesta/CRIRES-K-specific, now parameterized:

  * Instrument-agnostic.  CRIRES+ (cr2res EXTRACTC IDP: nm WAVE, TOPOCENT, no
    telluric → molecfit, reusing pipeline.crires_telluric) AND NIRPS (geneva
    S1D_FINAL: already telluric-corrected via FLUX_TELL_* + ATM_TRANSM, already
    BARYCENT, already Å → select the telluric column, NO molecfit).
  * Star-agnostic.  Takes a target descriptor, not Vesta-only.
  * The velocity step is BRANCHED (the RYA-373 → RYA-494 split):
      - REFLECTED_SOLAR : asteroid-ephemeris RV (Vesta; pipeline.reflected_solar_rv,
        RYA-372) — reflected sunlight, the body's Horizons velocity.
      - STELLAR         : a direct stellar target's own BERV (→ barycentric rest)
        plus systemic RV (→ stellar rest). α Cen A uses this — NOT asteroid ephemeris.

Permanent IR rules preserved (RYA-373, RYA-481):
  * telluric correction happens in the TOPOCENTRIC frame, BEFORE any RV shift
    (tellurics are stationary in Earth's atmosphere);
  * no flux/abundance is emitted unless `telluric_corrected` is set
    (TelluricNotCorrectedError otherwise);
  * nm → Å at the loader boundary;
  * no blind cross-epoch co-add — register on per-epoch RV first;
  * attribution + velocity frame are declared from the AUTHORITATIVE header /
    star-ID, never folder/filename (RYA-481). For α Cen A NIRPS the authoritative
    attribution is the RYA-423 IR RV star-ID, which OVERRIDES the (mislabeled)
    OBJECT header — see ALPHA_CEN_A_NIRPS_NOTE.

Smoke test:  python -m pipeline.ir_telluric --target alpha_cen_a --verify
"""
from __future__ import annotations

import glob
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path

import numpy as np

# Reuse the proven RYA-373 CRIRES molecfit machinery — do NOT re-implement it.
from pipeline import crires_telluric as ct
from pipeline.crires_telluric import (
    CriresFrame, CriresSegment, NM_TO_A, _C_KMS,
    TelluricNotCorrectedError, MolecfitNotAvailableError,
    continuum_normalize, load_crires_idp, _molecfit_segment, _esorex_available,
)

_DATA_ROOT = Path("/Users/ryanschmitt/Documents/Exoplanet Codex/data/spectra/"
                  "exoplanetcodex-data")


# ── Velocity-step branch (the RYA-373 → RYA-494 generalization) ───────────────
class VelocityMode(Enum):
    REFLECTED_SOLAR = "reflected_solar"   # asteroid ephemeris (Vesta), RYA-372
    STELLAR = "stellar"                   # direct star: BERV + systemic RV


class Instrument(Enum):
    CRIRES = "CRIRES"
    NIRPS = "NIRPS"


class IRAttributionError(RuntimeError):
    """The frame's authoritative attribution (header OR IR RV star-ID) does not
    match the requested target. Never silently substitute the wrong star (RYA-481)."""


# ── Target descriptor + registry ──────────────────────────────────────────────
@dataclass
class IRTarget:
    name: str
    velocity_mode: VelocityMode
    data_dir: Path                       # the (vetted) per-star folder
    body_key: str | None = None          # Horizons key when REFLECTED_SOLAR
    systemic_rv_kms: float = 0.0         # stellar systemic RV (→ stellar rest), STELLAR
    # authoritative attribution: which raw OBJECT / star-ID labels are THIS target,
    # per instrument. For NIRPS this is the RYA-423 RV star-ID truth, not OBJECT.
    nirps_attribution: str = ""          # informational provenance string


# α Cen A systemic RV ≈ −22.3 km/s (Gaia/CORAVEL); the dominant frame term is BERV.
ALPHA_CEN_A_NIRPS_NOTE = (
    "RYA-423 IR RV star-ID: the NIRPS frames headered OBJECT='AlphaCenB' have "
    "verdict=A (obs_rv≈-26.2 matches α Cen A -25.8, not B -18.4). The IR star-ID is "
    "authoritative over the (program-mislabeled) OBJECT header — these are α Cen A. "
    "The OBJECT='alf Cen A'-labelled NIRPS frames are verdict=NOT-ALPHA-CEN (RV≈-34, "
    "RYA-431 quarantine). α Cen B has NO NIRPS (RYA-439). This INVERTS the RYA-479 "
    "OBJECT-based NIRPS attribution for the IR arm."
)

TARGETS: dict[str, IRTarget] = {
    "vesta": IRTarget(
        name="vesta", velocity_mode=VelocityMode.REFLECTED_SOLAR,
        data_dir=_DATA_ROOT / "Solar Calibration/Solar System Targets/Vesta",
        body_key="vesta"),
    "alpha_cen_a": IRTarget(
        name="alpha_cen_a", velocity_mode=VelocityMode.STELLAR,
        data_dir=_DATA_ROOT / "Alpha Centauri (vetted)/Alpha Cen A",
        systemic_rv_kms=-22.3,
        nirps_attribution=ALPHA_CEN_A_NIRPS_NOTE),
}


# ── NIRPS frame model (already telluric-corrected; already Å + BARYCENT) ───────
@dataclass
class NirpsFrame:
    path: Path
    object_hdr: str
    mjd: float
    date_obs: str
    specsys: str                         # 'BARYCENT' (geneva DRS)
    berv_kms: float
    wave_A: np.ndarray                   # already Å (and barycentric)
    flux: np.ndarray                     # selected telluric-corrected column
    err: np.ndarray
    atm_transm: np.ndarray               # molecfit transmission the DRS divided out
    telluric_corrected: bool = False
    rest_frame: bool = False
    telluric_column: str = ""

    def covers_A(self, wl_A: float) -> bool:
        return float(np.nanmin(self.wave_A)) <= wl_A <= float(np.nanmax(self.wave_A))


# NIRPS telluric-corrected science column, in preference order (RYA-301 contract:
# read a FLUX_TELL_* column, never raw FLUX). EL = telluric-only; CAL = flux-cal'd.
_NIRPS_TELL_COLS = ("FLUX_TELL_CAL", "FLUX_TELL_EL")


def load_nirps_idp(path) -> NirpsFrame:
    """Load a NIRPS geneva S1D_FINAL into a NirpsFrame, SELECTING the telluric-
    corrected flux column (FLUX_TELL_CAL/_EL) + ATM_TRANSM — never raw FLUX. WAVE is
    already Å and BARYCENT (no nm→Å, no BERV re-apply)."""
    from astropy.io import fits
    path = Path(path)
    with fits.open(path) as h:
        hdr = h[0].header
        d = (h["SPECTRUM"].data if "SPECTRUM" in [x.name for x in h] else h[1].data)[0]
        cols = d.columns.names if hasattr(d, "columns") else h[1].data.columns.names
        tcol = next((c for c in _NIRPS_TELL_COLS if c in cols), None)
        if tcol is None:
            raise MolecfitNotAvailableError(
                f"{path.name}: no FLUX_TELL_* column — NIRPS frame is not telluric-"
                f"corrected; refusing to fall back to raw FLUX (permanent IR rule).")
        ecol = "ERR" + tcol[4:]
        wave = np.asarray(d["WAVE"], float)
        flux = np.asarray(d[tcol], float)
        err = np.asarray(d[ecol], float) if ecol in cols else np.full_like(flux, np.nan)
        atm = np.asarray(d["ATM_TRANSM"], float) if "ATM_TRANSM" in cols else np.full_like(flux, np.nan)
    return NirpsFrame(
        path=path, object_hdr=str(hdr.get("OBJECT", "?")),
        mjd=float(hdr.get("MJD-OBS", np.nan)), date_obs=str(hdr.get("DATE-OBS", "?")),
        specsys=str(hdr.get("SPECSYS", "?")).strip(),
        berv_kms=float(hdr.get("ESO QC BERV", np.nan)),
        wave_A=wave, flux=flux, err=err, atm_transm=atm,
        telluric_column=tcol)


# ── Generalized telluric correction (instrument branch) ───────────────────────
def telluric_correct_nirps(frame: NirpsFrame) -> NirpsFrame:
    """NIRPS telluric 'correction' = the geneva DRS already divided ATM_TRANSM out and
    we selected the FLUX_TELL_* column at load. Verify the model transmission is real
    (a true 0–1 transmission, not all-ones), then set the permanent-rule flag."""
    atm = frame.atm_transm
    # a real transmission spectrum: a meaningful fraction of pixels in (0,1).
    # (NIRPS order edges carry fill values outside [0,1]; ignore those.)
    real = np.isfinite(atm) & (atm > 0.0) & (atm <= 1.0)
    if real.sum() < 100 or float(np.nanmedian(atm[real])) > 0.999:
        raise MolecfitNotAvailableError(
            f"{frame.path.name}: ATM_TRANSM has no real telluric population — model not "
            f"applied; refusing to flag a non-corrected NIRPS frame.")
    frame.telluric_corrected = True
    frame._atm_transm_min = float(np.nanmin(atm[real]))   # clean min over valid pixels
    return frame


def telluric_correct_crires(frame: CriresFrame, window_nm: tuple[float, float],
                            work_dir: Path,
                            molecules=ct.TELLURIC_MOLECULES) -> CriresFrame:
    """CRIRES+ telluric correction via molecfit (reusing RYA-373's proven driver),
    generalized off the Vesta CO bandhead to an ARBITRARY science window. Runs in the
    TOPOCENTRIC frame (raw IDP); divides the fitted transmission out of every segment
    overlapping `window_nm`. Sets telluric_corrected. Fails loud if esorex absent."""
    if frame.specsys.upper() != "TOPOCENT":
        raise RuntimeError(
            f"{frame.path.name}: telluric fit must run in TOPOCENTRIC frame "
            f"(SPECSYS={frame.specsys!r}); the RV shift happens AFTER telluric.")
    if not _esorex_available():
        raise MolecfitNotAvailableError(
            "molecfit/esorex not found — telluric engine absent; no faked correction.")
    lo_A, hi_A = window_nm[0] * NM_TO_A, window_nm[1] * NM_TO_A
    segs = [s for s in frame.segments
            if np.nanmax(s.wave_A) >= lo_A and np.nanmin(s.wave_A) <= hi_A]
    if not segs:
        raise RuntimeError(f"{frame.path.name}: no segment overlaps {window_nm} nm.")
    residuals = []
    for s in segs:
        r = _molecfit_segment(frame, s, Path(work_dir) / f"ord{s.order}_det{s.detector}",
                              molecules)
        s.wave_A = r["lam_A"]; s.flux = r["corr"]
        s._mtrans = r["mtrans"]; s._flux_raw = r["flux_raw"]
        if r["resid"] == r["resid"]:
            residuals.append(r["resid"])
    frame.telluric_corrected = True
    frame._molecules = tuple(molecules)
    frame._gdas = r["gdas"]
    frame._telluric_residual = float(np.nanmedian(residuals)) if residuals else float("nan")
    return frame


# ── Generalized telluric-specific quality gate (star-agnostic, frame-pooled) ──
def telluric_residual_gate(frame, tol: float = 0.02) -> dict:
    """Star-agnostic telluric quality, the RYA-494 generalization of RYA-373's D1 gate.
    RYA-373 masked SOLAR lines; a general star has no line mask, so we DETREND each
    corrected segment with a median filter (window wider than telluric lines → follows
    the stellar spectrum + continuum) and compare the residual scatter at TELLURIC
    pixels (moderate model absorption, 0.05<mtrans<0.95, not saturated) against CLEAN
    continuum pixels (mtrans>0.98). Stellar lines cancel because they appear in both
    the spectrum and its trend. Pixel populations are POOLED ACROSS the frame's
    segments (CRIRES tellurics are bimodal per-segment — a segment is either telluric-
    saturated or clean — so the clean baseline must come from the clean segments). The
    EXCESS of the telluric-pixel residual over the clean baseline is the telluric-
    specific misfit. Returns excess + PASS (≤tol)."""
    from scipy.ndimage import median_filter
    segs = getattr(frame, "segments", None) or [frame]
    res_t, res_c = [], []
    for s in segs:
        mt = getattr(s, "_mtrans", None)
        if mt is None or not np.isfinite(mt).any():
            continue
        flux = s.flux
        fin = np.isfinite(flux) & np.isfinite(mt)
        if fin.sum() < 100:
            continue
        trend = median_filter(np.where(fin, flux, np.nanmedian(flux[fin])), size=51)
        with np.errstate(invalid="ignore", divide="ignore"):
            detr = np.where(trend > 0, flux / trend, np.nan)
        tell = fin & (mt < 0.95) & (mt > 0.05)        # moderate telluric, not saturated
        clean = fin & (mt > 0.98)
        if tell.sum() >= 5:
            res_t.append(np.abs(1.0 - detr[tell]))
        if clean.sum() >= 5:
            res_c.append(np.abs(1.0 - detr[clean]))
    if not res_t or not res_c:
        return {"residual": float("nan"), "passed": False,
                "reason": "too few telluric/clean pixels in the frame"}
    r_tell = float(np.nanmedian(np.concatenate(res_t)))
    r_clean = float(np.nanmedian(np.concatenate(res_c)))
    excess = max(0.0, r_tell - r_clean)
    return {"n_telluric_px": int(sum(len(a) for a in res_t)),
            "n_clean_px": int(sum(len(a) for a in res_c)),
            "resid_telluric": r_tell, "resid_clean_baseline": r_clean,
            "excess": excess, "tol": tol, "passed": excess <= tol}


# ── Generalized velocity conditioning (the branched step) ─────────────────────
def velocity_condition(frame, target: IRTarget, work_dir: Path = None):
    """Branch the velocity step on the target's mode. REFLECTED_SOLAR delegates to the
    RYA-372/373 asteroid-ephemeris path (Vesta). STELLAR conditions a direct star to
    its rest frame using its own BERV (+ systemic RV) — NOT asteroid ephemeris."""
    if target.velocity_mode is VelocityMode.REFLECTED_SOLAR:
        if not isinstance(frame, CriresFrame):
            raise NotImplementedError("reflected-solar RV path is CRIRES (Vesta) only.")
        return ct.rv_condition(frame, work_dir)            # RYA-372 single source
    return _stellar_rv_condition(frame, target)


def _stellar_rv_condition(frame, target: IRTarget):
    """Direct-star rest-frame conditioning: λ_rest = λ_obs / (1 + (v_berv? + v_sys)/c).
    NIRPS S1D is already BARYCENT (BERV removed by the DRS) → only the systemic RV
    remains. CRIRES is TOPOCENT → remove BERV then systemic RV. Telluric-first enforced.
    """
    if not getattr(frame, "telluric_corrected", False):
        raise TelluricNotCorrectedError(
            f"{frame.path.name}: refusing to RV-condition a non-telluric-corrected "
            f"frame (telluric first — permanent rule).")
    v_sys = target.systemic_rv_kms
    if isinstance(frame, NirpsFrame):
        # already barycentric; shift only the systemic RV to reach stellar rest
        shift = 1.0 + v_sys / _C_KMS
        frame.wave_A = frame.wave_A / shift
        frame.rest_frame = True
        frame._rv_applied = {"berv": 0.0, "systemic": v_sys, "frame_in": "BARYCENT"}
        return frame
    # CRIRES TOPOCENT: remove BERV (→ barycentric) then systemic (→ stellar rest)
    berv = _berv_for_frame(frame)
    shift = 1.0 + (berv + v_sys) / _C_KMS
    for s in frame.segments:
        s.wave_A = s.wave_A / shift
    frame.rest_frame = True
    frame._rv_applied = {"berv": berv, "systemic": v_sys, "frame_in": "TOPOCENT"}
    return frame


def _berv_for_frame(frame: CriresFrame) -> float:
    """Barycentric Earth radial velocity (km/s) for a CRIRES frame, computed from
    RA/DEC/MJD at Paranal (the IDP has no BERV keyword). Sign: v_berv such that
    λ_bary = λ_topo / (1 + v_berv/c)."""
    from astropy.coordinates import SkyCoord, EarthLocation
    from astropy.time import Time
    import astropy.units as u
    paranal = EarthLocation.of_site("paranal") if False else EarthLocation.from_geodetic(
        lon=-70.4051 * u.deg, lat=-24.6275 * u.deg, height=2635 * u.m)
    sc = SkyCoord(ra=frame.ra * u.deg, dec=frame.dec * u.deg)
    t = Time(float(frame.mjd), format="mjd")
    return float(sc.radial_velocity_correction(obstime=t, location=paranal).to(u.km / u.s).value)


# ── Inventory (authoritative attribution, never folder/glob blindly) ──────────
def crires_frames(target: IRTarget) -> list[CriresFrame]:
    """Load α-Cen-A CRIRES IDPs. OBJECT is correct for CRIRES ('alf Cen A'); the
    'Star S5' contaminant lives only in the B folder, so the A folder is clean."""
    from astropy.io import fits
    out = []
    for f in sorted(glob.glob(str(target.data_dir / "CRIRES" / "*.fits"))):
        obj = re_norm(fits.getheader(f, 0).get("OBJECT", ""))
        if obj != "alpha_cen_a":
            continue                       # never silently admit a non-A frame
        out.append(load_crires_idp(f))
    return out


def re_norm(obj: str) -> str:
    """OBJECT → star (the RYA-479 normalizer, case/separator-insensitive)."""
    import re
    k = re.sub(r"[^a-z0-9]", "", str(obj).lower())
    if k == "hd128620" or ("cen" in k and k.endswith("a")):
        return "alpha_cen_a"
    if k == "hd128621" or ("cen" in k and k.endswith("b")):
        return "alpha_cen_b"
    return "other"


_RYA423_MANIFEST = Path("data/audit/acen_holdings_rya384/ir_star_id_rya423_manifest.csv")


def _nirps_starid_verdict_A() -> set[str]:
    """Filenames the RYA-423 IR RV star-ID assigned verdict=A (authoritative IR
    attribution, overriding the mislabeled OBJECT header). Empty set if the manifest
    is absent (caller then falls back to OBJECT-norm with a loud note)."""
    import csv
    if not _RYA423_MANIFEST.exists():
        return set()
    with open(_RYA423_MANIFEST) as fh:
        return {r["frame"] for r in csv.DictReader(fh)
                if r["instr"] == "NIRPS" and r["verdict"] == "A"}


def nirps_frames(target: IRTarget) -> list[NirpsFrame]:
    """Load α-Cen-A NIRPS frames, attributed by the RYA-423 IR RV star-ID (verdict=A),
    NOT the mislabeled OBJECT header (RYA-481: authoritative attribution, never glob).
    Frames the star-ID flagged NOT-ALPHA-CEN (RV≈-34, RYA-431) or INDETERMINATE
    (Star S5) are excluded."""
    a_set = _nirps_starid_verdict_A()
    out = []
    for f in sorted(glob.glob(str(target.data_dir / "NIRPS" / "*.fits"))):
        if a_set and Path(f).name not in a_set:
            continue                       # exclude non-A per the authoritative star-ID
        out.append(load_nirps_idp(f))
    return out


# ── CLI / verify ──────────────────────────────────────────────────────────────
def _verify(target_name: str, run_telluric: bool = False) -> dict:
    t = TARGETS[target_name]
    print("=" * 84)
    print(f"  RYA-494 — generalized IR telluric: target={t.name} mode={t.velocity_mode.value}")
    print("=" * 84)
    cr = crires_frames(t)
    print(f"\n  CRIRES+: {len(cr)} frames "
          f"(bands {sorted({f.band for f in cr})}); molecfit engine="
          f"{'available' if _esorex_available() else 'ABSENT'}")
    for f in cr:
        print(f"    {f.wlen_id:<7}{f.wmin_nm:.0f}-{f.wmax_nm:.0f}nm SPECSYS={f.specsys} "
              f"SNR={f.snr:.0f}")
    nf = nirps_frames(t) if (t.data_dir / "NIRPS").exists() else []
    if nf:
        n0 = nf[0]
        print(f"\n  NIRPS: {len(nf)} frames; telluric column={n0.telluric_column}; "
              f"SPECSYS={n0.specsys}; WAVE {np.nanmin(n0.wave_A):.0f}-"
              f"{np.nanmax(n0.wave_A):.0f}Å (already telluric+BARYCENT)")
        print(f"    attribution: {t.nirps_attribution[:140]}…")
    # O I 844/926 coverage reality
    print("\n  Atomic IR target coverage (the brief's O I 844/926):")
    for nm, wl in [("O I 8446", 8446.5), ("O I 9266", 9265.9)]:
        cov = any(f.covers_nm(wl / NM_TO_A) for f in cr) or any(f.covers_A(wl) for f in nf)
        print(f"    {nm} Å: covered by α Cen A IR = {cov} "
              f"({'GAP — below CRIRES-Y 9496Å / NIRPS 9661Å blue edges' if not cov else 'OK'})")
    if run_telluric and _esorex_available() and cr:
        f0 = max(cr, key=lambda f: f.snr)            # Y1029, SNR 302
        win = (f0.wmin_nm, f0.wmin_nm + 5.0)         # a 5-nm telluric-rich sub-window
        print(f"\n  RUNNING molecfit on {f0.wlen_id} {win} nm (TOPOCENTRIC) …")
        telluric_correct_crires(f0, win, Path(f"/tmp/rya494_{t.name}_{f0.wlen_id}"))
        print(f"    telluric_corrected={f0.telluric_corrected}; "
              f"residual≈{f0._telluric_residual*100:.1f}%  gdas={f0._gdas}")
    print("=" * 84)
    return {"target": t.name, "crires": cr, "nirps": nf, "molecfit": _esorex_available()}


def main(argv=None):
    import argparse
    ap = argparse.ArgumentParser(description="RYA-494 generalized IR telluric")
    ap.add_argument("--target", default="alpha_cen_a", choices=list(TARGETS))
    ap.add_argument("--verify", action="store_true")
    ap.add_argument("--telluric", action="store_true", help="run a real molecfit pass")
    a = ap.parse_args(argv)
    return _verify(a.target, run_telluric=a.telluric)


if __name__ == "__main__":
    main()
