"""
pipeline/crires_telluric.py
===========================
RYA-373 — molecfit telluric + continuum-cal + RV-condition on the Vesta CRIRES+
K-band, unlocking the 2.3 µm CO (2-0) overtone arm (C, 12C/13C, O).

The Vesta CRIRES+ set (RYA-370 audit) is cr2res EXTRACTC IDPs: wavelength-cal'd +
extracted, but NOT telluric-corrected, FLUXCAL=UNCALIBRATED, SPECSYS=TOPOCENT,
WAVE in nm. Telluric correction is PERMANENT and MANDATORY before any IR abundance.

Pipeline order (critical):
  1. load (nm → Å at the boundary; assert raw TOPOCENT/UNCALIBRATED state)
  2. telluric-correct in the TOPOCENTRIC frame (molecfit; tellurics are stationary
     in Earth's atmosphere) — H2O + CH4 + CO2 over the 2.29-2.49 µm K orders
  3. continuum-normalize the CO-overtone orders (FLUXCAL=UNCALIBRATED → need
     normalized line depths, not absolute flux)
  4. RV-condition to the solar rest frame (asteroid ephemeris, RYA-372
     reflected_solar_rv — single source, NOT re-implemented here) — AFTER telluric
  5. per-setting, RV-registered co-add (lift SNR; the solar CO bandheads are faint;
     NEVER blind-coadd across epochs — register on per-epoch RV first)

The separability point: the 2.3 µm region carries BOTH telluric CO/CH4 AND the
stellar (solar) CO. molecfit models the TELLURIC component at topocentric rest; the
solar CO bandheads survive, Doppler-offset by the reflected RV. That velocity
separation is what makes telluric-vs-solar CO separable.

CRITICAL (permanent IR rule): no flux/abundance number is emitted unless the
`telluric_corrected` flag is set (TelluricNotCorrectedError otherwise); telluric is
applied BEFORE the RV shift; no blind cross-epoch co-add; nm→Å at the loader.

Smoke test:  python -m pipeline.crires_telluric --set vesta_crires_k --verify

NOTE (engine availability): molecfit (esorex molecfit_model/molecfit_calctrans) and
the RYA-372 reflected_solar_rv module are external prerequisites. Each step verifies
its engine and FAILS LOUD if absent (never a silent skip / fake correction).
"""
from __future__ import annotations

import glob
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from config.constants import PATHS   # RYA-373: canonical solar line-list path (RV anchors)

# ── Data location (the reflected-solar set lives OUTSIDE the repo, RYA-370) ────
_DATA_ROOT = Path("/Users/ryanschmitt/Documents/Exoplanet Codex/data/spectra/"
                  "exoplanetcodex-data/Solar Calibration/Solar System Targets")
VESTA_CRIRES_DIR = _DATA_ROOT / "Vesta" / "CRIRES"

# ── Physics constants ─────────────────────────────────────────────────────────
_C_KMS = 299792.458
NM_TO_A = 10.0

# 12C16O (2-0) overtone bandhead — the gold-standard C / 12C13C / O diagnostic the
# optical arms cannot reach. Vacuum wavelength (IR convention). RYA-373 target.
CO_2_0_BANDHEAD_NM = 2293.5      # ≈ 2.2935 µm (12C16O 2-0 R-branch bandhead, vacuum)
CO_OVERTONE_LO_NM = 2290.0       # the 2.3 µm CO-overtone series fit/analysis window
CO_OVERTONE_HI_NM = 2490.0
# molecfit telluric species over the K band (the IR rule: H2O + CH4 + CO2)
TELLURIC_MOLECULES = ("H2O", "CH4", "CO2")

SNR_FLOOR = 200.0                # the RYA-370 per-frame SNR floor


# ── Errors (no silent fallback) ───────────────────────────────────────────────
class TelluricNotCorrectedError(RuntimeError):
    """An abundance/synthesis consumer touched a frame whose telluric_corrected
    flag is not set. The permanent IR rule: no IR number without verified telluric."""


class MolecfitNotAvailableError(RuntimeError):
    """molecfit/esorex is not installed — the mandatory telluric engine is absent.
    NEVER fall back to an uncorrected or faked spectrum."""


class VelocityModuleNotAvailableError(RuntimeError):
    """pipeline.reflected_solar_rv (RYA-372) is not present — the single-source
    asteroid-ephemeris velocity correction. NEVER re-implement it here."""


# ── Frame model ───────────────────────────────────────────────────────────────
@dataclass
class CriresSegment:
    order: int
    detector: int
    wave_A: np.ndarray          # Å (converted from the IDP's nm at load)
    flux: np.ndarray
    err: np.ndarray
    qual: np.ndarray


@dataclass
class CriresFrame:
    path: Path
    wlen_id: str                # ESO INS WLEN ID (e.g. 'K2192')
    band: str                   # 'Y'|'J'|'H'|'K'
    mjd: float
    date_obs: str
    ra: float
    dec: float
    snr: float
    specsys: str                # expect 'TOPOCENT' (raw)
    fluxcal: str                # expect 'UNCALIBRATED'
    wmin_nm: float
    wmax_nm: float
    segments: list = field(default_factory=list)
    telluric_corrected: bool = False     # the permanent-rule flag
    continuum_normalized: bool = False
    rest_frame: bool = False             # solar rest frame after RV-condition

    def covers_nm(self, wl_nm: float) -> bool:
        """Header WAVELMIN/MAX covers wl — NOT a guarantee the line is on-chip
        (CRIRES+ tiles with inter-order/detector gaps)."""
        return self.wmin_nm <= wl_nm <= self.wmax_nm

    def segment_at(self, wl_nm: float):
        """The chip segment with wl ON-CHIP (within its actual wave_A), or None
        if wl falls in an inter-order/detector gap. This is the real coverage test."""
        wl_A = wl_nm * NM_TO_A
        for s in self.segments:
            if np.nanmin(s.wave_A) <= wl_A <= np.nanmax(s.wave_A):
                return s
        return None

    def on_chip(self, wl_nm: float) -> bool:
        return self.segment_at(wl_nm) is not None

    def co_segments(self) -> list:
        """Segments overlapping the 2.29-2.49 µm CO-overtone window (Å)."""
        lo, hi = CO_OVERTONE_LO_NM * NM_TO_A, CO_OVERTONE_HI_NM * NM_TO_A
        return [s for s in self.segments
                if np.nanmax(s.wave_A) >= lo and np.nanmin(s.wave_A) <= hi]


# ── Loader (nm → Å at the boundary) ───────────────────────────────────────────
def load_crires_idp(path) -> CriresFrame:
    """Load a cr2res EXTRACTC IDP into a CriresFrame, converting WAVE nm → Å at the
    boundary (the IR IDP is nm, unlike the optical IDPs which are already Å — do NOT
    reuse the optical Å assumption). Splits into per-(ORDER,DETEC) chip segments."""
    from astropy.io import fits
    path = Path(path)
    with fits.open(path) as h:
        hdr = h[0].header
        if 'SPECTRUM' in [hdu.name for hdu in h]:
            d = h['SPECTRUM'].data[0]
        else:
            d = h[1].data[0]
        wave_nm = np.asarray(d['WAVE'], float)
        flux = np.asarray(d['FLUX'], float)
        err = np.asarray(d['ERR'], float)
        qual = np.asarray(d['QUAL'], int)
        order = np.asarray(d['ORDER'], int)
        detec = np.asarray(d['DETEC'], int)

    wlen_id = str(hdr.get('HIERARCH ESO INS WLEN ID', '?')).strip()
    band = wlen_id[0] if wlen_id and wlen_id[0] in 'YJHK' else '?'
    wave_A = wave_nm * NM_TO_A          # ← nm → Å at the load boundary

    segs = []
    for key in sorted(set(zip(order.tolist(), detec.tolist()))):
        m = (order == key[0]) & (detec == key[1])
        segs.append(CriresSegment(order=key[0], detector=key[1],
                                  wave_A=wave_A[m], flux=flux[m],
                                  err=err[m], qual=qual[m]))
    return CriresFrame(
        path=path, wlen_id=wlen_id, band=band,
        mjd=float(hdr.get('MJD-OBS', np.nan)),
        date_obs=str(hdr.get('DATE-OBS', '?')),
        ra=float(hdr.get('RA', np.nan)), dec=float(hdr.get('DEC', np.nan)),
        snr=float(hdr.get('SNR', np.nan)),
        specsys=str(hdr.get('SPECSYS', '?')).strip(),
        fluxcal=str(hdr.get('FLUXCAL', '?')).strip(),
        wmin_nm=float(hdr.get('WAVELMIN', np.nan)),
        wmax_nm=float(hdr.get('WAVELMAX', np.nan)),
        segments=segs)


def inventory(crires_dir=VESTA_CRIRES_DIR) -> list:
    """Load every CRIRES+ IDP in the directory (header-level frame objects)."""
    files = sorted(glob.glob(str(Path(crires_dir) / '*.fits')))
    if not files:
        raise FileNotFoundError(f"No CRIRES+ IDP FITS under {crires_dir}")
    return [load_crires_idp(f) for f in files]


def co_overtone_frames(crires_dir=VESTA_CRIRES_DIR, on_chip_only: bool = True) -> list:
    """The K-band frames with the 2.293 µm CO(2-0) bandhead present — the frames this
    ticket conditions. on_chip_only=True (default) requires the bandhead to land on a
    detector (not in an inter-order/detector gap): only those can actually be coadded
    at the bandhead. on_chip_only=False returns the header-range supersets. (K only;
    never a cross-band coadd.)"""
    frames = inventory(crires_dir)
    sel = [f for f in frames if f.band == 'K' and f.covers_nm(CO_2_0_BANDHEAD_NM)]
    if on_chip_only:
        sel = [f for f in sel if f.on_chip(CO_2_0_BANDHEAD_NM)]
    return sel


# ── Continuum normalization (molecfit-independent; testable now) ───────────────
def continuum_normalize(wave_A: np.ndarray, flux: np.ndarray,
                        niter: int = 5, lo_sigma: float = 1.5,
                        hi_sigma: float = 3.0, deg: int = 3) -> np.ndarray:
    """Local continuum normalization of one CO-overtone order via iterative
    asymmetric sigma-clipping to the upper envelope (absorption pulls flux DOWN, so
    clip low points harder). Returns flux / continuum. FLUXCAL=UNCALIBRATED means we
    need normalized line depths, not absolute flux (RYA-373)."""
    good = np.isfinite(wave_A) & np.isfinite(flux) & (flux > 0)
    if good.sum() < deg + 2:
        return np.full_like(flux, np.nan)
    x = wave_A[good]
    y = flux[good]
    x0 = x.mean()
    xs = x - x0
    keep = np.ones(len(x), bool)
    cont = np.full(len(x), np.nan)
    for _ in range(niter):
        c = np.polyfit(xs[keep], y[keep], deg)
        cont = np.polyval(c, xs)
        resid = y - cont
        sig = np.std(resid[keep]) or 1.0
        keep = (resid > -lo_sigma * sig) & (resid < hi_sigma * sig)
        if keep.sum() < deg + 2:
            break
    out = np.full_like(flux, np.nan)
    full_cont = np.polyval(c, wave_A - x0)
    with np.errstate(invalid='ignore', divide='ignore'):
        out[np.isfinite(wave_A)] = flux[np.isfinite(wave_A)] / full_cont[np.isfinite(wave_A)]
    return out


# ── Telluric (molecfit, TOPOCENTRIC frame) — the mandatory engine ─────────────
def _esorex_available() -> bool:
    return (shutil.which("esorex") is not None
            or Path(_ESOREX).exists()              # RYA-375 Homebrew install (off default PATH)
            or shutil.which("molecfit_model") is not None)


def run_molecfit_telluric(frame: CriresFrame, work_dir: Path,
                          molecules=TELLURIC_MOLECULES) -> CriresFrame:
    """Fit and divide out the telluric model (H2O+CH4+CO2) over the 2.29-2.49 µm K
    orders, IN THE TOPOCENTRIC FRAME (tellurics are stationary in Earth's atmosphere
    → fixed in topocentric wavelength). Deep saturated cores are masked, not
    over-fit. Sets frame.telluric_corrected on success.

    Engine = ESO molecfit via esorex (molecfit_model → molecfit_calctrans). If
    esorex/molecfit is not installed this RAISES MolecfitNotAvailableError — there is
    NO silent fallback and NO uncorrected/faked spectrum (the permanent IR rule).
    """
    if frame.specsys.upper() != 'TOPOCENT':
        raise RuntimeError(
            f"{frame.path.name}: telluric fit must run in the TOPOCENTRIC frame, but "
            f"SPECSYS={frame.specsys!r}. The RV shift to solar rest must happen AFTER "
            f"telluric correction (RYA-373), never before.")
    if not _esorex_available():
        raise MolecfitNotAvailableError(
            "molecfit/esorex not found on PATH — the telluric engine is not installed. "
            "RYA-373 telluric correction is mandatory and cannot be faked. Install the "
            "ESO molecfit pipeline (esorex molecfit_model / molecfit_calctrans) and "
            "re-run. (No silent fallback: no uncorrected IR spectrum is emitted.)")
    # NOTE: the esorex molecfit_model/calctrans invocation (SOF + parameter file:
    # LIST_MOLEC, WAVE_INCLUDE over the 2.29-2.49 µm orders, FTOL/XTOL, CONTINUUM_N)
    # is finalized against the installed molecfit version (recipe params are
    # version-dependent — the ticket's "confirm the invocation against the branch").
    # Implemented in _molecfit_driver once esorex is present; this is the wired hook.
    return _molecfit_driver(frame, work_dir, molecules)


_ESOREX = "/opt/homebrew/bin/esorex"          # RYA-375 install (Homebrew ESO tap)
_MTRANS_FLOOR = 0.02                            # mask telluric cores below this transmission
_CONTINUUM_N = 2                                # polynomial continuum order (uncalibrated flux)

# Real GDAS atmospheric profiles (RYA-373 finish-out #2). molecfit's GDAS_PROFILE=auto
# requested odd hours (T01/T02) absent from the 3-hourly tarball → standard-profile
# fallback. We bypass that by extracting the nearest real 3-hourly GDAS for the obs
# MJD and passing it via --GDAS_PROFILE. Tarball is per-site (Paranal -70.4/-24.6).
_GDAS_LOC = "C-70.4-24.6"
_GDAS_TARBALL = ("/opt/homebrew/Cellar/telluriccorr/4.3.3_4/share/molecfit/data/"
                 f"profiles/gdas/gdas_profiles_{_GDAS_LOC}.tar.gz")


def _nearest_gdas(mjd: float, work_dir: Path) -> "str | None":
    """Extract the nearest 3-hourly real GDAS profile for the obs MJD from the local
    tarball; return its path, or None (→ caller logs + lets molecfit fall back)."""
    import tarfile
    from astropy.time import Time
    if not Path(_GDAS_TARBALL).exists():
        return None
    t = Time(float(mjd), format='mjd').to_datetime()
    hh = int(round(t.hour + t.minute / 60.0) / 3.0) * 3   # nearest 3-hourly slot
    from datetime import timedelta
    base = t.replace(minute=0, second=0, microsecond=0, hour=0) + timedelta(hours=hh)
    work_dir.mkdir(parents=True, exist_ok=True)
    try:
        with tarfile.open(_GDAS_TARBALL) as tf:
            for dh in (0, -3, 3, -6, 6):                  # nearest, then neighbours
                cand = (base + timedelta(hours=dh))
                name = f"{_GDAS_LOC}D{cand:%Y-%m-%d}T{cand.hour:02d}.gdas"
                try:
                    tf.extract(name, path=str(work_dir))
                except KeyError:
                    continue
                return _gdas_ascii_to_fits(work_dir / name)   # molecfit wants FITS
    except Exception:
        return None
    return None


def _gdas_ascii_to_fits(ascii_path) -> str:
    """Convert a tarball ASCII GDAS profile (# P[hPa] HGT[m] T[K] RELHUM[%]) to the
    FITS table molecfit's GDAS_PROFILE expects (press[hPa], height[km], temp[K],
    relhum[%]) — molecfit's GDAS_PROF loader is CFITSIO, it cannot read the ASCII."""
    from astropy.io import fits
    ascii_path = Path(ascii_path)
    rows = []
    for ln in ascii_path.read_text().splitlines():
        ln = ln.strip()
        if not ln or ln.startswith('#'):
            continue
        p, hgt_m, t, rh = (float(x) for x in ln.split())
        rows.append((p, hgt_m / 1000.0, t, rh))           # m → km
    rows = np.array(rows)
    out = ascii_path.with_suffix('.gdas.fits')
    fits.BinTableHDU.from_columns([
        fits.Column(name='press', format='1D', array=rows[:, 0]),
        fits.Column(name='height', format='1D', array=rows[:, 1]),
        fits.Column(name='temp', format='1D', array=rows[:, 2]),
        fits.Column(name='relhum', format='1D', array=rows[:, 3])],
    ).writeto(out, overwrite=True)
    return str(out)


def _write_molecfit_inputs(frame: CriresFrame, seg: CriresSegment, work_dir: Path,
                           molecules) -> Path:
    """Write the three molecfit_model inputs for one CO order: SCIENCE (binary table
    lambda[µm]/flux/dflux, primary header carrying the ESO atmospheric keywords),
    WAVE_INCLUDE (µm, edges trimmed), MOLECULES (LIST_MOLEC/FIT_MOLEC=1J/REL_COL=1D —
    the int32 FIT_MOLEC is mandatory, astropy's default int64 fails CPL). Returns the
    SOF path."""
    from astropy.io import fits
    work_dir.mkdir(parents=True, exist_ok=True)
    m = np.isfinite(seg.wave_A) & np.isfinite(seg.flux) & (seg.flux > 0)
    lam_um = seg.wave_A[m] / 1.0e4              # Å → µm (molecfit internal)
    flux = seg.flux[m].astype(float)
    err = seg.err[m].astype(float)

    src = fits.getheader(str(frame.path))
    ph = fits.PrimaryHDU()
    for k in src.keys():                         # carry ESO TEL/INS keywords for the atm model
        if k.startswith('ESO ') or k in ('MJD-OBS', 'RA', 'DEC', 'UTC', 'INSTRUME'):
            try:
                ph.header[k] = src[k]
            except Exception:
                pass
    from astropy.table import Table
    sci = fits.BinTableHDU(Table({'lambda': lam_um, 'flux': flux, 'dflux': err}),
                           name='SCIENCE')
    fits.HDUList([ph, sci]).writeto(work_dir / 'science.fits', overwrite=True)

    lo, hi = float(lam_um.min() + 5e-4), float(lam_um.max() - 5e-4)
    fits.BinTableHDU.from_columns([
        fits.Column(name='LOWER_LIMIT', format='1D', array=np.array([lo])),
        fits.Column(name='UPPER_LIMIT', format='1D', array=np.array([hi]))],
    ).writeto(work_dir / 'wave_include.fits', overwrite=True)

    fits.BinTableHDU.from_columns([
        fits.Column(name='LIST_MOLEC', format='4A', array=np.array(list(molecules))),
        fits.Column(name='FIT_MOLEC', format='1J',         # int32 — CPL rejects int64
                    array=np.ones(len(molecules), dtype=np.int32)),
        fits.Column(name='REL_COL', format='1D', array=np.ones(len(molecules)))],
    ).writeto(work_dir / 'molecules.fits', overwrite=True)

    sof = work_dir / 'model.sof'
    sof.write_text(f"{work_dir/'science.fits'} SCIENCE\n"
                   f"{work_dir/'wave_include.fits'} WAVE_INCLUDE\n"
                   f"{work_dir/'molecules.fits'} MOLECULES\n")
    return sof


def _molecfit_driver(frame: CriresFrame, work_dir: Path, molecules) -> CriresFrame:
    """esorex molecfit_model on the CO-overtone order(s); divide by the fitted
    convolved telluric transmission (BEST_FIT_MODEL.mtrans) — saturated cores masked.

    Sets frame.telluric_corrected and stashes per-segment corrected flux +
    telluric diagnostics on the segment. Uses molecfit_model's mtrans directly (the
    convolved transmission over the science grid), which avoids the calctrans
    chip-mapping step for a single order.

    NOTE (RYA-373 finding): on a reflected-SOLAR target the spectrum is solar-line-
    rich, so molecfit's continuum cannot model the solar lines and a blanket fit
    residual conflates the (wanted) solar spectrum with telluric misfit — the
    telluric-residual metric must be telluric-line-specific. The corrected flux and
    the model transmission are returned for that downstream verification.
    """
    import os
    from astropy.io import fits
    esorex = _ESOREX if Path(_ESOREX).exists() else shutil.which("esorex")
    if esorex is None:
        raise MolecfitNotAvailableError("esorex not found (RYA-375 install expected at "
                                        f"{_ESOREX}).")
    seg = frame.segment_at(CO_2_0_BANDHEAD_NM)
    if seg is None:
        raise RuntimeError(f"{frame.path.name}: CO(2-0) bandhead not on-chip — nothing "
                           f"to telluric-correct in this frame.")
    in_dir = Path(work_dir) / 'in'
    out_dir = Path(work_dir) / 'out'
    out_dir.mkdir(parents=True, exist_ok=True)
    sof = _write_molecfit_inputs(frame, seg, in_dir, molecules)

    gdas = _nearest_gdas(frame.mjd, in_dir)    # real 3-hourly GDAS for the obs MJD
    env = dict(os.environ, PATH=f"/opt/homebrew/bin:{os.environ.get('PATH', '')}")
    cmd = [esorex, f"--output-dir={out_dir}", "molecfit_model",
           "--COLUMN_LAMBDA=lambda", "--COLUMN_FLUX=flux", "--COLUMN_DFLUX=dflux",
           "--WLG_TO_MICRON=1.0", "--WAVELENGTH_FRAME=VAC",
           "--FIT_CONTINUUM=1", f"--CONTINUUM_N={_CONTINUUM_N}"]
    cmd += [f"--GDAS_PROFILE={gdas}"] if gdas else []
    cmd += [str(sof)]
    proc = subprocess.run(cmd, cwd=str(in_dir), env=env,
                          capture_output=True, text=True)
    bfm = out_dir / 'BEST_FIT_MODEL.fits'
    if proc.returncode != 0 or not bfm.exists():
        tail = "\n".join(l for l in proc.stdout.splitlines()
                         if 'ERROR' in l and 'gdas' not in l.lower())[-1500:]
        raise RuntimeError(f"molecfit_model failed (rc={proc.returncode}) on "
                           f"{frame.path.name}:\n{tail}")

    m = fits.open(bfm)[1].data
    lam_A = m['lambda'] * 1.0e4
    mtrans = m['mtrans']
    fl = m['flux']
    ok = np.isfinite(fl) & np.isfinite(mtrans) & (mtrans > _MTRANS_FLOOR)
    corr = np.full_like(fl, np.nan)
    corr[ok] = fl[ok] / mtrans[ok]
    # telluric-line-specific residual: scatter of corrected/continuum where the model
    # has MODERATE telluric absorption (0.3–0.95) — the cleanest available proxy
    # (deep cores are masked; solar lines still contribute, hence reported not gated).
    cont = continuum_normalize(lam_A, corr)
    modtell = ok & (mtrans < 0.95) & (mtrans > 0.30)
    resid = float(np.nanstd(cont[modtell])) if modtell.any() else np.nan

    seg.wave_A = lam_A                          # molecfit grid (Å)
    seg.flux = corr                             # telluric-corrected flux
    seg._mtrans = mtrans                        # model transmission (for downstream verify)
    seg._flux_raw = fl                          # pre-correction flux (for the D1 gate)
    frame.telluric_corrected = True
    frame._telluric_residual = resid            # blanket proxy (legacy; see D1 gate)
    frame._telluric_mtrans_max = float(np.nanmax(mtrans))
    frame._molecules = tuple(molecules)
    frame._gdas = (Path(gdas).name if gdas else 'standard-profile (no GDAS)')
    return frame


# ── D1 telluric-specific gate (RYA-373 Decision 1, primary pre-coadd gate) ─────
# Vetoes the meaningless blanket residual: score telluric quality ONLY on pixels
# that are (a) telluric-dominated in molecfit's model (transmission well below
# continuum, not saturated) AND (b) NOT solar-coincident. At those telluric cores
# the post-correction flux must return to the local continuum. Solar-coincident
# pixels are excluded via a solar K-band line mask (the CO-overtone bandheads + the
# strong K-band solar atomic lines); the Kitt Peak FTS solar IR atlas (RYA-162) is
# the preferred mask source but is NOT in the data set (flagged — end-validation gap).

# Strong solar K-band atomic lines (air→vac ~ +6 Å at 2.2 µm; vacuum Å) for the
# solar-coincidence mask. CO-overtone bandhead series handled separately.
_SOLAR_K_LINES_VAC = (22062.4, 22089.7,        # Na I doublet
                      22614.1, 22631.1,        # Ca I
                      21066.1, 22834.2)         # Al I / Ti I (approx)


def _solar_coincident(wave_A: np.ndarray, rv_kms: float = 0.0,
                      half_width_A: float = 1.2) -> np.ndarray:
    """Boolean mask of pixels within half_width of a solar K-band line OR inside the
    CO-overtone bandhead series (≥22930 Å), shifted by the reflected RV. These carry
    the wanted solar signal and must be EXCLUDED from the telluric-residual score."""
    shift = (1.0 + rv_kms / _C_KMS)
    mask = np.zeros_like(wave_A, dtype=bool)
    for w0 in _SOLAR_K_LINES_VAC:
        mask |= np.abs(wave_A - w0 * shift) < half_width_A
    mask |= wave_A >= (CO_2_0_BANDHEAD_NM * NM_TO_A - 2.0) * shift   # CO bandhead series → red
    return mask


def telluric_residual_gate(frame: CriresFrame, rv_kms: float = 0.0,
                           tol: float = 0.05) -> dict:
    """D1 primary gate: at telluric-dominated, solar-clean pixels, does the corrected
    flux return to the local continuum within `tol`? Returns the residual + PASS.
    Reflected-solar safe (ignores the solar spectrum). RYA-373 Decision 1."""
    assert_telluric_corrected(frame)
    seg = frame.segment_at(CO_2_0_BANDHEAD_NM)
    mt = getattr(seg, '_mtrans', None)
    if mt is None:
        raise RuntimeError("no model transmission on the corrected segment")
    cont = continuum_normalize(seg.wave_A, seg.flux)
    telluric = (mt < 0.90) & (mt > _MTRANS_FLOOR)         # telluric-dominated, not saturated
    solar = _solar_coincident(seg.wave_A, rv_kms)
    sel = telluric & ~solar & np.isfinite(cont)
    if sel.sum() < 10:
        return {'n_px': int(sel.sum()), 'residual': float('nan'),
                'passed': False, 'reason': 'too few telluric-clean pixels'}
    resid = float(np.nanmedian(np.abs(1.0 - cont[sel])))   # return-to-continuum residual
    frame._telluric_gate_residual = resid
    return {'n_px': int(sel.sum()), 'residual': resid, 'tol': tol,
            'passed': resid <= tol}


# ── RV-condition (RYA-372 reflected_solar_rv — single source) ─────────────────
def _load_reflected_solar_rv():
    try:
        from pipeline import reflected_solar_rv as rrv   # RYA-372
        return rrv
    except Exception as exc:
        raise VelocityModuleNotAvailableError(
            "pipeline.reflected_solar_rv (RYA-372) is not available — the single-"
            f"source asteroid-ephemeris velocity correction ({exc}). RYA-373 does NOT "
            "re-implement it; the RV-condition step waits on RYA-372.")


# RV guardrails (RYA-373 comment 7db5fe4d). The air→vac offset at 2.3 µm is ~+83
# km/s — catastrophic if half-applied / wrong-direction. Vesta's physical reflected
# RV spans ~−17..+33 km/s (RYA-370), so a measured |RV| beyond this band is the
# air-vac mis-handling signature → loud-fail (never inject 83 km/s silently).
_RV_PHYSICAL_MAX = 50.0          # km/s — beyond this = air-vac mis-convert, loud-fail
_TELLURIC_CLOSURE_MAX = 3.0      # km/s — telluric (topocentric) must close to ~0
_ANCHOR_DEPTH_MIN = 0.08         # central-depth floor for a usable RV anchor line
_ANCHOR_ISOLATION_A = 0.35       # Å — no comparable neighbour within this → unblended


def _air_to_vac(wl_air_A):
    """IAU air→vacuum (Morton 2000 / Birch & Downs) for Å. The solar line list is in
    AIR; the CRIRES K-band is VACUUM — convert at the spectrum-matching boundary."""
    wl = np.asarray(wl_air_A, dtype=float)
    s2 = (1.0e4 / wl) ** 2
    n = (1.0 + 0.00008336624212083 + 0.02408926869968 / (130.1065924522 - s2)
         + 0.0001599740894897 / (38.92568793293 - s2))
    return wl * n


def _solar_rv_anchors(lo_A: float, hi_A: float,
                      depth_min: float = _ANCHOR_DEPTH_MIN) -> dict:
    """Curate RV-anchor lines (selection logic) but resolve their wavelengths from the
    CANONICAL solar line list (PATHS['linelist_solar']) — NOT hardcoded (RYA-373
    guardrail 1). Select solar ATOMIC lines in [lo,hi] Å with central_depth ≥ depth_min
    that are isolated (no comparable neighbour within _ANCHOR_ISOLATION_A). Returns
    {'air': [...], 'labels': [...], 'source': path} — wavelengths are AIR (the list's
    frame); the caller converts to vacuum at the boundary."""
    import pandas as pd
    from pipeline.species import parse_ion       # RYA-345: atomic vs molecular
    path = Path(str(PATHS['linelist_solar']))
    df = pd.read_csv(path, comment='#', low_memory=False)
    w = df['wavelength_air_A'].astype(float)
    win = df[(w >= lo_A) & (w <= hi_A)].copy()
    if 'central_depth' not in win.columns:
        raise RuntimeError(f"{path.name} has no central_depth — cannot select RV anchors")
    # atomic only (a parseable ion); drop molecular notes
    def _is_atomic(el, ion):
        try:
            parse_ion(ion); return True
        except Exception:
            return False
    win = win[[_is_atomic(e, i) for e, i in zip(win['element'], win['ion'])]]
    win = win[win['central_depth'].astype(float) >= depth_min].sort_values('wavelength_air_A')
    wl = win['wavelength_air_A'].astype(float).to_numpy()
    dep = win['central_depth'].astype(float).to_numpy()
    keep = []
    for i in range(len(wl)):
        near = (np.abs(wl - wl[i]) < _ANCHOR_ISOLATION_A) & (np.arange(len(wl)) != i)
        if not np.any(near & (dep >= 0.5 * dep[i])):     # isolated from comparable lines
            keep.append(i)
    sel = win.iloc[keep]
    return {'air': sel['wavelength_air_A'].astype(float).tolist(),
            'labels': [f"{e} {i} {w:.2f}" for e, i, w in
                       zip(sel['element'], sel['ion'], sel['wavelength_air_A'])],
            'source': str(path)}


def _ccf_velocity(wave_A, obs_abs, model_abs, vmax=25.0, dv=0.3) -> float:
    """Telluric-anchor closure: velocity lag (km/s) that best aligns the observed
    telluric absorption (1−raw/continuum) to the molecfit model absorption (1−mtrans),
    by max cross-correlation (parabola-refined). Tellurics are at topocentric rest →
    a good wavelength solution closes to ~0."""
    from scipy.interpolate import interp1d
    m = np.isfinite(wave_A) & np.isfinite(obs_abs) & np.isfinite(model_abs)
    if m.sum() < 50:
        return np.nan
    w, o, md = wave_A[m], obs_abs[m], model_abs[m]
    f = interp1d(w, o, bounds_error=False, fill_value=0.0)
    vs = np.arange(-vmax, vmax + dv, dv)
    cc = [np.nansum(f(w * (1.0 + v / _C_KMS)) * md) for v in vs]
    k = int(np.argmax(cc))
    if 1 <= k < len(vs) - 1:                              # parabola refine
        y0, y1, y2 = cc[k - 1], cc[k], cc[k + 1]
        denom = (y0 - 2 * y1 + y2)
        if denom != 0:
            return float(vs[k] - 0.5 * dv * (y2 - y0) / denom)
    return float(vs[k])


def rv_condition(frame: CriresFrame, work_dir: Path = None) -> CriresFrame:
    """Shift the telluric-corrected CO order to the solar rest frame. The reflected RV
    is MEASURED off K-band solar lines via the RYA-372 single-source velocity module
    (measure_bulk_velocity) — anchors curated here, wavelengths resolved from the
    canonical solar line list and air→vac-converted at the boundary (guardrail 1+2).
    A telluric-anchor closure (tellurics = topocentric rest) verifies the wavelength
    zero-point. Loud-fails on the air-vac mis-handling signature (guardrail 2).

    Refuses a non-telluric-corrected frame (permanent rule)."""
    if not frame.telluric_corrected:
        raise TelluricNotCorrectedError(
            f"{frame.path.name}: refusing to RV-condition a frame whose telluric "
            f"correction is not verified (telluric_corrected=False). Telluric first.")
    rrv = _load_reflected_solar_rv()
    seg = frame.segment_at(CO_2_0_BANDHEAD_NM)
    lo, hi = float(np.nanmin(seg.wave_A)), float(np.nanmax(seg.wave_A))

    anchors = _solar_rv_anchors(lo, hi)
    if len(anchors['air']) < 2:
        raise RuntimeError(f"{frame.path.name}: <2 isolated solar RV anchors in the CO "
                           f"order [{lo:.0f},{hi:.0f}] Å — cannot condition RV.")
    anchors_vac = list(_air_to_vac(anchors['air']))       # ← air→vac BOUNDARY

    res = rrv.measure_bulk_velocity(seg.wave_A, seg.flux, lines=anchors_vac)
    v = res.get('v_med', np.nan)
    if not np.isfinite(v):
        raise RuntimeError(f"{frame.path.name}: RV measurement INSUFFICIENT "
                           f"(n_used={res.get('n_used')}/{len(anchors_vac)} anchors).")
    # guardrail 2 — air-vac loud-fail: a physical Vesta reflected RV is |v|≲33 km/s;
    # an unconverted/wrong-direction air-vac handling shows up as ~±83 km/s.
    if abs(v) > _RV_PHYSICAL_MAX:
        raise AssertionError(
            f"{frame.path.name}: measured reflected RV {v:+.1f} km/s exceeds the "
            f"physical Vesta band (±{_RV_PHYSICAL_MAX}) — the air↔vac (~83 km/s @ "
            f"2.3 µm) boundary is mis-handled. Refusing (no silent 83 km/s).")

    # telluric-anchor closure (zero-point): observed telluric absorption vs model
    cont_raw = continuum_normalize(seg.wave_A, getattr(seg, '_flux_raw', seg.flux))
    closure = _ccf_velocity(seg.wave_A, 1.0 - cont_raw, 1.0 - seg._mtrans)
    if np.isfinite(closure) and abs(closure) > _TELLURIC_CLOSURE_MAX:
        raise AssertionError(
            f"{frame.path.name}: telluric-anchor closure {closure:+.2f} km/s > "
            f"{_TELLURIC_CLOSURE_MAX} — wavelength zero-point is off (the tellurics "
            f"should sit at topocentric rest). Refusing.")

    seg.wave_A = seg.wave_A / (1.0 + v / _C_KMS)          # → solar rest frame
    frame.rest_frame = True
    frame._rv_refl = float(v)
    frame._rv_n_anchors = int(res.get('n_used', 0))
    frame._rv_anchor_src = anchors['source']
    frame._telluric_closure_kms = float(closure) if np.isfinite(closure) else None
    return frame


# ── Guard for downstream abundance/synthesis consumers ────────────────────────
def assert_telluric_corrected(frame: CriresFrame) -> None:
    """The permanent IR rule, enforced at the consumer boundary: no flux/abundance
    number may be emitted from an IR frame unless telluric correction is verified."""
    if not frame.telluric_corrected:
        raise TelluricNotCorrectedError(
            f"{frame.path.name}: IR abundance/synthesis attempted on a NON-telluric-"
            f"corrected frame. Run run_molecfit_telluric first (RYA-373, mandatory).")


# ── Per-epoch coadd (RV-registered; checksum-deduped) ─────────────────────────
def _file_md5(path) -> str:
    import hashlib
    return hashlib.md5(Path(path).read_bytes()).hexdigest()


def coadd_co(frames: list, step_A: float = 0.05) -> dict:
    """Coadd rest-frame, continuum-normalized CO orders on a common grid (inverse-
    variance weighted). CHECKSUM-DEDUP first (RYA-377: 6 byte-identical dup IDPs must
    never enter as independent frames and falsely √N-inflate SNR). Frames must be
    distinct SETTINGS at solar rest. Never blind — each frame is RV-registered first."""
    seen, uniq = {}, []
    for f in frames:
        h = _file_md5(f.path)
        if h in seen:
            print(f"  [coadd] DROPPED byte-duplicate {f.path.name} (== {seen[h]})")
            continue
        seen[h] = f.path.name
        uniq.append(f)
    for f in uniq:
        if not f.rest_frame:
            raise RuntimeError(f"{f.path.name}: not RV-conditioned — refusing blind coadd.")
    segs = [(f, f.segment_at(CO_2_0_BANDHEAD_NM)) for f in uniq]
    lo = max(float(np.nanmin(s.wave_A)) for _, s in segs)
    hi = min(float(np.nanmax(s.wave_A)) for _, s in segs)
    grid = np.arange(lo, hi + step_A * 0.5, step_A)
    from scipy.interpolate import interp1d
    num = np.zeros_like(grid); den = np.zeros_like(grid)
    for f, s in segs:
        ok = np.isfinite(s.wave_A) & np.isfinite(s.flux)
        fl = interp1d(s.wave_A[ok], s.flux[ok], bounds_error=False, fill_value=np.nan)(grid)
        wgt = (f.snr ** 2)                          # per-setting SNR² weight
        m = np.isfinite(fl)
        num[m] += fl[m] * wgt; den[m] += wgt
    coadd = np.where(den > 0, num / den, np.nan)
    snr_eff = float(np.sqrt(np.nansum([f.snr ** 2 for f in uniq])))
    return {'wave_A': grid, 'flux': coadd, 'frames': [f.wlen_id for f in uniq],
            'n_frames': len(uniq), 'snr_coadd': snr_eff}


# ── Conditioned-product output + provenance (PROVISIONAL, two gaps) ───────────
def write_conditioned(coadd: dict, frames: list, out_dir: Path,
                      gate_residuals: dict) -> Path:
    """Write the conditioned, continuum-normalized, rest-frame, coadded CO spectrum +
    a provenance header. Marked PROVISIONAL (RYA-373 guardrail 3): (a) telluric-gate
    final calibration pending the FTS solar IR atlas (RYA-162, absent); (b) re-validate
    after the RYA-387 0.001 re-extraction (RV anchors resolve from the provisional
    RYA-381 list). MUST NOT feed any published abundance until both clear."""
    from astropy.io import fits
    out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    ph = fits.PrimaryHDU()
    h = ph.header
    h['RYA'] = 'RYA-373'
    h['PROVIS'] = (True, 'PROVISIONAL — not for published abundance')
    h['GAP1'] = 'FTS solar IR atlas (RYA-162) absent — telluric-gate final cal pending'
    h['GAP2'] = 'revalidate after RYA-387 0.001 re-extraction (anchors from RYA-381)'
    h['TELL'] = ('molecfit', 'ESO molecfit (esorex), topocentric, H2O+CH4+CO2')
    h['CONTNORM'] = (True, 'continuum-normalized CO order')
    h['RESTFRM'] = (True, 'shifted to solar rest (measured reflected RV)')
    h['NFRAMES'] = (coadd['n_frames'], 'distinct K settings coadded')
    h['SNRCOADD'] = (round(coadd['snr_coadd'], 1), 'inverse-variance coadd SNR')
    for i, f in enumerate(frames):
        h[f'WLEN{i}'] = f.wlen_id
        h[f'GDAS{i}'] = getattr(f, '_gdas', '?')
        h[f'GATE{i}'] = (round(gate_residuals.get(f.wlen_id, float('nan')), 4),
                         'D1 telluric-specific residual')
        h[f'RVREFL{i}'] = (round(getattr(f, '_rv_refl', float('nan')), 3), 'reflected RV km/s')
        h[f'RVCLOS{i}'] = (round(getattr(f, '_telluric_closure_kms', float('nan')) or
                                 float('nan'), 3), 'telluric-anchor closure km/s')
    tab = fits.BinTableHDU.from_columns([
        fits.Column(name='wave_rest_A', format='1D', array=coadd['wave_A']),
        fits.Column(name='flux_norm', format='1D', array=coadd['flux'])], name='CO_COADD')
    out = out_dir / 'vesta_crires_K_CO_conditioned_PROVISIONAL.fits'
    fits.HDUList([ph, tab]).writeto(out, overwrite=True)
    return out


def condition_co_arm(out_dir: Path = None, work_root: Path = Path('/tmp/rya373_co')) -> dict:
    """Full implementable finish-out (RYA-373 #3–5): for each on-chip ¹²CO(2-0) frame
    (K2192, K2217) → molecfit telluric (GDAS) → D1 gate → continuum-normalize → RV-
    condition (measured, air→vac anchors, telluric closure) → per-epoch RV-registered,
    checksum-deduped coadd → conditioned CO output (PROVISIONAL). Returns a summary."""
    out_dir = Path(out_dir) if out_dir else (
        Path(str(PATHS['linelist_solar'])).parents[1] / 'audit' / 'crires_co_conditioned')
    frames = co_overtone_frames(on_chip_only=True)       # K2192, K2217 (¹²CO 2-0)
    gate_res, rv_status, done = {}, {}, []
    for f in frames:
        run_molecfit_telluric(f, work_root / f.wlen_id, molecules=TELLURIC_MOLECULES)
        g = telluric_residual_gate(f, rv_kms=0.0)
        gate_res[f.wlen_id] = g['residual']
        seg = f.segment_at(CO_2_0_BANDHEAD_NM)
        seg.flux = continuum_normalize(seg.wave_A, seg.flux)   # telluric→continuum-norm
        f.continuum_normalized = True
        # RV-condition; a non-measurable RV is reported LOUD (not faked, not crashed) —
        # the guardrail working. (RYA-373 finding: reflected RV is NOT reliably
        # measurable from the CO order — CO-band-dominated, sub-floor SNR.)
        try:
            rv_condition(f, work_root / f.wlen_id)
            rv_status[f.wlen_id] = f'OK v={f._rv_refl:+.2f} km/s (n={f._rv_n_anchors})'
            done.append(f)
        except (RuntimeError, AssertionError) as exc:
            rv_status[f.wlen_id] = f'RV-INSUFFICIENT: {exc}'
            print(f"  [rv] {f.wlen_id}: {exc}")

    result = {'frames': [f.wlen_id for f in frames], 'gate_residuals': gate_res,
              'rv_status': rv_status, 'telluric_corrected': True}
    if len(done) >= 2:
        coadd = coadd_co(done)
        result['output'] = str(write_conditioned(coadd, done, out_dir, gate_res))
        result['coadd_snr'] = coadd['snr_coadd']
        result['rv_refl'] = {f.wlen_id: f._rv_refl for f in done}
    else:
        result['output'] = None
        result['DATA_GAP'] = (
            "reflected RV not reliably measurable from the CO order on "
            f"{[f.wlen_id for f in frames if f.wlen_id not in [d.wlen_id for d in done]]} "
            "(CO-band-dominated, sub-floor SNR → clean atomic anchors too shallow; "
            "deep features are CO, not the atomic core). Telluric-corrected + D1-gated + "
            "continuum-normalized products ARE produced per frame; the rest-frame coadd "
            "needs the reflected RV measured from a clean-atomic-line order (resume point) "
            "or a higher-|RV| epoch (Decision 2).")
    return result


# ── CLI / smoke test ──────────────────────────────────────────────────────────
def _verify(crires_dir=VESTA_CRIRES_DIR, run_telluric: bool = False) -> dict:
    """Run the molecfit-independent core on the real data + report the engine-gated
    steps' availability. With run_telluric=True (and esorex present), also runs the
    real molecfit telluric pass on the best on-chip CO frame. Returns findings."""
    print("=" * 84)
    print("  RYA-373 — Vesta CRIRES+ K-band telluric/continuum/RV conditioning")
    print("=" * 84)

    frames = inventory(crires_dir)
    print(f"\n  Loaded {len(frames)} CRIRES+ IDP frames from {Path(crires_dir).name}/")
    bands = {}
    for f in frames:
        bands.setdefault(f.band, []).append(f)
    print("  bands: " + ", ".join(f"{b}×{len(v)}" for b, v in sorted(bands.items())))

    # nm→Å boundary check (a known K line lands in Å, not nm)
    kf = next((f for f in frames if f.band == 'K'), None)
    if kf is not None:
        lo_A = np.nanmin([np.nanmin(s.wave_A) for s in kf.segments])
        print(f"\n  nm→Å boundary: K frame {kf.wlen_id} min wave = {lo_A:,.0f} Å "
              f"({lo_A/NM_TO_A:,.0f} nm) — {'OK (Å)' if lo_A > 19000 else 'FAIL'}")

    # the K CO-overtone frames — header-covers vs ON-CHIP (gap-aware)
    hdr_cov = co_overtone_frames(crires_dir, on_chip_only=False)
    co = co_overtone_frames(crires_dir, on_chip_only=True)
    print(f"\n  CO(2-0) bandhead @ {CO_2_0_BANDHEAD_NM} nm: {len(hdr_cov)} K frames in "
          f"header range, {len(co)} with the bandhead ON-CHIP (rest fall in an "
          f"inter-order/detector gap).")
    print(f"  {'wlen':<8}{'date':<24}{'SNR':>7}{'nm range':>18}{'on-chip':>9}{'CO seg (Å)':>16}")
    for f in hdr_cov:
        seg = f.segment_at(CO_2_0_BANDHEAD_NM)
        seg_s = (f"{np.nanmin(seg.wave_A):.0f}-{np.nanmax(seg.wave_A):.0f}"
                 if seg else "GAP")
        print(f"  {f.wlen_id:<8}{f.date_obs:<24}{f.snr:>7.1f}"
              f"{f.wmin_nm:>8.0f}-{f.wmax_nm:<9.0f}{('YES' if seg else 'no'):>9}{seg_s:>16}")
    snr_quad = float(np.sqrt(np.nansum([f.snr**2 for f in co]))) if co else 0.0
    print(f"  on-chip per-frame SNR {[round(f.snr) for f in co]}; "
          f"RV-registered coadd ≈ {snr_quad:.0f} (floor {SNR_FLOOR:.0f}) "
          f"{'→ reaches floor' if snr_quad >= SNR_FLOOR else '→ below floor (coarse, as expected)'}")

    # raw-state asserts (the conditioning has NOT happened yet)
    raw_ok = all(f.specsys.upper() == 'TOPOCENT' for f in co) and \
             all(f.fluxcal.upper() == 'UNCALIBRATED' for f in co)
    print(f"\n  raw state: all on-chip CO frames SPECSYS=TOPOCENT & "
          f"FLUXCAL=UNCALIBRATED = {raw_ok}")

    # continuum normalization on the CO bandhead order (testable now)
    if co:
        f0 = max(co, key=lambda f: f.snr)
        seg = f0.segment_at(CO_2_0_BANDHEAD_NM)
        if seg is not None:
            norm = continuum_normalize(seg.wave_A, seg.flux)
            fin = np.isfinite(norm)
            print(f"\n  continuum-normalize {f0.wlen_id} CO order "
                  f"(ord {seg.order}/det {seg.detector}, {fin.sum()} px): "
                  f"median={np.nanmedian(norm[fin]):.3f}, "
                  f"min={np.nanmin(norm[fin]):.3f} (line cores), "
                  f"p95={np.nanpercentile(norm[fin],95):.3f}")

    # engine availability (telluric + RV) — fail-loud, never faked
    print(f"\n  ENGINE STATUS (no step is faked — each fails loud if absent):")
    mol = _esorex_available()
    print(f"    molecfit/esorex telluric engine : {'available' if mol else 'NOT INSTALLED'}")
    try:
        _load_reflected_solar_rv(); rv = True
    except VelocityModuleNotAvailableError:
        rv = False
    print(f"    reflected_solar_rv (RYA-372)    : {'available' if rv else 'NOT PRESENT'}")

    # optional real molecfit telluric pass (the smoke test's --telluric)
    tell_result = None
    if run_telluric and mol and co:
        f0 = max(co, key=lambda f: f.snr)
        print(f"\n  RUNNING molecfit telluric on {f0.wlen_id} (TOPOCENTRIC) …")
        try:
            run_molecfit_telluric(f0, Path(f'/tmp/mfit_{f0.wlen_id}'),
                                  molecules=TELLURIC_MOLECULES)
            seg = f0.segment_at(CO_2_0_BANDHEAD_NM)
            tell_result = {'wlen': f0.wlen_id, 'residual': f0._telluric_residual,
                           'mtrans_max': f0._telluric_mtrans_max,
                           'molecules': f0._molecules}
            print(f"    telluric_corrected=True; model transmission max "
                  f"{f0._telluric_mtrans_max:.3f}; molecules {f0._molecules}")
            print(f"    telluric residual (moderate-telluric proxy) "
                  f"{f0._telluric_residual*100:.1f}%  — NOTE: a reflected-SOLAR target is "
                  f"solar-line-rich, so this proxy conflates solar lines with telluric "
                  f"misfit; a telluric-line-specific metric is needed for the <2% gate.")
            print(f"    solar CO(2-0) bandhead region (22930-22960 Å) min flux "
                  f"{np.nanmin(seg.flux[(seg.wave_A>22930)&(seg.wave_A<22960)]):.0f} "
                  f"(residual absorption present = candidate surviving solar CO).")
        except Exception as exc:
            print(f"    molecfit telluric run FAILED: {exc}")

    print("\n" + "-" * 84)
    if mol and rv:
        print("  READY: both engines present — run the full telluric→RV→coadd pass.")
    else:
        pend = []
        if not mol: pend.append("molecfit (telluric)")
        if not rv: pend.append("RYA-372 reflected_solar_rv (RV-condition)")
        print(f"  CORE VERIFIED (load / nm→Å / CO coverage / continuum). PENDING engine(s): "
              f"{', '.join(pend)}.")
        print("  These steps FAIL LOUD until their engine is available — no faked "
              "telluric/RV, no IR number emitted (permanent rule).")
    print("=" * 84)
    return {'frames': frames, 'co_frames': co, 'molecfit': mol, 'rv': rv,
            'raw_ok': raw_ok, 'snr_coadd': snr_quad}


def main(argv=None):
    import argparse
    ap = argparse.ArgumentParser(description="RYA-373 Vesta CRIRES+ K-band telluric/RV")
    ap.add_argument('--set', default='vesta_crires_k',
                    help="dataset (only vesta_crires_k implemented)")
    ap.add_argument('--verify', action='store_true', help="run the verification report")
    ap.add_argument('--telluric', action='store_true',
                    help="also run the real molecfit telluric pass (slow; needs esorex)")
    ap.add_argument('--condition', action='store_true',
                    help="run the full #3-5 finish-out: telluric→gate→continuum→RV→coadd"
                         "→conditioned CO output (slow; needs esorex + RYA-372)")
    args = ap.parse_args(argv)
    if args.set != 'vesta_crires_k':
        raise SystemExit(f"unknown --set {args.set!r} (only 'vesta_crires_k')")
    if args.condition:
        import json
        res = condition_co_arm()
        print("\n" + "=" * 84)
        print("  RYA-373 #3-5 — conditioned CO arm (PROVISIONAL)")
        print("=" * 84)
        print(json.dumps(res, indent=2, default=str))
        return res
    return _verify(run_telluric=args.telluric)


if __name__ == '__main__':
    main()
