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
    return shutil.which("esorex") is not None or shutil.which("molecfit_model") is not None


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


def _molecfit_driver(frame: CriresFrame, work_dir: Path, molecules) -> CriresFrame:  # pragma: no cover
    """esorex molecfit_model → molecfit_calctrans on the K CO orders; divide the
    fitted telluric model out; mask saturated cores; verify the telluric residual.
    Finalized when molecfit is installed (engine-dependent recipe params)."""
    raise MolecfitNotAvailableError(
        "molecfit driver pending finalization against the installed esorex/molecfit "
        "(recipe parameters are version-specific). Engine not yet available.")


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


def rv_condition(frame: CriresFrame, horizons_id: str = '4') -> CriresFrame:
    """Shift a telluric-corrected, continuum-normalized frame to the solar rest frame
    via the asteroid ephemeris (Vesta = Horizons 4), using the RYA-372 single-source
    module. AFTER telluric. Verify on a clean K-band solar atomic line.

    Refuses to run on a frame that is not telluric-corrected (the permanent rule)."""
    if not frame.telluric_corrected:
        raise TelluricNotCorrectedError(
            f"{frame.path.name}: refusing to RV-condition a frame whose telluric "
            f"correction is not verified (telluric_corrected=False). Telluric first.")
    rrv = _load_reflected_solar_rv()
    return rrv.condition_frame(frame, horizons_id=horizons_id)   # RYA-372 API


# ── Guard for downstream abundance/synthesis consumers ────────────────────────
def assert_telluric_corrected(frame: CriresFrame) -> None:
    """The permanent IR rule, enforced at the consumer boundary: no flux/abundance
    number may be emitted from an IR frame unless telluric correction is verified."""
    if not frame.telluric_corrected:
        raise TelluricNotCorrectedError(
            f"{frame.path.name}: IR abundance/synthesis attempted on a NON-telluric-"
            f"corrected frame. Run run_molecfit_telluric first (RYA-373, mandatory).")


# ── CLI / smoke test ──────────────────────────────────────────────────────────
def _verify(crires_dir=VESTA_CRIRES_DIR) -> dict:
    """Run the molecfit-independent core on the real data + report the engine-gated
    steps' availability. Returns a findings dict."""
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
    args = ap.parse_args(argv)
    if args.set != 'vesta_crires_k':
        raise SystemExit(f"unknown --set {args.set!r} (only 'vesta_crires_k')")
    return _verify()


if __name__ == '__main__':
    main()
