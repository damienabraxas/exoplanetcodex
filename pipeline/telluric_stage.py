"""
pipeline/telluric_stage.py
==========================
RYA-424 — Telluric correction as a STANDING, instrument-aware data-input stage.

This promotes the per-dataset one-offs (RYA-373 molecfit driver on Vesta CRIRES+,
RYA-380 per-night-GDAS recipe) into a single gate that EVERY red-optical / IR
dataset passes through on data input, with a hard verification check before a
spectrum can be marked analysis-ready. It is upstream of the solar full-spectrum
run's Phase B (the 2.3 µm 12C/13C CO measurement), of Procyon, and of alpha Cen —
not parallel to them.

The four moving parts (RYA-424 §3):

  1. WAVELENGTH GATE (not instrument gate). On ingest of any spectrum, the
     [λ_lo, λ_hi] span is classified:
        * λ ≳ 6800 Å (red-optical + IR): sharp telluric forest → correction MANDATORY.
        * λ ≲ 3800 Å (blue) and space-UV (HST/STIS/COS): NO telluric step.
        * 3800–6800 Å (mid-optical): largely clean; case-by-case (not auto-required).
     See `telluric_regime` / `requires_telluric`.

  2. ENGINE ROUTING (single source of truth, no silent default). The instrument keys
     the engine — molecfit+GDAS for CRIRES+/UVES-red/ESPRESSO-red/FEROS-red/CHIRON/
     NIRPS; APERO+Wapiti for SPIRou (permanent rule — SPIRou is NOT molecfit). An
     unknown instrument LOUD-FAILS (`select_engine`).

  3. PER-NIGHT GDAS, LOUD-FAIL ON SILENT FALLBACK. molecfit's GDAS_PROFILE=auto
     silently falls back to a generic standard atmosphere when the real 3-hourly
     observation-night profile is missing (the RYA-373 failure mode). This stage
     resolves the real per-night profile up front and refuses to proceed without it,
     and re-asserts the actual profile used post-run (`resolve_gdas_profile`,
     `assert_real_gdas`).

  4. VERIFICATION + analysis_ready FLAG. Residuals in known telluric windows are
     scored against a tolerance (`telluric_residual_metric`, `TELLURIC_RESIDUAL_TOL`).
     A frame that does not pass is FLAGGED, not silently passed. The verdict
     (telluric_verified true/false + residual metric + GDAS provenance + engine) is
     written to a conditioning manifest (`TelluricManifest`, `write_manifest`). The
     abundance path refuses red-optical/IR input that is not telluric_verified
     (`require_telluric_verified`).

Permanent rules carried from RYA-373/380: loud failure over silent fallback; never
coadd across settings/resolutions; Ångström throughout; single source of truth for
engine selection + GDAS provenance.

Smoke test:  python -m pipeline.telluric_stage --vesta-crires
"""
from __future__ import annotations

import json
import shutil
from dataclasses import asdict, dataclass, field
from pathlib import Path

import numpy as np

from config import constants as C

# ── Errors (every one is a LOUD failure — never a silent skip / fake correction) ──
class EngineNotSelectedError(RuntimeError):
    """The instrument is not in the engine registry. The telluric engine is chosen
    explicitly per instrument (single source of truth) — there is NO silent default;
    an unknown instrument stops the stage rather than guessing molecfit-vs-Wapiti."""


class TelluricEngineNotInstalledError(RuntimeError):
    """The routed engine (molecfit / APERO+Wapiti) is not installed. The correction is
    mandatory and cannot be faked — install the engine and re-run."""


class GDASUnavailableError(RuntimeError):
    """The real per-night GDAS atmospheric profile could not be obtained (or molecfit
    silently fell back to a generic standard atmosphere). The RYA-373 failure mode —
    LOUD-FAIL, never silently pass a standard-profile correction as per-night."""


class TelluricNotVerifiedError(RuntimeError):
    """A red-optical/IR spectrum reached the abundance path without a verified telluric
    correction. The permanent rule: no red-optical/IR number without telluric_verified."""


# ── 1. Wavelength gate (the rule is wavelength-gated, not instrument-gated) ────────
def telluric_regime(wave_lo_A: float, wave_hi_A: float) -> str:
    """Classify a spectral span [lo, hi] Å into a telluric regime (RYA-380 rule):

        ir          — entirely ≳ TELLURIC_LAMBDA_MIN_A (6800 Å): sharp telluric forest.
        red_optical — crosses into / sits in the 6800–~10000 Å telluric forest.
        mid_optical — 3800–6800 Å: largely clean, case-by-case (NOT auto-required).
        blue        — entirely ≲ TELLURIC_BLUE_MAX_A (3800 Å), ground-based.

    A span that straddles the 6800 Å threshold is `red_optical` (it has telluric forest
    in its red end and must be corrected). Space-UV is NOT inferred from wavelength
    (HST/STIS/COS sit blueward but are above the atmosphere) — pass
    `space_uv=True` to `requires_telluric` for those.
    """
    lo, hi = float(wave_lo_A), float(wave_hi_A)
    if hi < lo:
        lo, hi = hi, lo
    if lo >= C.TELLURIC_LAMBDA_MIN_A:
        return C.TELLURIC_REGIME_IR
    if hi > C.TELLURIC_LAMBDA_MIN_A:
        return C.TELLURIC_REGIME_RED_OPTICAL
    if hi <= C.TELLURIC_BLUE_MAX_A:
        return C.TELLURIC_REGIME_BLUE
    return C.TELLURIC_REGIME_MID_OPTICAL


def requires_telluric(wave_lo_A: float, wave_hi_A: float,
                      space_uv: bool = False) -> bool:
    """True iff this span MUST pass telluric correction before analysis-ready (RYA-380):
    the IR and red-optical regimes (any forest above 6800 Å). Blue/mid-optical are not
    auto-required; space-UV (above the atmosphere) is never required."""
    if space_uv:
        return False
    return telluric_regime(wave_lo_A, wave_hi_A) in (
        C.TELLURIC_REGIME_IR, C.TELLURIC_REGIME_RED_OPTICAL)


# ── 2. Engine routing (single source of truth; no silent default) ─────────────────
def select_engine(instrument: str) -> str:
    """The telluric engine for an instrument (`constants.TELLURIC_ENGINES`). LOUD-FAIL
    on an unknown instrument — the engine is never guessed (molecfit vs APERO+Wapiti is
    a physics decision, not a default)."""
    eng = C.TELLURIC_ENGINES.get(instrument)
    if eng is None:
        raise EngineNotSelectedError(
            f"no telluric engine registered for instrument {instrument!r}. The engine is "
            f"instrument-keyed and explicit (no silent default). Known: "
            f"{sorted(C.TELLURIC_ENGINES)}. Add it to constants.TELLURIC_ENGINES "
            f"(molecfit+GDAS or apero_wapiti) before conditioning this data.")
    return eng


def site_for_instrument(instrument: str) -> dict:
    """The observatory site record (lat/lon/elev/GDAS loc) for an instrument's telluric
    engine. LOUD-FAIL if unregistered — the per-night GDAS profile is per-site."""
    site_key = C.TELLURIC_INSTRUMENT_SITE.get(instrument)
    if site_key is None or site_key not in C.SITES:
        raise EngineNotSelectedError(
            f"no observatory site registered for instrument {instrument!r} — cannot "
            f"locate the per-night GDAS profile. Add it to "
            f"constants.TELLURIC_INSTRUMENT_SITE + constants.SITES.")
    return dict(C.SITES[site_key], key=site_key)


def engine_available(engine: str) -> bool:
    """Is the routed engine installed? molecfit = esorex/molecfit on PATH; apero_wapiti
    = the APERO + Wapiti Python stack importable. Used by the stage to fail loud with a
    clear 'install X' message rather than a downstream crash."""
    if engine == C.TELLURIC_ENGINE_MOLECFIT:
        return (shutil.which('esorex') is not None
                or Path('/opt/homebrew/bin/esorex').exists()
                or shutil.which('molecfit_model') is not None)
    if engine == C.TELLURIC_ENGINE_APERO_WAPITI:
        try:
            import apero  # noqa: F401
            import wapiti  # noqa: F401
            return True
        except Exception:
            return False
    return False


# ── 3. Per-night GDAS resolver (LOUD-FAIL on silent fallback) ──────────────────────
def _telluriccorr_share() -> "Path | None":
    """The installed telluriccorr GDAS profile directory (Homebrew layout), or None."""
    cellar = Path('/opt/homebrew/Cellar/telluriccorr')
    if not cellar.exists():
        return None
    for ver in sorted(cellar.iterdir(), reverse=True):
        d = ver / 'share' / 'molecfit' / 'data' / 'profiles' / 'gdas'
        if d.exists():
            return d
    return None


def gdas_tarball_for_site(site: dict) -> Path:
    """Path to the per-site GDAS tarball (molecfit ships per-location archives keyed by
    the site's `gdas_loc`, e.g. 'C-70.4-24.6' for Paranal). LOUD-FAIL if the molecfit
    GDAS install or the site tarball is absent — never fall back to a standard atmosphere."""
    loc = site.get('gdas_loc')
    if loc is None:
        raise GDASUnavailableError(
            f"site {site.get('key')!r} has no molecfit GDAS location — its engine is not "
            f"molecfit (per-night GDAS does not apply). This is a routing error.")
    share = _telluriccorr_share()
    if share is None:
        raise GDASUnavailableError(
            "telluriccorr GDAS profile directory not found (expected under "
            "/opt/homebrew/Cellar/telluriccorr/*/share/molecfit/data/profiles/gdas). "
            "Install the ESO molecfit/telluriccorr GDAS data — no standard-profile fallback.")
    tb = share / f'gdas_profiles_{loc}.tar.gz'
    if not tb.exists():
        raise GDASUnavailableError(
            f"per-site GDAS tarball missing for {site.get('name')} ({loc}): {tb}. The "
            f"per-night profile cannot be resolved — LOUD-FAIL (no standard-atm fallback).")
    return tb


def resolve_gdas_profile(mjd: float, instrument: str, work_dir: Path) -> str:
    """Resolve the REAL nearest-3-hourly GDAS profile for the observation MJD at the
    instrument's site, as a molecfit-ready FITS path. LOUD-FAIL (GDASUnavailableError)
    if no real per-night profile is available — the whole point of RYA-380 is that a
    silent standard-atmosphere fallback is never accepted as 'per-night'.

    Implementation reuses the RYA-373 tarball-extract + ASCII→FITS mechanic, but its
    None-return (which let molecfit fall back) is upgraded here to a hard failure."""
    if np.isnan(mjd):
        raise GDASUnavailableError(
            f"{instrument}: NaN observation MJD — cannot locate a per-night GDAS profile "
            f"(BERV/GDAS both depend on a real header time). Fix the frame time first.")
    site = site_for_instrument(instrument)
    tarball = gdas_tarball_for_site(site)            # loud-fails if absent
    from pipeline import crires_telluric as ct
    # Point the RYA-373 extractor at the per-site tarball + loc, then resolve.
    orig_tb, orig_loc = ct._GDAS_TARBALL, ct._GDAS_LOC
    try:
        ct._GDAS_TARBALL, ct._GDAS_LOC = str(tarball), site['gdas_loc']
        gdas = ct._nearest_gdas(float(mjd), Path(work_dir))
    finally:
        ct._GDAS_TARBALL, ct._GDAS_LOC = orig_tb, orig_loc
    if not gdas:
        raise GDASUnavailableError(
            f"{instrument} MJD {mjd:.5f}: no real 3-hourly GDAS profile within ±6 h in "
            f"{tarball.name}. RYA-380 forbids the silent standard-atmosphere fallback — "
            f"LOUD-FAIL. (Fetch the observation-night GDAS for this epoch and re-run.)")
    return gdas


_STANDARD_PROFILE_SENTINELS = ('standard-profile', 'no gdas', 'standard atm')


def assert_real_gdas(provenance: str, instrument: str = '') -> None:
    """Post-run assertion: the GDAS provenance string recorded by the engine must name a
    real per-night profile, not the standard-atmosphere sentinel. Catches the case where
    the engine silently fell back AFTER resolution (RYA-373 critical failure mode)."""
    p = (provenance or '').lower()
    if (not provenance) or any(s in p for s in _STANDARD_PROFILE_SENTINELS):
        raise GDASUnavailableError(
            f"{instrument}: telluric correction used GDAS provenance {provenance!r} — a "
            f"generic standard atmosphere, NOT the observation-night GDAS. RYA-380 LOUD-"
            f"FAIL: this frame is not analysis-ready on a standard-profile correction.")


# ── 4. Verification metric (residuals in known telluric windows) ──────────────────
def telluric_residual_metric(wave_A: np.ndarray, corrected_flux: np.ndarray,
                             model_transmission: np.ndarray,
                             science_mask: "np.ndarray | None" = None,
                             continuum: "np.ndarray | None" = None,
                             mtrans_hi: float = 0.90, mtrans_floor: float = 0.02,
                             tol: float = None,
                             window_A: "tuple | None" = None) -> dict:
    """The single telluric VERIFICATION metric (RYA-424 §3.3, generalizing the RYA-373
    D1 gate). At pixels that are (a) telluric-DOMINATED in the engine's model
    (mtrans_floor < transmission < mtrans_hi — absorbing but not saturated) and (b) NOT
    science-coincident (`science_mask` masks the wanted stellar lines so a line-rich
    target is not scored as telluric misfit), the telluric-corrected flux must return
    to the local continuum. The metric is the median |1 − corrected/continuum| there.

    `window_A=(lo,hi)` (RYA-437 Part B): restrict the score to a wavelength window —
    telluric quality is wavelength-dependent, so for CO (12C/13C, C/O) science the gate
    scores LOCAL to the CO bandheads (TELLURIC_CO_LOCAL_WINDOW_A) while the global
    median (window_A=None) is kept as a secondary report.

    Returns {n_px, residual, tol, passed, window_A}. `passed` requires enough clean
    telluric pixels (≥10) AND residual ≤ tol; too-few-pixels is reported (not passed)."""
    tol = C.TELLURIC_RESIDUAL_TOL if tol is None else tol
    wave_A = np.asarray(wave_A, float)
    corrected_flux = np.asarray(corrected_flux, float)
    mt = np.asarray(model_transmission, float)
    if continuum is None:
        cont = _continuum_normalize(wave_A, corrected_flux)
    else:
        with np.errstate(invalid='ignore', divide='ignore'):
            cont = corrected_flux / np.asarray(continuum, float)
    telluric = (mt > mtrans_floor) & (mt < mtrans_hi)
    sel = telluric & np.isfinite(cont)
    if science_mask is not None:
        sel &= ~np.asarray(science_mask, bool)
    if window_A is not None:
        lo, hi = float(min(window_A)), float(max(window_A))
        sel &= (wave_A >= lo) & (wave_A <= hi)
    n = int(sel.sum())
    if n < 10:
        return {'n_px': n, 'residual': float('nan'), 'tol': float(tol),
                'passed': False, 'window_A': window_A,
                'reason': 'too few telluric-clean pixels to verify'}
    resid = float(np.nanmedian(np.abs(1.0 - cont[sel])))
    return {'n_px': n, 'residual': resid, 'tol': float(tol), 'passed': resid <= tol,
            'window_A': window_A}


def co_local_residual_metric(wave_A: np.ndarray, corrected_flux: np.ndarray,
                             model_transmission: np.ndarray,
                             science_mask: "np.ndarray | None" = None,
                             continuum: "np.ndarray | None" = None,
                             tol: float = None) -> dict:
    """RYA-437 Part B: the CO-science verification metric. Reports the residual scored
    LOCAL to the CO (2-0) bandheads (TELLURIC_CO_LOCAL_WINDOW_A, 2.2935→2.3448 µm) as the
    PRIMARY verdict, plus the GLOBAL median as a secondary. The verdict (`passed`,
    `residual`) is the CO-local one — telluric quality at the CO bandheads is what biases
    12C/13C, not the segment-wide median (which can pass while the CO region is dirty)."""
    co = telluric_residual_metric(wave_A, corrected_flux, model_transmission,
                                  science_mask=science_mask, continuum=continuum,
                                  tol=tol, window_A=C.TELLURIC_CO_LOCAL_WINDOW_A)
    glob = telluric_residual_metric(wave_A, corrected_flux, model_transmission,
                                    science_mask=science_mask, continuum=continuum,
                                    tol=tol, window_A=None)
    return {'residual': co['residual'], 'passed': co['passed'], 'tol': co['tol'],
            'n_px': co['n_px'], 'reason': co.get('reason'),
            'metric_used': 'co_local', 'window_A': C.TELLURIC_CO_LOCAL_WINDOW_A,
            'residual_co_local': co['residual'], 'residual_global': glob['residual'],
            'n_px_co_local': co['n_px'], 'n_px_global': glob['n_px']}


def _continuum_normalize(wave_A, flux):
    """Thin re-export of the RYA-373 continuum normalizer (one source). Imported lazily
    so the gate/engine/manifest logic stays importable without astropy/scipy."""
    from pipeline.crires_telluric import continuum_normalize
    return continuum_normalize(np.asarray(wave_A, float), np.asarray(flux, float))


# ── Conditioning manifest (the analysis_ready record) ─────────────────────────────
@dataclass
class TelluricManifest:
    """The per-frame telluric conditioning record written at data input. `telluric_verified`
    + `residual` + `gdas_profile` + `engine` is the analysis_ready contract: the abundance
    path reads this (require_telluric_verified) and refuses red-optical/IR input that is
    not verified."""
    dataset: str                     # logical dataset id (e.g. 'vesta_crires_k')
    source_path: str                 # the conditioned frame's source file
    instrument: str
    engine: str                      # molecfit | apero_wapiti
    site: str                        # observatory key (paranal/la_silla/ctio/cfht)
    wave_lo_A: float
    wave_hi_A: float
    regime: str                      # ir | red_optical | mid_optical | blue | space_uv
    telluric_required: bool          # did the wavelength gate require correction?
    telluric_verified: bool          # passed correction + GDAS + residual gate?
    residual: "float | None"         # the verification metric the verdict used (CO-local for CO science)
    tolerance: float                 # TELLURIC_RESIDUAL_TOL (RYA-437: 13C/CO-derived)
    n_verify_px: int                 # telluric-clean pixels the residual was scored on
    gdas_profile: str                # per-night GDAS provenance (real profile name)
    mjd: float = float('nan')
    # RYA-437: which metric drove the verdict + both numbers (CO-local primary, global secondary).
    metric_used: str = 'global'      # 'co_local' for CO science, else 'global'
    residual_co_local: "float | None" = None
    residual_global: "float | None" = None
    rya: str = 'RYA-424'
    provisional: bool = False        # carry RYA-373 PROVISIONAL where applicable
    notes: list = field(default_factory=list)

    @property
    def analysis_ready(self) -> bool:
        """A frame is analysis-ready iff telluric is not required, OR it is required AND
        verified. (Blue/space-UV/mid-optical that don't require telluric pass through.)"""
        return (not self.telluric_required) or self.telluric_verified


def write_manifest(manifest: TelluricManifest, path: Path) -> Path:
    """Persist one conditioning manifest as JSON (the conditioning record). Sibling
    `<dataset>.telluric.json` by convention."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    d = asdict(manifest)
    d['analysis_ready'] = manifest.analysis_ready
    path.write_text(json.dumps(d, indent=2, default=str))
    return path


def load_manifest(path: Path) -> dict:
    """Read a conditioning manifest JSON back to a dict (the abundance-path guard reads
    this; kept as a plain dict so a manifest written by a future engine still loads)."""
    return json.loads(Path(path).read_text())


# ── Abundance-path guard (the permanent rule, enforced at the consumer boundary) ──
def require_telluric_verified(manifest, wave_lo_A: float = None,
                              wave_hi_A: float = None, space_uv: bool = False) -> None:
    """Called by the abundance/synthesis path before consuming a red-optical/IR frame.
    Accepts a `TelluricManifest`, a manifest dict, or a path to a manifest JSON. If the
    span requires telluric (or the manifest says it was required) and the manifest is
    not telluric_verified, LOUD-FAIL (TelluricNotVerifiedError). A frame that does not
    require telluric (blue / space-UV / mid-optical) passes through."""
    if isinstance(manifest, (str, Path)):
        manifest = load_manifest(manifest)
    if isinstance(manifest, TelluricManifest):
        required = manifest.telluric_required
        verified = manifest.telluric_verified
        name = manifest.source_path or manifest.dataset
    elif isinstance(manifest, dict):
        required = bool(manifest.get('telluric_required'))
        verified = bool(manifest.get('telluric_verified'))
        name = manifest.get('source_path') or manifest.get('dataset', '<unknown>')
    else:
        raise TypeError(f"require_telluric_verified: unsupported manifest {type(manifest)!r}")
    # If the caller supplies the span, recompute the gate as a cross-check (defends
    # against a stale/wrong telluric_required in the manifest).
    if wave_lo_A is not None and wave_hi_A is not None:
        required = required or requires_telluric(wave_lo_A, wave_hi_A, space_uv=space_uv)
    if required and not verified:
        raise TelluricNotVerifiedError(
            f"{name}: red-optical/IR abundance input is NOT telluric_verified. The "
            f"permanent rule (RYA-380/424): no red-optical/IR number without a verified "
            f"molecfit/Wapiti correction against the observation-night GDAS. Run the "
            f"telluric stage (condition_crires_frame / the engine driver) first.")


# ── The data-input gate, wired end-to-end for the CRIRES+ (molecfit) engine ───────
def condition_crires_frame(frame, work_dir: Path, dataset: str = 'crires',
                           instrument: str = 'CRIRES+', rv_kms: float = 0.0,
                           molecules=None) -> TelluricManifest:
    """Run the full standing stage on ONE CRIRES+ frame (the RYA-373 first case), the
    molecfit engine end-to-end: route → resolve per-night GDAS (loud-fail) → molecfit
    telluric (topocentric) → assert real GDAS post-run → verify residual → build the
    conditioning manifest. Returns the `TelluricManifest` (NOT yet written — caller
    persists with `write_manifest`).

    Delegates the actual molecfit correction to the single-source RYA-373 driver
    (pipeline.crires_telluric) rather than re-implementing it."""
    from pipeline import crires_telluric as ct
    lo_A = float(np.nanmin([np.nanmin(s.wave_A) for s in frame.segments]))
    hi_A = float(np.nanmax([np.nanmax(s.wave_A) for s in frame.segments]))
    regime = telluric_regime(lo_A, hi_A)
    required = requires_telluric(lo_A, hi_A)
    engine = select_engine(instrument)
    site = site_for_instrument(instrument)
    if engine != C.TELLURIC_ENGINE_MOLECFIT:
        raise EngineNotSelectedError(
            f"{instrument}: condition_crires_frame is the molecfit path, but {instrument} "
            f"routes to {engine!r}. Use the {engine} driver.")
    if not engine_available(engine):
        raise TelluricEngineNotInstalledError(
            f"molecfit/esorex not installed — the {instrument} telluric engine is absent. "
            f"Install the ESO molecfit pipeline and re-run (no silent fallback).")
    molecules = tuple(molecules) if molecules else ct.TELLURIC_MOLECULES

    notes = []
    # Pre-flight: a real per-night GDAS profile MUST exist before we run (loud-fail).
    resolve_gdas_profile(frame.mjd, instrument, Path(work_dir) / 'gdas_preflight')

    # Correction via the single-source RYA-373 molecfit driver (topocentric).
    ct.run_molecfit_telluric(frame, Path(work_dir), molecules=molecules)
    gdas_prov = getattr(frame, '_gdas', '')
    assert_real_gdas(gdas_prov, instrument)          # no silent standard-profile fallback

    # Verification (RYA-437): for CO science the verdict is the CO-region-LOCAL residual
    # (telluric quality at the 2-0 bandheads, what biases 12C/13C), with the global median
    # kept as a secondary report. Both scored at telluric-dominated, solar-clean pixels.
    seg = frame.segment_at(ct.CO_2_0_BANDHEAD_NM)
    mt = getattr(seg, '_mtrans', None)
    science_mask = ct._solar_coincident(seg.wave_A, rv_kms) if mt is not None else None
    if mt is None:
        verify = {'n_px': 0, 'residual': float('nan'),
                  'tol': C.TELLURIC_RESIDUAL_TOL, 'passed': False,
                  'metric_used': 'co_local', 'residual_co_local': float('nan'),
                  'residual_global': float('nan'),
                  'reason': 'no model transmission on corrected segment'}
    else:
        verify = co_local_residual_metric(seg.wave_A, seg.flux, mt,
                                          science_mask=science_mask)
    if not verify['passed']:
        notes.append(f"VERIFY-FAIL: {verify.get('reason', 'CO-local residual above tolerance')} "
                     f"(residual_co_local={verify['residual_co_local']}, "
                     f"residual_global={verify['residual_global']}, tol={verify['tol']}, "
                     f"n_px={verify['n_px']})")
    notes.append('RYA-437: tol is 13C/CO-derived (binding 13CO 2-0); verdict on the CO-'
                 'region-local residual. RYA-373 PROVISIONAL: FTS solar IR atlas (RYA-162) '
                 'absent; re-validate after RYA-387 0.001 re-extract.')

    def _f(x):
        return None if (x is None or not np.isfinite(x)) else float(x)
    return TelluricManifest(
        dataset=dataset, source_path=str(frame.path), instrument=instrument,
        engine=engine, site=site['key'], wave_lo_A=lo_A, wave_hi_A=hi_A,
        regime=regime, telluric_required=required,
        telluric_verified=bool(verify['passed']),
        residual=_f(verify['residual']),
        tolerance=float(verify['tol']), n_verify_px=int(verify['n_px']),
        gdas_profile=str(gdas_prov), mjd=float(frame.mjd),
        metric_used=verify.get('metric_used', 'co_local'),
        residual_co_local=_f(verify.get('residual_co_local')),
        residual_global=_f(verify.get('residual_global')),
        provisional=True, notes=notes)


def condition_vesta_crires(out_dir: Path = None,
                           work_root: Path = Path('/tmp/rya424_telluric')) -> dict:
    """Smoke-driver: run the standing stage on the Vesta CRIRES+ on-chip CO frames end
    to end (the RYA-373 first case) and write one conditioning manifest per frame. This
    is the acceptance smoke test: extracted → telluric-corrected → verified → flagged."""
    from pipeline import crires_telluric as ct
    out_dir = Path(out_dir) if out_dir else (
        Path(str(C.PATHS['linelist_solar'])).parents[1] / 'audit' / 'telluric_stage_rya424')
    frames = ct.co_overtone_frames(on_chip_only=True)
    results = {}
    for f in frames:
        man = condition_crires_frame(f, work_root / f.wlen_id,
                                     dataset='vesta_crires_k', instrument='CRIRES+')
        path = write_manifest(man, out_dir / f'vesta_crires_{f.wlen_id}.telluric.json')
        results[f.wlen_id] = {'manifest': str(path), 'telluric_verified': man.telluric_verified,
                              'residual': man.residual, 'gdas': man.gdas_profile,
                              'analysis_ready': man.analysis_ready}
    return {'out_dir': str(out_dir), 'frames': results,
            'n_frames': len(frames)}


# ── CLI / smoke test ──────────────────────────────────────────────────────────────
def _availability_report() -> dict:
    """Engine + per-night-GDAS availability snapshot (the acceptance install check)."""
    rep = {}
    for inst, eng in sorted(C.TELLURIC_ENGINES.items()):
        rep[inst] = {'engine': eng, 'installed': engine_available(eng),
                     'site': C.TELLURIC_INSTRUMENT_SITE.get(inst)}
    return rep


def main(argv=None):
    import argparse
    ap = argparse.ArgumentParser(description='RYA-424 telluric data-input stage')
    ap.add_argument('--availability', action='store_true',
                    help='report engine + GDAS availability by instrument')
    ap.add_argument('--vesta-crires', action='store_true',
                    help='run the stage end-to-end on the Vesta CRIRES+ CO frames')
    args = ap.parse_args(argv)
    print('=' * 84)
    print('  RYA-424 — telluric correction as a standing data-input stage')
    print('=' * 84)
    rep = _availability_report()
    print('\n  ENGINE ROUTING (instrument-keyed, single source of truth):')
    print(f"    {'instrument':<14}{'engine':<16}{'site':<12}{'installed'}")
    for inst, r in rep.items():
        print(f"    {inst:<14}{r['engine']:<16}{str(r['site']):<12}"
              f"{'yes' if r['installed'] else 'NOT INSTALLED'}")
    gdas_dir = _telluriccorr_share()
    print(f"\n  per-night GDAS install: {gdas_dir if gdas_dir else 'NOT FOUND'}")
    if args.vesta_crires:
        print('\n  RUNNING the stage on the Vesta CRIRES+ CO frames …')
        res = condition_vesta_crires()
        print(json.dumps(res, indent=2, default=str))
    print('=' * 84)
    return rep


if __name__ == '__main__':
    main()
