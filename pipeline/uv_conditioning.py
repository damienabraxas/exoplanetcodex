#!/usr/bin/env python3
"""
pipeline/uv_conditioning.py — UV (FUV/NUV) conditioning as a standing data-input stage
(RYA-426). The no-telluric sibling of the RYA-424 telluric stage: space UV sits above the
atmosphere, so there is NO telluric step, but it has its own hard problems that silently
corrupt abundances if unhandled. This stage conditions a UV spectrum to the same standard
as RYA-424 — explicit gates, an analysis_ready flag, LOUD failure over silent fallback —
and writes the result to a conditioning manifest (the shared schema RYA-424 aligns to).

SCOPE BOUNDARY (deliberate): this is the CONDITIONING stage (is the spectrum analysis-ready)
— it does NOT rule which lines are scientifically usable nor apply NLTE values; that is
RYA-190 (line selection + NLTE policy). The diagnostic anchors here are CITED coverage
positions for the window map only, never a usability ruling.

THE GATES (RYA-426 §3):
  1. Star identity        — alpha Cen only (A/B header-mixed); single stars (Sun/Procyon) skip.
  2. Vacuum -> air         — at the loader boundary. UV is vacuum; pipeline is air. BUT air is
                             only defined for lambda >= 2000 A (IAU/VALD): the FUV (< 2000 A)
                             STAYS vacuum, and only NUV (>= 2000 A) is converted. Converting a
                             FUV line to "air" is meaningless and lands it at the wrong place.
  3. Scattered light       — STIS echelle inter-order scattered light is significant in the FUV.
                             VERIFY the pipeline background handling (do not assume; do not re-subtract).
  4. Chromospheric masking — Mg II h&k etc. have emission-filled cores that are NOT photospheric.
                             Mask the cores; a filled core fed to an abundance is a silent overestimate.
  5. FUV synthesis-not-EW  — the FUV has no true continuum (line haze). EW is INVALID there and is
                             REFUSED; abundances must come from line-by-line synthesis.
  6. ISM check             — resonance lines can carry narrow interstellar components; flag by distance.
  7. NLTE/3D coverage flag — UV lines are often strongly NLTE-sensitive and our NLTE grids do not
                             extend to UV transitions. LTE-flag LOUDLY where no UV grid exists.
  8. Window-to-diagnostic  — map which STIS/COS windows cover which diagnostics; flag
                             lines-without-spectra and spectra-without-lines.

analysis_ready is True only after gates 2,3,4,5,7 pass (single star skips 1; 6/8 are advisory flags).

Permanent rules (RYA-426): loud failure over silent fallback; single source of truth
(grating/LSF/window config cited below); never coadd across resolutions; synthesis mandatory
in the FUV; Angstrom (air) for lambda >= 2000 A, vacuum below; vacuum conversion at the loader.
"""
from __future__ import annotations

import json
import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

import config.constants as const
from pipeline.wavelength_util import AIR_VACUUM_BOUNDARY_A

ROOT = Path(str(const.ROOT))
CONDITIONING_DIR = ROOT / 'data' / 'conditioning'

# ── air/vacuum + FUV/NUV regime boundary ─────────────────────────────────────
# Air refractive index is only defined for lambda >= 2000 A (IAU 1991; Morton 2000;
# the VALD air/vacuum convention used by linelist_solar.csv). RYA-303 recorded the same:
# STIS arrays are VACUUM; the pipeline is air for lambda >= 2000 A. So 2000 A is BOTH the
# air/vacuum boundary AND the STIS FUV-MAMA / NUV-MAMA detector split — one number.
# The constant + the vac<->air converter now live in the shared wavelength_util SSOT
# (RYA-264); this is the single boundary value reused here for the regime split.
UV_AIR_VACUUM_BOUNDARY_A = AIR_VACUUM_BOUNDARY_A

# ── STIS / COS UV gratings (single source of truth; STIS IHB + RYA-384/303/426) ──
# 'frame' = native wavelength frame of the science product (always vacuum for HST UV).
# 'cov_A' is nominal full coverage; per-exposure windows come from the RYA-222 audit.
STIS_UV_GRATINGS = {
    'E140H': {'regime': 'FUV', 'detector': 'FUV-MAMA', 'R': 114000, 'cov_A': (1140.0, 1700.0)},
    'E140M': {'regime': 'FUV', 'detector': 'FUV-MAMA', 'R': 45800,  'cov_A': (1140.0, 1730.0)},
    'E230H': {'regime': 'NUV', 'detector': 'NUV-MAMA', 'R': 114000, 'cov_A': (1620.0, 3150.0)},
    'E230M': {'regime': 'NUV', 'detector': 'NUV-MAMA', 'R': 30000,  'cov_A': (1570.0, 3110.0)},
    'G230MB': {'regime': 'NUV', 'detector': 'CCD',     'R': 10000,  'cov_A': (1990.0, 3160.0)},
    # COS — FUV only (RYA-426 §2)
    'G130M': {'regime': 'FUV', 'detector': 'COS-FUV', 'R': 16000, 'cov_A': (1150.0, 1450.0)},
    'G160M': {'regime': 'FUV', 'detector': 'COS-FUV', 'R': 16000, 'cov_A': (1405.0, 1775.0)},
}
EXCLUDED_GRATINGS = {'WFC3'}      # permanently excluded (RYA-426 §2)

# ── chromospheric emission cores to MASK (cited; cores are NOT photospheric) ──
# lambda is the line centre in the stated frame; a window of +-mask_hw_A around it is masked.
UV_CHROMOSPHERIC_LINES = (
    {'name': 'Mg II k', 'lambda_A': 2796.3543, 'frame': 'vacuum', 'mask_hw_A': 0.7, 'ref': 'Morton 2003 ApJS 149,205'},
    {'name': 'Mg II h', 'lambda_A': 2803.5324, 'frame': 'vacuum', 'mask_hw_A': 0.7, 'ref': 'Morton 2003 ApJS 149,205'},
    {'name': 'C II 1334', 'lambda_A': 1334.532, 'frame': 'vacuum', 'mask_hw_A': 0.3, 'ref': 'Morton 2003'},
    {'name': 'C II* 1335', 'lambda_A': 1335.708, 'frame': 'vacuum', 'mask_hw_A': 0.3, 'ref': 'Morton 2003'},
    {'name': 'Si IV 1393', 'lambda_A': 1393.760, 'frame': 'vacuum', 'mask_hw_A': 0.3, 'ref': 'Morton 2003'},
    {'name': 'Si IV 1402', 'lambda_A': 1402.773, 'frame': 'vacuum', 'mask_hw_A': 0.3, 'ref': 'Morton 2003'},
    {'name': 'C IV 1548', 'lambda_A': 1548.187, 'frame': 'vacuum', 'mask_hw_A': 0.3, 'ref': 'Morton 2003'},
    {'name': 'C IV 1550', 'lambda_A': 1550.772, 'frame': 'vacuum', 'mask_hw_A': 0.3, 'ref': 'Morton 2003'},
    {'name': 'He II 1640', 'lambda_A': 1640.42, 'frame': 'vacuum', 'mask_hw_A': 0.4, 'ref': 'Morton 2003'},
    # NUV (>= 2000 A) chromospheric cores are quoted in AIR (their pipeline frame)
    {'name': 'Ca II K', 'lambda_A': 3933.66, 'frame': 'air', 'mask_hw_A': 0.8, 'ref': 'NIST ASD'},
    {'name': 'Ca II H', 'lambda_A': 3968.47, 'frame': 'air', 'mask_hw_A': 0.8, 'ref': 'NIST ASD'},
)

# ── diagnostic ANCHORS for the window map (CITED coverage positions ONLY) ──
# NOT a usability ruling and NOT an NLTE policy — RYA-190 owns line selection + NLTE.
# 'usable_ref' cites RYA-190's published verdict for context; this stage never decides it.
UV_DIAGNOSTIC_ANCHORS = (
    {'species': 'C I', 'element': 'C', 'lambda_A': 1657.38, 'frame': 'vacuum',
     'usable_ref': 'RYA-190: usable, NLTE ~+0.10 (Amarsi 2020); synthesis-not-EW (FUV)'},
    {'species': 'O I', 'element': 'O', 'lambda_A': 1355.60, 'frame': 'vacuum',
     'usable_ref': 'RYA-190: usable (semi-forbidden); NLTE-corrected'},
    {'species': 'O I', 'element': 'O', 'lambda_A': 1302.17, 'frame': 'vacuum',
     'usable_ref': 'RYA-190: DO NOT USE (resonance triplet saturated)'},
    {'species': 'N I', 'element': 'N', 'lambda_A': 1199.55, 'frame': 'vacuum',
     'usable_ref': 'RYA-190: DO NOT USE (resonance, NLTE 0.3-0.5 dex)'},
    {'species': 'S I', 'element': 'S', 'lambda_A': 1473.99, 'frame': 'vacuum',
     'usable_ref': 'RYA-190: high value (UV S access)'},
    {'species': 'C I', 'element': 'C', 'lambda_A': 1930.90, 'frame': 'vacuum',
     'usable_ref': 'RYA-190: C I UV (secondary multiplet)'},
)

# UV NLTE grids present in the repo (none yet — RYA-190/165 owe them). Used by gate 7 to
# LTE-flag loudly. This is a registry, not a ruling: it lists grids that EXIST.
_UV_NLTE_GRIDS_PRESENT: dict[str, str] = {}     # element -> grid file; empty == none on disk


class FUVSynthesisRequired(RuntimeError):
    """Raised when an EW measurement is attempted in the FUV (gate 5 — EW is invalid there)."""


class UVConditioningError(RuntimeError):
    """A UV conditioning gate failed hard (loud failure over silent fallback)."""


# ── gate 2: vacuum <-> air — the shared SSOT converter (RYA-264/426) ─────────
# Birch & Downs 1994 (VALD3/NIST). The implementation lives once in wavelength_util;
# re-exported here so existing `uvc.vac_to_air` / `uvc.air_to_vac` importers are
# unchanged. Below 2000 A both are the identity → FUV stays vacuum, NUV converts.
from pipeline.wavelength_util import vac_to_air, air_to_vac   # noqa: E402,F401


def to_pipeline_frame(wave_vac_A) -> np.ndarray:
    """Apply the loader-boundary convention to a native-vacuum UV grid: NUV -> air,
    FUV stays vacuum. This is what gate 2 hands downstream."""
    return vac_to_air(wave_vac_A)


# ── regime + gate 5 (FUV synthesis-not-EW) ───────────────────────────────────
def classify_regime(lambda_A: float) -> str:
    return 'FUV' if float(lambda_A) < UV_AIR_VACUUM_BOUNDARY_A else 'NUV'


def is_fuv(lambda_A) -> np.ndarray:
    return np.asarray(lambda_A, float) < UV_AIR_VACUUM_BOUNDARY_A


def refuse_ew_in_fuv(lambda_A) -> None:
    """Gate 5: EW is invalid in the FUV (no true continuum). Refuse it LOUDLY."""
    arr = np.atleast_1d(np.asarray(lambda_A, float))
    if np.any(arr < UV_AIR_VACUUM_BOUNDARY_A):
        bad = arr[arr < UV_AIR_VACUUM_BOUNDARY_A]
        raise FUVSynthesisRequired(
            f"EW measurement refused for FUV lambda {np.round(bad[:3], 2).tolist()} "
            f"(< {UV_AIR_VACUUM_BOUNDARY_A} A): the FUV has no true continuum — use spectral "
            f"synthesis against a pseudo-continuum (RYA-426 gate 5).")


# ── gate 3: scattered-light / background verification ────────────────────────
def scattered_light_check(flux, err=None, neg_frac_max: float = 0.05) -> dict:
    """VERIFY (not re-subtract) the pipeline background handling. HASP/x1d cspec products
    are already background-subtracted; over-subtracted FUV scattered light shows up as a
    population of negative-flux pixels. Flag the negative fraction; require it small for
    analysis_ready. Loud, explicit — never silently 'fix' it here."""
    f = np.asarray(flux, float)
    finite = np.isfinite(f)
    n = int(finite.sum())
    neg_frac = float((f[finite] < 0).sum() / n) if n else 1.0
    ok = neg_frac <= neg_frac_max
    return {'gate': 'scattered_light', 'method': 'verify cspec/x1d background (no re-subtraction)',
            'neg_flux_fraction': round(neg_frac, 4), 'neg_frac_max': neg_frac_max,
            'passed': bool(ok),
            'note': ('background handling verified (low negative-flux fraction)' if ok else
                     'HIGH negative-flux fraction — scattered-light/over-subtraction REVIEW required')}


# ── gate 4: chromospheric-core masking ───────────────────────────────────────
def chromospheric_core_mask(wave_pipeline_A) -> np.ndarray:
    """Boolean mask (True == masked emission core) over a pipeline-frame grid. Lines are
    compared in the SAME frame as the grid: NUV cores (>=2000) are converted vac->air to
    match, FUV cores stay vacuum. A filled core fed to an abundance is a silent overestimate."""
    w = np.asarray(wave_pipeline_A, float)
    mask = np.zeros(w.shape, dtype=bool)
    for line in UV_CHROMOSPHERIC_LINES:
        centre = line['lambda_A']
        # bring the line centre into the pipeline frame (NUV cores given in vacuum -> air)
        if line['frame'] == 'vacuum':
            centre = float(to_pipeline_frame(np.array([centre]))[0])
        mask |= np.abs(w - centre) <= line['mask_hw_A']
    return mask


def chromospheric_lines_in_range(lo_A: float, hi_A: float) -> list:
    out = []
    for line in UV_CHROMOSPHERIC_LINES:
        c = line['lambda_A']
        if line['frame'] == 'vacuum':
            c = float(to_pipeline_frame(np.array([c]))[0])
        if lo_A <= c <= hi_A:
            out.append({'name': line['name'], 'lambda_pipeline_A': round(c, 3),
                        'mask_hw_A': line['mask_hw_A'], 'ref': line['ref']})
    return out


# ── gate 6: ISM review flag ──────────────────────────────────────────────────
def ism_review_flag(distance_pc: float, anchors=UV_DIAGNOSTIC_ANCHORS,
                    distance_threshold_pc: float = 50.0) -> dict:
    """Resonance UV lines can carry narrow interstellar components. Negligible for nearby
    targets; flag all UV resonance anchors for ISM review beyond the threshold (RYA-190)."""
    flag = float(distance_pc) > distance_threshold_pc
    return {'gate': 'ism_review', 'distance_pc': float(distance_pc),
            'threshold_pc': distance_threshold_pc, 'review_required': bool(flag),
            'note': (f'distance > {distance_threshold_pc} pc — flag UV resonance lines for ISM review'
                     if flag else 'nearby target — ISM column negligible')}


# ── gate 7: NLTE/3D coverage flag ────────────────────────────────────────────
def nlte_coverage_flags(anchors=UV_DIAGNOSTIC_ANCHORS) -> list:
    """Per anchor: does a UV NLTE grid exist? None do yet (RYA-190/165 owe them) — so every
    UV diagnostic is LTE-flagged LOUDLY. This FLAGS coverage; it never applies a value
    (that is RYA-190's policy)."""
    out = []
    for a in anchors:
        has = a['element'] in _UV_NLTE_GRIDS_PRESENT
        out.append({'species': a['species'], 'lambda_A': a['lambda_A'],
                    'uv_nlte_grid': _UV_NLTE_GRIDS_PRESENT.get(a['element'], 'NONE'),
                    'lte_flag': 'LTE_ASSUMED_LOUD' if not has else 'NLTE_GRID_AVAILABLE',
                    'note': ('NO UV NLTE grid on disk — LTE-flagged loudly (RYA-190/165 owed)'
                             if not has else 'UV NLTE grid present')})
    return out


# ── gate 8: window-to-diagnostic map ─────────────────────────────────────────
def window_to_diagnostic_map(windows, anchors=UV_DIAGNOSTIC_ANCHORS) -> dict:
    """windows: iterable of dicts {grating, lo_A, hi_A} (pipeline-frame coverage). Map which
    windows cover which anchors; flag lines-without-spectra and spectra-without-lines."""
    windows = list(windows)
    covered, lines_without_spectra = [], []
    for a in anchors:
        # anchors are given in vacuum; bring to pipeline frame for the comparison
        c = float(to_pipeline_frame(np.array([a['lambda_A']]))[0])
        hits = [w['grating'] for w in windows if w['lo_A'] <= c <= w['hi_A']]
        rec = {'species': a['species'], 'lambda_A': a['lambda_A'],
               'lambda_pipeline_A': round(c, 3), 'covered_by': hits}
        (covered if hits else lines_without_spectra).append(rec)
    anchor_centres = [float(to_pipeline_frame(np.array([a['lambda_A']]))[0]) for a in anchors]
    spectra_without_lines = [w['grating'] for w in windows
                             if not any(w['lo_A'] <= c <= w['hi_A'] for c in anchor_centres)]
    return {'gate': 'window_to_diagnostic', 'n_windows': len(windows),
            'covered': covered, 'lines_without_spectra': lines_without_spectra,
            'spectra_without_lines': sorted(set(spectra_without_lines))}


# ── the conditioning manifest (shared schema; RYA-424 aligns to this) ────────
_MANIFEST_VERSION = 'uv_conditioning/v1'


@dataclass
class ConditioningResult:
    star: str
    grating: str
    regime: str                       # FUV | NUV
    instrument: str
    wave_range_A: tuple
    n_pixels: int
    analysis_ready: bool
    gates: dict = field(default_factory=dict)
    masked_fraction: float = 0.0
    telluric_applied: bool = False    # ALWAYS False for space UV (recorded explicitly)
    frame_out: str = 'air>=2000/vac<2000'
    schema: str = _MANIFEST_VERSION

    def to_manifest(self) -> dict:
        d = self.__dict__.copy()
        d['wave_range_A'] = [round(float(self.wave_range_A[0]), 3), round(float(self.wave_range_A[1]), 3)]
        return d


def conditioning_dir(star: str, *, create: bool = True) -> Path:
    d = CONDITIONING_DIR / str(star).strip().lower().replace(' ', '_')
    if create:
        d.mkdir(parents=True, exist_ok=True)
    return d


def condition_uv_spectrum(wave_vac_A, flux, err=None, *, star: str, grating: str,
                          instrument: str = 'STIS', distance_pc: float = 0.0,
                          is_binary: bool = False, identity_confirmed: bool = True,
                          windows=None, write: bool = True) -> ConditioningResult:
    """Run the full UV conditioning chain on one (native-vacuum) UV spectrum and return the
    ConditioningResult (+ write the manifest). Loud failure over silent fallback."""
    grating = str(grating).upper()
    if grating in EXCLUDED_GRATINGS:
        raise UVConditioningError(f"{grating} is permanently excluded (RYA-426 §2)")
    if grating not in STIS_UV_GRATINGS:
        raise UVConditioningError(f"unknown UV grating {grating!r} — add it to STIS_UV_GRATINGS "
                                  f"(single source of truth) before conditioning")
    gcfg = STIS_UV_GRATINGS[grating]
    wave_vac_A = np.asarray(wave_vac_A, float)
    flux = np.asarray(flux, float)
    gates: dict = {}

    # gate 1 — star identity (binary only)
    if is_binary and not identity_confirmed:
        raise UVConditioningError(
            f"{star}: binary UV identity INDETERMINATE — never route STIS frames by header "
            f"for an A/B pair (RYA-426 gate 1). Resolve via program/epoch cross-match first.")
    gates['identity'] = {'gate': 'identity', 'is_binary': bool(is_binary),
                         'confirmed': bool(identity_confirmed),
                         'note': 'single star — gate skipped' if not is_binary else 'confirmed'}

    # gate 2 — vacuum -> air (NUV) / stay vacuum (FUV)
    wave_out = to_pipeline_frame(wave_vac_A)
    n_converted = int(np.count_nonzero(wave_vac_A >= UV_AIR_VACUUM_BOUNDARY_A))
    gates['vacuum_to_air'] = {'gate': 'vacuum_to_air', 'standard': 'Birch & Downs 1994 (VALD/NIST)',
                              'boundary_A': UV_AIR_VACUUM_BOUNDARY_A, 'n_nuv_converted': n_converted,
                              'n_fuv_kept_vacuum': int(wave_out.size - n_converted), 'passed': True}

    # gate 3 — scattered light / background
    gates['scattered_light'] = scattered_light_check(flux, err)

    # gate 4 — chromospheric core masking
    core_mask = chromospheric_core_mask(wave_out)
    masked_frac = float(core_mask.sum() / core_mask.size) if core_mask.size else 0.0
    gates['chromospheric_mask'] = {
        'gate': 'chromospheric_mask', 'rule': 'mask +-hw around emission cores; cores are '
        'chromospheric (not photospheric) -> unusable for abundance',
        'lines_masked_in_range': chromospheric_lines_in_range(float(wave_out.min()), float(wave_out.max())),
        'masked_pixel_fraction': round(masked_frac, 4), 'passed': True}

    # gate 5 — FUV synthesis-not-EW
    regime = gcfg['regime']
    fuv_present = bool(np.any(wave_out < UV_AIR_VACUUM_BOUNDARY_A))
    gates['synthesis_gate'] = {'gate': 'synthesis_not_ew', 'regime': regime,
                               'fuv_pixels_present': fuv_present,
                               'ew_allowed': (not fuv_present),
                               'rule': 'FUV (<2000 A) abundances MUST use synthesis; EW refused',
                               'passed': True}

    # gate 6 — ISM
    gates['ism_review'] = ism_review_flag(distance_pc)

    # gate 7 — NLTE coverage
    nl = nlte_coverage_flags()
    gates['nlte_coverage'] = {'gate': 'nlte_coverage', 'per_anchor': nl,
                              'all_lte_flagged_loud': all(x['lte_flag'] == 'LTE_ASSUMED_LOUD' for x in nl),
                              'passed': True}

    # gate 8 — window-to-diagnostic (uses provided windows, else this spectrum's own span)
    if windows is None:
        windows = [{'grating': grating, 'lo_A': float(wave_out.min()), 'hi_A': float(wave_out.max())}]
    gates['window_map'] = window_to_diagnostic_map(windows)

    # analysis_ready: gates 2,3,4,5,7 must be clean (1 handled above; 6,8 advisory)
    analysis_ready = bool(gates['vacuum_to_air']['passed']
                          and gates['scattered_light']['passed']
                          and gates['chromospheric_mask']['passed']
                          and gates['synthesis_gate']['passed']
                          and gates['nlte_coverage']['passed'])

    result = ConditioningResult(
        star=star, grating=grating, regime=regime, instrument=instrument,
        wave_range_A=(float(wave_out.min()), float(wave_out.max())), n_pixels=int(wave_out.size),
        analysis_ready=analysis_ready, gates=gates, masked_fraction=round(masked_frac, 4),
        telluric_applied=False)

    if write:
        out = conditioning_dir(star) / f'{star.strip().lower().replace(" ", "_")}_{grating}_uv_conditioning.json'
        out.write_text(json.dumps(result.to_manifest(), indent=2))
    return result


# ── smoke harness (NOT the production loader — that is RYA-471) ───────────────
def _read_stis_cspec(path: str):
    """Thin smoke-test reader for one HASP/STIS _cspec.fits (vacuum WAVELENGTH/FLUX/ERROR).
    The production loader is RYA-471; this is only enough to exercise the stage on real data."""
    from astropy.io import fits
    with fits.open(path) as h:
        sci = h['SCI'].data if 'SCI' in [e.name for e in h] else h[1].data
        w = np.asarray(sci['WAVELENGTH']).ravel().astype(float)
        f = np.asarray(sci['FLUX']).ravel().astype(float)
        e = np.asarray(sci['ERROR']).ravel().astype(float) if 'ERROR' in sci.columns.names else None
        targ = h[0].header.get('TARGNAME', '')
    order = np.argsort(w)
    return w[order], f[order], (e[order] if e is not None else None), targ


def _smoke(star='procyon'):
    """End-to-end on real Procyon STIS (RYA-222 science-ready whitelist): vacuum->air ->
    scattered-light -> core-mask -> synthesis-ready -> flagged. Verifies vac->air on Mg II."""
    import csv
    # vacuum->air verification on a KNOWN line: Mg II k 2796.3543 vac -> 2795.528 air (Morton/NIST)
    air = float(vac_to_air(np.array([2796.3543, 2803.5324]))[0])
    air2 = float(vac_to_air(np.array([2803.5324]))[0])
    ok_h = abs(air - 2795.528) < 0.01
    ok_k = abs(air2 - 2802.705) < 0.01
    print(f"[vacuum->air] Mg II k 2796.3543(vac) -> {air:.3f}(air)  [known 2795.528]  {'OK' if ok_h else 'FAIL'}")
    print(f"[vacuum->air] Mg II h 2803.5324(vac) -> {air2:.3f}(air) [known 2802.705]  {'OK' if ok_k else 'FAIL'}")
    print(f"[boundary]    FUV C I 1657.38 stays vacuum: {float(vac_to_air(np.array([1657.38]))[0]):.3f} "
          f"(unchanged below {UV_AIR_VACUUM_BOUNDARY_A} A)")

    sr = ROOT / 'data' / 'audit' / 'procyon_hst' / 'procyon_hst_science_ready.csv'
    if not sr.exists():
        print(f"[smoke] science-ready whitelist absent ({sr}) — vacuum->air verified above; data step skipped")
        return ok_h and ok_k
    # RYA-222 whitelist ONLY (the 'Procyon HST' tree also holds 55 Cnc frames — never glob)
    rows = [r for r in csv.DictReader(open(sr))
            if str(r.get('is_science_grating', '')).strip() in ('True', 'true', '1')
            and Path(r['filepath']).exists()]
    windows = [{'grating': r['opt_elem'], 'lo_A': float(vac_to_air(np.array([float(r['wl_min_A'])]))[0]),
                'hi_A': float(vac_to_air(np.array([float(r['wl_max_A'])]))[0])} for r in rows]
    for gr in ('E230H', 'E140H'):                          # one NUV + one FUV
        cand = [r for r in rows if r['opt_elem'] == gr]
        if not cand:
            continue
        w, f, e, targ = _read_stis_cspec(cand[0]['filepath'])
        res = condition_uv_spectrum(w, f, e, star=star, grating=gr, instrument='STIS',
                                    distance_pc=3.51, windows=windows, write=True)
        print(f"\n[{gr}] {Path(cand[0]['filepath']).name}  targ={targ}")
        print(f"  regime={res.regime}  range={res.wave_range_A} A  npix={res.n_pixels}  "
              f"telluric_applied={res.telluric_applied}")
        print(f"  scattered-light neg-frac={res.gates['scattered_light']['neg_flux_fraction']} "
              f"pass={res.gates['scattered_light']['passed']}")
        print(f"  chromospheric masked-frac={res.masked_fraction} "
              f"(cores: {[l['name'] for l in res.gates['chromospheric_mask']['lines_masked_in_range']]})")
        print(f"  synthesis-not-EW: ew_allowed={res.gates['synthesis_gate']['ew_allowed']} "
              f"(FUV pixels present={res.gates['synthesis_gate']['fuv_pixels_present']})")
        print(f"  ANALYSIS_READY={res.analysis_ready}")
    wm = window_to_diagnostic_map(windows)
    print(f"\n[window->diagnostic] covered={[c['species']+' '+str(c['lambda_A']) for c in wm['covered']]}")
    print(f"  lines-without-spectra={[c['species']+' '+str(c['lambda_A']) for c in wm['lines_without_spectra']]}")
    return ok_h and ok_k


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser(description='RYA-426 UV conditioning stage')
    ap.add_argument('--smoke', action='store_true', help='end-to-end on real Procyon STIS')
    ap.add_argument('--star', default='procyon')
    args = ap.parse_args()
    if args.smoke:
        ok = _smoke(args.star)
        sys.exit(0 if ok else 1)
    ap.print_help()
