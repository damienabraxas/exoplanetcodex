"""
pipeline/reflected_solar_rv.py
==============================
Asteroid-ephemeris reflected-solar velocity-frame conditioning (RYA-372).

The RYA-370 audit cleared the Vesta reflected-solar **optical** set for science
but found it is **not run-ready as downloaded**: the standard velocity frames do
NOT land solar lines at rest (ESPRESSO BARYCENT residual −4.6…+33 km/s; UVES
TOPOCENT −17…+31 km/s). The residual is the **asteroid reflected radial velocity**
(Sun → Vesta → Earth), which the *stellar* BERV keyword cannot remove. This module
is the single source of truth for conditioning each frame to the solar rest frame —
the Phase-A blocker for RYA-371 (the optical CNO arms cannot synthesize against
velocity-shifted spectra). Arm-agnostic; reused by the CRIRES+ IR conditioning
(RYA-373). Upstream data conditioning ONLY — does not touch the EW/abundance modules.

Method (Molaro et al. 2013 / Lanza reflected-solar treatment)
------------------------------------------------------------
The observed solar-line shift is the sum of two legs:
  leg 1  Sun → Vesta : heliocentric range rate   (Horizons `r_rate`)
  leg 2  Vesta → observer : topocentric range rate (Horizons `delta_rate`)
`reflected_solar_rv()` returns this two-leg ephemeris model.

**But the sum + sign is asserted empirically, never trusted** — that is the entire
point of this ticket. Two facts make the raw two-leg sum insufficient on its own
(verified on these data, see RYA-372 comment):
  1. The IDP wavelength frames are NOT header-honest (the RYA-271 rule): an
     ESPRESSO BARYCENT frame has had a BERV applied, a UVES "TOPOCENT" frame
     behaves as if a barycentric-scale shift was applied — so the residual that
     must be removed is frame-dependent and the header cannot be trusted.
  2. Leg 2 in a barycentric frame needs Vesta's barycentric velocity projected on
     the *Earth→Vesta* line of sight, which differs from `delta_rate(@SSB)` by the
     parallactic angle (~7–9 km/s) — exactly where reflected-solar work goes wrong.

So the conditioning is **anchored to the measured photospheric bulk velocity**
(robust median over many clean Fe I cores; lands lines at rest by construction),
with the Horizons two-leg as the independent physical cross-check, and the result
**verified on a held-out line set** (CRITICAL loud-fail if it does not land at rest).
Per-frame only — Vesta's 5.3 h rotation makes the velocity epoch-variable; never
coadd before shifting. ESPRESSO: S1D_FINAL_A only (STACK_A is velocity-smeared —
rejected with a logged reason), `WAVE_AIR` (air); UVES: per-dichroic, air `WAVE`.

Smoke test
----------
  python -m pipeline.reflected_solar_rv --set vesta_espresso --verify
  python -m pipeline.reflected_solar_rv --set vesta_uves --verify
"""
from __future__ import annotations

import argparse
import sys
import warnings
from pathlib import Path

import numpy as np
from astropy.io import fits

from config.constants import REFLECTED_SOLAR_BODIES, codex_path   # RYA-394 single-source body registry

warnings.simplefilter('ignore')

C_KMS = 299792.458
# Default reflected-solar body key (the Vesta optical/IR set). RYA-394: this is a
# REGISTRY KEY, never a raw Horizons id — a bare id '4' resolves to Mars barycenter.
DEFAULT_BODY = 'vesta'


class BodyIDError(RuntimeError):
    """Horizons resolved a body to something other than the intended target (e.g. the
    bare-'4'⇒Mars-barycenter trap), or an unknown body key was passed. Never silently
    return the wrong body's velocity (RYA-394)."""


class RestFrameError(RuntimeError):
    """A frame's lines do NOT land at rest after RV conditioning — the applied velocity
    (hence the body ID or sign) is wrong. The closed-loop assert RYA-372 specced but left
    report-only, which let the Mars-for-Vesta swap ride four tickets (RYA-394). Carries
    the partial conditioning record on `.rec` for survey tools."""

    def __init__(self, message, rec=None):
        super().__init__(message)
        self.rec = rec
PARANAL = '309'            # JPL Horizons observatory code for ESO Paranal (VLT)

# Reflected-solar data root (outside the repo — RYA-370). Override with --root.
DEFAULT_ROOT = (codex_path('data.spectra_local')
                / 'Solar Calibration' / 'Solar System Targets' / 'Vesta')

# Chromospheric reference lines (deep cores form high → carry a small velocity
# offset vs the photosphere; reported, not used to anchor).
NAD2_AIR = 5889.951
HALPHA_AIR = 6562.797

# Disk-integrated solar gravitational redshift (Sun→observer), km/s — a known
# constant offset present in every solar spectrum; reported for context.
SOLAR_GRAVITATIONAL_REDSHIFT = 0.633

# Clean, relatively unblended photospheric lines (mostly Fe I), air Å, spanning the
# optical so at least one dichroic-arm subset is always in range. Split TRAIN/TEST
# by index parity so verification is on lines NOT used to derive the correction.
VEL_LINES_AIR = (
    3719.935, 3737.132, 3820.425, 3859.911, 4063.594, 4271.760, 4383.545,
    4404.750, 4459.118, 4528.614, 4602.941, 4647.434, 4736.773, 4871.318,
    4890.755, 4920.502, 4957.596, 5006.119, 5051.634, 5079.740, 5151.911,
    5191.455, 5232.940, 5247.050, 5269.537, 5328.039, 5371.489, 5397.128,
    5405.775, 5429.696, 5434.524, 5446.917, 5497.516, 5506.779, 5569.618,
    5576.089, 5586.756, 6065.482, 6137.692, 6173.334, 6191.558, 6213.430,
    6230.723, 6232.641, 6252.555, 6265.134, 6322.685, 6335.330, 6393.601,
    6411.649, 6421.350, 6430.846, 6494.980, 6592.914, 6677.987,
)
TRAIN_LINES = VEL_LINES_AIR[0::2]
TEST_LINES = VEL_LINES_AIR[1::2]

# Strong, deep, isolated lines used ONLY to seed the coarse bulk-shift estimate on
# low-SNR frames where the photospheric cores alone don't lock on (a wide search
# finds them reliably; their small chromospheric offset is negligible at the coarse
# stage and is refined away by the photospheric fine pass).
STRONG_ANCHORS_AIR = (4226.728, 4861.350, 5183.604, 5889.951, 5895.924, 6162.170, 6562.797)

PASS_TOL_KMS = 0.5        # held-out photospheric lines must land within this of rest
MIN_LINES = 5             # minimum clean lines for a valid bulk-velocity measurement

# Telluric-CCF wavelength-zero-point closure (RYA-372 review item 2). The held-out
# photospheric check CANNOT catch a uniform wavelength zero-point offset: a common-mode
# shift δ moves the anchor lines AND the held-out lines together, so zeroing the anchor
# lands everything at rest by construction. The telluric lines are the one header-
# independent reference — they sit at rest in the TOPOCENTRIC frame regardless of what
# the pipeline did. So we measure the telluric velocity (iSpec CCF) and assert it equals
# the pipeline's *claimed* applied frame shift: BERV for a BARYCENT frame, 0 for a
# TOPOCENT frame. A mismatch = a common-mode wavelength offset (or a frame mislabel) →
# loud-fail. Observed closure (RYA-372): ESPRESSO v_tell−BERV = +0.13 ± 0.03 km/s;
# UVES v_tell−0 = −0.04 ± 0.52 km/s (→ UVES is genuinely TOPOCENT, telluric at rest).
TELLURIC_CCF_MASK = ('input', 'linelists', 'CCF', 'Synth.Tellurics.500_1100nm', 'mask.lst')
TELL_CLOSURE_TOL = 1.5    # start (review §2); tighten per-instrument once distribution seen
BERV_KEYWORD = 'HIERARCH ESO QC BERV'   # confirmed against ESPRESSO headers (review note)

# "rest = mean Fe I core" carries the solar convective-blueshift zero-point
# (~−0.3 to −0.4 km/s, line-depth-dependent) + the +0.633 km/s gravitational redshift.
# Harmless for EW/synthesis (shift-invariant) and it cancels differentially, but it is
# NOT the laboratory rest frame — stamped in the --write provenance so a conditioned
# frame is not mistaken for it.
CONVECTIVE_BLUESHIFT_NOTE = (
    'rest frame = mean photospheric Fe I core (NOT laboratory rest): carries the solar '
    'convective blueshift (~-0.3..-0.4 km/s, depth-dependent) + gravitational redshift '
    '(+0.633 km/s). Shift-invariant for EW/synthesis; cancels in differential abundances.'
)


# ── Horizons two-leg ephemeris model ──────────────────────────────────────────
_HZ_CACHE: dict = {}


def reflected_solar_rv(mjd: float, body_key: str = DEFAULT_BODY,
                       obs_code: str = PARANAL) -> dict:
    """Molaro/Lanza two-leg reflected-solar velocity model from JPL Horizons, km/s.

    leg1 Sun→body  = `r_rate` (heliocentric range rate, v_helio);
    leg2 body→observer = `delta_rate` (topocentric range rate, v_obs);
    v_total = r_rate + delta_rate (the leading-order topocentric reflected-solar
    line-of-sight velocity; sign convention POSITIVE = receding / redshift).

    `body_key` indexes the single-source REFLECTED_SOLAR_BODIES registry (RYA-394) —
    raw Horizons ids are NO LONGER accepted, so the bare-'4'⇒Mars-barycenter trap cannot
    recur. The resolved Horizons targetname is asserted against the registry `match`
    (RYA-374 discipline): a body resolving to a non-matching target raises BodyIDError,
    never silently returns the wrong body's velocity. Returns the v_helio/v_obs/v_total
    breakdown plus the resolved `targetname`.

    This is the PHYSICAL MODEL / cross-check. The empirical photospheric anchor in
    `condition_frame` is what conditions to rest; this closes the loop against it.
    """
    from astroquery.jplhorizons import Horizons
    try:
        spec = REFLECTED_SOLAR_BODIES[body_key]
    except (KeyError, TypeError):
        raise BodyIDError(
            f"Unknown reflected-solar body key {body_key!r}. Pass a key from "
            f"REFLECTED_SOLAR_BODIES ({sorted(REFLECTED_SOLAR_BODIES)}), never a raw "
            f"Horizons id — a bare id like '4' resolves to Mars barycenter (RYA-394).")
    key = (round(float(mjd), 6), body_key, obs_code)
    if key not in _HZ_CACHE:
        eph = Horizons(id=spec['id'], location=obs_code, epochs=float(mjd) + 2400000.5,
                       id_type=spec['id_type']).ephemerides()
        targ = str(eph['targetname'][0])
        # GUARD (RYA-374): resolved target MUST be the intended body — not Mars-for-Vesta.
        if spec['match'].lower() not in targ.lower():
            raise BodyIDError(
                f"Horizons resolved body_key={body_key!r} (id={spec['id']!r}, "
                f"id_type={spec['id_type']!r}) to {targ!r}, which does not match "
                f"{spec['match']!r}. Refusing — the bare-'4'⇒Mars-barycenter trap.")
        _HZ_CACHE[key] = (float(eph['r_rate'][0]), float(eph['delta_rate'][0]), targ)
    r_rate, delta_rate, targ = _HZ_CACHE[key]
    return {'v_helio': r_rate, 'v_obs': delta_rate, 'v_total': r_rate + delta_rate,
            'targetname': targ}


def assert_rest_frame(measured_residual_kms: float, line_id: str,
                      tol_kms: float = PASS_TOL_KMS) -> None:
    """CLOSED LOOP — the empirical assert RYA-372 specced but left report-only. A nonzero
    residual after RV correction means the applied velocity (hence the body ID or sign)
    is wrong → loud-fail RestFrameError, never report-and-continue. A non-finite residual
    is NOT a rest failure (it is an unmeasurable / insufficient case handled upstream)."""
    if np.isfinite(measured_residual_kms) and abs(measured_residual_kms) > tol_kms:
        raise RestFrameError(
            f"{line_id} lands at {measured_residual_kms:+.3f} km/s after RV correction "
            f"(tol ±{tol_kms} km/s) — NOT at rest. The reflected-solar RV is wrong "
            f"(body-ID? sign?). This is the open assert that let the Mars-for-Vesta swap "
            f"ride four tickets (RYA-394).")


# ── Telluric velocity (iSpec CCF) — the header-independent frame reference ────
_ISPEC = {}


def _ispec():
    """Lazily bootstrap iSpec (heavy import); cache the module + telluric mask."""
    if 'mod' not in _ISPEC:
        from config.constants import ISPEC_DIR
        sys.path.insert(0, str(ISPEC_DIR))
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            import ispec
        _ISPEC['mod'] = ispec
        _ISPEC['mask'] = ispec.read_cross_correlation_mask(str(ISPEC_DIR.joinpath(*TELLURIC_CCF_MASK)))
    return _ISPEC['mod'], _ISPEC['mask']


def telluric_velocity(wave_air: np.ndarray, flux: np.ndarray) -> float:
    """Telluric line-of-sight velocity (km/s) via iSpec cross-correlation against the
    Synth.Tellurics mask. This is the *applied frame shift* of the spectrum: telluric
    absorption is at rest in the topocentric frame, so its measured velocity is exactly
    how far the pipeline moved the wavelength scale. NaN if iSpec/the mask is
    unavailable or the CCF finds no telluric peak (e.g. a blue arm with no telluric
    band) — the caller flags 'closure unavailable', never a silent pass.
    """
    try:
        ispec, mask = _ispec()
        spec = ispec.create_spectrum_structure(wave_air / 10.0)   # iSpec waveobs in nm
        spec['flux'] = flux
        spec['err'] = 1.0
        models, _ccf = ispec.cross_correlate_with_mask(
            spec, mask, lower_velocity_limit=-120, upper_velocity_limit=120,
            velocity_step=0.5, mask_depth=0.01, fourier=False)
        return float(models[0].mu()) if models else np.nan
    except Exception:
        return np.nan


# ── Line-core velocity measurement ────────────────────────────────────────────

def _core_velocity(wave_air: np.ndarray, flux: np.ndarray, rest: float,
                   win: float = 0.22, depth_min: float = 0.04,
                   center: float = None) -> float:
    """Velocity (km/s) of the absorption core relative to `rest`, via a 5-point
    parabola fit to the deepest pixel within ±`win` of `center` (defaults to `rest`).

    `center` ≠ `rest` lets a coarse bulk shift be applied before a narrow refine,
    so lines shifted by tens of km/s (Vesta reflected RV → up to ~0.6 Å) are still
    found. NaN if the window is empty, the line is too shallow (not present), or the
    minimum sits on the window edge (an incomplete/blended profile).
    """
    c0 = rest if center is None else center
    m = (wave_air > c0 - win) & (wave_air < c0 + win)
    if m.sum() < 5:
        return np.nan
    w, f = wave_air[m], flux[m]
    cont = np.nanmedian(np.concatenate([f[:2], f[-2:]]))
    i = int(np.argmin(f))
    if i < 2 or i > len(f) - 3:
        return np.nan
    if not np.isfinite(cont) or cont <= 0 or (cont - f[i]) / cont < depth_min:
        return np.nan
    c = np.polyfit(w[i - 2:i + 3], f[i - 2:i + 3], 2)
    if c[0] <= 0:
        return np.nan
    lam = -c[1] / (2 * c[0])
    if not (w[i - 2] < lam < w[i + 2]):
        return np.nan
    return (lam - rest) / rest * C_KMS


def _robust_median(vs: np.ndarray, clip: float = 2.5) -> tuple:
    """(median, std, n) after one MAD sigma-clip; (nan,nan,0) if empty."""
    v = vs[np.isfinite(vs)]
    if v.size == 0:
        return np.nan, np.nan, 0
    if v.size >= 3:
        med = np.median(v)
        mad = np.median(np.abs(v - med)) * 1.4826
        if mad > 0:
            v = v[np.abs(v - med) <= clip * mad]
    return float(np.median(v)), (float(np.std(v)) if v.size >= 2 else np.nan), int(v.size)


def measure_bulk_velocity(wave_air: np.ndarray, flux: np.ndarray, lines=VEL_LINES_AIR,
                          clip: float = 2.5) -> dict:
    """Robust bulk LOS velocity (km/s) from clean photospheric cores — two-pass.

    Pass 1 (coarse): a WIDE ±0.8 Å search (~±40 km/s) catches lines shifted by the
    large, epoch-variable Vesta reflected RV. Pass 2 (fine): re-measure each line in
    a NARROW window centred on the coarse-predicted position, so the precise core is
    fit without a neighbouring line stealing the minimum. Each pass is MAD
    sigma-clipped, then median. v_med is NaN if fewer than MIN_LINES clean lines fall
    in range (caller flags INSUFFICIENT — never a silent pass).
    """
    in_range = [x for x in lines if wave_air.min() < x < wave_air.max()]
    # Pass 1 — coarse bulk shift. Photospheric cores (±0.8 Å) plus strong anchors
    # (±1.8 Å) so a frame with a large shift and few clean photospheric cores still
    # seeds correctly. Strong anchors only feed the coarse median.
    anchors = [a for a in STRONG_ANCHORS_AIR if wave_air.min() < a < wave_air.max()]
    vc = np.array([_core_velocity(wave_air, flux, x, win=0.8) for x in in_range]
                  + [_core_velocity(wave_air, flux, a, win=1.8, depth_min=0.06) for a in anchors],
                  dtype=float)
    v_coarse, _, n_coarse = _robust_median(vc, clip)
    if not np.isfinite(v_coarse):
        return {'v_med': np.nan, 'v_std': np.nan, 'n_used': 0,
                'n_in_range': len(in_range), 'v_coarse': np.nan, 'per_line': {}}
    # Pass 2 — fine, narrow window centred on the coarse-shifted position.
    vf = np.array([_core_velocity(wave_air, flux, x, win=0.22,
                                  center=x * (1.0 + v_coarse / C_KMS))
                   for x in in_range], dtype=float)
    v_med, v_std, n_used = _robust_median(vf, clip)
    return {
        'v_med': v_med if n_used >= MIN_LINES else np.nan,
        'v_std': v_std, 'n_used': n_used, 'n_in_range': len(in_range),
        'v_coarse': round(v_coarse, 3),
        'per_line': {x: float(vv) for x, vv in zip(in_range, vf)},
    }


# ── Frame loading ─────────────────────────────────────────────────────────────

def _load_frame(path: Path) -> dict:
    """Load one IDP frame to the AIR wavelength scale, respecting the RYA-271 traps.

    ESPRESSO: `WAVE_AIR` (the `WAVE` column is vacuum); reject S1D_STACK_A
    (velocity-smeared). UVES: single `WAVE` is air. Wavelengths are read from the
    data COLUMN (Å), never the nm header keyword. Returns a dict with the spectrum
    + metadata, or {'reject': reason} for a frame that must not be conditioned.
    """
    with fits.open(path) as h:
        ph = h[0].header
        inst = ph.get('INSTRUME', '?')
        pro_catg = str(ph.get('HIERARCH ESO PRO CATG', '?'))
        cols = [c.name for c in h[1].columns] if len(h) > 1 and hasattr(h[1], 'columns') else []
        meta = {
            'file': path.name, 'instrument': inst, 'pro_catg': pro_catg,
            'specsys': ph.get('SPECSYS', '?'),
            'ins_mode': ph.get('HIERARCH ESO INS MODE', ph.get('HIERARCH ESO INS PATH', '?')),
            'mjd_mid': float(ph.get('MJD-OBS', np.nan)) + float(ph.get('EXPTIME', 0.0)) / 2 / 86400.0,
            'date_obs': str(ph.get('DATE-OBS', '?'))[:19],
            'berv': ph.get('HIERARCH ESO QC BERV', None),
            'snr': float(ph.get('SNR', np.nan)),
        }
        # ESPRESSO STACK_A → velocity-smeared multi-epoch coadd (RYA-370). Reject loud.
        if inst == 'ESPRESSO' and 'STACK' in pro_catg:
            meta['reject'] = f'STACK_A velocity-smeared (use S1D_FINAL_A); pro_catg={pro_catg}'
            return meta
        if 'WAVE_AIR' in cols:                 # ESPRESSO: vacuum WAVE + air WAVE_AIR
            wave_air = np.asarray(h[1].data['WAVE_AIR']).ravel().astype(float)
        elif 'WAVE' in cols:                   # UVES: air WAVE
            wave_air = np.asarray(h[1].data['WAVE']).ravel().astype(float)
        else:
            meta['reject'] = 'no WAVE/WAVE_AIR column'
            return meta
        meta['wave_air'] = wave_air
        meta['flux'] = np.asarray(h[1].data['FLUX']).ravel().astype(float)
    return meta


def discover_set(set_name: str, root: Path = DEFAULT_ROOT) -> list:
    """Frame paths for a named set. 'vesta_espresso' → ESPRESSO optical IDPs;
    'vesta_uves' → UVES optical IDPs (per-dichroic). STACK_A is included here and
    rejected-with-reason downstream (so the rejection is logged, not hidden)."""
    files = sorted(root.glob('ADP.*.fits'))
    out = []
    for f in files:
        inst = fits.getheader(f).get('INSTRUME', '?')
        if set_name == 'vesta_espresso' and inst == 'ESPRESSO':
            out.append(f)
        elif set_name == 'vesta_uves' and inst == 'UVES':
            out.append(f)
    return out


# ── Per-frame conditioning ────────────────────────────────────────────────────

def condition_frame(path: Path, obs_code: str = PARANAL,
                    body_key: str = DEFAULT_BODY, closure: bool = True) -> dict:
    """Condition one frame to the solar rest frame and verify on held-out lines.

    Anchor velocity = robust median of TRAIN photospheric cores. Apply
    λ_rest = λ_obs / (1 + v_anchor/C). Verify on the disjoint TEST set: the held-out
    photospheric residual must be ≤ PASS_TOL_KMS or the frame loud-fails RestFrameError
    (RYA-394 — the closed loop; never silently passed or report-only). Na D2 / Hα
    residuals are reported for context (they carry a small chromospheric core offset, so
    they are NOT asserted; the held-out PHOTOSPHERIC set is the rigorous rest proof). The
    Horizons two-leg model (`body_key`, single-source registry) is the independent
    cross-check — now Vesta, not the bare-'4'⇒Mars value (RYA-394).

    `closure` adds the telluric-CCF wavelength-zero-point check (review item 2): the
    measured telluric velocity must equal the pipeline's claimed applied frame shift
    (BERV for BARYCENT, 0 for TOPOCENT) within TELL_CLOSURE_TOL, else the frame is
    downgraded to CLOSURE_FAIL — this is the only guard against a common-mode offset
    the held-out check cannot see.
    """
    fr = _load_frame(Path(path))
    rec = {k: fr.get(k) for k in ('file', 'instrument', 'pro_catg', 'ins_mode',
                                  'specsys', 'date_obs', 'mjd_mid', 'snr', 'berv')}
    if 'reject' in fr:
        rec['status'] = 'REJECTED'
        rec['reason'] = fr['reject']
        return rec

    wave, flux = fr['wave_air'], fr['flux']

    # Anchor on TRAIN lines; the applied correction is the measured bulk velocity.
    train = measure_bulk_velocity(wave, flux, TRAIN_LINES)
    rec['v_anchor'] = train['v_med']
    rec['n_train'] = train['n_used']

    # Independent physical cross-check (Molaro two-leg) — body from the single-source
    # registry (RYA-394); BodyIDError (wrong-target) must surface, never be swallowed.
    try:
        eph = reflected_solar_rv(fr['mjd_mid'], body_key, obs_code=obs_code)
        rec.update({'v_helio': round(eph['v_helio'], 3), 'v_obs': round(eph['v_obs'], 3),
                    'v_total_eph': round(eph['v_total'], 3), 'eph_target': eph['targetname']})
    except BodyIDError:
        raise                                                 # body-ID bug is never silent
    except Exception as exc:                                   # transient Horizons/network
        rec.update({'v_helio': np.nan, 'v_obs': np.nan, 'v_total_eph': np.nan,
                    'eph_error': str(exc)[:80]})

    if not np.isfinite(train['v_med']):
        rec['status'] = 'INSUFFICIENT'
        rec['reason'] = (f'only {train["n_used"]} clean photospheric lines in range '
                         f'(need ≥{MIN_LINES}); SNR={fr["snr"]:.0f}. Not conditioned — flagged.')
        return rec

    v = train['v_med']
    wave_rest = wave / (1.0 + v / C_KMS)

    # Verify on the held-out TEST set + the named chromospheric lines.
    test = measure_bulk_velocity(wave_rest, flux, TEST_LINES)
    rec['n_test'] = test['n_used']
    rec['resid_test'] = round(test['v_med'], 3) if np.isfinite(test['v_med']) else np.nan
    rec['resid_NaD2'] = round(_core_velocity(wave_rest, flux, NAD2_AIR, 0.7), 3) \
        if wave_rest.min() < NAD2_AIR < wave_rest.max() else np.nan
    rec['resid_Halpha'] = round(_core_velocity(wave_rest, flux, HALPHA_AIR, 0.7), 3) \
        if wave_rest.min() < HALPHA_AIR < wave_rest.max() else np.nan
    rec['eph_minus_meas'] = round(rec['v_total_eph'] - v, 3) if np.isfinite(rec.get('v_total_eph', np.nan)) else np.nan

    rec['wave_rest'] = wave_rest          # for an optional writer (RYA-371)
    if test['n_used'] < MIN_LINES or not np.isfinite(test['v_med']):
        # Cannot VERIFY (data gap) — distinct from an off-rest failure; honestly flagged,
        # not a silent pass (the frame is not handed downstream as conditioned).
        rec['status'] = 'INSUFFICIENT'
        rec['reason'] = f'held-out set has {test["n_used"]} clean lines (<{MIN_LINES}) — cannot verify.'
        return rec
    # CLOSED-LOOP rest-frame assert (RYA-394): the held-out PHOTOSPHERIC set must land at
    # rest, else the applied velocity (body-ID? sign?) is wrong. Loud-fail — no report-only
    # path. verify_set catches RestFrameError (via `.rec`) to keep surveying.
    if abs(test['v_med']) > PASS_TOL_KMS:
        rec['status'] = 'CRITICAL'
        rec['reason'] = (f'held-out photospheric lines land at {test["v_med"]:+.3f} km/s '
                         f'(> {PASS_TOL_KMS} km/s) after conditioning — NOT at rest.')
        try:
            assert_rest_frame(test['v_med'], f"{rec['file']}: held-out photospheric set")
        except RestFrameError as e:
            e.rec = rec
            raise
    else:
        rec['status'] = 'PASS'

    # ── Telluric-CCF wavelength-zero-point closure (review item 2) ────────────
    # Measured on the AS-DELIVERED wave (the telluric tells us the applied frame
    # shift). Expected = the pipeline's claimed shift: BERV for BARYCENT, 0 for
    # TOPOCENT. Mismatch > tol = common-mode offset / frame mislabel → loud-fail.
    # (Reached only for PASS frames; an off-rest frame already raised RestFrameError.)
    if closure:
        v_tell = telluric_velocity(wave, flux)
        rec['v_telluric'] = round(v_tell, 3) if np.isfinite(v_tell) else np.nan
        specsys = str(fr.get('specsys', ''))
        berv = fr.get('berv', None)
        if specsys == 'BARYCENT' and berv is not None:
            expected = float(berv)
        elif specsys == 'TOPOCENT':
            expected = 0.0                      # validated: UVES telluric ≈ 0 (RYA-372)
        else:
            expected = np.nan
        rec['tell_expected'] = round(expected, 3) if np.isfinite(expected) else np.nan
        if np.isfinite(v_tell) and np.isfinite(expected):
            rec['closure_resid'] = round(v_tell - expected, 3)
            # Reflected-solar geometry cross-check (NOT gated): the measured
            # topocentric reflected velocity (v_anch − v_telluric) vs the Horizons
            # two-leg. The large/variable gap is the leading-order model limit.
            if np.isfinite(rec.get('v_total_eph', np.nan)):
                rec['geom_gap'] = round((v - v_tell) - rec['v_total_eph'], 3)
            if abs(rec['closure_resid']) > TELL_CLOSURE_TOL and rec['status'] == 'PASS':
                rec['status'] = 'CLOSURE_FAIL'
                rec['reason'] = (f'telluric at {v_tell:+.2f} km/s vs claimed frame shift '
                                 f'{expected:+.2f} ({specsys}); closure {v_tell-expected:+.2f} '
                                 f'> {TELL_CLOSURE_TOL} km/s — common-mode wavelength offset.')
        else:
            rec['closure_resid'] = np.nan
            rec['closure_note'] = 'telluric CCF unavailable (no telluric band / iSpec) — closure not checked'

    return rec


# ── Verify a whole set (smoke test) ───────────────────────────────────────────

_FAIL_STATES = ('CRITICAL', 'CLOSURE_FAIL')


def verify_set(set_name: str, root: Path = DEFAULT_ROOT, obs_code: str = PARANAL,
               body_key: str = DEFAULT_BODY, closure: bool = True) -> list:
    paths = discover_set(set_name, root)
    spec = REFLECTED_SOLAR_BODIES[body_key]
    print(f"\n{'='*120}\n  RYA-372 reflected-solar RV conditioning — set '{set_name}'  "
          f"({len(paths)} frames; root={root.name})\n{'='*120}")
    print(f"  ASSERTION: asteroid-ephemeris path used (JPL Horizons, body_key={body_key!r} "
          f"→ id={spec['id']!r}/id_type={spec['id_type']!r} @Paranal={obs_code}); NO bare-'4' "
          f"(=Mars) trap, NO stellar-BERV shortcut.")
    print(f"  Anchor = measured photospheric bulk velocity; CLOSED-LOOP rest assert on "
          f"held-out lines (±{PASS_TOL_KMS} km/s, RestFrameError).")
    print(f"  CLOSURE (teeth): telluric CCF must match the claimed frame shift "
          f"(BERV/0) within ±{TELL_CLOSURE_TOL} km/s — guards the common-mode offset the held-out check cannot.\n")
    hdr = (f"  {'file':32s} {'mode':9s} {'SPECSYS':9s} {'v_anch':>7s} {'v_eph':>7s} "
           f"{'v_tell':>7s} {'expect':>7s} {'closure':>8s} {'geomgap':>8s} {'test':>6s} {'status':12s}")
    print(hdr + "\n  " + "-" * 118)
    recs = []
    for p in paths:
        # condition_frame loud-fails RestFrameError on an off-rest frame (no report-only
        # path); the survey catches it (via .rec) to tabulate every frame as CRITICAL.
        try:
            r = condition_frame(p, obs_code=obs_code, body_key=body_key, closure=closure)
        except RestFrameError as e:
            r = e.rec if e.rec is not None else {'file': Path(p).name, 'status': 'CRITICAL',
                                                 'reason': str(e)[:60]}
        recs.append(r)
        if r['status'] == 'REJECTED':
            print(f"  {r['file']:32s} {str(r.get('ins_mode','')):9s} {str(r.get('specsys','')):9s} "
                  f"{'':>7s} {'':>7s} {'':>7s} {'':>7s} {'':>8s} {'':>8s} {'':>6s} "
                  f"REJECTED ({r['reason'][:34]})")
            continue

        def g(k, w=7):
            v = r.get(k)
            return f"{v:+.2f}" if isinstance(v, (int, float)) and np.isfinite(v) else '—'
        note = ''
        if r['status'] in ('CRITICAL', 'INSUFFICIENT', 'CLOSURE_FAIL'):
            note = ' ' + str(r.get('reason', ''))[:38]
        elif r.get('closure_note'):
            note = ' (closure n/a)'
        print(f"  {r['file']:32s} {str(r.get('ins_mode','')):9s} {str(r.get('specsys','')):9s} "
              f"{g('v_anchor'):>7s} {g('v_total_eph'):>7s} {g('v_telluric'):>7s} {g('tell_expected'):>7s} "
              f"{g('closure_resid'):>8s} {g('geom_gap'):>8s} {g('resid_test'):>6s} {r['status']:12s}{note}")

    from collections import Counter
    c = Counter(r['status'] for r in recs)
    print("  " + "-" * 118)
    print("  summary: " + "  ".join(f"{k}={v}" for k, v in sorted(c.items())))
    clo = np.array([r['closure_resid'] for r in recs
                    if r.get('status') == 'PASS'
                    and np.isfinite(r.get('closure_resid', np.nan))])
    if clo.size:
        print(f"  closure residual (telluric − claimed frame shift): median={np.median(clo):+.3f}  "
              f"std={np.std(clo):.3f}  max|·|={np.max(np.abs(clo)):.3f} km/s  (n={clo.size}; tol ±{TELL_CLOSURE_TOL})")
    fails = [r for r in recs if r['status'] in _FAIL_STATES]
    if fails:
        print(f"  ⚠️  {len(fails)} loud-fail frame(s) ({', '.join(sorted({r['status'] for r in fails}))}) "
              f"— NOT silently passed.")
    return recs


# ── --write: emit rest-frame conditioned spectra for RYA-371 ──────────────────

def write_set(set_name: str, out_dir: Path, recs: list = None, root: Path = DEFAULT_ROOT,
              obs_code: str = PARANAL) -> list:
    """Write per-frame rest-frame spectra for PASS frames (CRITICAL / CLOSURE_FAIL /
    INSUFFICIENT / REJECTED excluded, logged in the manifest). Per-frame only — never
    coadded (Vesta 5.3 h rotation). Each CSV carries a provenance header (applied
    v_anch, anchor set, held-out + closure residuals, the convective-blueshift note);
    columns are `wavelength_air_A,flux` (air per instrument convention). Reuses
    pre-computed `recs` (from verify_set) when given, to avoid a second conditioning pass."""
    import csv
    out_dir = Path(out_dir)
    if recs is None:
        recs = [condition_frame(p, obs_code=obs_code, closure=True)
                for p in discover_set(set_name, root)]
    manifest = []
    written = 0
    for r in recs:
        p = root / r['file']
        if r['status'] != 'PASS':
            manifest.append({**{k: r.get(k) for k in ('file', 'instrument', 'ins_mode',
                                 'specsys', 'status', 'v_anchor', 'closure_resid')},
                             'output': '', 'reason': r.get('reason', r.get('closure_note', ''))})
            continue
        inst_dir = out_dir / f'vesta_{str(r["instrument"]).lower()}'
        inst_dir.mkdir(parents=True, exist_ok=True)
        out_path = inst_dir / f'{Path(r["file"]).stem}_rest.csv'
        wave_rest = r['wave_rest']
        flux = _load_frame(Path(p))['flux']
        prov = [
            f'# RYA-372 reflected-solar rest-frame conditioned spectrum',
            f'# source: {r["file"]}  ({r["instrument"]} {r.get("ins_mode","")}, {r.get("specsys","")})',
            f'# date-obs: {r.get("date_obs","")}  mjd_mid: {r.get("mjd_mid","")}',
            f'# applied v_anch (removed): {r["v_anchor"]:+.4f} km/s  (lambda_rest = lambda_obs/(1+v/c))',
            f'# anchor: median of {r.get("n_train","?")} TRAIN photospheric Fe I cores (two-pass coarse->fine)',
            f'# held-out residual: {r.get("resid_test")} km/s   NaD2: {r.get("resid_NaD2")}   Ha: {r.get("resid_Halpha")}',
            f'# telluric closure: v_tell={r.get("v_telluric")} expected={r.get("tell_expected")} resid={r.get("closure_resid")} km/s',
            f'# Horizons two-leg cross-check: r+delta={r.get("v_total_eph")}  geom_gap={r.get("geom_gap")} km/s',
            f'# {CONVECTIVE_BLUESHIFT_NOTE}',
            f'# wavelength: AIR Angstrom (ESPRESSO WAVE_AIR / UVES WAVE); flux as delivered (FLUXCAL per header)',
        ]
        with open(out_path, 'w', newline='') as fh:
            fh.write('\n'.join(prov) + '\n')
            w = csv.writer(fh)
            w.writerow(['wavelength_air_A', 'flux'])
            for wl, fx in zip(wave_rest, flux):
                w.writerow([f'{wl:.6f}', f'{fx:.6g}'])
        written += 1
        manifest.append({'file': r['file'], 'instrument': r['instrument'],
                         'ins_mode': r.get('ins_mode'), 'specsys': r.get('specsys'),
                         'status': 'PASS', 'v_anchor': r['v_anchor'],
                         'closure_resid': r.get('closure_resid'),
                         'output': str(out_path.relative_to(out_dir)), 'reason': ''})
    out_dir.mkdir(parents=True, exist_ok=True)
    man_path = out_dir / f'{set_name}_manifest.csv'
    with open(man_path, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=['file', 'instrument', 'ins_mode', 'specsys',
                                           'status', 'v_anchor', 'closure_resid', 'output', 'reason'])
        w.writeheader()
        w.writerows(manifest)
    print(f"\n  [--write] {written} rest-frame spectra → {out_dir}/vesta_*/   "
          f"manifest → {man_path.name}  ({len(manifest)-written} excluded, logged)")
    return recs


def main(argv=None):
    ap = argparse.ArgumentParser(description='Reflected-solar RV conditioning (RYA-372)')
    ap.add_argument('--set', dest='set_name', required=True,
                    choices=['vesta_espresso', 'vesta_uves'])
    ap.add_argument('--verify', action='store_true', help='per-frame RV + closure table')
    ap.add_argument('--write', action='store_true', help='emit rest-frame spectra (PASS frames) + manifest')
    ap.add_argument('--out', default=None, help='output dir for --write (default data/processed/reflected_solar)')
    ap.add_argument('--no-closure', action='store_true', help='skip the telluric-CCF closure check')
    ap.add_argument('--root', default=str(DEFAULT_ROOT))
    ap.add_argument('--obs-code', default=PARANAL)
    args = ap.parse_args(argv)
    recs = verify_set(args.set_name, root=Path(args.root), obs_code=args.obs_code,
                      closure=not args.no_closure)
    if args.write:
        out = Path(args.out) if args.out else (
            Path(__file__).resolve().parent.parent / 'data' / 'processed' / 'reflected_solar')
        write_set(args.set_name, out, recs=recs, root=Path(args.root), obs_code=args.obs_code)
    # Non-zero exit on any loud-fail (CRITICAL / CLOSURE_FAIL) for CI / RYA-371 gating.
    return 1 if any(r['status'] in _FAIL_STATES for r in recs) else 0


if __name__ == '__main__':
    sys.exit(main())
