"""
pipeline/frame_object_contract.py
=================================
The standing loader convention, made executable — the single source of truth for
the things **every** loader must declare and verify about a frame before its
photons reach any fit, EW, or synthesis:

  A.  **OBJECT attribution** (RYA-481) — which star the frame is, by authoritative
      FITS header, NEVER by folder / filename / glob.
  B.  **Velocity frame** (RYA-481) — the frame's velocity reference (SPECSYS,
      applied BERV, systemic RV), with the corrections it needs applied *exactly*
      — no more, no less — and verified against a known line.
  C.  **Wavelength scale** (RYA-264) — the frame's native wavelength UNIT (nm/µm→Å)
      and SCALE (vacuum→air), converted to the pipeline's air-Å at the loader
      boundary from a CITED per-instrument convention, with a unit-sanity band gate
      (the ×10 / ×10000 catch) and a known-line air-rest verification.

D is the principle shared by all three: the correct response to *any* unverifiable
attribution, velocity frame, or wavelength scale is **raise/flag, never silently
proceed on a default**. A wrong-star, wrong-frame, or wrong-scale result is worse
than no result, because it looks right (a plausible, wrong abundance that passes
every check except an explicit one). Silence is the bug; loudness is the fix. All
three raise a common base, ``FrameContractError``.

Why this module exists (the receipts — five disguises of one disease)
--------------------------------------------------------------------
  1. Vesta-as-solar          — wrong SOURCE (reflected sunlight ≠ direct solar)
  2. glob-loading            — wrong STAR (globbed a tree, picked another star)
  3. arm-registry default    — wrong SOURCE (silent default resolved wrong arm)
  4. Procyon tree held 55Cnc — wrong STAR (srho01cnc frames in "Procyon HST")
  5. α Cen folders mislabeled— wrong STAR ("A/HARPS" was 63% B by OBJECT)
  6. UVES O I 777 no sys-RV  — wrong VELOCITY FRAME (BERV applied, systemic not)
     (+ HARPS double-BERV, CRIRES TOPOCENT→bary, UVES BERV sign-check)

Two root classes, one disease: wrong-source/wrong-star (attribution) and
wrong-wavelength (velocity frame). Both put good-looking photons in the wrong
place and let the fit absorb the error into the science number.

How loaders use it
------------------
  from pipeline.frame_object_contract import (
      assert_object, corrections_for_specsys, VelocityFrame, verify_line_position)

  star = assert_object(h0['OBJECT'], expected='alpha_cen_a', context=f"{fn}: ")
  pol  = corrections_for_specsys(h0['SPECSYS'])     # raises on unknown frame
  frame = VelocityFrame(specsys=h0['SPECSYS'], berv_applied=True,
                        systemic_rv_applied=False, wave_units='air',
                        wave_scale='angstrom')
  frame.validate()                                  # raises on the double-BERV trap
  verify_line_position(wave_A, flux, 6562.79, tol_kms=8.0)   # BERV sign-check

References (the incident trail): RYA-272 (BERV sign-check), RYA-426 (UV vac/air),
RYA-431 (NIRPS RV zero-point), RYA-464 (arm registry default), RYA-471 (no-glob
whitelist), RYA-478 (systemic-RV fix), RYA-479 (α Cen OBJECT re-split + no double
BERV), RYA-480 (CRIRES TOPOCENT→bary), RYA-288 (fail-loud guard precedent),
RYA-301/303 (the audits that surfaced the attribution disease). RYA-495 layers a
star-attribution *authority ranking* (RV star-ID > OBJECT header > folder) on top
of this map — register that authority via ``register_aliases`` when it lands.
"""
from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Optional

import numpy as np

from config.constants import PHYSICS
from pipeline.wavelength_util import vac_to_air   # the ONE vac→air impl (RYA-264/426)

_C_KMS = PHYSICS['c_kms']


class FrameContractError(RuntimeError):
    """Base for every loader-contract failure — OBJECT attribution, velocity frame,
    or wavelength scale. Catch this to handle any contract violation; the axis
    subclasses (``ObjectAttributionError``, ``VelocityFrameError``) narrow it. The
    wavelength-scale axis raises ``FrameContractError`` directly (RYA-264)."""


# ════════════════════════════════════════════════════════════════════════════
# A.  OBJECT attribution — central canonical map (never per-loader, never folder)
# ════════════════════════════════════════════════════════════════════════════

class ObjectAttributionError(FrameContractError):
    """Raised when a frame's OBJECT cannot be attributed to the expected star —
    unresolved, ambiguous, or a confident mismatch. Names the raw OBJECT and the
    reason; never a silent include or a folder-name fallback (RYA-481 §A)."""


# Canonical star ids match config/stars.yaml keys. Each maps to the set of
# normalized (lower-case, alphanumeric-only) OBJECT aliases that designate it:
# HD numbers, Bayer/proper names, HR numbers, catalogue strings. Maintained
# centrally here so a new alias is fixed once for every loader and audit.
_STAR_ALIASES: dict[str, set[str]] = {
    'solar':       {'sun', 'sol', 'thesun', 'solar'},
    'procyon':     {'procyon', 'hd61421', 'hd061421', 'alfcmi', 'alphacmi',
                    'alfacmi', 'hr2943', 'gj280'},
    'alpha_cen_a': {'hd128620', 'alfcena', 'alphacena', 'acena', 'rigilkentaurus',
                    'rigilkent', 'hr5459', 'gj559a'},
    'alpha_cen_b': {'hd128621', 'alfcenb', 'alphacenb', 'acenb', 'toliman',
                    'hr5460', 'gj559b'},
    '55cnc_a':     {'55cnc', '55cnca', '55cancri', 'rhocnc', 'rho1cnc', 'rho01cnc',
                    'srho01cnc', 'hd75732', 'hr3522', 'gj324a'},
}

# Reverse index, built once; refuses to register the same alias to two stars.
_ALIAS_TO_STAR: dict[str, str] = {}


def _index_aliases() -> None:
    _ALIAS_TO_STAR.clear()
    for star, aliases in _STAR_ALIASES.items():
        for a in aliases:
            prev = _ALIAS_TO_STAR.get(a)
            if prev is not None and prev != star:
                raise ObjectAttributionError(
                    f"alias {a!r} is ambiguous between {prev!r} and {star!r} — the "
                    f"canonical map must be one-to-one (RYA-481 §A).")
            _ALIAS_TO_STAR[a] = star


_index_aliases()


def register_aliases(star: str, aliases) -> None:
    """Add OBJECT aliases for a canonical star (extends the central map at import
    time — e.g. a new benchmark, or RYA-495's authority layer). Re-indexes and
    refuses a collision with another star, rolling back the new aliases so a
    rejected call never pollutes the map."""
    bucket = _STAR_ALIASES.setdefault(star, set())
    added = {normalize_object(a) for a in aliases} - bucket
    bucket.update(added)
    try:
        _index_aliases()
    except ObjectAttributionError:
        bucket.difference_update(added)
        _index_aliases()
        raise


def normalize_object(obj) -> str:
    """Collapse an OBJECT/target string to a case- and separator-insensitive key
    (lower-case, alphanumeric only): ``'alf Cen A' -> 'alfcena'``, ``'HD 128620'
    -> 'hd128620'``."""
    return re.sub(r'[^a-z0-9]', '', str(obj).lower())


def object_tokens(obj) -> frozenset:
    """Normalized tokens of a composite OBJECT (split on separators). Handles the
    multi-fibre HELIOS solar header ``'SUN,FP,G2V'`` -> ``{'sun','fp','g2v'}``."""
    return frozenset(normalize_object(t) for t in re.split(r'[^A-Za-z0-9]+', str(obj)) if t)


@dataclass(frozen=True)
class ObjectResolution:
    """Result of attributing an OBJECT string to a canonical star."""
    raw: str
    canonical: Optional[str]   # canonical star id, or None if unresolved/ambiguous
    resolved: bool
    reason: str


def resolve_object(obj) -> ObjectResolution:
    """Attribute a raw OBJECT header to a canonical star id by the central map.
    Resolves exact aliases, then per-token aliases (composite headers), then a
    documented ``cen``-suffix fallback for α Cen name variants. Returns an
    unresolved result (never a guess) on a missing, ambiguous, or unrecognized
    OBJECT — the caller decides whether that quarantines or raises."""
    raw = '' if obj is None else str(obj)
    if not raw.strip():
        return ObjectResolution(raw, None, False, 'missing OBJECT')

    full = normalize_object(raw)
    # 1) exact full-string alias
    if full in _ALIAS_TO_STAR:
        star = _ALIAS_TO_STAR[full]
        return ObjectResolution(raw, star, True, f"alias {full!r} -> {star}")

    # 2) per-token aliases (composite OBJECT like 'SUN,FP,G2V'); agreement required
    matches = {_ALIAS_TO_STAR[t] for t in object_tokens(raw) if t in _ALIAS_TO_STAR}
    if len(matches) == 1:
        star = next(iter(matches))
        return ObjectResolution(raw, star, True, f"token -> {star}")
    if len(matches) > 1:
        return ObjectResolution(raw, None, False,
                                f"ambiguous: tokens resolve to {sorted(matches)} ({raw!r})")

    # 3) documented α Cen name-suffix fallback (variants not in the explicit set)
    if 'cen' in full and 'cancri' not in full:
        if full.endswith('a'):
            return ObjectResolution(raw, 'alpha_cen_a', True, 'cen-name suffix -> A')
        if full.endswith('b'):
            return ObjectResolution(raw, 'alpha_cen_b', True, 'cen-name suffix -> B')
        return ObjectResolution(raw, None, False, f"ambiguous cen-name ({raw!r})")

    return ObjectResolution(raw, None, False, f"unrecognized OBJECT ({raw!r})")


def canonical_star_id(name: str) -> str:
    """Resolve an *expected*-star argument to a canonical id. Accepts a canonical
    id (``'alpha_cen_a'``) or any OBJECT alias; raises if it is not a known star."""
    if name in _STAR_ALIASES:
        return name
    keyed = re.sub(r'[\s\-]+', '_', str(name).strip().lower())
    if keyed in _STAR_ALIASES:
        return keyed
    res = resolve_object(name)
    if res.resolved:
        return res.canonical
    raise ObjectAttributionError(
        f"expected star {name!r} is not a known canonical id or alias "
        f"({sorted(_STAR_ALIASES)}) — add it to the central map (RYA-481 §A).")


def assert_object(obj, expected, *, context: str = '',
                  exc: type = ObjectAttributionError) -> str:
    """Verify a frame's OBJECT attributes to ``expected`` by header, and return the
    canonical star id. Raises ``exc`` (default ``ObjectAttributionError``; pass a
    loader-specific subclass to keep one error vocabulary) when the OBJECT is
    unresolvable/ambiguous or resolves to a *different* star. This is rule §A in
    one call — a loader that can't confirm the star fails loud, it does not
    default, glob, or trust the path."""
    exp = canonical_star_id(expected)
    res = resolve_object(obj)
    if not res.resolved:
        raise exc(
            f"{context}OBJECT={obj!r} could not be attributed to a known star "
            f"({res.reason}). Refusing a folder/filename fallback — quarantine + "
            f"verify by header (RYA-481 §A).")
    if res.canonical != exp:
        raise exc(
            f"{context}OBJECT={obj!r} resolves to {res.canonical!r}, expected {exp!r} "
            f"({res.reason}). A frame's star is its header, not its path (RYA-481 §A).")
    return res.canonical


# ════════════════════════════════════════════════════════════════════════════
# B.  Velocity frame — declared, corrected exactly, verified
# ════════════════════════════════════════════════════════════════════════════

class VelocityFrameError(FrameContractError):
    """Raised on an unrecognized velocity frame, a contradictory correction (the
    double-BERV trap), or a known line that lands off its rest wavelength after
    correction (a missing/sign-wrong term) — RYA-481 §B."""


@dataclass(frozen=True)
class FrameCorrection:
    """The correction a SPECSYS frame needs, from the project's incident record."""
    berv_needed: bool
    note: str


# The per-SPECSYS policy. ``berv_needed`` is whether the *loader* must apply the
# barycentric correction (True) or the flux is already barycentric and re-applying
# it would double-correct (False). Systemic RV is separate and frame-independent
# (see ``VelocityFrame``): BERV ≠ systemic RV; needing one does not mean the other
# is done (RYA-478).
_SPECSYS_POLICY: dict[str, FrameCorrection] = {
    'TOPOCENT': FrameCorrection(
        berv_needed=True,
        note='observer frame — loader MUST apply BERV (UVES IDP RYA-272; CRIRES IDP RYA-480)'),
    'BARYCENT': FrameCorrection(
        berv_needed=False,
        note='already barycentric — loader must NOT re-apply BERV (HARPS S1D; double-BERV trap RYA-479)'),
    'HELIOCEN': FrameCorrection(
        berv_needed=False,
        note='heliocentric ≈ barycentric — do NOT apply BERV (note: GES HELIOCEN UVES are quarantined, RYA-272)'),
}


def corrections_for_specsys(specsys, *, exc: type = VelocityFrameError) -> FrameCorrection:
    """Return the correction policy for a SPECSYS value, raising on an unrecognized
    frame rather than guessing whether BERV is needed (RYA-481 §B)."""
    key = str(specsys).strip().upper()
    if key not in _SPECSYS_POLICY:
        raise exc(
            f"SPECSYS={specsys!r} is not a recognized velocity reference frame. "
            f"Declare it explicitly before loading — refusing to guess whether BERV "
            f"is needed (RYA-481 §B). Known: {sorted(_SPECSYS_POLICY)}.")
    return _SPECSYS_POLICY[key]


_WAVE_UNITS = ('air', 'vacuum')
_WAVE_SCALES = ('angstrom', 'nm')


@dataclass(frozen=True)
class VelocityFrame:
    """A loader's explicit declaration of the frame it is returning. Stamp it into
    ``meta['velocity_frame']`` (via :meth:`declare`) so every spectrum carries its
    frame provenance, and call :meth:`validate` to catch a correction that
    contradicts the SPECSYS policy (the double-BERV / under-correction trap)."""
    specsys: str                 # raw header SPECSYS (or the loader's declared frame)
    berv_applied: bool           # has the loader applied the barycentric correction?
    systemic_rv_applied: bool    # has it been shifted to the stellar rest frame?
    wave_units: str              # 'air' | 'vacuum'
    wave_scale: str = 'angstrom'  # 'angstrom' | 'nm'
    note: str = ''

    def __post_init__(self):
        if self.wave_units not in _WAVE_UNITS:
            raise VelocityFrameError(
                f"wave_units={self.wave_units!r} must be one of {_WAVE_UNITS} "
                f"(declare vacuum vs air at the loader boundary — RYA-481 §B.4).")
        if self.wave_scale not in _WAVE_SCALES:
            raise VelocityFrameError(
                f"wave_scale={self.wave_scale!r} must be one of {_WAVE_SCALES}.")

    def validate(self) -> 'VelocityFrame':
        """Raise if the applied corrections contradict the SPECSYS policy: BERV
        applied to an already-barycentric frame (double-correction), or BERV not
        applied to a topocentric frame (left in the observer frame). Returns self
        so it chains. Loaders that *intend* to leave a frame uncorrected (and warn)
        simply don't call this."""
        pol = corrections_for_specsys(self.specsys)
        if self.berv_applied and not pol.berv_needed:
            raise VelocityFrameError(
                f"double-correction: BERV applied to a {self.specsys!r} frame that is "
                f"already barycentric ({pol.note}). This blueshifts/doubles the "
                f"correction — exactly the RYA-479 trap (RYA-481 §B.2).")
        if (not self.berv_applied) and pol.berv_needed:
            raise VelocityFrameError(
                f"under-correction: a {self.specsys!r} frame needs BERV but none was "
                f"applied ({pol.note}) — wavelengths are still observer-frame; do NOT "
                f"co-add or fit lines (RYA-481 §B.2).")
        return self

    def declare(self) -> str:
        """One-line human declaration for ``meta`` / docstrings."""
        rest = 'stellar-rest' if self.systemic_rv_applied else (
            'barycentric' if self.berv_applied else self.specsys.lower())
        return (f"{self.specsys} → {rest}; BERV {'applied' if self.berv_applied else 'not applied'}; "
                f"systemic-RV {'applied' if self.systemic_rv_applied else 'not applied'}; "
                f"{self.wave_units} {self.wave_scale}"
                f"{' — ' + self.note if self.note else ''}")


def line_velocity_kms(wave_A, flux, expected_A: float, *, window_A: float = 3.0,
                      absorption: bool = True) -> tuple[float, float]:
    """Measure where a known line actually lands. Returns ``(v_kms, centroid_A)``
    where ``v_kms`` is the velocity offset of the flux-weighted line centroid from
    its rest wavelength. The contrast weighting mirrors the RYA-272 Hα sign-check.
    Raises ``VelocityFrameError`` if there is too little signal to centroid."""
    w = np.asarray(wave_A, dtype=float).ravel()
    f = np.asarray(flux, dtype=float).ravel()
    m = np.isfinite(w) & np.isfinite(f) & (w > expected_A - window_A) & (w < expected_A + window_A)
    if int(m.sum()) < 3:
        raise VelocityFrameError(
            f"only {int(m.sum())} finite pixels within ±{window_A} Å of {expected_A} Å — "
            f"cannot locate the verification line (wrong coverage or all-NaN window).")
    w, f = w[m], f[m]
    if absorption:
        weight = np.clip(np.nanpercentile(f, 90) - f, 0.0, None)   # depth below pseudo-continuum
    else:
        weight = np.clip(f - np.nanpercentile(f, 10), 0.0, None)   # height above baseline
    if not np.any(weight > 0):
        raise VelocityFrameError(
            f"no line contrast within ±{window_A} Å of {expected_A} Å — cannot centroid "
            f"(is the line actually present in this frame?).")
    centroid = float(np.sum(w * weight) / np.sum(weight))
    return _C_KMS * (centroid - expected_A) / expected_A, centroid


def verify_line_position(wave_A, flux, expected_A: float, *, tol_kms: float = 5.0,
                         window_A: float = 3.0, absorption: bool = True,
                         context: str = '', exc: type = VelocityFrameError) -> float:
    """Confirm a known line lands at (or correctly shifted to) its expected
    wavelength after correction — the BERV sign-check pattern (RYA-272: Hα
    57.4→4.0 km/s). Returns the measured offset in km/s; raises ``exc`` when it
    exceeds ``tol_kms`` (a fixed wavelength offset is the signature of a missing or
    sign-wrong velocity term, and shows up downstream as an inflated χ²ᵣ)."""
    v, centroid = line_velocity_kms(wave_A, flux, expected_A,
                                    window_A=window_A, absorption=absorption)
    if abs(v) > tol_kms:
        raise exc(
            f"{context}line expected at {expected_A} Å lands at {centroid:.3f} Å "
            f"({v:+.2f} km/s, tol ±{tol_kms} km/s). A fixed offset means a missing or "
            f"sign-wrong velocity term (BERV or systemic RV) — verify, do not assume "
            f"(RYA-481 §B.3; RYA-272 sign-check, RYA-478 systemic-RV).")
    return v


# ════════════════════════════════════════════════════════════════════════════
# C.  Wavelength scale — native unit (nm/µm→Å) + scale (vacuum→air), CITED
# ════════════════════════════════════════════════════════════════════════════
# The third declared+verified axis (RYA-264), extending the same loader contract.
# The SPIRou/APERO vacuum-nm "vac→air OWED" that the RYA-481 build loud-declared is
# the concrete driver: reading vacuum IR wavelengths against an air line list
# mis-identifies lines by ~4–6 Å (~15–25 resolution elements at R~70k) — the same
# disease on a third axis.
#
# OPERATION ORDER vs. the §B velocity frame (canonical, documented, deterministic):
# wavelength-scale conversion and the velocity shift BOTH transform the wave axis.
# Convert UNIT+SCALE to air-Å in the OBSERVED frame FIRST (the refractive index n(λ)
# is evaluated at the observed wavelength), THEN apply the §B velocity shift. The
# loaders follow this: to_air_angstrom() runs on the raw grid, then ×(1+v/c). The
# effect of the order is sub-mÅ, but the contract must not be ambiguous.

_UNIT_TO_ANGSTROM = {'A': 1.0, 'nm': 10.0, 'um': 1.0e4}
_NATIVE_SCALES = ('air', 'vacuum')


@dataclass(frozen=True)
class WavelengthConvention:
    """One instrument's CITED native wavelength convention — the single source of
    truth for unit + scale. ``citation`` is the instrument pipeline/DRS doc the
    entry was VERIFIED against (an unsourced entry is a defect — never from memory).
    ``band_A`` is the plausible air-Å coverage for the unit-sanity gate (the ×10 /
    ×10000 catch); ``None`` skips the gate for multi-regime instruments."""
    native_unit: str          # 'A' | 'nm' | 'um'
    native_scale: str         # 'air' | 'vacuum'
    citation: str
    band_A: Optional[tuple] = None     # (lo, hi) air-Å, or None
    band_ref: str = ''


# instrument id (as the loader reports it) -> CITED convention. Unknown -> RAISE.
# VERIFIED 2026-06 against each instrument's DRS/pipeline doc (sources inline). The
# RYA-264 starting-hypothesis table was a guess to confirm, NOT authority — ESPRESSO
# (vacuum, not air) and NIRPS (vacuum) corrected from it on verification.
WAVELENGTH_CONVENTION: dict[str, WavelengthConvention] = {
    'HARPS': WavelengthConvention(
        'A', 'air', "ESO HARPS DRS User Manual (3M6-MAN-HAR-33110-0016): wavelengths "
        "in air, Solar-System barycentric rest frame",
        (3780.0, 6910.0), "HARPS 378–691 nm (Mayor et al. 2003 Msngr 114,20)"),
    'HARPS-N': WavelengthConvention(
        'A', 'air', "HARPS-N DRS (TNG; same DRS lineage as HARPS) — air, barycentric",
        (3830.0, 6930.0), "HARPS-N 383–693 nm (Cosentino et al. 2012 SPIE 8446)"),
    'UVES': WavelengthConvention(
        'A', 'air', "ESO UVES Phase-3 SDP: TUCD1=em.wl;obs.atmos ⇒ air (RYA-271 audit "
        "confirmed empirically)",
        (3000.0, 11000.0), "UVES 300–1100 nm (Dekker et al. 2000 SPIE 4008); matches "
        "uves_loader guard band"),
    'ESPRESSO': WavelengthConvention(
        'A', 'vacuum', "ESPRESSO DRS: native VACUUM, barycentric; an air WAVE_MATRIX "
        "extension is also provided — a loader may prefer that air extension (verified "
        "2026-06; corrects the RYA-264 'air' hypothesis)",
        (3780.0, 7880.0), "ESPRESSO 378–788 nm (Pepe et al. 2021 A&A 645,A96)"),
    'NIRPS': WavelengthConvention(
        'A', 'vacuum', "NIRPS-DRS (adapted from the ESPRESSO DRS ⇒ vacuum), BERV-applied; "
        "Mercier et al. 2025 A&A 700,A8 (arXiv:2507.21290). The α Cen S1D ADP carries a "
        "WAVE_AIR column — the loader PREFERS WAVE_AIR where present (RYA-481 NIRPS comment)",
        (9724.0, 19196.0), "NIRPS 972.4–1919.6 nm (Mercier et al. 2025)"),
    'SPIRou': WavelengthConvention(
        'nm', 'vacuum', "APERO (Cook et al. 2022 PASP 134,114509): vacuum nm. The vac→air "
        "conversion was declared OWED by the RYA-481 SPIRou wiring and is delivered here",
        (9558.0, 25158.0), "SPIRou 0.98–2.35 µm; APERO-extracted 955.8–2515.8 nm "
        "(Donati et al. 2020 MNRAS 498,5684); matches spirou_loader coverage"),
    'CRIRES+': WavelengthConvention(
        'nm', 'vacuum', "ESO cr2res pipeline: nm, vacuum (RYA-480 CRIRES nm→Å, "
        "TOPOCENT→bary). INSTRUME header is 'CRIRES'",
        None, "CRIRES+ 0.95–5.3 µm across settings — band per-setting, gate deferred"),
    'STIS': WavelengthConvention(
        'A', 'vacuum', "HST STIS Instrument Handbook: UV/optical wavelengths are vacuum "
        "(RYA-303/426; the FUV-stays-vacuum boundary at 2000 Å)",
        None, "STIS spans FUV→optical (multi-regime) — gate deferred to per-grating"),
    'COS': WavelengthConvention(
        'A', 'vacuum', "HST COS Instrument Handbook: vacuum wavelengths",
        None, "COS FUV/NUV (multi-segment) — gate deferred"),
    'APOGEE': WavelengthConvention(
        'A', 'vacuum', "SDSS APOGEE: vacuum Å (Nidever et al. 2015 AJ 150,173; "
        "Holtzman et al. 2015 AJ 150,148)",
        (15100.0, 17000.0), "APOGEE H-band 1.51–1.70 µm (Wilson et al. 2019 PASP 131,055001)"),
    'iSHELL': WavelengthConvention(
        'um', 'vacuum', "IRTF iSHELL / Spextool: µm, vacuum (Rayner et al. 2016 SPIE 9908; "
        "Cushing et al. 2004 PASP 116,362)",
        None, "iSHELL 1.06–5.3 µm across modes — gate deferred"),
}


def wavelength_convention(instrument: str, *, exc: type = FrameContractError) -> WavelengthConvention:
    """Return the CITED wavelength convention for an instrument, raising on an
    undeclared one rather than guessing unit/scale (RYA-264 §C). No default."""
    if instrument not in WAVELENGTH_CONVENTION:
        raise exc(
            f"[LOUD] wavelength convention undeclared for instrument {instrument!r} — no "
            f"default unit/scale (RYA-264). Add a CITED entry to WAVELENGTH_CONVENTION "
            f"(verify against its DRS doc) before loading. Known: {sorted(WAVELENGTH_CONVENTION)}.")
    return WAVELENGTH_CONVENTION[instrument]


def _band_gate(w: np.ndarray, conv: WavelengthConvention, instrument: str, exc: type) -> None:
    """Unit-sanity gate: catch a ×10 / ×10000 unit mismatch (the RYA-263 zero class).
    Skipped when the convention declares no band."""
    if conv.band_A is None:
        return
    finite = w[np.isfinite(w)]
    if finite.size == 0:
        raise exc(f"[LOUD] {instrument}: all wavelengths non-finite after conversion.")
    lo, hi = conv.band_A
    wmin, wmax = float(finite.min()), float(finite.max())
    if not (lo * 0.9 <= wmin and wmax <= hi * 1.1):
        raise exc(
            f"[LOUD] {instrument} wavelengths {wmin:.1f}–{wmax:.1f} Å outside the plausible "
            f"air band {lo:.0f}–{hi:.0f} Å ({conv.band_ref}) — unit mismatch (nm/µm read as "
            f"Å)? Refusing to proceed on ×10/×10000-off wavelengths (RYA-264 §C; RYA-263).")


def to_air_angstrom(wave, instrument: str, *, header=None, exc: type = FrameContractError) -> np.ndarray:
    """Convert a native instrument wavelength grid to the pipeline frame (air Å) from
    the CITED per-instrument convention: native unit → Å, then vacuum → air (the ONE
    shared Birch & Downs converter), then the unit-sanity band gate. Raises on an
    undeclared instrument or a band-gate failure — never a silent default or a ×10-off
    array. ``header`` is accepted for forward use (e.g. a per-product WAVUNIT/scale
    keyword); the DRS-provided air-column shortcut (NIRPS WAVE_AIR / ESPRESSO air ext)
    is a loader choice expressed via ``WavelengthScale(provided_air=True)``."""
    conv = wavelength_convention(instrument, exc=exc)
    w = np.asarray(wave, dtype=float) * _UNIT_TO_ANGSTROM[conv.native_unit]
    if conv.native_scale == 'vacuum':
        w = vac_to_air(w)                       # single shared SSOT impl (RYA-264/426)
    elif conv.native_scale != 'air':
        raise exc(f"[LOUD] unknown native scale {conv.native_scale!r} for {instrument}.")
    _band_gate(w, conv, instrument, exc)
    return w


@dataclass(frozen=True)
class WavelengthScale:
    """A loader's explicit declaration of the wavelength scale it converted from —
    the §C analogue of ``VelocityFrame``. Stamp ``.declare()`` into
    ``meta['wavelength_scale']`` so every spectrum carries its scale provenance; call
    ``.validate()`` to cross-check the declared native unit/scale against the CITED
    registry (a loader claiming the wrong convention fails loud); use
    ``.to_air_angstrom(wave)`` to perform the conversion the same way every time."""
    instrument: str
    native_unit: str            # 'A' | 'nm' | 'um'
    native_scale: str           # 'air' | 'vacuum'
    provided_air: bool = False  # loader read a DRS-provided air column (no vac→air needed)
    note: str = ''

    def __post_init__(self):
        if self.native_unit not in _UNIT_TO_ANGSTROM:
            raise FrameContractError(
                f"native_unit={self.native_unit!r} must be one of {sorted(_UNIT_TO_ANGSTROM)}.")
        if self.native_scale not in _NATIVE_SCALES:
            raise FrameContractError(
                f"native_scale={self.native_scale!r} must be one of {_NATIVE_SCALES}.")

    def validate(self) -> 'WavelengthScale':
        """Raise if the declared native unit/scale contradicts the CITED registry for
        this instrument (a loader must not silently disagree with the documented
        convention). Returns self so it chains."""
        conv = wavelength_convention(self.instrument)
        if (self.native_unit, self.native_scale) != (conv.native_unit, conv.native_scale):
            raise FrameContractError(
                f"[LOUD] {self.instrument}: loader declares native "
                f"{self.native_unit}/{self.native_scale} but the cited convention is "
                f"{conv.native_unit}/{conv.native_scale} ({conv.citation}). Resolve the "
                f"disagreement — do not load on a contradicted scale (RYA-264 §C).")
        return self

    def to_air_angstrom(self, wave, *, header=None, exc: type = FrameContractError) -> np.ndarray:
        """Convert ``wave`` to air Å. If ``provided_air`` (the loader read a DRS air
        column, e.g. NIRPS WAVE_AIR), only the unit factor + band gate are applied —
        the vacuum→air step is skipped because the DRS already did it."""
        if self.provided_air:
            conv = wavelength_convention(self.instrument, exc=exc)
            w = np.asarray(wave, dtype=float) * _UNIT_TO_ANGSTROM[self.native_unit]
            _band_gate(w, conv, self.instrument, exc)
            return w
        return to_air_angstrom(wave, self.instrument, header=header, exc=exc)

    def declare(self) -> str:
        """One-line human declaration for ``meta`` / docstrings."""
        path = (f"DRS-provided air column (no vac→air)" if self.provided_air
                else (f"{self.native_scale}→air" if self.native_scale == 'vacuum' else 'air (no conversion)'))
        unit = '' if self.native_unit == 'A' else f"{self.native_unit}→Å; "
        return (f"{self.instrument}: native {self.native_unit}/{self.native_scale} → air Å; "
                f"{unit}{path}{' — ' + self.note if self.note else ''}")


def verify_vac_to_air(known_vac_A: float, known_air_A: float, *, tol_A: float = 0.01,
                      context: str = '', exc: type = FrameContractError) -> float:
    """Scale-axis sign-check (the §C analogue of ``verify_line_position``): confirm the
    shared vac→air converter sends a KNOWN vacuum line to its tabulated air rest
    wavelength (e.g. Mg II k 2796.3543 vac → 2795.528 air). Returns the residual in Å;
    raises if it exceeds ``tol_A`` — a wrong sign or wrong formula shows up here as a
    fixed offset, the deterministic data-free counterpart to the BERV sign-check."""
    got = float(vac_to_air(np.array([known_vac_A]))[0])
    resid = got - known_air_A
    if abs(resid) > tol_A:
        raise exc(
            f"{context}vac→air check failed: {known_vac_A} Å(vac) → {got:.4f} Å, expected "
            f"{known_air_A} Å(air) (residual {resid:+.4f} Å > {tol_A} Å). The shared "
            f"converter is wrong or mis-applied (RYA-264 §C).")
    return resid
