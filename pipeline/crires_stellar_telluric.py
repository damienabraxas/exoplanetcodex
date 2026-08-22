"""
pipeline/crires_stellar_telluric.py
===================================
RYA-963 — molecfit telluric correction for DIRECT stellar CRIRES+ targets, starting with
α Cen A (2022-04-15, 6 IDPs, one per Y/J/H/H/K/K setting).

This is the RYA-373 driver generalized, not a second copy of it. Everything about the
telluric fit itself lives in `pipeline.crires_telluric` and is called from here: the IDP
loader (nm→Å at the boundary), the per-band molecule sets, the molecfit invocation, the
real per-night GDAS retrieval, the continuum normalizer. What is genuinely different for
a direct stellar target, and therefore lives here, is:

* **The velocity frame.** Vesta is reflected sunlight and needs the two-leg asteroid
  ephemeris RV (RYA-372). α Cen A is a direct target: its photosphere moves at the
  stellar radial velocity plus the barycentric correction, and there is no ephemeris leg
  at all. The telluric FIT is unchanged — tellurics sit at topocentric rest whatever the
  target is — so what changes is only where the STELLAR lines are, which is what the
  exclusion mask and the residual gate need to know.
* **The stellar-line mask.** RYA-373 hand-listed six K-band solar lines and the CO
  bandhead series because it fitted one order in one band. This set spans 950-2490 nm
  across four bands; six hand-picked wavelengths do not generalize, and hand-listing 24
  bands' worth would be the same defect at larger scale. The mask is therefore DERIVED
  from `PATHS['linelist_solar']` (RYA-381, 1150-25000 Å) — α Cen A is a G2V solar
  analogue, so the solar list is the right photospheric mask — air→vac converted and
  shifted to the frame's stellar velocity.
* **The star-ID gate.** α Cen A and B are 4-8" apart and astrometry cannot split them
  (RYA-952). See `identify_star` below for why this runs AFTER correction, not before.

Order of operations (each step's engine fails loud; nothing is faked):
  1. GDAS gate — a REAL per-night profile or `GDASUnavailable` (no standard atmosphere).
  2. telluric-correct every chip segment in the TOPOCENTRIC frame.
  3. per-frame RV from the corrected spectrum + the barycentric correction.
  4. star-ID against the α Cen AB orbit (RYA-423 `verdict`); anything not confirmed A is
     quarantined, never registered.
  5. D1 residual gate per segment; register the corrected holding with provenance.
"""
from __future__ import annotations

import glob
import hashlib
import json
import os
import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from config.constants import PATHS, codex_path
from pipeline.crires_telluric import (
    NM_TO_A, CriresFrame, CriresSegment, _C_KMS, _CONTINUUM_N, _CRIRES_R, _MTRANS_FLOOR,
    _air_to_vac, continuum_normalize, load_crires_idp, molecules_for_band,
    observatory_position)
from pipeline.telluric.esorex_runtime import (SUPPRESS_PREFIX, esorex_env,
                                              resolve_esorex)

# ── The datasets this driver serves ───────────────────────────────────────────
# A set names WHERE the frames are and WHICH star they are claimed to be. The claim is
# what step 4 tests; it is never assumed.
_VET = codex_path('data.spectra_local') / 'Alpha Centauri (vetted)'
# Corrected products are FITS, and *.fits is gitignored repo-wide, so they live on the
# DATA DRIVE beside the other staged spectra (the RYA-952 pattern for tau Ceti) — not in
# a per-ticket clone that gets deleted. What the repo carries is the SHA-256 manifest.
_SPECTRA = codex_path('data.spectra')

CRIRES_STELLAR_SETS = {
    'alpha_cen_a_crires': {
        'holding_id': 'alpha_cen_a_crires_plus',
        'dir': _VET / 'Alpha Cen A' / 'CRIRES',
        'claimed_star': 'A',
        'id_gate': 'acen_ab',          # RYA-423 orbit ID
        'epoch': '2022-04-15',
        'product_dir': _SPECTRA / 'alpha_cen_a' / 'CRIRESPlus_molecfit',
    },
    # The POSITIVE CONTROL for the α Cen A star-ID (RYA-963). One frame, K2192, the same
    # setting as A's K2192, from the same night 16 minutes later, through the same
    # reduction. Any systematic in the RV method — a wavelength zero-point, a BERV sign,
    # an air/vacuum slip — moves both frames by the same amount, so the pair
    # DISCRIMINATES where either alone only asserts: if the method is sound the two
    # frames must land on OPPOSITE branches of the AB orbit, 6.75 km/s apart at this
    # epoch. Both landing on the same branch means the labels are wrong; both landing
    # off-orbit means the method is.
    'alpha_cen_b_crires': {
        'holding_id': 'alpha_cen_b_crires_plus',
        'dir': _VET / 'Alpha Cen B' / 'CRIRES',
        'claimed_star': 'B',
        'id_gate': 'acen_ab',
        'epoch': '2022-04-15',
        'product_dir': _SPECTRA / 'alpha_cen_b' / 'CRIRESPlus_molecfit',
    },
}


def resolve_set(name: str) -> dict:
    try:
        return CRIRES_STELLAR_SETS[name]
    except KeyError:
        raise ValueError(f"unknown CRIRES+ stellar set {name!r}; declared: "
                         f"{sorted(CRIRES_STELLAR_SETS)}")


# ── Stellar-line mask (DERIVED from the canonical solar list, never hand-listed) ──
_MASK_DEPTH_MIN = 0.05        # central_depth floor for a line worth masking
_MASK_HALF_A = 0.30           # Å half-width of the exclusion around each line at rest


def _solar_list():
    """The canonical solar line list (air Å). Cached on the module — it is 354k rows and
    every frame needs it."""
    global _SOLAR_CACHE
    try:
        return _SOLAR_CACHE
    except NameError:
        pass
    import pandas as pd
    df = pd.read_csv(PATHS['linelist_solar'], low_memory=False)
    df = df[np.isfinite(df['wavelength_air_A']) & np.isfinite(df['central_depth'])]
    _SOLAR_CACHE = df[['wavelength_air_A', 'central_depth']].to_numpy()
    return _SOLAR_CACHE


def stellar_line_intervals(lo_A: float, hi_A: float, rv_kms: float = 0.0,
                           depth_min: float = _MASK_DEPTH_MIN,
                           half_A: float = _MASK_HALF_A) -> list:
    """Merged [lo, hi] intervals (Å, VACUUM, RV-shifted) covering the star's own
    photospheric lines over [lo_A, hi_A].

    molecfit fits continuum + telluric only. A stellar line left inside the fit region is
    absorbed into a telluric column density — the absorber gets credit for photospheric
    iron. Handing molecfit these intervals as WAVE_EXCLUDE is what keeps the fitted
    column an atmospheric quantity. The same intervals mark, for the residual gate, which
    pixels carry the wanted signal and must not be scored as telluric residual."""
    tbl = _solar_list()
    w_air, depth = tbl[:, 0], tbl[:, 1]
    shift = 1.0 + rv_kms / _C_KMS
    w_vac = _air_to_vac(w_air) * shift
    sel = (depth >= depth_min) & (w_vac >= lo_A - half_A) & (w_vac <= hi_A + half_A)
    w = np.sort(w_vac[sel])
    if w.size == 0:
        return []
    lo = w - half_A
    hi = w + half_A
    out = [[float(lo[0]), float(hi[0])]]
    for a, b in zip(lo[1:], hi[1:]):
        if a <= out[-1][1]:
            out[-1][1] = max(out[-1][1], float(b))
        else:
            out.append([float(a), float(b)])
    return out


def _intervals_to_um(intervals) -> list:
    return [(a / 1.0e4, b / 1.0e4) for a, b in intervals]


def stellar_mask(wave_A: np.ndarray, intervals) -> np.ndarray:
    """Boolean: pixel falls inside one of the stellar-line intervals."""
    m = np.zeros(wave_A.shape, dtype=bool)
    for a, b in intervals:
        m |= (wave_A >= a) & (wave_A <= b)
    return m


# ── Velocity frame ────────────────────────────────────────────────────────────
def barycentric_correction_kms(frame: CriresFrame) -> float:
    """The barycentric velocity correction for this exposure, km/s, computed from the
    frame's own RA/DEC/MJD and the observatory position in its header.

    The cr2res IDP carries NO BERV keyword — unlike the HARPS/ESPRESSO DRS products,
    which is exactly the kind of per-instrument difference that has to be checked rather
    than assumed. Sign convention: ADD this to a topocentric radial velocity to get the
    barycentric one (astropy's `radial_velocity_correction`)."""
    from astropy import units as u
    from astropy.coordinates import EarthLocation, SkyCoord
    from astropy.time import Time
    pos = observatory_position(frame)
    loc = EarthLocation.from_geodetic(lon=pos['lon'] * u.deg, lat=pos['lat'] * u.deg,
                                      height=pos['elevation_m'] * u.m)
    sc = SkyCoord(ra=frame.ra * u.deg, dec=frame.dec * u.deg)
    t = Time(frame.mjd, format='mjd', scale='utc')      # location goes on the call, not
    return float(sc.radial_velocity_correction(          # the Time — astropy refuses both
        'barycentric', obstime=t, location=loc).to(u.km / u.s).value)


# ── Fit-window selection (DERIVED from the frame, not hand-listed) ────────────
_WINDOW_A = 60.0            # Å per fit window
_MIN_TELLURIC_FRAC = 0.02   # a window must be at least this absorbed by non-stellar lines
# …and must also carry at least this fraction of the BEST window's information. An
# absolute floor alone is not enough: Y1029's four windows scored 0.267, 0.128, 0.062,
# 0.058, and the bottom two are regions where the "non-stellar absorption" is mostly
# stellar residual the mask did not catch. Fitting there is fitting noise — and it is not
# free. LBLRTM's cost is a WAVENUMBER quantity, so a fixed-Å window is 4.4x wider in
# cm^-1 at 1.05 um (54) than at 2.19 um (12.5); Y1029 was still at chi2 2.4e7 after 193
# mpfit calls (K2192 converged in 34) when esorex reached 13.7 GB and the kernel killed
# it. Weak windows cost the most exactly where they help least.
_RELATIVE_TELLURIC_FLOOR = 0.4


_INFORMATIVE_LO, _INFORMATIVE_HI = 0.15, 0.90


def _absorbed_fraction(wave_A, flux, stellar, lo=_INFORMATIVE_LO, hi=_INFORMATIVE_HI):
    """Fraction of non-stellar pixels whose normalized flux lies in the INFORMATIVE
    depth band — absorbed enough to constrain a column density, not so absorbed that the
    line is saturated.

    Scoring on "below 0.97" alone picks the most absorbed chip, which in an H2O band is
    the SATURATED one: a saturated core carries almost no information about the column
    (dI/dN → 0 there) while dominating chi2, and RYA-373 already had to mask such cores
    out after the fact. Scoring the unsaturated depth band picks the window where the
    column density is actually measurable."""
    cont = continuum_normalize(wave_A, flux)
    ok = np.isfinite(cont) & ~stellar
    if ok.sum() < 50:
        return 0.0
    c = cont[ok]
    return float(np.mean((c > lo) & (c < hi)))


def _score_chips(frame: CriresFrame, rv_kms: float, width_A: float) -> list:
    """Per chip, the best `width_A` window and its informative-absorption score.
    Returns [{score, order, detector, lo_A, hi_A}, …] sorted best-first."""
    scored = []
    for seg in frame.segments:
        ok = np.isfinite(seg.wave_A) & np.isfinite(seg.flux) & (seg.flux > 0)
        if ok.sum() < 500:
            continue
        w, f = seg.wave_A[ok], seg.flux[ok]
        iv = stellar_line_intervals(float(w.min()), float(w.max()), rv_kms)
        stellar = stellar_mask(w, iv)
        best, best_lohi = 0.0, None
        lo, step = float(w.min()) + 2.0, width_A / 2.0
        while lo + width_A <= float(w.max()) - 2.0:
            m = (w >= lo) & (w <= lo + width_A)
            if m.sum() >= 200:
                frac = _absorbed_fraction(w[m], f[m], stellar[m])
                if frac > best:
                    best, best_lohi = frac, (lo, lo + width_A)
            lo += step
        if best_lohi is not None:
            scored.append({'absorbed_frac': best, 'order': seg.order,
                           'detector': seg.detector,
                           'lo_A': best_lohi[0], 'hi_A': best_lohi[1]})
    scored.sort(key=lambda d: -d['absorbed_frac'])
    return scored


def select_fit_windows(frame: CriresFrame, rv_kms: float = 0.0, n_windows: int = 4,
                       width_A: float = _WINDOW_A) -> list:
    """Choose the WAVE_INCLUDE windows molecfit fits this frame's atmosphere on.

    Derived, not tabulated. A hand-written window list per band is the thing that does
    not survive the next setting (RYA-373's six K-band wavelengths do not reach Y), and
    a window list is a physics claim that should be checkable against the data in hand.
    So: continuum-normalize each chip, mark the pixels that are absorbed but are NOT at a
    known stellar line, and take the windows where the INFORMATIVE (unsaturated) part of
    that absorption is largest — the definition of "telluric-rich and stellar-clean"
    applied to the frame itself.

    Returns the `n_windows` best windows, one per chip. Raises if none clears
    `_MIN_TELLURIC_FRAC`: a frame with no telluric signal must say so rather than hand
    molecfit an unconstrained atmosphere."""
    scored = _score_chips(frame, rv_kms, width_A)
    picked = [d for d in scored if d['absorbed_frac'] >= _MIN_TELLURIC_FRAC][:n_windows]
    if not picked:
        top = f"{scored[0]['absorbed_frac']:.4f}" if scored else "no scorable chip"
        raise RuntimeError(
            f"{frame.path.name} ({frame.wlen_id}): no fit window reaches the "
            f"{_MIN_TELLURIC_FRAC:.0%} informative non-stellar absorbed fraction (best "
            f"{top}). There is no telluric signal here to constrain a fit; refusing to "
            f"hand molecfit an unconstrained atmosphere.")
    return picked


# ── Window-matched molecules (a molecule with no band in the window is not free) ──
# Significant near-IR telluric absorption bands, vacuum µm (HITRAN band systems; the
# same set molecfit's own documentation lists for the NIR). A molecule is fitted only
# where one of the frame's fit windows lands inside one of its bands: the RYA-963
# prototype fitted CO2 and CO over 2185-2199 nm, where neither has a band, and the
# optimiser duly drove rel_mol_col_CO2 to its 1e-5 floor and left CO pinned at 1.0 with
# zero uncertainty. Two free parameters, no signal, and a floor-pegged column that a
# later reader could mistake for a measurement of a CO2-free atmosphere.
MOLECULE_BANDS_UM = {
    'H2O': ((0.92, 0.99), (1.09, 1.17), (1.30, 1.51), (1.72, 2.00), (2.45, 2.90)),
    'O2':  ((0.75, 0.78), (1.24, 1.30)),
    'CO2': ((1.55, 1.66), (1.93, 2.12), (2.63, 2.92)),
    'CH4': ((1.60, 1.78), (2.15, 2.50)),
    'CO':  ((2.28, 2.50),),
}


def _overlaps(lo_um, hi_um, bands) -> bool:
    return any(hi_um > a and lo_um < b for a, b in bands)


def plan_fit(frame: CriresFrame, rv_kms: float = 0.0, n_windows: int = 4) -> dict:
    """Choose the fit windows AND the molecule table for one frame.

    Two different questions, which the first cut of this driver conflated:

    * **Which molecules must the MODEL contain?** Every molecule with a band anywhere in
      the frame's wavelength range. calctrans evaluates the fitted atmosphere over the
      whole frame, and a molecule left out of MOLECULES is not modelled at all — the
      transmission simply lacks it, and those pixels come out uncorrected while
      everything about the run says "corrected". K2148's windows land at 1.99-2.07 µm,
      so a window-derived molecule list dropped CH4 and CO and would have shipped an
      unmodelled 2.3 µm CH4/CO forest inside a product registered `telluric_applied`.
    * **Which molecules may be FITTED?** Only those with a band inside a fit window; a
      column with no signal in the fit region is a free parameter with nothing to
      determine it (the prototype pegged rel_mol_col_CO2 at its 1e-5 floor). The rest
      are modelled at their profile column (REL_COL=1, FIT_MOLEC=0) — which is the right
      prior for the well-mixed gases (CO2/CH4/CO), the variable one being H2O.

    So: pick the best windows; then, for any frame molecule not yet covered, add the best
    remaining window that lands in one of its bands, so a fittable molecule is not held
    fixed merely because the top-N ranking missed it."""
    width = _WINDOW_A
    scored = _score_chips(frame, rv_kms, width)
    best = scored[0]['absorbed_frac'] if scored else 0.0
    floor = max(_MIN_TELLURIC_FRAC, _RELATIVE_TELLURIC_FLOOR * best)
    picked = [d for d in scored if d['absorbed_frac'] >= floor][:n_windows]
    if not picked:
        top = f"{best:.4f}" if scored else "no scorable chip"
        raise RuntimeError(
            f"{frame.path.name} ({frame.wlen_id}): no fit window reaches the "
            f"{_MIN_TELLURIC_FRAC:.0%} informative non-stellar absorbed fraction (best "
            f"{top}). There is no telluric signal here to constrain a fit; refusing to "
            f"hand molecfit an unconstrained atmosphere.")

    lo_um = float(np.nanmin([np.nanmin(s.wave_A) for s in frame.segments])) / 1.0e4
    hi_um = float(np.nanmax([np.nanmax(s.wave_A) for s in frame.segments])) / 1.0e4
    frame_molecules = [m for m in molecules_for_band(frame.band)
                       if _overlaps(lo_um, hi_um, MOLECULE_BANDS_UM.get(m, ()))]
    if not frame_molecules:
        raise RuntimeError(
            f"{frame.wlen_id}: none of {molecules_for_band(frame.band)} has a band in "
            f"{lo_um:.3f}-{hi_um:.3f} µm — the band's molecule set does not match its "
            f"own wavelength range.")

    def covered(mol, wins):
        bands = MOLECULE_BANDS_UM.get(mol, ())
        return any(_overlaps(w['lo_A'] / 1e4, w['hi_A'] / 1e4, bands) for w in wins)

    for mol in frame_molecules:
        if covered(mol, picked):
            continue
        bands = MOLECULE_BANDS_UM.get(mol, ())
        for cand in scored:
            # A coverage window is exempt from the RELATIVE floor: it is added because a
            # molecule has no other window at all, so "weaker than the best" is beside
            # the point — the alternative is not fitting that molecule.
            if cand in picked or cand['absorbed_frac'] < _MIN_TELLURIC_FRAC:
                continue
            if _overlaps(cand['lo_A'] / 1e4, cand['hi_A'] / 1e4, bands):
                picked.append(dict(cand, added_for=mol))
                break
    picked.sort(key=lambda d: d['lo_A'])
    fit_flags = {m: bool(covered(m, picked)) for m in frame_molecules}
    return {'windows': picked, 'molecules': tuple(frame_molecules),
            'fit_molec': fit_flags,
            'frame_range_um': (lo_um, hi_um)}


# ── Per-frame molecfit: model (on windows) → calctrans (over the whole frame) ──
@dataclass
class FrameCorrection:
    frame: CriresFrame
    wave_A: np.ndarray                 # full frame, sorted by wavelength
    flux_raw: np.ndarray
    err: np.ndarray
    mtrans: np.ndarray                 # molecfit transmission, same grid
    flux_corr: np.ndarray              # flux_raw / mtrans, deep cores masked
    seg_index: np.ndarray              # row -> index into frame.segments
    molecules: tuple = ()
    fit_molec: dict = field(default_factory=dict)
    moved: dict = field(default_factory=dict)
    err_usable: bool = False
    windows: list = field(default_factory=list)
    fit: dict = field(default_factory=dict)     # BEST_FIT_PARAMETERS, as a dict
    gdas: str = ''
    gdas_md5: str = ''
    rv_kms: float = 0.0
    berv_kms: float = 0.0


def _frame_table(frame: CriresFrame):
    """The frame as one array of chips laid end to end in wavelength order, plus the
    row→segment index that puts it back.

    Chips are ordered by their starting wavelength and each keeps its own ascending
    pixel order; the concatenation is NOT globally sorted, because **CRIRES+ echelle
    orders overlap**. Y1029's ord9/det3 [9629.76, 9687.52] and ord8/det1 [9659.52,
    9724.42] share 28 Å — normal cross-dispersed behaviour, increasing toward the blue,
    and simply absent from the K band that RYA-373 worked in. A global sort would
    interleave two chips' pixels, and the transmission model would then be mapped back
    across the wrong chip boundary in the overlap.

    What IS required, and asserted, is that each chip's own wavelengths increase
    monotonically — that is what makes a chip one contiguous row block, which the FITS
    extensions and the per-chip mtrans mapping both depend on."""
    order = sorted(range(len(frame.segments)),
                   key=lambda i: float(np.nanmin(frame.segments[i].wave_A)))
    waves, fluxes, errs, idx = [], [], [], []
    for i in order:
        s = frame.segments[i]
        ok = np.isfinite(s.wave_A) & np.isfinite(s.flux)
        if ok.sum() < 2:
            continue
        w = s.wave_A[ok]
        f = s.flux[ok]
        e = (s.err[ok] if len(s.err) == len(s.wave_A) else np.full(int(ok.sum()), np.nan))
        if np.any(np.diff(w) <= 0):
            o = np.argsort(w)
            w, f, e = w[o], f[o], e[o]
            if np.any(np.diff(w) <= 0):
                raise RuntimeError(
                    f"{frame.path.name}: chip ord{s.order}/det{s.detector} has duplicate "
                    f"or non-increasing wavelengths even after sorting — it is not a "
                    f"single monotonic chip and cannot be one FITS extension.")
        waves.append(w); fluxes.append(f); errs.append(e)
        idx.append(np.full(len(w), i))
    if not waves:
        raise RuntimeError(f"{frame.path.name}: no chip has usable pixels")
    return (np.concatenate(waves), np.concatenate(fluxes),
            np.concatenate(errs), np.concatenate(idx))

def _md5(path) -> str:
    h = hashlib.md5()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def _write_science(path: Path, frame: CriresFrame, wave_A, flux, err, seg_index) -> dict:
    """SCIENCE for molecfit: ONE FITS EXTENSION PER CRIRES+ CHIP, in wavelength order,
    each carrying lambda[µm]/flux/dflux. Read with `--CHIP_EXTENSIONS=TRUE`, which is
    molecfit's declared form for "chip-specific subranges of a single observation to be
    fitted as a single combined spectrum".

    🔴 This structure is not cosmetic, and getting it wrong produced a run that reported
    success while fitting nothing. Handing the whole 1946-2486 nm setting over as ONE
    chip makes molecfit's `WLC_CONST` default catastrophic: it is documented as a "shift
    relative to half wavelength range", i.e. a FRACTION, so the default -0.05 is -3.4 Å
    on a 134 Å chip (harmless, and duly fitted away) and **-135 Å on a 540 nm one**. The
    model then lands nowhere near the data — `mflux` came out identically 0 across all
    40704 pixels — so chi2 was constant at 2.3229e7 for every mpfit call, the Jacobian
    was zero, and mpfit exited at status 4 with EVERY parameter still at its initial
    value and uncertainty exactly 0: `rel_mol_col_H2O = 1 ± 0`. A "converged" fit that
    never moved, which would have shipped a transmission of 1.0 everywhere as a
    correction. Per-chip extensions keep each WLC polynomial over its own ~140 Å, and
    `WLC_CONST=0` starts from the cr2res wavelength solution rather than assuming a 5%
    half-range error in it.

    Returns {'err_usable': …, 'n_chips': …, 'chip_rows': [(start, stop), …]} — the row
    slices that put the per-chip products back together in the input row order."""
    from astropy.io import fits
    src = fits.getheader(str(frame.path))
    ph = fits.PrimaryHDU()
    for k in src.keys():
        if k.startswith('ESO ') or k in ('MJD-OBS', 'RA', 'DEC', 'UTC', 'INSTRUME'):
            try:
                ph.header[k] = src[k]
            except Exception:
                pass
    err_usable = bool(np.all(np.isfinite(err)) and np.all(err > 0))
    hdus, slices, chip_of, chip_span = [ph], [], {}, []
    # seg_index is already wavelength-sorted, so each chip is one contiguous row block.
    edges = np.flatnonzero(np.diff(seg_index) != 0) + 1
    for a, b in zip(np.r_[0, edges], np.r_[edges, len(seg_index)]):
        cols = [fits.Column(name='lambda', format='1D', unit='um',
                            array=wave_A[a:b] / 1.0e4),
                fits.Column(name='flux', format='1D', array=flux[a:b])]
        if err_usable:
            cols.append(fits.Column(name='dflux', format='1D', array=err[a:b]))
        n = len(slices) + 1
        hdus.append(fits.BinTableHDU.from_columns(cols, name=f'CHIP{n}'))
        slices.append((int(a), int(b)))
        seg = frame.segments[int(seg_index[a])]
        chip_of[(seg.order, seg.detector)] = n
        chip_span.append((float(wave_A[a]), float(wave_A[b - 1])))
    fits.HDUList(hdus).writeto(path, overwrite=True)
    return {'err_usable': err_usable, 'n_chips': len(slices), 'chip_rows': slices,
            'chip_of': chip_of, 'chip_span_A': chip_span}


def _write_ranges(path: Path, intervals_um) -> None:
    from astropy.io import fits
    fits.BinTableHDU.from_columns([
        fits.Column(name='LOWER_LIMIT', format='1D',
                    array=np.array([a for a, _ in intervals_um], float)),
        fits.Column(name='UPPER_LIMIT', format='1D',
                    array=np.array([b for _, b in intervals_um], float))],
    ).writeto(path, overwrite=True)


def _write_molecules(path: Path, molecules, fit_flags=None) -> None:
    """MOLECULES table. FIT_MOLEC must be int32 — astropy's default int64 fails CPL with
    "Type mismatch" (RYA-373). A 0 means modelled at the profile column but not fitted."""
    from astropy.io import fits
    flags = [1 if (fit_flags is None or fit_flags.get(m, True)) else 0 for m in molecules]
    fits.BinTableHDU.from_columns([
        fits.Column(name='LIST_MOLEC', format='4A', array=np.array(list(molecules))),
        fits.Column(name='FIT_MOLEC', format='1J',      # int32 — CPL rejects int64
                    array=np.array(flags, dtype=np.int32)),
        fits.Column(name='REL_COL', format='1D', array=np.ones(len(molecules)))],
    ).writeto(path, overwrite=True)


def _best_fit_dict(out_dir: Path) -> dict:
    """BEST_FIT_PARAMETERS as {name: (value, uncertainty)}. Scans every extension for
    the parameter table rather than assuming extension 1: with CHIP_EXTENSIONS the file
    also carries WAVE_INCLUDE/WAVE_EXCLUDE extensions, and a run that fitted nothing
    leaves the CHIP1 extension empty — which must read as "the fit produced no
    parameters", not as a subscripting crash."""
    from astropy.io import fits
    out = {}
    with fits.open(Path(out_dir) / 'BEST_FIT_PARAMETERS.fits') as h:
        for hdu in h[1:]:
            d = getattr(hdu, 'data', None)
            if d is None or getattr(d, 'columns', None) is None:
                continue
            if 'parameter' not in d.columns.names:
                continue
            for name, val, unc in zip(d['parameter'], d['value'], d['uncertainty']):
                key = str(name).strip().strip('\x00')
                if key:
                    out[key] = (float(val), float(unc))
            break
    if not out:
        raise RuntimeError(
            f"BEST_FIT_PARAMETERS in {out_dir} carries no parameter table — molecfit "
            f"wrote the product but fitted nothing.")
    return out


def assert_fit_moved(best: dict, molecules) -> dict:
    """Refuse a fit that reports success without having moved.

    `converged is not correct`: molecfit exits 0 and writes every product even when mpfit
    took the exit at status 4 with a zero Jacobian, leaving each column at its prior with
    uncertainty exactly 0. That is what a -135 A WLC_CONST produced here, and the
    resulting transmission is 1.0 everywhere — an uncorrected spectrum wearing a
    correction's provenance. The invariant declared in advance: at least one FITTED
    molecular column must carry a non-zero uncertainty, and best_chi2 must be below
    initial_chi2."""
    init = best.get('initial_chi2', (float('nan'), -1))[0]
    bestc = best.get('best_chi2', (float('nan'), -1))[0]
    moved = [m for m in molecules
             if best.get(f'rel_mol_col_{m}', (0.0, 0.0))[1] > 0]
    if not moved:
        raise RuntimeError(
            f"molecfit reported success but no molecular column was actually fitted: "
            f"{ {m: best.get(f'rel_mol_col_{m}') for m in molecules} }. Every "
            f"uncertainty is 0, so mpfit never varied them (initial_chi2={init:.6g}, "
            f"best_chi2={bestc:.6g}). Refusing to emit a transmission that is the prior, "
            f"not a fit.")
    if np.isfinite(init) and np.isfinite(bestc) and not (bestc < init):
        raise RuntimeError(
            f"molecfit reported success but chi2 never improved "
            f"(initial={init:.6g}, best={bestc:.6g}) — the model did not respond to any "
            f"parameter. Refusing to emit this as a correction.")
    return {'fitted_columns': moved, 'initial_chi2': init, 'best_chi2': bestc}


#: Address-space cap for one esorex process, GiB. Sirius has 15 GB and is SHARED with
#: other tickets' runs. On 2026-08-22 the Y1029 model grew to 13.7 GB and the kernel's
#: OOM killer fired — and the kernel chooses the victim, so the process it killed could
#: just as easily have been the concurrent 2.5-hour synthesis another session was running.
#: Capping our own address space converts "the machine kills something" into "this run
#: fails, loudly, with its own name on it". Override with RYA963_MEM_GIB.
_MEM_CAP_GIB = float(os.environ.get('RYA963_MEM_GIB', 8))


def _mem_limiter():
    """preexec_fn that caps the child's address space. Returns None where RLIMIT_AS is
    not meaningful (macOS reports it but does not enforce it usefully for this)."""
    import resource
    cap = int(_MEM_CAP_GIB * (1 << 30))

    def _apply():
        try:
            soft, hard = resource.getrlimit(resource.RLIMIT_AS)
            new_hard = hard if hard != resource.RLIM_INFINITY else cap
            resource.setrlimit(resource.RLIMIT_AS, (min(cap, new_hard), new_hard))
        except (ValueError, OSError):
            pass          # a cap we cannot set is not a reason to refuse to run
    return _apply


def _run(cmd, cwd, esorex, log_stem: Path):
    import subprocess
    kw = {}
    if sys.platform.startswith('linux'):
        kw['preexec_fn'] = _mem_limiter()
    proc = subprocess.run(cmd, cwd=str(cwd), env=esorex_env(esorex),
                          capture_output=True, text=True, **kw)
    log_stem.with_suffix('.stdout.txt').write_text(proc.stdout)
    log_stem.with_suffix('.stderr.txt').write_text(proc.stderr)
    return proc


def correct_frame(frame: CriresFrame, work_dir, rv_kms: float = 0.0,
                  n_windows: int = 4, gdas_path=None) -> FrameCorrection:
    """Telluric-correct ONE CRIRES+ frame end to end, in the TOPOCENTRIC frame.

    `molecfit_model` fits ONE atmosphere per frame on the derived fit windows —
    one exposure saw one sky, so one PWV and one set of column densities is the
    physically coherent thing to fit; per-chip independent fits would let the water
    column disagree with itself across a single 4-second exposure. `molecfit_calctrans`
    then evaluates that atmosphere over every pixel of the frame, which is what makes
    the 20-27 chips affordable: the optimiser runs on ~4 windows, not on 27 chips.

    `rv_kms` is the star's TOPOCENTRIC radial velocity, used only to place the stellar
    exclusion mask. Pass 0 on the first pass and re-run once the RV is measured; the
    telluric fit itself does not depend on it."""
    from astropy.io import fits
    from pipeline.crires_telluric import _resolve_gdas
    work_dir = Path(work_dir)
    in_dir, out_dir, ct_dir = work_dir / 'in', work_dir / 'model', work_dir / 'calctrans'
    for d in (in_dir, out_dir, ct_dir):
        d.mkdir(parents=True, exist_ok=True)
    if frame.specsys.upper() != 'TOPOCENT':
        raise RuntimeError(
            f"{frame.path.name}: telluric fit must run in the TOPOCENTRIC frame, but "
            f"SPECSYS={frame.specsys!r}. The RV shift happens AFTER correction, never "
            f"before (RYA-373).")

    wave_A, flux, err, seg_index = _frame_table(frame)
    plan = plan_fit(frame, rv_kms=rv_kms, n_windows=n_windows)
    windows, molecules, fit_flags = plan['windows'], plan['molecules'], plan['fit_molec']
    exclude_A = []
    for w in windows:
        exclude_A += stellar_line_intervals(w['lo_A'], w['hi_A'], rv_kms)

    sci = _write_science(in_dir / 'science.fits', frame, wave_A, flux, err, seg_index)
    err_usable = sci['err_usable']
    region_chips = ",".join(str(sci['chip_of'][(w['order'], w['detector'])])
                            for w in windows)
    _write_ranges(in_dir / 'wave_include.fits',
                  [(w['lo_A'] / 1.0e4, w['hi_A'] / 1.0e4) for w in windows])
    _write_molecules(in_dir / 'molecules.fits', molecules, fit_flags)
    sof_lines = ['science.fits SCIENCE', 'wave_include.fits WAVE_INCLUDE',
                 'molecules.fits MOLECULES']
    if exclude_A:
        _write_ranges(in_dir / 'wave_exclude.fits', _intervals_to_um(exclude_A))
        sof_lines.insert(2, 'wave_exclude.fits WAVE_EXCLUDE')
    (in_dir / 'model.sof').write_text("\n".join(sof_lines) + "\n")

    gdas = str(gdas_path) if gdas_path else _resolve_gdas(frame, in_dir)
    esorex = resolve_esorex()
    pos = observatory_position(frame)
    px_A = float(np.median(np.diff(wave_A)))
    centre_A = 0.5 * (windows[0]['lo_A'] + windows[0]['hi_A'])
    gauss_px = (centre_A / _CRIRES_R) / px_A

    cmd = [esorex, f"--output-dir={out_dir}", SUPPRESS_PREFIX, "molecfit_model",
           "--COLUMN_LAMBDA=lambda", "--COLUMN_FLUX=flux",
           f"--COLUMN_DFLUX={'dflux' if err_usable else 'NULL'}",
           "--DEFAULT_ERROR=0.01", "--WLG_TO_MICRON=1.0", "--WAVELENGTH_FRAME=VAC",
           "--CHIP_EXTENSIONS=TRUE",     # extensions are chips of ONE observation
           # Which chip each WAVE_INCLUDE range belongs to. Omitting it does not fail:
           # molecfit prints "Assuming that all regions are mapped to Chip 1" and puts
           # every range on the first chip, so three of four ranges then sit outside
           # their chip's wavelength span and the fit has nothing to work with.
           f"--MAP_REGIONS_TO_CHIP={region_chips}",
           "--FIT_CONTINUUM=1", f"--CONTINUUM_N={_CONTINUUM_N}",
           # WLC_CONST is a fraction of the chip's HALF WAVELENGTH RANGE, so its physical
           # size depends on how much spectrum a chip covers. Start from zero: cr2res
           # delivers a wavelength-calibrated product, and the default -0.05 asserts a 5%
           # half-range error in it that nothing here has measured.
           "--FIT_WLC=1", "--WLC_N=1", "--WLC_CONST=0.0",
           "--FIT_RES_BOX=FALSE", "--RES_BOX=0.0",
           "--FIT_RES_LORENTZ=FALSE", "--RES_LORENTZ=0.0",
           "--FIT_RES_GAUSS=TRUE", f"--RES_GAUSS={gauss_px:.4f}",
           # VARKERN was added here without a reason and is not in RYA-931's proven
           # configuration; a kernel allowed to vary across a 27-chip range is extra
           # cost and extra freedom that nothing in this data constrains.
           "--KERNMODE=FALSE", "--KERNFAC=5.0",
           "--MIRROR_TEMPERATURE_KEYWORD=NONE",
           "--SLIT_WIDTH_KEYWORD=NONE", "--SLIT_WIDTH_VALUE=0.2",
           f"--ELEVATION_VALUE={pos['elevation_m']}",
           f"--LONGITUDE_VALUE={pos['lon']}", f"--LATITUDE_VALUE={pos['lat']}",
           f"--GDAS_PROFILE={gdas}", "model.sof"]
    proc = _run(cmd, in_dir, esorex, work_dir / 'molecfit_model')
    _require_product(proc, out_dir, 'BEST_FIT_MODEL.fits', f"{frame.path.name} model")
    best = _best_fit_dict(out_dir)
    moved = assert_fit_moved(best, [m for m in molecules if fit_flags.get(m)])

    # ---- calctrans: the fitted atmosphere, evaluated over every pixel of the frame ----
    for name in ('ATM_PARAMETERS.fits', 'BEST_FIT_PARAMETERS.fits', 'MODEL_MOLECULES.fits'):
        (ct_dir / name).write_bytes((out_dir / name).read_bytes())
    # MODEL_MOLECULES comes from molecfit_model (it carries the FITTED columns) — do
    # not regenerate it here or calctrans would evaluate the prior, not the fit.
    _write_science(ct_dir / 'science.fits', frame, wave_A, flux, err, seg_index)
    from astropy.io import fits as _f
    # The mappings go in as RECIPE PARAMETERS, not as SOF frames. Supplying them as
    # MAPPING_ATMOSPHERIC / MAPPING_CONVOLVE tables — the form the man page lists first —
    # makes calctrans die with "Access beyond boundaries" in cpl_column_get_int followed
    # by "Cannot find atm_parameters. Please check the incorrectly set
    # MAPPING_ATMOSPHERIC parameter", which reads as a wrong index and is really a
    # rejected input form. One science extension maps to atmosphere 1 (one exposure, one
    # atmosphere), so the parameter form says everything the table would have.
    (ct_dir / 'calctrans.sof').write_text(
        "science.fits SCIENCE\n"
        "MODEL_MOLECULES.fits MODEL_MOLECULES\n"
        "ATM_PARAMETERS.fits ATM_PARAMETERS\n"
        "BEST_FIT_PARAMETERS.fits BEST_FIT_PARAMETERS\n")
    ct_out = ct_dir / 'out'
    ct_out.mkdir(parents=True, exist_ok=True)
    # calctrans takes NO column/WLG_TO_MICRON parameters — it reads them from the
    # molecfit_model products. Passing them makes esorex refuse the command line outright
    # ("Command line parameter '--COLUMN_LAMBDA=lambda' not recognized"), which the
    # recipe never even sees.
    mapping = ",".join(['1'] * sci['n_chips'])   # every chip -> the one fitted atmosphere
    cmd2 = [esorex, f"--output-dir={ct_out}", SUPPRESS_PREFIX, "molecfit_calctrans",
            "--CHIP_EXTENSIONS=TRUE",
            f"--MAPPING_ATMOSPHERIC={mapping}", f"--MAPPING_CONVOLVE={mapping}",
            "--USE_ONLY_INPUT_PRIMARY_DATA=FALSE", "calctrans.sof"]
    proc2 = _run(cmd2, ct_dir, esorex, work_dir / 'molecfit_calctrans')
    _require_product(proc2, ct_out, 'TELLURIC_DATA.fits', f"{frame.path.name} calctrans")

    # TELLURIC_DATA carries one extension per chip. Concatenate them and interpolate
    # onto the frame grid; chips do not overlap in wavelength (asserted in _frame_table),
    # so the concatenation is single-valued.
    # MAP BY ROW, and PROVE the mapping against the flux column.
    #
    # Under CHIP_EXTENSIONS=TRUE calctrans COMBINES the chips: TELLURIC_CORR comes back
    # as a single 1D image of exactly as many samples as we handed in, in the order we
    # handed them in (CALCTRANS_CHIPS_COMBINED is the combined input it corresponds to).
    # So no interpolation is needed at all — and interpolating anyway is not merely
    # redundant, it is unsafe wherever two echelle orders overlap, because "the model at
    # this wavelength" is ambiguous there.
    #
    # The row correspondence is then ASSERTED, not assumed, by checking calctrans's own
    # copy of the flux against the flux we supplied. RYA-931 learned this the hard way in
    # the other direction: BEST_FIT_MODEL.lambda is VACUUM, so mapping transmission back
    # BY WAVELENGTH mis-registered by ~1.9 A = 190 pixels. Map by row; check the flux.
    mt_parts = []
    with fits.open(ct_out / 'TELLURIC_CORR.fits') as h:
        for hdu in h[1:]:
            if getattr(hdu, 'data', None) is not None and hdu.data.size:
                mt_parts.append(np.asarray(hdu.data, float).ravel())
    if not mt_parts:
        raise RuntimeError(f"{frame.path.name}: TELLURIC_CORR carries no transmission")
    mt = np.concatenate(mt_parts)
    if mt.size != wave_A.size:
        raise RuntimeError(
            f"{frame.path.name}: calctrans returned {mt.size} transmission samples for "
            f"{wave_A.size} input pixels. The row mapping back to the frame is only "
            f"defined when they correspond one to one.")
    with fits.open(ct_out / 'TELLURIC_DATA.fits') as h:
        col = None
        for hdu in h[1:]:
            d = getattr(hdu, 'data', None)
            if d is not None and getattr(d, 'columns', None) is not None \
                    and 'flux' in d.columns.names:
                col = np.asarray(d['flux'], float)
                break
    if col is None or col.size != flux.size:
        raise RuntimeError(
            f"{frame.path.name}: TELLURIC_DATA carries no flux column of the input "
            f"length, so the row mapping cannot be verified — refusing to assume it.")
    scale = np.nanmedian(np.abs(flux)) or 1.0
    misaligned = int(np.count_nonzero(np.abs(col - flux) > 1e-6 * scale))
    if misaligned:
        raise RuntimeError(
            f"{frame.path.name}: calctrans's flux column differs from the input flux in "
            f"{misaligned} of {flux.size} rows, so its output is NOT in input row order "
            f"and the transmission would be applied to the wrong pixels.")

    ok = np.isfinite(mt) & (mt > _MTRANS_FLOOR) & np.isfinite(flux)
    corr = np.full_like(flux, np.nan)
    corr[ok] = flux[ok] / mt[ok]

    return FrameCorrection(
        frame=frame, wave_A=wave_A, flux_raw=flux, err=err, mtrans=mt, flux_corr=corr,
        seg_index=seg_index, molecules=molecules, windows=windows, fit=best,
        fit_molec=dict(fit_flags), err_usable=err_usable, moved=moved,
        gdas=Path(gdas).name, gdas_md5=_md5(gdas), rv_kms=rv_kms,
        berv_kms=barycentric_correction_kms(frame))


def _require_product(proc, out_dir: Path, name: str, tag: str) -> None:
    """A recipe that returned 0 but wrote no product is NOT a failed fit, and reporting
    it as one sent RYA-939 hunting a problem that did not exist. Name what is wrong."""
    if proc.returncode < 0:
        # A NEGATIVE return code is a SIGNAL, not a fit failure, and reporting it as
        # "esorex FAILED" sends the next reader looking for a problem in the fit. -9 on
        # this box is the OOM killer (confirmed in dmesg for the Y1029 model at 13.7 GB).
        import signal
        try:
            name = signal.Signals(-proc.returncode).name
        except ValueError:
            name = f'signal {-proc.returncode}'
        hint = (' — the kernel OOM killer, or the RLIMIT_AS cap this driver sets '
                f'({_MEM_CAP_GIB:g} GiB). Not a fit failure: the process never got to '
                f'finish. Narrow the fit windows or raise RYA963_MEM_GIB.'
                if -proc.returncode == 9 else '')
        raise RuntimeError(f"{tag}: esorex was KILLED by {name}{hint}")
    if proc.returncode != 0:
        tail = "\n".join(l for l in proc.stdout.splitlines() if 'ERROR' in l)[-2000:]
        raise RuntimeError(f"{tag}: esorex FAILED (rc={proc.returncode}):\n{tail}")
    if not (Path(out_dir) / name).exists():
        produced = sorted(f.name for f in Path(out_dir).glob('*.fits'))
        raise RuntimeError(
            f"{tag}: esorex SUCCEEDED (rc=0) but {name} is missing from {out_dir}. "
            f"Produced: {produced or 'nothing'}. Generic out_NNNN.fits names mean "
            f"--suppress-prefix did not take effect.")


# ── D1 residual gate, generalized off the K band ──────────────────────────────
_GATE_TOL = 0.05


def telluric_residual_gate(fc: FrameCorrection, tol: float = _GATE_TOL) -> dict:
    """RYA-373 Decision 1, generalized: at pixels molecfit says are telluric-DOMINATED
    but not saturated, and that are NOT at one of the star's own lines, does the
    corrected flux return to the local continuum?

    A blanket residual over the whole frame is meaningless for a G2V star — most of the
    frame is photosphere, and the score would be dominated by stellar line depth, which
    is the signal. Restricting to telluric-dominated, stellar-clean pixels is what makes
    the number a statement about the CORRECTION.

    Scored per chip so a single bad order cannot hide inside a frame-wide median, and
    reported both ways. `passed` is the FRAME verdict; `chips` carries the detail."""
    cont = np.full_like(fc.flux_corr, np.nan)
    raw_cont = np.full_like(fc.flux_raw, np.nan)
    for i in range(len(fc.frame.segments)):
        m = fc.seg_index == i
        if m.sum() < 200:
            continue
        cont[m] = continuum_normalize(fc.wave_A[m], fc.flux_corr[m])
        raw_cont[m] = continuum_normalize(fc.wave_A[m], fc.flux_raw[m])
    iv = stellar_line_intervals(float(fc.wave_A.min()), float(fc.wave_A.max()), fc.rv_kms)
    stellar = stellar_mask(fc.wave_A, iv)
    telluric = np.isfinite(fc.mtrans) & (fc.mtrans < 0.90) & (fc.mtrans > _MTRANS_FLOOR)
    sel = telluric & ~stellar & np.isfinite(cont)

    chips = []
    for i, seg in enumerate(fc.frame.segments):
        m = sel & (fc.seg_index == i)
        if m.sum() < 20:
            continue
        chips.append({
            'order': seg.order, 'detector': seg.detector, 'n_px': int(m.sum()),
            'before': float(np.nanmedian(np.abs(1.0 - raw_cont[m]))),
            'after': float(np.nanmedian(np.abs(1.0 - cont[m])))})
    if not chips:
        return {'n_px': int(sel.sum()), 'residual_before': float('nan'),
                'residual_after': float('nan'), 'passed': False, 'chips': [],
                'reason': 'no telluric-dominated, stellar-clean pixel in this frame'}
    before = float(np.nanmedian(np.abs(1.0 - raw_cont[sel])))
    after = float(np.nanmedian(np.abs(1.0 - cont[sel])))
    return {'n_px': int(sel.sum()), 'n_chips': len(chips),
            'residual_before': before, 'residual_after': after,
            'improvement': (before - after) / before if before > 0 else float('nan'),
            'tol': tol, 'passed': bool(after <= tol), 'chips': chips}


# ── Radial velocity from the CORRECTED spectrum ───────────────────────────────
_RV_MASK_DEPTH_MIN = 0.20      # only reasonably strong lines carry CCF weight
_RV_VMAX = 60.0
_RV_DV = 0.25


def _ccf(wave_A, absorption, mask_w, mask_weight, vmax=_RV_VMAX, dv=_RV_DV):
    """Weighted-mask CCF. Returns (velocity grid km/s, CCF). `absorption` is
    1 - normalized flux, so the CCF peaks (maximum) at the line velocity."""
    v = np.arange(-vmax, vmax + dv, dv)
    good = np.isfinite(absorption)
    w, a = wave_A[good], absorption[good]
    out = np.empty_like(v)
    for i, vv in enumerate(v):
        shifted = mask_w * (1.0 + vv / _C_KMS)
        vals = np.interp(shifted, w, a, left=np.nan, right=np.nan)
        ok = np.isfinite(vals)
        out[i] = float(np.sum(vals[ok] * mask_weight[ok]) / max(mask_weight[ok].sum(), 1e-9))
    return v, out


def _ccf_peak(v, ccf) -> dict:
    """Parabolic interpolation about the CCF maximum, WITH the two ways it can fail.

    A peak sitting on the first or last velocity sample is not a measurement — the true
    maximum is outside the search grid and the returned number is just the edge. That is
    exactly what a frame with mostly-NaN transmission produced here: -80.0 km/s, the
    grid edge, reported as if it were an RV. `railed` says so instead. `contrast` is the
    peak's height above the CCF floor in units of its own scatter; a flat CCF has
    nothing to interpolate and no velocity to report."""
    if not np.any(np.isfinite(ccf)):
        return {'rv_kms': float('nan'), 'railed': True, 'contrast': 0.0,
                'reason': 'CCF is entirely non-finite'}
    i = int(np.nanargmax(ccf))
    floor = float(np.nanmedian(ccf))
    scatter = float(np.nanstd(ccf))
    contrast = (float(ccf[i]) - floor) / scatter if scatter > 0 else 0.0
    if i in (0, len(v) - 1):
        return {'rv_kms': float(v[i]), 'railed': True, 'contrast': contrast,
                'reason': f'CCF peak is at the edge of the +/-{abs(v[0]):.0f} km/s '
                          f'search grid — the true maximum is outside it'}
    y0, y1, y2 = ccf[i - 1], ccf[i], ccf[i + 1]
    denom = (y0 - 2 * y1 + y2)
    peak = (float(v[i]) if denom == 0
            else float(v[i] - 0.5 * (v[i] - v[i - 1]) * (y2 - y0) / denom))
    return {'rv_kms': peak, 'railed': False, 'contrast': contrast}


def telluric_zero_point(fc: FrameCorrection) -> dict:
    """The wavelength-solution zero-point, measured against the atmosphere itself.

    Tellurics are at rest in the topocentric frame, so molecfit's transmission model and
    the observed telluric absorption must line up at v=0. Whatever offset the CCF finds
    is the frame's wavelength zero-point error, and it must be subtracted from the
    stellar RV before that RV is compared to an ephemeris. This is also the air-vs-vacuum
    tripwire: cr2res is vacuum, and reading it as air would show up here as ≈ -83 km/s at
    2.2 µm, not as a subtle bias (RYA-373)."""
    telluric = np.isfinite(fc.mtrans) & (fc.mtrans < 0.90) & (fc.mtrans > 0.30)
    if telluric.sum() < 200:
        return {'rv_kms': float('nan'), 'n_lines': 0,
                'reason': 'too few unsaturated telluric pixels'}
    cont = np.full_like(fc.flux_raw, np.nan)
    for i in range(len(fc.frame.segments)):
        m = fc.seg_index == i
        if m.sum() >= 200:
            cont[m] = continuum_normalize(fc.wave_A[m], fc.flux_raw[m])
    obs_abs = 1.0 - cont
    mod_abs = 1.0 - fc.mtrans
    # mask = the model's own telluric line centres
    idx = np.where(telluric)[0]
    deep = idx[(mod_abs[idx] > 0.10)]
    if deep.size < 50:
        return {'rv_kms': float('nan'), 'n_lines': int(deep.size),
                'reason': 'too few deep telluric pixels'}
    v, c = _ccf(fc.wave_A, obs_abs, fc.wave_A[deep], mod_abs[deep], vmax=20.0, dv=0.1)
    pk = _ccf_peak(v, c)
    return dict(pk, n_lines=int(deep.size))


def measure_rv(fc: FrameCorrection, zero_point_kms: float = 0.0) -> dict:
    """Topocentric and barycentric stellar RV from the CORRECTED spectrum.

    Uses the canonical solar line list as the CCF mask (α Cen A is G2V), restricted to
    pixels molecfit says are telluric-clean so a residual telluric line cannot pull the
    peak toward v=0 — which would bias the answer toward the systemic velocity and make
    the star-ID look better than it is."""
    cont = np.full_like(fc.flux_corr, np.nan)
    for i in range(len(fc.frame.segments)):
        m = fc.seg_index == i
        if m.sum() >= 200:
            cont[m] = continuum_normalize(fc.wave_A[m], fc.flux_corr[m])
    clean = np.isfinite(cont) & (np.isfinite(fc.mtrans) & (fc.mtrans > 0.95))
    absorption = np.where(clean, 1.0 - cont, np.nan)

    tbl = _solar_list()
    w_air, depth = tbl[:, 0], tbl[:, 1]
    sel = (depth >= _RV_MASK_DEPTH_MIN)
    w_vac = _air_to_vac(w_air[sel])
    d = depth[sel]
    inrange = (w_vac >= fc.wave_A.min()) & (w_vac <= fc.wave_A.max())
    w_vac, d = w_vac[inrange], d[inrange]
    if w_vac.size < 20:
        return {'rv_topo_kms': float('nan'), 'n_mask_lines': int(w_vac.size),
                'reason': 'too few mask lines in range'}
    if clean.sum() < 500:
        return {'rv_topo_kms': float('nan'), 'railed': True, 'n_clean_px': int(clean.sum()),
                'reason': f'only {int(clean.sum())} telluric-clean pixels — not enough to '
                          f'cross-correlate; an RV from this would be noise'}
    v, c = _ccf(fc.wave_A, absorption, w_vac, d)
    pk = _ccf_peak(v, c)
    rv_topo = pk['rv_kms'] - zero_point_kms
    return {'rv_topo_kms': rv_topo, 'rv_bary_kms': rv_topo + fc.berv_kms,
            'berv_kms': fc.berv_kms, 'zero_point_kms': zero_point_kms,
            'n_mask_lines': int(w_vac.size), 'ccf_contrast_sigma': pk['contrast'],
            'railed': pk['railed'], 'reason': pk.get('reason'),
            'n_clean_px': int(clean.sum())}


# ── Star ID: the A/B split, run on the CORRECTED spectrum ─────────────────────
def _rya423_verdict():
    """RYA-423's `verdict` — imported, never re-implemented. It lives in a script rather
    than a package, so load it by path off the repo root."""
    import importlib.util
    from config.constants import codex_root
    path = Path(codex_root('repo')) / 'scripts' / 'ir_star_id_rya423.py'
    spec = importlib.util.spec_from_file_location('_rya423', path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def identify_star(fc: FrameCorrection, rv: dict, id_gate: str = 'acen_ab') -> dict:
    """Which α Cen component is this frame, by RYA-423's rule?

    **Why this runs AFTER correction, not before.** RYA-423 is explicit that for CRIRES
    it can return only INDETERMINATE: its PRIMARY discriminator is an absolute RV against
    the AB orbit, CRIRES ships no pipeline CCF RV, and "the reduced spectra are
    telluric-dominated -> RV not recoverable per-frame without telluric correction". Its
    own CRIRES branch ends "needs telluric correction + IR templates (downstream gate)".
    Running it on the uncorrected frames, as a gate BEFORE this ticket's work, therefore
    cannot confirm anything — it is a circular dependency, and satisfying it in that
    order would mean quarantining all six frames for want of the very correction they
    are being quarantined from receiving.

    Telluric correction is star-agnostic: molecfit fits Earth's atmosphere at topocentric
    rest and does not care which star is behind it. So correcting first costs nothing and
    breaks the circle — and then RYA-423's PRIMARY becomes available, on this frame's own
    measured RV. Anything the rule does not confirm as A is quarantined, never registered.

    Note also (RYA-423 defect, found here): that script's CRIRES loop globs
    `Alpha Centauri A/CHIRES` — a directory that does not exist under any spectra root,
    with `CRIRES` misspelt. It has therefore always matched ZERO files, so the CRIRES
    branch has never executed on a frame at all. The INDETERMINATE rows it would have
    produced were never produced either."""
    if id_gate != 'acen_ab':
        # RYA-965 will point this same driver at tau Ceti / eps Eri / 55 Cnc, which are
        # singletons: there is no close pair to split, so the id gate there is RYA-964's
        # alias lookup at intake, not an orbit. Dispatch on the set's declared gate
        # rather than growing an alpha-Cen branch inside a generic driver — "ONE recipe;
        # a star that fails it is a finding, not a knob" (RYA-965).
        raise NotImplementedError(
            f"id_gate={id_gate!r} is declared but not implemented here. The alpha Cen AB "
            f"orbit gate is the only close-pair split this driver carries; singleton "
            f"targets resolve their identity through RYA-964 at intake.")
    from pipeline.acen_orbit import GAMMA, K_A, K_B, SOURCE, predicted_rv, rv_bounds
    mod = _rya423_verdict()
    p = predicted_rv(fc.frame.mjd)
    rv_bary = rv.get('rv_bary_kms', float('nan'))
    if rv.get('railed') or not np.isfinite(rv_bary):
        # A railed or absent CCF peak is NOT a velocity. Handing the edge of the search
        # grid to the orbit test would read as "outside the orbit bounds" and return the
        # confident, wrong verdict NOT-ALPHA-CEN on a frame we simply could not measure.
        rv_bary = float('nan')
    # jd/contrast are the NIRPS-calibrated SECONDARY. There is no CRIRES equivalent on
    # the same scale — the CCF contrast measured here is in sigma, RYA-423's threshold is
    # a NIRPS percentage — so they are passed as absent rather than as a wrong-unit
    # number that would silently trip the low-contrast INDETERMINATE branch.
    verdict, evidence = mod.verdict(rv_bary, p['rv_A'], p['rv_B'], None, None)
    lo, hi = rv_bounds()
    return {'verdict': verdict, 'evidence': evidence, 'rv_bary_kms': rv_bary,
            'pred_A_kms': p['rv_A'], 'pred_B_kms': p['rv_B'],
            'delta_AB_kms': p['delta_AB'], 'orbit_bounds_kms': (lo, hi),
            'rv_tol_kms': mod.RV_TOL, 'orbit_source': SOURCE,
            'gamma_kms': GAMMA, 'K_A': K_A, 'K_B': K_B,
            'secondary_available': False,
            'secondary_note': ('RYA-423 SECONDARY (J-band depth / CCF contrast) is '
                               'NIRPS-calibrated and has no CRIRES equivalent on the '
                               'same scale; passed as absent, not as a wrong-unit '
                               'number.')}


# ── Corrected-product output ──────────────────────────────────────────────────
def _num(header, key, value, comment: str, nd: int = 4) -> None:
    """Write a numeric card, or the string UNMEASURED when the value is not finite.

    FITS headers cannot hold NaN at all (astropy refuses outright), and the tempting
    workarounds are both wrong: a sentinel like -999 reads as a measurement, and
    dropping the card makes an unmeasured quantity indistinguishable from one nobody
    thought to record. UNMEASURED says which it is."""
    try:
        v = float(value)
    except (TypeError, ValueError):
        header[key] = ('UNMEASURED', comment)
        return
    header[key] = ((round(v, nd), comment) if np.isfinite(v)
                   else ('UNMEASURED', comment))


def write_corrected(fc: FrameCorrection, out_dir, gate: dict, rv: dict, ident: dict,
                    ticket: str = 'RYA-963') -> Path:
    """Persist ONE telluric-corrected frame: per-chip wavelength, raw flux, corrected
    flux, error and the molecfit transmission model, plus the provenance that makes the
    correction auditable — which GDAS profile (name + md5), which molecules were fitted
    vs held, the fit statistics, the residual gate before/after, and the star ID.

    The uncorrected EXTRACTC IDP remains the archival base; this is a derived product and
    says so."""
    from astropy.io import fits
    out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    f = fc.frame
    ph = fits.PrimaryHDU(); h = ph.header
    h['RYA'] = ticket
    h['ORIGIN'] = ('exoplanetcodex', 'derived product; the EXTRACTC IDP is the base')
    h['BASEFILE'] = (f.path.name, 'uncorrected cr2res EXTRACTC IDP')
    h['BASEMD5'] = _md5(f.path)
    # The raw header label, read from the base frame — never a constant written in here.
    # It is recorded as a LABEL, not an identity: `OBJECT` is not an identifier (RYA-952
    # found tau Ceti hiding under OBJECT='STD'), which is exactly why STARID below
    # carries the measured verdict instead.
    try:
        h['RAWOBJ'] = (str(fits.getheader(str(f.path)).get('OBJECT', '?')).strip(),
                       'raw header label, NOT an identity')
    except Exception:
        h['RAWOBJ'] = ('?', 'raw header label unavailable')
    h['WLEN'] = (f.wlen_id, 'ESO INS WLEN ID')
    h['BAND'] = f.band
    h['MJD-OBS'] = f.mjd
    h['DATE-OBS'] = f.date_obs
    h['SPECSYS'] = ('TOPOCENT', 'telluric fit + product are topocentric')
    h['TELLAPP'] = (True, 'telluric_applied')
    h['TELLENG'] = ('molecfit', 'esorex molecfit_model + molecfit_calctrans')
    h['GDAS'] = (fc.gdas, 'real per-night GDAS profile; no standard atmosphere')
    h['GDASMD5'] = fc.gdas_md5
    h['MOLEC'] = (','.join(fc.molecules), 'modelled molecules')
    h['MOLFIT'] = (','.join(m for m in fc.molecules if fc.fit_molec.get(m)),
                   'fitted; the rest held at the profile column')
    h['NWINDOW'] = (len(fc.windows), 'molecfit fit windows (derived from the frame)')
    h['ERRUSED'] = (fc.err_usable, 'IDP ERR column usable as molecfit dflux')
    for i, w in enumerate(fc.windows, 1):
        h[f'WIN{i}'] = (f"{w['lo_A']:.1f}-{w['hi_A']:.1f} A ord{w['order']}/det{w['detector']}",
                        f"absorbed_frac={w['absorbed_frac']:.3f}")
    for key, card in (('reduced_chi2', 'REDCHI2'), ('rms_rel_to_mean', 'RMSREL'),
                      ('gaussfwhm', 'LSFGAUSS'), ('h2o_col_mm', 'H2OCOLMM')):
        if key in fc.fit:
            _num(h, card, fc.fit[key][0], key, 6)
    for mol in fc.molecules:
        k = f'rel_mol_col_{mol}'
        if k in fc.fit:
            _num(h, f'RCOL{mol}'[:8], fc.fit[k][0],
                 f'{k} +/- {fc.fit[k][1]:.4g}', 6)
    _num(h, 'GATEBEF', gate['residual_before'], 'D1 residual BEFORE correction', 5)
    _num(h, 'GATEAFT', gate['residual_after'], 'D1 residual AFTER correction', 5)
    h['GATENPX'] = gate['n_px']
    h['GATEPASS'] = (gate['passed'], f"tol={gate.get('tol')}")
    _num(h, 'BERV', fc.berv_kms, 'km/s, add to topocentric RV', 5)
    _num(h, 'RVTOPO', rv.get('rv_topo_kms'), 'km/s, measured')
    _num(h, 'RVBARY', rv.get('rv_bary_kms'), 'km/s')
    _num(h, 'RVZP', rv.get('zero_point_kms'), 'km/s telluric zero-point, subtracted')
    h['RVRAILED'] = (bool(rv.get('railed')), 'CCF peak hit the search-grid edge')
    if rv.get('reason'):
        h['RVWHY'] = (str(rv['reason'])[:68].encode('ascii', 'replace').decode(),
                      'why the RV is not a measurement')
    h['STARID'] = (ident['verdict'], 'RYA-423 rule on the measured RV')
    _num(h, 'IDPRED_A', ident['pred_A_kms'], 'km/s predicted alpha Cen A', 3)
    _num(h, 'IDPRED_B', ident['pred_B_kms'], 'km/s predicted alpha Cen B', 3)

    order = np.array([fc.frame.segments[i].order for i in fc.seg_index], dtype=np.int32)
    detec = np.array([fc.frame.segments[i].detector for i in fc.seg_index], dtype=np.int32)
    tab = fits.BinTableHDU.from_columns([
        fits.Column(name='WAVE', format='1D', unit='Angstrom', array=fc.wave_A),
        fits.Column(name='FLUX', format='1D', array=fc.flux_corr),
        fits.Column(name='FLUX_RAW', format='1D', array=fc.flux_raw),
        fits.Column(name='ERR', format='1D', array=fc.err),
        fits.Column(name='MTRANS', format='1D', array=fc.mtrans),
        fits.Column(name='ORDER', format='1J', array=order),
        fits.Column(name='DETEC', format='1J', array=detec)], name='SPECTRUM')
    # A separate MTRANS extension, not only the column. pipeline.telluric_intake derives
    # `telluric_applied` from the product itself rather than taking our word for it, and
    # its rule 2 checks that a transmission extension is NOT all-unity. That gives the
    # frozen-fit failure a SECOND, independent detector: a run whose columns never moved
    # produces transmission 1.0 everywhere, and intake would then refuse to call this
    # holding corrected no matter what the header claims.
    mt = fits.ImageHDU(data=fc.mtrans.astype(np.float64), name='MTRANS')
    mt.header['COMMENT'] = 'molecfit model transmission, aligned row-for-row with SPECTRUM'
    out = out_dir / f"alpha_cen_a_crires_{f.wlen_id}_{f.date_obs[:10]}_telluric.fits"
    fits.HDUList([ph, tab, mt]).writeto(out, overwrite=True)
    return out


# ── Orchestrator ──────────────────────────────────────────────────────────────
def gdas_gate(night: str, site: str = 'paranal', mjds=()) -> dict:
    """STEP 0. Resolve the REAL per-night GDAS profile for every exposure, or raise
    GDASUnavailable. Reported before any molecfit runs, because the one external
    dependency of this whole leg is data availability, not compute."""
    from pipeline.telluric.gdas_fetch import fetch_gdas, nearest_3hourly
    slots, paths = {}, {}
    for mjd in (mjds or ()):
        slot = nearest_3hourly(night, float(mjd))
        slots[float(mjd)] = f"{slot:%Y-%m-%dT%H}"
        paths[f"{slot:%Y-%m-%dT%H}"] = fetch_gdas(site, night=night, mjd=float(mjd))
    if not mjds:
        slot = nearest_3hourly(night)
        slots['night'] = f"{slot:%Y-%m-%dT%H}"
        paths[f"{slot:%Y-%m-%dT%H}"] = fetch_gdas(site, night=night)
    return {'site': site, 'night': night,
            'slots': slots,
            'profiles': {k: str(v) for k, v in paths.items()},
            'md5': {k: _md5(v) for k, v in paths.items()},
            'n_distinct_profiles': len(paths),
            'standard_atmosphere_fallback': False}


def run_set(name: str = 'alpha_cen_a_crires', work_root=None, out_dir=None,
            n_windows: int = 4, limit=None) -> dict:
    """The whole ticket for one declared set: GDAS gate → per-frame telluric correction
    → RV → star ID → residual gate → product. Frames the star-ID gate does not confirm
    are QUARANTINED: their correction is kept (it is real work and star-agnostic) but
    they are not registered under the star's holding."""
    rec = resolve_set(name)
    files = sorted(glob.glob(str(Path(rec['dir']) / '*.fits')))
    if not files:
        raise FileNotFoundError(f"no CRIRES+ IDP under {rec['dir']}")
    all_frames = [load_crires_idp(f) for f in files]
    # A set is ONE EPOCH. The night fixes the GDAS profile, and a frame from a different
    # night needs its own — quietly folding it in would correct it against the wrong
    # night's atmosphere, which is the RYA-373 failure mode wearing a different hat. The
    # α Cen B directory holds two 2025-03-11 'Star S5' frames (RYA-423 quarantine) that
    # this excludes by date rather than by name.
    frames, other_epoch = [], []
    for fr in all_frames:
        (frames if str(fr.date_obs)[:10] == rec['epoch'] else other_epoch).append(fr)
    if not frames:
        raise RuntimeError(
            f"{name}: none of {len(all_frames)} frames under {rec['dir']} is from the "
            f"declared epoch {rec['epoch']} "
            f"(found {sorted({str(f.date_obs)[:10] for f in all_frames})}).")
    frames = frames[:limit]
    gate0 = gdas_gate(rec['epoch'], mjds=[f.mjd for f in frames])

    from config.constants import codex_root
    work_root = Path(work_root) if work_root else Path(codex_root('work')) / 'rya963'
    out_dir = Path(out_dir) if out_dir else Path(rec['product_dir'])
    results, confirmed, quarantined, failed = [], [], [], []
    for fr in frames:
      try:
        # Pass 1 places the stellar mask at the EXPECTED velocity; the RV it then
        # measures is what pass 2 would use. The mask is 0.6 A wide and the expected and
        # measured velocities differ by far less than that, so one pass suffices — but
        # the measured offset is reported so that assumption stays checkable.
        from pipeline.acen_orbit import predicted_rv
        berv = barycentric_correction_kms(fr)
        rv_expected_topo = predicted_rv(fr.mjd)['rv_A'] - berv
        fc = correct_frame(fr, work_root / fr.wlen_id, rv_kms=rv_expected_topo,
                           n_windows=n_windows,
                           gdas_path=gate0['profiles'][gate0['slots'][fr.mjd]])
        zp = telluric_zero_point(fc)
        # A railed zero-point is not a zero-point; correcting the stellar RV by the edge
        # of a search grid would inject a fabricated tens-of-km/s shift.
        zp_use = (zp['rv_kms'] if (not zp.get('railed')
                                   and np.isfinite(zp.get('rv_kms', np.nan)))
                  else 0.0)
        rv = measure_rv(fc, zero_point_kms=zp_use)
        ident = identify_star(fc, rv, id_gate=rec['id_gate'])
        gate = telluric_residual_gate(fc)
        product = write_corrected(fc, out_dir, gate, rv, ident)
        row = {'frame': fr.path.name, 'wlen_id': fr.wlen_id, 'band': fr.band,
               'mjd': fr.mjd, 'date_obs': fr.date_obs, 'snr_hdr': fr.snr,
               'molecules': fc.molecules, 'fitted': fc.moved.get('fitted_columns'),
               'initial_chi2': fc.moved.get('initial_chi2'),
               'best_chi2': fc.moved.get('best_chi2'),
               'reduced_chi2': fc.fit.get('reduced_chi2', (None,))[0],
               'h2o_col_mm': fc.fit.get('h2o_col_mm', (None,))[0],
               'lsf_gauss_px': fc.fit.get('gaussfwhm', (None,))[0],
               'windows': [(round(w['lo_A'], 1), round(w['hi_A'], 1)) for w in fc.windows],
               'gdas': fc.gdas, 'gdas_md5': fc.gdas_md5,
               'berv_kms': fc.berv_kms, 'telluric_zero_point_kms': zp.get('rv_kms'),
               'telluric_zero_point_railed': bool(zp.get('railed')),
               'telluric_zero_point_applied': zp_use,
               'rv_railed': bool(rv.get('railed')),
               'rv_reason': rv.get('reason'),
               'ccf_contrast_sigma': rv.get('ccf_contrast_sigma'),
               'rv_topo_kms': rv.get('rv_topo_kms'), 'rv_bary_kms': rv.get('rv_bary_kms'),
               'rv_expected_topo_kms': rv_expected_topo,
               'star_id': ident['verdict'], 'star_id_evidence': ident['evidence'],
               'pred_A_kms': ident['pred_A_kms'], 'pred_B_kms': ident['pred_B_kms'],
               'gate_before': gate['residual_before'], 'gate_after': gate['residual_after'],
               'gate_n_px': gate['n_px'], 'gate_passed': gate['passed'],
               'product': str(product), 'product_sha256': _sha256(product)}
        results.append(row)
        (confirmed if ident['verdict'] == rec['claimed_star'] else quarantined).append(row)
      except Exception as exc:
        # One frame failing must not cost the other five. A failure is a per-frame
        # FINDING, recorded with its reason and carried through to the report; it is
        # never a silent skip, and the set's verdict names it.
        failed.append({'frame': fr.path.name, 'wlen_id': fr.wlen_id,
                       'error': f"{type(exc).__name__}: {exc}"})
        print(f"  [FAIL] {fr.wlen_id}: {type(exc).__name__}: {exc}")
    return {'set': name, 'holding_id': rec['holding_id'], 'gdas_gate': gate0,
            'epoch': rec['epoch'], 'claimed_star': rec['claimed_star'],
            'other_epoch_frames': [{'frame': f.path.name, 'date_obs': f.date_obs,
                                    'wlen_id': f.wlen_id} for f in other_epoch],
            'frames': results, 'n_confirmed': len(confirmed),
            'n_quarantined': len(quarantined), 'n_failed': len(failed),
            'quarantined': [r['frame'] for r in quarantined],
            'failed': failed,
            'out_dir': str(out_dir)}


def _sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def main(argv=None):
    import argparse
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[2])
    ap.add_argument('--set', default='alpha_cen_a_crires', dest='set_name')
    ap.add_argument('--gdas-only', action='store_true',
                    help='STEP 0 only: report the GDAS verdict and stop.')
    ap.add_argument('--plan-only', action='store_true',
                    help='report the derived fit windows + molecule table, no molecfit.')
    ap.add_argument('--limit', type=int, default=None)
    ap.add_argument('--n-windows', type=int, default=4)
    ap.add_argument('--work-root', default=None)
    ap.add_argument('--out-dir', default=None)
    ap.add_argument('--json', default=None, help='write the run record here')
    a = ap.parse_args(argv)

    rec = resolve_set(a.set_name)
    files = sorted(glob.glob(str(Path(rec['dir']) / '*.fits')))[:a.limit]
    if not files:
        # A set whose frames are not reachable must SAY so. Returning quietly here would
        # report "nothing wrong" on a machine that simply cannot see the spectra — and
        # the CRIRES+ holdings live on Sirius, so that is the normal case on the Mac.
        raise SystemExit(
            f"no CRIRES+ IDP FITS under {rec['dir']} — this set's frames are not "
            f"reachable from here. The spectra root resolves to "
            f"{codex_path('data.spectra_local')}; on the Mac that is a Sirius path "
            f"(the holdings migrated 2026-08-16), so run this on Sirius or point "
            f"CODEX_SPECTRA_LOCAL at a reachable copy.")
    frames = [load_crires_idp(f) for f in files]
    if a.gdas_only:
        g = gdas_gate(rec['epoch'], mjds=[f.mjd for f in frames])
        print(f"STEP 0 GDAS gate: {rec['epoch']} @ {g['site']}")
        for slot, path in g['profiles'].items():
            print(f"  slot {slot}  {Path(path).name}  md5={g['md5'][slot]}")
        print(f"  {len(frames)} frames -> {g['n_distinct_profiles']} distinct profile(s); "
              f"standard-atmosphere fallback: {g['standard_atmosphere_fallback']}")
        return 0
    if a.plan_only:
        from pipeline.acen_orbit import predicted_rv
        for fr in frames:
            rv = predicted_rv(fr.mjd)['rv_A'] - barycentric_correction_kms(fr)
            p = plan_fit(fr, rv_kms=rv, n_windows=a.n_windows)
            print(f"{fr.wlen_id} {p['frame_range_um'][0]:.3f}-{p['frame_range_um'][1]:.3f} um "
                  f"molecules={p['molecules']} fit={p['fit_molec']}")
            for w in p['windows']:
                print(f"    ord{w['order']}/det{w['detector']} "
                      f"{w['lo_A']:.1f}-{w['hi_A']:.1f} A  f={w['absorbed_frac']:.3f}"
                      + (f"  [added for {w['added_for']}]" if 'added_for' in w else ""))
        return 0

    out = run_set(a.set_name, work_root=a.work_root, out_dir=a.out_dir,
                  n_windows=a.n_windows, limit=a.limit)
    for r in out['frames']:
        print(f"{r['wlen_id']:6s} chi2 {r['initial_chi2']:.4g} -> {r['best_chi2']:.4g}  "
              f"PWV {r['h2o_col_mm']:.3f} mm  gate {r['gate_before']:.4f} -> "
              f"{r['gate_after']:.4f} {'PASS' if r['gate_passed'] else 'FAIL'}  "
              f"RV_bary {r['rv_bary_kms']:+.2f}  ID {r['star_id']}")
    for f in out.get('failed', []):
        print(f"{f['wlen_id']:6s} FAILED — {f['error']}")
    print(f"\nconfirmed {out['n_confirmed']} / quarantined {out['n_quarantined']} / "
          f"failed {out.get('n_failed', 0)}")
    if a.json:
        Path(a.json).write_text(json.dumps(out, indent=2, default=str) + "\n")
        print(f"[record] {a.json}")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
