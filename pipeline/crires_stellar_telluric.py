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

CRIRES_STELLAR_SETS = {
    'alpha_cen_a_crires': {
        'holding_id': 'alpha_cen_a_crires_plus',
        'dir': _VET / 'Alpha Cen A' / 'CRIRES',
        'claimed_star': 'A',
        'id_gate': 'acen_ab',          # RYA-423 orbit ID
        'epoch': '2022-04-15',
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
    picked = [d for d in scored if d['absorbed_frac'] >= _MIN_TELLURIC_FRAC][:n_windows]
    if not picked:
        top = f"{scored[0]['absorbed_frac']:.4f}" if scored else "no scorable chip"
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
    err_usable: bool = False
    windows: list = field(default_factory=list)
    fit: dict = field(default_factory=dict)     # BEST_FIT_PARAMETERS, as a dict
    gdas: str = ''
    gdas_md5: str = ''
    rv_kms: float = 0.0
    berv_kms: float = 0.0


def _frame_table(frame: CriresFrame):
    """The whole frame as one wavelength-sorted spectrum, plus the row→segment index that
    puts it back. CRIRES+ orders tile without overlap, so a global sort is well-defined;
    it is asserted rather than assumed."""
    waves, fluxes, errs, idx = [], [], [], []
    for i, s in enumerate(frame.segments):
        ok = np.isfinite(s.wave_A) & np.isfinite(s.flux)
        waves.append(s.wave_A[ok]); fluxes.append(s.flux[ok])
        errs.append(s.err[ok] if len(s.err) == len(s.wave_A) else np.full(ok.sum(), np.nan))
        idx.append(np.full(int(ok.sum()), i))
    w = np.concatenate(waves); f = np.concatenate(fluxes)
    e = np.concatenate(errs); ix = np.concatenate(idx)
    o = np.argsort(w)
    w, f, e, ix = w[o], f[o], e[o], ix[o]
    dup = int(np.count_nonzero(np.diff(w) <= 0))
    if dup:
        raise RuntimeError(
            f"{frame.path.name}: {dup} non-increasing wavelength steps after sorting the "
            f"{len(frame.segments)} chips together — the orders overlap, so a single "
            f"wavelength-sorted spectrum is not a faithful representation of this frame.")
    return w, f, e, ix


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
    hdus, slices = [ph], []
    # seg_index is already wavelength-sorted, so each chip is one contiguous row block.
    edges = np.flatnonzero(np.diff(seg_index) != 0) + 1
    for a, b in zip(np.r_[0, edges], np.r_[edges, len(seg_index)]):
        cols = [fits.Column(name='lambda', format='1D', unit='um',
                            array=wave_A[a:b] / 1.0e4),
                fits.Column(name='flux', format='1D', array=flux[a:b])]
        if err_usable:
            cols.append(fits.Column(name='dflux', format='1D', array=err[a:b]))
        hdus.append(fits.BinTableHDU.from_columns(cols, name=f'CHIP{len(slices) + 1}'))
        slices.append((int(a), int(b)))
    fits.HDUList(hdus).writeto(path, overwrite=True)
    return {'err_usable': err_usable, 'n_chips': len(slices), 'chip_rows': slices}


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
    from astropy.io import fits
    with fits.open(Path(out_dir) / 'BEST_FIT_PARAMETERS.fits') as h:
        d = h[1].data
        out = {}
        for name, val, unc in zip(d['parameter'], d['value'], d['uncertainty']):
            key = str(name).strip().strip('\x00')
            if key:
                out[key] = (float(val), float(unc))
    return out


def _run(cmd, cwd, esorex, log_stem: Path):
    import subprocess
    proc = subprocess.run(cmd, cwd=str(cwd), env=esorex_env(esorex),
                          capture_output=True, text=True)
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
           "--FIT_CONTINUUM=1", f"--CONTINUUM_N={_CONTINUUM_N}",
           # WLC_CONST is a fraction of the chip's HALF WAVELENGTH RANGE, so its physical
           # size depends on how much spectrum a chip covers. Start from zero: cr2res
           # delivers a wavelength-calibrated product, and the default -0.05 asserts a 5%
           # half-range error in it that nothing here has measured.
           "--FIT_WLC=1", "--WLC_N=1", "--WLC_CONST=0.0",
           "--FIT_RES_BOX=FALSE", "--RES_BOX=0.0",
           "--FIT_RES_LORENTZ=FALSE", "--RES_LORENTZ=0.0",
           "--FIT_RES_GAUSS=TRUE", f"--RES_GAUSS={gauss_px:.4f}",
           "--KERNMODE=FALSE", "--KERNFAC=5.0", "--VARKERN=TRUE",
           "--MIRROR_TEMPERATURE_KEYWORD=NONE",
           "--SLIT_WIDTH_KEYWORD=NONE", "--SLIT_WIDTH_VALUE=0.2",
           f"--ELEVATION_VALUE={pos['elevation_m']}",
           f"--LONGITUDE_VALUE={pos['lon']}", f"--LATITUDE_VALUE={pos['lat']}",
           f"--GDAS_PROFILE={gdas}", "model.sof"]
    proc = _run(cmd, in_dir, esorex, work_dir / 'molecfit_model')
    _require_product(proc, out_dir, 'BEST_FIT_MODEL.fits', f"{frame.path.name} model")
    best = _best_fit_dict(out_dir)

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
    lam, mtr = [], []
    with fits.open(ct_out / 'TELLURIC_DATA.fits') as h:
        for hdu in h[1:]:
            d = hdu.data
            if d is None or not getattr(d, 'columns', None):
                continue
            names = d.columns.names
            if 'mtrans' not in names:
                continue
            lam.append(np.asarray(d['mlambda' if 'mlambda' in names else 'lambda'], float))
            mtr.append(np.asarray(d['mtrans'], float))
    if not lam:
        raise RuntimeError(f"{frame.path.name}: TELLURIC_DATA carries no mtrans column")
    lam_um = np.concatenate(lam); mtrans = np.concatenate(mtr)
    o = np.argsort(lam_um)
    mt = np.interp(wave_A, lam_um[o] * 1.0e4, mtrans[o], left=np.nan, right=np.nan)
    ok = np.isfinite(mt) & (mt > _MTRANS_FLOOR) & np.isfinite(flux)
    corr = np.full_like(flux, np.nan)
    corr[ok] = flux[ok] / mt[ok]

    return FrameCorrection(
        frame=frame, wave_A=wave_A, flux_raw=flux, err=err, mtrans=mt, flux_corr=corr,
        seg_index=seg_index, molecules=molecules, windows=windows, fit=best,
        fit_molec=dict(fit_flags), err_usable=err_usable,
        gdas=Path(gdas).name, gdas_md5=_md5(gdas), rv_kms=rv_kms,
        berv_kms=barycentric_correction_kms(frame))


def _require_product(proc, out_dir: Path, name: str, tag: str) -> None:
    """A recipe that returned 0 but wrote no product is NOT a failed fit, and reporting
    it as one sent RYA-939 hunting a problem that did not exist. Name what is wrong."""
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


def _ccf_peak(v, ccf) -> float:
    """Parabolic interpolation about the maximum."""
    i = int(np.nanargmax(ccf))
    if i in (0, len(v) - 1):
        return float(v[i])
    y0, y1, y2 = ccf[i - 1], ccf[i], ccf[i + 1]
    denom = (y0 - 2 * y1 + y2)
    if denom == 0:
        return float(v[i])
    return float(v[i] - 0.5 * (v[i] - v[i - 1]) * (y2 - y0) / denom)


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
    return {'rv_kms': _ccf_peak(v, c), 'n_lines': int(deep.size)}


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
    v, c = _ccf(fc.wave_A, absorption, w_vac, d)
    peak = _ccf_peak(v, c)
    contrast = float((np.nanmax(c) - np.nanmedian(c)) / max(np.nanstd(c), 1e-9))
    rv_topo = peak - zero_point_kms
    return {'rv_topo_kms': rv_topo, 'rv_bary_kms': rv_topo + fc.berv_kms,
            'berv_kms': fc.berv_kms, 'zero_point_kms': zero_point_kms,
            'n_mask_lines': int(w_vac.size), 'ccf_contrast_sigma': contrast,
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


def identify_star(fc: FrameCorrection, rv: dict) -> dict:
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
    from pipeline.acen_orbit import GAMMA, K_A, K_B, SOURCE, predicted_rv, rv_bounds
    mod = _rya423_verdict()
    p = predicted_rv(fc.frame.mjd)
    rv_bary = rv.get('rv_bary_kms', float('nan'))
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
    h['OBJECT'] = str(getattr(f, 'object_name', 'alf Cen A'))
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
            h[card] = (round(fc.fit[key][0], 6), key)
    for mol in fc.molecules:
        k = f'rel_mol_col_{mol}'
        if k in fc.fit:
            h[f'RCOL{mol}'[:8]] = (round(fc.fit[k][0], 6), f'{k} +/- {fc.fit[k][1]:.4g}')
    h['GATEBEF'] = (round(gate['residual_before'], 5), 'D1 residual BEFORE correction')
    h['GATEAFT'] = (round(gate['residual_after'], 5), 'D1 residual AFTER correction')
    h['GATENPX'] = gate['n_px']
    h['GATEPASS'] = (gate['passed'], f"tol={gate.get('tol')}")
    h['BERV'] = (round(fc.berv_kms, 5), 'km/s, add to topocentric RV')
    h['RVTOPO'] = (round(rv.get('rv_topo_kms', float('nan')), 4), 'km/s, measured')
    h['RVBARY'] = (round(rv.get('rv_bary_kms', float('nan')), 4), 'km/s')
    h['RVZP'] = (round(rv.get('zero_point_kms', float('nan')), 4),
                 'km/s telluric wavelength zero-point, subtracted')
    h['STARID'] = (ident['verdict'], 'RYA-423 rule on the measured RV')
    h['IDPRED_A'] = (round(ident['pred_A_kms'], 3), 'km/s predicted alpha Cen A')
    h['IDPRED_B'] = (round(ident['pred_B_kms'], 3), 'km/s predicted alpha Cen B')

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
    out = out_dir / f"alpha_cen_a_crires_{f.wlen_id}_{f.date_obs[:10]}_telluric.fits"
    fits.HDUList([ph, tab]).writeto(out, overwrite=True)
    return out
