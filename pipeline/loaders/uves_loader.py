"""
pipeline/loaders/uves_loader.py
================================
Loader for ESO UVES Phase-3 IDP spectra (Procyon; RYA-271 audit → GO WITH CAVEATS).

File structure (ESO SDP 1D):
  PRIMARY    — observation metadata header
  SPECTRUM   — single-row BINTABLE, each cell a full array:
               WAVE (air Å), FLUX, ERR, QUAL, SNR, CONTRIB (+ *_REDUCED variants)

Key audit facts (RYA-271, confirmed empirically):
  - SPECSYS = TOPOCENT on every valid IDP → the loader MUST apply BERV.
    Verified via Hα 6562.79 vs per-epoch BERV (19–28 km/s discrimination margin).
  - Wavelengths are AIR Angstrom (TUCD1 = em.wl;obs.atmos).
  - 4 Gaia-ESO FLAMES-UVES files (SPECSYS=HELIOCEN, nm, R≈48.5k, OBSTECH=MOS) are
    quarantined — they must NEVER reach this loader. The guards below reject them
    loudly (frame / units) so a stray one can never be double-BERV-corrected.
  - UVES IDPs are NOT telluric-corrected → telluric_corrected=False; the per-epoch
    O I telluric verdict lives in the registry (RYA-272 step 3), enforced downstream.

BERV correction (the classic sign bug — guarded by the smoke test):
  v_corr = SkyCoord(RA,DEC).radial_velocity_correction('barycentric', MJD-OBS @ Paranal)
  λ_bary = λ_topo × (1 + v_corr / c)

References:
  UVES: Dekker et al. 2000, SPIE 4008, 534.  ESO Phase-3 SDP standard.
  Linear: RYA-271 (audit) / RYA-272 (this loader).
"""
from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord, EarthLocation
from astropy.time import Time
import astropy.units as u

from config.constants import PHYSICS
from .base_loader import SpectrumData
from pipeline import frame_object_contract as foc

# Paranal (UT2/Kueyen), hardcoded so the loader never depends on an astropy site
# download (mirrors the SPIRou loader's hardcoded CFHT location).
_PARANAL = EarthLocation(lat=-24.62743941 * u.deg, lon=-70.40498688 * u.deg, height=2669.0 * u.m)
_C_KMS = PHYSICS['c_kms']

_WAVE_MIN_A, _WAVE_MAX_A = 3000.0, 11000.0   # plausible UVES air-wavelength span

# RYA-264 §C: UVES Phase-3 IDP is native AIR Å — declared explicitly (an explicit
# "air, no conversion" is now a stated fact, not an assumption). Cross-checked against
# the cited registry; to_air_angstrom is a no-op (×1) + the unit-sanity band gate.
_UVES_SCALE = foc.WavelengthScale('UVES', native_unit='A', native_scale='air',
                                  note='UVES IDP air Å — explicit no-op (RYA-264)')


class UVESProductError(RuntimeError):
    """Raised when a file is not a loadable topocentric UVES IDP (frame / instrument /
    units / wavelength guard failure). Names the file and the specific failed guard."""


class UVESLoader:
    """Load a single ESO UVES Phase-3 IDP into a barycentric-frame SpectrumData.

        spec = UVESLoader('ADP.....fits').load()
        # spec.wave_A — air Å, barycentric (BERV applied)
        # spec.flux   — flux-calibrated (NaN where QUAL != 0)
        # spec.meta   — provenance incl. berv_kms, setting, R, original SPECSYS
    """

    def __init__(self, filepath: str | Path, *, expected_star: str | None = None):
        self.filepath = Path(filepath)
        # RYA-481 §A: when the caller knows which star this file should be, the
        # loader verifies it by the OBJECT header (never by folder/filename) and
        # fails loud on a mismatch. Default None preserves prior behaviour.
        self.expected_star = expected_star

    def load(self) -> SpectrumData:
        with fits.open(self.filepath) as hdl:
            h0 = hdl[0].header
            self._guard(h0, hdl)
            if self.expected_star is not None:
                foc.assert_object(h0.get('OBJECT'), expected=self.expected_star,
                                  context=f"{self.filepath.name}: ", exc=UVESProductError)
            tab = hdl[1].data
            wave = np.asarray(tab['WAVE'][0], dtype=np.float64)            # air Å, topocentric
            # RYA-264 §C: explicit air-Å declaration (no-op conversion) + unit-sanity
            # gate, in the OBSERVED frame before the velocity shift (canonical order).
            wave = _UVES_SCALE.validate().to_air_angstrom(wave)
            flux = self._first_present(tab, ('FLUX', 'FLUX_REDUCED')).astype(np.float64)
            err = self._first_present(tab, ('ERR', 'ERR_REDUCED'), optional=True)
            qual = self._first_present(tab, ('QUAL', 'QUAL_REDUCED'), optional=True)

        # mask bad pixels (QUAL != 0) → NaN; never crash on NaN
        if qual is not None:
            bad = np.asarray(qual, dtype=np.int64) != 0
            flux = flux.copy()
            flux[bad] = np.nan
            if err is not None:
                err = err.copy().astype(np.float64)
                err[bad] = np.nan

        berv = self._compute_berv(h0)
        wave_bary = wave * (1.0 + berv / _C_KMS)

        return SpectrumData(
            wave_A=wave_bary,
            flux=flux,
            err=(np.asarray(err, dtype=np.float64) if err is not None else None),
            meta=self._build_meta(h0, berv, wave),
        )

    # ── guards (fail loudly, name the file + the failed guard) ────────────────
    def _guard(self, h0: fits.Header, hdl: fits.HDUList) -> None:
        fn = self.filepath.name
        inst = str(h0.get('INSTRUME', ''))
        if not inst.upper().startswith('UVES'):
            raise UVESProductError(f"{fn}: INSTRUME={inst!r} is not UVES — wrong instrument.")
        specsys = str(h0.get('SPECSYS', ''))
        if specsys != 'TOPOCENT':
            raise UVESProductError(
                f"{fn}: SPECSYS={specsys!r}, expected TOPOCENT. A non-topocentric product "
                f"(e.g. a Gaia-ESO HELIOCEN FLAMES-UVES file) must NOT be loaded — applying "
                f"BERV would double-correct. Quarantine it (RYA-271/272).")
        # units guard: BINTABLE wavelength must be Angstrom (GES files are nm)
        unit = str(hdl[1].header.get('TUNIT1', '')).lower()
        if unit and unit not in ('angstrom', 'angstroms', 'a', 'aa'):
            raise UVESProductError(
                f"{fn}: wavelength TUNIT1={unit!r}, expected angstrom (GES products are nm).")
        wave0 = np.asarray(hdl[1].data['WAVE'][0], dtype=np.float64)
        lo, hi = float(np.nanmin(wave0)), float(np.nanmax(wave0))
        if not (_WAVE_MIN_A <= lo and hi <= _WAVE_MAX_A and hi > lo):
            raise UVESProductError(
                f"{fn}: wavelength span [{lo:.1f}, {hi:.1f}] Å outside the plausible UVES "
                f"air range [{_WAVE_MIN_A:.0f}, {_WAVE_MAX_A:.0f}] — units/product check failed.")

    @staticmethod
    def _first_present(tab, names, optional=False):
        for n in names:
            if n in tab.columns.names:
                return np.asarray(tab[n][0])
        if optional:
            return None
        raise UVESProductError(f"BINTABLE missing all of {names}")

    def _compute_berv(self, h0: fits.Header) -> float:
        try:
            target = SkyCoord(ra=float(h0['RA']) * u.deg, dec=float(h0['DEC']) * u.deg, frame='icrs')
            obstime = Time(float(h0['MJD-OBS']), format='mjd', scale='utc', location=_PARANAL)
            return float(target.radial_velocity_correction(kind='barycentric', obstime=obstime)
                         .to(u.km / u.s).value)
        except Exception as exc:
            warnings.warn(f"BERV computation failed for {self.filepath.name}: {exc}. "
                          "Wavelengths left topocentric — do NOT co-add or fit lines.",
                          RuntimeWarning, stacklevel=3)
            return 0.0

    def _setting(self, h0: fits.Header) -> str:
        path = str(h0.get('HIERARCH ESO INS PATH', '?'))
        wlen = h0.get('HIERARCH ESO INS GRAT2 WLEN', h0.get('HIERARCH ESO INS GRAT1 WLEN'))
        return f"{path}{int(wlen)}" if wlen else path

    def _build_meta(self, h0: fits.Header, berv: float, wave_topo: np.ndarray) -> dict:
        # RYA-481 §B: declare + self-check the velocity frame. TOPOCENT IDP → BERV
        # applied here, systemic RV NOT (rest-frame shift is the fit's job, RYA-478).
        # validate() raises on the double-BERV / under-correction contradiction.
        vframe = foc.VelocityFrame(
            specsys=str(h0.get('SPECSYS', 'TOPOCENT')), berv_applied=True,
            systemic_rv_applied=False, wave_units='air', wave_scale='angstrom',
            note='UVES IDP — BERV applied by loader (RYA-272)').validate()
        return {
            'instrument'         : 'UVES',
            'object'             : h0.get('OBJECT', 'UNKNOWN'),
            'date_obs'           : h0.get('DATE-OBS', 'UNKNOWN'),
            'mjd_obs'            : float(h0.get('MJD-OBS', 0.0)),
            'exptime_s'          : float(h0.get('EXPTIME', 0.0)),
            'snr_summary'        : float(h0.get('SNR', -1.0)),
            'resolution_R'       : float(h0.get('SPEC_RES', -1.0)),
            'setting'            : self._setting(h0),
            'program_id'         : h0.get('HIERARCH ESO OBS PROG ID', 'UNKNOWN'),
            'airmass'            : float(h0.get('HIERARCH ESO TEL AIRM START', -1.0)),
            'fluxcal'            : h0.get('HIERARCH ESO PRO SCIENCE', h0.get('FLUXCAL', 'UNKNOWN')),
            # reduction software (the base contract's 'apero_version' slot, generalised)
            'apero_version'      : h0.get('HIERARCH ESO PRO REC1 PIPE ID', 'uves'),
            'berv_kms'           : berv,
            'frame'              : 'barycentric (BERV applied)',
            'velocity_frame'     : vframe.declare(),      # RYA-481 §B declaration
            'wavelength_scale'   : _UVES_SCALE.declare(),  # RYA-264 §C declaration
            'object_canonical'   : foc.resolve_object(h0.get('OBJECT')).canonical,
            'specsys_original'   : h0.get('SPECSYS', 'UNKNOWN'),
            'wave_units_raw'     : 'angstrom',
            'wave_coverage_topo_A': (round(float(np.nanmin(wave_topo)), 1),
                                     round(float(np.nanmax(wave_topo)), 1)),
            'telluric_corrected' : False,   # UVES IDPs are NOT telluric-corrected
            'wapiti_applied'     : False,   # SPIRou-specific; N/A for UVES
            'filepath'           : str(self.filepath),
        }
