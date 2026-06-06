"""
pipeline/loaders/spirou_loader.py
==================================
Loader for SPIRou APERO-reduced _t.fits files (telluric-corrected).

File structure (12 HDUs):
  PRIMARY  — observation metadata header
  FluxAB   — telluric-corrected science flux, fiber A+B, (49, 4088)
  WaveAB   — per-order wavelength solution in nm, (49, 4088)
  BlazeAB  — blaze function, (49, 4088)
  Recon    — APERO telluric model (transmission 0–1), (49, 4088)
  OHLine   — OH sky emission model, (49, 4088)
  FluxA/WaveA/BlazeA — fiber A only (not used)
  FluxB/WaveB/BlazeB — fiber B only (not used)

Key header facts (confirmed 2026-06-06, APERO 0.7.295, RYA-183):
  - Wavelengths in nm (no CUNIT1 keyword — inferred from 955–2515 nm range)
  - BERV = NaN in header — observer-frame wavelengths, correction applied here
  - Wapiti NOT applied to this dataset — documented as known limitation
  - NaN fraction ~25% in FluxAB — normal for deep telluric bands + order gaps
  - QSOGRADE = 1 confirmed on sample files

BERV correction:
  astropy.coordinates computes correction from MJD-OBS + RA_DEG/DEC_DEG in header.
  Applied as: wave_bary = wave_obs * (1 + berv_kms / c_kms)

Telluric rule (PERMANENT — no exceptions):
  Raises TelluricNotCorrectedError if file is not _t.fits.

References:
  APERO: Cook et al. 2022, PASP 134, 1041
  Wapiti: Cook et al. 2022, A&A — NOT applied to this dataset
  SPIRou: Donati et al. 2020, MNRAS 498, 5684

Linear: RYA-183 / RYA-184
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


_CFHT_LOCATION = EarthLocation(
    lat=19.825252  * u.deg,
    lon=-155.468876 * u.deg,
    height=4204.0  * u.m,
)

_C_KMS = PHYSICS['c_kms']


class TelluricNotCorrectedError(RuntimeError):
    """Raised when a non-telluric-corrected SPIRou file is passed to this loader."""
    pass


class SPIRouLoader:
    """
    Load a single SPIRou APERO _t.fits file into a SpectrumData object.

    Parameters
    ----------
    filepath : str or Path
    fiber : str — 'AB' (default), 'A', or 'B'

    Usage
    -----
        spec = SPIRouLoader('2374860t.fits').load()
        # spec.wave_A  — (49, 4088) Å, barycentric
        # spec.flux    — (49, 4088) telluric-corrected
        # spec.meta    — provenance dict
    """

    def __init__(self, filepath: str | Path, fiber: str = 'AB'):
        self.filepath = Path(filepath)
        if fiber not in ('AB', 'A', 'B'):
            raise ValueError(f"fiber must be 'AB', 'A', or 'B', got {fiber!r}")
        self.fiber = fiber

    def load(self) -> SpectrumData:
        with fits.open(self.filepath) as hdl:
            self._check_telluric_product(hdl)
            h0      = hdl[0].header
            flux    = hdl[f'Flux{self.fiber}'].data.astype(np.float64)
            wave_nm = hdl[f'Wave{self.fiber}'].data.astype(np.float64)

        # nm → Å (no CUNIT1 in header — confirmed nm from value range 955–2515)
        wave_A_obs = wave_nm * 10.0

        # Barycentric correction
        berv_kms   = self._compute_berv(h0)
        wave_A_bary = wave_A_obs * (1.0 + berv_kms / _C_KMS)

        # Replace zeros with NaN (APERO uses 0 for masked pixels in some orders)
        flux[flux == 0] = np.nan

        return SpectrumData(
            wave_A = wave_A_bary,
            flux   = flux,
            err    = None,
            meta   = self._build_meta(h0, berv_kms),
        )

    def _check_telluric_product(self, hdl: fits.HDUList) -> None:
        """Enforce telluric rule. Raises TelluricNotCorrectedError on non-_t.fits."""
        h0       = hdl[0].header
        filename = self.filepath.name
        drsoutid = h0.get('DRSOUTID', '')

        if not filename.endswith('t.fits') and 'TELLU' not in str(drsoutid).upper():
            raise TelluricNotCorrectedError(
                f"{filename} does not appear to be a telluric-corrected SPIRou "
                f"product (_t.fits). DRSOUTID={drsoutid!r}. "
                "Telluric rule is absolute — use APERO _t.fits only."
            )

        if 'WAPITI' not in h0:
            warnings.warn(
                f"{filename}: Wapiti systematics correction not applied "
                "(no WAPITI keyword). Known limitation of the Q57 dataset. "
                "Acceptable for EW abundance work — document in data report.",
                UserWarning,
                stacklevel=3,
            )

    def _compute_berv(self, h0: fits.Header) -> float:
        """Compute barycentric correction (km/s) from header coords + MJD."""
        try:
            ra_deg  = float(h0['RA_DEG'])
            dec_deg = float(h0['DEC_DEG'])
            mjd_mid = float(h0.get('MJDMID', h0['MJD-OBS']))

            target  = SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg, frame='icrs')
            # location set on obstime only — astropy raises if also passed to
            # radial_velocity_correction when obstime already carries a location.
            obstime = Time(mjd_mid, format='mjd', scale='utc', location=_CFHT_LOCATION)

            berv = target.radial_velocity_correction(
                kind='barycentric', obstime=obstime,
            ).to(u.km / u.s).value

        except Exception as exc:
            warnings.warn(
                f"BERV computation failed for {self.filepath.name}: {exc}. "
                "Wavelengths will be in observer frame. "
                "Do NOT co-add or fit lines without resolving this.",
                RuntimeWarning,
                stacklevel=3,
            )
            berv = 0.0

        return float(berv)

    def _build_meta(self, h0: fits.Header, berv_kms: float) -> dict:
        return {
            'instrument'         : 'SPIRou',
            'object'             : h0.get('OBJECT', 'UNKNOWN'),
            'date_obs'           : h0.get('DATE-OBS', 'UNKNOWN'),
            'exptime_s'          : float(h0.get('EXPTIME', 0.0)),
            'snr_summary'        : float(h0.get('SPEMSNR', -1.0)),
            'apero_version'      : h0.get('VERSION', 'UNKNOWN'),
            'apero_drspdate'     : h0.get('DRSPDATE', 'UNKNOWN'),
            'apero_drsvdate'     : h0.get('DRSVDATE', 'UNKNOWN'),
            'berv_kms'           : berv_kms,
            'wave_units_raw'     : 'nm',
            'wave_coverage_A'    : (9558.0, 25158.0),
            'n_orders'           : 49,
            'n_pixels_per_order' : 4088,
            'telluric_corrected' : True,
            'wapiti_applied'     : 'WAPITI' in h0,
            'qso_grade'          : h0.get('QSOGRADE', -1),
            'pi_name'            : h0.get('PI_NAME', 'UNKNOWN'),
            'run_id'             : h0.get('RUNID', 'UNKNOWN'),
            'program_id'         : h0.get('QRUNID', 'UNKNOWN'),
            'airmass'            : float(h0.get('AIRMASS', -1.0)),
            'mjd_mid'            : float(h0.get('MJDMID', h0.get('MJD-OBS', 0.0))),
            'fiber'              : self.fiber,
            'filepath'           : str(self.filepath),
        }
