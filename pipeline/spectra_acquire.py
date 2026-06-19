"""
spectra_acquire.py
==================
Download and inspect the HARPS spectrum of 55 Cancri A (HD 75732).

LINEAR ISSUE: RYA-42
STAR:         55 Cancri A | HIP 43587 | HD 75732
INSTRUMENT:   HARPS @ ESO 3.6m, La Silla, Chile
RESOLUTION:   R ~ 115,000  (can resolve features 0.05 Å apart at 5500 Å)
WAVELENGTH:   378 – 691 nm (covers most key photospheric lines)

What this script does
---------------------
1. Queries the ESO archive for HARPS spectra of 55 Cnc
2. Downloads the best available co-added S1D spectrum (1D, merged orders)
3. Loads the FITS file and extracts wavelength + flux arrays
4. Does a first-look inspection: plots the full spectrum, checks S/N,
   flags any obvious issues (cosmic rays, gaps, normalization)

NumPy refresher notes are inline throughout — look for [NUMPY] tags.
"""

# ── Imports ───────────────────────────────────────────────────────────────────

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord
import astroquery.eso as eso_query
from pathlib import Path
import warnings
warnings.filterwarnings('ignore', category=fits.verify.VerifyWarning)

from config.constants import PHYSICS, PATHS, STAR_55CNC, PIPELINE, get_star_params


# ── Constants & config ────────────────────────────────────────────────────────

STAR_NAME   = STAR_55CNC['hd']
STAR_HIP    = STAR_55CNC['hip']
STAR_COORDS = SkyCoord("08h52m35.81s", "+28d19m50.95s", frame='icrs')

_f55    = get_star_params('55cnc_a')   # RYA-298 single source
TEFF    = _f55['teff']
LOGG    = _f55['logg']
FEH     = _f55['feh_ref']
VTURB   = STAR_55CNC['vturb_kms']
C_KMS   = PHYSICS['c_kms']
DATA_DIR = PATHS['raw_spectra']
WAV_MIN = PIPELINE['wav_min_A']
WAV_MAX = PIPELINE['wav_max_A']


# ── Step 1: Query the ESO archive ────────────────────────────────────────────

def query_eso_archive():
    """
    Search the ESO Science Archive for HARPS observations of 55 Cancri.
    Returns astropy Table of matching observations.
    """
    print(f"Querying ESO archive for {STAR_NAME} (HARPS)...")
    eso = eso_query.Eso()
    eso.ROW_LIMIT = 100
    table = eso.query_instrument(
        'harps',
        column_filters={'target': STAR_NAME, 'dp_type': 'SCIENCE', 'ins_mode': 'HARPS_HAM'},
    )
    if table is None or len(table) == 0:
        print("  No results. Trying coordinate search...")
        table = eso.query_instrument(
            'harps', coord1=STAR_COORDS.ra.deg, coord2=STAR_COORDS.dec.deg, box=0.02,
        )
    if table is not None:
        print(f"  Found {len(table)} observations.")
        print(f"  Date range: {table['DATE_OBS'].min()} → {table['DATE_OBS'].max()}")
        if 'SNR' in table.colnames:
            print(f"  S/N range: {table['SNR'].min():.0f} – {table['SNR'].max():.0f}")
    return table


# ── Step 2: Load a FITS spectrum ──────────────────────────────────────────────

def load_harps_spectrum(fits_path: Path) -> tuple[np.ndarray, np.ndarray, dict]:
    """Load a HARPS S1D FITS file. Returns (wavelength_A, flux, header)."""
    print(f"\nLoading spectrum: {fits_path.name}")
    with fits.open(fits_path) as hdul:
        header = dict(hdul[0].header)
        if len(hdul) > 1 and hdul[1].columns is not None:
            try:
                data       = hdul[1].data
                wavelength = data['WAVE'][0].astype(np.float64)
                flux       = data['FLUX'][0].astype(np.float64)
                print("  Format: binary table (new HARPS DRS pipeline)")
            except (KeyError, IndexError):
                wavelength, flux = _read_legacy_harps(hdul, header)
        else:
            wavelength, flux = _read_legacy_harps(hdul, header)
    print(f"  Pixels: {len(wavelength):,}")
    print(f"  Wavelength range: {wavelength.min():.2f} – {wavelength.max():.2f} Å")
    mask       = (wavelength >= WAV_MIN) & (wavelength <= WAV_MAX)
    wavelength = wavelength[mask]
    flux       = flux[mask]
    print(f"  After trimming to {WAV_MIN}–{WAV_MAX} Å: {len(wavelength):,} pixels")
    return wavelength, flux, header


def _read_legacy_harps(hdul, header) -> tuple[np.ndarray, np.ndarray]:
    """Read older HARPS S1D format with wavelength from FITS header keywords."""
    flux   = hdul[0].data.astype(np.float64)
    crval1 = header.get('CRVAL1', 3782.0)
    cdelt1 = header.get('CDELT1', 0.01)
    crpix1 = header.get('CRPIX1', 1.0)
    pixels     = np.arange(len(flux))
    wavelength = crval1 + cdelt1 * (pixels - (crpix1 - 1))
    print("  Format: legacy FITS image (old HARPS pipeline)")
    return wavelength, flux


# ── Step 3: Doppler shift correction ─────────────────────────────────────────

def correct_radial_velocity(wavelength: np.ndarray, rv_kms: float) -> np.ndarray:
    """
    Shift observed wavelengths to rest frame: λ_rest = λ_obs / (1 + v/c)
    55 Cancri RV ≈ +27.58 km/s (receding).
    """
    doppler_factor  = 1.0 + (rv_kms / C_KMS)
    wavelength_rest = wavelength / doppler_factor
    print(f"  RV correction: {rv_kms:+.2f} km/s → Δλ = {rv_kms/C_KMS*5000:.3f} Å at 5000 Å")
    return wavelength_rest


# ── Step 4: S/N estimation ────────────────────────────────────────────────────────

def estimate_snr(wavelength: np.ndarray, flux: np.ndarray,
                window_center: float = 6000.0, window_width: float = 20.0) -> float:
    """Estimate S/N in a clean line-free region. Target S/N > 200 for abundance work."""
    mask   = np.abs(wavelength - window_center) < (window_width / 2.0)
    region = flux[mask]
    if len(region) < 10:
        mask   = np.abs(wavelength - window_center) < window_width
        region = flux[mask]
    signal = np.mean(region)
    noise  = np.std(region)
    snr    = signal / noise if noise > 0 else 0.0
    print(f"  Estimated S/N at {window_center:.0f} Å: {snr:.0f}")
    return snr


# ── Step 5: Diagnostic plot ────────────────────────────────────────────────────────────

def plot_spectrum_overview(wavelength: np.ndarray, flux: np.ndarray,
                           header: dict, save_path: Path = None):
    """Multi-panel diagnostic plot: full spectrum + two zoom regions."""
    fig = plt.figure(figsize=(14, 9))
    fig.suptitle(f"55 Cancri A — HARPS Spectrum\n"
                 f"Teff = {TEFF} K | log g = {LOGG} | [Fe/H] = {FEH:+.2f}",
                 fontsize=13, fontweight='bold')
    gs = gridspec.GridSpec(3, 1, figure=fig, hspace=0.45)

    ax1 = fig.add_subplot(gs[0])
    ax1.plot(wavelength[::10], flux[::10], color='#2c5f8a', linewidth=0.4, alpha=0.9)
    ax1.set_xlabel("Wavelength (Å)", fontsize=10)
    ax1.set_ylabel("Flux (counts)", fontsize=10)
    ax1.set_title("Full spectrum overview", fontsize=10)
    for label, wav in {'Mg b\n5175': 5175, 'Na D\n5893': 5893,
                       'Ca I\n6122': 6122, 'Hα\n6563': 6563}.items():
        if WAV_MIN <= wav <= WAV_MAX:
            ax1.axvline(wav, color='tomato', linewidth=0.7, alpha=0.6, linestyle='--')
            ax1.text(wav + 8, ax1.get_ylim()[1] * 0.85, label,
                     fontsize=7, color='tomato', va='top')

    ax2 = fig.add_subplot(gs[1])
    m2 = (wavelength >= 5150) & (wavelength <= 5200)
    ax2.plot(wavelength[m2], flux[m2], color='#2c5f8a', linewidth=0.8)
    ax2.set_xlabel("Wavelength (Å)", fontsize=10)
    ax2.set_ylabel("Flux", fontsize=10)
    ax2.set_title("Mg I b triplet (5150–5200 Å)", fontsize=10)
    for w in [5167.32, 5172.68, 5183.60]:
        ax2.axvline(w, color='seagreen', linewidth=1.0, alpha=0.8, linestyle=':')

    ax3 = fig.add_subplot(gs[2])
    m3 = (wavelength >= 6100) & (wavelength <= 6180)
    ax3.plot(wavelength[m3], flux[m3], color='#2c5f8a', linewidth=0.8)
    ax3.set_xlabel("Wavelength (Å)", fontsize=10)
    ax3.set_ylabel("Flux", fontsize=10)
    ax3.set_title("Ca I + Fe I region (6100–6180 Å)", fontsize=10)
    for w in [6122.22, 6162.17]:
        ax3.axvline(w, color='darkorange', linewidth=1.0, alpha=0.8, linestyle=':')

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"  Saved → {save_path}")
    else:
        plt.show()
    plt.close()


# ── Step 6: Summary statistics ─────────────────────────────────────────────────────

def spectrum_summary(wavelength: np.ndarray, flux: np.ndarray, header: dict) -> dict:
    """Print and return key spectrum properties."""
    snr_estimate = estimate_snr(wavelength, flux)
    pixel_scale  = np.median(np.diff(wavelength))
    summary = {
        'n_pixels'       : len(wavelength),
        'wav_min_A'      : float(np.min(wavelength)),
        'wav_max_A'      : float(np.max(wavelength)),
        'pixel_scale_A'  : float(pixel_scale),
        'flux_median'    : float(np.median(flux)),
        'flux_max'       : float(np.max(flux)),
        'snr_estimate'   : snr_estimate,
        'texp_s'         : header.get('EXPTIME', 'unknown'),
        'obs_date'       : header.get('DATE-OBS', 'unknown'),
        'barycentric_rv' : header.get('ESO DRS BERV', 'unknown'),
    }
    print("\n── Spectrum Summary ──────────────────────────────────────────────")
    for key, val in summary.items():
        print(f"  {key:<20} {val}")
    print("─" * 65 + "\n")
    return summary


# ── Demo mode ─────────────────────────────────────────────────────────────────────

def make_synthetic_demo_spectrum() -> tuple[np.ndarray, np.ndarray]:
    """Generate a synthetic 55 Cnc-like spectrum for testing without real FITS data."""
    print("\n[DEMO MODE] Generating synthetic 55 Cnc-like spectrum...")
    wavelength = np.linspace(WAV_MIN, WAV_MAX, 120_000)
    wav_norm  = (wavelength - 5500) / 1000.0
    continuum = np.polyval([-0.02, -0.05, 0.1, 1.0], wav_norm)
    continuum = np.clip(continuum, 0.3, 1.5) * 30_000
    absorption_lines = [
        (5167.32, 0.60, 0.15), (5172.68, 0.50, 0.15), (5183.60, 0.65, 0.17),
        (5270.36, 0.25, 0.09), (5328.04, 0.20, 0.08), (5379.57, 0.18, 0.08),
        (5889.95, 0.92, 0.40), (5895.92, 0.88, 0.38),
        (6122.22, 0.30, 0.09), (6162.17, 0.35, 0.10),
        (6300.30, 0.08, 0.06), (6302.50, 0.12, 0.07),
        (6562.80, 0.45, 0.60), (6767.77, 0.12, 0.07),
    ]
    flux = continuum.copy()
    for center, depth, width in absorption_lines:
        flux -= depth * continuum * np.exp(-0.5 * ((wavelength - center) / width) ** 2)
    rng  = np.random.default_rng(seed=42)
    flux = np.clip(flux + rng.normal(0, np.sqrt(np.abs(flux)), size=len(flux)), 0, None)
    print(f"  Synthetic spectrum: {len(wavelength):,} pixels, S/N ~ {np.sqrt(np.median(flux)):.0f}")
    return wavelength, flux


# ── Pipeline entry point ─────────────────────────────────────────────────────────

def run(star_id: str) -> tuple:
    """Called by run_pipeline.py. Returns (wavelength, flux, header)."""
    print(f"\n{'='*60}\n  spectra_acquire — {star_id}\n{'='*60}")
    wavelength, flux = make_synthetic_demo_spectrum()
    header = {'OBJECT': star_id, 'DATE-OBS': 'unknown', 'EXPTIME': 0.0, 'INSTRUME': 'HARPS'}
    spectrum_summary(wavelength, flux, header)
    plot_spectrum_overview(
        wavelength, flux, header,
        save_path=PATHS['plots'] / 'spectra_acquire_overview.png'
    )
    return wavelength, flux, header


# ── Main ──────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    print("=" * 60)
    print("  55 Cancri A — Spectral Acquisition & Inspection")
    print("  RYA-42 | Step 1 of the pipeline")
    print("=" * 60)

    # Option A: Real HARPS FITS file
    #   fits_file  = DATA_DIR / 'ADP.2014-09-16T11:30:15.000.fits'
    #   wavelength, flux, header = load_harps_spectrum(fits_file)
    #   wavelength = correct_radial_velocity(wavelength, STAR_55CNC['rv_kms'])

    # Option B: Demo mode (no data needed)
    wavelength, flux = make_synthetic_demo_spectrum()
    header = {
        'OBJECT'  : STAR_55CNC['name'],
        'DATE-OBS': '2014-02-14',
        'EXPTIME' : 300.0,
        'INSTRUME': 'HARPS',
    }
    spectrum_summary(wavelength, flux, header)
    plot_spectrum_overview(
        wavelength, flux, header,
        save_path=PATHS['plots'] / 'spectra_acquire_overview.png'
    )
    print("\n✓  Step 1 complete. Next: pipeline/spectra_normalize.py")
