"""
01_acquire.py
=============
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

import numpy as np                        # The core numerical library.
                                          # Everything we do with arrays lives here.

import matplotlib.pyplot as plt           # Plotting. We'll use this constantly.
import matplotlib.gridspec as gridspec    # For multi-panel figure layouts.

from astropy.io import fits               # Read FITS files (the standard astronomy format).
                                          # FITS = Flexible Image Transport System.
                                          # Every spectrum, image, table from a telescope
                                          # comes as a FITS file.

from astropy import units as u            # Physical units (Angstrom, km/s, etc.)
from astropy.coordinates import SkyCoord  # Celestial coordinates (RA, Dec)
import astroquery.eso as eso_query        # Query the ESO archive programmatically.

from pathlib import Path

import warnings
warnings.filterwarnings('ignore', category=fits.verify.VerifyWarning)

from config.constants import PHYSICS, PATHS, STAR_55CNC, PIPELINE


# ── Constants & config ────────────────────────────────────────────────────────

STAR_NAME   = STAR_55CNC['hd']           # Primary identifier in ESO archive
STAR_HIP    = STAR_55CNC['hip']          # Hipparcos catalog ID
STAR_COORDS = SkyCoord("08h52m35.81s", "+28d19m50.95s", frame='icrs')

# Stellar parameters (CHARA interferometry, von Braun et al. 2011)
# Used to set up the model atmosphere in 03_line_fitting.py
TEFF    = STAR_55CNC['teff_K']           # K  — effective temperature
LOGG    = STAR_55CNC['logg']             # cgs — surface gravity
FEH     = STAR_55CNC['feh']              # dex — [Fe/H]
VTURB   = STAR_55CNC['vturb_kms']        # km/s — microturbulence (initial estimate)

C_KMS   = PHYSICS['c_kms']              # km/s — speed of light (NIST CODATA 2022)

DATA_DIR = PATHS['raw_spectra']          # Where downloaded spectra are saved
                                          # (directory auto-created on import)

# Wavelength range of interest (Angstroms)
# HARPS covers 3780–6910 Å across 72 echelle orders
WAV_MIN = PIPELINE['wav_min_A']          # Å  — blue limit
WAV_MAX = PIPELINE['wav_max_A']          # Å  — red limit


# ── Step 1: Query the ESO archive ────────────────────────────────────────────

def query_eso_archive():
    """
    Search the ESO Science Archive for HARPS observations of 55 Cancri.

    The ESO archive is public — no account needed for querying.
    Downloading requires a free ESO user portal account.

    Returns
    -------
    table : astropy Table of matching observations
    """
    print(f"Querying ESO archive for {STAR_NAME} (HARPS)...")

    eso = eso_query.Eso()
    eso.ROW_LIMIT = 100

    table = eso.query_instrument(
        'harps',
        column_filters={
            'target': STAR_NAME,
            'dp_type': 'SCIENCE',
            'ins_mode': 'HARPS_HAM',
        }
    )

    if table is None or len(table) == 0:
        print("  No results. Trying coordinate search...")
        table = eso.query_instrument(
            'harps',
            coord1=STAR_COORDS.ra.deg,
            coord2=STAR_COORDS.dec.deg,
            box=0.02,
        )

    if table is not None:
        print(f"  Found {len(table)} observations.")
        print(f"  Date range: {table['DATE_OBS'].min()} → {table['DATE_OBS'].max()}")
        if 'SNR' in table.colnames:
            print(f"  S/N range: {table['SNR'].min():.0f} – {table['SNR'].max():.0f}")

    return table


# ── Step 2: Load a FITS spectrum ──────────────────────────────────────────────

def load_harps_spectrum(fits_path: Path) -> tuple[np.ndarray, np.ndarray, dict]:
    """
    Load a HARPS 1D merged spectrum (S1D format) from a FITS file.

    HARPS produces several data products per observation:
      - S2D: 2D echelle orders (72 separate stripes) — raw, highest detail
      - S1D: 1D merged & rebinned spectrum — easiest to work with
      - CCF: Cross-Correlation Function — used for radial velocity

    We want S1D for our abundance work: it's already wavelength-calibrated,
    barycentric-corrected, and merged into a single continuous array.

    Parameters
    ----------
    fits_path : Path to the FITS file

    Returns
    -------
    wavelength : np.ndarray  — wavelengths in Angstroms
    flux       : np.ndarray  — normalized flux (roughly, raw counts or e⁻)
    header     : dict        — FITS header metadata (Texp, date, S/N, etc.)
    """
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

    # ── [NUMPY] np.ndarray basics ──────────────────────────────────────────
    # wavelength and flux are 1D NumPy arrays — think of them as Python lists
    # on steroids. Each index i corresponds to one wavelength pixel:
    #   wavelength[i] = the wavelength at pixel i  (in Å)
    #   flux[i]       = the measured photon count at that wavelength
    #
    # The key superpower: array operations work element-wise WITHOUT loops.
    #   flux * 2          → doubles every element (no for loop needed)
    #   flux[flux < 0]    → boolean mask: grab only negative values
    #   np.mean(flux)     → mean of all elements in one call
    # ──────────────────────────────────────────────────────────────────────

    print(f"  Pixels: {len(wavelength):,}")
    print(f"  Wavelength range: {wavelength.min():.2f} – {wavelength.max():.2f} Å")
    print(f"  Flux range: {flux.min():.1f} – {flux.max():.1f}")

    # ── [NUMPY] Boolean indexing ───────────────────────────────────────────
    # (wavelength >= WAV_MIN) creates a boolean array: [False, False, True, ...]
    # We combine two conditions with & (AND) — both must be True.
    # This "boolean masking" filters the array without a for loop.
    # ──────────────────────────────────────────────────────────────────────
    mask       = (wavelength >= WAV_MIN) & (wavelength <= WAV_MAX)
    wavelength = wavelength[mask]
    flux       = flux[mask]

    print(f"  After trimming to {WAV_MIN}–{WAV_MAX} Å: {len(wavelength):,} pixels")

    return wavelength, flux, header


def _read_legacy_harps(hdul, header) -> tuple[np.ndarray, np.ndarray]:
    """
    Read older HARPS S1D FITS format where the spectrum is encoded
    as a 1D array in HDU 0, with wavelength reconstructed from header keywords.
    """
    flux = hdul[0].data.astype(np.float64)

    # FITS linear wavelength solution: WAV = CRVAL1 + CDELT1 * (pixel - CRPIX1)
    crval1 = header.get('CRVAL1', 3782.0)
    cdelt1 = header.get('CDELT1', 0.01)
    crpix1 = header.get('CRPIX1', 1.0)

    # ── [NUMPY] np.arange — building a pixel index array ──────────────────
    # np.arange(n) creates [0, 1, 2, ..., n-1] — like Python's range() but
    # returns a NumPy array, so we can do math on it directly.
    # ──────────────────────────────────────────────────────────────────────
    pixels     = np.arange(len(flux))
    wavelength = crval1 + cdelt1 * (pixels - (crpix1 - 1))

    print("  Format: legacy FITS image (old HARPS pipeline)")
    return wavelength, flux


# ── Step 3: Doppler shift correction ─────────────────────────────────────────

def correct_radial_velocity(wavelength: np.ndarray,
                             rv_kms: float) -> np.ndarray:
    """
    Correct the observed spectrum for the star's radial velocity.

    When 55 Cancri moves toward or away from us, all its spectral lines
    shift in wavelength (Doppler effect). We need to shift them back to
    the rest-frame lab wavelengths so they match our line list.

    The non-relativistic Doppler formula:
        Δλ/λ = v/c
        λ_rest = λ_obs / (1 + v/c)

    55 Cancri's systemic radial velocity ≈ +27.58 km/s (moving away from us).
    HARPS spectra are usually already barycentric-corrected (Earth's motion
    removed), but the stellar RV still needs to be applied.

    Parameters
    ----------
    wavelength : observed wavelength array (Å)
    rv_kms     : radial velocity in km/s (positive = receding)

    Returns
    -------
    wavelength_rest : rest-frame wavelength array (Å)
    """
    # ── [NUMPY] Broadcasting — scalar operation on an array ────────────────
    # When you divide an array by a scalar, NumPy applies it to EVERY element
    # simultaneously. No loop needed.
    # ──────────────────────────────────────────────────────────────────────
    doppler_factor  = 1.0 + (rv_kms / C_KMS)
    wavelength_rest = wavelength / doppler_factor

    print(f"  RV correction: {rv_kms:+.2f} km/s → shifted wavelengths by "
          f"{rv_kms/C_KMS * 5000:.3f} Å at 5000 Å")

    return wavelength_rest


# ── Step 4: Signal-to-noise estimation ───────────────────────────────────────

def estimate_snr(wavelength: np.ndarray,
                 flux: np.ndarray,
                 window_center: float = 6000.0,
                 window_width: float = 20.0) -> float:
    """
    Estimate the signal-to-noise ratio in a clean, line-free spectral region.

    S/N = signal / noise = mean(flux) / std(flux)  [in a region with no lines]

    We use a region near 6000 Å which is relatively clean for K dwarfs.
    Target S/N > 200 for abundance work (PIPELINE['snr_min_science']).
    55 Cnc (V=5.95) achieves S/N > 400 in a single HARPS exposure.

    Parameters
    ----------
    wavelength    : wavelength array (Å)
    flux          : flux array
    window_center : center of the clean region (Å)
    window_width  : width of the window (Å)

    Returns
    -------
    snr : float
    """
    # ── [NUMPY] Boolean masking for window selection ───────────────────────
    # np.abs computes absolute value element-wise.
    # We select pixels within ±(window_width/2) Angstroms of window_center.
    # ──────────────────────────────────────────────────────────────────────
    mask   = np.abs(wavelength - window_center) < (window_width / 2.0)
    region = flux[mask]

    if len(region) < 10:
        print(f"  Warning: fewer than 10 pixels in SNR window. Using wider region.")
        mask   = np.abs(wavelength - window_center) < window_width
        region = flux[mask]

    # ── [NUMPY] np.mean and np.std ─────────────────────────────────────────
    # np.mean(arr) = sum / count
    # np.std(arr)  = standard deviation = sqrt( mean( (x - mean(x))² ) )
    # For photon noise in a flat region: std ≈ sqrt(mean) (Poisson stats)
    # ──────────────────────────────────────────────────────────────────────
    signal = np.mean(region)
    noise  = np.std(region)
    snr    = signal / noise if noise > 0 else 0.0

    print(f"  Estimated S/N at {window_center:.0f} Å: {snr:.0f}")
    return snr


# ── Step 5: First-look diagnostic plot ───────────────────────────────────────

def plot_spectrum_overview(wavelength: np.ndarray,
                            flux: np.ndarray,
                            header: dict,
                            save_path: Path = None):
    """
    Make a multi-panel diagnostic plot of the full spectrum.
    Three panels:
      1. Full spectrum overview (WAV_MIN – WAV_MAX)
      2. Zoom: Fe + Mg region (5150 – 5200 Å) — forest of lines
      3. Zoom: Ca I + Fe I region (6100 – 6180 Å) — cleaner region
    """
    fig = plt.figure(figsize=(14, 9))
    fig.suptitle(f"55 Cancri A — HARPS Spectrum\n"
                 f"Teff = {TEFF} K | log g = {LOGG} | [Fe/H] = {FEH:+.2f}",
                 fontsize=13, fontweight='bold')

    gs = gridspec.GridSpec(3, 1, figure=fig, hspace=0.45)

    # ── Panel 1: Full spectrum ──────────────────────────────────────────────
    ax1 = fig.add_subplot(gs[0])
    # ── [NUMPY] Downsampling with slicing ─────────────────────────────────
    # arr[::10] takes every 10th element — fast visual downsampling of
    # a 200k-pixel array without losing the line structure visually.
    # ──────────────────────────────────────────────────────────────────────
    ax1.plot(wavelength[::10], flux[::10], color='#2c5f8a', linewidth=0.4, alpha=0.9)
    ax1.set_xlabel("Wavelength (Å)", fontsize=10)
    ax1.set_ylabel("Flux (counts)", fontsize=10)
    ax1.set_title("Full spectrum overview", fontsize=10)

    key_lines = {'Mg b\n5175': 5175, 'Na D\n5893': 5893,
                 'Ca I\n6122': 6122, 'Hα\n6563': 6563}
    for label, wav in key_lines.items():
        if WAV_MIN <= wav <= WAV_MAX:
            ax1.axvline(wav, color='tomato', linewidth=0.7, alpha=0.6, linestyle='--')
            ax1.text(wav + 8, ax1.get_ylim()[1] * 0.85, label,
                     fontsize=7, color='tomato', va='top')

    # ── Panel 2: Mg I b triplet ─────────────────────────────────────────────
    ax2 = fig.add_subplot(gs[1])
    zoom1_mask = (wavelength >= 5150) & (wavelength <= 5200)
    ax2.plot(wavelength[zoom1_mask], flux[zoom1_mask], color='#2c5f8a', linewidth=0.8)
    ax2.set_xlabel("Wavelength (Å)", fontsize=10)
    ax2.set_ylabel("Flux", fontsize=10)
    ax2.set_title("Mg I b triplet (5150–5200 Å) — note high line density", fontsize=10)
    for w in [5167.32, 5172.68, 5183.60]:
        ax2.axvline(w, color='seagreen', linewidth=1.0, alpha=0.8, linestyle=':')
        ax2.text(w + 0.3, ax2.get_ylim()[1] * 0.95,
                 f'Mg I\n{w:.2f}', fontsize=7, color='seagreen', va='top')

    # ── Panel 3: Ca I + Fe I region ─────────────────────────────────────────
    ax3 = fig.add_subplot(gs[2])
    zoom2_mask = (wavelength >= 6100) & (wavelength <= 6180)
    ax3.plot(wavelength[zoom2_mask], flux[zoom2_mask], color='#2c5f8a', linewidth=0.8)
    ax3.set_xlabel("Wavelength (Å)", fontsize=10)
    ax3.set_ylabel("Flux", fontsize=10)
    ax3.set_title("Ca I + Fe I region (6100–6180 Å) — cleaner for EW measurement",
                  fontsize=10)
    for w in [6122.22, 6162.17]:
        ax3.axvline(w, color='darkorange', linewidth=1.0, alpha=0.8, linestyle=':')
        ax3.text(w + 0.3, ax3.get_ylim()[1] * 0.95,
                 f'Ca I\n{w:.2f}', fontsize=7, color='darkorange', va='top')

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"  Saved overview plot → {save_path}")
    else:
        plt.show()
    plt.close()


# ── Step 6: Summary statistics ────────────────────────────────────────────────

def spectrum_summary(wavelength: np.ndarray,
                     flux: np.ndarray,
                     header: dict) -> dict:
    """Print and return a summary of key spectrum properties."""
    # ── [NUMPY] Aggregation functions ─────────────────────────────────────
    # np.min/max/mean/median/std all operate on the entire array in one call,
    # orders of magnitude faster than a Python for loop.
    # ──────────────────────────────────────────────────────────────────────
    snr_estimate = estimate_snr(wavelength, flux)

    # ── [NUMPY] np.diff — differences between adjacent elements ───────────
    # np.diff(arr) returns arr[1:] - arr[:-1]
    # np.median of wavelength diffs gives a robust pixel scale estimate.
    # ──────────────────────────────────────────────────────────────────────
    pixel_scale = np.median(np.diff(wavelength))  # Å/pixel

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


# ── Demo mode: generate a synthetic spectrum for testing ──────────────────────

def make_synthetic_demo_spectrum() -> tuple[np.ndarray, np.ndarray]:
    """
    Generate a fake spectrum for testing the pipeline before we have
    real HARPS data. Replace with load_harps_spectrum() when data is available.

    Returns
    -------
    wavelength, flux : np.ndarrays matching what load_harps_spectrum returns
    """
    print("\n[DEMO MODE] Generating synthetic 55 Cnc-like spectrum...")
    print("  (Replace with real HARPS FITS file when available)\n")

    # ── [NUMPY] np.linspace — evenly spaced array ─────────────────────────
    # np.linspace(start, stop, num) creates 'num' values from start to stop
    # inclusive. For wavelength grids this is ideal — we specify the range
    # and pixel count directly.
    # ──────────────────────────────────────────────────────────────────────
    wavelength = np.linspace(WAV_MIN, WAV_MAX, 120_000)

    wav_norm  = (wavelength - 5500) / 1000.0
    # ── [NUMPY] np.polyval([a3,a2,a1,a0], x) = a3x³ + a2x² + a1x + a0 ───
    continuum = np.polyval([-0.02, -0.05, 0.1, 1.0], wav_norm)
    continuum = np.clip(continuum, 0.3, 1.5) * 30_000

    absorption_lines = [
        (5167.32, 0.60, 0.15),   # Mg I b (strong)
        (5172.68, 0.50, 0.15),   # Mg I b
        (5183.60, 0.65, 0.17),   # Mg I b (strongest)
        (5270.36, 0.25, 0.09),   # Fe I
        (5328.04, 0.20, 0.08),   # Fe I
        (5379.57, 0.18, 0.08),   # Fe I
        (5889.95, 0.92, 0.40),   # Na D1 (very strong — broad wings)
        (5895.92, 0.88, 0.38),   # Na D2
        (6122.22, 0.30, 0.09),   # Ca I
        (6162.17, 0.35, 0.10),   # Ca I
        (6300.30, 0.08, 0.06),   # O I forbidden (weak but critical!)
        (6302.50, 0.12, 0.07),   # Fe I (blend near O I line — tricky!)
        (6562.80, 0.45, 0.60),   # Hα (broad — chromospheric)
        (6767.77, 0.12, 0.07),   # Ni I
    ]

    flux = continuum.copy()
    for center, depth, width in absorption_lines:
        # ── [NUMPY] Gaussian applied element-wise to entire array ─────────
        # No loop — NumPy computes exp() for all 120k pixels at once.
        # ──────────────────────────────────────────────────────────────────
        flux -= depth * continuum * np.exp(-0.5 * ((wavelength - center) / width) ** 2)

    # ── [NUMPY] Reproducible random noise ─────────────────────────────────
    # Fixed seed → same "random" numbers every run. Essential for
    # reproducible science.
    # ──────────────────────────────────────────────────────────────────────
    rng          = np.random.default_rng(seed=42)
    photon_noise = rng.normal(0, np.sqrt(np.abs(flux)), size=len(flux))
    flux         = np.clip(flux + photon_noise, 0, None)

    print(f"  Synthetic spectrum: {len(wavelength):,} pixels, "
          f"{len(absorption_lines)} absorption lines, "
          f"S/N ~ {np.sqrt(np.median(flux)):.0f}")

    return wavelength, flux


# ── Main ──────────────────────────────────────────────────────────────────────

if __name__ == "__main__":

    print("=" * 60)
    print("  55 Cancri A — Spectral Acquisition & Inspection")
    print("  RYA-42 | Step 1 of the pipeline")
    print("=" * 60)

    # ── Option A: Load a real HARPS FITS file ─────────────────────────────
    # Once you've downloaded from ESO archive, point this at your file:
    #
    #   fits_file  = DATA_DIR / 'ADP.2014-09-16T11:30:15.000.fits'
    #   wavelength, flux, header = load_harps_spectrum(fits_file)
    #   wavelength = correct_radial_velocity(wavelength, STAR_55CNC['rv_kms'])
    #
    # ── Option B: Demo mode (runs right now, no data needed) ─────────────
    wavelength, flux = make_synthetic_demo_spectrum()
    header = {
        'OBJECT'  : STAR_55CNC['name'],
        'DATE-OBS': '2014-02-14',
        'EXPTIME' : 300.0,
        'INSTRUME': 'HARPS',
    }

    summary = spectrum_summary(wavelength, flux, header)
    plot_spectrum_overview(
        wavelength, flux, header,
        save_path=PATHS['plots'] / '01_spectrum_overview.png'
    )

    print("\n✓  Step 1 complete.")
    print("   Next: run 02_normalize.py to continuum-normalize the spectrum.")
    print(f"   Key output saved: {PATHS['plots'] / '01_spectrum_overview.png'}")
