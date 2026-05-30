"""
pipeline/spectra_normalize.py
==============================
Solar continuum normalization for the HARPS Ceres calibration spectra.

Loads all S1D FITS exposures, applies BERV correction, co-adds with
sigma-clipping, then fits an iterative sigma-clipping spline continuum
and normalizes to continuum = 1.0.

Linear issue : RYA-96
Data         : Solar Calibration/archive/ (HARPS S1D binary table format)
Citation     : Dumusque et al. 2021, arXiv:2009.01945, ESO 1102.D-0954(A)
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from astropy.io import fits
from scipy.interpolate import CubicSpline
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

from config.constants import PIPELINE, PATHS, PHYSICS, STAR_SOLAR


# ── FITS loading ──────────────────────────────────────────────────────────────

def _load_fits(fits_path: Path) -> tuple:
    """Load one HARPS S1D binary table FITS → (wavelength, flux, err, header)."""
    with fits.open(fits_path) as hdul:
        hdr = hdul[0].header
        # astropy drops the 'HIERARCH' prefix when converting to dict, so read
        # BERV from the Header object directly before converting.
        header = dict(hdr)
        header['HIERARCH ESO DRS BERV'] = float(hdr.get('HIERARCH ESO DRS BERV', 0.0))
        data = hdul[1].data
        wav  = data['WAVE'][0].astype(np.float64)
        flux = data['FLUX'][0].astype(np.float64)
        err  = data['ERR'][0].astype(np.float64)
    return wav, flux, err, header


def _apply_berv(wavelength: np.ndarray, berv_kms: float) -> np.ndarray:
    """Shift observed wavelengths to barycentric rest frame: λ_bary = λ_obs*(1 + BERV/c)."""
    return wavelength * (1.0 + berv_kms / PHYSICS['c_kms'])


# ── Stacking ──────────────────────────────────────────────────────────────────

def stack_spectra(fits_dir: Path) -> tuple:
    """
    Load all FITS exposures, apply per-exposure BERV correction, interpolate
    to a common wavelength grid, sigma-clip across exposures, and co-add.

    Returns (wavelength, flux_master, n_exposures).
    """
    fits_files = sorted(fits_dir.glob('*.fits'))
    if not fits_files:
        raise FileNotFoundError(f"No FITS files found in {fits_dir}")

    wav_min = PIPELINE['wav_min_A']
    wav_max = PIPELINE['wav_max_A']
    sigma   = PIPELINE['sigma_clip']

    print(f"  Loading {len(fits_files)} exposures from {fits_dir.name}/")
    spectra = []
    for f in fits_files:
        wav, flux, err, hdr = _load_fits(f)
        berv = float(hdr.get('HIERARCH ESO DRS BERV', 0.0))
        wav_bary = _apply_berv(wav, berv)
        mask = (wav_bary >= wav_min) & (wav_bary <= wav_max)
        spectra.append((wav_bary[mask], flux[mask]))
        print(f"    {f.name}: BERV={berv:+.4f} km/s  pixels={mask.sum():,}")

    # Reference grid: wavelength array of the first exposure
    wav_ref = spectra[0][0]

    # Interpolate all exposures onto the reference grid
    flux_stack = np.zeros((len(spectra), len(wav_ref)))
    for i, (wav, flux) in enumerate(spectra):
        flux_stack[i] = np.interp(wav_ref, wav, flux, left=np.nan, right=np.nan)

    # Sigma-clip across the exposure stack at each pixel, then take mean
    med = np.nanmedian(flux_stack, axis=0)
    std = np.nanstd(flux_stack, axis=0)
    clipped = np.where(np.abs(flux_stack - med) <= sigma * std, flux_stack, np.nan)
    flux_master = np.nanmean(clipped, axis=0)

    print(f"  Co-added {len(spectra)} exposures → {len(wav_ref):,} pixels")
    return wav_ref, flux_master, len(spectra)


# ── Continuum fitting ─────────────────────────────────────────────────────────

def fit_continuum(wavelength: np.ndarray, flux: np.ndarray) -> np.ndarray:
    """
    Iterative sigma-clipping spline continuum fit.

    Each iteration:
      1. Divide spectrum into windows of width continuum_knot_spacing_A
      2. Take the continuum_upper_percentile of flux in each window as anchor
      3. Fit CubicSpline through anchors
      4. Reject pixels more than sigma_clip * std BELOW the fit (absorption lines)
    Repeat for continuum_n_iter iterations. All params from PIPELINE config.
    """
    sigma    = PIPELINE['sigma_clip']
    n_iter   = PIPELINE['continuum_n_iter']
    spacing  = PIPELINE['continuum_knot_spacing_A']
    pct      = PIPELINE['continuum_upper_percentile']

    mask = np.ones(len(flux), dtype=bool)
    continuum = np.ones_like(flux)

    for iteration in range(n_iter):
        edges = np.arange(wavelength[0], wavelength[-1] + spacing, spacing)
        w_anchors, f_anchors = [], []
        for w0, w1 in zip(edges[:-1], edges[1:]):
            seg = mask & (wavelength >= w0) & (wavelength < w1)
            if seg.sum() > 10:
                w_anchors.append(0.5 * (w0 + w1))
                f_anchors.append(np.percentile(flux[seg], pct))

        if len(w_anchors) < 4:
            break

        spl = CubicSpline(w_anchors, f_anchors, extrapolate=True)
        continuum = np.clip(spl(wavelength), 1e-6, None)

        residuals = flux - continuum
        std = np.std(residuals[mask])
        new_mask = residuals > -sigma * std
        n_clipped = mask.sum() - new_mask.sum()
        mask = new_mask
        print(f"    iter {iteration + 1}/{n_iter}: "
              f"{mask.sum():,} pixels retained, {n_clipped:,} clipped")

    return continuum


# ── Output ────────────────────────────────────────────────────────────────────

def _save_normalized(wavelength: np.ndarray, flux_raw: np.ndarray,
                     continuum: np.ndarray, out_path: Path) -> np.ndarray:
    """Write (wavelength, flux_raw, continuum, flux_normalized) to CSV."""
    flux_norm = flux_raw / continuum
    pd.DataFrame({
        'wavelength_air_A': wavelength,
        'flux_raw'        : flux_raw,
        'continuum'       : continuum,
        'flux_normalized' : flux_norm,
    }).to_csv(out_path, index=False)
    print(f"  Saved → {out_path.name}  ({len(wavelength):,} rows)")
    return flux_norm


def _plot_diagnostic(wavelength: np.ndarray, flux_raw: np.ndarray,
                     continuum: np.ndarray, flux_norm: np.ndarray,
                     n_exp: int, out_path: Path) -> None:
    """3-panel diagnostic: raw+continuum, full normalized, Mg b zoom."""
    fig = plt.figure(figsize=(14, 10))
    fig.suptitle(
        f"Solar continuum normalization — HARPS Ceres calibration\n"
        f"Teff = {STAR_SOLAR['teff_K']:.0f} K   "
        f"log g = {STAR_SOLAR['logg']}   "
        f"ξ = {STAR_SOLAR['vturb_kms']} km/s   "
        f"{n_exp} exposures co-added\n"
        f"{STAR_SOLAR['citation']}  |  ESO {STAR_SOLAR['program_id']}",
        fontsize=10,
    )
    gs = gridspec.GridSpec(3, 1, figure=fig, hspace=0.52)

    # Panel 1 — raw + continuum
    ax1 = fig.add_subplot(gs[0])
    step = max(1, len(wavelength) // 8000)
    ax1.plot(wavelength[::step], flux_raw[::step],
             color='#4a90d9', lw=0.3, alpha=0.8, label='Co-added raw')
    ax1.plot(wavelength[::step], continuum[::step],
             color='tomato', lw=1.0, label='Continuum fit')
    ax1.set_ylabel('Flux (ADU)', fontsize=9)
    ax1.set_title('Co-added raw spectrum + iterative spline continuum', fontsize=9)
    ax1.legend(fontsize=8, loc='upper right')
    ax1.set_xlabel('Wavelength (Å)', fontsize=9)

    # Panel 2 — full normalized spectrum
    ax2 = fig.add_subplot(gs[1])
    ax2.plot(wavelength[::step], flux_norm[::step],
             color='#2c5f2e', lw=0.3, alpha=0.9)
    ax2.axhline(1.0, color='tomato', lw=0.7, linestyle='--', alpha=0.7)
    ax2.set_ylim(0.0, 1.15)
    ax2.set_ylabel('Normalized flux', fontsize=9)
    ax2.set_title('Normalized spectrum (continuum = 1.0)', fontsize=9)
    ax2.set_xlabel('Wavelength (Å)', fontsize=9)

    # Panel 3 — Mg b triplet zoom (normalization QA)
    ax3 = fig.add_subplot(gs[2])
    zoom = (wavelength >= 5150) & (wavelength <= 5205)
    ax3.plot(wavelength[zoom], flux_norm[zoom], color='#2c5f2e', lw=0.7)
    ax3.axhline(1.0, color='tomato', lw=0.7, linestyle='--', alpha=0.7)
    ax3.set_ylim(0.0, 1.12)
    for w, lbl in [(5167.322, 'Mg Ib₁'), (5172.684, 'Mg Ib₂'), (5183.604, 'Mg Ib₃')]:
        ax3.axvline(w, color='darkorange', lw=0.8, linestyle=':', alpha=0.85)
        ax3.text(w + 0.4, 1.06, lbl, fontsize=7, color='darkorange')
    ax3.set_xlabel('Wavelength (Å)', fontsize=9)
    ax3.set_ylabel('Normalized flux', fontsize=9)
    ax3.set_title('Mg I b triplet zoom — normalization QA (5150–5205 Å)', fontsize=9)

    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved → {out_path.name}")


# ── Pipeline entry point ──────────────────────────────────────────────────────

def run(star_id: str = 'solar') -> tuple:
    """
    Called by run_pipeline.py. Returns (wavelength, flux_normalized).

    star_id: any value containing 'solar' routes to the solar calibration path.
    """
    print(f"\n{'='*60}\n  spectra_normalize — {star_id}\n{'='*60}")

    if 'solar' not in star_id.lower():
        raise NotImplementedError(
            f"spectra_normalize not yet implemented for '{star_id}'. "
            "Only star_id containing 'solar' is supported."
        )

    fits_dir  = PATHS['solar_spectra']
    out_csv   = PATHS['solar_normalized']
    out_plot  = PATHS['solar_diagnostic']

    print(f"\n[1/4] Stacking exposures")
    wavelength, flux_raw, n_exp = stack_spectra(fits_dir)

    print(f"\n[2/4] Fitting continuum")
    continuum = fit_continuum(wavelength, flux_raw)

    print(f"\n[3/4] Saving normalized spectrum → {out_csv.relative_to(out_csv.parent.parent)}")
    flux_norm = _save_normalized(wavelength, flux_raw, continuum, out_csv)

    print(f"\n[4/4] Saving diagnostic plot → {out_plot.relative_to(out_plot.parent.parent)}")
    _plot_diagnostic(wavelength, flux_raw, continuum, flux_norm, n_exp, out_plot)

    # Quick QA report
    near_unity = np.abs(flux_norm - 1.0) < 0.05
    print(f"\n── Normalization QA ─────────────────────────────────────")
    print(f"  Pixels within 5% of continuum : {near_unity.sum():,} / {len(flux_norm):,} "
          f"({100 * near_unity.mean():.1f}%)")
    print(f"  Median normalized flux         : {np.median(flux_norm):.4f}")
    print(f"  Flux range                     : "
          f"{flux_norm.min():.3f} – {flux_norm.max():.3f}")
    print(f"────────────────────────────────────────────────────────\n")
    print(f"✓  spectra_normalize complete.")

    return wavelength, flux_norm


# ── Main ──────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    run(star_id='solar')
