"""
pipeline/spectra_normalize.py
==============================
Continuum normalization for HARPS S1D/ADP spectra.

Supports solar (Dumusque feed) and Procyon (20 HARPS Phase 3 ADPs).
Generalised via star routing — solar behaviour is bit-for-bit unchanged.

Solar  : RYA-96   — Dumusque HARPS direct solar feed (BERV correction applied)
Procyon: RYA-274  — 20 HARPS Phase 3 ADPs, SPECSYS=BARYCENT (no BERV applied)
"""

import argparse
import datetime
import subprocess
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

from config.constants import (
    PIPELINE, PATHS, PHYSICS, CONTINUUM_PARAMS,
    STAR_SOLAR, STAR_PROCYON,
)

# ── FITS loading ──────────────────────────────────────────────────────────────

def _load_fits(fits_path: Path, require_barycentric: bool = False) -> tuple:
    """Load one HARPS S1D/ADP binary table FITS → (wavelength, flux, err, header_dict).

    require_barycentric=True: asserts SPECSYS=BARYCENT and skips BERV correction.
    require_barycentric=False (solar default): applies BERV correction.
    """
    with fits.open(fits_path) as hdul:
        hdr = hdul[0].header
        header = dict(hdr)
        # Read BERV before dict conversion (HIERARCH prefix may be dropped)
        header['HIERARCH ESO DRS BERV'] = float(hdr.get('HIERARCH ESO DRS BERV', 0.0))
        # Record SPECSYS for provenance
        header['SPECSYS'] = str(hdr.get('SPECSYS', 'UNKNOWN'))
        data = hdul[1].data
        wav  = data['WAVE'][0].astype(np.float64)
        flux = data['FLUX'][0].astype(np.float64)
        err  = data['ERR'][0].astype(np.float64)

    if require_barycentric:
        if header['SPECSYS'] != 'BARYCENT':
            raise ValueError(
                f"{fits_path.name}: expected SPECSYS=BARYCENT, got {header['SPECSYS']}. "
                "Applying BERV to an already-corrected spectrum would double-correct."
            )

    return wav, flux, err, header


def _apply_berv(wavelength: np.ndarray, berv_kms: float) -> np.ndarray:
    """Shift observed wavelengths to barycentric rest frame: λ_bary = λ_obs*(1 + BERV/c)."""
    return wavelength * (1.0 + berv_kms / PHYSICS['c_kms'])


# ── Stacking ──────────────────────────────────────────────────────────────────

def stack_spectra(fits_dir: Path, star_id: str = 'solar') -> tuple:
    """
    Load all FITS exposures, co-add with sigma-clipping.

    Solar path  : applies per-exposure BERV correction.
    Procyon path: confirms SPECSYS=BARYCENT per file; no BERV applied.

    Returns (wavelength, flux_master, n_exposures, per_file_info).
    per_file_info is a list of dicts with filename/snr/mjd/date_obs/specsys/berv_kms.
    """
    fits_files = sorted(fits_dir.glob('*.fits'))
    if not fits_files:
        raise FileNotFoundError(f"No FITS files found in {fits_dir}")

    is_procyon = 'procyon' in star_id.lower()
    wav_min    = PIPELINE['wav_min_A']
    wav_max    = PIPELINE['wav_max_A']
    cparams    = CONTINUUM_PARAMS.get(star_id.lower(), CONTINUUM_PARAMS['solar'])
    sigma      = cparams['sigma_clip']

    print(f"  Loading {len(fits_files)} exposures from {fits_dir.name}/")
    spectra    = []
    file_info  = []

    for f in fits_files:
        wav, flux, err, hdr = _load_fits(f, require_barycentric=is_procyon)
        berv_kms = float(hdr.get('HIERARCH ESO DRS BERV', 0.0))
        snr      = float(hdr.get('SNR', 0.0))
        mjd      = float(hdr.get('MJD-OBS', 0.0))
        date_obs = str(hdr.get('DATE-OBS', '?'))
        prog_id  = str(hdr.get('HIERARCH ESO OBS PROG ID', hdr.get('PROG_ID', '?')))
        specsys  = hdr.get('SPECSYS', '?')

        if is_procyon:
            # Already barycentric — record BERV for provenance but do NOT apply
            wav_final = wav
            print(f"    {f.name}: SPECSYS={specsys}  BERV={berv_kms:+.4f} km/s (not applied)  SNR={snr:.1f}")
        else:
            wav_final = _apply_berv(wav, berv_kms)
            print(f"    {f.name}: BERV={berv_kms:+.4f} km/s  SNR={snr:.1f}")

        mask = (wav_final >= wav_min) & (wav_final <= wav_max)
        spectra.append((wav_final[mask], flux[mask]))
        file_info.append({
            'filename': f.name, 'snr': snr, 'mjd_obs': mjd,
            'date_obs': date_obs, 'specsys': specsys,
            'berv_kms': berv_kms, 'prog_id': prog_id,
        })
        print(f"      pixels retained: {mask.sum():,}")

    # Reference grid: wavelength array of the first exposure (native HARPS sampling)
    wav_ref = spectra[0][0]

    # Interpolate all exposures onto the reference grid
    flux_stack = np.zeros((len(spectra), len(wav_ref)))
    for i, (wav, flux) in enumerate(spectra):
        flux_stack[i] = np.interp(wav_ref, wav, flux, left=np.nan, right=np.nan)

    # Sigma-clip across the exposure stack at each pixel, then take mean
    med     = np.nanmedian(flux_stack, axis=0)
    std     = np.nanstd(flux_stack, axis=0)
    clipped = np.where(np.abs(flux_stack - med) <= sigma * std, flux_stack, np.nan)
    n_rejected = int(np.sum(np.isnan(clipped)) - np.sum(np.isnan(flux_stack)))
    flux_master = np.nanmean(clipped, axis=0)

    # Combined SNR ~ quadrature sum
    combined_snr = float(np.sqrt(np.sum([fi['snr']**2 for fi in file_info])))
    mjd_min = min(fi['mjd_obs'] for fi in file_info)
    mjd_max = max(fi['mjd_obs'] for fi in file_info)

    print(f"  Co-added {len(spectra)} exposures → {len(wav_ref):,} pixels")
    print(f"  Sigma-clip ({sigma}σ) rejected {n_rejected:,} pixel-epoch samples")
    print(f"  Grid: native HARPS S1D sampling  "
          f"{wav_ref[0]:.2f}–{wav_ref[-1]:.2f} Å  step≈{np.median(np.diff(wav_ref)):.4f} Å")
    print(f"  MJD range: {mjd_min:.5f} – {mjd_max:.5f}")
    print(f"  Combined SNR (quadrature): {combined_snr:.1f}")

    for fi in file_info:
        fi['combined_snr'] = combined_snr
        fi['n_spectra']    = len(spectra)
        fi['n_rejected_pixel_samples'] = n_rejected

    return wav_ref, flux_master, len(spectra), file_info


# ── Continuum fitting ─────────────────────────────────────────────────────────

def fit_continuum(wavelength: np.ndarray, flux: np.ndarray,
                  star_id: str = 'solar') -> np.ndarray:
    """
    Iterative sigma-clipping spline continuum fit.

    Parameters come from CONTINUUM_PARAMS[star_id]; solar values are identical
    to the PIPELINE dict to maintain bit-for-bit backward compatibility.
    """
    cparams = CONTINUUM_PARAMS.get(star_id.lower(), CONTINUUM_PARAMS['solar'])
    sigma   = cparams['sigma_clip']
    n_iter  = cparams['n_iter']
    spacing = cparams['knot_spacing_A']
    pct     = cparams['upper_percentile']

    print(f"  Continuum params [{star_id}]: "
          f"spacing={spacing} Å  pct={pct}  n_iter={n_iter}  sigma={sigma}")

    mask     = np.ones(len(flux), dtype=bool)
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
            print(f"    iter {iteration + 1}: only {len(w_anchors)} anchor windows — stopping early")
            break

        spl      = CubicSpline(w_anchors, f_anchors, extrapolate=True)
        continuum = np.clip(spl(wavelength), 1e-6, None)

        residuals = flux - continuum
        std       = np.std(residuals[mask])
        new_mask  = residuals > -sigma * std
        n_clipped = mask.sum() - new_mask.sum()
        mask      = new_mask
        print(f"    iter {iteration + 1}/{n_iter}: "
              f"{mask.sum():,} pixels retained, {n_clipped:,} clipped")

    return continuum


# ── Output ────────────────────────────────────────────────────────────────────

def _save_normalized(wavelength: np.ndarray, flux_raw: np.ndarray,
                     continuum: np.ndarray, out_path: Path,
                     provenance: str = '') -> np.ndarray:
    """Write (wavelength, flux_raw, continuum, flux_normalized) to CSV.

    If provenance is provided, writes it as '#'-prefixed comment lines before
    the column header so the file is self-describing. pandas reads it back with
    read_csv(..., comment='#').  Solar path passes provenance='' to keep the
    existing file format unchanged.
    """
    flux_norm = flux_raw / continuum
    df = pd.DataFrame({
        'wavelength_air_A': wavelength,
        'flux_raw'        : flux_raw,
        'continuum'       : continuum,
        'flux_normalized' : flux_norm,
    })
    if provenance:
        with open(out_path, 'w') as fh:
            fh.write(provenance + '\n')
            fh.write(df.to_csv(index=False))
    else:
        df.to_csv(out_path, index=False)
    print(f"  Saved → {out_path.name}  ({len(wavelength):,} rows)")
    return flux_norm


def _build_provenance(star_id: str, file_info: list, continuum_params: dict,
                      wavelength: np.ndarray, out_csv: Path) -> str:
    """Build the provenance header block for the output file."""
    try:
        commit = subprocess.check_output(
            ['git', 'rev-parse', '--short', 'HEAD'],
            cwd=out_csv.parent, stderr=subprocess.DEVNULL
        ).decode().strip()
    except Exception:
        commit = 'unknown'

    filenames  = [fi['filename'] for fi in file_info]
    programs   = sorted(set(fi['prog_id'] for fi in file_info))
    mjds       = [fi['mjd_obs'] for fi in file_info]
    snrs       = [fi['snr'] for fi in file_info]
    combined_snr = file_info[0]['combined_snr'] if file_info else 0.0
    n_rejected   = file_info[0].get('n_rejected_pixel_samples', 0) if file_info else 0

    lines = [
        "# ─────────────────────────────────────────────────────────────────────",
        f"# PROVENANCE: {out_csv.name}",
        f"# Generated : {datetime.datetime.utcnow().strftime('%Y-%m-%dT%H:%M:%S')} UTC",
        f"# Script    : pipeline/spectra_normalize.py  commit {commit}",
        f"# Star      : {star_id}",
        "# ─────────────────────────────────────────────────────────────────────",
        f"# Source files ({len(filenames)}):",
    ]
    for fn in filenames:
        lines.append(f"#   {fn}")
    lines += [
        f"# Programs   : {', '.join(programs)}",
        f"# MJD range  : {min(mjds):.5f} – {max(mjds):.5f}",
        f"# Dates      : {sorted(set(fi['date_obs'][:10] for fi in file_info))}",
        f"# Per-file SNR: {[round(s,1) for s in snrs]}",
        f"# Combined SNR (quadrature): {combined_snr:.1f}",
        f"# SPECSYS    : {file_info[0]['specsys'] if file_info else '?'}  "
         "(BERV applied: " + ("No — already BARYCENT" if 'procyon' in star_id.lower() else "Yes — per-exposure") + ")",
        "# ─────────────────────────────────────────────────────────────────────",
        "# Co-add method  : sigma-clipping stack → nanmean",
        f"# Stack sigma    : {continuum_params['sigma_clip']}",
        f"# Rejected pixel-epoch samples: {n_rejected:,}",
        f"# Grid           : native HARPS S1D sampling (first-exposure reference)",
        f"# Grid range     : {wavelength[0]:.2f}–{wavelength[-1]:.2f} Å",
        f"# Grid step (med): {np.median(np.diff(wavelength)):.4f} Å",
        "# ─────────────────────────────────────────────────────────────────────",
        "# Continuum method: iterative sigma-clipping CubicSpline",
        f"# knot_spacing_A : {continuum_params['knot_spacing_A']}",
        f"# upper_pct      : {continuum_params['upper_percentile']}",
        f"# n_iter         : {continuum_params['n_iter']}",
        f"# sigma_clip     : {continuum_params['sigma_clip']}",
        "# ─────────────────────────────────────────────────────────────────────",
    ]
    return '\n'.join(lines)


def _plot_diagnostic(wavelength: np.ndarray, flux_raw: np.ndarray,
                     continuum: np.ndarray, flux_norm: np.ndarray,
                     n_exp: int, star_id: str, out_path: Path) -> None:
    """4-panel diagnostic: raw+continuum, full normalized, Hβ Balmer wing, Fe-rich zoom."""
    is_procyon = 'procyon' in star_id.lower()
    star_params = STAR_PROCYON if is_procyon else STAR_SOLAR
    star_name   = star_params.get('name', star_id)

    fig = plt.figure(figsize=(14, 12))
    fig.suptitle(
        f"{star_name} continuum normalization — HARPS\n"
        f"Teff = {star_params['teff_K']:.0f} K   "
        f"log g = {star_params['logg']}   "
        f"ξ = {star_params['vturb_kms']} km/s   "
        f"{n_exp} exposures co-added",
        fontsize=10,
    )
    gs = gridspec.GridSpec(4, 1, figure=fig, hspace=0.58)

    step = max(1, len(wavelength) // 8000)

    # Panel 1 — raw + continuum
    ax1 = fig.add_subplot(gs[0])
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

    if is_procyon:
        # Panel 3 — Hβ Balmer wing (key normalization diagnostic for F stars)
        ax3 = fig.add_subplot(gs[2])
        zoom3 = (wavelength >= 4830) & (wavelength <= 4910)
        ax3.plot(wavelength[zoom3], flux_norm[zoom3], color='#2c5f2e', lw=0.7)
        ax3.axhline(1.0, color='tomato', lw=0.7, linestyle='--', alpha=0.7)
        ax3.axvline(4861.33, color='royalblue', lw=0.8, linestyle=':', alpha=0.9)
        ax3.text(4861.8, ax3.get_ylim()[1] * 0.85 if ax3.get_ylim()[1] > 0 else 0.9,
                 'Hβ 4861', fontsize=7, color='royalblue')
        ax3.set_xlabel('Wavelength (Å)', fontsize=9)
        ax3.set_ylabel('Normalized flux', fontsize=9)
        ax3.set_title('Hβ Balmer wing — continuum placement QA (4830–4910 Å)', fontsize=9)

        # Panel 4 — Fe-rich region (Procyon optical Fe pool core)
        ax4 = fig.add_subplot(gs[3])
        zoom4 = (wavelength >= 5000) & (wavelength <= 5100)
        ax4.plot(wavelength[zoom4], flux_norm[zoom4], color='#2c5f2e', lw=0.7)
        ax4.axhline(1.0, color='tomato', lw=0.7, linestyle='--', alpha=0.7)
        for w, lbl in [(5041.072, 'Si I'), (5055.984, 'Si I'), (5079.740, 'Fe I')]:
            ax4.axvline(w, color='darkorange', lw=0.8, linestyle=':', alpha=0.85)
            ax4.text(w + 0.3, 1.04, lbl, fontsize=6, color='darkorange', rotation=90)
        ax4.set_xlabel('Wavelength (Å)', fontsize=9)
        ax4.set_ylabel('Normalized flux', fontsize=9)
        ax4.set_title('Fe-rich region zoom — normalization QA (5000–5100 Å)', fontsize=9)
    else:
        # Solar: Mg b zoom (original Panel 3)
        ax3 = fig.add_subplot(gs[2])
        zoom3 = (wavelength >= 5150) & (wavelength <= 5205)
        ax3.plot(wavelength[zoom3], flux_norm[zoom3], color='#2c5f2e', lw=0.7)
        ax3.axhline(1.0, color='tomato', lw=0.7, linestyle='--', alpha=0.7)
        for w, lbl in [(5167.322, 'Mg Ib₁'), (5172.684, 'Mg Ib₂'), (5183.604, 'Mg Ib₃')]:
            ax3.axvline(w, color='darkorange', lw=0.8, linestyle=':', alpha=0.85)
            ax3.text(w + 0.4, 1.06, lbl, fontsize=7, color='darkorange')
        ax3.set_xlabel('Wavelength (Å)', fontsize=9)
        ax3.set_ylabel('Normalized flux', fontsize=9)
        ax3.set_title('Mg I b triplet zoom — normalization QA (5150–5205 Å)', fontsize=9)

        # Panel 4 — Fe I dense region
        ax4 = fig.add_subplot(gs[3])
        zoom4 = (wavelength >= 5200) & (wavelength <= 5250)
        ax4.plot(wavelength[zoom4], flux_norm[zoom4], color='#2c5f2e', lw=0.7)
        ax4.axhline(1.0, color='tomato', lw=0.7, linestyle='--', alpha=0.7)
        ax4.set_xlabel('Wavelength (Å)', fontsize=9)
        ax4.set_ylabel('Normalized flux', fontsize=9)
        ax4.set_title('5200–5250 Å zoom', fontsize=9)

    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved → {out_path.name}")


# ── Pipeline entry point ──────────────────────────────────────────────────────

def run(star_id: str = 'solar') -> tuple:
    """
    Called by run_pipeline.py or directly via --star flag.
    Returns (wavelength, flux_normalized).
    """
    print(f"\n{'='*60}\n  spectra_normalize — {star_id}\n{'='*60}")

    is_procyon = 'procyon' in star_id.lower()
    is_solar   = 'solar'   in star_id.lower()

    if not is_solar and not is_procyon:
        raise NotImplementedError(
            f"spectra_normalize not yet implemented for '{star_id}'. "
            "Supported: 'solar', 'procyon'."
        )

    if is_solar:
        fits_dir  = PATHS['solar_spectra']
        out_csv   = PATHS['solar_normalized']
        out_plot  = PATHS['solar_diagnostic']
    else:
        fits_dir  = PATHS['procyon_spectra']
        out_csv   = PATHS['procyon_normalized']
        out_plot  = PATHS['procyon_diagnostic']

    print(f"\n[1/4] Stacking exposures")
    wavelength, flux_raw, n_exp, file_info = stack_spectra(fits_dir, star_id=star_id)

    print(f"\n[2/4] Fitting continuum")
    continuum = fit_continuum(wavelength, flux_raw, star_id=star_id)

    # Build provenance before writing — embedded in CSV for Procyon, printed for both
    cparams    = CONTINUUM_PARAMS.get(star_id.lower(), CONTINUUM_PARAMS['solar'])
    provenance = _build_provenance(star_id, file_info, cparams, wavelength, out_csv)

    print(f"\n[3/4] Saving normalized spectrum → {out_csv.relative_to(out_csv.parent.parent)}")
    # Solar: no provenance header in CSV (preserves existing file format for regression)
    # Procyon: embed provenance as '#' comment lines at top of CSV
    csv_provenance = '' if is_solar else provenance
    flux_norm = _save_normalized(wavelength, flux_raw, continuum, out_csv, csv_provenance)

    print(f"\n[4/4] Saving diagnostic plot → {out_plot.relative_to(out_plot.parent.parent)}")
    _plot_diagnostic(wavelength, flux_raw, continuum, flux_norm, n_exp, star_id, out_plot)

    print(f"\n── Provenance ───────────────────────────────────────────")
    print(provenance)
    print(f"────────────────────────────────────────────────────────")

    # Quick QA report
    near_unity = np.abs(flux_norm - 1.0) < 0.05

    # Continuum sanity: clean windows away from strong lines
    if is_procyon:
        # Fe-continuum windows in Procyon optical range (avoid Balmer + known blends)
        win1 = (wavelength >= 5150) & (wavelength <= 5155)
        win2 = (wavelength >= 5380) & (wavelength <= 5385)
        win3 = (wavelength >= 5570) & (wavelength <= 5575)
    else:
        win1 = (wavelength >= 5070) & (wavelength <= 5075)
        win2 = (wavelength >= 5380) & (wavelength <= 5385)
        win3 = (wavelength >= 5570) & (wavelength <= 5575)
    clean_windows = win1 | win2 | win3
    cont_median = float(np.median(flux_norm[clean_windows])) if clean_windows.sum() > 5 else float('nan')
    # Continuum-PLACEMENT metric: the upper envelope of the normalized flux must
    # sit at ~1.0 (the true continuum is where the line-free peaks reach). This
    # is the correct continuum-quality gate. The clean-window MEDIAN, by
    # contrast, is biased low for line-rich stars (F dwarfs have no truly
    # line-free optical windows — even "clean" 5 Å windows carry weak-line haze),
    # so it measures line blanketing, not continuum error. Report both.
    cont_placement = float(np.percentile(flux_norm, 99.0))
    cont_win_upper = (float(np.percentile(flux_norm[clean_windows], 95))
                      if clean_windows.sum() > 5 else float('nan'))

    print(f"\n── Normalization QA ─────────────────────────────────────")
    print(f"  Pixels within 5% of continuum : {near_unity.sum():,} / {len(flux_norm):,} "
          f"({100 * near_unity.mean():.1f}%)")
    print(f"  Median normalized flux         : {np.median(flux_norm):.4f}")
    print(f"  Flux range                     : "
          f"{flux_norm.min():.3f} – {flux_norm.max():.3f}")
    print(f"  Continuum placement (p99)      : {cont_placement:.4f}  (target 1.0 ± 0.01) ← gate")
    print(f"  Clean-window upper env (p95)   : {cont_win_upper:.4f}")
    print(f"  Clean-window median            : {cont_median:.4f}  "
          f"(line-haze-dominated for F stars — not a continuum-error gate)")
    print(f"────────────────────────────────────────────────────────\n")
    print(f"✓  spectra_normalize [{star_id}] complete.")

    return wavelength, flux_norm


# ── Main ──────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Continuum normalization')
    parser.add_argument('--star', default='solar',
                        help="Star ID: 'solar' (default) or 'procyon'")
    args = parser.parse_args()
    run(star_id=args.star)
