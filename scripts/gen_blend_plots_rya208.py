# =============================================================================
# THE EXOPLANET CODEX
# exoplanetcodex.org  |  github.com/damienabraxas/exoplanetcodex
# =============================================================================
#
# File:         gen_blend_plots_rya208.py
# Module:       scripts (one-off diagnostic tools)
# Description:  Generate RYA-224-standard problem line diagnostic plots for
#               the two Fe I vetted exclusions confirmed in RYA-208.
#               Fe I 4918.994 Å and Fe I 4970.496 Å — both blend_flag=True.
#               Outputs attached to Linear RYA-208.
#
# Author:       Ryan Schmitt
# Contributors: Claude (Anthropic) via Claude Code
# Created:      2026-06-11
# Last modified: 2026-06-11
# Linear issue: RYA-209 (vetted exclusion documentation)
#               RYA-224 (plot standard)
#
# Usage:
#   python scripts/gen_blend_plots_rya208.py
#   Outputs: plots/problem_lines/solar/Fe_I_4919A_2026-06-11.png
#            plots/problem_lines/solar/Fe_I_4970A_2026-06-11.png
# =============================================================================

from pathlib import Path
from astropy.io import fits
import numpy as np
# Standalone-script bootstrap: these run under foreign interpreters (e.g. venv_pysme)
# and from arbitrary cwds, so put the REPO ROOT on sys.path BEFORE importing anything
# from `pipeline`. Derived from __file__, never from cwd. (RYA-313)
import os as _os_boot, sys as _sys_boot
_sys_boot.path.insert(0, _os_boot.path.dirname(_os_boot.path.dirname(
    _os_boot.path.abspath(__file__))))
from pipeline._numcompat import trapezoid as _trapezoid  # numpy>=2 removed np.trapz (RYA-313)
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.optimize import curve_fit
import warnings

SPECTRA_DIR = Path(
    '/Users/ryanschmitt/Documents/Exoplanet Codex/data/spectra/'
    'exoplanetcodex-data/Solar Calibration/archive'
)
OUT_DIR = Path('plots/problem_lines/solar')
OUT_DIR.mkdir(parents=True, exist_ok=True)

RUN_DATE   = '2026-06-11'
FONT_MONO  = 'DejaVu Sans Mono'
BG_COLOR   = '#0d1117'
TEXT_COLOR = '#ffffff'
ACCENT     = '#4af0c4'
RED        = '#eb5757'
GRID_COLOR = '#1e2a3a'


# ---------------------------------------------------------------------------
# Blend definitions — key contaminating species at each target wavelength
# ---------------------------------------------------------------------------
BLEND_SPECS = {
    4918.994: {
        'label':        'Fe I 4918.994 Å',
        'primary_wl':   4918.99349,  # stronger of the two Fe I entries
        'secondary_wl': 4918.95379,  # 2nd Fe I entry only 0.04 Å away
        'excpot':       2.865,
        'log_gf':       -0.342,
        'depth':        0.860,
        'blend_labels': [
            ('Fe I',  4918.9535, 'seagreen'),
            ('Ni I',  4918.7019, RED),
            ('V I*',  4918.9800, RED),   # centroid of V I forest
        ],
        'blend_note':   ('Fe I 4918.954 Å (Δ=0.04 Å, inseparable at HARPS R)\n'
                         'Ni I 4918.702 Å (Δ=0.29 Å, depth=0.24)\n'
                         'V I forest: 40+ lines within 0.05 Å'),
        'flag_reason':  'Two Fe I lines 0.04 Å apart; V I forest; inseparable at R~115k',
        'int_wl':       4919,
    },
    4970.496: {
        'label':        'Fe I 4970.496 Å',
        'primary_wl':   4970.49575,
        'secondary_wl': None,
        'excpot':       3.635,
        'log_gf':       -1.740,
        'depth':        0.553,
        'blend_labels': [
            ('La II',  4970.393, RED),
            ('Fe I',   4970.646, RED),
        ],
        'blend_note':   ('La II cluster at 4970.39 Å (7 lines, Δ=0.10 Å, depth~0.08)\n'
                         'Fe I 4970.646 Å (Δ=0.15 Å, depth=0.35)\n'
                         'Nd II 4970.91 Å + Co I crowd'),
        'flag_reason':  'La II cluster 0.10 Å blue; Fe I contamination 0.15 Å red',
        'int_wl':       4970,
    },
}


def load_best_solar_spectrum():
    """Load the highest-quality HARPS solar spectrum (longest exposure)."""
    fits_files = sorted(SPECTRA_DIR.glob('ADP*.fits'))
    if not fits_files:
        raise FileNotFoundError(f"No solar FITS files in {SPECTRA_DIR}")

    best_file, best_exptime = None, 0.0
    for fp in fits_files:
        with fits.open(str(fp)) as hdul:
            et = float(hdul[0].header.get('EXPTIME', 0))
            if et > best_exptime:
                best_exptime, best_file = et, fp

    with fits.open(str(best_file)) as hdul:
        spec = hdul[1].data[0]
        wl   = spec['WAVE'].astype(np.float64)
        flux = spec['FLUX'].astype(np.float64)

    print(f"Loaded: {best_file.name}  exptime={best_exptime:.0f}s  "
          f"range={wl.min():.1f}–{wl.max():.1f} Å  n={len(wl)}")
    return wl, flux


def extract_window(wl, flux, center, half_width=2.2):
    """Extract a wavelength window around the target line."""
    mask = (wl >= center - half_width) & (wl <= center + half_width)
    return wl[mask].copy(), flux[mask].copy()


def normalize_spectrum(wl_win, flux_win, center, line_half_width=0.6):
    """
    Normalize spectrum by fitting a polynomial to continuum regions
    flanking the line window. Excludes ±line_half_width Å from centre.
    """
    cont_mask = (np.abs(wl_win - center) > line_half_width)
    if cont_mask.sum() < 10:
        # Fallback: use 80th percentile
        return flux_win / np.percentile(flux_win, 80)

    from numpy.polynomial import polynomial as P
    coeffs = P.polyfit(wl_win[cont_mask], flux_win[cont_mask], deg=3)
    continuum = P.polyval(wl_win, coeffs)
    continuum = np.where(continuum > 0, continuum, np.median(flux_win))
    return flux_win / continuum


def gaussian_plus_continuum(x, amp, cen, sig, slope, offset):
    """Gaussian absorption profile on a linear continuum."""
    gauss = -amp * np.exp(-0.5 * ((x - cen) / sig) ** 2)
    return 1.0 + gauss + slope * (x - cen) + offset


def fit_gaussian(wl_norm, flux_norm, center, fit_half_width=0.5):
    """Fit a Gaussian absorption to the normalised line profile."""
    mask = np.abs(wl_norm - center) < fit_half_width
    if mask.sum() < 10:
        return None, None

    x = wl_norm[mask]
    y = flux_norm[mask]
    depth_guess = 1.0 - np.min(y)
    p0 = [max(depth_guess, 0.05), center, 0.08, 0.0, 0.0]
    bounds = ([0, center - 0.3, 0.01, -5, -0.5], [1.5, center + 0.3, 0.5, 5, 0.5])
    try:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            popt, _ = curve_fit(gaussian_plus_continuum, x, y, p0=p0, bounds=bounds,
                                maxfev=10000)
        return popt, x
    except Exception:
        return None, None


def estimate_ew(wl_norm, flux_norm, center, ew_half_width=0.4):
    """Estimate EW by direct integration of (1 - flux) over the line window."""
    mask = np.abs(wl_norm - center) < ew_half_width
    if mask.sum() < 5:
        return float('nan')
    x = wl_norm[mask]
    y = np.clip(1.0 - flux_norm[mask], 0, None)
    return float(_trapezoid(y, x) * 1000)  # Å → mÅ


def make_problem_line_plot(wl_all, flux_all, target_wl, spec_info):
    """Generate one RYA-224-standard problem line diagnostic plot."""
    center    = spec_info['primary_wl']
    int_wl    = spec_info['int_wl']

    # Extract and normalise window
    wl_win, flux_win = extract_window(wl_all, flux_all, center, half_width=2.5)
    flux_norm = normalize_spectrum(wl_win, flux_win, center)

    # Gaussian fit
    popt, fit_x = fit_gaussian(wl_win, flux_norm, center)
    ew_mA = estimate_ew(wl_win, flux_norm, center)

    # Residuals
    if popt is not None:
        fit_mask = np.abs(wl_win - center) < 0.5
        flux_model = gaussian_plus_continuum(wl_win[fit_mask], *popt)
        residuals  = flux_norm[fit_mask] - flux_model
        noise_est  = np.std(flux_norm[(np.abs(wl_win - center) > 1.5)])
    else:
        fit_mask  = np.zeros(len(wl_win), dtype=bool)
        residuals = np.array([])
        noise_est = 0.02

    # -----------------------------------------------------------------------
    # Figure setup
    # -----------------------------------------------------------------------
    fig = plt.figure(figsize=(10, 8), facecolor=BG_COLOR)
    gs  = gridspec.GridSpec(2, 1, height_ratios=[0.70, 0.30],
                            hspace=0.05, top=0.88, bottom=0.12,
                            left=0.10, right=0.97)
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1], sharex=ax1)

    for ax in (ax1, ax2):
        ax.set_facecolor(BG_COLOR)
        ax.tick_params(colors=TEXT_COLOR, labelsize=9)
        ax.spines[:].set_color('#3d4a5c')
        ax.grid(True, color=GRID_COLOR, alpha=0.4, linewidth=0.6)

    plt.setp(ax1.get_xticklabels(), visible=False)

    # -----------------------------------------------------------------------
    # Panel 1 — spectrum
    # -----------------------------------------------------------------------
    plot_mask = np.abs(wl_win - center) <= 2.0
    ax1.plot(wl_win[plot_mask], flux_norm[plot_mask],
             color='#cccccc', linewidth=0.9, alpha=0.9, label='Observed')

    # Gaussian fit overlay
    if popt is not None:
        x_dense = np.linspace(center - 0.5, center + 0.5, 500)
        y_fit   = gaussian_plus_continuum(x_dense, *popt)
        ax1.plot(x_dense, y_fit, color=ACCENT, linewidth=1.5,
                 alpha=0.9, label='Gaussian fit')

    # Continuum reference
    ax1.axhline(1.0, color='#6d7a8a', linestyle='--', linewidth=0.8, alpha=0.7)

    # Target line
    ax1.axvline(center, color=ACCENT, linestyle='--', linewidth=1.0, alpha=0.85)

    # Secondary Fe I line (if present)
    if spec_info['secondary_wl'] is not None:
        ax1.axvline(spec_info['secondary_wl'], color='#7bd4f0',
                    linestyle=':', linewidth=1.0, alpha=0.7,
                    label=f"Fe I {spec_info['secondary_wl']:.3f}")

    # Blend species markers
    blend_x_used = set()
    for (species, blend_wl, color) in spec_info['blend_labels']:
        if abs(blend_wl - center) <= 2.0:
            if blend_wl not in blend_x_used:
                ax1.axvline(blend_wl, color=color, linestyle='--',
                            linewidth=1.0, alpha=0.7)
                ax1.text(blend_wl, 0.97, f"{species}\n{blend_wl:.3f}",
                         color=color, fontsize=7, fontfamily=FONT_MONO,
                         ha='center', va='top', rotation=90,
                         transform=ax1.get_xaxis_transform())
                blend_x_used.add(blend_wl)

    ax1.set_ylabel('Normalised Flux', color=TEXT_COLOR, fontsize=10,
                   fontfamily=FONT_MONO)
    ax1.set_ylim(0.0, 1.15)
    ax1.yaxis.label.set_color(TEXT_COLOR)
    ax1.legend(loc='lower left', fontsize=8, framealpha=0.3,
               facecolor='#1a2233', edgecolor='#3d4a5c', labelcolor=TEXT_COLOR)

    # -----------------------------------------------------------------------
    # Annotation box
    # -----------------------------------------------------------------------
    ew_str   = f"{ew_mA:.1f} mÅ (estimated)" if np.isfinite(ew_mA) else "N/A"
    snr_est  = (1.0 / noise_est) if noise_est > 0 else 0.0
    chi2_str = "N/A"
    if popt is not None and len(residuals) > 1:
        chi2_str = f"{np.sum(residuals**2 / noise_est**2):.1f}"

    ann_lines = [
        f"Star:        solar (Sun)",
        f"Line:        {spec_info['label']}",
        f"EP / log gf: {spec_info['excpot']:.3f} eV / {spec_info['log_gf']:.3f}",
        f"Central dep: {spec_info['depth']:.3f}",
        f"EW (est.):   {ew_str}",
        f"A(X) 1DLTE:  excluded — blend",
        f"A(X) NLTE:   excluded — blend",
        f"Grade:       excluded (blend_flag=True)",
        f"SNR:         ~{snr_est:.0f}",
        f"χ²:          {chi2_str}",
        f"Flag reason: confirmed blend (RYA-208)",
    ]
    ann_text = '\n'.join(ann_lines)
    ax1.text(0.98, 0.97, ann_text,
             transform=ax1.transAxes,
             fontsize=7.5, fontfamily=FONT_MONO,
             color=TEXT_COLOR, va='top', ha='right',
             bbox=dict(boxstyle='round,pad=0.5', facecolor='#0d1a2d',
                       edgecolor=ACCENT, alpha=0.85))

    # -----------------------------------------------------------------------
    # Panel 2 — residuals
    # -----------------------------------------------------------------------
    if len(residuals) > 0:
        res_wl = wl_win[fit_mask]
        ax2.plot(res_wl, residuals, color='#cccccc', linewidth=0.8, alpha=0.9)
        ax2.axhline(0, color='#6d7a8a', linewidth=0.8)
        ax2.fill_between(res_wl, -noise_est, noise_est,
                         color='#4d5a6a', alpha=0.3)
        # Flag large residuals
        large = np.abs(residuals) > 3 * noise_est
        if large.sum() > 0:
            ax2.scatter(res_wl[large], residuals[large], color=RED, s=20, zorder=5)
    else:
        ax2.text(0.5, 0.5, 'No fit — residuals unavailable',
                 transform=ax2.transAxes, ha='center', va='center',
                 color='#6d7a8a', fontsize=9, fontfamily=FONT_MONO)

    ax2.set_xlabel('Wavelength (Å, air)', color=TEXT_COLOR, fontsize=10,
                   fontfamily=FONT_MONO)
    ax2.set_ylabel('Residual', color=TEXT_COLOR, fontsize=9, fontfamily=FONT_MONO)
    ax2.xaxis.label.set_color(TEXT_COLOR)

    # -----------------------------------------------------------------------
    # Title
    # -----------------------------------------------------------------------
    fig.suptitle(
        f"PROBLEM LINE DIAGNOSTIC — {spec_info['label']}\n"
        f"blend_flag=True — confirmed non-separable blend at HARPS R~115,000",
        color=TEXT_COLOR, fontsize=11, fontfamily=FONT_MONO, y=0.98,
        fontweight='bold',
    )

    # -----------------------------------------------------------------------
    # Blend note panel below title
    # -----------------------------------------------------------------------
    fig.text(0.10, 0.905,
             f"Blend context:  {spec_info['blend_note'].splitlines()[0]}",
             color='#f0c04a', fontsize=8, fontfamily=FONT_MONO)

    # -----------------------------------------------------------------------
    # Footer
    # -----------------------------------------------------------------------
    fig.text(0.5, 0.02,
             f"The Exoplanet Codex | exoplanetcodex.org  |  "
             f"Run date: {RUN_DATE}  |  Pipeline: iSpec + SPECTRUM  |  "
             f"RYA-208 / RYA-209",
             ha='center', color='#6d7a8a', fontsize=7.5, fontfamily=FONT_MONO)

    # -----------------------------------------------------------------------
    # Save
    # -----------------------------------------------------------------------
    out_path = OUT_DIR / f"Fe_I_{int_wl}A_{RUN_DATE}.png"
    fig.savefig(str(out_path), dpi=150, facecolor=BG_COLOR,
                bbox_inches='tight')
    plt.close(fig)
    print(f"Saved: {out_path}")
    return out_path


def main():
    print("Loading solar spectrum...")
    wl_all, flux_all = load_best_solar_spectrum()

    output_paths = []
    for target_wl, spec_info in BLEND_SPECS.items():
        print(f"\nGenerating plot for Fe I {target_wl} Å...")
        out_path = make_problem_line_plot(wl_all, flux_all, target_wl, spec_info)
        output_paths.append(out_path)

    print(f"\nGenerated {len(output_paths)} plots:")
    for p in output_paths:
        print(f"  {p}")
    return output_paths


if __name__ == '__main__':
    main()
