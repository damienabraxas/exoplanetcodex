"""
pipeline/diagnostics_abundance.py
==================================
Solar Fe I abundance diagnostics from measured EWs.

Loads solar_ew.csv + linelist_solar.csv, merges on wavelength, derives a
rough COG abundance A(Fe I) per line, and produces four diagnostic plots
used to judge the quality of the stellar parameter solution before
proceeding to 55 Cnc A.

Plots
-----
1  A(Fe I) vs Excitation Potential  — slope → Teff sensitivity
2  A(Fe I) vs Reduced EW            — slope → vturb sensitivity
3  A(Fe I) vs Wavelength            — trends → blends / continuum errors
4  Per-line A(Fe I) histogram       — mean, std, N

Usage
-----
    python -m pipeline.diagnostics_abundance
    # or
    from pipeline.diagnostics_abundance import run
    run()

Linear issue: RYA-160
"""

import warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

from config.constants import PATHS, STAR_SOLAR

TITLE_SUFFIX = "Solar HARPS — Dumusque ESO 1102.D-0954"
ASPLUND_FE   = 7.46
THETA_SOLAR  = 5040.0 / STAR_SOLAR['teff_K']   # ≈ 0.8721

DIAG_DIR = PATHS['plots'] / 'diagnostics'


# ── Data loading ──────────────────────────────────────────────────────────────

def _load_and_merge() -> pd.DataFrame:
    """
    Load solar_ew.csv and linelist_solar.csv; merge on wavelength (±0.05 Å).
    Returns a DataFrame with EW columns plus excitation_potential_eV and log_gf.
    """
    ew_path = PATHS['solar_ew']
    ll_path = PATHS['linelist_solar']

    if not ew_path.exists():
        raise FileNotFoundError(
            f"solar_ew.csv not found at {ew_path}\nRun lines_fit.py first."
        )
    if not ll_path.exists():
        raise FileNotFoundError(f"linelist_solar.csv not found at {ll_path}")

    ew = pd.read_csv(ew_path)
    ll = pd.read_csv(ll_path, comment='#', low_memory=False)

    ll_fe = ll[ll['element'] == 'Fe'][
        ['wavelength_air_A', 'excitation_potential_eV', 'log_gf']
    ].copy()

    # Nearest-neighbour wavelength join within ±0.05 Å
    merged = pd.merge_asof(
        ew.sort_values('wavelength_air_A'),
        ll_fe.sort_values('wavelength_air_A'),
        on='wavelength_air_A',
        direction='nearest',
        tolerance=0.05,
    )
    return merged


def _filter_fe(df: pd.DataFrame) -> pd.DataFrame:
    """Apply science-quality cuts; derive A(Fe I) via rough COG inversion."""
    fe = df[
        (df['element'] == 'Fe') &
        (df['ion'] == 'I') &
        (df['blend_flag'] == False) &
        (df['ew_mA'].between(10.0, 150.0)) &
        (df['chi2'] < 0.05) &
        df['excitation_potential_eV'].notna() &
        df['log_gf'].notna()
    ].copy()

    fe['A_Fe_COG'] = (
        np.log10(fe['ew_mA'] / fe['wavelength_air_A'])
        - fe['log_gf']
        + fe['excitation_potential_eV'] * THETA_SOLAR
        + ASPLUND_FE
    )
    fe = fe[np.isfinite(fe['A_Fe_COG'])].reset_index(drop=True)

    fe['ew_regime'] = pd.cut(
        fe['ew_mA'],
        bins=[0, 30, 90, np.inf],
        labels=['weak (<30 mÅ)', 'medium (30–90 mÅ)', 'strong (>90 mÅ)'],
    )
    return fe


# ── Plot helpers ──────────────────────────────────────────────────────────────

_REGIME_COLOURS = {
    'weak (<30 mÅ)'     : '#5b9bd5',
    'medium (30–90 mÅ)' : '#70ad47',
    'strong (>90 mÅ)'   : '#e05c5c',
}


def _scatter_by_regime(ax, x, y, regime):
    for label, colour in _REGIME_COLOURS.items():
        mask = regime == label
        if mask.any():
            ax.scatter(x[mask], y[mask], c=colour, s=18, alpha=0.75,
                       label=label, zorder=3)


def _fit_line(x, y):
    coeffs = np.polyfit(x, y, 1)
    x_line = np.linspace(x.min(), x.max(), 200)
    y_line = np.polyval(coeffs, x_line)
    return coeffs[0], coeffs[1], x_line, y_line


def _apply_common_style(ax, title, xlabel, ylabel):
    ax.set_title(f"{title}\n{TITLE_SUFFIX}", fontsize=11, pad=8)
    ax.set_xlabel(xlabel, fontsize=10)
    ax.set_ylabel(ylabel, fontsize=10)
    ax.legend(fontsize=8, framealpha=0.3)
    ax.tick_params(labelsize=9)


# ── Individual plots ──────────────────────────────────────────────────────────

def _plot_excitation(fe: pd.DataFrame, out_dir: Path) -> None:
    """Plot 1 — A(Fe I) vs Excitation Potential."""
    x = fe['excitation_potential_eV'].values
    y = fe['A_Fe_COG'].values
    slope, _, x_line, y_line = _fit_line(x, y)

    with plt.style.context('dark_background'):
        fig, ax = plt.subplots(figsize=(10, 6))
        _scatter_by_regime(ax, x, y, fe['ew_regime'])
        ax.plot(x_line, y_line, '--', color='gold', lw=1.5, label='linear fit', zorder=4)
        ax.text(0.03, 0.06,
                f"slope = {slope:+.3f} dex/eV  |  target: |slope| < 0.02",
                transform=ax.transAxes, fontsize=9, color='gold',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='#333333', alpha=0.7))
        _apply_common_style(ax, 'A(Fe I) vs Excitation Potential',
                            'Excitation Potential (eV)', 'A(Fe I)  [COG]')
        fig.tight_layout()
        fig.savefig(out_dir / 'fe_excitation_potential.png', dpi=150)
        plt.close(fig)

    print(f"  Plot 1 — excitation slope: {slope:+.4f} dex/eV  (target |slope| < 0.02)")


def _plot_reduced_ew(fe: pd.DataFrame, out_dir: Path) -> None:
    """Plot 2 — A(Fe I) vs Reduced EW."""
    x = np.log10(fe['ew_mA'] / fe['wavelength_air_A']).values
    y = fe['A_Fe_COG'].values
    slope, _, x_line, y_line = _fit_line(x, y)

    with plt.style.context('dark_background'):
        fig, ax = plt.subplots(figsize=(10, 6))
        _scatter_by_regime(ax, x, y, fe['ew_regime'])
        ax.plot(x_line, y_line, '--', color='gold', lw=1.5, label='linear fit', zorder=4)
        ax.text(0.03, 0.06,
                f"slope = {slope:+.3f} dex  |  target: flat (|slope| < 0.05)",
                transform=ax.transAxes, fontsize=9, color='gold',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='#333333', alpha=0.7))
        _apply_common_style(ax, 'A(Fe I) vs Reduced EW',
                            'log₁₀(EW / λ)  [reduced EW]', 'A(Fe I)  [COG]')
        fig.tight_layout()
        fig.savefig(out_dir / 'fe_reduced_ew.png', dpi=150)
        plt.close(fig)

    print(f"  Plot 2 — reduced-EW slope: {slope:+.4f}  (target flat)")


def _plot_wavelength(fe: pd.DataFrame, out_dir: Path) -> None:
    """Plot 3 — A(Fe I) vs Wavelength."""
    x = fe['wavelength_air_A'].values
    y = fe['A_Fe_COG'].values
    BLUE_EDGE = 4500.0

    with plt.style.context('dark_background'):
        fig, ax = plt.subplots(figsize=(10, 6))
        _scatter_by_regime(ax, x, y, fe['ew_regime'])
        ax.axvline(BLUE_EDGE, color='tomato', lw=1.2, ls='--', zorder=4)
        ax.axvspan(x.min() - 50, BLUE_EDGE, color='tomato', alpha=0.08,
                   label='blue edge — exclude', zorder=2)
        _apply_common_style(ax, 'A(Fe I) vs Wavelength',
                            'Wavelength (Å)', 'A(Fe I)  [COG]')
        fig.tight_layout()
        fig.savefig(out_dir / 'fe_wavelength.png', dpi=150)
        plt.close(fig)

    print(f"  Plot 3 — wavelength coverage: {x.min():.1f}–{x.max():.1f} Å")


def _plot_histogram(fe: pd.DataFrame, out_dir: Path) -> None:
    """Plot 4 — Per-line A(Fe I) histogram."""
    a_fe = fe['A_Fe_COG'].values
    mean = float(np.mean(a_fe))
    std  = float(np.std(a_fe))
    n    = len(a_fe)

    with plt.style.context('dark_background'):
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.hist(a_fe, bins=40, color='#5b9bd5', edgecolor='#3a7abd', alpha=0.85)
        ax.axvline(ASPLUND_FE, color='gold', lw=1.8, ls='--',
                   label=f'Asplund 2021: A(Fe) = {ASPLUND_FE}')
        ax.axvline(mean, color='white', lw=1.2, ls=':',
                   label=f'This run: {mean:.3f} ± {std:.3f}  (N={n})')
        ax.text(0.97, 0.95,
                f"mean = {mean:.3f}\nstd  = {std:.3f}\nN    = {n}",
                transform=ax.transAxes, fontsize=9, color='white',
                va='top', ha='right', family='monospace',
                bbox=dict(boxstyle='round,pad=0.4', facecolor='#222222', alpha=0.8))
        _apply_common_style(ax, 'Fe I Abundance Distribution',
                            'A(Fe I)  [COG]', 'N lines')
        fig.tight_layout()
        fig.savefig(out_dir / 'fe_histogram.png', dpi=150)
        plt.close(fig)

    print(f"  Plot 4 — histogram: mean={mean:.3f}, std={std:.3f}, N={n}")


# ── Entry point ───────────────────────────────────────────────────────────────

def run() -> pd.DataFrame:
    """Build all four diagnostic plots. Returns the clean Fe I sample used."""
    DIAG_DIR.mkdir(parents=True, exist_ok=True)

    print("── Solar abundance diagnostics (RYA-160) ─────────────────")
    print(f"  Loading {PATHS['solar_ew'].name} …")

    merged = _load_and_merge()
    fe     = _filter_fe(merged)

    if len(fe) == 0:
        raise RuntimeError(
            "No Fe I lines passed quality cuts. "
            "Check solar_ew.csv and linelist_solar.csv are populated."
        )

    print(f"  Fe I lines after cuts: {len(fe)}  "
          f"(from {len(merged[merged['element'] == 'Fe'])} total Fe rows)")

    _plot_excitation(fe, DIAG_DIR)
    _plot_reduced_ew(fe, DIAG_DIR)
    _plot_wavelength(fe, DIAG_DIR)
    _plot_histogram(fe, DIAG_DIR)

    print(f"\n  Plots saved → {DIAG_DIR}/")
    print("──────────────────────────────────────────────────────────")
    return fe


if __name__ == '__main__':
    run()
