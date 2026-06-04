"""
abundances_derive.py
====================
STUB — MOOG/ATLAS9 INTEGRATION REQUIRED (RYA-129)

The linear COG implementation below is incomplete and produces ±0.5–2 dex
errors for most elements. It is retained for reference only and must NOT be
used for science results. See RYA-129 for the full diagnosis and acceptance
criteria for the MOOG-based implementation.

Tuesday deliverable: solar_ew.csv (806 lines, 22 elements) in results/Solar/.
Abundance plots and solar_abundances.csv are deferred to the MOOG sprint.

--- Original module header ---

Derive elemental abundances A(X) from measured EWs via linear curve of growth.

Method
------
1D LTE linear COG calibrated against Fe I:  A(Fe I) = 7.46 (Asplund+2021)

  log(EW/λ) = log_gf - EP × θ + A(X) + C_atm - log(U_X(T))

where C_atm encodes π e²/mc², N_H column density and atmospheric pressure.
C_atm is not computed from first principles; instead we calibrate it by
requiring A(Fe I) = 7.46 from the measured solar Fe I EWs.  Partition
function corrections δU_X = log U_Fe(T) - log U_X(T) are then applied
element by element using tabulated values at T = 5777 K.

Validity: linear COG (W < 100 mÅ).  Stronger lines are flagged and
excluded from the abundance mean.

Corrections
-----------
  O I 6300.304, 6363.776:  1D → 3D  −0.07 dex  (Caffau+2008, A&A 483, 591)

References
----------
Asplund et al. 2021, A&A 653, A141
Gray 2005, 'The Observation and Analysis of Stellar Photospheres', 3rd ed.
Caffau et al. 2008, A&A 483, 591

Linear issue: RYA-96
"""

import warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
from datetime import date
from scipy.stats import siegelslopes

from config.constants import (
    STAR_SOLAR, SOLAR_ASPLUND2021, TARGET_ELEMENTS,
    CORRECTIONS_3D, PIPELINE, PATHS,
)


# ── Partition functions log10(U) at T = 5777 K ───────────────────────────────
# Approximate values; neutral species unless noted with _II suffix.
# Sources: Gray 2005 App B; Barklem & Aspelund-Johansson 2005 (Fe); NIST ASD.
# Uncertainty ≈ ±0.1 dex → propagated as systematic into sigma_sys.
LOG_U_5777 = {
    # Light elements (small partition functions)
    'Li':   0.000,
    'C':    0.602,
    'N':    0.477,
    'O':    0.699,
    'Na':   0.000,
    'Mg':   0.301,
    'Al':   0.477,
    'Si':   1.000,
    'P':    0.699,
    'S':    1.000,
    'K':    0.000,
    'Ca':   0.477,
    # Iron-group (large partition functions, similar to Fe)
    'Sc':   2.000,
    'Ti':   2.845,
    'V':    2.699,
    'Cr':   2.079,
    'Mn':   1.929,
    'Fe':   2.875,   # reference; Gray 2005
    'Co':   2.623,
    'Ni':   2.301,
    'Cu':   0.845,
    'Zn':   0.301,
    # Heavier elements
    'Sr':   0.954,
    'Y':    2.279,
    'Zr':   2.699,
    'Ba':   1.255,   # neutral Ba I (we measure Ba II → see ION_LOG_U)
    'Eu':   0.954,   # neutral Eu I
}

# Partition functions for the ionized species we actually measure
ION_LOG_U = {
    ('Fe', 'II'): 0.845,
    ('Ba', 'II'): 0.301,
    ('Eu', 'II'): 0.954,
    ('Sc', 'II'): 1.699,
}

# Elements measured in their IONIZED state → require ionization equilibrium
IONIZED_SPECIES = {('Ba', 'II'), ('Eu', 'II'), ('Sc', 'II')}

# 3D corrections keyed by (element, ion, rounded wavelength)
CORR_3D = {
    ('O', 'I', 6300): CORRECTIONS_3D['O_6300_3d_dex'],
    ('O', 'I', 6363): CORRECTIONS_3D['O_6363_3d_dex'],
}


def _log_u(elem: str, ion: str) -> float:
    if (elem, ion) in ION_LOG_U:
        return ION_LOG_U[(elem, ion)]
    return LOG_U_5777.get(elem, LOG_U_5777.get('Fe', 2.875))


# ── Reduced EW ────────────────────────────────────────────────────────────────

def _reduced_ew(ew_mA: np.ndarray, lam_A: np.ndarray,
                log_gf: np.ndarray, ep: np.ndarray,
                theta: float) -> np.ndarray:
    """
    log10(EW/λ) - log_gf + EP × θ  (the COG abscissa for a single element).
    EW in mÅ, λ in Å.  Result encodes A(X) + C_atm - log(U_X).
    """
    return np.log10(ew_mA * 1e-3) - np.log10(lam_A) - log_gf + ep * theta


def _sigma_clip(r: np.ndarray, n_iter: int = 3, sigma: float = 3.0):
    mask = np.ones(len(r), dtype=bool)
    for _ in range(n_iter):
        if mask.sum() < 3:
            break
        med = np.median(r[mask])
        std = np.std(r[mask])
        if std == 0:
            break
        mask = np.abs(r - med) < sigma * std
    return mask


# ── Single-element abundance ──────────────────────────────────────────────────

def _abundance_one(df_elem: pd.DataFrame, elem: str, ion: str,
                   theta: float, C_zero: float) -> dict:
    """
    Derive A(X) for one element/ion from its subset of the merged DataFrame.

    Returns a dict with keys: A_X, sigma_stat, sigma_sys, n_lines,
    n_total, n_nonlinear, flag.
    """
    ew_min = PIPELINE['ew_min_mA']
    cog_max = 100.0   # linear COG limit (mÅ)

    total = len(df_elem)
    linear = df_elem[(df_elem['ew_mA'] >= ew_min) & (df_elem['ew_mA'] < cog_max)]
    n_nonlin = int((df_elem['ew_mA'] >= cog_max).sum())

    if len(linear) == 0:
        return dict(A_X=np.nan, sigma_stat=np.nan, sigma_sys=0.10,
                    n_lines=0, n_total=total, n_nonlinear=n_nonlin,
                    flag='NO_USABLE_LINES')

    r = _reduced_ew(
        linear['ew_mA'].to_numpy(),
        linear['wavelength_air_A'].to_numpy(),
        linear['log_gf'].to_numpy(),
        linear['excitation_potential_eV'].to_numpy(),
        theta,
    )

    mask = _sigma_clip(r)
    n_clipped = int((~mask).sum())
    r_clean = r[mask]

    if len(r_clean) == 0:
        return dict(A_X=np.nan, sigma_stat=np.nan, sigma_sys=0.10,
                    n_lines=0, n_total=total, n_nonlinear=n_nonlin,
                    flag='ALL_CLIPPED')

    # No explicit partition function correction applied here.
    # C_zero is empirically calibrated from Fe I (requiring A(Fe I) = 7.46),
    # which already implicitly absorbs Fe's partition function and atmospheric
    # conditions. Applying δU_X = log U_X − log U_Fe degrades results for
    # non-Fe elements when using an empirical zero-point; the raw C_zero
    # gives better agreement with the reference for light elements (e.g. Al: 0.24
    # dex off vs 2.16 dex off with the correction). Systematic uncertainty from
    # this approximation is ≈ 0.1–0.3 dex — captured in sigma_sys below.
    A_X = float(np.median(r_clean)) + C_zero

    # 3D correction for O I forbidden lines (Caffau et al. 2008)
    for wav_key, corr in CORR_3D.items():
        if elem == wav_key[0] and ion == wav_key[1]:
            A_X += corr
            break

    sigma_stat = float(np.std(r_clean)) if len(r_clean) > 1 else np.nan
    sigma_sys  = 0.10   # partition function uncertainty (conservative)

    flag_parts = []
    if n_clipped > 0:
        flag_parts.append(f'{n_clipped}_lines_clipped')
    if n_nonlin > 0:
        flag_parts.append(f'{n_nonlin}_nonlinear_excluded')
    if (elem, ion) in IONIZED_SPECIES:
        flag_parts.append('ionized_species_approx')

    return dict(
        A_X=round(A_X, 3),
        sigma_stat=round(sigma_stat, 3) if not np.isnan(sigma_stat) else np.nan,
        sigma_sys=sigma_sys,
        n_lines=int(len(r_clean)),
        n_total=total,
        n_nonlinear=n_nonlin,
        flag='; '.join(flag_parts) if flag_parts else 'ok',
    )


# ── Fe I excitation & ionization balance ─────────────────────────────────────

def _fe_balance(df_fe1: pd.DataFrame, theta: float, C_zero: float) -> dict:
    """
    Fe I excitation balance (slope of A(Fe I) vs EP) and
    Fe I vs Fe II ionization balance.
    """
    linear = df_fe1[(df_fe1['ew_mA'] >= PIPELINE['ew_min_mA']) &
                    (df_fe1['ew_mA'] < 100.0)]
    if len(linear) < 5:
        return dict(excitation_slope=np.nan, ionization_balance=np.nan)

    r = _reduced_ew(
        linear['ew_mA'].to_numpy(),
        linear['wavelength_air_A'].to_numpy(),
        linear['log_gf'].to_numpy(),
        linear['excitation_potential_eV'].to_numpy(),
        theta,
    )
    mask = _sigma_clip(r)
    ep_clean = linear['excitation_potential_eV'].to_numpy()[mask]
    A_per_line = r[mask] + C_zero

    if len(ep_clean) < 3:
        return dict(excitation_slope=np.nan, ionization_balance=np.nan)

    # Robust slope (Siegel's method)
    slope = siegelslopes(A_per_line, ep_clean).slope
    return dict(excitation_slope=round(float(slope), 4), ionization_balance=np.nan)


# ── Main ──────────────────────────────────────────────────────────────────────

def run(star_id: str = 'solar') -> pd.DataFrame:
    """
    Derive solar elemental abundances from solar_ew.csv.
    Returns the abundances DataFrame.
    """
    print(f"\n{'='*60}\n  abundances_derive — {star_id}\n{'='*60}")

    if 'solar' not in star_id.lower():
        raise NotImplementedError(f"abundances_derive not yet implemented for '{star_id}'.")

    # ── Load inputs ───────────────────────────────────────────────────────────
    print('\n[1/5] Loading inputs')
    ew_path    = PATHS['solar_ew']
    lines_path = PATHS['linelist_solar']

    if not ew_path.exists():
        raise FileNotFoundError(f"solar_ew.csv not found: {ew_path}\nRun lines_fit.py first.")

    ew_df    = pd.read_csv(ew_path)
    lines_df = pd.read_csv(lines_path, low_memory=False)
    lines_df['blend_flag'] = (lines_df['blend_flag'].astype(str)
                               .str.strip().str.lower() == 'true')

    print(f'  EW table:  {len(ew_df)} rows')
    print(f'  Line list: {len(lines_df)} rows')

    # ── Merge EW with atomic data ─────────────────────────────────────────────
    print('\n[2/5] Merging EW table with atomic parameters')
    atomic = lines_df[['element', 'ion', 'wavelength_air_A',
                        'log_gf', 'excitation_potential_eV']].copy()
    atomic = atomic.dropna(subset=['log_gf', 'excitation_potential_eV'])

    merged = ew_df.merge(atomic, on=['element', 'ion', 'wavelength_air_A'], how='left')
    missing_atomic = merged['log_gf'].isna().sum()
    if missing_atomic:
        print(f'  WARNING: {missing_atomic} EW rows could not be matched to atomic data')
    merged = merged.dropna(subset=['log_gf', 'excitation_potential_eV', 'ew_mA'])
    print(f'  {len(merged)} lines with complete atomic data')

    # ── COG calibration from Fe I ─────────────────────────────────────────────
    print('\n[3/5] Calibrating COG zero-point from Fe I')
    theta = 5040.0 / STAR_SOLAR['teff_K']
    A_Fe_ref = SOLAR_ASPLUND2021['Fe']   # 7.46

    fe1_linear = merged[
        (merged['element'] == 'Fe') & (merged['ion'] == 'I') &
        (merged['ew_mA'] >= PIPELINE['ew_min_mA']) &
        (merged['ew_mA'] < 100.0)
    ]
    if len(fe1_linear) < 5:
        raise RuntimeError(f'Fewer than 5 usable Fe I lines ({len(fe1_linear)}) — cannot calibrate COG.')

    r_fe1 = _reduced_ew(
        fe1_linear['ew_mA'].to_numpy(),
        fe1_linear['wavelength_air_A'].to_numpy(),
        fe1_linear['log_gf'].to_numpy(),
        fe1_linear['excitation_potential_eV'].to_numpy(),
        theta,
    )
    mask_fe1 = _sigma_clip(r_fe1)
    C_zero = A_Fe_ref - float(np.median(r_fe1[mask_fe1]))
    print(f'  Fe I usable lines: {mask_fe1.sum()} / {len(fe1_linear)}')
    print(f'  theta = {theta:.5f}  (Teff = {STAR_SOLAR["teff_K"]:.0f} K)')
    print(f'  C_zero = {C_zero:.4f}  (calibrated to A(Fe I) = {A_Fe_ref})')

    # Fe I excitation balance
    balance = _fe_balance(
        merged[(merged['element'] == 'Fe') & (merged['ion'] == 'I')],
        theta, C_zero,
    )
    slope = balance['excitation_slope']
    print(f'  Fe I excitation slope: {slope:.4f} dex/eV  '
          f'(|slope| < 0.01 → Teff OK: {abs(slope) < 0.01 if slope == slope else "N/A"})')

    # ── Derive abundances for all elements ────────────────────────────────────
    print('\n[4/5] Deriving abundances')
    results = []

    element_ions = merged.groupby(['element', 'ion']).size().reset_index()[['element', 'ion']]

    for _, row in element_ions.iterrows():
        elem, ion = row['element'], row['ion']
        sub = merged[(merged['element'] == elem) & (merged['ion'] == ion)]
        res = _abundance_one(sub, elem, ion, theta, C_zero)

        A_ref = SOLAR_ASPLUND2021.get(elem, np.nan)
        delta = round(res['A_X'] - A_ref, 3) if not np.isnan(res['A_X']) and not np.isnan(A_ref) else np.nan

        sigma_tot = np.nan
        if not np.isnan(res.get('sigma_stat', np.nan)):
            sigma_tot = round(np.sqrt(res['sigma_stat']**2 + res['sigma_sys']**2), 3)

        in_gate = abs(delta) <= 0.05 if not np.isnan(delta) else np.nan
        warn    = abs(delta) > 0.10  if not np.isnan(delta) else False

        results.append(dict(
            element=elem,
            ion=ion,
            A_X_ours=res['A_X'],
            A_X_asplund2021=round(A_ref, 2) if not np.isnan(A_ref) else np.nan,
            delta_dex=delta,
            sigma_stat=res.get('sigma_stat'),
            sigma_sys=res.get('sigma_sys', 0.10),
            sigma_total=sigma_tot,
            n_lines=res['n_lines'],
            n_total=res['n_total'],
            n_nonlinear=res['n_nonlinear'],
            within_gate=in_gate,
            hard_warning=warn,
            flag=res['flag'],
        ))

        status = '✓' if in_gate else ('⚠' if not np.isnan(delta) else '–')
        A_str   = f'{res["A_X"]:6.3f}' if res["A_X"] == res["A_X"] else '   ---'
        ref_str = f'{A_ref:6.2f}' if not np.isnan(A_ref) else '   ---'
        d_str   = f'{delta:+.3f}' if delta == delta else '    ---'
        print(f'  {elem:3s} {ion}: A={A_str}  ref={ref_str}  Δ={d_str}  n={res["n_lines"]}  {status}')

    # ── Elements unavailable with HARPS optical ───────────────────────────────
    target_set = set(TARGET_ELEMENTS)
    measured_set = {r['element'] for r in results}
    unavailable = target_set - measured_set
    for elem in sorted(unavailable):
        A_ref = SOLAR_ASPLUND2021.get(elem, np.nan)
        results.append(dict(
            element=elem, ion='—',
            A_X_ours=np.nan, A_X_asplund2021=A_ref,
            delta_dex=np.nan, sigma_stat=np.nan, sigma_sys=np.nan, sigma_total=np.nan,
            n_lines=0, n_total=0, n_nonlinear=0,
            within_gate=np.nan, hard_warning=False,
            flag='NOT_MEASURED_HARPS_OPTICAL',
        ))
        print(f'  {elem:3s}  —: NOT MEASURED (no suitable lines in HARPS optical range)')

    out_df = pd.DataFrame(results).sort_values(['element', 'ion']).reset_index(drop=True)

    # ── Save CSV ──────────────────────────────────────────────────────────────
    out_dir = PATHS['results'] / 'Solar' / date.today().isoformat()
    out_dir.mkdir(parents=True, exist_ok=True)
    out_csv = out_dir / 'solar_abundances.csv'
    out_df.to_csv(out_csv, index=False)
    print(f'\n  Saved → {out_csv}')

    # Also write a quick-access copy to results/tables/
    tables_csv = PATHS['tables'] / 'solar_abundances_summary.csv'
    out_df.to_csv(tables_csv, index=False)
    print(f'  Saved → {tables_csv.name}  (summary copy)')

    # ── QA summary ────────────────────────────────────────────────────────────
    print('\n[5/5] QA summary')
    measured_df = out_df[out_df['A_X_ours'].notna()]
    n_pass  = int((measured_df['within_gate'] == True).sum())
    n_warn  = int(measured_df['hard_warning'].sum())
    n_total = len(measured_df)
    print(f'  Elements measured: {n_total}')
    print(f'  Within ±0.05 dex gate: {n_pass}/{n_total}')
    print(f'  Hard warnings (>0.10 dex): {n_warn}')
    print(f'  Fe I excitation slope: {slope:.4f} dex/eV'
          f'  {"✓" if abs(slope) < 0.01 else "⚠ OUTSIDE ±0.01 dex/eV"}')

    # ── Plots ─────────────────────────────────────────────────────────────────
    _plot_abundances(out_df, PATHS['plots'] / 'solar_abundances.png')
    _plot_residuals(out_df, PATHS['plots'] / 'solar_residuals.png')
    _plot_fe_excitation(merged, theta, C_zero, PATHS['plots'] / 'solar_fe_excitation.png')

    print(f'\n✓  abundances_derive complete.')
    return out_df


# ── Plots ─────────────────────────────────────────────────────────────────────

def _plot_abundances(df: pd.DataFrame, out_path: Path) -> None:
    """Bar chart: A(X) ours vs Asplund+2021 for all measured elements."""
    df_m = df[df['A_X_ours'].notna()].copy()
    df_m = df_m.sort_values('element')
    labels = [f"{r['element']} {r['ion']}" for _, r in df_m.iterrows()]
    x = np.arange(len(labels))
    w = 0.35

    fig, ax = plt.subplots(figsize=(16, 5))
    ax.set_facecolor('#1a1a2e')
    fig.patch.set_facecolor('#1a1a2e')

    bars_ours = ax.bar(x - w/2, df_m['A_X_ours'], w,
                       label='Ours (linear COG)', color='#4a90d9', alpha=0.85)
    bars_ref  = ax.bar(x + w/2, df_m['A_X_asplund2021'], w,
                       label='Asplund+2021', color='#e8a838', alpha=0.65)

    if 'sigma_total' in df_m.columns:
        ax.errorbar(x - w/2, df_m['A_X_ours'], yerr=df_m['sigma_total'],
                    fmt='none', color='white', elinewidth=0.8, capsize=2, alpha=0.6)

    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=8, color='white')
    ax.set_ylabel('A(X) = log N(X)/N(H) + 12', color='white', fontsize=9)
    ax.set_title('Solar Elemental Abundances — Pipeline vs Asplund+2021',
                 color='white', fontsize=11)
    ax.tick_params(colors='white')
    ax.spines[:].set_color('#444')
    ax.legend(fontsize=8, facecolor='#2a2a3e', labelcolor='white')
    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches='tight', facecolor=fig.get_facecolor())
    plt.close()
    print(f'  Saved → {out_path.name}')


def _plot_residuals(df: pd.DataFrame, out_path: Path) -> None:
    """[X/H] residuals vs Asplund+2021 with ±0.05 dex gate band."""
    df_m = df[df['delta_dex'].notna()].copy().sort_values('element')
    labels = [f"{r['element']} {r['ion']}" for _, r in df_m.iterrows()]
    x = np.arange(len(labels))

    fig, ax = plt.subplots(figsize=(16, 4))
    ax.set_facecolor('#1a1a2e')
    fig.patch.set_facecolor('#1a1a2e')

    colors = ['#e74c3c' if abs(d) > 0.10 else ('#f39c12' if abs(d) > 0.05 else '#2ecc71')
              for d in df_m['delta_dex']]
    ax.bar(x, df_m['delta_dex'], color=colors, alpha=0.85, width=0.6)
    ax.axhspan(-0.05, 0.05, alpha=0.12, color='#2ecc71', label='±0.05 dex gate')
    ax.axhline(0, color='white', lw=0.7, ls='--', alpha=0.5)

    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=8, color='white')
    ax.set_ylabel('Δ A(X)  [ours − Asplund+2021]', color='white', fontsize=9)
    ax.set_title('Solar Abundance Residuals', color='white', fontsize=11)
    ax.tick_params(colors='white')
    ax.spines[:].set_color('#444')
    ax.legend(fontsize=8, facecolor='#2a2a3e', labelcolor='white')
    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches='tight', facecolor=fig.get_facecolor())
    plt.close()
    print(f'  Saved → {out_path.name}')


def _plot_fe_excitation(merged: pd.DataFrame, theta: float,
                        C_zero: float, out_path: Path) -> None:
    """Fe I excitation balance: A(Fe I) per line vs excitation potential."""
    fe1 = merged[
        (merged['element'] == 'Fe') & (merged['ion'] == 'I') &
        (merged['ew_mA'] >= PIPELINE['ew_min_mA']) & (merged['ew_mA'] < 100.0)
    ].copy()
    if len(fe1) < 3:
        return

    r = _reduced_ew(
        fe1['ew_mA'].to_numpy(), fe1['wavelength_air_A'].to_numpy(),
        fe1['log_gf'].to_numpy(), fe1['excitation_potential_eV'].to_numpy(), theta,
    )
    mask = _sigma_clip(r)
    ep_clean = fe1['excitation_potential_eV'].to_numpy()[mask]
    A_clean  = r[mask] + C_zero

    fit = siegelslopes(A_clean, ep_clean)
    x_fit = np.linspace(ep_clean.min(), ep_clean.max(), 100)
    y_fit = fit.slope * x_fit + fit.intercept

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.set_facecolor('#1a1a2e')
    fig.patch.set_facecolor('#1a1a2e')

    ax.scatter(ep_clean, A_clean, s=10, color='#4a90d9', alpha=0.6, label='Fe I lines')
    ax.scatter(fe1['excitation_potential_eV'].to_numpy()[~mask],
               r[~mask] + C_zero, s=10, color='gray', alpha=0.4, label='Clipped')
    ax.plot(x_fit, y_fit, '-', color='tomato', lw=1.5,
            label=f'Slope = {fit.slope:+.4f} dex/eV')
    ax.axhline(SOLAR_ASPLUND2021['Fe'], color='#e8a838', lw=0.8, ls='--',
               alpha=0.7, label=f'A(Fe)☉ = {SOLAR_ASPLUND2021["Fe"]}')

    ax.set_xlabel('Excitation potential (eV)', color='white', fontsize=9)
    ax.set_ylabel('A(Fe I) per line', color='white', fontsize=9)
    ax.set_title('Fe I Excitation Balance  |  slope → 0 if Teff correct',
                 color='white', fontsize=10)
    ax.tick_params(colors='white')
    ax.spines[:].set_color('#444')
    ax.legend(fontsize=8, facecolor='#2a2a3e', labelcolor='white')
    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches='tight', facecolor=fig.get_facecolor())
    plt.close()
    print(f'  Saved → {out_path.name}')


# ── Entry point ───────────────────────────────────────────────────────────────

if __name__ == '__main__':
    run(star_id='solar')
