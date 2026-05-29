"""
config/constants.py
===================
Single source of truth for all constants, paths, and parameters.
NO hardcoded values anywhere else in the pipeline.

Physical constants: NIST CODATA 2022 recommended values
Stellar parameters: CHARA interferometry (von Braun et al. 2011)
Solar reference:    Asplund et al. 2021 (A&A 653, A141)
"""

__version__ = "2026-05-28"

from pathlib import Path

# ── Repository root ───────────────────────────────────────────────────────────
ROOT = Path(__file__).parent.parent

# ── Physical constants (NIST CODATA 2022) ────────────────────────────────────
PHYSICS = {
    'c_kms'        : 299792.458,        # Speed of light (km/s) — exact
    'c_ms'         : 299792458.0,       # Speed of light (m/s) — exact
    'h_eV'         : 4.135667696e-15,   # Planck constant (eV·s)
    'k_eV'         : 8.617333262e-5,    # Boltzmann constant (eV/K)
    'k_cgs'        : 1.380649e-16,      # Boltzmann constant (erg/K)
    'sigma_sb'     : 5.670374419e-5,    # Stefan-Boltzmann (erg/cm²/s/K⁴)
    'G_cgs'        : 6.67430e-8,        # Gravitational constant (cm³/g/s²)
    'amu_g'        : 1.66053906660e-24, # Atomic mass unit (g)
}

# ── Astronomical constants ────────────────────────────────────────────────────
ASTRO = {
    'Msun_g'       : 1.989e33,          # Solar mass (g)
    'Rsun_cm'      : 6.957e10,          # Solar radius (cm)
    'Lsun_erg'     : 3.828e33,          # Solar luminosity (erg/s)
    'Teff_sun'     : 5778.0,            # Solar effective temperature (K)
    'logg_sun'     : 4.438,             # Solar surface gravity (cgs)
    'pc_cm'        : 3.085677581e18,    # Parsec (cm)
    'AU_cm'        : 1.495978707e13,    # Astronomical unit (cm)
}

# ── Solar abundances: Asplund et al. 2021 (A&A 653, A141) ────────────────────
# A(X) = log(N_X/N_H) + 12
# This is the reference — every [X/H] in the Codex is defined as A(X)_star - A(X)_sun
SOLAR_ASPLUND2021 = {
    'H'  : 12.00,   # by definition
    'He' : 10.914,
    'Li' :  1.05,
    'C'  :  8.46,
    'N'  :  7.83,
    'O'  :  8.69,
    'Na' :  6.24,
    'Mg' :  7.55,
    'Al' :  6.43,
    'Si' :  7.51,
    'P'  :  5.41,
    'S'  :  7.12,
    'K'  :  5.07,
    'Ca' :  6.30,
    'Sc' :  3.14,
    'Ti' :  4.97,
    'V'  :  3.90,
    'Cr' :  5.62,
    'Mn' :  5.42,
    'Fe' :  7.46,   # NOTE: was 7.50 in Lodders 2003 used in 2010 thesis
    'Co' :  4.94,
    'Ni' :  6.20,
    'Ba' :  2.27,
    'Y'  :  2.21,
    'Zr' :  2.59,
    'Eu' :  0.52,
}

# ── 55 Cancri A stellar parameters ───────────────────────────────────────────
# Primary source: CHARA interferometry, von Braun et al. 2011 (ApJ 729, 63)
# Secondary: Valenti & Fischer 2005, Teske et al. 2013
STAR_55CNC = {
    'name'         : '55 Cancri A',
    'iau_name'     : 'Copernicus',
    'hd'           : 'HD 75732',
    'hip'          : 'HIP 43587',
    'teff_K'       : 5196.0,    # ± 24 K (CHARA)
    'teff_err'     : 24.0,
    'logg'         : 4.41,      # ± 0.02 (CHARA)
    'logg_err'     : 0.02,
    'feh'          : 0.32,      # ± 0.02 (spectroscopic)
    'feh_err'      : 0.02,
    'vturb_kms'    : 0.9,       # km/s initial estimate — refined from Fe lines
    'mass_msun'    : 0.905,     # ± 0.015 (CHARA)
    'radius_rsun'  : 0.943,     # ± 0.010 (CHARA)
    'age_gyr'      : 10.2,      # ± 2.5
    'distance_pc'  : 12.55,
    'vmag'         : 5.95,
    'vsini_kms'    : 2.3,       # Very slow rotator — no rotational broadening issues
    'rv_kms'       : 27.58,     # Systemic radial velocity — for Doppler correction
}

# ── Pipeline settings ─────────────────────────────────────────────────────────
PIPELINE = {
    # Wavelength range (full HARPS optical range)
    'wav_min_A'         : 3780.0,
    'wav_max_A'         : 6910.0,

    # S/N requirements
    'snr_min_science'   : 200,    # Minimum S/N for abundance work
    'snr_target'        : 400,    # Target S/N after co-adding

    # Gaussian fitting
    'fit_window_A'      : 2.0,    # ± Angstroms around line center for fitting
    'voigt_threshold_mA': 150,    # EW threshold above which to use Voigt profile

    # Continuum normalization
    'continuum_poly_deg': 3,      # Polynomial degree for continuum fit
    'sigma_clip'        : 3.0,    # Sigma clipping for continuum anchor points

    # Line quality thresholds
    'ew_min_mA'         : 5.0,    # Minimum detectable EW (mÅ)
    'ew_max_mA'         : 300.0,  # Maximum EW before saturation concerns
    'snr_line_min'      : 5.0,    # Minimum line S/N for reliable EW

    # NIST grade policy
    'min_nist_grade'    : 'B',    # Minimum acceptable grade for science
}

# ── File paths ────────────────────────────────────────────────────────────────
PATHS = {
    'data_root'         : ROOT / 'data',
    'raw_spectra'       : ROOT / 'data' / 'raw',
    'processed_spectra' : ROOT / 'data' / 'processed',
    'linelists'         : ROOT / 'data' / 'linelists',
    'model_atmospheres' : ROOT / 'data' / 'model_atmospheres',
    'results'           : ROOT / 'results',
    'plots'             : ROOT / 'results' / 'plots',
    'tables'            : ROOT / 'results' / 'tables',

    # Key data files
    'linelist_master'   : ROOT / 'data' / 'linelists' / 'linelist_master.csv',
    'solar_reference'   : ROOT / 'data' / 'linelists' / 'solar_asplund2021.csv',
    'atlas9_grid'       : ROOT / 'data' / 'model_atmospheres' / 'atlas9',
}

# Auto-create all output directories on import
for _key in ('data_root', 'raw_spectra', 'processed_spectra', 'linelists',
             'model_atmospheres', 'results', 'plots', 'tables'):
    PATHS[_key].mkdir(parents=True, exist_ok=True)

# ── Model atmosphere parameters ───────────────────────────────────────────────
MODEL = {
    'type'              : 'ATLAS9',          # Castelli-Kurucz
    'geometry'          : 'plane-parallel',  # Valid for dwarfs (log g > 3.5)
    'radiative_transfer': 'LTE',
    'grid_source'       : 'https://www.stsci.edu/hst/instrumentation/reference-data-for-calibration-and-tools/astronomical-catalogs/castelli-and-kurucz-atlas',
    'nlte_corrections'  : 'Amarsi et al. 2020',  # Applied post-LTE for O, Na, C
}


# ── Validation ────────────────────────────────────────────────────────────────

_REQUIRED_KEYS = {
    'PHYSICS'          : ['c_kms', 'c_ms', 'h_eV', 'k_eV', 'k_cgs', 'sigma_sb', 'G_cgs', 'amu_g'],
    'ASTRO'            : ['Msun_g', 'Rsun_cm', 'Lsun_erg', 'Teff_sun', 'logg_sun', 'pc_cm', 'AU_cm'],
    'SOLAR_ASPLUND2021': ['H', 'C', 'N', 'O', 'Na', 'Mg', 'Al', 'Si', 'Ca', 'Fe', 'Ni'],
    'STAR_55CNC'       : ['teff_K', 'logg', 'feh', 'vturb_kms', 'rv_kms', 'hd', 'hip'],
    'PIPELINE'         : ['wav_min_A', 'wav_max_A', 'snr_min_science', 'ew_min_mA', 'ew_max_mA'],
    'PATHS'            : ['data_root', 'raw_spectra', 'results', 'plots'],
}

_DICTS = {
    'PHYSICS'          : PHYSICS,
    'ASTRO'            : ASTRO,
    'SOLAR_ASPLUND2021': SOLAR_ASPLUND2021,
    'STAR_55CNC'       : STAR_55CNC,
    'PIPELINE'         : PIPELINE,
    'PATHS'            : PATHS,
}


def validate_constants():
    """Check all required keys exist and critical values are sane. Raises on failure."""
    missing = []
    for dict_name, required_keys in _REQUIRED_KEYS.items():
        d = _DICTS[dict_name]
        for key in required_keys:
            if key not in d:
                missing.append(f"{dict_name}['{key}']")
    if missing:
        raise KeyError(f"Missing constants: {', '.join(missing)}")

    assert abs(PHYSICS['c_kms'] - 299792.458) < 0.001, "Speed of light mismatch"
    assert abs(SOLAR_ASPLUND2021['Fe'] - 7.46) < 0.01, "Solar Fe abundance mismatch"
    assert 4000 < STAR_55CNC['teff_K'] < 7000, "Teff out of range"

    print("Constants validated ✓")
