"""
data/linelists/loader.py
========================
Loads and filters the NIST/VALD-traceable master line list.
All pipeline scripts import from here — never hardcode wavelengths.

Usage:
    from data.linelists.loader import load_linelist, get_element_lines

    # Load all Priority 1 lines for abundance work
    lines = load_linelist(elements=['Fe', 'Ca', 'Mg', 'O', 'C'], min_nist_grade='B')

    # Get Fe I lines only
    fe_lines = get_element_lines('Fe', ion='I')
"""

import pandas as pd
import numpy as np
from pathlib import Path
from config.constants import PATHS, PIPELINE

# NIST grade ordering — lower number = better quality
GRADE_ORDER = {'A+': 0, 'A': 1, 'B': 2, 'C': 3, 'D': 4, 'E': 5}


def load_linelist(
    elements=None,
    ion=None,
    min_nist_grade=None,
    wav_min=None,
    wav_max=None,
    exclude_blends=True,
    priority=None,
) -> pd.DataFrame:
    """
    Load the master line list with optional filtering.

    Parameters
    ----------
    elements      : list of str or None — e.g. ['Fe', 'Ca', 'Mg']
                    None = return all elements
    ion           : str or None — 'I' for neutral, 'II' for singly ionized
    min_nist_grade: str or None — minimum quality grade ('A+', 'A', 'B', 'C')
                    None = use PIPELINE['min_nist_grade'] from constants.py
    wav_min       : float or None — minimum wavelength (Å); None = PIPELINE default
    wav_max       : float or None — maximum wavelength (Å); None = PIPELINE default
    exclude_blends: bool — if True, exclude lines flagged as blended
    priority      : int or None — filter to lines with priority <= this value

    Returns
    -------
    df : pandas DataFrame with columns:
         element, ion, wavelength_air_A, excitation_potential_eV,
         log_gf, loggf_source, nist_grade, blend_flag, priority, notes
    """
    linelist_path = PATHS['linelist_master']
    if not linelist_path.exists():
        raise FileNotFoundError(
            f"Master line list not found at {linelist_path}\n"
            f"Run RYA-64 to generate it from VALD3 + NIST."
        )

    df = pd.read_csv(linelist_path, comment='#')


# ── Canonical line reliability model ─────────────────────────────

base_blend = df['blend_flag'].astype(float)

depth_penalty = 0.0
if 'central_depth' in df.columns:
    depth_penalty = np.clip(df['central_depth'], 0, 1)

# NEW: unified blend score (0 clean → 1 unreliable)
df['blend_score'] = np.clip(
    base_blend * 0.7 + depth_penalty * 0.3,
    0.0,
    1.0
)

# Optional interpretability split (useful later)
df['crowding_score'] = depth_penalty
df['catalog_blend']  = base_blend


    # ── LINE FEATURE ENRICHMENT (catalog-level) ─────────────────────────

# Simple physics-informed blend score
# 0 = clean, 1 = heavily blended / unreliable
df['blend_score'] = np.where(df['blend_flag'], 1.0, 0.0)

# Optional refinement using central depth (proxy for crowding)
if 'central_depth' in df.columns:
    df['blend_score'] += np.clip(df['central_depth'], 0, 1) * 0.3

df['blend_score'] = np.clip(df['blend_score'], 0.0, 1.0)

#end new add by Ryan
    # Wavelength range
    wmin = wav_min if wav_min is not None else PIPELINE['wav_min_A']
    wmax = wav_max if wav_max is not None else PIPELINE['wav_max_A']
    df = df[(df['wavelength_air_A'] >= wmin) & (df['wavelength_air_A'] <= wmax)]

    # Element filter
    if elements is not None:
        df = df[df['element'].isin(elements)]

    # Ionization state filter
    if ion is not None:
        df = df[df['ion'] == ion]

    # NIST grade filter
    grade = min_nist_grade if min_nist_grade is not None else PIPELINE['min_nist_grade']
    max_grade_num = GRADE_ORDER.get(grade, 2)
    df = df[df['nist_grade'].map(lambda g: GRADE_ORDER.get(g, 99)) <= max_grade_num]

    # Blend filter
    if exclude_blends:
        df = df[df['blend_flag'] == False]

   # Priority filter
if priority is not None:
    df = df[df['priority'] <= priority]

# ─────────────────────────────────────────────
# LINE PHYSICS ENRICHMENT LAYER (STEP 1)
# ─────────────────────────────────────────────

# Saturation penalty (continuous, not binary)
if 'central_depth' in df.columns:
    depth = df['central_depth'].clip(0, 1)

    # soft exponential penalty past ~0.6
    df['saturation_penalty'] = np.exp(-5 * (depth - 0.6).clip(lower=0))
else:
    df['saturation_penalty'] = 1.0
# Reduced EW placeholder (will be filled later after EW step)
df['reduced_ew'] = np.nan

# Physics-based quality weight
grade_weight = df['nist_grade'].map({
    'A+': 1.0,
    'A': 0.9,
    'B': 0.75,
    'C': 0.5,
    'D': 0.2
}).fillna(0.1)

#blend_penalty = df['blend_flag'].map({False: 1.0, True: 0.3})
blend_penalty = np.exp(-df['blend_score'])
#blend_penalty = exp(-blend_score)
#blend_score =
   # Σ (neighbor_line_strength × overlap_function(Δλ))

priority_weight = df['priority'].map({
    1: 1.0,
    2: 0.7,
    3: 0.4
}).fillna(0.5)

df['line_quality_weight'] = (grade_weight * blend_penalty * priority_weight * df[df['saturation_penalty'])

return df.reset_index(drop=True)

    return df.reset_index(drop=True)


def get_element_lines(element: str, ion: str = 'I', **kwargs) -> pd.DataFrame:
    """Convenience wrapper — get all lines for one element and ionization state."""
    return load_linelist(elements=[element], ion=ion, **kwargs)


def get_fe_excitation_range(df: pd.DataFrame = None) -> pd.DataFrame:
    """
    Return Fe I lines sorted by excitation potential.

    Used for the Teff self-consistency check: Fe I abundance must NOT
    correlate with excitation potential. Need lines spanning ~0.5–5.0 eV
    to get a reliable slope for the excitation equilibrium test.
    """
    if df is None:
        df = load_linelist(elements=['Fe'], ion='I', min_nist_grade='B')
    return df.sort_values('excitation_potential_eV').reset_index(drop=True)


def summarize_linelist(df: pd.DataFrame) -> None:
    """Print a QA summary of a loaded line list."""
    print(f"\n── Line list summary ───────────────────────────────────")
    print(f"  Total lines : {len(df)}")
    print(f"  Elements    : {sorted(df['element'].unique())}")
    print(f"  Wav range   : {df['wavelength_air_A'].min():.1f} – "
          f"{df['wavelength_air_A'].max():.1f} Å")
    print(f"  NIST grades : {df['nist_grade'].value_counts().to_dict()}")
    print(f"  Priorities  : {df['priority'].value_counts().to_dict()}")
    print(f"─────────────────────────────────────────────────────\n")
