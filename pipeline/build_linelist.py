"""
pipeline/build_linelist.py
==========================
Rebuild linelist CSVs with the RYA-209 column schema:

  vald_proximity_flag (float 0.0–1.0)
      Continuous proximity contamination estimate, computed from VALD neighbour
      density and relative central depths. Replaces the old proximity-based
      binary blend_flag. Used as an input to the RYA-220 quality scoring system.

  blend_flag (bool, default False)
      Vetted spectroscopic exclusion flag. Set True only when a line is confirmed
      non-separable at HARPS R~115,000 via synthesis or literature. This is the
      authoritative source consumed by abundances_derive.py for hard exclusions.

Run this script whenever adding new vetted exclusions or rebuilding from VALD.

RYA-209
"""

import numpy as np
import pandas as pd
from pathlib import Path

from config.constants import PATHS

# Vetted spectroscopic exclusions — RYA-208 confirmed blends at HARPS R~115,000.
# Add new entries here (with literature reference in the comment) as lines are vetted.
VETTED_BLENDS = [
    ('Fe', 'I', 4918.994),  # RYA-208: confirmed non-separable blend
    ('Fe', 'I', 4970.496),  # RYA-208: confirmed non-separable blend
]
VETTING_TOLERANCE = 0.05  # Å — tight; these are exact vetted wavelengths

# vald_proximity_flag formula parameters (RYA-209 spec)
PROXIMITY_WINDOW = 0.5  # Å — neighbours within this window contribute
PROXIMITY_SCALE  = 0.1  # Å — exponential decay scale


def compute_vald_proximity_flag(df: pd.DataFrame) -> np.ndarray:
    """
    Compute continuous proximity contamination score for each line.

    Formula (RYA-209):
        penalty_i = Σ_j  min(depth_j / depth_i, 1.0) * exp(-|Δλ| / 0.1)
        vald_proximity_flag_i = min(penalty_i, 1.0)

    where j ranges over all neighbours within PROXIMITY_WINDOW Å of line i.
    Uses scipy.spatial.cKDTree for O(n + k) pair enumeration.
    """
    try:
        from scipy.spatial import cKDTree
    except ImportError:
        raise ImportError("scipy is required for vald_proximity_flag computation: pip install scipy")

    wl    = df['wavelength_air_A'].values.astype(float)
    depth = df['central_depth'].values.astype(float)
    n     = len(wl)

    tree  = cKDTree(wl.reshape(-1, 1))
    pairs = tree.query_pairs(PROXIMITY_WINDOW, output_type='ndarray')  # shape (k, 2)

    penalty = np.zeros(n)
    if len(pairs) > 0:
        i_idx, j_idx = pairs[:, 0], pairs[:, 1]
        delta = np.abs(wl[i_idx] - wl[j_idx])

        # Contribution of j to i: depth_j / depth_i
        contrib_ji_to_i = (
            np.minimum(depth[j_idx] / np.maximum(depth[i_idx], 1e-9), 1.0)
            * np.exp(-delta / PROXIMITY_SCALE)
        )
        # Contribution of i to j: depth_i / depth_j
        contrib_ij_to_j = (
            np.minimum(depth[i_idx] / np.maximum(depth[j_idx], 1e-9), 1.0)
            * np.exp(-delta / PROXIMITY_SCALE)
        )

        np.add.at(penalty, i_idx, contrib_ji_to_i)
        np.add.at(penalty, j_idx, contrib_ij_to_j)

    return np.minimum(penalty, 1.0)


def build_vetted_blend_flag(df: pd.DataFrame) -> pd.Series:
    """
    Build the new vetted blend_flag column.
    Default False; set True for lines in VETTED_BLENDS (within VETTING_TOLERANCE Å).
    """
    flags = pd.Series(False, index=df.index, dtype=bool)
    for elem, ion, wl in VETTED_BLENDS:
        mask = (
            (df['element'] == elem) &
            (df['ion'] == ion) &
            (np.abs(df['wavelength_air_A'] - wl) < VETTING_TOLERANCE)
        )
        n_matched = int(mask.sum())
        if n_matched > 0:
            flags[mask] = True
            print(f"    Vetted blend: {elem} {ion} {wl:.3f} Å — {n_matched} row(s) flagged")
        else:
            print(f"    WARNING: vetted blend {elem} {ion} {wl:.3f} Å NOT FOUND in {df.attrs.get('name', 'linelist')}")
    return flags


def rebuild_linelist(path: Path) -> None:
    print(f"\n  Rebuilding {path.name} ({path})...")
    df = pd.read_csv(str(path), low_memory=False)
    df.attrs['name'] = path.name
    print(f"    Loaded {len(df)} rows, columns: {list(df.columns)}")

    # Validate required columns
    for col in ('wavelength_air_A', 'central_depth'):
        if col not in df.columns:
            raise ValueError(f"{path.name} missing required column '{col}'")

    # Compute vald_proximity_flag
    print(f"    Computing vald_proximity_flag (window={PROXIMITY_WINDOW} Å, scale={PROXIMITY_SCALE} Å)...")
    vpf = compute_vald_proximity_flag(df)
    print(f"    vald_proximity_flag: mean={vpf.mean():.4f}  max={vpf.max():.4f}  "
          f">=0.5: {(vpf >= 0.5).sum()}  >=0.9: {(vpf >= 0.9).sum()}")

    # Build vetted blend_flag
    print(f"    Building vetted blend_flag...")
    new_blend_flag = build_vetted_blend_flag(df)
    print(f"    blend_flag True={new_blend_flag.sum()}  (all others False)")

    # Drop old blend_flag (proximity-based bool), insert new columns in its place
    old_cols = list(df.columns)
    if 'blend_flag' in old_cols:
        insert_pos = old_cols.index('blend_flag')
        df = df.drop(columns=['blend_flag'])
    else:
        insert_pos = old_cols.index('central_depth') + 1

    df.insert(insert_pos,     'vald_proximity_flag', np.round(vpf, 6))
    df.insert(insert_pos + 1, 'blend_flag',          new_blend_flag)

    df.to_csv(str(path), index=False)
    print(f"    Saved → {path.name}  columns: {list(df.columns)}")


def run():
    print(f"\n{'='*60}")
    print(f"  build_linelist  |  RYA-209 schema migration")
    print(f"{'='*60}")

    targets = [
        Path(str(PATHS['linelist_solar'])),
        Path(str(PATHS['linelist_master'])),
    ]
    for path in targets:
        if not path.exists():
            print(f"\n  SKIP: {path} not found")
            continue
        rebuild_linelist(path)

    print(f"\n{'='*60}")
    print(f"  build_linelist complete.")
    print(f"{'='*60}\n")


if __name__ == '__main__':
    run()
