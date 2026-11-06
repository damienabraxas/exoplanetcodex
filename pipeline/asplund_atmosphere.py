"""RYA-1315 -- Asplund, Amarsi & Grevesse 2021 Sect. 2.1's two atmosphere terms, per line.

    temperature structure = 1/2 |mean(<3D> - HM74)|     inhomogeneities = 1/2 |mean(3D - <3D>)|

read from Scott+2015 I/II's per-line LTE abundances in the 3D, <3D> and Holweger & Mueller
(1974) models (data/reference/scott2015). Used by the budget assembler, the live-product
re-stamp and the legacy 3D-NLTE route, so all three price the same thing.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SCOTT2015 = ROOT / "data/reference/scott2015/scott2015_line_models.csv"


def scott_atmosphere(element: str, ion: str, pool_w) -> dict | None:
    """Asplund+2021 Sect. 2.1's two atmosphere terms for a pool, from Scott+2015 I/II's per-line
    LTE abundances in the 3D, <3D> and Holweger & Mueller (1974) models:
        temperature structure = 1/2 |mean(<3D> - HM)|,  inhomogeneities = 1/2 |mean(3D - <3D>)|
    averaged over the pool's lines Scott tabulates (0.05 A); a pool with none of them takes the
    species' own Scott lines (Asplund's per-element value). None when Scott has no such species."""
    t = pd.read_csv(SCOTT2015)
    t = t[(t.element == element) & (t.ion == ion)]
    if t.empty:
        return None
    wa = t.wavelength_nm.astype(float) * 10.0
    pool = [float(w) for w in pool_w]
    hit = t[[any(abs(w - x) <= 0.05 for x in pool) for w in wa]]
    use, basis = (hit, "pool lines") if len(hit) else (t, "species mean (no pool line in Scott)")
    temp = float(np.mean(use.A_mean3D - use.A_HM))
    inhom = float(np.mean(use.A_3D - use.A_mean3D))
    return {"temperature_structure_dex": abs(temp) / 2.0, "inhomogeneity_dex": abs(inhom) / 2.0,
            "mean_3Dmean_minus_HM_dex": temp, "mean_3D_minus_3Dmean_dex": inhom,
            "n_scott_lines": int(len(use)), "basis": basis,
            "scott_wavelengths_A": sorted(round(float(w) * 10.0, 3) for w in use.wavelength_nm)}
