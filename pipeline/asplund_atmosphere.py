"""RYA-1315 -- Asplund, Amarsi & Grevesse 2021 Sect. 2.1's two atmosphere terms, per line.

    temperature structure = 1/2 |mean(<3D> - HM74)|     inhomogeneities = 1/2 |mean(3D - <3D>)|

read from Scott+2015 I/II's per-line LTE abundances in the 3D, <3D> and Holweger & Mueller
(1974) models (data/reference/scott2015). Used by the budget assembler, the live-product
re-stamp and the legacy 3D-NLTE route, so all three price the same thing. Atomic C, N, O (not
in Scott) take Asplund+2021 Table 3's own per-indicator <3D>/HM/3D results.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SCOTT2015 = ROOT / "data/reference/scott2015/scott2015_line_models.csv"
#: Asplund+2021 Table 3's own <3D>/HM/3D results for the atomic C, N, O indicators, which Scott
#: does not tabulate (non-LTE columns, Asplund's adopted ones).
ASPLUND2021_CNO = ROOT / "data/reference/asplund2021_cno/asplund2021_table3_atomic.csv"
#: the forbidden lines Table 3 lists separately: [C I] 872.7 nm, [O I] 630.0 / 636.4 nm
FORBIDDEN_A = {"C": (8727.13,), "O": (6300.30, 6363.78)}


def scott_atmosphere(element: str, ion: str, pool_w) -> dict | None:
    """Asplund+2021 Sect. 2.1's two atmosphere terms for a pool, from Scott+2015 I/II's per-line
    LTE abundances in the 3D, <3D> and Holweger & Mueller (1974) models:
        temperature structure = 1/2 |mean(<3D> - HM)|,  inhomogeneities = 1/2 |mean(3D - <3D>)|
    averaged over the pool's lines Scott tabulates (0.05 A); a pool with none of them takes the
    species' own Scott lines (Asplund's per-element value). None when Scott has no such species."""
    t = pd.read_csv(SCOTT2015)
    t = t[(t.element == element) & (t.ion == ion)]
    if t.empty:
        return asplund2021_cno_atmosphere(element, ion, pool_w)
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


def asplund2021_cno_atmosphere(element: str, ion: str, pool_w) -> dict | None:
    """The same two terms for atomic C, N, O from Asplund+2021 Table 3 itself (per indicator,
    non-LTE): each pool line takes its indicator's row ([C I]/[O I] for the forbidden lines,
    within 0.05 A; else C I/N I/O I), averaged over the pool as for Scott. None otherwise."""
    t = pd.read_csv(ASPLUND2021_CNO)
    t = t[(t.element == element) & (t.ion == ion)]
    pool = [float(w) for w in pool_w]
    if t.empty or not pool:
        return None
    forb = FORBIDDEN_A.get(element, ())
    rows = [t[t.forbidden == int(any(abs(w - f) <= 0.05 for f in forb))].iloc[0] for w in pool]
    temp = float(np.mean([r.nlte_mean3D - r.nlte_HM for r in rows]))
    inhom = float(np.mean([r.nlte_3D - r.nlte_mean3D for r in rows]))
    return {"temperature_structure_dex": abs(temp) / 2.0, "inhomogeneity_dex": abs(inhom) / 2.0,
            "mean_3Dmean_minus_HM_dex": temp, "mean_3D_minus_3Dmean_dex": inhom,
            "n_scott_lines": len(rows),
            "basis": "Asplund+2021 Table 3 non-LTE, per indicator (" +
                     ", ".join(sorted({r.indicator for r in rows})) + ")",
            "scott_wavelengths_A": sorted(round(w, 3) for w in pool)}
