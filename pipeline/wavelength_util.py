"""
pipeline/wavelength_util.py
===========================
The **one** vacuum↔air wavelength converter in the codebase (SSOT) — Birch & Downs
(1994), as adopted by VALD3 / NIST / iSpec, which is the standard the project's air
line lists are on.

Promoted here from RYA-426 `uv_conditioning` (by RYA-264) so a single implementation
backs **both** the UV conditioning path (`uv_conditioning.to_pipeline_frame`) and the
loader wavelength-scale axis (`frame_object_contract.to_air_angstrom`). A second copy
of a vac↔air formula is exactly the SSOT defect this project forbids — import these,
never re-derive Morton/Birch-Downs locally.

Convention boundary: the air refractive index is only defined for λ ≥ 2000 Å (IAU
1991; Morton 2000). Below that, "air" is undefined and the VALD/IAU convention keeps
wavelengths in **vacuum** (the FUV-stays-vacuum boundary, RYA-303/426). So both
converters are the **identity** below 2000 Å and only transform λ ≥ 2000 Å. (The IR
instruments this serves — SPIRou, NIRPS, CRIRES+ — are all far above the boundary, so
they always convert.)

Formula: n(λ_vac) = 1 + 8.34254e-5 + 2.406147e-2/(130 − s²) + 1.5998e-4/(38.9 − s²),
with s = 1e4/λ_vac(Å) in µm⁻¹. λ_air = λ_vac / n.  Reference: Birch & Downs 1994,
Metrologia 31, 315 (revising Edlén 1966).
"""
from __future__ import annotations

import numpy as np

# Below this wavelength air is undefined → vac↔air is the identity (IAU/VALD).
AIR_VACUUM_BOUNDARY_A = 2000.0


def _bd1994_n(wave_vac_A) -> np.ndarray:
    """Birch & Downs (1994) refractive index of standard air at vacuum wavelength λ (Å)."""
    s2 = (1.0e4 / np.asarray(wave_vac_A, float)) ** 2          # (1/λ_µm)²
    return (1.0 + 8.34254e-5 + 2.406147e-2 / (130.0 - s2)
            + 1.5998e-4 / (38.9 - s2))


def vac_to_air(wave_vac_A) -> np.ndarray:
    """Vacuum → air (Birch & Downs 1994). Identity below 2000 Å (air undefined there;
    VALD/IAU keep vacuum). λ_air = λ_vac / n(λ_vac)."""
    w = np.asarray(wave_vac_A, float)
    out = w.copy()
    conv = w >= AIR_VACUUM_BOUNDARY_A
    out[conv] = w[conv] / _bd1994_n(w[conv])
    return out


def air_to_vac(wave_air_A) -> np.ndarray:
    """Air → vacuum inverse (one Birch & Downs step; n varies negligibly over the
    air↔vac offset, so evaluating n at the air wavelength is accurate to ≪1 mÅ).
    Identity below 2000 Å. λ_vac = λ_air · n(λ_air)."""
    w = np.asarray(wave_air_A, float)
    out = w.copy()
    conv = w >= AIR_VACUUM_BOUNDARY_A
    out[conv] = w[conv] * _bd1994_n(w[conv])
    return out
