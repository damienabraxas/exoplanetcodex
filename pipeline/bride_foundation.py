"""Read-only boundary checks for experimental Bride columns; no RT or abundance solve.

The Multi3D reader returns SI arrays in (z, y, x) order. Lightweaver 0.16.1
``Atmosphere.make_1d`` expects decreasing geometric height and logG in m/s^2,
whereas STAR_PARAMS logg is cgs. See the installed atmosphere.py docstring and
https://goobley.github.io/Lightweaver/lightweaver-api.html . These checks establish
an input contract, not a validated atmosphere, atom, solver, or science product.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

import numpy as np


class FoundationInputError(ValueError):
    """A required input is invalid; never repair, clip, or infer physics here."""


def _finite_array(value, name: str, shape: tuple[int, ...]) -> np.ndarray:
    array = np.asarray(value, dtype=np.float64)
    if array.shape != shape or not np.isfinite(array).all():
        raise FoundationInputError(f"{name}: expected finite array of shape {shape}")
    return array


def validate_snapshot(snapshot: Mapping) -> tuple[int, int, int]:
    """Validate shapes, coordinates, and signs, without asserting solar realism.

    x/y increase, z decreases (top to bottom); never reorder a cube to make it
    pass. Hydrogen entries are the reader's neutral and proton densities, not
    a multi-level NLTE hydrogen population vector. Zero in either is allowed;
    their sum must be positive. All number densities are m^-3, velocities m/s.
    """
    sizes = []
    for name in ("nz", "ny", "nx"):
        value = snapshot[name]
        if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)) or value < 1:
            raise FoundationInputError(f"{name}: expected a positive integer")
        sizes.append(int(value))
    shape = tuple(sizes)
    if shape[0] < 2:
        raise FoundationInputError("z: a column needs at least two depth points")
    for name, length, direction in (("z", shape[0], -1), ("y", shape[1], 1), ("x", shape[2], 1)):
        axis = _finite_array(snapshot[name], name, (length,))
        if not np.all(direction * np.diff(axis) > 0):
            raise FoundationInputError(f"{name}: coordinates are not strictly ordered")
    for name in ("temperature", "electron_density", "hydrogen_density", "proton_density", "vx", "vy", "vz"):
        array = _finite_array(snapshot[name], name, shape)
        if name in ("temperature", "electron_density") and np.any(array <= 0):
            raise FoundationInputError(f"{name}: must be positive")
        if name in ("hydrogen_density", "proton_density") and np.any(array < 0):
            raise FoundationInputError(f"{name}: must be nonnegative")
    total = np.asarray(snapshot["hydrogen_density"]) + np.asarray(snapshot["proton_density"])
    if not np.isfinite(total).all() or np.any(total <= 0):
        raise FoundationInputError("total hydrogen density: must be finite and positive")
    return shape


@dataclass(frozen=True)
class ColumnInputs:
    """Detached SI arrays plus explicit source identity for one column, not a fit."""

    star_id: str
    ix: int
    iy: int
    height_m: np.ndarray
    temperature_K: np.ndarray
    ne_m3: np.ndarray
    nH_total_m3: np.ndarray
    vx_m_s: np.ndarray
    vy_m_s: np.ndarray
    vz_m_s: np.ndarray
    vturb_m_s: np.ndarray
    logg_cgs: float
    logg_si: float


def column_inputs(snapshot: Mapping, *, star_id: str, ix: int, iy: int,
                  vturb_m_s: np.ndarray) -> ColumnInputs:
    """Select a column and resolve gravity; require the microturbulence explicitly.

    No default xi is borrowed from a 1D analysis: whether unresolved broadening
    is appropriate for a resolved-velocity cube is a benchmark physics choice.
    Preserve all three velocity components so a later oblique-ray calculation
    cannot mistake a vertical projection for complete 3D transport.
    """
    nz, ny, nx = validate_snapshot(snapshot)
    for name, index, limit in (("ix", ix, nx), ("iy", iy, ny)):
        if isinstance(index, (bool, np.bool_)) or not isinstance(index, (int, np.integer)) or not 0 <= index < limit:
            raise FoundationInputError(f"{name}: column index outside snapshot")
    vturb = _finite_array(vturb_m_s, "vturb_m_s", (nz,))
    if np.any(vturb < 0):
        raise FoundationInputError("vturb_m_s: must be nonnegative")
    from config.constants import get_star_params

    logg = float(get_star_params(star_id)["logg"])
    if not np.isfinite(logg):
        raise FoundationInputError("canonical logg: must be finite")
    # g_SI = g_cgs / 100; this is a unit conversion, not a solar parameter.
    logg_si = logg - np.log10(100.0)

    def take(name):
        return np.ascontiguousarray(np.asarray(snapshot[name])[:, iy, ix], dtype=np.float64).copy()

    return ColumnInputs(star_id, ix, iy, np.asarray(snapshot["z"]).copy(),
                        take("temperature"), take("electron_density"),
                        take("hydrogen_density") + take("proton_density"),
                        take("vx"), take("vy"), take("vz"), vturb.copy(), logg, logg_si)


def make_vertical_atmosphere(column: ColumnInputs, *, lightweaver, lower_bc, upper_bc):
    """Construct an atmosphere only, never Context/SE/RT; explicit boundary objects.

    This adapter is for the inspected 0.16.1 camelCase API. A different signature
    fails loudly rather than substituting defaults. This is a vertical 1D column;
    the retained horizontal velocities are deliberately not projected into it.
    Constructing it does not authorize a production solve or assert H/EOS parity.
    """
    if lower_bc is None or upper_bc is None:
        raise FoundationInputError("both boundary conditions must be supplied explicitly")
    return lightweaver.Atmosphere.make_1d(
        scale=lightweaver.ScaleType.Geometric,
        depthScale=column.height_m.copy(), temperature=column.temperature_K.copy(),
        ne=column.ne_m3.copy(), nHTot=column.nH_total_m3.copy(),
        vlos=column.vz_m_s.copy(), vturb=column.vturb_m_s.copy(),
        logG=column.logg_si, lowerBc=lower_bc, upperBc=upper_bc,
    )
