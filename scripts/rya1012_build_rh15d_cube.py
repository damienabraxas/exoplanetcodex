"""
scripts/rya1012_build_rh15d_cube.py
===================================
RYA-1012 — build an RH 1.5D HDF5 model atmosphere that VoronoiRT can read.

WHY THIS EXISTS. VoronoiRT's `get_atmos()` (src/atmosphere.jl) reads the **RH 1.5D
HDF5 atmosphere format** — Pereira's own format, and he co-authored VoronoiRT. The
match is exact on every axis: dataset names (`x`, `y`, `z`, `temperature`,
`electron_density`, `hydrogen_populations`, `velocity_x/y/z`), SI units
(m, m/s, K, m^-3), and the otherwise-cryptic `[:,:,:,1,1]` indexing on
`hydrogen_populations`, which lands on RH's `(nt, nhydr, nx, ny, nz)` layout
(level 1, snapshot 1). Muspel.jl's `read_atmos_rh` reads the identical layout.

So a cube built once in this format is readable by VoronoiRT, Muspel.jl, RH 1.5D
and helita alike. THIS IS THE POINT: we are not inventing a format, we are writing
the ecosystem's lingua franca.

TWO NON-OBVIOUS EDITS ARE REQUIRED (both cost a run to find):

  (1) `velocity_x` and `velocity_y` MUST be present. `get_atmos` reads them
      unconditionally, but the RH 1.5D format itself does not require them —
      RH is a 1.5D (column-by-column) code and only needs the line-of-sight
      `velocity_z`. Muspel's own bundled FALC.hdf5 omits them.

  (2) `x` and `y` MUST have real, non-zero spacing. Muspel's FALC.hdf5 carries
      `x = [0,0,0]`, `y = [0,0,0]` — a 1D model tiled with ZERO horizontal
      extent. That is correct for 1.5D RT and FATAL for a 3D solver, which
      divides by dx/dy. The failure is silent and nasty: VoronoiRT exits 0 and
      writes a well-formed, entirely-NaN intensity map. Populations, opacity and
      the source function are all clean; the NaN is born inside
      `short_characteristics_up`. ALWAYS read the values, never the exit code.

NOTE ON `dxy`: a horizontally-tiled 1D atmosphere has no intrinsic horizontal
scale, so `dxy` is a FREE PARAMETER set by the caller. It does not affect the
1D-consistency check (a homogeneous atmosphere gives the same emergent intensity
whatever the spacing), but it DOES set the Voronoi cell geometry and therefore
the irregular-vs-regular comparison. Record whatever you choose.

Usage
-----
    # Tile a 1D RH atmosphere (e.g. Muspel's FALC) into a VoronoiRT-readable cube
    python scripts/rya1012_build_rh15d_cube.py from-rh \
        --src /path/to/FALC.hdf5 --out cube.hdf5 --dxy 48000

Verification (what a correct cube produces, from RYA-1012):
    regular-grid continuum, 500 nm, disk centre = 42.29 kW nm^-1 m^-2
    (solar reference ~40.8 -> within 3.6%), uniform across the map to 6e-14.
"""

import argparse
import shutil

import h5py
import numpy as np

# Datasets VoronoiRT's get_atmos() reads. Anything missing here is a hard failure.
REQUIRED = (
    "x",
    "y",
    "z",
    "temperature",
    "electron_density",
    "hydrogen_populations",
    "velocity_x",
    "velocity_y",
    "velocity_z",
)


def _horizontal_axes(f, dxy):
    """Replace degenerate x/y with real spacing -- see edit (2) in the docstring."""
    for name in ("x", "y"):
        n = len(f[name][...])
        del f[name]
        f.create_dataset(name, data=(np.arange(n) * dxy).astype("f4"))


def from_rh(src, out, dxy):
    """Make an existing RH 1.5D atmosphere readable by VoronoiRT.

    The source is copied first: we never modify the upstream file.
    """
    shutil.copy(src, out)
    with h5py.File(out, "a") as f:
        vz = f["velocity_z"][...]
        for name in ("velocity_x", "velocity_y"):
            if name not in f:
                # Edit (1): get_atmos reads these unconditionally.
                f.create_dataset(name, data=np.zeros_like(vz))
        _horizontal_axes(f, dxy)
        verify(f)
    return out


def verify(f):
    """Fail loudly on anything get_atmos would trip over."""
    missing = [d for d in REQUIRED if d not in f]
    if missing:
        raise SystemExit(f"FAIL: missing datasets {missing}")

    for name in ("x", "y"):
        axis = np.asarray(f[name][...], dtype=float)
        if axis.size > 1 and np.ptp(axis) == 0.0:
            raise SystemExit(
                f"FAIL: '{name}' has zero extent -- a 3D solver divides by d{name} "
                "and will emit an all-NaN map while exiting 0. Set --dxy."
            )

    hpop = f["hydrogen_populations"]
    if hpop.ndim != 5:
        raise SystemExit(
            f"FAIL: hydrogen_populations is {hpop.ndim}-D, expected 5-D "
            "(nt, nhydr, nx, ny, nz) -- get_atmos indexes [:,:,:,1,1]."
        )

    for name in ("temperature", "electron_density", "hydrogen_populations"):
        arr = np.asarray(f[name][...], dtype=float)
        if not np.all(np.isfinite(arr)):
            raise SystemExit(f"FAIL: '{name}' contains non-finite values")
        if arr.min() <= 0:
            raise SystemExit(f"FAIL: '{name}' has non-positive values")

    temp = np.asarray(f["temperature"][...], dtype=float)
    print("OK  datasets  :", sorted(f.keys()))
    print("OK  shape     :", f["temperature"].shape, "(nt, nx, ny, nz)")
    print("OK  hpops     :", hpop.shape, "(nt, nhydr, nx, ny, nz)")
    print(f"OK  T   [K]   : {temp.min():.1f} .. {temp.max():.1f}")
    print("OK  x   [m]   :", f["x"][...])
    print("OK  y   [m]   :", f["y"][...])


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("from-rh", help="tile/repair an existing RH 1.5D atmosphere")
    p.add_argument("--src", required=True, help="source RH 1.5D HDF5 (not modified)")
    p.add_argument("--out", required=True, help="destination cube")
    p.add_argument(
        "--dxy",
        type=float,
        required=True,
        help="horizontal grid spacing in metres (FREE PARAMETER -- record it)",
    )

    v = sub.add_parser("verify", help="check an existing cube against get_atmos")
    v.add_argument("--cube", required=True)

    args = ap.parse_args()
    if args.cmd == "from-rh":
        out = from_rh(args.src, args.out, args.dxy)
        print(f"\nwrote {out}  (dxy = {args.dxy:g} m)")
    else:
        with h5py.File(args.cube, "r") as f:
            verify(f)


if __name__ == "__main__":
    main()
