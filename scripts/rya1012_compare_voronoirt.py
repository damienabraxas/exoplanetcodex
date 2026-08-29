"""
scripts/rya1012_compare_voronoirt.py
====================================
RYA-1012 — score VoronoiRT's output against its own published reference arrays
(searchlight) and against its own regular-grid arm (continuum).

TWO DIFFERENT KINDS OF COMPARISON, and the distinction matters:

  searchlight -- an EXTERNAL reference. data/searchlight_data/*.npy are the
      author's published arrays. Only the REGULAR-grid ones are reproducible:
      searchlight_regular() has no RNG at all, while searchlight_irregular()
      draws rand(3, n_sites) with NO Random.seed! (seeds exist in the sibling
      scripts, so the omission looks accidental). Measured: two runs of the
      irregular arm on the SAME machine minutes apart differ by 30.6% / 15.1%
      in total flux -- as large as the gap to the published arrays. So the
      Voronoi arm is not a reference for anyone, including the author. Always
      re-run a generator and measure its own spread before calling a mismatch
      a portability failure.

  continuum -- an INTERNAL reference. LTE_regular() and LTE_compare() both emit
      emergent intensity at theta=180 deg, lambda=500 nm; the criterion is
      whether the irregular grid reproduces the regular one. That IS the paper's
      claim, and it needs no published data. This is why the missing Bifrost
      cube does not block a meaningful continuum test.

Usage
-----
    python scripts/rya1012_compare_voronoirt.py searchlight \
        --ref  /path/VoronoiRT/data/searchlight_data \
        --ours /path/VoronoiRT/data/rya1012_out_run1

    python scripts/rya1012_compare_voronoirt.py continuum \
        --lte /path/VoronoiRT/data/LTE --n-sites 10000

Reference results (RYA-1012, 2026-08-29):
    searchlight  I_160_45_regular  bitwise identical
                 I_20_15_regular   9 ULP, max rel 1.16e-15   (at phi=195, not 15)
    continuum    regular 42.29 kW nm^-1 m^-2 (solar ~40.8, within 3.6%)
                 irregular vs regular: max rel diff 0.50%
"""

import argparse
import os

import numpy as np

# Solar disk-centre continuum intensity at 500 nm, ~4.08e13 W m^-3 sr^-1,
# expressed in the kW nm^-1 m^-2 that VoronoiRT writes.
SOLAR_I500_KW_NM_M2 = 4.08e13 * 1e-9 / 1e3


def _ulps(a, b):
    """Distance in units-in-the-last-place between two float64 arrays."""
    ia = np.frombuffer(np.ascontiguousarray(a, dtype=np.float64).tobytes(), dtype=np.int64)
    ib = np.frombuffer(np.ascontiguousarray(b, dtype=np.float64).tobytes(), dtype=np.int64)
    return int(np.abs(ia - ib).max())


def searchlight(ref_dir, ours_dir):
    names = [
        "x_regular.npy",
        "y_regular.npy",
        "I_20_15_regular.npy",
        "I_160_45_regular.npy",
        "I_20_15_voronoi.npy",
        "I_160_45_voronoi.npy",
    ]
    for name in names:
        rp, op = os.path.join(ref_dir, name), os.path.join(ours_dir, name)
        if not (os.path.exists(rp) and os.path.exists(op)):
            print(f"{name:26s} SKIP (absent)")
            continue
        ref, ours = np.load(rp), np.load(op)
        print("=" * 72)
        print(name)
        if ref.shape != ours.shape:
            print(f"  SHAPE MISMATCH ref={ref.shape} ours={ours.shape}")
            continue
        diff = np.abs(ref - ours)
        exact = int((ref == ours).sum())
        print(f"  shape={ref.shape}  ref.sum={ref.sum():.12g}  ours.sum={ours.sum():.12g}")
        print(f"  bitwise identical elements: {exact} / {ref.size}")
        print(f"  max|diff| = {diff.max():.6e}   max ULPs = {_ulps(ref, ours)}")
        nz = ref != 0
        if nz.any():
            print(f"  max rel diff = {(diff[nz] / np.abs(ref[nz])).max():.6e}")
        if "voronoi" in name:
            print("  NOTE: unseeded RNG -- NOT reproducible, not a valid reference.")


def continuum(lte_dir, n_sites):
    reg = os.path.join(lte_dir, "I_regular_full.npy")
    irr = os.path.join(lte_dir, f"I_irregular_{n_sites}_extinction.npy")
    for path in (reg, irr):
        if not os.path.exists(path):
            raise SystemExit(f"missing {path}")
    r, v = np.load(reg), np.load(irr)

    print("=" * 72)
    print("CONTINUUM -- internal reference (irregular vs regular), 500 nm, disk centre")
    for label, arr in (("regular", r), ("irregular", v)):
        if np.isnan(arr).any():
            # An all-NaN map is what a zero-extent x/y axis produces, at exit 0.
            raise SystemExit(
                f"FAIL: {label} map contains NaN -- check the cube's x/y spacing "
                "(see scripts/rya1012_build_rh15d_cube.py)"
            )
        print(f"  {label:10s} mean = {arr.mean():.4f} kW nm^-1 m^-2   shape={arr.shape}")

    print(f"  regular spread across map = {r.max() - r.min():.3e} "
          "(~0 expected iff the cube is horizontally homogeneous)")
    print(f"  solar reference at 500 nm = {SOLAR_I500_KW_NM_M2:.2f} kW nm^-1 m^-2")
    print(f"  regular / solar           = {r.mean() / SOLAR_I500_KW_NM_M2:.3f}")

    d = np.abs(v - r)
    print(f"  max|diff| irregular-regular = {d.max():.6e}")
    print(f"  max rel diff                = {100 * (d / np.abs(r)).max():.4f} %  <-- the paper's criterion")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("searchlight", help="score against the published arrays")
    s.add_argument("--ref", required=True, help="VoronoiRT data/searchlight_data")
    s.add_argument("--ours", required=True, help="our output dir")

    c = sub.add_parser("continuum", help="score irregular against regular")
    c.add_argument("--lte", required=True, help="VoronoiRT data/LTE")
    c.add_argument("--n-sites", type=int, default=10000)

    args = ap.parse_args()
    if args.cmd == "searchlight":
        searchlight(args.ref, args.ours)
    else:
        continuum(args.lte, args.n_sites)


if __name__ == "__main__":
    main()
