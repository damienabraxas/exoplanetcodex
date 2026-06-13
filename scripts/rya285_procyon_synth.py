"""
scripts/rya285_procyon_synth.py
RYA-285 Step 2 — Procyon Turbospectrum synthesis-EW run at the Gaia FGK
benchmark anchor (Heiter et al. 2015: Teff 6550, log g 3.99, [Fe/H] +0.01),
at both ξ = 1.66 km/s (current fixed value) and ξ = 2.00 km/s (benchmark
literature value / RYA-284 target).

Params are PINNED (skip_convergence=True) so the engine diff and the ξ effect
stay separable — convergence is not allowed to walk Teff/log g.

Writes:
  procyon_xi166_abundances_synth.csv  (+ _per_line_synth.csv)
  procyon_xi200_abundances_synth.csv  (+ _per_line_synth.csv)
  procyon_abundances_synth.csv        (= the ξ=1.66 anchor, for rya_engine_diff)

Run from repo root: python scripts/rya285_procyon_synth.py
"""
import sys
import time
import shutil
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from config.constants import PATHS
from pipeline.abundances_derive import run

PROCYON_BENCHMARK = {'teff_K': 6550.0, 'logg': 3.99, 'feh': 0.01}
XI_VALUES = [1.66, 2.00]
OUT_DIR = Path(str(PATHS['solar_ew'])).parent


def _tag(xi: float) -> str:
    return f"xi{int(round(xi * 100)):03d}"


def main():
    for xi in XI_VALUES:
        tag = _tag(xi)
        print("\n" + "#" * 68)
        print(f"#  Procyon synthesis  |  benchmark anchor  |  ξ = {xi:.2f} km/s  ({tag})")
        print("#" * 68)

        params = dict(PROCYON_BENCHMARK, vturb_kms=float(xi))
        t0 = time.time()
        run('procyon', model_grid='ATLAS9.Castelli',
            stellar_params_override=params,
            engine='synthesis', skip_convergence=True)
        wall = time.time() - t0
        print(f"\n  [{tag}] wall-clock: {wall:.1f} s ({wall/60:.1f} min)")

        # Tag the outputs so both ξ runs are preserved
        for base in ('abundances_synth', 'per_line_synth'):
            src = OUT_DIR / f'procyon_{base}.csv'
            dst = OUT_DIR / f'procyon_{tag}_{base}.csv'
            if src.exists():
                shutil.copy(src, dst)
                print(f"  [{tag}] saved → {dst.name}")

    # Leave the ξ=1.66 anchor as the canonical procyon_abundances_synth.csv
    anchor = OUT_DIR / f'procyon_{_tag(1.66)}_abundances_synth.csv'
    canon  = OUT_DIR / 'procyon_abundances_synth.csv'
    if anchor.exists():
        shutil.copy(anchor, canon)
        print(f"\n  Canonical (ξ=1.66) → {canon.name} for rya_engine_diff")


if __name__ == '__main__':
    main()
