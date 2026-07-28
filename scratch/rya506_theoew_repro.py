#!/usr/bin/env python3
"""
scratch/rya506_theoew_repro.py — RYA-506 diagnostic (analysis-only, no production edits).

Reproduces + instruments the Procyon Fe theo-EW collapse: a clean-from-raw regen quarantines
100% of the Fe pool (theo < 5 mÅ, 106 → 0) → "MOOG: No abundances". This harness monkeypatches
`abundances_derive._fe2_theoretical_ew` to LOG the theoretical-EW array + the atmosphere it
consumes at the quarantine step (it does NOT change the computation — it wraps the original).

ROOT CAUSE (confirmed): iSpec `calculate_theoretical_ew_and_depth` runs the SPECTRUM synthesis
in a child multiprocessing.Process that imports the C `synthesizer` extension and expects FORK
semantics, and initialises `output_ew = np.zeros(num_lines)`. On macOS the Python 3.8+ default
start method is `spawn`, under which that child dies without enqueueing → iSpec silently
returns the zeros → every Fe line fails the theo<5 quarantine. FIX (applied in
pipeline/abundances_derive.py): force the start method to `fork` (eliminates it, deterministic)
+ a loud all-zero-batch guard. Post-fix the plain run reproduces A(FeI)=7.557/7.571. This
harness remains a probe of the theo-EW values; it is no longer needed to make the run succeed.

Prereq: stage Procyon data first (fast):
    python -m pipeline.spectra_normalize --star procyon
    python -m pipeline.lines_fit --star procyon
Then: python scratch/rya506_theoew_repro.py   (run from the repo root)
"""
import sys
from pathlib import Path

import numpy as np

_REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_REPO))
sys.path.insert(0, str(_REPO.parent / 'ispec'))
import warnings; warnings.filterwarnings('ignore')
from pipeline import abundances_derive as ad   # noqa: E402

_orig = ad._fe2_theoretical_ew
_calls = []


def _patched(fe_linemasks, stellar_params, atmosphere):
    theo = _orig(fe_linemasks, stellar_params, atmosphere)
    if len(_calls) < 3:
        t = np.asarray(theo, dtype=float)
        obs = np.array([float(x) for x in fe_linemasks['ew']])
        print(f"\n[THEO-EW CALL {len(_calls)}] n_lines={len(t)}  "
              f"params teff={stellar_params['teff_K']} logg={stellar_params['logg']} "
              f"feh={stellar_params['feh']} vturb={stellar_params['vturb_kms']}")
        print(f"  theo_ew mÅ: min={np.nanmin(t):.3f} median={np.nanmedian(t):.3f} "
              f"max={np.nanmax(t):.3f}  n(theo<5)={int((t < 5).sum())}/{len(t)}")
        print(f"  obs_ew  mÅ: min={np.nanmin(obs):.1f} median={np.nanmedian(obs):.1f} "
              f"max={np.nanmax(obs):.1f}")
        print("  --> theo all-~0 here == the iSpec dead-child silent-zeros bug; "
              "sane theo == child survived.")
        _calls.append(1)
    return theo


ad._fe2_theoretical_ew = _patched
try:
    cp, res = ad.run('procyon', engine='spectrum')
    print("\n=== RUN COMPLETED (theo-EW child survived) ===")
    print(res[['element', 'ion', 'A_X', 'A_X_nlte', 'n_lines']].to_string(index=False))
except Exception as e:
    print(f"\n=== RUN FAILED (theo-EW child died → silent zeros): {type(e).__name__}: {str(e)[:120]}")
