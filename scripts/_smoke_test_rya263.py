"""Smoke test for RYA-263 fix: confirm per-line a_1dlte is non-zero for Procyon."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import pandas as pd
from pipeline.abundances_derive import run

# Procyon params (hardcoded — STAR_PROCYON not yet on this branch)
PROCYON = {'teff_K': 6530, 'logg': 3.96, 'feh': -0.04, 'vturb_kms': 1.66}

if __name__ == '__main__':
    print("RYA-263 smoke test: Procyon Fe I per-line abundances (post fix)\n")
    converged_params, results = run('procyon', stellar_params_override=PROCYON)

    # run() returns (converged_params dict, results DataFrame)
    print("\n--- Gate table (all elements) ---")
    print(results.to_string(index=False))

    fe1_row = results[(results['element'] == 'Fe') & (results['ion'] == 'I')]

    if not fe1_row.empty:
        a_fe = float(fe1_row['A_X'].iloc[0])
        std  = float(fe1_row['A_X_std'].iloc[0])
        n    = int(fe1_row['n_lines'].iloc[0])
        print(f"\nFe I: A(Fe)={a_fe:.4f}  std={std:.4f}  n={n}")
        if a_fe != 0.0:
            print("SMOKE TEST: PASS")
        else:
            print("SMOKE TEST: FAIL — A(Fe I) still 0.0")
            sys.exit(1)
    else:
        print("SMOKE TEST: FAIL — no Fe I row in results")
        sys.exit(1)
