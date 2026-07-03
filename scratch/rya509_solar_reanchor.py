"""RYA-509: solar clean-from-raw re-validation through the RYA-506-fixed pipeline.
Compares reproduced solar A(X) to the BANKED v1 reference (not a memorized value)."""
import os
for v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[v] = '1'
import sys
sys.path.insert(0, '/Users/ryanschmitt/Documents/Exoplanet Codex/exoplanetcodex-rya509')
sys.path.insert(0, '/Users/ryanschmitt/Documents/Exoplanet Codex/ispec')
import warnings; warnings.filterwarnings('ignore')
import numpy as np
import pandas as pd
from pipeline import abundances_derive as ad
from pipeline import data_namespace as ns
import multiprocessing as mp
print("start_method =", mp.get_start_method())

cp, res = ad.run('solar', engine='spectrum')
ref, ver = ns.read_solar_reference('CURRENT')

print("\n=== CONVERGED:", {k: (round(v, 3) if isinstance(v, float) else v)
                          for k, v in cp.items() if k in ('teff_K', 'logg', 'feh', 'vturb_kms')})

# join reproduced vs banked on (element, ion)
key = ['element', 'ion']
a = res[['element', 'ion', 'A_X', 'A_X_nlte', 'n_lines']].copy()
b = ref[['element', 'ion', 'A_X', 'A_X_nlte', 'n_lines']].copy()
m = a.merge(b, on=key, how='outer', suffixes=('_repro', '_bank'))
m['dA_LTE'] = (m['A_X_repro'] - m['A_X_bank']).round(3)
m['dA_NLTE'] = (m['A_X_nlte_repro'] - m['A_X_nlte_bank']).round(3)
cols = ['element', 'ion', 'A_X_bank', 'A_X_repro', 'dA_LTE',
        'A_X_nlte_bank', 'A_X_nlte_repro', 'dA_NLTE', 'n_lines_bank', 'n_lines_repro']
print(f"\n=== REPRODUCED vs BANKED (ref {ver}) ===")
print(m[cols].to_string(index=False))

# headline: Fe + max abs delta
fe1 = m[(m.element == 'Fe') & (m.ion == 'I')]
if len(fe1):
    r = fe1.iloc[0]
    print(f"\nFe I: banked A_nlte={r['A_X_nlte_bank']} repro={r['A_X_nlte_repro']} "
          f"Δ={r['dA_NLTE']}  (banked n={r['n_lines_bank']} repro n={r['n_lines_repro']})")
maxd = m['dA_NLTE'].abs().max()
print(f"max |ΔA_NLTE| across species = {maxd:.3f} dex")
print("REPRODUCES within 0.05 dex on all species" if maxd is not None and maxd <= 0.05
      else "SOME species exceed 0.05 dex — inspect")
m[cols].to_csv('/private/tmp/claude-501/-Users-ryanschmitt/23af2ec8-2a58-48ac-a59f-e836f81767d5/scratchpad/rya509_compare.csv', index=False)
