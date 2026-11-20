"""Bride (default physics) vs our MULTI 2.3 run with the model's own ABUND. argv: record-glob multi_run_dir [variant k=v,...]"""
import sys, os, json, glob
sys.path.insert(0, '.'); sys.path.insert(0, 'scripts')
import numpy as np
from dataclasses import replace
from pathlib import Path
from config.constants import require_codex_path
from pipeline.bride_reference import read_marcs_structure
from pipeline.bride_columns import column_context, marcs_column, solve_column
from pipeline.multi_atom import read_multi_atom
import bride_qa_benchmark as Q
rec_path = sorted(glob.glob(os.path.expanduser('~/codex/bride_inputs/rya1282/bench_gerber2023_O_MARCS_zladder/' + sys.argv[1])))[0]
rec = json.loads(Path(rec_path).read_text()); ref = rec['reference']
atom = replace(read_multi_atom(require_codex_path('experiments.bride_oxygen_atom')), abundance=float(ref['abundance']))
S = read_marcs_structure(require_codex_path('grids.marcs_standard') / rec['atmosphere_file']['name'])
bm = np.load(Path(sys.argv[2]) / 'b.npy')
kw = dict(abundances_dex=S['abundances_dex'], logg_cgs=S['logg_cgs'], vturb_m_s=S['vturb_m_s'], rays=3)
for item in (sys.argv[3].split(',') if len(sys.argv) > 3 else []):
    k, v = item.split('=')
    if k == 'extras':
        import pipeline.bride_ts_background as TB
        _f = TB.ts_metal_factory
        TB.ts_metal_factory = lambda *a, **q: _f(*a, **{**q, 'extras': bool(int(v))}); continue
    if k == 'haze':
        import pickle; kw['line_haze'] = dict(lines=pickle.load(open('/private/tmp/claude-501/-Users-ryanschmitt/8eefa0eb-2b5d-45e6-98b2-baa84a382ba2/scratchpad/vald_atoms.pkl', 'rb')), abundances_dex=S['abundances_dex']); continue
    kw[k] = type(kw[k])(v) if k in kw else v
ctx = column_context(marcs_column(S), atom, ng=True, **kw)
n, hist, _, _ = solve_column(ctx, maxiter=6000, seconds=3000.)
b = n / ctx.activeAtoms[0].nStar
lt = np.asarray(S['logtau_500']); phot = (lt >= -3) & (lt <= 0)
TL = list(Q.TRIPLET_LEVELS)
d = np.log10(b) - np.log10(bm)
dt = d[TL]
i = int(np.argmax(np.abs(np.where(phot, dt, 0)).max(0)))
print(f'{Path(rec_path).stem[:22]} {sys.argv[3] if len(sys.argv)>3 else "default"} Bride-vs-MULTI(model ABUND) triplet phot max|d| {np.abs(dt[:, phot]).max():.4f} at logtau {lt[i]:+.2f}; median {np.median(np.abs(dt[:, phot])):.4f}')
worst = np.argsort(-np.abs(d[:, phot]).max(1))[:6]
print('   worst levels (idx, max|dlogb| phot):', [(int(k), round(float(np.abs(d[k, phot]).max()), 3)) for k in worst])
for tau in (-3, -2, -1, 0):
    j = int(np.argmin(abs(lt - tau)))
    print(f'   logtau {lt[j]:+.2f}  ground b bride/multi {b[0,j]:.4f}/{bm[0,j]:.4f}  triplet low {b[TL[0],j]:.4f}/{bm[TL[0],j]:.4f}  up {b[TL[1],j]:.4f}/{bm[TL[1],j]:.4f}  OII {b[-1,j]:.3f}/{bm[-1,j]:.3f}')
