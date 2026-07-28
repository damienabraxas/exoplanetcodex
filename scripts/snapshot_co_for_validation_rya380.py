#!/usr/bin/env python3
"""
scripts/snapshot_co_for_validation_rya380.py
============================================
RYA-380 — snapshot RYA-373's per-frame conditioned CRIRES+ CO products (gitignored
*.fits) to the committed CSV the RYA-390 validation harness consumes, NOW CARRYING the
persisted molecfit transmission model (`mtrans`). This is what unblocks RYA-390
telluric-MODEL check 2 (mtrans vs Wallace) — previously BLOCKED because only the
corrected flux was saved.

Columns: wave_vac_A (VACUUM Å; CRIRES IDP / molecfit WAVELENGTH_FRAME=VAC),
flux_norm (continuum-normalized), err, mtrans (molecfit BEST_FIT_MODEL.mtrans). The
FITS provenance cards are copied into the `# KEY = value` comment block the harness
parses (RESTFRM/SPECSYS/RVSTATUS/GAP1/…/GDAS/MTRANS).

Usage:  python -m scripts.snapshot_co_for_validation_rya380
"""
from __future__ import annotations

import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits

from config.constants import PATHS

COND = Path(PATHS['data_root']) / 'audit' / 'crires_co_conditioned'
# header cards copied into the CSV provenance block (in this order)
_PROV_KEYS = ('RYA', 'PROVIS', 'TELL', 'CONTNORM', 'RESTFRM', 'SPECSYS', 'WLEN',
              'GDAS', 'GATE', 'SNRFRAME', 'RVSTATUS', 'MTRANS', 'GAP1', 'GAP2')


def _git_desc() -> str:
    try:
        return subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD'],
                                       cwd=str(COND), text=True).strip()
    except Exception:
        return 'unknown'


def snapshot(fits_path: Path) -> Path:
    hdul = fits.open(fits_path)
    h = hdul[0].header
    d = hdul['CO_ORDER'].data
    cols = {'wave_vac_A': np.asarray(d['wave_A'], float),
            'flux_norm': np.asarray(d['flux_norm'], float),
            'err': np.asarray(d['err'], float)}
    if 'mtrans' in [c.name for c in d.columns]:
        cols['mtrans'] = np.asarray(d['mtrans'], float)
    df = pd.DataFrame(cols)

    out = fits_path.with_suffix('.csv')
    head = [f"# RYA-380 conditioned CRIRES+ Vesta K-band CO (PROVISIONAL) — snapshot @ "
            f"build/rya-380 {_git_desc()}; input to RYA-390 Part B validation.",
            "# wave_vac_A: VACUUM Å (CRIRES IDP / molecfit WAVELENGTH_FRAME=VAC). "
            "flux_norm: continuum-normalized. err. mtrans: molecfit transmission model."]
    for k in _PROV_KEYS:
        if k in h:
            head.append(f"# {k} = {h[k]}")
    with open(out, 'w') as fh:
        fh.write("\n".join(head) + "\n")
        df.to_csv(fh, index=False)
    has_mt = 'mtrans' in df.columns and np.isfinite(df['mtrans']).any()
    print(f"  {fits_path.name} -> {out.name}  ({len(df)} rows; "
          f"mtrans {'persisted' if has_mt else 'ABSENT'})")
    return out


def main() -> int:
    fits_products = sorted(COND.glob('vesta_crires_K_CO_*_topocent_PROVISIONAL.fits'))
    if not fits_products:
        raise SystemExit(f"no conditioned FITS products in {COND} — run "
                         f"`python -m pipeline.crires_telluric --condition` first.")
    print(f"RYA-380 snapshot — {len(fits_products)} conditioned CO product(s):")
    for p in fits_products:
        snapshot(p)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
