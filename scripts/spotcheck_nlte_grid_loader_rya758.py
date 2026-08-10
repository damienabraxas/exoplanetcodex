#!/usr/bin/env python3
"""
scripts/spotcheck_nlte_grid_loader_rya758.py
============================================
RYA-758 smoke-test step 8 — prove the audit's grid-tier claims sit on a REAL
integration path, not on a paper URL.

The ticket asks for "the ``pipeline.nlte_grid_loader`` interface".  No such module
exists; the live entry point is ``pipeline.nlte_corrections._load_mpia_element_grid``
(per-element, registry-driven, RYA-413-guarded), reached in production through
``apply_element_nlte_corrections`` / ``element_grid_in_bounds``.  That is what this
script exercises, on already-registered CONTROL elements only.

Controls are deliberately elements we already carry (Ca = Mashonkina+2017 via MPIA
Spectrum Tools, Ti = Mallinson+2024 via the Amarsi PySME family).  No new grid is
downloaded: grids install on Sirius only, never the Mac.

    python scripts/spotcheck_nlte_grid_loader_rya758.py
"""
from __future__ import annotations

import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

# (element, a wavelength that IS in the grid, solar teff/logg/feh)
CONTROLS = [
    ('Ca', 6122.220),
    ('Ti', 5689.460),
]


def main() -> int:
    import pandas as pd
    from config.constants import NLTE_CORRECTION_ELEMENTS, committed_grid_artifact
    from pipeline import nlte_corrections as nc

    teff, logg, feh = 5772.0, 4.438, 0.0            # solar, per config/stars.yaml
    failures = []

    print('entry point: pipeline.nlte_corrections._load_mpia_element_grid '
          '(there is no pipeline.nlte_grid_loader)')
    print(f'grid dir   : {committed_grid_artifact("nlte_grids")}')
    print()

    for element, wave in CONTROLS:
        spec = NLTE_CORRECTION_ELEMENTS[element]
        path = committed_grid_artifact('nlte_grids', spec['grid'])
        try:
            raw = pd.read_csv(path)
            grid = nc._load_mpia_element_grid(element)
            delta = nc._mpia_element_delta(element, wave, teff, logg, feh)
            in_box = nc.element_grid_in_bounds(element, teff, logg, feh)
        except Exception as exc:                    # LOUD: never swallow
            failures.append(f'{element}: {type(exc).__name__}: {exc}')
            print(f'[FAIL] {element}: {type(exc).__name__}: {exc}')
            continue

        print(f'[OK] {element}  file={spec["grid"]}')
        print(f'     raw table shape      : {raw.shape}  columns={list(raw.columns)}')
        print(f'     interpolated lines   : {len(grid["waves"])}  {list(grid["waves"])}')
        print(f'     ion / reference      : {grid["ion"]} / {grid["ref"][:70]}...')
        print(f'     (teff,logg,feh) hull : {grid["bounds"]}')
        print(f'     solar inside hull    : {in_box}')
        print(f'     delta_nlte @ {wave:9.3f} A, solar = {delta:+.4f} dex')
        print()

        if not in_box:
            failures.append(f'{element}: solar node outside the grid hull')
        if delta != delta:                          # NaN
            failures.append(f'{element}: delta_nlte is NaN at {wave} A')

    if failures:
        print('GRID LOADER SPOT-CHECK: FAIL')
        for f in failures:
            print('  -', f)
        return 1
    print('GRID LOADER SPOT-CHECK: PASS '
          f'({len(CONTROLS)} registered grids loaded and interpolated)')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
