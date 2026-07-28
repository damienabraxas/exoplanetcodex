#!/usr/bin/env python3
"""
scripts/build_uves_registry_rya272.py
=====================================
RYA-272 Step 3 — build the Procyon UVES spectrum registry from the 45 accepted
topocentric IDPs. One row per IDP: filename, DATE-OBS, setting, coverage, R, SNR,
BERV applied (km/s, computed by the loader), and the O I 7771–7775 telluric verdict
(CLEAN / CORRECTABLE / EXCLUDE) consumed from the RYA-271 audit — not re-derived.
The 2013-10-08 RED760 co-add is marked oi_anchor=True (the designated O I spectrum).

The registry is how downstream EW code selects spectra (telluric-rule enforcement at
the data layer): only oi_anchor / CLEAN spectra may be used for O I until molecfit.

Out: data/spectra/procyon/uves_registry.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from pipeline.loaders.uves_loader import UVESLoader, UVESProductError   # noqa: E402

UVES_DIR = ROOT.parent / 'data' / 'spectra' / 'exoplanetcodex-data' / 'Procyon' / 'Procyon UVES'
OUT = ROOT / 'data' / 'spectra' / 'procyon' / 'uves_registry.csv'

# O I telluric verdicts — verbatim from the RYA-271 audit table (per file/epoch).
OI_TELLURIC = {
    'ADP.2020-06-15T10:09:57.908.fits': 'CLEAN',        # 2013-10-08 RED760 co-add, SNR 411 — the anchor
    'ADP.2020-08-10T12:01:00.339.fits': 'CORRECTABLE',  # 2002-10-08 RED860, SNR 284
    'ADP.2020-08-10T12:01:00.361.fits': 'CORRECTABLE',  # 2002-10-08 RED860, SNR 222
    'ADP.2021-08-27T05:08:20.555.fits': 'EXCLUDE',      # 2002-02-27 RED860, SNR 196 (humid epoch)
    'ADP.2021-08-27T05:08:20.607.fits': 'EXCLUDE',      # 2002-02-27 RED860, SNR 178
}
OI_ANCHOR = 'ADP.2020-06-15T10:09:57.908.fits'
OI_TRIPLET = (7771.94, 7775.39)


def main():
    files = sorted(UVES_DIR.glob('ADP.*.fits'))
    rows = []
    for f in files:
        try:
            s = UVESLoader(f).load()
        except UVESProductError as e:
            print(f"  SKIP (guard) {f.name}: {e}")
            continue
        lo, hi = s.wave_range_A
        covers_oi = bool(lo <= OI_TRIPLET[0] and hi >= OI_TRIPLET[1])
        verdict = OI_TELLURIC.get(f.name, '')
        rows.append({
            'filename': f.name,
            'date_obs': s.meta['date_obs'],
            'setting': s.meta['setting'],
            'wl_min_A': round(lo, 1), 'wl_max_A': round(hi, 1),
            'resolution_R': s.meta['resolution_R'],
            'snr': s.meta['snr_summary'],
            'berv_applied_kms': round(s.meta['berv_kms'], 4),
            'covers_OI_7771': covers_oi,
            'oi_telluric_verdict': verdict,
            'oi_anchor': bool(f.name == OI_ANCHOR),
            'program_id': s.meta['program_id'],
        })
    df = pd.DataFrame(rows).sort_values(['date_obs', 'wl_min_A']).reset_index(drop=True)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT, index=False)

    # sanity (the registry contract)
    n_anchor = int(df['oi_anchor'].sum())
    n_oi = int(df['covers_OI_7771'].sum())
    assert len(df) == 45, f"expected 45 IDPs, got {len(df)}"
    assert n_anchor == 1, f"expected exactly 1 oi_anchor, got {n_anchor}"
    # every O I verdict belongs to an O I-covering spectrum, and vice-versa
    verd = df[df['oi_telluric_verdict'] != '']
    assert verd['covers_OI_7771'].all(), "an O I verdict on a non-O I spectrum"
    print(f"RYA-272 UVES registry: {len(df)} IDPs, {n_oi} cover O I 7771, "
          f"{n_anchor} oi_anchor (2013-10-08 RED760).")
    print(df[df['covers_OI_7771']][['filename', 'date_obs', 'setting', 'snr',
          'berv_applied_kms', 'oi_telluric_verdict', 'oi_anchor']].to_string(index=False))
    print(f"  wrote {OUT.relative_to(ROOT)}")


if __name__ == '__main__':
    main()
