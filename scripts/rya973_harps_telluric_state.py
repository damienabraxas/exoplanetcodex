#!/usr/bin/env python3
"""RYA-973 — DETERMINE the HARPS telluric state for a star, and gate. Do not heavy-lift.

The optical is telluric-LIGHT, so the answer here is a disposition, not a molecfit run
(RYA-460): inside HARPS's range the only registered telluric band is **O2 B, 6867-6884 A**.
What this script produces is the state, measured two ways, and the red-edge interval a
consumer must then exclude or correct.

TWO ROUTES, because they fail differently:

* **HEADERS** (`pipeline.telluric_intake`, RYA-806) — authoritative when the product says
  something. HARPS Phase-3 `SCIENCE.SPECTRUM` says NOTHING: no PRO REC chain, no
  transmission extension, no telluric column. That is `unknown`, and `telluric_policy`
  REFUSES unknown rather than defaulting it either way.
* **THE FLUX** — a header can be silent; the O2 band cannot. O2 is purely terrestrial with
  no photospheric counterpart, so the fraction of band pixels driven below 0.7 of the
  local continuum, measured against a telluric-FREE control shoulder on the same spectrum,
  says what the header would not. RYA-931 calibrated the scale on solar HARPS: 22.06%
  before correction, 0.00% after.

    python3 scripts/rya973_harps_telluric_state.py --star tau_ceti
"""
from __future__ import annotations

import argparse
import glob
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from config.constants import codex_path                       # noqa: E402
from pipeline.telluric_intake import from_many                # noqa: E402

#: The only registered telluric band inside HARPS's range (RYA-931/460).
O2_B = (6867.0, 6884.0)
#: A telluric-free shoulder on the SAME spectrum. Without a control the O2 number is
#: uncalibrated: it measures the star's own line density as much as the atmosphere's.
CONTROL = (6800.0, 6840.0)
DEPTH = 0.7


def _frac_below(w, f, lo, hi, thr=DEPTH):
    m = (w >= lo) & (w <= hi) & np.isfinite(f)
    if m.sum() < 50:
        return float('nan'), 0
    cont = np.nanpercentile(f[m], 95)
    if not np.isfinite(cont) or cont <= 0:
        return float('nan'), int(m.sum())
    return float(np.mean(f[m] / cont < thr)), int(m.sum())


def measure(files) -> dict:
    from astropy.io import fits
    rows = []
    for path in files:
        try:
            with fits.open(path) as h:
                hdr = h[0].header
                d = h[1].data[0]
                w = np.asarray(d['WAVE'], float)
                f = np.asarray(d['FLUX'], float)
        except Exception as exc:
            rows.append({'file': Path(path).name, 'error': f"{type(exc).__name__}: {exc}"})
            continue
        if not np.isfinite(w).any() or np.nanmax(w) < O2_B[1]:
            rows.append({'file': Path(path).name, 'error': 'does not reach the O2 B band'})
            continue
        o2, n_o2 = _frac_below(w, f, *O2_B)
        ctl, _ = _frac_below(w, f, *CONTROL)
        rows.append({'file': Path(path).name, 'date_obs': str(hdr.get('DATE-OBS', ''))[:10],
                     'o2b_frac_below_0.7': o2, 'control_frac_below_0.7': ctl,
                     'n_px_o2b': n_o2,
                     # A control that is itself absorbed means the frame is unusable for
                     # this test -- reported, never quietly averaged in.
                     'control_usable': bool(np.isfinite(ctl) and ctl < 0.05)})
    good = [r for r in rows if r.get('control_usable')]
    return {'rows': rows, 'n': len(rows), 'n_usable': len(good),
            'o2b_median': float(np.median([r['o2b_frac_below_0.7'] for r in good])) if good else float('nan'),
            'control_median': float(np.median([r['control_frac_below_0.7'] for r in good])) if good else float('nan')}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument('--star', default='tau_ceti')
    ap.add_argument('--dir', default=None)
    ap.add_argument('--limit', type=int, default=None)
    ap.add_argument('--out', default=None)
    a = ap.parse_args()

    d = Path(a.dir) if a.dir else codex_path('data.spectra') / 'tau_cet' / 'HARPS'
    files = sorted(set(glob.glob(str(d / '**' / '*.fits'), recursive=True)))[:a.limit]
    if not files:
        raise SystemExit(f"no HARPS products under {d}")

    header_state, evidence = from_many(files)
    m = measure(files)
    ratio = (m['o2b_median'] / m['control_median']) if m['control_median'] else float('inf')
    # The verdict is the FLUX where the header is silent -- and the header being silent is
    # itself recorded, not smoothed over.
    flux_state = ('not-applied' if np.isfinite(ratio) and ratio > 10
                  else 'applied' if np.isfinite(ratio) and ratio < 3
                  else 'indeterminate')
    out = {
        'ticket': 'RYA-973', 'star': a.star, 'dir': str(d), 'n_products': len(files),
        'header_state': header_state,
        'header_reason': evidence[0].reasons[0] if evidence and evidence[0].reasons else '',
        'o2b_window_A': list(O2_B), 'control_window_A': list(CONTROL),
        'depth_threshold': DEPTH,
        'o2b_median_frac_below': m['o2b_median'],
        'control_median_frac_below': m['control_median'],
        'excess_ratio': ratio,
        'flux_state': flux_state,
        'n_frames_measured': m['n'], 'n_frames_control_usable': m['n_usable'],
        'calibration': ('RYA-931 solar HARPS: 22.06% of O2 B pixels below 0.7 before '
                        'correction, 0.00% after'),
        'disposition': ('DETERMINE-AND-GATE (RYA-460): the optical is telluric-light, so '
                        'the owed action is a red-edge disposition, NOT a molecfit run. '
                        'Any line inside 6867-6884 A needs exclusion or correction; lines '
                        'blueward of it are CLEAN_OUTSIDE_REGISTERED_BANDS.'),
        'per_frame': m['rows'],
    }
    dest = Path(a.out) if a.out else (ROOT / 'data' / 'audit' / 'rya973_harps_telluric'
                                      / f'{a.star}_harps_telluric_state.json')
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(json.dumps(out, indent=2) + '\n')
    print(f"{a.star} HARPS: {len(files)} products")
    print(f"  header state : {header_state}")
    print(f"  O2 B median  : {m['o2b_median']*100:.2f}%  vs control {m['control_median']*100:.2f}%"
          f"  (x{ratio:.0f})")
    print(f"  FLUX state   : {flux_state}")
    print(f"  usable frames: {m['n_usable']}/{m['n']}")
    print(f"[out] {dest}")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
