#!/usr/bin/env python3
"""RYA-963 — is `pipeline.acen_orbit`'s A/B branch assignment inverted?

The α Cen A CRIRES+ frames measure at the RV the module predicts for **B**. Two
explanations fit that: the holdings are mislabelled, or the module's A/B assignment is
swapped. Neither is decidable from the RV alone, so this script runs the two independent
tests that separate them and prints both.

**Test 1 — photometry, model-free.** α Cen A is the brighter component: K = -1.49 vs
-0.60, a flux ratio of 2.27. The A and B directories each hold a K2192 frame from the
same night, 16 minutes apart, through the same reduction, so their count RATES are
directly comparable. This identifies the components without any orbit at all.

**Test 2 — the ω convention.** Kervella Table 1's ω is documented in `acen_orbit` itself
as "component B relative to A" — the visual-orbit convention. A secondary's barycentric
orbit is parallel to the relative orbit, so ω_B = ω_rel; the primary's is ANTI-parallel,
so ω_A = ω_rel + 180°. The module assigns the opposite. This script evaluates the
predicted RVs both ways and compares each to the measurement.

    python3 scripts/rya963_acen_branch_check.py --measured-rv -19.17

Nothing here modifies `acen_orbit`. Flipping it would change the star-ID verdict of every
frame RYA-423 has ever judged, and the module's docstring says the present assignment was
pinned by RYA-384's spectral-type call — so it is a science decision to ratify, not a
one-line patch to apply in passing.
"""
from __future__ import annotations

import argparse
import glob
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from config.constants import codex_path                      # noqa: E402
from pipeline import acen_orbit as ao                        # noqa: E402
from pipeline.crires_telluric import load_crires_idp         # noqa: E402

#: 2MASS/Bessell K magnitudes. The only external numbers here, and they are photometry,
#: not a model: the brighter component IS alpha Cen A.
K_MAG_A, K_MAG_B = -1.49, -0.60
EXPECTED_RATIO = 10 ** (0.4 * (K_MAG_B - K_MAG_A))

VET = codex_path('data.spectra_local') / 'Alpha Centauri (vetted)'


def _k2192(subdir: str, epoch: str = '2022-04-15'):
    for f in sorted(glob.glob(str(VET / subdir / 'CRIRES' / '*.fits'))):
        fr = load_crires_idp(f)
        if fr.wlen_id == 'K2192' and str(fr.date_obs)[:10] == epoch:
            return fr
    raise SystemExit(f"no K2192 frame from {epoch} under {VET / subdir / 'CRIRES'}")


def brightness_test(epoch: str = '2022-04-15') -> dict:
    from astropy.io import fits
    a, b = _k2192('Alpha Cen A', epoch), _k2192('Alpha Cen B', epoch)
    out = {}
    for tag, fr in (('A', a), ('B', b)):
        h = fits.getheader(str(fr.path))
        rates = {}
        for s in fr.segments:
            m = np.isfinite(s.flux) & (s.flux > 0)
            if m.sum() > 500:
                rates[(s.order, s.detector)] = float(np.median(s.flux[m]))
        out[tag] = {'exptime': float(h['EXPTIME']), 'rates': rates, 'snr': fr.snr,
                    'seeing': next((float(h[k]) for k in h.keys()
                                    if 'AMBI FWHM START' in str(k)), float('nan')),
                    'n_saturated': next((int(h[k]) for k in h.keys()
                                         if 'QC NUMSAT' in str(k)), -1),
                    'dit': next((float(h[k]) for k in h.keys()
                                 if 'DET SEQ1 DIT' in str(k)), float('nan'))}
    common = sorted(set(out['A']['rates']) & set(out['B']['rates']))
    ratios = np.array([(out['A']['rates'][k] / out['A']['exptime'])
                       / (out['B']['rates'][k] / out['B']['exptime']) for k in common])
    return {'ratio_median': float(np.median(ratios)), 'ratio_min': float(ratios.min()),
            'ratio_max': float(ratios.max()), 'n_chips': len(common),
            'snr2_over_t': ((out['A']['snr'] ** 2 / out['A']['exptime'])
                            / (out['B']['snr'] ** 2 / out['B']['exptime'])),
            'expected_if_labels_correct': EXPECTED_RATIO,
            'expected_if_labels_swapped': 1.0 / EXPECTED_RATIO,
            'A': out['A'], 'B': out['B'], 'mjd': a.mjd}


def convention_test(mjd: float, measured_rv: float) -> dict:
    from astropy.time import Time
    by = Time(mjd, format='mjd').byear
    as_written = {'A': ao._rv(by, ao.K_A, ao._OMEGA_A),
                  'B': ao._rv(by, ao.K_B, ao._OMEGA_B)}
    # primary anti-parallel to the relative orbit
    visual = {'A': ao._rv(by, ao.K_A, ao._OMEGA_B),
              'B': ao._rv(by, ao.K_B, ao._OMEGA_A)}
    return {'byear': float(by), 'as_written': as_written, 'visual_convention': visual,
            'resid_as_written_A': abs(measured_rv - as_written['A']),
            'resid_visual_A': abs(measured_rv - visual['A'])}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument('--measured-rv', type=float, required=True,
                    help='barycentric RV measured for the A-directory frames, km/s')
    ap.add_argument('--epoch', default='2022-04-15')
    a = ap.parse_args()

    b = brightness_test(a.epoch)
    print("=" * 78)
    print("TEST 1 — photometry (model-free): which directory holds the BRIGHTER star?")
    print("=" * 78)
    print(f"  matched K2192 pair, {a.epoch}, {b['n_chips']} common chips")
    print(f"  A-dir seeing {b['A']['seeing']}\"  DIT {b['A']['dit']}s  "
          f"saturated px {b['A']['n_saturated']}")
    print(f"  B-dir seeing {b['B']['seeing']}\"  DIT {b['B']['dit']}s  "
          f"saturated px {b['B']['n_saturated']}")
    print(f"  count-rate ratio A/B = {b['ratio_median']:.3f} "
          f"(per-chip {b['ratio_min']:.3f}-{b['ratio_max']:.3f})")
    print(f"  SNR^2/exptime ratio  = {b['snr2_over_t']:.3f}")
    print(f"  expected, labels CORRECT : {b['expected_if_labels_correct']:.3f}")
    print(f"  expected, labels SWAPPED : {b['expected_if_labels_swapped']:.3f}")
    correct = (abs(b['ratio_median'] - b['expected_if_labels_correct'])
               < abs(b['ratio_median'] - b['expected_if_labels_swapped']))
    print(f"  --> the A directory holds the {'BRIGHTER' if correct else 'FAINTER'} star "
          f"=> labels are {'CORRECT' if correct else 'SWAPPED'}")
    print("  NOTE: the A frame has the most saturated pixels, and saturation SUPPRESSES")
    print("        its measured counts, so the true ratio is if anything LARGER.")

    c = convention_test(b['mjd'], a.measured_rv)
    print()
    print("=" * 78)
    print("TEST 2 — the omega convention in pipeline/acen_orbit.py")
    print("=" * 78)
    print(f"  measured RV of the A-directory frames: {a.measured_rv:+.2f} km/s")
    print(f"  as written        : rv_A={c['as_written']['A']:+.3f}  "
          f"rv_B={c['as_written']['B']:+.3f}   |resid_A| = {c['resid_as_written_A']:.2f}")
    print(f"  visual convention : rv_A={c['visual_convention']['A']:+.3f}  "
          f"rv_B={c['visual_convention']['B']:+.3f}   |resid_A| = {c['resid_visual_A']:.2f}")
    print()
    if correct and c['resid_visual_A'] < c['resid_as_written_A']:
        print("  VERDICT: the holdings are correctly labelled AND the measured RV matches")
        print("  alpha Cen A only under the visual-orbit convention. _OMEGA_A and _OMEGA_B")
        print("  in pipeline/acen_orbit.py are swapped.")
        print()
        print("  Scope: rv_bounds()/consistent_with_orbit() are SYMMETRIC about gamma and")
        print("  are UNAFFECTED, so RYA-431's off-orbit NOT-ALPHA-CEN finding stands. What")
        print("  inverts is every A-vs-B call RYA-423 has made.")
        print()
        print("  Proposed one-line change (NOT applied here - ratify first):")
        print("      _OMEGA_A = np.radians(OMEGA_DEG) + np.pi")
        print("      _OMEGA_B = np.radians(OMEGA_DEG)")
        return 0
    print("  VERDICT: the two tests do not agree - do not act on either alone.")
    return 1


if __name__ == '__main__':
    raise SystemExit(main())
