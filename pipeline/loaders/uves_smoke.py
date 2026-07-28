"""
pipeline/loaders/uves_smoke.py
===============================
RYA-272 Step 4 — UVES loader smoke test. Four checks:
  1. Quarantine listing (8 quarantined, 45 science IDPs remain).
  2. Guard test — a quarantined GES (HELIOCEN) file must be REJECTED loudly.
  3. BERV SIGN CHECK (CRITICAL) — load ≥3 epochs spanning the BERV range; measure the
     Hα core in the BERV-CORRECTED spectrum. The epoch-dependent BERV signature must be
     GONE: corrected velocities cluster at Procyon's systemic RV, spread << the raw BERV
     spread. Wrong sign would DOUBLE the spread → fail.
  4. Registry — 45 rows, exactly one oi_anchor, O I verdicts present.

Run:  python -m pipeline.loaders.uves_smoke
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from config.constants import PHYSICS                                       # noqa: E402
from pipeline.loaders.uves_loader import UVESLoader, UVESProductError      # noqa: E402

PROC = ROOT.parent / 'data' / 'spectra' / 'exoplanetcodex-data' / 'Procyon'
UVES_DIR = PROC / 'Procyon UVES'
QUAR = PROC / 'quarantine'
REGISTRY = ROOT / 'data' / 'spectra' / 'procyon' / 'uves_registry.csv'
_C = PHYSICS['c_kms']
HALPHA = 6562.79          # air Å

# one Hα-covering IDP per epoch, spanning the BERV range (−25 … +28 km/s)
BERV_CHECK_FILES = [
    'ADP.2020-07-17T12:58:31.921.fits',   # 2005-03-18  BERV ≈ −25.2
    'ADP.2021-08-27T05:08:19.964.fits',   # 2002-02-27  BERV ≈ −19.3
    'ADP.2020-08-10T12:01:00.198.fits',   # 2002-10-08  BERV ≈ +28.4
    'ADP.2020-06-15T10:09:57.908.fits',   # 2013-10-08  BERV ≈ +28.5  (anchor)
]


def _halpha_velocity(spec) -> float:
    """Flux-weighted centroid of the Hα absorption core in the (corrected) spectrum →
    radial velocity (km/s). Robust for the broad F-star Balmer core (~3 km/s scatter)."""
    w, f = spec.wave_A, spec.flux
    m = (w > HALPHA - 3.0) & (w < HALPHA + 3.0) & np.isfinite(f)
    w, f = w[m], f[m]
    cont = np.nanpercentile(f, 90)
    absorp = np.clip(cont - f, 0.0, None)
    centroid = float(np.sum(w * absorp) / np.sum(absorp))
    return _C * (centroid - HALPHA) / HALPHA


def main() -> int:
    ok = True
    print(f"\n{'='*70}\n  RYA-272 UVES loader smoke test\n{'='*70}")

    # ── 1. quarantine listing ──
    qfits = sorted(QUAR.glob('ADP.*.fits'))
    sci = sorted(UVES_DIR.glob('ADP.*.fits'))
    print(f"\n[1] Quarantine: {len(qfits)} files quarantined, {len(sci)} science IDPs remain.")
    for q in qfits:
        print(f"      quarantined: {q.name}")
    ok &= (len(qfits) == 8 and len(sci) == 45)
    print(f"    -> {'PASS' if len(qfits) == 8 and len(sci) == 45 else 'FAIL'} "
          f"(expect 8 quarantined / 45 science)")

    # ── 2. guard test ──
    ges = next((q for q in qfits if '2020-12-07' in q.name), None)
    print(f"\n[2] Guard test — loading a quarantined GES file ({ges.name}):")
    try:
        UVESLoader(ges).load()
        print("    -> FAIL: GES file loaded (guard did not fire)")
        ok = False
    except UVESProductError as e:
        print(f"    REJECTED: {e}")
        print("    -> PASS")

    # ── 3. BERV sign check (CRITICAL) ──
    print(f"\n[3] BERV SIGN CHECK (critical):")
    print(f"    {'file':38s} {'date':11s} {'BERV':>8s} {'v_raw':>8s} {'v_corr':>8s}")
    raw_v, corr_v, bervs = [], [], []
    for fn in BERV_CHECK_FILES:
        spec = UVESLoader(UVES_DIR / fn).load()
        berv = spec.meta['berv_kms']
        v_corr = _halpha_velocity(spec)          # corrected (loader applied BERV)
        v_raw = v_corr - berv                     # back out to topocentric for display
        raw_v.append(v_raw); corr_v.append(v_corr); bervs.append(berv)
        print(f"    {fn:38s} {spec.meta['date_obs'][:10]:11s} {berv:+8.2f} {v_raw:+8.2f} {v_corr:+8.2f}")
    raw_spread = float(np.max(raw_v) - np.min(raw_v))
    corr_spread = float(np.max(corr_v) - np.min(corr_v))
    sys_rv = float(np.mean(corr_v))
    print(f"    raw (topocentric) Hα spread = {raw_spread:.1f} km/s  (= the BERV signature)")
    print(f"    corrected Hα spread         = {corr_spread:.1f} km/s  → systemic RV ~ {sys_rv:+.1f} km/s")
    sign_ok = corr_spread < 0.4 * raw_spread and corr_spread < 8.0
    print(f"    -> {'PASS' if sign_ok else 'FAIL — sign wrong, STOP'} "
          f"(corrected spread must collapse vs raw {raw_spread:.0f} km/s)")
    ok &= sign_ok

    # ── 4. registry ──
    print(f"\n[4] Registry ({REGISTRY.relative_to(ROOT)}):")
    if not REGISTRY.exists():
        print("    -> FAIL: registry not built (run scripts/build_uves_registry_rya272.py)")
        return 1
    df = pd.read_csv(REGISTRY)
    n_anchor = int(df['oi_anchor'].sum())
    print(f"    {len(df)} rows, {int(df['covers_OI_7771'].sum())} cover O I, {n_anchor} oi_anchor")
    print(df[df['covers_OI_7771']][['filename', 'date_obs', 'setting', 'snr',
          'berv_applied_kms', 'oi_telluric_verdict', 'oi_anchor']].to_string(index=False))
    verds = set(df['oi_telluric_verdict'].dropna().astype(str).str.strip()) - {''}
    reg_ok = (len(df) == 45 and n_anchor == 1 and
              verds <= {'CLEAN', 'CORRECTABLE', 'EXCLUDE'})
    ok &= reg_ok
    print(f"    -> {'PASS' if reg_ok else 'FAIL'} (45 rows / 1 anchor / verdicts ∈ audit set)")

    print(f"\n{'='*70}\n  SMOKE TEST: {'PASS' if ok else 'FAIL'}\n{'='*70}\n")
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
