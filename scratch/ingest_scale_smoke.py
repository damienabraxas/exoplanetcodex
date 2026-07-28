#!/usr/bin/env python3
"""
scratch/ingest_scale_smoke.py
=============================
RYA-264 §3 smoke — the wavelength unit+scale axis end-to-end. Demonstrates that
`frame_object_contract.to_air_angstrom` (the 3rd loader-contract axis) converts each
instrument's native unit/scale to air Å from the CITED registry, runs the
unit-sanity band gate (the ×10/×10000 catch), and RAISES — never returns silent
zeros or a defaulted scale — on an unknown instrument or a mis-united array.

Deterministic checks run with NO data. Pass real frames to also exercise the loaders:

    python scratch/ingest_scale_smoke.py [--spirou <a SPIRou _t.fits>] [--nirps <a NIRPS S1D>]
"""
import argparse
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from pipeline import frame_object_contract as foc          # noqa: E402


def _ok(cond):
    return 'OK' if cond else 'FAIL'


def deterministic():
    ok = True
    print("== RYA-264 wavelength unit+scale axis — deterministic checks ==\n")

    # SPIRou: native vacuum nm -> air Å (×10 + vac→air)
    w_nm = np.array([980.0, 1500.0, 2400.0])
    air = foc.to_air_angstrom(w_nm, 'SPIRou')
    off = air - w_nm * 10.0
    spirou_ok = bool(np.all(off < 0) and np.all(np.abs(off) > 2.0) and np.all(np.abs(off) < 8.0))
    ok &= spirou_ok
    print(f"[SPIRou]  vacuum nm -> air Å : {np.round(air, 2)}")
    print(f"          vac–air offset Å  : {np.round(off, 2)} (air < vacuum, ~3–7 Å)  {_ok(spirou_ok)}")

    # HARPS / UVES: declared air, no conversion
    uves = foc.to_air_angstrom(np.array([4000.0, 6000.0, 9000.0]), 'UVES')
    harps = foc.to_air_angstrom(np.array([4000.0, 6000.0]), 'HARPS')
    noop_ok = bool(np.allclose(uves, [4000, 6000, 9000]) and np.allclose(harps, [4000, 6000]))
    ok &= noop_ok
    print(f"[UVES]    air Å (declared, no conversion) : {uves}  {_ok(noop_ok)}")
    print(f"[HARPS]   air Å (declared, no conversion) : {harps}")

    # ESPRESSO / NIRPS: native vacuum (registry-corrected from the hypothesis table)
    for inst in ('ESPRESSO', 'NIRPS'):
        c = foc.WAVELENGTH_CONVENTION[inst]
        print(f"[{inst}] native {c.native_unit}/{c.native_scale}  «{c.citation[:70]}…»")

    # scale sign-check: Mg II k 2796.3543 vac -> 2795.528 air
    resid = foc.verify_vac_to_air(2796.3543, 2795.528, tol_A=0.01)
    sign_ok = abs(resid) < 0.01
    ok &= sign_ok
    print(f"\n[sign]    Mg II k 2796.3543(vac) -> air, residual {resid:+.4f} Å (tol 0.01)  {_ok(sign_ok)}")

    # LOUD: unknown instrument raises (no default)
    try:
        foc.to_air_angstrom(w_nm, 'BOGUS'); raised = False
    except foc.FrameContractError:
        raised = True
    ok &= raised
    print(f"[loud]    unknown instrument -> raises (no default)  {_ok(raised)}")

    # LOUD: unit-sanity gate catches nm-as-Å (×10 too small)
    try:
        foc.to_air_angstrom(w_nm / 10.0, 'SPIRou'); gate = False
    except foc.FrameContractError:
        gate = True
    ok &= gate
    print(f"[loud]    nm-array read as Å -> band gate raises  {_ok(gate)}")

    print(f"\nDeterministic: {'ALL OK' if ok else 'FAILURES PRESENT'}")
    return ok


def real_spirou(path):
    from pipeline.loaders.spirou_loader import SPIRouLoader
    print(f"\n== SPIRou loader on {Path(path).name} ==")
    s = SPIRouLoader(path).load()
    lo, hi = s.wave_range_A
    print(f"  range = {lo:.1f}–{hi:.1f} Å   wavelength_scale = {s.meta['wavelength_scale']!r}")
    band_ok = 8800.0 <= lo and hi <= 26000.0
    print(f"  air-Å band ~9800–24000 (with margin): {_ok(band_ok)}")
    return band_ok


def real_nirps(path):
    print(f"\n== NIRPS frame {Path(path).name} ==")
    # No NIRPS loader yet (RYA-498 builds it); demonstrate the axis on its WAVE grid.
    from astropy.io import fits
    with fits.open(path) as h:
        wave = None
        for hdu in h:
            if hdu.data is not None and getattr(hdu, 'columns', None) is not None:
                for col in ('WAVE_AIR', 'WAVE', 'wavelength'):
                    if col in hdu.columns.names:
                        wave = np.asarray(hdu.data[col]).ravel().astype(float)
                        provided_air = (col == 'WAVE_AIR')
                        print(f"  using column {col!r} (provided_air={provided_air})")
                        break
            if wave is not None:
                break
    if wave is None:
        print("  no WAVE/WAVE_AIR column found — skipped")
        return True
    ws = foc.WavelengthScale('NIRPS', native_unit='A', native_scale='vacuum',
                             provided_air=provided_air)
    air = ws.to_air_angstrom(wave[np.isfinite(wave)])
    print(f"  air-Å range = {air.min():.1f}–{air.max():.1f} (YJH ~9724–19196)")
    return True


def main():
    ap = argparse.ArgumentParser(description='RYA-264 wavelength unit+scale smoke')
    ap.add_argument('--spirou', help='a SPIRou _t.fits')
    ap.add_argument('--nirps', help='a NIRPS S1D')
    args = ap.parse_args()
    ok = deterministic()
    if args.spirou:
        ok &= real_spirou(args.spirou)
    if args.nirps:
        ok &= real_nirps(args.nirps)
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
