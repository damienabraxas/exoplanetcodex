#!/usr/bin/env python3
"""RYA-770 — RCA for the SynthesisHandler adapter failing the optical control.

    python3 scripts/rya770_adapter_rca.py --n 12

DISPOSABLE DIAGNOSTIC. Reads only; changes no production path and writes no science
artifact. It exists to answer one question with numbers instead of argument:

    the adapter measured 3 of 12 lines (7 CHI2-GATE with red_chi2 16-302, 2
    NOT-IN-SYNTH-LINELIST) while the SAME core fit, on the SAME spectrum, is banked at
    med red_chi2 = 2.976. What differs?

The prime hypothesis on the ticket is iSpec's edge-zeroing: it sets the first and last
pixel of every synthesis window to ~1e-10. That is checked FIRST and explicitly, because
if true it is a two-pixel artefact producing a 300x chi2 and everything else is noise.

Per control line this reports:
  * EDGE       — are sf[0] / sf[-1] ~1e-10 on the core's fine grid, and (the part that
                 matters) do any OBSERVED pixels land in the zeroed span?
  * ALIGN      — offset between the observed and synthetic line-core wavelength. A
                 rest-frame slip shows up here and inflates chi2 without bounding it.
  * CONTINUUM  — observed flux at the window edges. A window sitting on a blend rather
                 than continuum inflates chi2 the same way and is NOT a defect.
  * CHI2       — recomputed with the edge pixels included vs excluded, which separates
                 "edges did it" from "edges are innocent" with one number.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.control_synthesis_handler import (          # noqa: E402
    build_context, CONTROL_LO, CONTROL_HI)

EDGE_ZERO = 1e-6          # anything below this is iSpec's ~1e-10 sentinel, not flux
SIG = 0.01                # the core's model-adequacy sigma


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--element', default='Fe')
    ap.add_argument('--ion', default='I')
    ap.add_argument('--n', type=int, default=12)
    ap.add_argument('--out', default=str(ROOT / 'data' / 'audit' / 'synthesis_control'))
    a = ap.parse_args()

    from pipeline.abundances_derive import (
        _fit_synth_flux, _synth_flux_at_abund, _wingwide_window_nm,
        _load_observed_spectrum, SYNTH_FIT_WSTEP_NM)
    from config.constants import HARPS_R

    # Same 12 lines the control uses -- same filter, same linspace subsample.
    pool = pd.read_csv(ROOT / 'data' / 'measured' / 'sol_ew_results_v1.csv')
    ref = pool[(pool.element == a.element) & (pool.ion == a.ion) &
               pool.wavelength_air_A.between(CONTROL_LO, CONTROL_HI) &
               pool.ew_mA.between(10, 90)].copy()
    ref = ref.sort_values('wavelength_air_A').reset_index(drop=True)
    ref = ref.iloc[np.unique(np.linspace(0, len(ref) - 1, a.n).astype(int))]
    ref = ref.reset_index(drop=True)

    ctx = build_context(a.element, a.ion, float(HARPS_R))
    ow_nm, oflux = _load_observed_spectrum('solar')
    ow_A = ow_nm * 10.0
    print(f"HARPS normalised solar: {ow_A.size} pts, {ow_A.min():.1f}-{ow_A.max():.1f} A")
    print(f"broadening: R={HARPS_R}, vmac={ctx['macroturbulence']}, vsini={ctx['vsini']}\n")

    kw = dict(atmosphere=ctx['atmosphere'], teff=ctx['teff'], logg=ctx['logg'],
              feh=ctx['feh'], vturb=ctx['vturb'], linelist=ctx['linelist'],
              isotopes=ctx['isotopes'], solar_abund=ctx['solar_abund'],
              element=a.element, atom_code=ctx['atom_code'], R=float(HARPS_R),
              macroturbulence=float(ctx['macroturbulence']),
              vsini=float(ctx['vsini']), tmp_dir='/tmp/ispec_codex_synth')

    rows = []
    for _, r in ref.iterrows():
        c = float(r.wavelength_air_A)
        ew_obs = float(r.ew_mA)
        wbase, wtop = _wingwide_window_nm(c / 10.0, ew_obs)

        fit = _fit_synth_flux(ow_nm, oflux, ctx['atmosphere'], ctx['teff'], ctx['logg'],
                              ctx['feh'], ctx['vturb'], ctx['linelist'], ctx['isotopes'],
                              ctx['solar_abund'], a.element, ctx['atom_code'],
                              wbase, wtop, 4.0, 10.0, float(HARPS_R),
                              float(ctx['macroturbulence']), float(ctx['vsini']),
                              tmp_dir='/tmp/ispec_codex_synth')
        a_best = fit.get('A_X')
        if a_best is None or not np.isfinite(a_best):
            rows.append(dict(wave=c, status=fit['status'], note=fit.get('reason', '')))
            print(f"{c:10.3f}  {fit['status']}: {fit.get('reason','')}")
            continue

        # Reproduce the core's internals so the residual can be dissected. The step comes
        # from the core, not a local literal: a diagnostic that measures the zeroed span
        # on a DIFFERENT grid than the fit uses would be measuring the wrong thing.
        wstep = float(SYNTH_FIT_WSTEP_NM)
        sw = np.arange(wbase, wtop + wstep * 0.5, wstep)
        sf = _synth_flux_at_abund(sw, trial_A=float(a_best), **kw)

        m = (ow_nm >= wbase) & (ow_nm <= wtop)
        ow_w, of_w = ow_nm[m], oflux[m]

        # ── EDGE: is iSpec zeroing, and does it reach any OBSERVED pixel? ──────
        n_zero = int(np.sum(sf < EDGE_ZERO))
        zero_idx = np.flatnonzero(sf < EDGE_ZERO)
        edges_zeroed = bool(n_zero and (0 in zero_idx or len(sf) - 1 in zero_idx))
        # the observed pixels only ever see sf through interp1d, so what matters is
        # whether an observed pixel falls inside a zeroed SPAN
        obs_in_zero = 0
        if n_zero:
            for zi in zero_idx:
                lo = sw[max(zi - 1, 0)]
                hi = sw[min(zi + 1, len(sw) - 1)]
                obs_in_zero += int(np.sum((ow_w >= lo) & (ow_w <= hi)))

        from scipy.interpolate import interp1d
        sf_i = interp1d(sw, sf, bounds_error=False, fill_value=1.0)(ow_w)
        resid = (of_w - sf_i) / SIG
        chi2_all = float(np.nansum(resid ** 2)) / max(of_w.size - 1, 1)

        # chi2 with any observed pixel that sees a zeroed span removed
        keep = np.ones(ow_w.size, dtype=bool)
        if n_zero:
            for zi in zero_idx:
                lo = sw[max(zi - 1, 0)]
                hi = sw[min(zi + 1, len(sw) - 1)]
                keep &= ~((ow_w >= lo) & (ow_w <= hi))
        chi2_trim = (float(np.nansum(resid[keep] ** 2)) / max(int(keep.sum()) - 1, 1)
                     if keep.sum() > 2 else float('nan'))

        # ── ALIGN: observed vs synthetic core position ────────────────────────
        align_mA = float('nan')
        try:
            align_mA = (ow_w[np.argmin(of_w)] - sw[np.argmin(sf)]) * 1e7  # nm -> mA
        except Exception:
            pass

        # ── CONTINUUM: what the window edges actually sit on ──────────────────
        cont_lo = float(np.median(of_w[:3])) if of_w.size >= 3 else float('nan')
        cont_hi = float(np.median(of_w[-3:])) if of_w.size >= 3 else float('nan')

        rows.append(dict(
            wave=c, ew_ref=ew_obs, A=a_best, red_chi2=fit['red_chi2'],
            n_pix=int(of_w.size), status=fit['status'],
            n_synth_zero=n_zero, edges_zeroed=edges_zeroed, obs_pix_in_zero=obs_in_zero,
            chi2_all=round(chi2_all, 3), chi2_edge_trimmed=round(chi2_trim, 3),
            align_mA=round(align_mA, 2), obs_min_flux=round(float(of_w.min()), 4),
            synth_min_flux=round(float(sf.min()), 4),
            cont_edge_lo=round(cont_lo, 4), cont_edge_hi=round(cont_hi, 4)))

        print(f"{c:10.3f}  A={a_best:6.3f}  chi2={fit['red_chi2']:8.2f}  "
              f"npix={of_w.size:3d}  synth_zeros={n_zero:3d} "
              f"edges={'Y' if edges_zeroed else 'n'}  obs_in_zero={obs_in_zero:2d}  "
              f"chi2_trim={chi2_trim:8.2f}  align={align_mA:+7.1f} mA  "
              f"obs_min={of_w.min():.3f} synth_min={sf.min():.3f}  "
              f"edgeflux={cont_lo:.3f}/{cont_hi:.3f}")

    d = pd.DataFrame(rows)
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    d.to_csv(out / 'rya770_adapter_rca.csv', index=False)

    print("\n" + "=" * 78)
    print("EDGE-ZEROING HYPOTHESIS")
    if 'obs_pix_in_zero' in d:
        tot = int(d.obs_pix_in_zero.fillna(0).sum())
        print(f"  synthetic pixels at ~1e-10, summed over lines : "
              f"{int(d.n_synth_zero.fillna(0).sum())}")
        print(f"  OBSERVED pixels that see a zeroed span        : {tot}")
        both = d.dropna(subset=['chi2_all', 'chi2_edge_trimmed'])
        if len(both):
            print(f"  median red_chi2 as-is                         : "
                  f"{both.chi2_all.median():.2f}")
            print(f"  median red_chi2 with those pixels removed     : "
                  f"{both.chi2_edge_trimmed.median():.2f}")
        print(f"  VERDICT: {'EDGES IMPLICATED' if tot else 'EDGES ARE NOT THE ROOT'}"
              f" (no observed pixel lands in a zeroed span)" if not tot else "")
    print("\nALIGNMENT (rest-frame slip would show as a consistent non-zero offset)")
    if 'align_mA' in d:
        print(d[['wave', 'align_mA', 'red_chi2']].to_string(index=False))
    print(f"\nwrote {out / 'rya770_adapter_rca.csv'}")


if __name__ == '__main__':
    main()
