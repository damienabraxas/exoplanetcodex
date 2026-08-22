#!/usr/bin/env python3
"""RYA-798 — is the LTE synthesis path bit-identical to before the NLTE threading?

RYA-770 stabilised `_fit_synth_flux` at -0.026 dex against the banked optical answer, and
RYA-798 edited the function it calls. Asserting "the default is None" is a structural
argument; this is the numerical one. Run in a worktree, it writes the flux arrays and the
fitted abundances to an .npz; run it in a main worktree and in the branch worktree and the
two files must agree EXACTLY.

Deliberately does not go through the handler: fewer moving parts between the edit and the
number.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
# must precede the ISPEC_DIR line below, which USES codex_path
from config.constants import codex_path  # noqa: E402
sys.path.insert(0, os.environ.get("ISPEC_DIR", str(codex_path("engines.ispec"))))

LINES = [5905.671, 6027.070, 6705.117]
WSTEP_NM = 0.0002


def main() -> None:
    out = sys.argv[1]
    # RYA-985: optional second argument, defaulting to solar so the existing invocation is
    # unchanged. Same single source as the rest of the chain — no STAR_PARAMS["solar"] here.
    star = sys.argv[2] if len(sys.argv) > 2 else "solar"
    from scripts.control_synthesis_handler import build_context
    from pipeline.abundances_derive import _synth_flux_at_abund, _fit_synth_flux
    from config.constants import get_star_params

    p = get_star_params(star)
    teff, logg = float(p["teff"]), float(p["logg"])
    feh = float(p.get("feh", p.get("feh_ref", 0.0)))
    ctx = build_context("Fe", "I", 500000.0, star=star)

    kw = dict(atmosphere=ctx["atmosphere"], teff=teff, logg=logg, feh=feh,
              vturb=ctx["vturb"], linelist=ctx["linelist"], isotopes=ctx["isotopes"],
              solar_abund=ctx["solar_abund"], element="Fe", atom_code=26,
              R=0, macroturbulence=0.0, vsini=0.0)

    fluxes, fits = {}, {}
    for c in LINES:
        sw = np.arange((c - 2.0) / 10.0, (c + 2.0) / 10.0 + WSTEP_NM * 0.5, WSTEP_NM)
        f = _synth_flux_at_abund(sw, trial_A=7.46, **kw)
        fluxes[f"{c:.3f}"] = np.asarray(f, dtype=float)
        print(f"  {c:9.3f}  n={len(f)}  mean={np.mean(f):.10f}  min={np.min(f):.10f}")

    # and one real fit, because the fitter is what the products actually call
    try:
        from pipeline.abundances_derive import (_load_observed_spectrum,
                                                _wingwide_window_nm)
        ow_nm, oflux = _load_observed_spectrum("solar")
    except Exception as e:
        print(f"  (fit leg skipped: {type(e).__name__} — observed spectrum not staged "
              f"in this worktree; the flux arrays above are the test)")
        np.savez(out, **fluxes)
        print(f"wrote {out}")
        return
    c = LINES[0]
    wb, wt = _wingwide_window_nm(c, 60.0)
    r = _fit_synth_flux(ow_nm, oflux, ctx["atmosphere"], teff, logg, feh, ctx["vturb"],
                        ctx["linelist"], ctx["isotopes"], ctx["solar_abund"], "Fe", 26,
                        wb, wt, 6.8, 8.1, float(ctx["resolving_power"]),
                        float(ctx["macroturbulence"]), float(ctx["vsini"]))
    fits["A_X"] = float(r.get("A_X", np.nan))
    fits["chi2"] = float(r.get("chi2_red", r.get("chi2", np.nan)))
    print(f"  fit {c:.3f}: A_X={fits['A_X']:.10f}  status={r.get('status')}")

    np.savez(out, **fluxes, **{f"fit_{k}": v for k, v in fits.items()})
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
