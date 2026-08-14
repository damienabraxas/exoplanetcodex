"""RYA-783 — decompose the Engine-B-NLTE failure THROUGH THE HANDLER, as the driver does.

The previous probe called `_fit_synth_flux` directly with a hand-built window and every arm
failed on "too few observed" — a fault in the probe, not a result. This one uses exactly
what derive_band_products uses, so the only differences between arms are the two the
gerber-nlte path actually changes: the atmosphere, and the departures.
"""
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
# Standalone-script bootstrap (RYA-313): repo root on sys.path BEFORE importing
# config/pipeline, so this runs from any cwd. Derived from __file__, never cwd.
import os as _os_boot, sys as _sys_boot
_sys_boot.path.insert(0, _os_boot.path.dirname(_os_boot.path.dirname(
    _os_boot.path.abspath(__file__))))
from config.constants import codex_path  # RYA-810 path register

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, os.environ.get("ISPEC_DIR", str(codex_path('engines.ispec'))))


def main() -> None:
    from scripts.control_synthesis_handler import build_context
    from scripts.measure_band_ew import kp_segments, load_kp_window
    from pipeline.abundances_derive import _load_atmosphere
    from pipeline.measure import resolve_handler
    from pipeline.band_policy import resolve as policy_for
    from config.constants import STAR_PARAMS

    p = STAR_PARAMS["solar"]
    teff, logg = float(p["teff"]), float(p["logg"])
    feh = float(p.get("feh", p.get("feh_ref", 0.0)))
    ctx = build_context("Fe", "II", 500000.0)
    marcs = _load_atmosphere(teff, logg, feh, ctx["vturb"], model_grid="MARCS.GES")
    pol = policy_for(7500.0)
    segs = kp_segments()

    ew = pd.read_csv(ROOT / "data" / "measured" / "band_ew" /
                     "FeII_6910_9199_kpno_solar_atlas_PROFILEFIT_ew.csv")
    ok = ew[ew.in_aggregate.fillna(False)]

    arms = [("ATLAS9+LTE", dict(ctx), None),
            ("MARCS+LTE", {**ctx, "atmosphere": marcs}, None),
            ("MARCS+NLTE", {**ctx, "atmosphere": marcs, "nlte_deck": "gerber"}, "gerber")]

    for name, c_ctx, _deck in arms:
        handler = resolve_handler(3400.0)
        handler.prepare(pol, {**c_ctx, "instrument": "kpno_solar_atlas"})
        print(f"\n=== {name} ===")
        for _, r in ok.iterrows():
            c = float(r.wavelength_air_A)
            try:
                w_obs, f_obs, _ = load_kp_window(segs, c, pad=1.4)
            except Exception as e:
                print(f"  {c:9.3f}  window: {type(e).__name__}")
                continue
            lb = handler.measure_line(w_obs, f_obs, element="Fe", ion="II",
                                      wavelength_A=c, instrument="kpno_solar_atlas",
                                      policy=pol, pre_normalised=True,
                                      context={**c_ctx, "ew_hint_mA": float(r.ew_mA)})
            if lb.in_aggregate and lb.abundance is not None:
                print(f"  {c:9.3f}  A={lb.abundance:.4f}")
            else:
                print(f"  {c:9.3f}  {str(lb.excluded_reason)[:110]}")


if __name__ == "__main__":
    main()
