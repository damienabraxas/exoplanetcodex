#!/usr/bin/env python3
"""RYA-925: visual diagnosis of the anomalous Kitt Peak Al I 3057.144 A fit."""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "scripts")]

from pipeline.nearuv_synth import build_solar_context, synthesize_band
from rya759_nearuv_fe_product import NEARUV_LINELIST
from rya759_nearuv_synth import _kp_segments, _load_kp_window

CENTRE = 3057.144
FIT_HALF_WIDTH = 0.40
PLOT_HALF_WIDTH = 0.72
BEST_A = 4.198
LITERATURE_A = 6.43
NO_AL_A = 1.43
OUT = ROOT / "data" / "results" / "rya925"


def main() -> None:
    ctx = build_solar_context("Al", 500_000, linelist_file=str(NEARUV_LINELIST),
                              apply_canonical_gf=False)
    ow, obs = _load_kp_window(_kp_segments(), CENTRE, PLOT_HALF_WIDTH)
    curves = {}
    for label, abundance in (("best_fit", BEST_A), ("literature", LITERATURE_A),
                             ("effectively_no_Al", NO_AL_A)):
        curves[label] = synthesize_band(
            ctx, CENTRE - PLOT_HALF_WIDTH, CENTRE + PLOT_HALF_WIDTH,
            element="Al", trial_A=abundance, step_A=0.002,
            check_sensitivity=False)

    wave = curves["best_fit"]["wave_A"]
    obs_i = np.interp(wave, ow, obs)
    table = pd.DataFrame({"wave_A": wave, "observed_kpno": obs_i})
    for label, curve in curves.items():
        table[f"synth_{label}"] = curve["flux"]
    table["al_signal_best_minus_no_al"] = (
        table["synth_best_fit"] - table["synth_effectively_no_Al"])
    table.to_csv(OUT / "AlI_3057_nearuv_diagnostic_spectra.csv", index=False)

    ll = ctx["linelist"]
    m = ((np.asarray(ll["wave_A"], float) >= CENTRE - PLOT_HALF_WIDTH) &
         (np.asarray(ll["wave_A"], float) <= CENTRE + PLOT_HALF_WIDTH))
    fields = [x for x in ("element", "wave_A", "lower_state_eV", "loggf",
                           "theoretical_depth", "reference_code")
              if x in ll.dtype.names]
    neighbours = pd.DataFrame({f: ll[f][m] for f in fields}).sort_values("wave_A")
    neighbours.to_csv(OUT / "AlI_3057_nearuv_neighbours.csv", index=False)

    fit = (wave >= CENTRE - FIT_HALF_WIDTH) & (wave <= CENTRE + FIT_HALF_WIDTH)
    chi = {}
    for label in curves:
        resid = (obs_i[fit] - table.loc[fit, f"synth_{label}"].to_numpy()) / 0.01
        chi[label] = float(np.sum(resid * resid) / max(resid.size - 1, 1))

    fig, (ax, rx) = plt.subplots(2, 1, figsize=(12, 7.2), sharex=True,
                                 gridspec_kw={"height_ratios": [3, 1]})
    ax.plot(wave, obs_i, color="black", lw=1.25, label="Kitt Peak observed")
    ax.plot(wave, table.synth_best_fit, color="#268bd2", lw=1.3,
            label=f"1D-LTE best A(Al)={BEST_A:.3f}  (chi2r={chi['best_fit']:.1f})")
    ax.plot(wave, table.synth_literature, color="#2aa198", lw=1.2,
            label=f"1D-LTE A(Al)={LITERATURE_A:.2f}  (chi2r={chi['literature']:.1f})")
    ax.plot(wave, table.synth_effectively_no_Al, color="#b58900", lw=1.1, ls="--",
            label=f"effectively no Al  (chi2r={chi['effectively_no_Al']:.1f})")
    ax.axvspan(CENTRE - FIT_HALF_WIDTH, CENTRE + FIT_HALF_WIDTH, color="#268bd2",
               alpha=.06, label="fitted window")
    ax.axvline(CENTRE, color="#dc322f", ls=":", lw=1.2)
    ax.set_ylabel("Residual flux")
    ax.set_title("RYA-925 · Kitt Peak Al I 3057.144 Å near-UV diagnostic")
    ax.legend(loc="lower left", fontsize=8, ncol=2)

    if not neighbours.empty:
        depth_col = "theoretical_depth" if "theoretical_depth" in neighbours else None
        ranked = neighbours.sort_values(depth_col, ascending=False).head(12) if depth_col else neighbours.head(12)
        y0, y1 = ax.get_ylim()
        for j, row in enumerate(ranked.itertuples(index=False)):
            w = float(getattr(row, "wave_A"))
            species = str(getattr(row, "element", "?"))
            ax.vlines(w, y1 - .045 * (y1 - y0), y1, color="#6c7a96", lw=.7)
            ax.text(w, y1 - (.055 + .055 * (j % 3)) * (y1 - y0), species,
                    rotation=90, va="top", ha="center", fontsize=6, color="#4f5b73")

    rx.axhline(0, color="#6c7a96", lw=.8)
    rx.plot(wave, obs_i - table.synth_best_fit, color="#268bd2", lw=1,
            label="obs - best")
    rx.plot(wave, obs_i - table.synth_literature, color="#2aa198", lw=1,
            label="obs - literature")
    rx.fill_between(wave, -0.01, 0.01, color="#6c7a96", alpha=.12,
                    label="±0.01 adequacy floor")
    rx.axvline(CENTRE, color="#dc322f", ls=":", lw=1.2)
    rx.set_ylabel("Residual")
    rx.set_xlabel("Air wavelength (Å)")
    rx.legend(loc="lower left", fontsize=8, ncol=3)
    fig.tight_layout()
    fig.savefig(OUT / "AlI_3057_nearuv_diagnostic.png", dpi=180)

    peak_al_signal = float(np.max(np.abs(
        table.synth_best_fit - table.synth_effectively_no_Al)))
    rms_best = float(np.sqrt(np.mean((obs_i[fit] - table.loc[fit, "synth_best_fit"])**2)))
    rms_lit = float(np.sqrt(np.mean((obs_i[fit] - table.loc[fit, "synth_literature"])**2)))
    print(f"neighbours={len(neighbours)} peak_Al_signal={peak_al_signal:.6f} "
          f"rms_best={rms_best:.6f} rms_literature={rms_lit:.6f} chi={chi}")


if __name__ == "__main__":
    main()
