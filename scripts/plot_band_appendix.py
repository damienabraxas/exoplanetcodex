#!/usr/bin/env python3
"""Appendix plots — one panel per line, so a reader can check the claim — RYA-711/713.

    python3 scripts/plot_band_appendix.py --element Fe --ion I --lo 6910 --hi 9199

Ryan, 2026-08-09: *"if that line failed, you would still have to prove it to me in the
Appendix, which is a check in itself. From a plot we can determine, hey do we have the
right line? Is it normalized correctly? Is there a blend? is this a gf ghost?"*

So the panel is not decoration — it is the evidence, and it is drawn for **every** line,
passed or quarantined. A quarantined line without a plot is an assertion; with one it is
a demonstration, and the reader can overrule us.

Each panel is built to answer Ryan's four questions directly:

  1. RIGHT LINE?     catalogued position (dashed) vs the deepest feature actually found
                     (marker + offset in mA of a Angstrom). If they diverge, the window
                     is dominated by something else.
  2. NORMALISED?     the continuum in use is drawn as a solid rule at 1.0, and the
                     side-band regions are shaded with their 95th percentile marked. If
                     that percentile sits well below the rule, the side-bands are
                     absorbed and any local re-normalisation would eat real line flux.
  3. BLEND?          catalogued neighbours are ticked. Absence of ticks is NOT absence of
                     blends -- our IR inventory is sparse -- so the profile itself is
                     shown wide enough that an asymmetric or double core is visible.
  4. GF GHOST?       the depth predicted from the line parameters is drawn as a short
                     horizontal bar against the observed depth. A line the catalogue says
                     is strong but the Sun does not show (or the reverse) is a ghost, and
                     the mismatch is the point of the bar.

Element, ion and band are arguments; nothing here knows which element it is drawing.
"""
from __future__ import annotations

import argparse
import sys
import textwrap
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.measure_band_ew import (  # noqa: E402
    kp_segments, load_kp_window, window_half_width, PRE_NORMALISED)
from pipeline.band_products import SIDEBAND_CLEAN_MIN  # noqa: E402

EW_DIR = ROOT / "data" / "measured" / "band_ew"
OUT = ROOT / "data" / "plots" / "band_appendix"

PANELS_PER_FIG = 12


def _text(v) -> str:
    """NaN-safe string. `float('nan') or ""` returns the NaN, not the empty string, so
    `str(x or "")` yields the literal 'nan' -- which is truthy and made three PASSING
    lines render as QUARANTINED in the prototype. Guard on isna, not on truthiness."""
    return "" if v is None or (isinstance(v, float) and pd.isna(v)) or pd.isna(v) else str(v)


def wrap(text: str, width: int = 62) -> str:
    """Wrap a title so a long root cause cannot overrun into the next panel."""
    return "\n".join(textwrap.wrap(text, width)) if text else ""


def short_reason(row) -> str:
    """One compact line of WHY, for the panel title."""
    r = _text(row.excluded_reason)
    if not r.strip():
        return "IN AGGREGATE — passed every check"
    if "FEATURE-VERIFICATION" in r:
        # Named for Ryan's four questions, so the panel title says WHICH check failed.
        if "GF-GHOST-ABSENT" in r:
            return "GF GHOST — catalogue promises a line, Sun shows none"
        if "GF-GHOST" in r:
            return "GF GHOST — right position, wrong strength"
        if "BLEND-DOMINATED" in r:
            return "BLEND — a different feature owns this window"
        return "QUARANTINED — failed feature verification"
    if "saturation ceiling" in r:
        return "QUARANTINED — saturated, inversion ill-conditioned"
    return "QUARANTINED"


def draw_panel(ax, w, f, row, hw, predicted_depth, pre_norm: bool) -> None:
    c = float(row.wavelength_air_A)
    cont = 1.0 if pre_norm else np.percentile(f, 95)

    ax.plot(w, f, lw=0.9, color="#1a1a1a", zorder=3)

    # (2) NORMALISED? -- the continuum actually used, and what the side-bands say.
    ax.axhline(cont, color="#0b7285", lw=1.1, zorder=2)
    d = np.abs(w - c)
    sb = (d > hw) & (d <= hw * 2.0)
    if sb.sum() >= 5:
        p95 = float(np.percentile(f[sb], 95))
        for sgn in (-1, 1):
            ax.axvspan(c + sgn * hw, c + sgn * hw * 2.0, color="#0b7285", alpha=0.07, zorder=1)
        ax.axhline(p95, color="#e8590c", lw=0.9, ls=":", zorder=2)
        if p95 < SIDEBAND_CLEAN_MIN:
            ax.text(0.02, 0.06, f"side-bands ABSORBED ({p95:.3f})", transform=ax.transAxes,
                    fontsize=6.5, color="#e8590c", weight="bold")

    # the integration window
    ax.axvspan(c - hw, c + hw, color="#ffd43b", alpha=0.16, zorder=1)

    # (1) RIGHT LINE? -- catalogued position vs the deepest feature actually present.
    ax.axvline(c, color="#c92a2a", lw=1.0, ls="--", zorder=4)
    m = d <= hw
    if m.sum():
        i = int(np.argmin(f[m]))
        peak, fmin = float(w[m][i]), float(f[m][i])
        ax.plot([peak], [fmin], marker="v", ms=5, color="#5f3dc4", zorder=5)
        off = peak - c
        if abs(off) > 0.02:
            ax.annotate(f"{off:+.3f} A", xy=(peak, fmin), xytext=(0, -11),
                        textcoords="offset points", ha="center", fontsize=6.5,
                        color="#5f3dc4", weight="bold")

        # (4) GF GHOST? -- predicted depth against what is AT THE LINE, which is the
        # quantity the ghost test uses. The deepest point in the window is a different
        # number and is shown separately; conflating them made panel 9179.742 read as
        # "obs 0.137" beside a verdict of "Sun shows none".
        if predicted_depth and predicted_depth > 0:
            ax.hlines(cont - predicted_depth * cont, c - hw * 0.42, c + hw * 0.42,
                      color="#2b8a3e", lw=2.0, zorder=5)
            j = int(np.argmin(np.abs(w - c)))
            at_line = 1.0 - float(f[j]) / cont
            deepest = 1.0 - fmin / cont
            ax.text(0.98, 0.06,
                    f"at line {at_line:.3f} / pred {predicted_depth:.3f}"
                    + (f"  (deepest {deepest:.3f})" if abs(deepest - at_line) > 0.02 else ""),
                    transform=ax.transAxes, fontsize=6.5, ha="right", color="#2b8a3e")

    ax.set_xlim(c - hw * 2.2, c + hw * 2.2)
    lo = min(f.min(), cont - 0.02)
    ax.set_ylim(lo - 0.04 * (cont - lo + 0.05), cont + 0.035)
    dom = _text(getattr(row, "fault_domain", ""))
    mech = _text(getattr(row, "mechanism", ""))
    cause = wrap(f"ROOT CAUSE [{dom}] {mech}") if dom else ""
    title = f"{row.element} {row.ion} {c:.3f}  ·  EW {row.ew_mA:.1f} mA\n" \
            f"{wrap(short_reason(row))}" + (f"\n{cause}" if cause else "")
    ax.set_title(title, fontsize=6.4, pad=4,
                 color=("#2b8a3e" if row.in_aggregate else "#c92a2a"))
    ax.tick_params(labelsize=6)
    # Absolute wavelengths -- an offset like "+7.296e3" hides the number that matters.
    ax.ticklabel_format(axis="x", useOffset=False, style="plain")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--element", required=True)
    ap.add_argument("--ion", default="I")
    ap.add_argument("--lo", type=float, required=True)
    ap.add_argument("--hi", type=float, required=True)
    ap.add_argument("--instrument", default="kpno_solar_atlas")
    ap.add_argument("--out", default=str(OUT))
    ap.add_argument("--waves", default=None,
                    help="comma-separated wavelengths for a focused figure (a PROTOTYPE "
                         "showing one example per failure mode). Omit for the full appendix.")
    ap.add_argument("--tag", default="", help="filename suffix for a focused figure")
    a = ap.parse_args()

    stem = f"{a.element}{a.ion}_{int(a.lo)}_{int(a.hi)}_{a.instrument}"
    df = pd.read_csv(EW_DIR / f"{stem}_ew.csv")
    acc = pd.read_csv(ROOT / "data" / "audit" / "line_accounting" / "per_line.csv")
    df = df.merge(acc[["wave_air_A", "predicted_depth", "log_gf"]],
                  left_on="wavelength_air_A", right_on="wave_air_A", how="left")
    # Root causes ride alongside, so a panel states WHERE the fault lives, not just
    # what the symptom looked like. QA wants the mechanism.
    rc_path = EW_DIR / f"{stem}_root_causes.csv"
    if rc_path.exists():
        rc = pd.read_csv(rc_path)
        if len(rc):
            df = df.merge(rc[["wave", "fault_domain", "mechanism"]],
                          left_on="wavelength_air_A", right_on="wave", how="left")
    allw = acc.wave_air_A.values
    segs = kp_segments()
    pre = PRE_NORMALISED[a.instrument]

    # Quarantined FIRST -- the failures are what the appendix exists to prove.
    df = df.sort_values(["in_aggregate", "wavelength_air_A"]).reset_index(drop=True)
    if a.waves:
        want = [float(x) for x in a.waves.split(",")]
        keep = [df.index[(df.wavelength_air_A - x).abs().argsort()[:1]][0] for x in want]
        missing = [x for x, k in zip(want, keep)
                   if abs(float(df.loc[k, "wavelength_air_A"]) - x) > 0.05]
        if missing:
            raise SystemExit(f"requested wavelengths not in this band's measurements: {missing}")
        df = df.loc[keep].reset_index(drop=True)

    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    written = []
    for start in range(0, len(df), PANELS_PER_FIG):
        chunk = df.iloc[start:start + PANELS_PER_FIG]
        nrow = int(np.ceil(len(chunk) / 3))
        fig, axes = plt.subplots(nrow, 3, figsize=(13, 3.35 * nrow), squeeze=False)
        for ax in axes.ravel():
            ax.axis("off")
        for k, (_, row) in enumerate(chunk.iterrows()):
            ax = axes.ravel()[k]; ax.axis("on")
            c = float(row.wavelength_air_A)
            hw = window_half_width(allw, c)
            try:
                w, f, _ = load_kp_window(segs, c, pad=hw * 3.0)
            except Exception as e:
                ax.text(0.5, 0.5, f"{c:.3f}\n{type(e).__name__}", ha="center",
                        va="center", fontsize=7, transform=ax.transAxes)
                continue
            draw_panel(ax, w, f, row, hw, float(row.predicted_depth or 0.0), pre)
        fig.suptitle(
            f"{a.element} {a.ion} · {a.instrument} · {a.lo:.0f}-{a.hi:.0f} A · "
            f"appendix panels {start + 1}-{start + len(chunk)} of {len(df)}\n"
            f"red dash = catalogued position · purple = deepest feature · yellow = integration "
            f"window · teal = continuum used · dotted orange = side-band 95th pct · green bar = predicted depth",
            fontsize=8.5)
        fig.tight_layout(rect=(0, 0, 1, 0.95), h_pad=2.4, w_pad=1.4)
        suffix = a.tag or f"appendix_{start // PANELS_PER_FIG + 1}"
        p = out / f"{stem}_{suffix}.png"
        fig.savefig(p, dpi=145); plt.close(fig)
        written.append(p)

    nq = int((~df.in_aggregate).sum())
    print(f"{a.element} {a.ion}: {len(df)} lines plotted ({nq} quarantined, shown first)")
    for p in written:
        print(f"  wrote {p}")


if __name__ == "__main__":
    main()
