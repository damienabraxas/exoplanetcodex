#!/usr/bin/env python3
"""
scripts/plot_unmeasured_lines_rya707.py — RYA-707
=================================================
A CANONICAL LINE THAT PRODUCES NO MEASUREMENT MUST BE ABLE TO SHOW WHY.

Ryan, 2026-08-08: *"I want a plot for any canonical line that doesn't give a
measurement. That way it is archived and proven."*

The gap this closes
-------------------
Every stage of the pipeline records a line's fate in words — `NO_EW_POOL`,
`GRADE`, `SATURATION_COG`, "line set exhausted". Those strings are conclusions.
Nothing in the repo let anyone check a conclusion against the spectrum, so a
wrong one survived every review it passed through.

Al is the worked case, and it was wrong twice. The RYA-673 wiring audit recorded
`NO_EW_POOL`; RYA-701 re-diagnosed it as a saturation artifact. Both are false.
Al has two lines in the committed pool, and this script's own numbers show why
neither is usable:

    window                     observed absorption   fitted EW      ratio
    6695.80-6696.60 A                66.48 mA        119.95 mA      1.8x
      (Al I 6696.020 + 6696.185 fitted independently, 0.165 A apart)
    6631.00-6631.45 A                 7.41 mA         34.53 mA      4.7x
      (Al I 6631.218, observed depth 2.2%)

The fits claim more absorbed light than the spectrum is missing. That is the
+0.98 dex, and no amount of NLTE, gf or 3D work would ever have touched it.

What this plots
---------------
For each element, every CANONICAL line that yields no measured value, with the
four facts that decide its fate, measured here and not read from any label:

  * is it inside our wavelength coverage at all;
  * its observed central depth;
  * the integrated absorption actually present in its window;
  * what the fitter claimed, and what curation then did with it.

`--element` restricts the sweep; the default is every element the tracker
carries at an owed / gf_floor / upper_limit tier. Figures land in
`results/plots/unmeasured/` and are copied to the artifact store, because a
gitignored proof that lives in one worktree is not archived (RYA-461).

Integrated absorption is the honest yardstick and is deliberately naive: the
trapezoid of (1 - normalised flux) over the window, continuum taken as 1.0. It
does not deblend and it does not model. It cannot tell you a line's abundance —
what it CAN do is put a ceiling on any EW claimed inside it, which is the only
question being asked here.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from config import constants as const           # noqa: E402

OUT_DIR = ROOT / "results" / "plots" / "unmeasured"

#: A linelist entry is a CANONICAL LINE candidate if its predicted core is at least
#: this deep. Below ~1% the feature cannot carry an abundance at solar S/N, so
#: plotting it would bury the real cases in noise. This is a plotting cutoff, never a
#: science one — nothing downstream reads it.
DEPTH_FLOOR = 0.010
#: Linelist entries closer than this are one feature (HFS/isotopic components, and
#: blends the fitter cannot separate). Grouped so a 6-component HFS line plots once.
GROUP_A = 0.05
#: Half-width of the window over which absorption is integrated and depth measured.
WINDOW_A = 0.40
#: Half-width over which absorption is integrated for the OVER-SUBSCRIPTION screen, and
#: the (smaller) half-width inside which a fit is counted against it.
#:
#: These differ ON PURPOSE, and the first version of this screen got it wrong. Summing
#: every fit centred within +/-0.4 A and comparing to the absorption in that same +/-0.4 A
#: fires on any strong line, because a line centred near the window edge has much of its
#: EW OUTSIDE the window -- the fit is honest and the comparison is not. That first pass
#: flagged 273 in-pool lines, nearly all of them false. Integrating over a window WIDER
#: than the one fits are drawn from keeps the wings inside the budget.
ABSORB_A = 1.00
FIT_A = 0.40
#: Ratio above which a window is called over-subscribed. Not a physical constant: a
#: screening threshold, deliberately blunt, and it is why this verdict names a window to
#: go and look at rather than asserting a line is wrong.
OVERSUB_RATIO = 1.5

# ── WHAT OVER-SUBSCRIBED IS AND IS NOT ───────────────────────────────────────
# It is a TRIAGE FLAG. It is not a verdict, and it must not be counted up and reported
# as "N bad lines" — I did exactly that on the first pass and it was wrong.
#
# It compares the sum of fitted EWs in a window against the light missing from a wider
# one. In an ISOLATED window that is close to airtight: Al I 6696.020 and 6696.185 sit
# 0.165 A apart on a flat continuum, and 119.95 mA of fitted EW against 66.48 mA of
# missing light is a real double-count, verified by hand.
#
# In a CROWDED window it over-fires, for a reason that is physics and not a bug: a
# strong line's damping wings carry EW far outside any window you draw, so a fit centred
# inside the window legitimately claims absorption from outside it. Ca I 6122 (78% deep)
# and the Na D-region lines trip it while being perfectly well measured.
#
# So: the flag says LOOK HERE. A window it flags is confirmed only by checking that the
# neighbourhood is clean and the fitted centres are well inside the integration span.
# Anything else is the over-claiming this script was written to catch.
MAX_PANELS = 8

INK, MUTE = "#1a1a1a", "#8a8a8a"
C_OK, C_BAD, C_SAT, C_NODATA = "#1a7f37", "#b3261e", "#8250df", "#8a8a8a"


def _spectrum() -> tuple[np.ndarray, np.ndarray]:
    p = Path(str(const.PATHS["solar_normalized"]))
    if not p.exists():
        raise SystemExit(
            f"normalised solar spectrum missing at {p} — this script proves things "
            f"against the spectrum and must not run without it (RYA-518).")
    d = pd.read_csv(p)
    return d.wavelength_air_A.values, d.flux_normalized.values


def _integrated_absorption(W, F, lo, hi) -> float:
    """mA of light actually missing between lo and hi. Continuum = 1.0."""
    m = (W > lo) & (W < hi)
    if m.sum() < 3:
        return float("nan")
    from pipeline._numcompat import trapezoid          # RYA numpy-2 shim
    return float(trapezoid(1.0 - F[m], W[m]) * 1000.0)


def _depth(W, F, w0, half=0.06) -> float:
    m = (W > w0 - half) & (W < w0 + half)
    return float(1.0 - F[m].min()) if m.sum() else float("nan")


def canonical_lines(ll: pd.DataFrame, element: str) -> pd.DataFrame:
    """Group the linelist into distinct features for one element."""
    a = ll[(ll.element == element)].copy()
    if a.empty:
        return a
    # Cluster on the GAP between consecutive wavelengths, never on a rounded key.
    # `round(w / 0.05)` splits two entries 0.019 A apart across a bin boundary -- which
    # is how Al I 6698.662 and 6698.681 first plotted as two separate "lines". Same
    # defect as RYA-703 (a rounded key that misses a match) and RYA-704 (one that merges
    # two things), filed the same day. A rounded wavelength is a tolerance, not an identity.
    a = a.sort_values("wavelength_air_A")
    gap = a.wavelength_air_A.diff().fillna(9e9)
    a["_key"] = (gap > GROUP_A).cumsum()
    g = a.groupby("_key").agg(
        wavelength_air_A=("wavelength_air_A", "mean"),
        ion=("ion", lambda s: s.mode().iat[0]),
        log_gf=("log_gf", "max"),
        excitation_potential_eV=("excitation_potential_eV", "min"),
        predicted_depth=("central_depth", "max"),
        n_components=("wavelength_air_A", "size")).reset_index(drop=True)
    return g[g.predicted_depth >= DEPTH_FLOOR].sort_values(
        "predicted_depth", ascending=False)


def synthesis_required() -> dict[str, str]:
    """Elements whose value comes from SYNTHESIS, not from an equivalent width.

    Ryan, 2026-08-08: *"remember we have to use synth on some of these as well"* and
    *"nothing should be suppressed."* Both are right, and the second is a correction to
    how this was first written.

    Nothing is suppressed. The EW is still measured for every line and still shown here.
    What RYA-520 changes for these elements is which channel REPORTS the value: a
    blended, saturated or HFS-split feature cannot be inverted from an equivalent width,
    so the synthesis fit supersedes the EW as the reported number. The EW remains a
    diagnostic, which is exactly why the ladder below computes it either way.

    Membership is read from the problem-children registry (`required_treatment`) plus the
    CNO set, exactly as `rya527_two_engine_run` computes it — never a hardcoded list here,
    so an element changing treatment needs no edit in this file.
    """
    import pipeline.problem_children as pc
    out = {el: "CNO — multi-arm synthesis (RYA-491/237)" for el in ("C", "N", "O")}
    for el in const.TARGET_ELEMENTS:
        d = pc.disposition_for(el) or {}
        t = d.get("required_treatment")
        if t in ("synthesis", "HFS_sum"):
            out[el] = f"{t} (problem-children registry)"
    return out


def fate(row, W, F, fitted: pd.DataFrame, pool_waves: set,
         synth_note: str | None = None) -> dict:
    """The four measured facts, plus what the pipeline did."""
    w0 = float(row.wavelength_air_A)
    in_cov = bool(W.min() <= w0 <= W.max())
    out = {"wave": w0, "ion": row.ion, "in_coverage": in_cov,
           "predicted_depth": float(row.predicted_depth),
           "n_components": int(row.n_components)}
    if not in_cov:
        out.update(depth=float("nan"), absorbed=float("nan"), fitted_ew=None,
                   in_pool=False, verdict="NO DATA", colour=C_NODATA,
                   why="outside our wavelength coverage")
        return out
    out["depth"] = _depth(W, F, w0)
    out["absorbed"] = _integrated_absorption(W, F, w0 - WINDOW_A, w0 + WINDOW_A)
    near = fitted[(fitted.wavelength_air_A - w0).abs() < GROUP_A]
    out["fitted_ew"] = float(near.ew_mA.sum()) if len(near) else None
    out["n_fits"] = int(len(near))
    out["in_pool"] = any(abs(pw - w0) < GROUP_A for pw in pool_waves)
    # OVER-SUBSCRIPTION: every fit whose centre falls in the window, summed, against the
    # light actually missing from it. A single fit can look reasonable while two fits of
    # one blended feature each claim the whole thing -- Al I 6696.020 and 6696.185 are
    # 0.165 A apart and were fitted independently at 60.30 and 59.65 mA, 119.95 mA in
    # total, from a window containing 66.48 mA. Comparing one fit at a time cannot see it.
    win = fitted[(fitted.wavelength_air_A > w0 - FIT_A)
                 & (fitted.wavelength_air_A < w0 + FIT_A)]
    out["window_fitted"] = float(win.ew_mA.sum()) if len(win) else 0.0
    out["window_n"] = int(len(win))
    out["window_absorbed"] = _integrated_absorption(W, F, w0 - ABSORB_A, w0 + ABSORB_A)

    # THE LADDER (Ryan, 2026-08-08): "check EW / 1D / LTE, if something wonky /
    # special / etc, then move onto our other models."
    #
    # So the EW verdict is computed for EVERY line, including synthesis-required ones.
    # The earlier version returned early on `synth_note` and that was wrong in a way
    # that matters before a freeze: it ASSUMED the escalation was warranted instead of
    # showing it. Escalation is a consequence of the EW check, not a substitute for it,
    # and the two cases worth catching are only visible when both are computed:
    #
    #   escalated, EW clean   -> why is this element on synthesis? possibly unnecessary
    #   not escalated, EW wonky -> it should have been escalated and was not
    #
    # `escalation` therefore rides ALONGSIDE the verdict and never replaces it.
    out["escalation"] = synth_note

    if out["depth"] >= 0.90:
        out.update(verdict="SATURATED", colour=C_SAT,
                   why=f"core is {out['depth']*100:.0f}% deep — an EW carries no "
                       f"abundance information here")
    elif out["fitted_ew"] is None:
        out.update(verdict="NEVER FITTED", colour=C_BAD,
                   why="in coverage, and no fit was attempted")
    elif (out["window_n"] > 1
          and out["window_fitted"] > out["window_absorbed"] * OVERSUB_RATIO):
        out.update(verdict="OVER-SUBSCRIBED", colour=C_BAD,
                   why=f"{out['window_n']} fits claim {out['window_fitted']:.1f} mA "
                       f"between them; only {out['window_absorbed']:.1f} mA is missing "
                       f"across +/-{ABSORB_A} A "
                       f"({out['window_fitted']/out['window_absorbed']:.1f}x) — "
                       f"look at this window")
    elif out["fitted_ew"] > out["window_absorbed"] * OVERSUB_RATIO:
        out.update(verdict="OVER-MEASURED", colour=C_BAD,
                   why=f"one fit claims {out['fitted_ew']:.1f} mA where "
                       f"{out['window_absorbed']:.1f} mA is missing across "
                       f"+/-{ABSORB_A} A ({out['fitted_ew']/out['window_absorbed']:.1f}x)")
    elif not out["in_pool"]:
        out.update(verdict="FITTED, CULLED", colour=C_BAD,
                   why=f"fitted {out['fitted_ew']:.1f} mA and dropped before the "
                       f"committed pool")
    else:
        out.update(verdict="IN POOL", colour=C_OK,
                   why=f"fitted {out['fitted_ew']:.1f} mA, kept")

    ew_wonky = out["verdict"] in ("SATURATED", "OVER-SUBSCRIBED", "OVER-MEASURED")
    if out["escalation"] and ew_wonky:
        out["ladder"] = "ESCALATED, warranted"          # the system working
    elif out["escalation"]:
        out["ladder"] = "ESCALATED, EW looked usable"   # check the escalation is needed
    elif ew_wonky:
        out["ladder"] = "NOT ESCALATED, EW is wonky"    # the gap that matters pre-freeze
        out["colour"] = C_BAD
    else:
        out["ladder"] = "EW route, clean"
    return out


def plot_element(element: str, W, F, ll, fitted, pool_waves, out_dir: Path,
                 synth_note: str | None = None) -> Path | None:
    lines = canonical_lines(ll, element)
    if lines.empty:
        print(f"  {element}: no canonical line above the {DEPTH_FLOOR:.0%} depth floor")
        return None
    facts = [fate(r, W, F, fitted, pool_waves, synth_note)
             for _, r in lines.iterrows()]
    # "IN POOL" is the only verdict that means the line is fine. OVER-MEASURED and
    # OVER-SUBSCRIBED lines may well BE in the pool -- that is worse than absence,
    # not better, so they stay in the figure.
    unmeasured = [f for f in facts
                  if f["verdict"] != "IN POOL" or f.get("ladder", "").startswith("NOT ESC")]
    # Rank IN-COVERAGE first. Sorting on predicted depth alone fills the figure with
    # far-UV lines we have never observed and never will from the ground — they are
    # the deepest entries in a VALD pull and the least informative panels here.
    covered = sorted([f for f in unmeasured if f["in_coverage"]],
                     key=lambda f: -f["predicted_depth"])
    # A couple of out-of-coverage panels still earn their place: for Al they ARE the
    # answer (7835/8772 are the best lines in the element). Keep the strongest few,
    # inside the visible range first, and never let them crowd out a real spectrum.
    missing = sorted([f for f in unmeasured if not f["in_coverage"]],
                     key=lambda f: (abs(f["wave"] - 0.5 * (W.min() + W.max())),))[:2]
    show = (covered + missing)[:MAX_PANELS]
    if not show:
        print(f"  {element}: every canonical line is in the pool — nothing to prove")
        return None

    n = len(show)
    ncol = 2 if n > 1 else 1
    nrow = (n + ncol - 1) // ncol
    fig, axes = plt.subplots(nrow, ncol, figsize=(6.0 * ncol, 2.5 * nrow), squeeze=False)
    for ax, f in zip(axes.ravel(), show):
        w0 = f["wave"]
        if not f["in_coverage"]:
            ax.text(0.5, 0.5, f"{element} {f['ion']}  {w0:.3f} Å\n\nNO DATA\n"
                              f"outside {W.min():.0f}–{W.max():.0f} Å",
                    ha="center", va="center", fontsize=10, color=C_NODATA,
                    fontweight="bold", transform=ax.transAxes)
            ax.set_xticks([]); ax.set_yticks([])
            for s in ax.spines.values():
                s.set_color(C_NODATA); s.set_linestyle(":")
            continue
        m = (W > w0 - 1.2) & (W < w0 + 1.2)
        ax.plot(W[m], F[m], color=INK, lw=0.85)
        ax.axhline(1.0, color=MUTE, lw=0.5, ls=":")
        ax.axvspan(w0 - WINDOW_A, w0 + WINDOW_A, color=f["colour"], alpha=0.07)
        ax.axvline(w0, color=f["colour"], lw=1.0, ls="--")
        lo = min(F[m].min() - 0.04, 0.96)
        ax.set_ylim(lo, 1.06); ax.set_xlim(w0 - 1.2, w0 + 1.2)
        hdr = f"{element} {f['ion']} {w0:.3f}  —  {f['verdict']}"
        if f.get("ladder") and f["ladder"] != "EW route, clean":
            hdr += f"   [{f['ladder']}]"
        ax.set_title(hdr, fontsize=8.5, color=f["colour"], fontweight="bold", loc="left")
        bits = [f"observed depth {f['depth']*100:.1f}%",
                f"absorbed {f['absorbed']:.1f} mÅ in ±{WINDOW_A} Å"]
        if f["fitted_ew"] is not None:
            bits.append(f"fitted {f['fitted_ew']:.1f} mÅ"
                        + (f" ({f['n_fits']} fits)" if f["n_fits"] > 1 else ""))
        ax.text(0.02, 0.06, "\n".join(bits) + f"\n{f['why']}", transform=ax.transAxes,
                fontsize=6.9, color=MUTE, va="bottom",
                bbox=dict(fc="white", ec="none", alpha=0.75, pad=1.5))
        ax.tick_params(labelsize=7)
    for ax in axes.ravel()[n:]:
        ax.axis("off")
    head = f"{element} — canonical lines that produce no measurement"
    if synth_note:
        head += "   [SYNTHESIS element: EW measured, synthesis reports]"
    fig.suptitle(head,
                 fontsize=12, fontweight="bold", color=INK, x=0.01, ha="left")
    fig.text(0.01, 0.005,
             "Absorption is the trapezoid of (1 − normalised flux); it does not deblend. "
             "It cannot give an abundance — it bounds any EW claimed inside the window.",
             fontsize=7.5, color=MUTE, ha="left")
    fig.tight_layout(rect=[0, 0.015, 1, 0.97])
    out_dir.mkdir(parents=True, exist_ok=True)
    p = out_dir / f"{element}_unmeasured.png"
    fig.savefig(p, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  {element}: {len(show)} line(s) plotted → {p.name}")
    return p


def owed_elements() -> list[str]:
    from pipeline import state_surfaces
    df = pd.read_csv(ROOT / state_surfaces.TRACKER, comment="#")
    keep = {"owed", "nlte_owed", "curation_owed", "gf_floor", "upper_limit"}
    return [str(r.element) for _, r in df.iterrows()
            if str(r.get("tier", "")).strip().lower() in keep]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--element", action="append",
                    help="restrict to this element (repeatable); default = every "
                         "owed / gf_floor / upper_limit element in the tracker")
    ap.add_argument("--out", default=str(OUT_DIR))
    ap.add_argument("--no-archive", action="store_true",
                    help="skip the artifact-store copy (RYA-461)")
    args = ap.parse_args(argv)

    W, F = _spectrum()
    ll = pd.read_csv(str(const.PATHS["linelist_solar"]), low_memory=False)
    pool = pd.read_csv(str(const.PATHS["solar_ew_canonical"]), comment="#")
    pool_waves = set(pool.wavelength_air_A.astype(float))

    staging = Path(str(const.PATHS["solar_ew"]))
    if staging.exists():
        fitted = pd.read_csv(staging, low_memory=False)
    else:
        fitted = pool[["element", "ion", "wavelength_air_A", "ew_mA"]].copy()
        print(f"note: staging fit output absent; 'NEVER FITTED' below means 'not in the "
              f"committed pool'. Regenerate {staging.name} for the full fit history.")

    synth = synthesis_required()
    els = args.element or owed_elements()
    ssel = sorted(set(els) & set(synth))
    print(f"synthesis-reported (EW still measured, RYA-520): {ssel or 'none'}")
    print(f"spectrum {W.min():.1f}–{W.max():.1f} Å; {len(els)} element(s)")
    out_dir, made = Path(args.out), []
    for el in els:
        p = plot_element(el, W, F, ll, fitted, pool_waves, out_dir, synth.get(el))
        if p:
            made.append(p)

    if made and not args.no_archive:
        from pipeline.artifact_store import save_artifact
        for p in made:
            save_artifact(str(p.relative_to(ROOT)), kind="plots")
        print(f"archived {len(made)} figure(s) to the artifact store")
    print(f"\n{len(made)} figure(s) in {out_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
