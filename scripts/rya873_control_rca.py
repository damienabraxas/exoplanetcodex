#!/usr/bin/env python3
"""
RYA-873 — is ANGLE 1 comparing like with like? The same-species neighbour test.
==============================================================================
    ISPEC_DIR=... python3 scripts/rya873_control_rca.py [--control <per-line csv>]

The SynthesisHandler's optical control PASSES on abundance (dex_offset −0.0100) and FAILS
on EW (median synth/banked ratio 1.433, +0.1562 dex), and reports itself as
*"compensating errors. Treat with MORE suspicion than a clean failure, because it looks
like success."* RYA-873 has to decide what harness residual that earns, and the answer
turns entirely on whether the EW disagreement is real.

WHY IT MIGHT NOT BE
-------------------
The ratios are not an offset. They span **0.50 to 7.60** with MAD 0.359 about a median of
1.433, while the abundances from those same fits scatter by only 0.065 dex. A handler that
mis-scaled EW would be wrong by A factor, not by a factor that varies 15x line to line.

`SynthesisHandler.measure_line` builds its EW as a DIFFERENCE OF TWO SYNTHESES — the
window at the fitted abundance, minus the window with the element depleted by 4 dex:

    ew = trapz(clip(without - with_el, 0, None), ew_w) * 1000        # ew_hw = max(3*hw_A, 1.0)

🔴 THAT CANCELS OTHER SPECIES AND DOES NOT CANCEL THIS ONE. The code says "blends cancel",
and for Ni, Ca, CH it is exactly right: they absorb identically in both syntheses, so they
subtract away. But ANOTHER Fe I LINE also vanishes when Fe is depleted 4 dex, so its
absorption lands in the difference and is counted into the target line's EW. The window is
at least +/-1.0 A wide, and in the solar optical Fe I is dense.

Meanwhile `ew_ref` is a banked PROFILE FIT to the target line alone. If the mechanism is
real, the two are simply not the same quantity for a crowded line — and the control's
"compensating errors" verdict would be an artifact of the comparison rather than evidence
about the handler.

This is the RYA-785 finding in a second place (`ew()` over +/-1.2 A: a same-species
neighbour does NOT cancel) and the RYA-679 shape (full-window red_chi2 measures the BLEND
LIST, not the element: the same Zr II lines scored 83.12 vs 0.39 deblended).

THE FALSIFIABLE PREDICTION
--------------------------
If the mechanism is the explanation, then:

1. lines with NO Fe I neighbour inside the EW window sit at ratio ~1;
2. for the rest, the ratio tracks how much of the window's Fe I absorption belongs to
   somebody else -- predicted ratio ~ (sum of Fe I strength in window) / (target strength).

If instead the isolated lines are ALSO scattered, the mechanism is refuted, the EW
disagreement is real, and the compensating-error caveat stands. Both outcomes are useful
and the script states which one it found rather than arguing for either.

⚠️ THE STRENGTH PROXY IS A PROXY. `10^(loggf - ep*theta)` is the curve-of-growth strength
ordering, not an EW: it is monotonic in line strength but not linear in it, and it ignores
saturation, damping and the fact that overlapping lines do not add linearly. So the
CORRELATION is the evidence and the predicted ratio is an ordering, not a fit. The
prediction that carries real weight is (1), which needs no proxy at all.
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
sys.path.insert(0, str(ROOT / "scripts"))

OUT = ROOT / "data" / "results" / "rya873"
CONTROL_CSV = ROOT / "data" / "audit" / "synthesis_control" / "synthesis_control_FeI.csv"

#: EW-window half-widths to score, in A. The handler uses `max(3 * hw_A, 1.0)`, where
#: `hw_A` is wing-wide from the line's own EW hint and is NOT recorded in the control
#: artifact. Rather than reconstruct it (and inherit any drift in that rule), the test is
#: run at a ladder of fixed windows: if the mechanism is real the correlation should
#: STRENGTHEN as the window grows toward the handler's, and that trend is itself evidence.
WINDOWS_A = (0.25, 0.5, 1.0, 1.5, 2.0)

#: Boltzmann factor for the solar photosphere, 5040/Teff at Teff = 5772 K (stars.yaml).
THETA_SOLAR = 5040.0 / 5772.0


def fe1_lines(linelist) -> pd.DataFrame:
    """Fe I rows of the loaded synthesis list, with a curve-of-growth strength ordering."""
    from pipeline import gf_rung
    d = gf_rung.linelist_frame(linelist)
    d = d[d.species == "Fe 1"].copy()
    d["strength"] = 10.0 ** (d.log_gf - d.ep_eV * THETA_SOLAR)
    return d.reset_index(drop=True)


def crowding(control: pd.DataFrame, fe1: pd.DataFrame, half_A: float) -> pd.DataFrame:
    """Per control line: how much of the window's Fe I absorption is NOT the target."""
    lw = fe1.wavelength_air_A.to_numpy()
    ls = fe1.strength.to_numpy()
    rows = []
    for c in control.wave.astype(float):
        m = np.abs(lw - c) <= half_A
        if not m.any():
            rows.append({"n_fe1_in_window": 0, "n_neighbours": 0,
                         "target_strength": np.nan, "window_strength": 0.0,
                         "predicted_ratio": np.nan})
            continue
        idx = np.flatnonzero(m)
        # The TARGET is the nearest Fe I row; every other Fe I row in the window is a
        # neighbour whose absorption the difference-of-syntheses cannot remove.
        t = idx[np.argmin(np.abs(lw[idx] - c))]
        tot = float(ls[idx].sum())
        tgt = float(ls[t])
        rows.append({"n_fe1_in_window": int(m.sum()),
                     "n_neighbours": int(m.sum()) - 1,
                     "target_strength": tgt, "window_strength": tot,
                     "predicted_ratio": (tot / tgt) if tgt > 0 else np.nan})
    return pd.DataFrame(rows)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--control", type=Path, default=CONTROL_CSV)
    ap.add_argument("--out", type=Path, default=OUT)
    a = ap.parse_args()

    d = pd.read_csv(a.control)
    ok = d[(d.status == "ok") & d.ew_synth.notna() & d.ew_ref.notna()].copy()
    ok["ratio"] = ok.ew_synth / ok.ew_ref
    print(f"=== the control's accepted lines ===")
    print(f"  {len(ok)} of {len(d)} accepted; ratio median {ok.ratio.median():.3f}, "
          f"range {ok.ratio.min():.2f}-{ok.ratio.max():.2f}")

    # RYA-847 part 2 is merged and the NON-MINIMUM check is always on. Report whether it
    # touched this control at all, so the numbers below are anchored to a stated line set.
    if "frac_rise_weaker" in d.columns:
        fr = pd.to_numeric(d.frac_rise_weaker, errors="coerce")
        n_nonmin = int((fr <= 0).sum())
        print(f"  RYA-847 non-minimum check: {n_nonmin} of {len(d)} lines have "
              f"frac_rise <= 0" + ("" if n_nonmin else "  (the gate does not bite here)"))

    from pipeline.abundances_derive import _load_synth_resources
    ll, _, _ = _load_synth_resources()
    fe1 = fe1_lines(ll)
    print(f"  synthesis list: {len(fe1)} Fe I rows")

    results, tables = [], {}
    for half in WINDOWS_A:
        cr = crowding(ok.reset_index(drop=True), fe1, half)
        t = pd.concat([ok.reset_index(drop=True), cr], axis=1)
        t["half_A"] = half
        isolated = t[t.n_neighbours == 0]
        crowded = t[t.n_neighbours > 0]
        v = t[t.predicted_ratio.notna() & np.isfinite(t.predicted_ratio)]
        rho = (float(np.corrcoef(np.log10(v.predicted_ratio), np.log10(v.ratio))[0, 1])
               if len(v) > 2 else np.nan)
        results.append({
            "half_A": half, "n": len(t),
            "n_isolated": len(isolated), "n_crowded": len(crowded),
            "isolated_ratio_median": (round(float(isolated.ratio.median()), 3)
                                      if len(isolated) else None),
            "crowded_ratio_median": (round(float(crowded.ratio.median()), 3)
                                     if len(crowded) else None),
            "corr_log_predicted_vs_log_observed": (None if not np.isfinite(rho)
                                                   else round(rho, 3)),
        })
        tables[half] = t

    r = pd.DataFrame(results)
    print(f"\n=== the prediction, at a ladder of EW windows ===")
    print(f"  (1) isolated lines should sit at ratio ~1   "
          f"(2) log-log correlation should STRENGTHEN with window")
    print(f"{'half_A':>8}{'n_iso':>7}{'iso median':>12}{'n_crowd':>9}"
          f"{'crowd median':>14}{'corr':>8}")
    for _, x in r.iterrows():
        print(f"{x.half_A:>8.2f}{x.n_isolated:>7}"
              f"{(x.isolated_ratio_median if x.isolated_ratio_median is not None else float('nan')):>12}"
              f"{x.n_crowded:>9}"
              f"{(x.crowded_ratio_median if x.crowded_ratio_median is not None else float('nan')):>14}"
              f"{(x.corr_log_predicted_vs_log_observed if x.corr_log_predicted_vs_log_observed is not None else float('nan')):>8}")

    # ── the verdict, stated as the thing that would refute it ────────────────────
    best = max((x for x in results
                if x["corr_log_predicted_vs_log_observed"] is not None),
               key=lambda x: x["corr_log_predicted_vs_log_observed"], default=None)
    print(f"\n=== verdict ===")
    if best is None:
        print("  the correlation is undefined — too few lines carry a prediction")
    else:
        print(f"  strongest at half_A = {best['half_A']} A: "
              f"corr(log predicted, log observed) = "
              f"{best['corr_log_predicted_vs_log_observed']}")
        iso = best["isolated_ratio_median"]
        print(f"  isolated-line median ratio there: {iso} "
              f"({best['n_isolated']} lines)")
        print("  A ratio near 1 on isolated lines AND a positive correlation on the rest")
        print("  supports: ew_synth counts SAME-SPECIES neighbours the difference cannot")
        print("  cancel, so ANGLE 1 is not comparing like with like and its +0.1562 dex")
        print("  is an artifact of the comparison, not a property of the handler.")
        print("  Scattered isolated lines, or no correlation, REFUTES that and leaves the")
        print("  compensating-error caveat standing.")

    a.out.mkdir(parents=True, exist_ok=True)
    pd.concat(tables.values(), ignore_index=True).to_csv(
        a.out / "rya873_crowding_by_line.csv", index=False)
    (a.out / "rya873_crowding_summary.json").write_text(json.dumps({
        "ticket": "RYA-873",
        "control": str(a.control),
        "n_accepted": int(len(ok)),
        "ratio_median": round(float(ok.ratio.median()), 4),
        "ratio_min": round(float(ok.ratio.min()), 4),
        "ratio_max": round(float(ok.ratio.max()), 4),
        "theta_solar": round(THETA_SOLAR, 5),
        "windows": results,
    }, indent=2) + "\n")
    print(f"\n  wrote {a.out.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
