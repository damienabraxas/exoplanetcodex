#!/usr/bin/env python3
"""
scripts/stage2_pool_ledger_rya706.py — RYA-706
==============================================
NO SILENT DROP, STAGE 2: FIT OUTPUT -> COMMITTED POOL.

RYA-429 built the no-silent-drop ledger and it covers stage 1 only.

  stage 1  line list -> fit    23371 -> 8950   14421 rejections ALL classified,
                                               n_unaccounted 0, invariant_holds true
  stage 2  fit -> committed pool 8951 ->  808   8102 dropped, NO recorded reason

`scripts/promote_solar_ew.py` says why in its own docstring: "lines new in staging are
reported but not auto-added (adding a line to the canonical is a curation decision)."
So the 808-line pool was hand-picked once, at the RYA-396 vendoring in June 2026, and
the reasoning was never written down.

That is not a bookkeeping complaint. Solar Al I 6698.671 is clean -- everything non-Al
within +/-0.35 A is 0.1-0.3% deep against a 16.7% feature, the four Al entries share one
excitation potential (hyperfine components of ONE line), and it fits at chi2 0.0001. It
is not in the pool, and nothing says why. Al then reports no value at all, and that
absence was explained twice, wrongly.

WHAT THIS DOES
--------------
Re-derives stage 2 against the pipeline's OWN thresholds -- never invented ones. Each
dropped line gets the FIRST reason that applies, in the order the pipeline would have
applied them, and anything no rule explains lands in `UNEXPLAINED`, which is the whole
point of the exercise.

Ryan ratified 2026-08-08: "run everything on every model, just to be sure" and "every
line you think is not working, we double check, we prove, and use all models to prove
it." A line's fate is never decided by an unrecorded selection step.

WHAT IT DOES NOT DO
-------------------
It does not add a line to the pool. `promote_solar_ew.py` is right that adding one is a
reviewed curation act; this makes the gap visible and reproducible so that review has
something to read. It changes no value anywhere.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from config import constants as const                      # noqa: E402
from config.constants import PIPELINE                      # noqa: E402
import pipeline.curate_nonfe_pools as cn                   # noqa: E402

OUT = ROOT / "data" / "audit" / "stage2_pool_ledger"
#: The June fit output, preserved to the artifact store before the stale copy was
#: cleared (RYA-461 save-before-clean). Without it stage 2 is unreconstructable.
STAGING_FALLBACK = Path("/Users/ryanschmitt/Documents/Exoplanet Codex/data/solar_ew.csv")

#: Contaminant tolerance for the "is this line clean" test. A neighbour shallower than
#: this fraction of the target's own depth cannot meaningfully corrupt its EW. 0.10 is a
#: STATED CHOICE, not a pipeline constant -- there is no existing rule for it, and
#: inventing one silently is the failure this ticket is about.
CLEAN_CONTAM_FRAC = 0.10
CLEAN_WINDOW_A = 0.35


def _spectrum():
    d = pd.read_csv(str(const.PATHS["solar_normalized"]))
    return d.wavelength_air_A.values, d.flux_normalized.values


def _depth(W, F, w0, half=0.06):
    m = (W > w0 - half) & (W < w0 + half)
    return float(1.0 - F[m].min()) if m.sum() else float("nan")


def _rew(ew_mA, lam):
    if not (ew_mA and ew_mA > 0 and lam):
        return float("nan")
    return math.log10(ew_mA * 1e-3 / lam)


def classify(row, W, F, ll, synth: dict) -> tuple[str, str]:
    """FIRST applicable reason, in pipeline order. Returns (code, detail)."""
    el, lam, ew = str(row.element), float(row.wavelength_air_A), float(row.ew_mA)
    err = float(row.ew_err_mA) if pd.notna(row.ew_err_mA) else float("nan")

    if el in synth:
        return "SYNTH_ROUTED", f"{el}: synthesis is the reporting channel ({synth[el]})"
    if not (W.min() <= lam <= W.max()):
        return "OUT_OF_COVERAGE", f"{lam:.3f} outside {W.min():.0f}-{W.max():.0f} A"
    if ew < PIPELINE["ew_min_mA"]:
        return "EW_BELOW_MIN", f"{ew:.2f} < {PIPELINE['ew_min_mA']} mA"
    if ew > PIPELINE["ew_max_mA"]:
        return "EW_ABOVE_MAX", f"{ew:.2f} > {PIPELINE['ew_max_mA']} mA"
    rew = _rew(ew, lam)
    if rew > cn.REW_LINEAR_CEILING:
        return "SATURATED_REW", f"REW {rew:.3f} > {cn.REW_LINEAR_CEILING}"
    d = _depth(W, F, lam)
    if d >= 0.90:
        return "SATURATED_CORE", f"observed depth {d*100:.1f}%"
    if np.isfinite(err) and ew > 0 and err / ew > cn.EW_ERR_FRAC_MAX:
        return "EW_ERR_EXCESS", f"err/EW {err/ew:.2f} > {cn.EW_ERR_FRAC_MAX}"
    return "UNEXPLAINED", "no pipeline rule accounts for this drop"


def contaminants(ll, el, lam, own_depth):
    """Neighbours of a DIFFERENT species within the clean window, deeper than the
    stated fraction of the line's own observed depth."""
    w = ll[(ll.wavelength_air_A > lam - CLEAN_WINDOW_A)
           & (ll.wavelength_air_A < lam + CLEAN_WINDOW_A)
           & (ll.element != el)]
    if w.empty or not np.isfinite(own_depth) or own_depth <= 0:
        return 0, 0.0
    worst = float(w.central_depth.max())
    n_bad = int((w.central_depth > own_depth * CLEAN_CONTAM_FRAC).sum())
    return n_bad, worst


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--staging", default=None)
    ap.add_argument("--out", default=str(OUT))
    args = ap.parse_args(argv)

    staging = Path(args.staging) if args.staging else Path(str(const.PATHS["solar_ew"]))
    if not staging.exists():
        staging = STAGING_FALLBACK
    if not staging.exists():
        raise SystemExit(
            f"stage-2 reconstruction needs the fit output and it is absent. Looked at "
            f"{const.PATHS['solar_ew']} and {STAGING_FALLBACK}. It is gitignored and "
            f"regenerable; the June copy lives in the artifact store (RYA-461).")

    W, F = _spectrum()
    fit = pd.read_csv(staging, low_memory=False)
    pool = pd.read_csv(str(const.PATHS["solar_ew_canonical"]), comment="#")
    ll = pd.read_csv(str(const.PATHS["linelist_solar"]), low_memory=False)

    import importlib.util as iu
    sp = iu.spec_from_file_location("g", ROOT / "scripts" / "plot_unmeasured_lines_rya707.py")
    gen = iu.module_from_spec(sp); sp.loader.exec_module(gen)
    synth = gen.synthesis_required()

    key = lambda d: set(zip(d.element.astype(str), d.ion.astype(str),
                            d.wavelength_air_A.astype(float).round(2)))
    in_pool = key(pool)

    rows = []
    for _, r in fit.iterrows():
        k = (str(r.element), str(r.ion), round(float(r.wavelength_air_A), 2))
        if k in in_pool:
            continue
        code, detail = classify(r, W, F, ll, synth)
        d = _depth(W, F, float(r.wavelength_air_A))
        n_bad, worst = contaminants(ll, str(r.element), float(r.wavelength_air_A), d)
        rows.append(dict(element=r.element, ion=r.ion,
                         wavelength_air_A=round(float(r.wavelength_air_A), 4),
                         ew_mA=r.ew_mA, ew_err_mA=r.ew_err_mA, chi2=r.get("chi2"),
                         observed_depth=round(d, 4) if np.isfinite(d) else None,
                         rew=round(_rew(float(r.ew_mA), float(r.wavelength_air_A)), 4),
                         n_contaminants=n_bad, worst_contaminant_depth=round(worst, 4),
                         drop_reason=code, detail=detail))
    led = pd.DataFrame(rows)

    n_unexplained = int((led.drop_reason == "UNEXPLAINED").sum())
    # RECOVERY CANDIDATE: unexplained, clean, unsaturated, well-fitted, in coverage.
    # Deliberately conservative -- this list is meant to be short and defensible, not
    # long. Every criterion is one an unrecovered line would visibly fail.
    rec = led[(led.drop_reason == "UNEXPLAINED")
              & (led.n_contaminants == 0)
              & (led.observed_depth.astype(float) >= 0.02)
              & (led.chi2.astype(float) <= 0.01)].copy()
    rec = rec.sort_values(["element", "wavelength_air_A"])

    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    led.to_csv(out / "stage2_drops.csv", index=False)
    rec.to_csv(out / "recovery_candidates.csv", index=False)
    summary = dict(
        ticket="RYA-706", stage="fit -> committed pool",
        staging=str(staging), staging_rows=int(len(fit)),
        pool_rows=int(len(pool)), dropped=int(len(led)),
        n_unexplained=n_unexplained,
        n_recovery_candidates=int(len(rec)),
        reason_counts=dict(Counter(led.drop_reason)),
        unexplained_by_element=dict(Counter(
            led.loc[led.drop_reason == "UNEXPLAINED", "element"]).most_common()),
        recovery_by_element=dict(Counter(rec.element).most_common()),
        invariant_note=("every dropped line carries a classified reason; UNEXPLAINED is "
                        "itself a classification and counts as an unrecorded drop, not "
                        "as accounted-for"))
    (out / "summary.json").write_text(json.dumps(summary, indent=2, default=str))

    print(f"stage 2: {len(fit)} fitted -> {len(pool)} pooled; {len(led)} dropped")
    for k, v in sorted(summary["reason_counts"].items(), key=lambda kv: -kv[1]):
        print(f"  {k:18s} {v}")
    print(f"\nUNEXPLAINED: {n_unexplained}")
    for e, n in list(summary["unexplained_by_element"].items())[:12]:
        print(f"    {e:4s} {n}")
    print(f"\nRECOVERY CANDIDATES: {len(rec)}")
    for e, n in summary["recovery_by_element"].items():
        print(f"    {e:4s} {n}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
