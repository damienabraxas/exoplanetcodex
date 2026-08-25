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
from config.constants import PIPELINE, codex_path          # noqa: E402
import pipeline.curate_nonfe_pools as cn                   # noqa: E402

OUT = ROOT / "data" / "audit" / "stage2_pool_ledger"
#: The June fit output, preserved to the artifact store before the stale copy was
#: cleared (RYA-461 save-before-clean). Without it stage 2 is unreconstructable.
STAGING_FALLBACK = codex_path('data.solar_ew_staging')

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


# ── PROMOTION GATE: if lines get promoted, CHECK gf (ratified 2026-08-08) ───────
# Ryan: "the gf thang has popped up before in NLTE. So I propose a decision node in our
# arch logic here, where if Lines get promoted, check gf."
#
# It is the fourth time gf has turned out to be the real answer, after three different
# framings pointed elsewhere first: Ti's inflated NLTE traced to the model-atom vintage
# (RYA-542/545/546), the RYA-398 graded cull, the RYA-161/697 zero-point cluster, and now
# Al. The pattern is that a line-set change is where a latent gf error becomes visible --
# promotion adds or swaps lines, and a RELATIVE gf error between them shows up instantly
# as line-to-line scatter that measurement error cannot explain.
#
# So promotion is a decision node, not a write. Two things are computed here and reported
# with every candidate, because a promotion decision made without them is the Al mistake
# repeated: the gf GRADE of each line, and whether the set's scatter exceeds what its EW
# errors can account for.
#
# Al is the worked case. Its two recovered lines share EP 3.1427 exactly -- same lower
# level -- and differ only in log gf (-1.886 vs -2.249), both ungraded with K75
# provenance. They give 6.601 and 6.262: a 0.339 dex spread against EW errors of 4.2% and
# 4.7% worth ~0.02-0.04 dex. The mean lands on the literature and means nothing, because
# the spread is larger than the Mg 5528-vs-5711 discordance (0.21 dex) that keeps Mg owed.
GF_UNGRADED = frozenset({"", "nan", "none", "?"})


#: THE adjudicated gf source (RYA-354). `linelist_solar` also carries a `nist_grade`
#: column and reading it here would UNDER-report: the linelist is a VALD pull, and the
#: adjudication that attaches a graded value lives in canonical_gf with its
#: `adjudication_status`. Matching is nearest-within-tolerance, never a rounded key --
#: RYA-703/704, and I made that exact mistake in this ticket's sibling script.
CANONICAL_GF = ROOT / "data" / "linelists" / "canonical_gf.csv"
GF_MATCH_A = 0.03
_gf_cache: dict = {}


def gf_status(ll, el, ion, lam) -> tuple[str, float]:
    """(grade, log_gf) from the adjudicated gf table. Grade absent -> UNGRADED, which
    HOLDS a promotion until adjudicated. It does not block the MEASUREMENT; it blocks
    the claim that the measurement sits on a trusted scale."""
    if "df" not in _gf_cache:
        g = pd.read_csv(CANONICAL_GF, low_memory=False)
        g["_el"] = g.species.astype(str).str.split().str[0]
        _gf_cache["df"] = g
    g = _gf_cache["df"]
    r = g[(g._el == el) & ((g.wavelength_air_A - lam).abs() <= GF_MATCH_A)]
    if r.empty:
        return "NOT_IN_CANONICAL_GF", float("nan")
    best = r.loc[(r.wavelength_air_A - lam).abs().idxmin()]
    grade = str(best.get("nist_grade") or "").strip()
    return ("UNGRADED" if grade.lower() in GF_UNGRADED else grade), float(best.log_gf)


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

    # 🔴 RYA-1033 — membership was tested on `round(wavelength_air_A, 2)`. A fitted line
    # whose wavelength straddles a rounding boundary then reads as ABSENT from the pool and
    # gets a drop_reason invented for it, when it is in the pool and was never dropped. The
    # ledger's whole claim is "every drop is accounted for", so a spurious drop is worse
    # here than anywhere else. Tolerance match, per (element, ion).
    from pipeline import line_match

    fit = fit.reset_index(drop=True)
    in_pool_mask = np.zeros(len(fit), dtype=bool)
    for k, grp in fit.groupby([fit.element.astype(str), fit.ion.astype(str)], sort=False):
        sub = pool[(pool.element.astype(str) == k[0]) & (pool.ion.astype(str) == k[1])]
        if sub.empty:
            continue
        res = line_match.match(grp.wavelength_air_A.to_numpy(float),
                               sub.wavelength_air_A.to_numpy(float))
        in_pool_mask[grp.index.to_numpy()] = res.index >= 0

    rows = []
    for pos, r in fit.iterrows():
        if in_pool_mask[pos]:
            continue
        code, detail = classify(r, W, F, ll, synth)
        d = _depth(W, F, float(r.wavelength_air_A))
        n_bad, worst = contaminants(ll, str(r.element), float(r.wavelength_air_A), d)
        grade, lgf = gf_status(ll, str(r.element), str(r.ion), float(r.wavelength_air_A))
        rows.append(dict(gf_grade=grade, log_gf=lgf,
                         element=r.element, ion=r.ion,
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

    # SYNTHESIS CANDIDATES -- and this class is LARGER than the recovery class, which is
    # the answer to "does anything else need Al's treatment". A line that fits well and
    # is deep enough, whose ONLY disqualifier is a contaminant, is not a bad measurement:
    # it is a measurement the EW rung cannot finish and the synthesis rung can. Al I
    # 6696.020 is the worked example -- Fe I 6696.315 at depth 0.105 against its own
    # 0.228. Counting these separately is what stops "clean recovery" from being read as
    # the whole opportunity.
    syn = led[(led.drop_reason == "UNEXPLAINED")
              & (led.n_contaminants > 0)
              & (led.observed_depth.astype(float) >= 0.02)
              & (led.chi2.astype(float) <= 0.01)].copy()
    syn = syn.sort_values(["element", "wavelength_air_A"])

    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    led.to_csv(out / "stage2_drops.csv", index=False)
    rec.to_csv(out / "recovery_candidates.csv", index=False)
    syn.to_csv(out / "synthesis_candidates.csv", index=False)
    summary = dict(
        ticket="RYA-706", stage="fit -> committed pool",
        staging=str(staging), staging_rows=int(len(fit)),
        pool_rows=int(len(pool)), dropped=int(len(led)),
        n_unexplained=n_unexplained,
        n_recovery_candidates=int(len(rec)),
        n_synthesis_candidates=int(len(syn)),
        synthesis_by_element=dict(Counter(syn.element).most_common()),
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
    print(f"\nRECOVERY CANDIDATES (clean, straight EW recovery): {len(rec)}")
    for e, n in summary["recovery_by_element"].items():
        print(f"    {e:4s} {n}")
    print(f"\nSYNTHESIS CANDIDATES (good fit, blocked ONLY by a blend): {len(syn)}")
    for e, n in summary["synthesis_by_element"].items():
        print(f"    {e:4s} {n}")

    # The promotion gate, reported per element over BOTH candidate classes.
    both = pd.concat([rec, syn])
    print("\nPROMOTION GATE — gf grade of every candidate (ratified 2026-08-08):")
    print(f"  {'el':4s}{'cands':>6s}{'ungraded':>10s}{'graded':>8s}   gate")
    for e in sorted(set(both.element)):
        b = both[both.element == e]
        ung = int((b.gf_grade == "UNGRADED").sum())
        print(f"  {e:4s}{len(b):6d}{ung:10d}{len(b)-ung:8d}   "
              + ("HOLD — adjudicate gf before promoting" if ung else "gf graded, clear"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
