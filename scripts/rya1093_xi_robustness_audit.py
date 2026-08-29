#!/usr/bin/env python3
"""RYA-1093 — does the 58-line pool constrain solar xi at all, or is 0.709 a fit artifact?

    python3 scripts/rya1093_xi_robustness_audit.py --star solar --element Fe --ion I

THE QUESTION, AND WHY THE SPREAD IS THE RED FLAG
------------------------------------------------
RYA-311 measured xi = 0.7088 +/- 0.0588 (FITEXY) on 58 Fe I lines against a pinned 1.0. The
same 58 lines give 0.909 under OLS. A 0.3 km/s spread decided by FIT METHOD, with a formal
error a fifth of it, is not a measurement — it is a pool that does not pin the quantity, and
this audit's job is to say so with numbers rather than adopt the prettiest of three.

🔴 FINDING 1 — THE "FITEXY vs OLS" GAP IS **WEIGHTING**, NOT BOTH-AXIS ERRORS.
This is the ticket's own premise, and it does not survive measurement. Mucciarelli's
argument is that sigma_REW belongs in the fit; run it and the abscissa errors turn out to be
irrelevant here. At the solution, per line:

    median sigma_y^2        = 8.70e-04
    median b^2 sigma_x^2    = 4.42e-09        <- the entire both-axis correction
    ratio                   = 5.1e-06

Driving sigma_x to zero (which makes FITEXY exactly a WEIGHTED OLS) moves the slope by
5.7e-04 dex/dex. The gap to the production number is 0.089 dex/dex — **157x larger**. The
whole difference is that production's `_compute_rew_slope` is UNWEIGHTED:

    unweighted OLS   +0.0834      <- production (`np.polyfit`, no weights)
    weighted OLS     -0.0053      <- same fit, 1/sigma_A^2
    FITEXY           -0.0058      <- weighted OLS plus a 5.7e-04 both-axis nudge

sigma_A ranges 0.0150-0.0665 across the pool, a factor 4.4, so ignoring it is not a small
approximation. The right production fix is therefore WEIGHTING (and FITEXY is the correct
weighted estimator to adopt) — but a reader who "switches to FITEXY" believing the both-axis
term is what was missing has misdiagnosed a real bug and will not understand why some other
pool does not respond the same way.

🔴 FINDING 2 — THE POOL IS TRUNCATED AT EXACTLY THE LINES THAT PIN xi.
xi is constrained by SATURATED lines: on the flat part of the curve of growth the abundance
depends on microturbulence and little else. The production ceiling cuts EW > 100 mA, and the
pool's maximum EW is 100.0 mA exactly — the distribution is piled against its own cut, with
18 of 58 lines in the topmost REW bin and only 4 above REW = -4.8.

The solve is correspondingly unstable to those few lines:

    all 58 lines            FITEXY -0.0058
    drop REW >= -4.80  (-4) FITEXY -0.0478
    drop REW >= -4.85 (-13) FITEXY +0.0801
    drop REW >= -4.90 (-20) FITEXY +0.1497

A 0.20 dex/dex swing from removing 20 lines at the strong end. The slope is a property of a
handful of near-ceiling points, not of the pool.

🔴 FINDING 3 — EVERY POOL THAT WOULD TEST IT HAS NO SOLUTION.
The two independent checks the ticket asks for both come back empty, and RYA-311 already ran
them:

  * REMOVE THE CEILING (n=78, the saturated lines restored) -> `root_found: false`. Adding
    the deep lines does not reconcile FITEXY with OLS; it destroys the root entirely.
  * GRADED-ONLY (n=9, cleanest laboratory gf) -> `root_found: false` over xi 0.2-3.0.

So the gf-slope hypothesis cannot even be tested on the graded subset — not because it was
refuted, but because 9 lab-gf lines do not solve xi. That is reported as a GAP, never as a
pass (RYA-805: an absence needs a control, and this one has none).

🔴 FINDING 4 — THE FORMAL ERROR IS NOT THE HONEST ERROR.
chi2_red = 27.9 at the solution: the stated point errors explain about a fifth of the
observed scatter (A scatters 0.149 dex against a median sigma_A of 0.030). The Delta-chi2 = 1
error assumes chi2_red = 1 and is therefore an underestimate by ~sqrt(27.9) = 5.3x, which is
exactly the gap between the 0.0588 formal error and the 0.27-0.31 method+selection spread.
`pipeline/fitexy.py` says this in its own docstring; this audit is where it bites.

THE VERDICT AND THE RULE IT COMES FROM
--------------------------------------
RYA-1093 spec 2B: "If the pool is strong-line-poor / xi under-constrained (formal error stays
<< selection spread even with deep lines): KEEP the pin at 1.0 -- 0.709 is a fit artifact,
not a measurement."

All four findings point the same way, so the pin is KEPT at 1.0 BY THAT FINDING. Ryan's
2026-08-28 amendment makes this the branch that needs no ratification; only an ADOPT would
wait for him, and there is nothing to adopt.

⚠️ THIS AUDIT IS READ-ONLY AND CANNOT MOVE A NUMBER. It re-derives from RYA-311's committed
per-line pool and its committed solves. It writes one report. It does not touch
`config/stars.yaml`, the EW pool, canonical_gf or any product, and nothing in it selects
among the three xi candidates by which one passes a gate (RYA-161).
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline.fitexy import fitexy                                # noqa: E402

RYA311 = ROOT / "data" / "results" / "rya311"
OUT_DEFAULT = ROOT / "data" / "results" / "rya1093"

XI_PINNED = 1.0
#: The REW above which a line is carrying real curve-of-growth saturation, i.e. is actually
#: informative about xi. Mucciarelli 2011 works the solar Fe I solve around this value.
REW_SATURATED = -4.8
#: The 0.05 dex solar gate this feeds (RYA-282). DISPLAY AND COMPARISON ONLY — no xi
#: candidate is chosen by whether it clears it (RYA-161).
SOLAR_GATE_DEX = 0.05


class AuditRefusal(SystemExit):
    """A question this audit cannot answer from committed material."""


def _load(name: str) -> dict:
    p = RYA311 / f"{name}.json"
    if not p.exists():
        raise AuditRefusal(
            f"RYA-311 solve artifact missing: {p}\n"
            f"  This audit re-derives RYA-311's result; it does not re-measure it (the EW "
            f"route needs MOOGSILENT, a Linux binary). Land or stage RYA-311 first.")
    return json.loads(p.read_text())


def per_line() -> pd.DataFrame:
    p = RYA311 / "rya311_xi_solar_FeI_all_ceil100_per_line.csv"
    if not p.exists():
        raise AuditRefusal(f"RYA-311 per-line pool missing: {p}")
    return pd.read_csv(p)


# ── A1: does the pool contain lines that constrain xi? ────────────────────────

def constraint_census(d: pd.DataFrame, solve: dict) -> dict:
    """The strong/weak split, and whether the pool is truncated at its own informative end."""
    rew = d["rew"].to_numpy(float)
    ew = d["ew_mA"].to_numpy(float)
    ceiling = solve.get("ew_ceiling_mA")
    n_sat = int((rew >= REW_SATURATED).sum())
    # ⚠️ THE TELL IS NOT THE COUNT, IT IS THE PILE-UP. A pool can be weak-line-dominated
    # honestly. This one has its MAXIMUM EW sitting exactly on the cut, which means the
    # informative lines were not absent — they were removed.
    at_ceiling = (None if ceiling is None
                  else int(np.isclose(ew, float(ceiling), rtol=0, atol=0.05).sum()))
    hist, edges = np.histogram(rew, bins=8)
    return {
        "n_lines": int(len(d)),
        "rew_min": round(float(rew.min()), 4), "rew_max": round(float(rew.max()), 4),
        "rew_median": round(float(np.median(rew)), 4),
        "ew_max_mA": round(float(ew.max()), 2),
        "ew_ceiling_mA": ceiling,
        "n_at_the_ceiling": at_ceiling,
        "n_saturated_rew_ge_-4.8": n_sat,
        "fraction_saturated": round(n_sat / len(d), 3),
        "n_in_top_rew_bin": int(hist[-1]),
        "histogram": [{"rew_lo": round(float(edges[i]), 3),
                       "rew_hi": round(float(edges[i + 1]), 3), "n": int(hist[i])}
                      for i in range(len(hist))],
        "verdict": (
            f"strong-line-POOR and TRUNCATED: only {n_sat} of {len(d)} lines reach "
            f"REW >= {REW_SATURATED}, and the pool's maximum EW is {ew.max():.1f} mA "
            f"against a {ceiling} mA ceiling — the lines that pin xi were CUT, not absent"
            if ceiling is not None and np.isclose(ew.max(), float(ceiling), atol=0.05)
            else f"{n_sat} of {len(d)} lines at REW >= {REW_SATURATED}"),
    }


# ── A2: what actually separates FITEXY from the production slope ─────────────

def method_decomposition(d: pd.DataFrame) -> dict:
    """🔴 Separate the BOTH-AXIS term from the WEIGHTING term. They are not the same thing.

    The ticket attributes the FITEXY-vs-OLS gap to Mucciarelli's both-axis argument. That is
    testable in one line — run FITEXY with sigma_x driven to zero, which IS a weighted OLS —
    and the answer is that the abscissa errors contribute essentially nothing here while the
    weights contribute everything.

    Getting this right matters beyond the solar number: "switch production to FITEXY" and
    "start weighting the production slope" are different changes, and only the second one is
    what moves this pool's answer.
    """
    x = d["rew"].to_numpy(float)
    y = d["A_at_measured_xi"].to_numpy(float)
    sx = d["sigma_rew"].to_numpy(float)
    sy = d["sigma_A_from_ew"].to_numpy(float)

    f = fitexy(x, y, sx, sy)
    w_only = fitexy(x, y, np.full_like(sx, 1e-12), sy)      # sigma_x -> 0 == weighted OLS
    unweighted = float(np.polyfit(x, y, 1)[0])

    b = float(f.slope)
    y_term = float(np.median(sy ** 2))
    x_term = float(np.median((b ** 2) * (sx ** 2)))
    return {
        "slope_fitexy": round(b, 6),
        "slope_weighted_ols": round(float(w_only.slope), 6),
        "slope_unweighted_ols_PRODUCTION": round(unweighted, 6),
        "sigma_slope_formal": round(float(f.sigma_slope), 6),
        "sigma_slope_chi2_scaled": round(float(f.sigma_slope_scaled), 6),
        "chi2_red": round(float(f.chi2_red), 3),
        "both_axis_term_dex_per_dex": round(abs(b - float(w_only.slope)), 8),
        "weighting_term_dex_per_dex": round(abs(unweighted - float(w_only.slope)), 8),
        "weighting_over_both_axis": round(
            abs(unweighted - float(w_only.slope)) / max(abs(b - float(w_only.slope)), 1e-12), 1),
        "median_sigma_y_squared": y_term,
        "median_b2_sigma_x_squared": x_term,
        "x_over_y_denominator_ratio": float(x_term / y_term),
        "sigma_A_min": round(float(sy.min()), 4), "sigma_A_max": round(float(sy.max()), 4),
        "sigma_A_dynamic_range": round(float(sy.max() / sy.min()), 2),
        "verdict": (
            "the gap is WEIGHTING, not both-axis errors: driving sigma_x to zero (a weighted "
            "OLS) reproduces FITEXY to within the both-axis term, which is "
            f"{abs(b - float(w_only.slope)):.1e} dex/dex, while the unweighted production "
            f"fit sits {abs(unweighted - float(w_only.slope)):.4f} away — "
            f"{abs(unweighted - float(w_only.slope)) / max(abs(b - float(w_only.slope)), 1e-12):.0f}x "
            "larger. sigma_A spans a factor "
            f"{sy.max() / sy.min():.1f} across the pool, so ignoring it is not a small "
            "approximation."),
    }


# ── A3: how much of the slope is a few near-ceiling lines? ───────────────────

def leverage(d: pd.DataFrame) -> dict:
    """Strip the strong end progressively. A constrained solve barely moves; this one flips.

    ⚠️ This is the test that decides the ticket, and it is deliberately a SENSITIVITY rather
    than a fit statistic: chi2 and the formal error both describe the fit you ran, and say
    nothing about how much of it rests on four points.
    """
    x = d["rew"].to_numpy(float)
    y = d["A_at_measured_xi"].to_numpy(float)
    sx = d["sigma_rew"].to_numpy(float)
    sy = d["sigma_A_from_ew"].to_numpy(float)
    rows = [{"cut": None, "n": int(len(x)),
             "slope_fitexy": round(float(fitexy(x, y, sx, sy).slope), 5),
             "slope_unweighted_ols": round(float(np.polyfit(x, y, 1)[0]), 5)}]
    for cut in (-4.80, -4.85, -4.90, -5.00):
        m = x < cut
        if m.sum() < 10:
            continue
        rows.append({"cut": cut, "n": int(m.sum()),
                     "slope_fitexy": round(float(fitexy(x[m], y[m], sx[m], sy[m]).slope), 5),
                     "slope_unweighted_ols": round(float(np.polyfit(x[m], y[m], 1)[0]), 5)})
    span = max(r["slope_fitexy"] for r in rows) - min(r["slope_fitexy"] for r in rows)
    return {
        "trials": rows, "fitexy_slope_span": round(span, 5),
        "verdict": (
            f"the FITEXY slope moves {span:.3f} dex/dex as the strong end is stripped — the "
            f"slope is a property of a handful of near-ceiling lines, not of the pool. A "
            f"solve this sensitive to its own truncation does not constrain xi."),
    }


# ── A4: the checks that have no answer at all ────────────────────────────────

def missing_solves() -> dict:
    """The deep-line and graded-only solves — both already run by RYA-311, both rootless.

    Reported as a GAP, loudly. "The gf-slope hypothesis was not confirmed" and "the gf-slope
    hypothesis could not be tested" are different statements and only the second is true.
    """
    noceil = _load("rya311_xi_solar_FeI_all_noceil")
    graded = _load("rya311_xi_solar_FeI_graded_ceil100_xi0.2-3")
    return {
        "deep_lines_ceiling_removed": {
            "artifact": "rya311_xi_solar_FeI_all_noceil.json",
            "n_fitted": noceil.get("n_fitted"), "root_found": bool(noceil.get("root_found")),
            "xi_measured": noceil.get("xi_measured"),
            "verdict": ("adding the saturated lines does NOT reconcile FITEXY with OLS — it "
                        "removes the root altogether. Spec 2B's ADOPT branch requires a "
                        "reconciliation that does not exist."),
        },
        "graded_only_lab_gf": {
            "artifact": "rya311_xi_solar_FeI_graded_ceil100_xi0.2-3.json",
            "n_fitted": graded.get("n_fitted"), "root_found": bool(graded.get("root_found")),
            "xi_measured": graded.get("xi_measured"),
            "verdict": ("the gf-slope hypothesis CANNOT BE TESTED on this pool: 9 lab-gf "
                        "lines give no root over xi 0.2-3.0. Reported as a GAP, not as a "
                        "clean bill of health (RYA-805). Ryan's 2026-08-28 amendment asked "
                        "for graded-only vs full — this is the graded-only answer."),
        },
    }


# ── the honest delta_xi and the gate ─────────────────────────────────────────

def honest_delta_xi(solve: dict, meth: dict, lev: dict) -> dict:
    """WHICH error bar is honest, argued from the audit rather than chosen.

    🔴 RYA-161: the candidate is NOT selected by which one clears the 0.05 gate. The rule is
    the ticket's own — the formal error is honest only if the solve is well constrained — and
    this audit has just shown it is not, on four independent grounds. So the honest bar is
    the method+selection spread, and it is the larger number. That it then FAILS the gate is
    a result, not a reason to revisit the choice.
    """
    formal = float(solve["delta_xi_formal"])
    chi2s = float(solve["delta_xi_chi2_scaled"])
    xi_fitexy = float(solve["xi_measured"])
    xi_ols = float(solve["xi_ols_root"])
    method_spread = abs(xi_ols - xi_fitexy)
    pin_spread = abs(XI_PINNED - xi_fitexy)
    selection = max(method_spread, pin_spread)
    return {
        "candidates": {
            "formal_delta_chi2_1": round(formal, 4),
            "chi2_scaled": round(chi2s, 4),
            "method_spread_fitexy_vs_ols": round(method_spread, 4),
            "pin_spread_fitexy_vs_1.0": round(pin_spread, 4),
        },
        "chosen": round(selection, 4),
        "chosen_name": "method+selection spread",
        "why": (
            f"the formal error would be honest only if the solve were well constrained. It "
            f"is not: chi2_red = {meth['chi2_red']} (the point errors explain ~1/"
            f"{round(meth['chi2_red'] ** 0.5)} of the scatter), the slope moves "
            f"{lev['fitexy_slope_span']:.3f} dex/dex when the strong end is stripped, the "
            f"ceiling-free and graded-only solves have no root at all, and the pool's "
            f"maximum EW sits exactly on its own cut. The formal "
            f"{formal:.4f} is {selection / formal:.1f}x smaller than the spread the method "
            f"and selection actually produce."),
        # ⚠️ THE GATE CONSTANT IS DELIBERATELY NOT REFERENCED IN THIS FUNCTION.
        # `test_nothing_selects_a_candidate_by_the_gate` reads this function's AST and
        # refuses any mention of `SOLAR_GATE_DEX` inside it — even in a message. A
        # function that chooses an error bar must not be able to see the threshold that
        # bar will be judged against, and "it was only in a string" is exactly the
        # argument that would let the reference back in (RYA-161).
        "explicitly_not_chosen": (
            f"{formal:.4f} (the formal error) is NOT chosen, and specifically not because "
            f"it would clear the solar gate. Choosing an error bar because it passes is "
            f"the RYA-161 failure this ticket names as CRITICAL; the gate is applied in "
            f"`gate_report`, which cannot influence this choice."),
    }


def gate_report(delta_xi: dict, solve: dict) -> dict:
    """sigma_reported against the 0.05 solar gate, with the honest bar. No tuning.

    ⚠️ The gate is REPORTED, not satisfied. If the honest xi allowance does not fit inside
    0.05 dex, the finding is about the gate — and the field's own solar Fe floor is 0.04-0.06
    dex, so a gate at 0.05 with an honestly budgeted xi term is tight by construction.
    """
    d_xi = float(delta_xi["chosen"])
    # dA/dxi from RYA-311's own sweep: the A shift between the measured xi and the pin,
    # divided by the xi difference. Read from the artifact rather than re-modelled.
    dA = float(solve["A_shift_measured_minus_pinned"])
    dxi = float(solve["xi_measured"]) - XI_PINNED
    sens = abs(dA / dxi) if dxi else float("nan")
    contrib = sens * d_xi
    # ⚠️ SAY WHAT THIS NUMBER IS. The chosen delta_xi is the FITEXY-to-pin spread, and the
    # sensitivity is measured over that same interval, so `contrib` reduces to |A(FITEXY) -
    # A(pin)| exactly. That is not a coincidence to hide behind two multiplications: the xi
    # allowance spans the two candidates, so its abundance cost IS the gap between them.
    # Stated, because a reader who does not notice would take it for an independent estimate.
    # Tolerance is 1e-4, not 1e-6, and the reason is arithmetic rather than slack:
    # `chosen` is published rounded to 4 dp (0.2912) while the sensitivity is computed
    # from the unrounded interval (0.29124...), so the identity holds to ~9e-6 and no
    # tighter. Stating the bound rather than loosening it silently.
    reduces_to_a_shift = bool(abs(contrib - abs(dA)) < 1e-4)
    return {
        "dA_dxi_dex_per_kms": round(sens, 4),
        "delta_xi_used": round(d_xi, 4),
        "xi_contribution_to_sigma_params_dex": round(contrib, 4),
        "solar_gate_dex": SOLAR_GATE_DEX,
        "xi_term_alone_exceeds_the_gate": bool(contrib > SOLAR_GATE_DEX),
        "reduces_to_the_A_shift_between_candidates": reduces_to_a_shift,
        "how_it_is_built": (
            "delta_xi is the FITEXY-to-pin spread and dA/dxi is measured over that same "
            "interval, so this term is identically |A(FITEXY) - A(pin)| = "
            f"{abs(dA):.4f} dex. It is the abundance cost of not knowing which of the two "
            "candidates is right — which is exactly what an honest xi allowance means when "
            "the pool cannot choose between them."),
        "verdict": (
            f"the xi term ALONE contributes {contrib:.4f} dex, against a {SOLAR_GATE_DEX} "
            f"dex gate — so with an honestly budgeted xi the solar gate cannot be met by "
            f"this pool, before any other systematic is added. Reported as a finding about "
            f"the GATE (RYA-161); the adopt-or-relax call is Ryan's."
            if contrib > SOLAR_GATE_DEX else
            f"the xi term contributes {contrib:.4f} dex, inside the {SOLAR_GATE_DEX} gate."),
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--star", default="solar", choices=["solar"])
    ap.add_argument("--element", default="Fe")
    ap.add_argument("--ion", default="I")
    ap.add_argument("--out", type=Path, default=OUT_DEFAULT)
    args = ap.parse_args(argv)

    solve = _load("rya311_xi_solar_FeI_all_ceil100")
    d = per_line()
    census = constraint_census(d, solve)
    meth = method_decomposition(d)
    lev = leverage(d)
    gaps = missing_solves()
    dxi = honest_delta_xi(solve, meth, lev)
    gate = gate_report(dxi, solve)

    constrained = not (
        census["fraction_saturated"] < 0.25
        or lev["fitexy_slope_span"] > 0.05
        or not gaps["deep_lines_ceiling_removed"]["root_found"]
        or meth["chi2_red"] > 4.0)
    verdict = ("ADOPT" if constrained else "KEEP-PIN-1.0")

    rep = {
        "ticket": "RYA-1093",
        "generated_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "read_only": True,
        "star": args.star, "species": f"{args.element} {args.ion}",
        "inputs": {"solve": "data/results/rya311/rya311_xi_solar_FeI_all_ceil100.json",
                   "per_line": "data/results/rya311/"
                               "rya311_xi_solar_FeI_all_ceil100_per_line.csv"},
        "xi_candidates": {"fitexy": solve["xi_measured"],
                          "ols": solve["xi_ols_root"], "pinned": XI_PINNED},
        "A_shift_if_adopted_dex": solve["A_shift_measured_minus_pinned"],
        "constraint_census": census,
        "method_decomposition": meth,
        "leverage": lev,
        "unavailable_solves": gaps,
        "honest_delta_xi": dxi,
        "gate": gate,
        "verdict": verdict,
        "verdict_reason": (
            "the 58-line pool does NOT constrain solar xi. It is truncated at exactly the "
            "saturated lines that pin it (max EW = its own 100 mA ceiling); the slope swings "
            f"{lev['fitexy_slope_span']:.3f} dex/dex when the strong end is stripped; "
            f"chi2_red = {meth['chi2_red']}; and BOTH independent checks — ceiling removed "
            "and graded-only — have no root at all. 0.709 is a fit artifact of the "
            "truncation, not a measurement. The pin is KEPT at 1.0 by RYA-1093 spec 2B, "
            "which Ryan's 2026-08-28 amendment makes the branch needing no ratification."),
        "firewall": ("RYA-161: no xi candidate and no delta_xi was selected by whether it "
                     "clears the 0.05 gate. The honest bar is the larger number and it "
                     "FAILS the gate; that is reported, not resolved."),
    }
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "xi_robustness_audit.json").write_text(json.dumps(rep, indent=2) + "\n")

    print("=" * 90)
    print(f"RYA-1093 — solar xi robustness audit   {args.element} {args.ion}  {args.star}")
    print("=" * 90)
    print(f"\n  candidates: FITEXY {solve['xi_measured']:.4f} | OLS "
          f"{solve['xi_ols_root']:.4f} | PIN {XI_PINNED:.2f}   "
          f"(adopting FITEXY would move A(Fe) by "
          f"{solve['A_shift_measured_minus_pinned']:+.4f} dex)")
    print(f"\nA1  CONSTRAINT — {census['verdict']}")
    print(f"      REW {census['rew_min']:.3f}..{census['rew_max']:.3f}; "
          f"{census['n_saturated_rew_ge_-4.8']}/{census['n_lines']} at REW >= -4.8; "
          f"{census['n_in_top_rew_bin']} of {census['n_lines']} in the TOPMOST REW bin "
          f"(max EW {census['ew_max_mA']} = the {census['ew_ceiling_mA']} mA cut)")
    print(f"\nA2  METHOD — {meth['verdict']}")
    print(f"      unweighted OLS {meth['slope_unweighted_ols_PRODUCTION']:+.5f}  |  "
          f"weighted OLS {meth['slope_weighted_ols']:+.5f}  |  "
          f"FITEXY {meth['slope_fitexy']:+.5f}")
    print(f"      chi2_red {meth['chi2_red']}  ->  formal sigma_slope "
          f"{meth['sigma_slope_formal']:.5f} vs chi2-scaled "
          f"{meth['sigma_slope_chi2_scaled']:.5f}")
    print(f"\nA3  LEVERAGE — {lev['verdict']}")
    for t in lev["trials"]:
        c = "all lines" if t["cut"] is None else f"drop REW >= {t['cut']:.2f}"
        print(f"      {c:<22} n={t['n']:<4} FITEXY {t['slope_fitexy']:+.5f}")
    print("\nA4  THE CHECKS WITH NO ANSWER")
    for k, g in gaps.items():
        print(f"      {k:<28} n={g['n_fitted']}  root_found={g['root_found']}")
        print(f"          {g['verdict']}")
    print(f"\nDELTA_XI — chosen {dxi['chosen']:.4f} ({dxi['chosen_name']})")
    print(f"      candidates: {dxi['candidates']}")
    print(f"      {dxi['why']}")
    print(f"      ⚠️ {dxi['explicitly_not_chosen']}")
    print(f"\nGATE — {gate['verdict']}")
    print(f"\nVERDICT: {verdict}")
    print(f"      {rep['verdict_reason']}")
    print(f"\nwrote {args.out / 'xi_robustness_audit.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
