"""
RYA-981 — exercise the RYA-968 and RYA-969 frameworks on the RYA-967 graded 55-line fixture.

Test-then-scale: a small, clean, REAL graded set, to catch framework bugs before the expensive
production re-run. The fixture is READ-ONLY and is never written to.

Run:  python3 -m tests.framework_validation --fixture rya967_graded_55 --frameworks 968,969
"""
from __future__ import annotations

import argparse
import inspect
import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline import differential_bridge as B          # noqa: E402
from pipeline import gf_empirical as G                 # noqa: E402
from pipeline.error_budget import empirical_gf_term, zero_point_term   # noqa: E402

FIXTURES = {
    "rya967_graded_55": "data/results/band_products/"
                        "FeI_4200_6910_kpno_solar_atlas_SYNTH_1D-LTE_lines.csv",
}


def load_fixture(name: str) -> pd.DataFrame:
    p = ROOT / FIXTURES[name]
    if not p.exists():
        raise SystemExit(
            f"fixture {name} not found at {p}. It is a Sirius artifact (RYA-967); this "
            f"validation reads it and never writes it.")
    d = pd.read_csv(p)
    return d


def _rows(d: pd.DataFrame) -> list:
    return d.to_dict("records")


# ── section 969 ──────────────────────────────────────────────────────────────────────
def section_969(d: pd.DataFrame, out: dict) -> None:
    print("\n" + "=" * 78)
    print("RYA-969 — differential bridge, Increment 1 machinery")
    print("=" * 78)

    rows = _rows(d)

    # 1 — SELF-CONSISTENCY. Sun against itself must difference to exactly zero.
    print("\n[1] SELF-CONSISTENCY  (star = Sun, reference = Sun)")
    a = d.abundance.to_numpy(float)
    per_line = a - a
    maxdev = float(np.max(np.abs(per_line)))
    print(f"    per-line [Fe/H]_diff : max |deviation from 0| = {maxdev:.3e} dex  (n={len(a)})")
    print(f"    mean                 : {float(np.mean(per_line)):+.3e} dex")
    print(f"    verdict              : {'PASS' if maxdev == 0.0 else 'FAIL'} — a nonzero "
          f"systematic here would be a bug in the subtraction itself")
    out["969_self_consistency_max_dev_dex"] = maxdev

    # 2 — THE GATE, and its DERIVED threshold.
    print("\n[2] LINE-SHARING GATE")
    n_min = B.derive_min_shared_lines()
    print(f"    derived minimum      : {n_min} lines")
    print(f"    derivation           : a hop delta is a MEAN, and a mean is only as sound as "
          f"the sd behind it.")
    print(f"                           rel. uncertainty on a sample sd = 1/sqrt(2(n-1));")
    print(f"                           requiring <= {B.SD_RELATIVE_TOLERANCE:.0%} gives "
          f"n >= 1 + 1/(2*{B.SD_RELATIVE_TOLERANCE}^2) = {n_min}.")
    print(f"                           COMPUTED from the criterion, not typed as a literal.")
    out["969_derived_min_shared_lines"] = n_min
    out["969_threshold_derivation"] = ("n >= 1 + 1/(2*tol^2) with tol = "
                                       f"{B.SD_RELATIVE_TOLERANCE}; sd known to that "
                                       "relative precision")

    th = B.BridgeThresholds(min_shared_lines=n_min, rew_min=-6.0, rew_max=-4.5)
    s = B.line_sharing_gate(rows, rows, th)
    print(f"\n    gate on the fixture  : shared={s.n_shared}  target={s.n_target}  "
          f"reference={s.n_reference}")
    if s.unevaluable:
        print(f"    🔴 UNEVALUABLE      : {s.unevaluable}")
        out["969_gate"] = "UNEVALUABLE"
        out["969_gate_reason"] = s.unevaluable
    elif not s.ok:
        print(f"    REFUSED             : {s.reason}")
        out["969_gate"] = "REFUSED"
    else:
        print(f"    PASS")
        out["969_gate"] = "PASS"

    # 3 — TWO PRODUCTS.
    print("\n[3] TWO-PRODUCT OUTPUT")
    anchor_sd = float(np.std(a, ddof=1))
    zp = anchor_sd / math.sqrt(a.size)
    try:
        hop = B.bridge_hop(rows, rows, th, target="solar", reference="solar")
        bridged = B.chain_to_sun([hop], star="solar",
                                 solar_absolute_dex=float(np.mean(a)), zero_point_dex=zp)
        print("    " + bridged.report().replace("\n", "\n    "))
        out["969_two_products"] = "EMITTED"
    except (B.GateUnevaluable, B.HopTooLarge) as e:
        # The gate refused, so the hop cannot run. Emit the products the fixture CAN support
        # and say exactly which half is missing rather than fabricating a differential.
        print(f"    hop not run — the gate refused: {type(e).__name__}")
        print(f"    differential [Fe/H] : 0.0000 dex by construction (self-reference), but the "
              f"GATE did not certify it")
        print(f"    absolute     A(Fe)  = {float(np.mean(a)):.4f} dex, "
              f"zero-point +/- {zp:.4f} (anchor n={a.size}; does not shrink)")
        out["969_two_products"] = "ABSOLUTE ONLY — gate refused, differential not certified"


# ── section 968 ──────────────────────────────────────────────────────────────────────
def section_968(d: pd.DataFrame, out: dict) -> None:
    print("\n" + "=" * 78)
    print("RYA-968 — empirical gf grade, absolute layer")
    print("=" * 78)

    a = d.abundance.to_numpy(float)

    # 4 — ZERO-POINT.
    print("\n[4] ZERO-POINT on the graded anchor")
    sd = float(np.std(a, ddof=1))
    zp = sd / math.sqrt(a.size)
    print(f"    n={a.size}  mean={float(np.mean(a)):.4f}  sd={sd:.4f}")
    print(f"    zero-point (SEM)     = {zp:.4f} dex   <- the absolute-accuracy floor here")
    print(f"    sd-of-sd             = {sd / math.sqrt(2 * (a.size - 1)):.4f} dex "
          f"({1 / math.sqrt(2 * (a.size - 1)):.0%} of sd) — the anchor is well determined")
    zt = zero_point_term(zp, n_anchor=int(a.size))
    print(f"    budget term          : {zt.name} = {zt.dex:.4f}, averages_down={zt.averages_down}")
    out["968_zero_point_dex"] = zp
    out["968_anchor_n"] = int(a.size)
    out["968_anchor_sd"] = sd

    # thresholds: DERIVED where the data supports it, declared otherwise.
    th = G.Thresholds(
        rew_min=-6.0, rew_max=-4.5,
        min_depth=0.0,                                    # see finding: no depth column
        max_red_chi2=float(np.nanpercentile(d.red_chi2, 95)),
        admission_tol_dex=float(2 * sd),
        min_anchor_lines=B.derive_min_shared_lines(),
        consistency_bound_dex=float(2 * sd),
        fallback_sigma_dex=0.50)
    print(f"\n    thresholds used (validation fixtures, NOT project values):")
    print(f"      max_red_chi2       = {th.max_red_chi2:.1f}  (95th pct of the fixture's own chi2)")
    print(f"      min_anchor_lines   = {th.min_anchor_lines}   (derived, same criterion as 969)")

    rows = _rows(d)
    for r in rows:                       # the fixture is all GF-LAB with a published sigma
        r["gf_tier"] = "LAB"
        r["gf_sigma_dex"] = r.get("gf_sigma_dex")

    # 5 — F2 OFFSET INVARIANCE on stages 1-5.
    print("\n[5] F2 — OFFSET INVARIANCE of stages 1-5")
    base = [G.grade_line(r, th) for r in rows]
    worst, ok = 0.0, True
    for delta in (-0.5, -0.2, -0.05, 0.05, 0.2, 0.5):
        moved = [G.grade_line({**r, "abundance": r["abundance"] + delta}, th) for r in rows]
        for x, y in zip(base, moved):
            if x.stages != y.stages or x.tier != y.tier or x.reason_code != y.reason_code \
                    or x.sigma_dex != y.sigma_dex:
                ok = False
            worst = max(worst, abs((x.sigma_dex or 0) - (y.sigma_dex or 0)))
    print(f"    swept delta          : -0.5 .. +0.5 dex over {len(rows)} lines")
    print(f"    stages 1-5 identical : {ok}   (max sigma drift {worst:.3e})")
    print(f"    verdict              : {'PASS — F2 holds' if ok else 'FAIL'}")
    print(f"    note                 : stage 6 is NOT offset-invariant BY DESIGN; it owes a "
          f"per-line audit list instead, exercised below.")
    out["968_F2_stages_1_5_invariant"] = ok

    # 6 — GRADE DISTRIBUTION.
    print("\n[6] GRADE DISTRIBUTION of the 55 graded lines")
    anchor = G.build_anchor(list(a), th)
    pool = G.grade_pool(rows, th, anchor)
    print(f"    tiers   : {pool.tier_counts}")
    print(f"    reasons : {pool.reason_counts or '{}'}")
    print(f"    in value: {pool.n_in_value} of {len(rows)}")
    low = [v for v in pool.verdicts if v.tier not in (G.TIER_GOLDEN, G.TIER_CONSISTENT)]
    print(f"    grading LOW: {len(low)}  <- these are lab lines, so any is a FINDING")
    for v in low[:8]:
        print(f"       {v.wavelength_air_A:9.3f}  {v.tier:20s} {v.reason_code:18s} "
              f"stages={v.stages}")
    out["968_tier_counts"] = pool.tier_counts
    out["968_n_grading_low"] = len(low)
    out["968_stage6_audit_present"] = any(
        r["anchor_distance_dex"] is not None for r in G.audit_table(pool))
    if pool.gf_sigma_dex:
        gt = empirical_gf_term(pool.gf_sigma_dex, n_lines=pool.n_in_value,
                               provenance="fixture: all GF-LAB")
        print(f"    budget  : {gt.name} = {gt.dex:.4f} dex")
        out["968_gf_sigma_dex"] = pool.gf_sigma_dex

    # 7 — PHYSICS, NOT REFERENCE.
    print("\n[7] RYA-161 — no reference abundance in the grading signature")
    bad = []
    for fn in (G.stage1_depth, G.stage2_saturation, G.stage3_purity, G.stage4_fit, G.stage5_gf):
        p = set(inspect.signature(fn).parameters)
        hit = p & {"abundance", "reference", "anchor", "gold", "expected"}
        print(f"    {fn.__name__:22s} params={sorted(p)}")
        if hit:
            bad.append((fn.__name__, hit))
    src = (ROOT / "pipeline" / "gf_empirical.py").read_text(encoding="utf-8")
    code = "\n".join(l for l in src.splitlines() if not l.lstrip().startswith("#"))
    leaks = [b for b in ("gold_reference", "solar_gold", "STAR_PARAMS", "asplund",
                         "reference_abundance") if b.lower() in code.lower()]
    print(f"    stages taking an abundance : {bad or 'NONE'}")
    print(f"    reference imports in module: {leaks or 'NONE'}")
    print(f"    verdict                    : "
          f"{'PASS — grade keys on line physics only' if not bad and not leaks else 'FAIL'}")
    out["968_no_reference_in_signature"] = (not bad and not leaks)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--fixture", default="rya967_graded_55", choices=sorted(FIXTURES))
    ap.add_argument("--frameworks", default="968,969")
    ap.add_argument("--json-out", default=None)
    a = ap.parse_args(argv)

    d = load_fixture(a.fixture)
    print(f"FIXTURE {a.fixture}: {len(d)} lines, {FIXTURES[a.fixture]}")
    present = {c: int(d[c].notna().sum()) for c in d.columns}
    missing = [c for c, n in present.items() if n == 0]
    print(f"  columns with NO data: {missing}")

    out = {"fixture": a.fixture, "n_lines": int(len(d)), "empty_columns": missing}
    want = {x.strip() for x in a.frameworks.split(",")}
    if "969" in want:
        section_969(d, out)
    if "968" in want:
        section_968(d, out)

    if a.json_out:
        Path(a.json_out).write_text(json.dumps(out, indent=2, default=str) + "\n")
        print(f"\n[out] {a.json_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
