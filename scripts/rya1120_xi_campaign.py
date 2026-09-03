"""
scripts/rya1120_xi_campaign.py
==============================
RYA-1120 / RYA-282 §2 — measure dA/dxi on EACH PRODUCT'S OWN POOL, then sigma_params.

WHY A CAMPAIGN AND NOT A NUMBER. RYA-1089 measured dA/dvmic = -0.24 dex/(km/s) on ONE
62-line pool. RYA-1093 showed xi-sensitivity is a STRONG-LINE phenomenon, so that
derivative is a property of the LINE SET, not of Fe: a DEEPGRADED pool (feature depth
> 0.60, saturated by construction) and a GRADED pool (<= 0.60) do not share one. The
only honest per-product term is a per-product measurement — RYA-282 §2's
perturb-and-re-derive, run on the pool that product actually used.

THE UNIT OF WORK IS NOT THE PRODUCT. One `derive_band_products` invocation emits every
treatment for a given (ion x holding x tier) at once, so 42 xi-applicable products
collapse to ~11 RUN UNITS. The campaign perturbs each unit, not each product.

🔴 THREE THINGS THAT WOULD MAKE THIS SILENTLY WRONG, AND WHAT STOPS EACH

  1. AN INERT OVERRIDE. `xi` reaches the run through `get_star_params`, and the
     ATMOSPHERE is loaded with it — so an override applied after the context is built
     changes nothing, every difference comes out 0.000, and the campaign reports "xi does
     not matter" with total confidence. `--control` runs one unit at two xi and REFUSES
     to proceed unless the answers actually differ. A guard whose two sides converge is
     not a guard (RYA-853).

  2. COLLIDING ARTIFACTS. xi is NOT in the artifact stem, so two perturbed runs write the
     same filename and the second silently overwrites the first — RYA-1099 lost its xi=0
     product exactly this way and it survives only in a log. Every run here gets its own
     output directory, and the harness refuses to read a directory it did not just write.

  3. A MOVING POOL. 🔴 The control caught this: LINE ACCEPTANCE ITSELF DEPENDS ON xi, so
     the two runs of a pair do not necessarily measure the same lines — ENGINE-B came back
     n=9 at xi=1.00 and n=8 at xi=1.10. Differencing those two AGGREGATES is the RYA-1083
     error wearing a different hat: it is a difference of aggregates over two different
     line sets, not a derivative. So dA/dxi is computed as a PER-LINE PAIRED DIFFERENTIAL
     through `pipeline.paired_differential` — the SAME implementation RYA-1083 already
     built for the NLTE effect, not a second copy of it — and the pool change is recorded.

  4. A PARTIAL DIFFERENCE. If one side of a pair fails, differencing what is left is not
     a smaller measurement, it is a different quantity. A unit missing either side is
     reported UNMEASURED, never half-differenced.

WHAT IT DOES NOT DO. It does not measure Teff. delta_Teff for the Sun is 1.0 K, and at
RYA-1089's measured dA/dTeff = 0.0665 dex/100 K the term is 0.000665 dex — it moves
sigma_params by 3.2e-06, about 30x BELOW the published rounding step, and dA/dTeff would
have to be 4.0x larger than measured to shift even the last published digit. That is
carried as a declared bound with its arithmetic, per the ticket's "or a declared, sourced
bound" — see `--teff-bound-report`. ⚠️ SOLAR ONLY: for any other star delta_Teff is tens
of kelvin and the term is real; the bound must be re-derived per star, never inherited.

Usage
-----
    python scripts/rya1120_xi_campaign.py --control          # prove the override BITES
    python scripts/rya1120_xi_campaign.py --plan             # the run units, no compute
    python scripts/rya1120_xi_campaign.py --teff-bound-report
"""
from __future__ import annotations

import argparse
import json
import math
import os
import pathlib
import subprocess
import sys

REPO = pathlib.Path(__file__).resolve().parents[1]
AUDIT = REPO / "data" / "results" / "rya1112" / "vis_fe_uncertainty_audit.json"

#: The derivative STEP (how far we perturb to get a slope) — `uncertainty_stack`'s own.
#: Deliberately NOT the same thing as delta_p, which is the star's UNCERTAINTY and is what
#: the slope gets multiplied by. Conflating them is a factor-of-three error.
XI_STEP_KMS = 0.10

#: Teff arithmetic, quoted from the measurement rather than re-derived here.
TEFF_DELTA_K = 1.0
TEFF_DA_DTEFF_PER_100K = 0.0665     # RYA-1089, Fe I
DEX_PLACES = 4


def run_units():
    """(ion, holding, tier) units that carry a xi-applicable product. Derived, not typed."""
    doc = json.loads(AUDIT.read_text(encoding="utf-8"))
    units = {}
    for p in doc["products"]:
        if str(p["xi_applicability"]).startswith("NOT APPLICABLE"):
            continue
        units.setdefault((p["ion"], p["holding"], p["tier"]), set()).add(p["n_lines"])
    return {k: sorted(v) for k, v in sorted(units.items())}


def teff_bound_report():
    """Show that the Teff term cannot reach the published precision — the BOUND."""
    xi_term = 0.0699                       # RYA-1089's stamped sigma_B_vmic, as a scale
    t = TEFF_DA_DTEFF_PER_100K * (TEFF_DELTA_K / 100.0)
    step = 10 ** -DEX_PLACES
    need = math.sqrt(max((xi_term + step / 2) ** 2 - xi_term ** 2, 0.0))
    return {
        "delta_Teff_K": TEFF_DELTA_K,
        "dA_dTeff_per_100K_measured": TEFF_DA_DTEFF_PER_100K,
        "teff_term_dex": t,
        "sigma_params_xi_only": xi_term,
        "sigma_params_with_teff": math.hypot(xi_term, t),
        "shift_dex": math.hypot(xi_term, t) - xi_term,
        "published_rounding_step_dex": step,
        "changes_published_bar": round(math.hypot(xi_term, t), DEX_PLACES) != round(xi_term, DEX_PLACES),
        "dA_dTeff_needed_to_shift_last_digit_per_100K": need / (TEFF_DELTA_K / 100.0),
        "factor_larger_than_measured": (need / (TEFF_DELTA_K / 100.0)) / TEFF_DA_DTEFF_PER_100K,
        "scope": ("SOLAR ONLY. delta_Teff = 1.0 K is a property of the Sun; for any other "
                  "star it is tens of kelvin and this term is real. Re-derive per star."),
        "invalidated_if": ("any measured dA/dTeff on a solar Fe pool exceeds the "
                           "'needed' value above — then the bound is void and Teff must "
                           "be measured, not bounded."),
    }


def leg_cmd(ion, holding, tier, *, element="Fe", lo=4200.0, hi=6910.0, instrument=None):
    tier_flag = "graded" if tier == "GRADED" else "all"
    inst = instrument or ("kpno_solar_atlas" if "kpno" in holding
                          else "harps" if "harps" in holding else "iag_fts_solar_atlas")
    return [sys.executable, "scripts/derive_band_products.py",
            "--element", element, "--ion", ion,
            "--lo", str(lo), "--hi", str(hi),
            "--instrument", inst, "--holding", holding,
            "--lines-tier", tier_flag]


def _read_products(outdir):
    """Every product CSV in a run's OWN output dir, keyed by treatment."""
    import csv
    out = {}
    for f in sorted(pathlib.Path(outdir).glob("*_products.csv")):
        with open(f, newline="", encoding="utf-8") as fh:
            for r in csv.DictReader(fh):
                out[(r["treatment"], r["ion"], r["band"])] = float(r["A"])
    return out


def dA_dxi(lines_hi, lines_lo, *, step_kms, minus_dir=None, plus_dir=None,
           xi_nominal=None):
    """dA/dxi as a PER-LINE PAIRED DIFFERENTIAL, divided by the central-difference span.

    Reuses `pipeline.paired_differential` rather than differencing two aggregates. The
    control proved why: the pool is not fixed across xi (ENGINE-B moved n=9 -> n=8), so a
    difference of aggregates compares two different line sets. `difference_of_aggregates`
    is carried alongside so the two can be seen to disagree — the same disclosure
    RYA-1083 makes for the NLTE effect.
    """
    sys.path.insert(0, str(REPO))
    from pipeline.paired_differential import paired_differential

    #: 🔴 RYA-1178 A — THE PAIRING MUST BE PROVEN BY THE ARTIFACTS, NOT BY ISOLATION.
    #: Before this, `lines_hi` and `lines_lo` were correct only because the caller wrote
    #: each leg into its own directory and remembered which was which. The near-UV
    #: worktree collision showed what that is worth: two runs wrote one stem, the second
    #: overwrote the first, and nothing downstream could have noticed. When the run dirs
    #: are supplied, their stamps must form the expected nominal -/+ step pair or this
    #: refuses outright — and the span is then taken FROM THE STAMPS, so the number the
    #: derivative divides by is the number the runs were actually made at.
    pair = None
    if minus_dir is not None or plus_dir is not None:
        from pipeline.xi_pairing import assert_pair, span_kms
        if minus_dir is None or plus_dir is None:
            from pipeline.xi_pairing import XiPairingError
            raise XiPairingError(
                "one run directory was supplied without the other — a derivative needs "
                "both legs; half a pair is a different quantity, not a smaller one.")
        lo_stamp, hi_stamp = assert_pair(minus_dir, plus_dir, xi_nominal=xi_nominal,
                                         step_kms=step_kms)
        pair = {"minus": lo_stamp, "plus": hi_stamp}
        span = span_kms(lo_stamp, hi_stamp)
    else:
        span = 2.0 * float(step_kms)

    pd_ = paired_differential(lines_hi, lines_lo)
    d = pd_.as_dict()
    d.update({
        "dA_dxi_paired": pd_.median / span,
        "dA_dxi_from_aggregates": pd_.difference_of_aggregates / span,
        "pool_moved": bool(len(lines_hi) != len(lines_lo)),
        "n_hi": int(len(lines_hi)), "n_lo": int(len(lines_lo)),
        "step_kms": float(step_kms),
        "span_kms": float(span),
        "xi_pairing": pair,
        "xi_pairing_verified": pair is not None,
    })
    return d


def sigma_xi(dA_dxi_value, delta_xi_kms):
    """|dA/dxi| * delta_xi -- the slope times the star's UNCERTAINTY, not the step."""
    return abs(float(dA_dxi_value)) * float(delta_xi_kms)


def control(xi_lo, xi_hi, workdir, outroot):
    """🔴 PROVE THE OVERRIDE BITES before trusting a single derivative.

    An override that silently does nothing yields dA/dxi = 0 everywhere and a campaign
    that confidently reports "xi does not matter". This runs ONE unit at two xi and
    requires the answers to differ.
    """
    got = {}
    for tag, xi in (("lo", xi_lo), ("hi", xi_hi)):
        od = pathlib.Path(outroot) / f"ctl_{tag}"
        od.mkdir(parents=True, exist_ok=True)
        env = dict(os.environ, CODEX_XI_OVERRIDE=str(xi))
        subprocess.run(leg_cmd("I", "solar_kpno_molecfit_corrected", "GRADED"),
                       cwd=workdir, env=env, check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
        src = pathlib.Path(workdir) / "data" / "results" / "band_products"
        for f in src.glob("*_products.csv"):
            f.replace(od / f.name)
        got[tag] = _read_products(od)

    common = set(got["lo"]) & set(got["hi"])
    if not common:
        raise SystemExit("CONTROL FAILED: the two xi runs share no product — cannot compare.")
    diffs = {k: got["hi"][k] - got["lo"][k] for k in sorted(common)}
    moved = {k: d for k, d in diffs.items() if abs(d) > 0}
    print(f"control: xi {xi_lo} -> {xi_hi} over {len(common)} product(s)")
    for k, d in diffs.items():
        print(f"   {k[0]:<34s} dA = {d:+.4f}")
    if not moved:
        raise SystemExit(
            "🔴 CONTROL FAILED: every product returned the IDENTICAL abundance at two "
            "different xi. The override is INERT — it is being applied after the "
            "atmosphere is loaded, or not at all. Measuring dA/dxi now would report 0.000 "
            "for every pool and look like a result. Fix the hook before running the "
            "campaign.")
    print(f"CONTROL PASSED — {len(moved)}/{len(common)} products moved with xi.")
    return diffs


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--plan", action="store_true")
    ap.add_argument("--teff-bound-report", action="store_true")
    ap.add_argument("--control", action="store_true")
    ap.add_argument("--xi-nominal", type=float, default=1.0)
    ap.add_argument("--workdir", default=".")
    ap.add_argument("--outroot", default="/tmp/rya1120_campaign")
    a = ap.parse_args()

    if a.teff_bound_report:
        print(json.dumps(teff_bound_report(), indent=2))
        return 0
    if a.plan:
        u = run_units()
        print(f"{len(u)} run units ({sum(max(v) for v in u.values())} lines at the widest pool each)")
        for (ion, holding, tier), ns in u.items():
            print(f"  Fe {ion:<3s} {holding:<38s} {tier:<12s} pools={ns}")
        print("\neach unit runs TWICE (xi -/+ step); one run emits every treatment for that pool")
        return 0
    if a.control:
        control(a.xi_nominal - XI_STEP_KMS, a.xi_nominal + XI_STEP_KMS,
                a.workdir, a.outroot)
        return 0
    ap.error("choose --plan, --control or --teff-bound-report")


if __name__ == "__main__":
    raise SystemExit(main())
