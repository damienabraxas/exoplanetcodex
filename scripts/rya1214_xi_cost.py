#!/usr/bin/env python3
"""RYA-1214 Step 5 — what a CNO xi campaign costs, measured rather than guessed.

WHY IT HAS TO BE A NEW CAMPAIGN AND NOT A NUMBER WE ALREADY HAVE
----------------------------------------------------------------
`dA/dxi` is a property of the LINE SET, not of the element (RYA-1093: xi-sensitivity is a
STRONG-LINE phenomenon, which is why a DEEPGRADED pool and a GRADED pool of the same ion
do not share one). Every measured derivative in the repo was taken on an Fe pool. The CNO
indicators are high-excitation permitted lines and two forbidden ones, and the forbidden
pair in particular is a different physical regime — so no Fe number may stand in for
them, in either direction. The honest state today is `UNMEASURED`, and that is what the
products carry.

WHAT THE COST ACTUALLY IS
-------------------------
The unit of work is NOT the product. One `derive_band_products` invocation emits every
treatment for a given (element x ion x band x holding) at once, so the campaign perturbs
RUN UNITS, and each unit needs TWO perturbed runs -- xi = 0.90 and xi = 1.10 -- because
`dA/dxi` is a per-line PAIRED differential (RYA-1120): line acceptance itself depends on
xi, so differencing the two aggregates would difference two different line sets.

    cost  =  2 x (run units)  x  (wall clock per run)

Both factors are read from artifacts here. The unit count comes from the products this
ticket actually emitted, and the wall clock from the pool's own per-cell logs -- not from
an Fe timing, which would be a different line count on a different band.

⚠️ AND EACH PERTURBED RUN NEEDS ITS OWN OUTPUT DIRECTORY. xi is not in the artifact stem,
so two perturbed runs of one unit write the same filename and the second silently
overwrites the first (RYA-1099 lost a product exactly this way; RYA-1168's first near-UV
run was thrown away because four concurrent jobs shared one worktree and one job's
directory came back holding another's line files). Isolation is not an optimisation here,
it is the only thing that makes the pairing TRUE -- nothing downstream can detect a
stolen file, because a stolen file could have come from either xi.
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
OUT = ROOT / "data" / "audit" / "rya1214_cno_products"

#: RYA-1120's own perturbation, reused so the CNO derivative is on the SAME scale as the
#: Fe ones and the two can be compared. Not re-chosen here.
XI_LO, XI_HI = 0.90, 1.10
XI_SPAN_KMS = XI_HI - XI_LO
#: The sourced solar microturbulence uncertainty (RYA-1089), which is what dA/dxi is
#: multiplied by to become sigma_xi.
DELTA_XI_KMS = 0.2912
#: Below this many PAIRED lines a slope is not publishable (RYA-907: unmeasured is not a
#: value). Same floor RYA-1163/1168 apply.
MIN_PAIRED = 3


#: Where this ticket's own pool logs were copied to. ⚠️ NOT `Path.home()/scratch` — the
#: pools ran on SIRIUS and `Path.home()` resolves on whichever box runs this script, so the
#: first version found no logs and reported "NO TIMED CELLS" while the timings existed. A
#: cost estimate that silently has no basis is worse than one that refuses.
TIMING_DIRS = (Path("/tmp/rya1214_timing"),
               Path.home() / "scratch" / "rya1214_logs",
               Path.home() / "scratch" / "rya1214_set_logs")


def measured_run_times(log_dir: Path) -> dict[str, float]:
    """Wall-clock seconds per cell, from the driver logs' own START/OK lines."""
    drv = None
    for cand in (log_dir / "_driver.log", log_dir):
        if cand.is_file():
            drv = cand
            break
    if drv is None:
        return {}
    starts, times = {}, {}
    for line in drv.read_text().splitlines():
        m = re.match(r"(START|OK|FAIL)\s+(\S+)\s+(\d\d):(\d\d):(\d\d)", line.strip())
        if not m:
            continue
        kind, name, h, mi, s = m.group(1), m.group(2), *map(int, m.groups()[2:])
        t = h * 3600 + mi * 60 + s
        if kind == "START":
            starts[name] = t
        elif name in starts:
            dt = t - starts[name]
            times[name] = dt + 86400 if dt < 0 else dt      # midnight rollover
    return times


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    # 🔴 THE PUBLISHED FEEDS ARE THE SOURCE, not a side artifact of the run. A cost
    # estimate read from run outputs would count cells that were never published, and
    # `publish_product.py` is what decides which those are (RYA-1034: "a product is
    # published HERE or it does not exist").
    cells = []
    for el in ("C", "N", "O"):
        feed = ROOT / "data" / "products" / "solar" / f"{el}.json"
        if not feed.exists():
            continue
        doc = json.loads(feed.read_text())
        # 🔴 QUARANTINED PRODUCTS COUNT. Every CNO product is in `quarantine[]` and not
        # `products[]`, on ONE code — NOT_YET_DEFENSIBLE, Ryan's REAL_SPECIES ruling that
        # everything outside {Fe I, Fe II, Al} is a pilot. Reading `products[]` alone would
        # size the campaign at ZERO run units and report "nothing to perturb", which is
        # the answer for a feed with no CNO in it, not for this one. The xi term is owed on
        # a product whichever list it sits in — and it is one of the things standing
        # between these and the live set.
        for sec in ("products", "quarantine"):
            for p in doc.get(sec, []):
                cells.append({"element": p["element"], "ion": p["ion"], "band": p["band"],
                              "holding": p["holding"], "treatment": p["treatment"],
                              "selector": p.get("selector") or "", "feed_section": sec})
    if not cells:
        print("no published C/N/O product — run and publish the pool before sizing its "
              "xi campaign; there is nothing to perturb.", file=sys.stderr)

    # The RUN UNIT includes the SELECTOR: a named-set run and a depth-ranked run on one
    # cell are two pools measured on different lines, and dA/dxi is a property of the LINE
    # SET (RYA-1093). Collapsing them would perturb one and serve the derivative to both.
    #
    # ⚠️ THE SELECTOR IS NORMALISED FIRST, because three historical records carry a
    # mis-parsed one. `derive_band_products` names its artifacts
    # <...>_SYNTH_<SELECTOR>_<TREATMENT>, and an earlier publish of the N ENGINE-A products
    # split that on `_SYNTH_` and kept the whole tail — so they went out with selector
    # `SET-AGSS21_ENGINE-A` while `treatment` already said ENGINE-A. Those are NOT separate
    # pools: an ENGINE-A product comes from the SAME invocation as its 1D-LTE sibling, so
    # counting them separately would size the campaign for three runs nobody needs to make.
    # Normalised with `band_products.TREATMENTS`, the same registry the publisher now uses.
    from pipeline.band_products import TREATMENTS
    def _norm(sel: str) -> str:
        for t in TREATMENTS:
            if sel.endswith("_" + t):
                return sel[: -(len(t) + 1)]
        return sel
    units = sorted({(c["element"], c["ion"], c["band"], c["holding"],
                     _norm(c["selector"])) for c in cells})
    treatments = sorted({c["treatment"] for c in cells})
    # A product whose treatment resolves the velocity field itself has no xi term to
    # measure -- full 3D is the case (RYA-1185), and there is no 3D CNO product here, so
    # the applicable set is every unit. Stated rather than assumed.
    n_units = len(units)

    times = {}
    for d in TIMING_DIRS:
        if d.is_dir():
            for f in sorted(d.glob("*driver*.log")):
                times.update(measured_run_times(f))
        times.update(measured_run_times(d))
    # A cell that FAILED tells us nothing about the cost of a successful perturbed run —
    # the 11 depth-gate refusals took ~60 s each because they never reached a synthesis.
    times = {k: v for k, v in times.items() if v > 120}
    vals = sorted(times.values())
    median_s = vals[len(vals) // 2] if vals else None

    report = {
        "ticket": "RYA-1214", "step": "5 — xi campaign cost for the CNO pools",
        "why_new": "dA/dxi is a property of the LINE SET (RYA-1093), and every measured "
                   "derivative in the repo was taken on an Fe pool. The CNO indicators "
                   "are high-excitation permitted lines plus two forbidden ones; no Fe "
                   "number may stand in for them in either direction.",
        "current_state": "UNMEASURED on every CNO product emitted by this ticket — "
                         "carried as a declared state, never as a borrowed value "
                         "(RYA-907).",
        "perturbation": {"xi_lo_kms": XI_LO, "xi_hi_kms": XI_HI,
                         "span_kms": XI_SPAN_KMS, "delta_xi_kms": DELTA_XI_KMS,
                         "min_paired_lines": MIN_PAIRED},
        "run_units": [{"element": e, "ion": i, "band": b, "holding": h, "selector": s}
                      for e, i, b, h, s in units],
        "n_run_units": n_units,
        "treatments_covered_per_unit": treatments,
        "runs_required": 2 * n_units,
        "measured_wall_clock_s": {
            "n_cells_timed": len(times),
            "median_s": median_s,
            "min_s": vals[0] if vals else None,
            "max_s": vals[-1] if vals else None,
            "basis": "this ticket's own pool logs — same bands, same line counts, same "
                     "box, two jobs in parallel against a box already running another "
                     "campaign. NOT an Fe timing.",
        },
        "estimated_cost": (
            None if median_s is None else {
                "serial_core_hours": round(2 * n_units * median_s / 3600.0, 2),
                "wall_clock_hours_at_P2": round(2 * n_units * median_s / 3600.0 / 2, 2),
                "note": "Sirius has 4 cores and CI shares them, so P2 is the sustainable "
                        "parallelism; P4 halves the wall clock and starves CI.",
            }),
        "isolation_requirement": (
            "EACH perturbed run needs its own output directory. xi is not in the artifact "
            "stem, so two runs of one unit collide on the filename and the second wins "
            "silently (RYA-1099). Nothing downstream can detect it, because a stolen file "
            "could have come from either xi — isolation is what makes the pairing true, "
            "not a speed-up (RYA-1168)."),
        "method": "per-line PAIRED differential via pipeline.paired_differential, never a "
                  "difference of aggregates: line acceptance depends on xi, so the two "
                  "sides of a pair do not necessarily measure the same lines (RYA-1083).",
        "what_it_buys": (
            "sigma_xi = |dA/dxi| x 0.2912 dex per product. On the Fe pools that term runs "
            "0.03-0.09 dex, comparable to or larger than the gf term, so a CNO budget "
            "without it is not 'the same standard as Fe' — it is short by its largest "
            "missing piece."),
        "what_it_does_not_cover": (
            "Teff and logg. RYA-1120's solar Teff bound (delta_Teff = 1.0 K) is SOLAR and "
            "was derived on an Fe pool's dA/dTeff; it must be re-derived for CNO before "
            "it can be cited, not inherited."),
    }
    (OUT / "xi_cost_estimate.json").write_text(json.dumps(report, indent=2) + "\n")

    print("=== RYA-1214 Step 5 — CNO xi campaign cost ===")
    print(f"  run units            : {n_units}   (element x ion x band x holding)")
    print(f"  perturbed runs       : {2 * n_units}   (xi={XI_LO} and xi={XI_HI} per unit)")
    if median_s:
        print(f"  measured wall clock  : median {median_s / 60:.1f} min/run "
              f"(min {vals[0] / 60:.1f}, max {vals[-1] / 60:.1f}, n={len(times)})")
        print(f"  estimated cost       : "
              f"{report['estimated_cost']['serial_core_hours']:.1f} serial core-hours, "
              f"~{report['estimated_cost']['wall_clock_hours_at_P2']:.1f} h wall clock at P2")
    else:
        print("  measured wall clock  : NO TIMED CELLS — run the pool first; this script "
              "refuses to quote an Fe timing for a CNO band")
    for u in units:
        print(f"    {u[0]} {u[1]:3s}  {u[2]:12s} {u[3]:32s} {u[4]}")
    print(f"\n  wrote {(OUT / 'xi_cost_estimate.json').relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
