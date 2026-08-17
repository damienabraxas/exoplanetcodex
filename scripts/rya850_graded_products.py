#!/usr/bin/env python3
"""
RYA-850 — promote the graded (lab-gf) pool to the reported value
================================================================
Two products per cell, GRADED primary, UNGRADED kept and labelled.

WHAT WAS ACTUALLY MISSING. Nothing here is a new measurement. The lab-gf pools were built
by RYA-824 (VIS/IR) and RYA-836 (near-UV), `error_budget.gf_term(graded=...)` already
carries both terms (0.041 vs 0.170), and the near-UV lab-gf cell already reaches the
matrix. What never happened is the wiring: `derive_band_products.py:770` passes
`gf_graded=False` UNCONDITIONALLY, so every EW-route cell is charged the Kurucz 0.17
regardless of what its lines are, and RYA-824's two pools were never emitted as cells at
all. This script closes that, and reports the pair.

🔴 ONLY THREE OF EIGHTEEN CELLS CAN HAVE A GRADED TWIN TODAY, and that is a finding, not a
shortfall of this script. A lab-gf pool exists for Fe I in the near-UV (RYA-836, n=59),
Fe I VIS (RYA-824, n=9) and Fe I red-optical (RYA-824, n=20). There is NO primary-lab-gf
pool for Fe II in any band, nor for any ENGINE-A/ENGINE-B product — nobody has measured
one. Those cells fall back to UNGRADED and say so, rather than being relabelled graded on
the strength of a term they are not entitled to.

⚠️ THE 0.041 TERM UNDERSTATES THESE POOLS. `GRADED_GF_SYSTEMATIC_DEX` is a generic NIST
grade-B bound. RYA-824 recorded what its lines' CITED laboratory sigmas actually are:
0.0524 dex (IR) and 0.0600 dex (VIS) — 28% and 46% larger. The spec asks for the generic
term, so that is what is wired; both are reported side by side so the choice is visible,
and the measured figure is the more defensible one.

⚠️ NO COMBINED HEADLINE IS COMPUTED. The ticket asks for a "combined/headline Fe value".
RYA-712 forbids a cross-band combined product and `pipeline.band_products` deliberately has
no `combine()` — the per-band values are separate measurements whose spread IS the result.
So this emits the per-band graded table and names the candidate headline cell; collapsing
them into one number is a reporting decision for a human, not arithmetic to do quietly here.

⚠️ THE LINE SETS ARE PRE-RYA-847. RYA-847 will remove unconstrained synthesis fits, and the
graded pool is then lab-gf lines INTERSECT constrained lines. It is unmerged and in flight,
so every number here is provisional and must be regenerated once it lands. The near-UV
cells are the exposed ones (synthesis route); the VIS/IR lab-gf pools are EW-route and are
outside 847's scope.

Usage:
    python3 scripts/rya850_graded_products.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline.band_policy import resolve as resolve_band          # noqa: E402
from pipeline.error_budget import build as build_budget           # noqa: E402
from pipeline.graded_reporting import (                           # noqa: E402
    annotate, element_table, format_value, pair_products, stat_from_scatter,
    total_sigma)

MATRIX = ROOT / "data" / "results" / "rya783" / "fe_product_matrix.csv"
RYA824 = ROOT / "data" / "results" / "rya824" / "rya824_lab_gf_per_line.csv"
OUT = ROOT / "data" / "results" / "rya850"

#: The EW/profile-fit handler's MEASURED residual, as `derive_band_products` charges it.
#: Kept identical so a graded cell differs from its ungraded twin in the gf term ALONE.
PROFILE_FIT_RESIDUAL_DEX = 0.0129

#: The pool's OWN cited laboratory sigma, derived from the lines it actually contains
#: rather than hardcoded. RYA-853 refereed this table against the source papers and found
#: Ruffoni 142/142 and Den Hartog 203/203 PERFECT on both value and cited sigma, so these
#: numbers are audited data now, not a plausible-looking alternative.
LAB_CSV = ROOT / "data" / "reference" / "fe_gf_lab" / "fe1_lab_loggf.csv"
RYA836_PER_LINE = (ROOT / "data" / "results" / "rya836"
                   / "rya836_nearuv_lab_gf_per_line.csv")
CITED_MATCH_TOL_A = 0.05


def cited_sigma(waves) -> tuple[int, float]:
    """RMS of the per-line cited laboratory sigma over a pool.

    RMS rather than the median because these combine in quadrature into the budget, and
    the median throws away exactly the tail that a quadrature sum is sensitive to: the
    IR pool's median is 0.030 while its RMS is 0.0524. Reproduces RYA-824's published
    0.0524 / 0.0600 exactly, which is the check that this is the same quantity 824 meant.
    """
    lab = pd.read_csv(LAB_CSV)
    got = []
    for w in waves:
        m = lab[np.abs(lab.wavelength_air_A - float(w)) <= CITED_MATCH_TOL_A]
        if len(m):
            got.append(float(m.iloc[0].e_loggf_dex))
    if not got:
        return 0, float("nan")
    a = np.asarray(got, dtype=float)
    return len(a), float(np.sqrt((a ** 2).mean()))

#: 824 named its bands by instrument regime; the matrix names them by band policy.
BAND_FROM_824 = {"IR": "red-optical", "VIS": "VIS"}


def graded_cells_from_824() -> pd.DataFrame:
    """Emit RYA-824's two lab-gf pools as matrix cells, on the graded gf term.

    The pools were measured and published in RYA-824 and then never became cells. This
    re-aggregates the SAME per-line abundances — no re-inversion — and charges the budget
    the graded term, which is the entire difference from the ungraded twin.
    """
    d = pd.read_csv(RYA824)
    rows = []
    for band824, grp in d.groupby("band"):
        band = BAND_FROM_824.get(str(band824))
        if band is None:
            continue
        ok = grp[np.isfinite(grp.a_labgf)]
        if len(ok) < 3:
            continue
        value = float(ok.a_labgf.median())
        scatter = float(ok.a_labgf.std(ddof=1))
        n = int(len(ok))

        pol = resolve_band(float(ok.wavelength_air_A.median()))
        b = build_budget("Fe", float(ok.wavelength_air_A.median()), n,
                         scatter_dex=scatter, gf_graded=True,
                         harness_residual_dex=PROFILE_FIT_RESIDUAL_DEX,
                         handler="ProfileFitHandler")
        stat, syst = b.total()
        rows.append({
            "element": "Fe", "ion": str(ok.ion.iloc[0]), "band": pol.name,
            "instrument": "kpno_solar_atlas", "treatment": "1D-LTE-LABGF",
            "A": round(value, 3), "n_lines": n, "n_excluded": int(len(grp) - len(ok)),
            "stat_dex": round(float(stat), 4), "syst_dex": round(float(syst), 4),
            "dominant": (b.dominant().name if b.dominant() else ""),
            "_src": "rya824 lab-gf pool (RYA-850 wiring)",
            "scatter_dex": round(scatter, 4),
        })
    return pd.DataFrame(rows)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    if not MATRIX.exists():
        raise SystemExit(f"Fe matrix absent: {MATRIX}")

    matrix = pd.read_csv(MATRIX)
    new = graded_cells_from_824()
    print(f"[wire] RYA-824 lab-gf pools -> {len(new)} new graded cells")
    for _, r in new.iterrows():
        print(f"       {r.band:<12} {r.ion:<3} n={r.n_lines:<3} "
              f"A={r.A:.3f}  stat={r.stat_dex:.4f} syst={r.syst_dex:.4f} "
              f"(scatter {r.scatter_dex:.4f})")

    combined = pd.concat([matrix, new.drop(columns=["scatter_dex"])], ignore_index=True)
    d = annotate(combined)

    # ── the pairing ───────────────────────────────────────────────────────────────
    pairs = pair_products(combined)
    n_graded = sum(1 for p in pairs if p.primary_is_graded)
    print(f"\n=== {len(pairs)} cells: {n_graded} GRADED, {len(pairs)-n_graded} "
          f"UNGRADED (no lab-gf pool measured) ===")
    print(f"{'band':<13}{'ion':<4}{'engine':<16}{'pool':<10}"
          f"{'A':>8}{'stat':>8}{'syst':>8}{'total':>8}   vs ungraded")
    for p in pairs:
        prim = p.primary
        if prim is None:
            continue
        pool = "GRADED" if p.primary_is_graded else "ungraded"
        beats = p.graded_beats_ungraded
        cmp_txt = ""
        if beats is not None:
            u = p.ungraded["total_dex"]
            g = p.graded["total_dex"]
            cmp_txt = (f"  {u:.4f} -> {g:.4f}  "
                       f"{'TIGHTER' if beats else '⚠️ WIDER — do not promote'}")
        print(f"{p.band:<13}{p.ion:<4}{p.base:<16}{pool:<10}"
              f"{prim['A']:>8.3f}{prim['stat_dex']:>8.4f}{prim['syst_dex']:>8.4f}"
              f"{prim['total_dex']:>8.4f}{cmp_txt}")

    # ── the main-table view ───────────────────────────────────────────────────────
    et = element_table(combined)
    print(f"\n=== main Solar table view (graded where it exists, +/- VISIBLE) ===")
    for _, r in et[et.primary_is_graded].iterrows():
        print(f"  A({r.element} {r.ion}; {r.band}, {r.engine}) = "
              f"{format_value(r.A, r.total_dex)}   "
              f"[stat {r.stat_dex:.3f}, sys {r.syst_dex:.3f}, n={int(r.n_lines)}]")

    # ── the generic term vs what RYA-824 actually cited ───────────────────────────
    # ── the generic term against the pools' own AUDITED cited sigma ───────────────
    per_line_waves = {}
    d824 = pd.read_csv(RYA824)
    per_line_waves["red-optical"] = d824[d824.band == "IR"].wavelength_air_A.tolist()
    per_line_waves["VIS"] = d824[d824.band == "VIS"].wavelength_air_A.tolist()
    if RYA836_PER_LINE.exists():
        d836 = pd.read_csv(RYA836_PER_LINE)
        per_line_waves["near-UV"] = d836[np.isfinite(
            d836.a_labgf)].wavelength_air_A.tolist()

    cited_by_band = {}
    for band, waves in per_line_waves.items():
        n, sig = cited_sigma(waves)
        if n:
            cited_by_band[band] = {"n": n, "cited_sigma": sig}

    print(f"\n=== ⚠️ the generic graded term vs the pools' OWN cited lab sigma ===")
    print(f"  (RYA-853 refereed this table against the source papers: Ruffoni 142/142 and")
    print(f"   Den Hartog 203/203 PERFECT on value AND cited sigma — this is audited data)")
    for band, c in sorted(cited_by_band.items()):
        print(f"  {band:<12} n={c['n']:3d}  cited sigma {c['cited_sigma']:.4f}  "
              f"vs generic {GRADED_GF_TERM if False else 0.041:.3f}  "
              f"({c['cited_sigma']/0.041:.1f}x)")

    # Every graded cell, not just the two wired from RYA-824 — the near-UV cell comes
    # from the matrix and is entitled to the same treatment.
    cited_cells = []
    for pr in pairs:
        if not pr.primary_is_graded:
            continue
        c = cited_by_band.get(pr.band, {}).get("cited_sigma")
        if c is None:
            continue
        g = pr.graded
        b = build_budget("Fe", {"VIS": 5000.0, "red-optical": 8000.0,
                                "near-UV": 3390.0}.get(pr.band, 5000.0),
                         int(g["n_lines"]), scatter_dex=0.4, gf_graded=True,
                         harness_residual_dex=(0.0 if pr.band == "near-UV"
                                               else PROFILE_FIT_RESIDUAL_DEX),
                         handler="SynthesisHandler")
        others = [t.dex for t in b.terms
                  if not t.averages_down and "gf" not in t.name.lower()]
        syst_cited = float(np.sqrt(np.sum(np.square(others)) + c ** 2))
        tot_generic = float(g["total_dex"])
        tot_cited = total_sigma(float(g["stat_dex"]), syst_cited)
        cited_cells.append({
            "band": pr.band, "engine": pr.base, "A": g["A"],
            "cited_sigma": c,
            "syst_generic": float(g["syst_dex"]), "syst_cited": syst_cited,
            "total_generic": tot_generic, "total_cited": tot_cited,
            "pct": 100.0 * (tot_cited / tot_generic - 1.0),
        })
        print(f"  {pr.band:<12} generic 0.041 -> total {tot_generic:.4f}   |   "
              f"cited {c:.4f} -> total {tot_cited:.4f}   ({tot_cited/tot_generic-1:+.1%})")

    # ── do the graded bands agree? ────────────────────────────────────────────────
    # If they disagree by more than their bars, "the graded value" is not well defined and
    # picking one is a choice rather than a measurement. RYA-712 keeps them separate; this
    # says whether that separation is hiding a conflict.
    g = [p for p in pairs if p.primary_is_graded]
    print(f"\n=== do the {len(g)} graded cells agree? (pairwise, in sigma) ===")
    consistency = []
    for i in range(len(g)):
        for j in range(i + 1, len(g)):
            a, b_ = g[i].graded, g[j].graded
            d = float(a["A"]) - float(b_["A"])
            sig = float(np.hypot(a["total_dex"], b_["total_dex"]))
            n_sig = abs(d) / sig if sig > 0 else float("nan")
            consistency.append({"a": f"{g[i].band}", "b": f"{g[j].band}",
                                "delta": d, "combined_sigma": sig, "n_sigma": n_sig})
            print(f"  {g[i].band:<12} vs {g[j].band:<12} "
                  f"delta {d:+.3f}  combined sigma {sig:.3f}  => {n_sig:.1f} sigma"
                  f"{'   ⚠️ DISCREPANT' if n_sig > 2 else ''}")
    spread = (max(float(x.graded["A"]) for x in g)
              - min(float(x.graded["A"]) for x in g)) if g else float("nan")
    print(f"  band-to-band spread of the graded values: {spread:.3f} dex")
    print(f"  (the 1D reference is 7.516 — RYA-669/783; gold's 7.466 is on the 3D scale)")

    combined.to_csv(OUT / "rya850_two_product_matrix.csv", index=False)
    et.to_csv(OUT / "rya850_element_table.csv", index=False)
    summary = {
        "ticket": "RYA-850",
        "n_cells": len(pairs),
        "n_graded_primary": n_graded,
        "n_ungraded_fallback": len(pairs) - n_graded,
        "graded_cells": [
            {"band": p.band, "ion": p.ion, "engine": p.base,
             "A": p.graded["A"], "n_lines": int(p.graded["n_lines"]),
             "stat_dex": p.graded["stat_dex"], "syst_dex": p.graded["syst_dex"],
             "total_dex": p.graded["total_dex"],
             "ungraded_total_dex": (p.ungraded or {}).get("total_dex"),
             "graded_beats_ungraded": p.graded_beats_ungraded}
            for p in pairs if p.primary_is_graded],
        "cited_gf_sigma_by_band": cited_by_band,
        "cells_with_cited_sigma": cited_cells,
        "graded_consistency": consistency,
        "graded_band_spread_dex": spread,
        "caveats": [
            "line sets are PRE-RYA-847; the graded pool is lab-gf INTERSECT constrained "
            "once 847 lands, so these numbers are provisional",
            "no combined cross-band headline is computed — RYA-712 forbids it and "
            "band_products has no combine(); the headline is a SELECTED cell",
            "the generic graded term 0.041 understates RYA-824's own cited lab sigma "
            "(0.0524 IR / 0.0600 VIS)",
            "15 of 18 cells have NO lab-gf pool measured and fall back to ungraded",
        ],
    }
    (OUT / "rya850_summary.json").write_text(json.dumps(summary, indent=2))
    print(f"\n[out] {OUT}/rya850_two_product_matrix.csv"
          f"\n[out] {OUT}/rya850_element_table.csv"
          f"\n[out] {OUT}/rya850_summary.json")


if __name__ == "__main__":
    main()
